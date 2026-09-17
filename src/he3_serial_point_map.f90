module he3_serial_point_map
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t
  use spinful_state_2d, only : spinful_state_2d_t
  use straight_trajectory_2d, only : straight_trajectory_2d_t, &
                                     build_straight_trajectory_2d
  use trajectory_self_energy_2d, only : trajectory_self_energy_2d_t, &
                                        sample_trajectory_self_energy_2d
  use legacy_riccati_trajectory_2d, only : &
    propagate_legacy_riccati_path_to_target
  use free_vortex_asymptotic_2d, only : free_vortex_endpoint_2d_t, &
    extend_with_free_vortex_asymptotic_2d
  use quasiclassical_propagator, only : quasiclassical_propagator_t, &
                                        reconstruct_legacy_quasiclassical_propagator
  use he3_self_consistency_integrand, only : he3_point_map_accumulator_t, &
                                             accumulate_he3_point_map_contribution, &
                                             finalize_he3_point_map
  use he3_quadrature, only : angular_quadrature_3d_t, ozaki_quadrature_t
  implicit none
  private

  complex(rk), parameter :: imaginary_unit = &
    cmplx(0.0_rk, 1.0_rk, kind=rk)

  type, public :: he3_point_map_diagnostics_t
    integer :: direction_count = 0
    integer :: pole_count = 0
    integer :: attempted_contributions = 0
    integer :: accepted_contributions = 0
    integer :: failed_propagators = 0
    real(rk) :: maximum_normalization_error = 0.0_rk
    real(rk) :: sampling_seconds = 0.0_rk
    real(rk) :: propagation_seconds = 0.0_rk
    real(rk) :: accumulation_seconds = 0.0_rk
  end type he3_point_map_diagnostics_t

  public :: evaluate_serial_he3_point_map

contains

  subroutine evaluate_serial_he3_point_map( &
      mesh, state, origin, angular_quadrature, energy_quadrature, &
      feedback_scale, trajectory_maximum_step, minimum_substep_count, &
      boundary_relaxation_distance, mapped_gap, mapped_current_mean_field, &
      diagnostics, succeeded, free_vortex_endpoint)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    real(rk), intent(in) :: origin(2)
    type(angular_quadrature_3d_t), intent(in) :: angular_quadrature
    type(ozaki_quadrature_t), intent(in) :: energy_quadrature
    real(rk), intent(in) :: feedback_scale, trajectory_maximum_step
    integer, intent(in) :: minimum_substep_count
    real(rk), intent(in) :: boundary_relaxation_distance
    complex(rk), intent(out) :: mapped_gap(3, 3)
    real(rk), intent(out) :: mapped_current_mean_field(3)
    type(he3_point_map_diagnostics_t), intent(out) :: diagnostics
    logical, intent(out) :: succeeded
    type(free_vortex_endpoint_2d_t), intent(in), optional :: &
      free_vortex_endpoint

    type(he3_point_map_accumulator_t) :: accumulator
    type(quasiclassical_propagator_t) :: propagator
    type(straight_trajectory_2d_t) :: trajectory
    type(trajectory_self_energy_2d_t) :: propagation_self_energy, self_energy
    complex(rk), allocatable :: coherence_from_entry(:, :)
    complex(rk), allocatable :: coherence_from_exit(:, :)
    complex(rk) :: spectral_energy
    real(rk), allocatable :: propagation_coordinate(:)
    real(rk) :: finish_time, normalization_error, start_time
    integer :: direction, pole, propagation_target, substep_count
    logical :: reconstruction_succeeded, use_free_vortex_endpoint

    if (.not. mesh%is_valid()) error stop "serial point map needs a valid mesh"
    if (.not. state%is_valid_for(mesh)) &
      error stop "serial point map state does not match its mesh"
    if (.not. angular_quadrature%is_valid()) &
      error stop "serial point map needs a valid angular quadrature"
    if (.not. energy_quadrature%is_valid()) &
      error stop "serial point map needs a valid Ozaki quadrature"
    if (trajectory_maximum_step <= 0.0_rk) &
      error stop "serial point-map trajectory step must be positive"
    if (minimum_substep_count < 1) &
      error stop "serial point map needs at least one Riccati substep"
    if (boundary_relaxation_distance < 0.0_rk) &
      error stop "serial point-map boundary distance cannot be negative"
    use_free_vortex_endpoint = .false.
    if (present(free_vortex_endpoint)) then
      use_free_vortex_endpoint = free_vortex_endpoint%enabled
      if (.not. free_vortex_endpoint%is_valid_for(mesh)) &
        error stop "serial point map has an invalid free-vortex endpoint"
    end if

    diagnostics = he3_point_map_diagnostics_t()
    diagnostics%direction_count = angular_quadrature%direction_count()
    diagnostics%pole_count = energy_quadrature%pole_count()
    call accumulator%reset()

    do direction = 1, angular_quadrature%direction_count()
      call cpu_time(start_time)
      call build_straight_trajectory_2d( &
        mesh, origin, angular_quadrature%momentum(:, direction), &
        trajectory_maximum_step, trajectory)
      call sample_trajectory_self_energy_2d( &
        state, trajectory, feedback_scale, self_energy)
      if (allocated(propagation_coordinate)) deallocate(propagation_coordinate)
      if (use_free_vortex_endpoint) then
        call extend_with_free_vortex_asymptotic_2d( &
          mesh, state, trajectory, self_energy, feedback_scale, &
          free_vortex_endpoint, propagation_coordinate, propagation_target, &
          propagation_self_energy)
      else
        allocate(propagation_coordinate(size(trajectory%path_coordinate)), &
          source=trajectory%path_coordinate)
        propagation_target = trajectory%target_sample
        propagation_self_energy = self_energy
      end if
      call cpu_time(finish_time)
      diagnostics%sampling_seconds = diagnostics%sampling_seconds + &
        max(0.0_rk, finish_time - start_time)
      do pole = 1, energy_quadrature%pole_count()
        diagnostics%attempted_contributions = &
          diagnostics%attempted_contributions + 1
        spectral_energy = imaginary_unit * &
          cmplx(energy_quadrature%pole(pole), 0.0_rk, kind=rk)
        ! new_src/riccati.f90 uses 1+int(abs(en)) substeps for each
        ! trajectory interval.  Retain that energy-dependent floor while
        ! allowing convergence studies to request a larger minimum.
        substep_count = max(minimum_substep_count, &
          1 + int(abs(energy_quadrature%pole(pole))))

        call cpu_time(start_time)
        call propagate_legacy_riccati_path_to_target( &
          propagation_coordinate, propagation_target, &
          propagation_self_energy, spectral_energy, substep_count, &
          boundary_relaxation_distance, coherence_from_entry, &
          coherence_from_exit)
        call reconstruct_legacy_quasiclassical_propagator( &
          coherence_from_entry(:, propagation_target), &
          coherence_from_exit(:, propagation_target), &
          propagator, reconstruction_succeeded)
        call cpu_time(finish_time)
        diagnostics%propagation_seconds = diagnostics%propagation_seconds + &
          max(0.0_rk, finish_time - start_time)

        if (.not. reconstruction_succeeded) then
          diagnostics%failed_propagators = diagnostics%failed_propagators + 1
          cycle
        end if
        normalization_error = propagator%normalization_error()
        diagnostics%maximum_normalization_error = max( &
          diagnostics%maximum_normalization_error, normalization_error)

        call cpu_time(start_time)
        call accumulate_he3_point_map_contribution( &
          accumulator, propagator, angular_quadrature%momentum(:, direction), &
          angular_quadrature%weight(direction), &
          energy_quadrature%residue(pole), energy_quadrature%temperature)
        call cpu_time(finish_time)
        diagnostics%accumulation_seconds = &
          diagnostics%accumulation_seconds + &
          max(0.0_rk, finish_time - start_time)
      end do
    end do

    diagnostics%accepted_contributions = accumulator%contribution_count
    call finalize_he3_point_map( &
      accumulator, energy_quadrature%gap_prefactor, mapped_gap, &
      mapped_current_mean_field)
    succeeded = diagnostics%failed_propagators == 0 .and. &
      diagnostics%accepted_contributions == diagnostics%attempted_contributions
  end subroutine evaluate_serial_he3_point_map

end module he3_serial_point_map
