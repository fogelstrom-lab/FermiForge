module he3_trajectory_point_map
  use he3_kinds, only : rk
  use spinful_state_2d, only : spinful_state_2d_t
  use straight_trajectory_2d, only : straight_trajectory_2d_t
  use trajectory_self_energy_2d, only : trajectory_self_energy_2d_t, &
                                        sample_trajectory_self_energy_2d
  use legacy_riccati_trajectory_2d, only : &
    propagate_legacy_riccati_to_target
  use quasiclassical_propagator, only : quasiclassical_propagator_t, &
                                        reconstruct_legacy_quasiclassical_propagator
  use he3_self_consistency_integrand, only : he3_point_map_accumulator_t, &
                                             accumulate_he3_point_map_contribution
  implicit none
  private

  public :: accumulate_legacy_trajectory_point_map

contains

  subroutine accumulate_legacy_trajectory_point_map( &
      state, trajectory, feedback_scale, spectral_energy, substep_count, &
      boundary_relaxation_distance, direction_weight, pole_weight, &
      temperature, accumulator, normalization_error, succeeded)
    type(spinful_state_2d_t), intent(in) :: state
    type(straight_trajectory_2d_t), intent(in) :: trajectory
    real(rk), intent(in) :: feedback_scale
    complex(rk), intent(in) :: spectral_energy
    integer, intent(in) :: substep_count
    real(rk), intent(in) :: boundary_relaxation_distance
    real(rk), intent(in) :: direction_weight, pole_weight, temperature
    type(he3_point_map_accumulator_t), intent(inout) :: accumulator
    real(rk), intent(out) :: normalization_error
    logical, intent(out) :: succeeded

    type(quasiclassical_propagator_t) :: propagator
    type(trajectory_self_energy_2d_t) :: self_energy
    complex(rk), allocatable :: coherence_from_entry(:, :)
    complex(rk), allocatable :: coherence_from_exit(:, :)
    integer :: target

    normalization_error = huge(1.0_rk)
    succeeded = .false.
    call sample_trajectory_self_energy_2d( &
      state, trajectory, feedback_scale, self_energy)
    call propagate_legacy_riccati_to_target( &
      trajectory, self_energy, spectral_energy, substep_count, &
      boundary_relaxation_distance, coherence_from_entry, coherence_from_exit)

    target = trajectory%target_sample
    call reconstruct_legacy_quasiclassical_propagator( &
      coherence_from_entry(:, target), coherence_from_exit(:, target), &
      propagator, succeeded)
    if (.not. succeeded) return

    normalization_error = propagator%normalization_error()
    call accumulate_he3_point_map_contribution( &
      accumulator, propagator, trajectory%momentum, direction_weight, &
      pole_weight, temperature)
  end subroutine accumulate_legacy_trajectory_point_map

end module he3_trajectory_point_map
