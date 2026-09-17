program test_he3_serial_point_map
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, make_uniform_cartesian_mesh
  use spinful_state_2d, only : spinful_state_2d_t, allocate_spinful_state_2d
  use straight_trajectory_2d, only : straight_trajectory_2d_t, &
                                     build_straight_trajectory_2d
  use he3_quadrature, only : angular_quadrature_3d_t, ozaki_quadrature_t, &
                             make_legacy_angular_quadrature, &
                             make_ozaki_quadrature, &
                             read_legacy_gauss_table, &
                             read_legacy_ozaki_table
  use he3_self_consistency_integrand, only : he3_point_map_accumulator_t, &
                                             finalize_he3_point_map
  use he3_trajectory_point_map, only : accumulate_legacy_trajectory_point_map
  use he3_serial_point_map, only : he3_point_map_diagnostics_t, &
                                   evaluate_serial_he3_point_map
  use free_vortex_asymptotic_2d, only : free_vortex_endpoint_2d_t
  implicit none

  character(len=1024) :: gauss_table_path, ozaki_table_path

  if (command_argument_count() /= 2) &
    error stop "test_he3_serial_point_map needs Gauss and Ozaki table paths"
  call get_command_argument(1, gauss_table_path)
  call get_command_argument(2, ozaki_table_path)

  call test_serial_loop_matches_expanded_sequence()
  call test_full_legacy_quadrature_on_uniform_b_phase( &
    trim(gauss_table_path), trim(ozaki_table_path))

  print '(a)', "3He serial point-map tests passed"

contains

  subroutine test_serial_loop_matches_expanded_sequence()
    type(angular_quadrature_3d_t) :: angular
    type(cartesian_mesh_2d_t) :: mesh
    type(he3_point_map_accumulator_t) :: expanded
    type(he3_point_map_diagnostics_t) :: diagnostics
    type(ozaki_quadrature_t) :: energy
    type(spinful_state_2d_t) :: state
    type(straight_trajectory_2d_t) :: trajectory
    complex(rk) :: expanded_gap(3, 3), mapped_gap(3, 3), spectral_energy
    real(rk) :: expanded_current(3), mapped_current(3)
    real(rk) :: normalization_error, x, y
    integer :: direction, orbital, point, pole, spin, x_node, y_node
    logical :: contribution_succeeded, map_succeeded

    call make_uniform_cartesian_mesh( &
      -0.4_rk, 0.4_rk, 4, -0.4_rk, 0.4_rk, 4, mesh)
    call allocate_spinful_state_2d(mesh, state)
    do y_node = 0, mesh%y_cell_count()
      y = mesh%y_coordinate(y_node)
      do x_node = 0, mesh%x_cell_count()
        x = mesh%x_coordinate(x_node)
        point = mesh%point_index(x_node, y_node)
        do spin = 1, 3
          do orbital = 1, 3
            state%order_parameter(spin, orbital, point) = cmplx( &
              0.04_rk * real(spin + orbital, rk) + &
                0.012_rk * x - 0.007_rk * y, &
              0.009_rk * real(spin - orbital, rk) - &
                0.005_rk * x + 0.011_rk * y, kind=rk)
          end do
        end do
        state%current_mean_field(:, point) = &
          [0.02_rk * x, -0.015_rk * y, 0.01_rk * (x - y)]
      end do
    end do

    call make_legacy_angular_quadrature( &
      2, [-0.35_rk, 0.35_rk], [1.0_rk, 1.0_rk], angular)
    call make_ozaki_quadrature( &
      0.3_rk, 12, [0.55_rk, 1.25_rk], [0.8_rk, 1.1_rk], energy)

    call evaluate_serial_he3_point_map( &
      mesh, state, [0.03_rk, -0.02_rk], angular, energy, 0.7_rk, &
      0.16_rk, 3, 0.08_rk, mapped_gap, mapped_current, diagnostics, &
      map_succeeded)
    call require(map_succeeded, "serial point map reported a failed propagator")

    call expanded%reset()
    do direction = 1, angular%direction_count()
      call build_straight_trajectory_2d( &
        mesh, [0.03_rk, -0.02_rk], angular%momentum(:, direction), &
        0.16_rk, trajectory)
      do pole = 1, energy%pole_count()
        spectral_energy = cmplx(0.0_rk, energy%pole(pole), kind=rk)
        call accumulate_legacy_trajectory_point_map( &
          state, trajectory, 0.7_rk, spectral_energy, 3, 0.08_rk, &
          angular%weight(direction), energy%residue(pole), &
          energy%temperature, expanded, normalization_error, &
          contribution_succeeded)
        call require(contribution_succeeded, &
                     "expanded point-map contribution failed")
      end do
    end do
    call finalize_he3_point_map( &
      expanded, energy%gap_prefactor, expanded_gap, expanded_current)

    call require(maxval(abs(mapped_gap - expanded_gap)) < 2.0e-14_rk, &
                 "serial point-map gap differs from expanded sequence")
    call require(maxval(abs(mapped_current - expanded_current)) < 2.0e-16_rk, &
                 "serial point-map current differs from expanded sequence")
    call require(diagnostics%direction_count == 4 .and. &
                 diagnostics%pole_count == 2 .and. &
                 diagnostics%attempted_contributions == 8 .and. &
                 diagnostics%accepted_contributions == 8 .and. &
                 diagnostics%failed_propagators == 0, &
                 "serial point-map diagnostic counts are wrong")
    call require(diagnostics%maximum_normalization_error < 1.0e-14_rk, &
                 "serial point-map normalization diagnostic is too large")
    call require(diagnostics%sampling_seconds >= 0.0_rk .and. &
                 diagnostics%propagation_seconds >= 0.0_rk .and. &
                 diagnostics%accumulation_seconds >= 0.0_rk, &
                 "serial point-map timing diagnostics are invalid")
  end subroutine test_serial_loop_matches_expanded_sequence


  subroutine test_full_legacy_quadrature_on_uniform_b_phase( &
      gauss_path, ozaki_path)
    character(len=*), intent(in) :: gauss_path, ozaki_path

    type(angular_quadrature_3d_t) :: angular
    type(cartesian_mesh_2d_t) :: mesh
    type(he3_point_map_diagnostics_t) :: diagnostics
    type(free_vortex_endpoint_2d_t) :: endpoint
    type(ozaki_quadrature_t) :: energy
    type(spinful_state_2d_t) :: state
    complex(rk) :: local_gap(3, 3), mapped_gap(3, 3)
    complex(rk) :: off_diagonal_gap(3, 3)
    real(rk) :: diagonal_scale, local_current(3), mapped_current(3)
    integer :: component, point
    logical :: map_succeeded

    call read_legacy_gauss_table(gauss_path, 48, 11, angular)
    call read_legacy_ozaki_table(ozaki_path, energy)
    call make_uniform_cartesian_mesh( &
      -1.0_rk, 1.0_rk, 6, -1.0_rk, 1.0_rk, 6, mesh)
    call allocate_spinful_state_2d(mesh, state)
    do point = 1, mesh%point_count()
      do component = 1, 3
        state%order_parameter(component, component, point) = &
          cmplx(0.22_rk, 0.0_rk, kind=rk)
      end do
    end do

    call evaluate_serial_he3_point_map( &
      mesh, state, [0.0_rk, 0.0_rk], angular, energy, 0.0_rk, &
      0.25_rk, 2, 0.5_rk, mapped_gap, mapped_current, diagnostics, &
      map_succeeded)
    call require(map_succeeded, &
                 "full legacy quadrature point map reported a failure")
    call require(diagnostics%direction_count == 528 .and. &
                 diagnostics%pole_count == 8 .and. &
                 diagnostics%attempted_contributions == 4224 .and. &
                 diagnostics%accepted_contributions == 4224 .and. &
                 diagnostics%failed_propagators == 0, &
                 "full legacy quadrature diagnostic counts are wrong")
    call require(diagnostics%maximum_normalization_error < 2.0e-14_rk, &
                 "full legacy quadrature normalization error is too large")

    off_diagonal_gap = mapped_gap
    do component = 1, 3
      off_diagonal_gap(component, component) = &
        cmplx(0.0_rk, 0.0_rk, kind=rk)
    end do
    diagonal_scale = max(1.0_rk, maxval(abs([ &
      mapped_gap(1, 1), mapped_gap(2, 2), mapped_gap(3, 3)])))
    call require(maxval(abs(off_diagonal_gap)) < 5.0e-13_rk * diagonal_scale, &
                 "uniform B-phase map produced off-diagonal gap components")
    call require(abs(mapped_gap(1, 1) - mapped_gap(2, 2)) < &
                   5.0e-13_rk * diagonal_scale .and. &
                 abs(mapped_gap(2, 2) - mapped_gap(3, 3)) < &
                   5.0e-13_rk * diagonal_scale, &
                 "uniform B-phase map broke rotational component symmetry")
    call require(maxval(abs(mapped_current)) < 5.0e-13_rk, &
                 "uniform equilibrium B-phase map produced a current")

    local_gap = mapped_gap
    local_current = mapped_current
    endpoint%enabled = .true.
    endpoint%bulk_gap = 0.22_rk
    endpoint%winding = 0.0_rk
    endpoint%outer_radius = 2.0_rk
    endpoint%maximum_step = 0.25_rk
    call evaluate_serial_he3_point_map( &
      mesh, state, [0.0_rk, 0.0_rk], angular, energy, 0.0_rk, &
      0.25_rk, 2, 0.5_rk, mapped_gap, mapped_current, diagnostics, &
      map_succeeded, endpoint)
    call require(map_succeeded, &
      "uniform B-phase map failed with the free-vortex endpoint")
    call require(maxval(abs(mapped_gap - local_gap)) < 5.0e-13_rk .and. &
                 maxval(abs(mapped_current - local_current)) < 5.0e-13_rk, &
      "free-vortex endpoint changes a spatially uniform bulk state")
  end subroutine test_full_legacy_quadrature_on_uniform_b_phase


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_he3_serial_point_map
