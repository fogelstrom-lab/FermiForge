program test_he3_field_map_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, &
                                make_uniform_cartesian_mesh
  use spinful_state_2d, only : spinful_state_2d_t, allocate_spinful_state_2d
  use he3_quadrature, only : angular_quadrature_3d_t, ozaki_quadrature_t, &
                             make_legacy_angular_quadrature, &
                             make_ozaki_quadrature
  use he3_serial_point_map, only : he3_point_map_diagnostics_t, &
                                   evaluate_serial_he3_point_map
  use he3_field_map_2d, only : he3_field_map_diagnostics_t, &
                               evaluate_serial_he3_field_map
  use new_src_iteration_layout_2d, only : &
    new_src_iteration_vector_size_2d
  use legacy_anderson_mixing, only : legacy_anderson_t, anderson_report_t
  use he3_nonlinear_iteration_2d, only : he3_field_residual_report_t, &
                                         compute_he3_field_residual, &
                                         update_he3_state_with_anderson
  implicit none

  call test_complete_field_map()
  call test_residual_and_anderson_adapter()

  print '(a)', "3He field-map and nonlinear-adapter tests passed"

contains

  subroutine test_complete_field_map()
    type(angular_quadrature_3d_t) :: angular
    type(cartesian_mesh_2d_t) :: mesh
    type(he3_field_map_diagnostics_t) :: diagnostics
    type(he3_point_map_diagnostics_t) :: point_diagnostics
    type(ozaki_quadrature_t) :: energy
    type(spinful_state_2d_t) :: mapped, state
    complex(rk) :: expected_gap(3, 3)
    real(rk) :: expected_current(3), x, y
    integer :: orbital, point, spin, x_node, y_node
    logical :: point_succeeded, succeeded

    call make_uniform_cartesian_mesh( &
      -0.3_rk, 0.3_rk, 2, -0.3_rk, 0.3_rk, 2, mesh)
    call allocate_spinful_state_2d(mesh, state)
    do y_node = 0, mesh%y_cell_count()
      y = mesh%y_coordinate(y_node)
      do x_node = 0, mesh%x_cell_count()
        x = mesh%x_coordinate(x_node)
        point = mesh%point_index(x_node, y_node)
        do spin = 1, 3
          do orbital = 1, 3
            state%order_parameter(spin, orbital, point) = cmplx( &
              0.05_rk * real(spin + orbital, rk) + 0.01_rk * x, &
              0.007_rk * real(spin - orbital, rk) - 0.008_rk * y, rk)
          end do
        end do
        state%current_mean_field(:, point) = &
          [0.01_rk * x, -0.02_rk * y, 0.005_rk * (x + y)]
      end do
    end do
    call make_legacy_angular_quadrature( &
      2, [-0.4_rk, 0.4_rk], [1.0_rk, 1.0_rk], angular)
    call make_ozaki_quadrature( &
      0.3_rk, 10, [0.45_rk], [1.0_rk], energy)

    call evaluate_serial_he3_field_map( &
      mesh, state, angular, energy, 0.5_rk, 0.15_rk, 2, 0.1_rk, &
      mapped, diagnostics, succeeded)
    call require(succeeded, "serial field map reported a failed point")
    call require(diagnostics%requested_point_count == 9 .and. &
                 diagnostics%completed_point_count == 9 .and. &
                 diagnostics%failed_point_count == 0, &
                 "serial field-map point counts are wrong")
    call require(diagnostics%attempted_contributions == 36 .and. &
                 diagnostics%accepted_contributions == 36 .and. &
                 diagnostics%failed_propagators == 0, &
                 "serial field-map contribution counts are wrong")

    do point = 1, mesh%point_count()
      call evaluate_serial_he3_point_map( &
        mesh, state, mesh%point_coordinate(point), angular, energy, &
        0.5_rk, 0.15_rk, 2, 0.1_rk, expected_gap, expected_current, &
        point_diagnostics, point_succeeded)
      call require(point_succeeded, "independent point map failed")
      call require(maxval(abs(mapped%order_parameter(:, :, point) - &
                                  expected_gap)) < 2.0e-14_rk, &
                   "field map differs from independent point evaluation")
      call require(maxval(abs(mapped%current_mean_field(:, point) - &
                                  expected_current)) < 2.0e-16_rk, &
                   "field-map mean field differs from point evaluation")
    end do
  end subroutine test_complete_field_map


  subroutine test_residual_and_anderson_adapter()
    type(anderson_report_t) :: anderson_report
    type(cartesian_mesh_2d_t) :: mesh
    type(he3_field_residual_report_t) :: residual_report
    type(legacy_anderson_t) :: accelerator
    type(spinful_state_2d_t) :: current, expected, mapped, next
    logical, allocatable :: active_point(:)
    integer :: point

    call make_uniform_cartesian_mesh( &
      -1.0_rk, 1.0_rk, 2, -1.0_rk, 1.0_rk, 1, mesh)
    call allocate_spinful_state_2d(mesh, current)
    do point = 1, mesh%point_count()
      current%order_parameter(1, 1, point) = &
        cmplx(0.2_rk + 0.01_rk * real(point, rk), 0.03_rk, rk)
      current%current_mean_field(2, point) = -0.02_rk * real(point, rk)
    end do
    mapped = current
    mapped%order_parameter(2, 3, 4) = cmplx(0.4_rk, -0.2_rk, rk)
    mapped%current_mean_field(1, 4) = 0.3_rk
    mapped%order_parameter(1, 1, 2) = &
      mapped%order_parameter(1, 1, 2) + cmplx(0.01_rk, -0.02_rk, rk)

    call compute_he3_field_residual(current, mapped, residual_report)
    call require(residual_report%value_count == 21 * mesh%point_count(), &
                 "field residual used the wrong vector length")
    call require(residual_report%maximum_point == 4, &
                 "field residual identified the wrong maximum point")
    call require(residual_report%maximum_absolute_residual > 0.39_rk .and. &
                 residual_report%rms_residual > 0.0_rk .and. &
                 residual_report%mean_absolute_residual > 0.0_rk, &
                 "field residual statistics are invalid")

    allocate(active_point(mesh%point_count()), source=.false.)
    active_point(2) = .true.
    call compute_he3_field_residual( &
      current, mapped, residual_report, active_point_mask=active_point)
    call require(residual_report%value_count == 21, &
                 "masked residual used the wrong active vector length")
    call require(residual_report%maximum_point == 2, &
                 "masked residual included a frozen halo point")
    call require(abs(residual_report%maximum_absolute_residual - 0.02_rk) < &
                 8.0_rk * epsilon(1.0_rk), &
                 "masked residual maximum is wrong")

    call accelerator%initialize( &
      new_src_iteration_vector_size_2d(current), 5, 0.1_rk, 5.0_rk)
    call update_he3_state_with_anderson( &
      mesh, accelerator, current, mapped, 0.0_rk, next, &
      anderson_report, residual_report)
    expected = current
    expected%order_parameter = current%order_parameter + &
      cmplx(0.01_rk, 0.0_rk, rk) * &
      (mapped%order_parameter - current%order_parameter)
    expected%current_mean_field = current%current_mean_field + &
      0.01_rk * (mapped%current_mean_field - current%current_mean_field)
    call require(maxval(abs(next%order_parameter - &
                                expected%order_parameter)) < 2.0e-15_rk, &
                 "first field Anderson update differs from simple mixing")
    call require(maxval(abs(next%current_mean_field - &
                                expected%current_mean_field)) < 2.0e-15_rk, &
                 "first mean-field Anderson update is wrong")
    call require(anderson_report%iteration == 1 .and. &
                 anderson_report%history_size == 1, &
                 "field Anderson adapter returned the wrong state")
  end subroutine test_residual_and_anderson_adapter


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_he3_field_map_2d
