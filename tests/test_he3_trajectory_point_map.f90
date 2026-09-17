program test_he3_trajectory_point_map
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, make_uniform_cartesian_mesh
  use spinful_state_2d, only : spinful_state_2d_t, allocate_spinful_state_2d
  use straight_trajectory_2d, only : straight_trajectory_2d_t, &
                                     build_straight_trajectory_2d
  use quasiclassical_propagator, only : quasiclassical_propagator_t, &
                                        reconstruct_legacy_quasiclassical_propagator
  use he3_self_consistency_integrand, only : he3_point_map_accumulator_t, &
                                             accumulate_he3_point_map_contribution
  use he3_trajectory_point_map, only : accumulate_legacy_trajectory_point_map
  implicit none

  call test_complete_trajectory_contribution()

  print '(a)', "3He trajectory point-map tests passed"

contains

  subroutine test_complete_trajectory_contribution()
    type(cartesian_mesh_2d_t) :: mesh
    type(spinful_state_2d_t) :: state
    type(straight_trajectory_2d_t) :: trajectory
    type(quasiclassical_propagator_t) :: frozen_propagator
    type(he3_point_map_accumulator_t) :: expected, mapped
    complex(rk) :: forward(0:3), reverse(0:3)
    real(rk) :: momentum(3), normalization_error, planar_magnitude
    real(rk) :: source_coordinate, x, y
    integer :: orbital, point, spin, x_node, y_node
    logical :: frozen_succeeded, mapped_succeeded

    planar_magnitude = sqrt(1.0_rk - 0.35_rk**2)
    momentum = [cos(0.41_rk) * planar_magnitude, &
                sin(0.41_rk) * planar_magnitude, 0.35_rk]
    call make_uniform_cartesian_mesh( &
      -0.4_rk * momentum(1), 0.4_rk * momentum(1), 4, &
      -0.4_rk * momentum(2), 0.4_rk * momentum(2), 4, mesh)
    call allocate_spinful_state_2d(mesh, state)

    do y_node = 0, mesh%y_cell_count()
      y = mesh%y_coordinate(y_node)
      do x_node = 0, mesh%x_cell_count()
        x = mesh%x_coordinate(x_node)
        point = mesh%point_index(x_node, y_node)
        source_coordinate = &
          (x * momentum(1) + y * momentum(2)) / &
          (momentum(1)**2 + momentum(2)**2)
        state%order_parameter(1, 1, point) = cmplx( &
          0.23_rk + 0.075_rk * source_coordinate, &
          0.04_rk - 0.06_rk * source_coordinate, kind=rk) / &
          cmplx(momentum(1), 0.0_rk, kind=rk)
        state%order_parameter(2, 1, point) = cmplx( &
          -0.08_rk + 0.045_rk * source_coordinate, &
          0.17_rk + 0.03_rk * source_coordinate, kind=rk) / &
          cmplx(momentum(1), 0.0_rk, kind=rk)
        state%order_parameter(3, 1, point) = cmplx( &
          0.11_rk - 0.035_rk * source_coordinate, &
          -0.13_rk + 0.05_rk * source_coordinate, kind=rk) / &
          cmplx(momentum(1), 0.0_rk, kind=rk)
        do spin = 1, 3
          do orbital = 2, 3
            state%order_parameter(spin, orbital, point) = &
              cmplx(0.0_rk, 0.0_rk, kind=rk)
          end do
        end do
        state%current_mean_field(:, point) = &
          0.09_rk * source_coordinate * momentum
      end do
    end do

    call build_straight_trajectory_2d( &
      mesh, [0.0_rk, 0.0_rk], momentum, 0.2_rk, trajectory)
    call require(trajectory%sample_count() == 5 .and. &
                 trajectory%target_sample == 3, &
                 "point-map regression trajectory has the wrong samples")

    call frozen_coherences(forward, reverse)
    call reconstruct_legacy_quasiclassical_propagator( &
      forward, reverse, frozen_propagator, frozen_succeeded)
    call require(frozen_succeeded, "frozen point-map propagator is singular")
    call accumulate_he3_point_map_contribution( &
      expected, frozen_propagator, momentum, 0.11_rk * 0.17_rk, &
      0.42_rk, 0.37_rk)

    call accumulate_legacy_trajectory_point_map( &
      state, trajectory, 1.0_rk, cmplx(0.0_rk, 0.83_rk, kind=rk), &
      3, 0.0_rk, 0.11_rk * 0.17_rk, 0.42_rk, 0.37_rk, mapped, &
      normalization_error, mapped_succeeded)
    call require(mapped_succeeded, "trajectory point-map reconstruction failed")
    call require(normalization_error < 8.0e-15_rk, &
                 "trajectory point-map propagator is not normalized")
    if (maxval(abs(mapped%gap - expected%gap)) >= 2.0e-15_rk) then
      print '(a,es24.16)', "maximum gap difference: ", &
        maxval(abs(mapped%gap - expected%gap))
      print '(a,3(2(1x,es18.10)))', "mapped row 1:  ", mapped%gap(1, :)
      print '(a,3(2(1x,es18.10)))', "expected row 1:", expected%gap(1, :)
    end if
    call require(maxval(abs(mapped%gap - expected%gap)) < 2.0e-15_rk, &
                 "complete trajectory gap contribution differs from source")
    call require(maxval(abs(mapped%current_mean_field - &
                 expected%current_mean_field)) < 2.0e-17_rk, &
                 "complete trajectory current contribution differs from source")
    call require(mapped%contribution_count == 1, &
                 "complete trajectory contribution count is wrong")
  end subroutine test_complete_trajectory_contribution


  subroutine frozen_coherences(forward, reverse)
    complex(rk), intent(out) :: forward(0:3), reverse(0:3)

    forward = cmplx(0.0_rk, 0.0_rk, kind=rk)
    forward(1:3) = [ &
      cmplx(-3.45416635611957e-2_rk, 1.17489705712979e-1_rk, kind=rk), &
      cmplx(-9.01737612952456e-2_rk, -5.58797007439699e-2_rk, kind=rk), &
      cmplx(8.02906423983453e-2_rk, 7.10955985562741e-2_rk, kind=rk)]
    reverse = cmplx(0.0_rk, 0.0_rk, kind=rk)
    reverse(1:3) = [ &
      cmplx(1.62659859791399e-2_rk, 1.42678313138293e-1_rk, kind=rk), &
      cmplx(9.95709656830779e-2_rk, -4.16922195641589e-2_rk, kind=rk), &
      cmplx(-6.37043298757304e-2_rk, 5.90495834745839e-2_rk, kind=rk)]
  end subroutine frozen_coherences


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_he3_trajectory_point_map
