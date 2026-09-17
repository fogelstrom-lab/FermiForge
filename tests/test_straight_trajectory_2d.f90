program test_straight_trajectory_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, make_uniform_cartesian_mesh
  use spinful_state_2d, only : spinful_state_2d_t, allocate_spinful_state_2d
  use straight_trajectory_2d, only : straight_trajectory_2d_t, &
                                     build_straight_trajectory_2d, &
                                     sample_pair_potential_along_trajectory_2d
  implicit none

  call test_oblique_domain_intersection()
  call test_nonzero_pz_projection()
  call test_invariant_axis_trajectory()
  call test_cached_affine_field_sampling()

  print '(a)', "straight 2D trajectory tests passed"

contains

  subroutine test_oblique_domain_intersection()
    type(cartesian_mesh_2d_t) :: mesh
    type(straight_trajectory_2d_t) :: trajectory
    real(rk), parameter :: origin(2) = [0.5_rk, 0.5_rk]
    real(rk), parameter :: momentum(3) = [0.6_rk, 0.8_rk, 0.0_rk]
    real(rk) :: entry(2), exit(2)

    call make_uniform_cartesian_mesh(-2.0_rk, 3.0_rk, 20, &
                                     -1.0_rk, 2.0_rk, 12, mesh)
    call build_straight_trajectory_2d(mesh, origin, momentum, 0.2_rk, trajectory)

    call require(trajectory%is_valid(), "oblique trajectory is invalid")
    call require(abs(trajectory%entry_coordinate + 1.875_rk) < 1.0e-13_rk, &
                 "oblique trajectory entry coordinate is wrong")
    call require(abs(trajectory%exit_coordinate - 1.875_rk) < 1.0e-13_rk, &
                 "oblique trajectory exit coordinate is wrong")
    entry = origin + trajectory%entry_coordinate * momentum(1:2)
    exit = origin + trajectory%exit_coordinate * momentum(1:2)
    call require(abs(entry(2) + 1.0_rk) < 1.0e-13_rk, &
                 "oblique trajectory does not enter on the lower boundary")
    call require(abs(exit(2) - 2.0_rk) < 1.0e-13_rk, &
                 "oblique trajectory does not exit on the upper boundary")
    call require(maxval(trajectory%path_coordinate(2:) - &
                        trajectory%path_coordinate(:trajectory%sample_count()-1)) &
                 <= 0.2_rk * (1.0_rk + 64.0_rk * epsilon(1.0_rk)), &
                 "trajectory contains a step larger than requested")
    call require(abs(trajectory%path_coordinate(trajectory%target_sample)) < &
                 tiny(1.0_rk), &
                 "trajectory does not contain its target exactly")
  end subroutine test_oblique_domain_intersection


  subroutine test_nonzero_pz_projection()
    type(cartesian_mesh_2d_t) :: mesh
    type(straight_trajectory_2d_t) :: trajectory
    real(rk), parameter :: origin(2) = [0.5_rk, 0.5_rk]
    real(rk), parameter :: momentum(3) = [0.6_rk, 0.0_rk, 0.8_rk]

    call make_uniform_cartesian_mesh(-2.0_rk, 3.0_rk, 20, &
                                     -1.0_rk, 2.0_rk, 12, mesh)
    call build_straight_trajectory_2d(mesh, origin, momentum, 0.25_rk, trajectory)

    call require(abs(trajectory%entry_coordinate + 2.5_rk / 0.6_rk) < 1.0e-13_rk, &
                 "p_z trajectory entry did not use projected x speed")
    call require(abs(trajectory%exit_coordinate - 2.5_rk / 0.6_rk) < 1.0e-13_rk, &
                 "p_z trajectory exit did not use projected x speed")
  end subroutine test_nonzero_pz_projection


  subroutine test_invariant_axis_trajectory()
    type(cartesian_mesh_2d_t) :: mesh
    type(straight_trajectory_2d_t) :: trajectory

    call make_uniform_cartesian_mesh(-1.0_rk, 1.0_rk, 4, &
                                     -1.0_rk, 1.0_rk, 4, mesh)
    call build_straight_trajectory_2d(mesh, [0.2_rk, -0.3_rk], &
                                      [0.0_rk, 0.0_rk, 1.0_rk], &
                                      0.1_rk, trajectory)

    call require(trajectory%parallel_to_invariant_axis, &
                 "z-parallel trajectory was not identified")
    call require(trajectory%sample_count() == 1, &
                 "z-parallel trajectory should have one spatial sample")
    call require(trajectory%target_sample == 1, &
                 "z-parallel trajectory target index is wrong")
  end subroutine test_invariant_axis_trajectory


  subroutine test_cached_affine_field_sampling()
    type(cartesian_mesh_2d_t) :: mesh
    type(spinful_state_2d_t) :: state
    type(straight_trajectory_2d_t) :: trajectory
    complex(rk), allocatable :: pair_potential(:, :)
    complex(rk) :: expected(3)
    real(rk), allocatable :: sampled_mean(:, :)
    real(rk), parameter :: origin(2) = [0.1_rk, -0.2_rk]
    real(rk), parameter :: momentum(3) = [0.36_rk, 0.48_rk, 0.8_rk]
    real(rk) :: x, y
    integer :: orbital, point, sample, spin, x_node, y_node

    call make_uniform_cartesian_mesh(-1.0_rk, 1.0_rk, 8, &
                                     -1.0_rk, 1.0_rk, 8, mesh)
    call allocate_spinful_state_2d(mesh, state)
    do y_node = 0, mesh%y_cell_count()
      y = mesh%y_coordinate(y_node)
      do x_node = 0, mesh%x_cell_count()
        x = mesh%x_coordinate(x_node)
        point = mesh%point_index(x_node, y_node)
        do spin = 1, 3
          do orbital = 1, 3
            state%order_parameter(spin, orbital, point) = cmplx( &
              real(spin + orbital, rk) + 0.2_rk * x - 0.3_rk * y, &
              real(spin - orbital, rk) - 0.1_rk * x + 0.4_rk * y, kind=rk)
          end do
        end do
        state%current_mean_field(:, point) = [x, y, x - y]
      end do
    end do

    call build_straight_trajectory_2d(mesh, origin, momentum, 0.15_rk, trajectory)
    call sample_pair_potential_along_trajectory_2d(state, trajectory, &
                                                    pair_potential, sampled_mean)
    do sample = 1, trajectory%sample_count()
      x = origin(1) + trajectory%path_coordinate(sample) * momentum(1)
      y = origin(2) + trajectory%path_coordinate(sample) * momentum(2)
      expected = cmplx(0.0_rk, 0.0_rk, kind=rk)
      do spin = 1, 3
        do orbital = 1, 3
          expected(spin) = expected(spin) + &
            cmplx(momentum(orbital), 0.0_rk, kind=rk) * cmplx( &
              real(spin + orbital, rk) + 0.2_rk * x - 0.3_rk * y, &
              real(spin - orbital, rk) - 0.1_rk * x + 0.4_rk * y, kind=rk)
        end do
      end do
      call require(maxval(abs(pair_potential(:, sample) - expected)) < 3.0e-13_rk, &
                   "cached trajectory does not reproduce affine pair potential")
      call require(maxval(abs(sampled_mean(:, sample) - [x, y, x - y])) < &
                   3.0e-13_rk, &
                   "cached trajectory does not reproduce affine mean field")
    end do
  end subroutine test_cached_affine_field_sampling


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_straight_trajectory_2d
