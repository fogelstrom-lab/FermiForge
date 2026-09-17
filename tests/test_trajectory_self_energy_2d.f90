program test_trajectory_self_energy_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, make_uniform_cartesian_mesh
  use spinful_state_2d, only : spinful_state_2d_t, allocate_spinful_state_2d
  use straight_trajectory_2d, only : straight_trajectory_2d_t, &
                                     build_straight_trajectory_2d
  use trajectory_self_energy_2d, only : trajectory_self_energy_2d_t, &
                                        legacy_forward_branch, &
                                        legacy_reverse_branch, &
                                        sample_trajectory_self_energy_2d, &
                                        make_legacy_riccati_coefficients
  implicit none

  call test_affine_self_energy_sampling()
  call test_legacy_branch_conventions()

  print '(a)', "2D trajectory self-energy tests passed"

contains

  subroutine test_affine_self_energy_sampling()
    type(cartesian_mesh_2d_t) :: mesh
    type(spinful_state_2d_t) :: state
    type(straight_trajectory_2d_t) :: trajectory
    type(trajectory_self_energy_2d_t) :: self_energy
    complex(rk) :: expected(3)
    real(rk), parameter :: momentum(3) = [0.36_rk, 0.48_rk, 0.8_rk]
    real(rk), parameter :: origin(2) = [0.1_rk, -0.2_rk]
    real(rk) :: expected_shift, x, y
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

    call build_straight_trajectory_2d(mesh, origin, momentum, 0.17_rk, trajectory)
    call sample_trajectory_self_energy_2d(state, trajectory, 1.7_rk, self_energy)
    call require(self_energy%is_valid(), "sampled trajectory self energy is invalid")
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
      expected_shift = 1.7_rk * dot_product(momentum, [x, y, x - y])
      call require(maxval(abs( &
        self_energy%triplet_pair_potential(:, sample) - expected)) < 3.0e-13_rk, &
        "trajectory pair self energy does not reproduce an affine field")
      call require(abs(self_energy%diagonal_shift(sample) - expected_shift) < &
                   3.0e-13_rk, &
                   "trajectory diagonal self energy uses the wrong projection")
    end do
  end subroutine test_affine_self_energy_sampling


  subroutine test_legacy_branch_conventions()
    type(trajectory_self_energy_2d_t) :: self_energy
    complex(rk) :: pair_a(0:3), pair_b(0:3), energy, exchange(3)
    complex(rk), parameter :: spectral = cmplx(0.0_rk, 0.7_rk, kind=rk)

    allocate(self_energy%triplet_pair_potential(3, 1))
    allocate(self_energy%diagonal_shift(1))
    self_energy%triplet_pair_potential(:, 1) = &
      [cmplx(0.1_rk, 0.2_rk, kind=rk), &
       cmplx(-0.3_rk, 0.4_rk, kind=rk), &
       cmplx(0.5_rk, -0.6_rk, kind=rk)]
    self_energy%diagonal_shift(1) = 0.08_rk

    call make_legacy_riccati_coefficients( &
      self_energy, 1, legacy_forward_branch, spectral, pair_a, pair_b, &
      energy, exchange)
    call require(maxval(abs(pair_a(1:3) - &
      self_energy%triplet_pair_potential(:, 1))) < tiny(1.0_rk), &
      "forward Riccati branch changed the triplet pair potential")
    call require(maxval(abs(pair_b - conjg(pair_a))) < tiny(1.0_rk), &
                 "forward Riccati branch has the wrong conjugate field")

    call make_legacy_riccati_coefficients( &
      self_energy, 1, legacy_reverse_branch, spectral, pair_a, pair_b, &
      energy, exchange)
    call require(maxval(abs(pair_b(1:3) - &
      self_energy%triplet_pair_potential(:, 1))) < tiny(1.0_rk), &
      "reverse Riccati branch changed the triplet pair potential")
    call require(maxval(abs(pair_a - conjg(pair_b))) < tiny(1.0_rk), &
                 "reverse Riccati branch has the wrong conjugate field")
    call require(abs(energy - cmplx(-0.08_rk, 0.7_rk, kind=rk)) < &
                 tiny(1.0_rk), &
                 "legacy Riccati effective energy has the wrong shift")
    call require(maxval(abs(exchange)) < tiny(1.0_rk), &
                 "inactive spin exchange field is not explicit zero")
  end subroutine test_legacy_branch_conventions


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_trajectory_self_energy_2d
