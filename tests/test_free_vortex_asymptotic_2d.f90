program test_free_vortex_asymptotic_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, &
                                make_uniform_cartesian_mesh
  use spinful_state_2d, only : spinful_state_2d_t, allocate_spinful_state_2d
  use radial_mesh, only : radial_mesh_t, make_uniform_radial_mesh
  use axial_radial_embedding, only : axial_radial_profile_t, &
    allocate_axial_radial_profile, sample_axial_radial_profile
  use straight_trajectory_2d, only : straight_trajectory_2d_t, &
                                     build_straight_trajectory_2d
  use trajectory_self_energy_2d, only : trajectory_self_energy_2d_t, &
                                        sample_trajectory_self_energy_2d
  use free_vortex_asymptotic_2d, only : free_vortex_endpoint_2d_t, &
    sample_free_vortex_asymptotic_state_2d, &
    extend_with_free_vortex_asymptotic_2d, &
    set_axisymmetric_radial_reference
  implicit none

  call test_source_decay_powers()
  call test_axisymmetric_radial_reference()
  call test_trajectory_extension()

  print '(a)', "free-vortex asymptotic endpoint tests passed"

contains

  subroutine test_source_decay_powers()
    type(cartesian_mesh_2d_t) :: mesh
    type(free_vortex_endpoint_2d_t) :: endpoint
    type(spinful_state_2d_t) :: state
    complex(rk) :: sampled(3, 3)
    real(rk) :: current(3)

    call make_test_state(mesh, state)
    endpoint = make_endpoint()
    call require(endpoint%is_valid_for(mesh), &
      "free-vortex endpoint was unexpectedly invalid")

    call sample_free_vortex_asymptotic_state_2d( &
      mesh, state, [4.0_rk, 0.0_rk], endpoint, sampled, current)
    call require(abs(sampled(1, 1) - cmplx(0.32_rk, 0.0_rk, rk)) < &
                   2.0e-15_rk, &
      "in-plane diagonal correction does not decay as inverse radius squared")
    call require(abs(sampled(1, 2) - cmplx(0.0075_rk, 0.0_rk, rk)) < &
                   2.0e-15_rk, &
      "in-plane off-diagonal correction has the wrong decay")
    call require(abs(sampled(1, 3) - cmplx(0.06_rk, 0.0_rk, rk)) < &
                   2.0e-15_rk, &
      "mixed z/in-plane correction does not decay as inverse radius")
    call require(maxval(abs(current - [0.02_rk, -0.01_rk, 0.005_rk])) < &
                   2.0e-15_rk, &
      "current-related asymptotic field has the wrong inverse-radius tail")

    call sample_free_vortex_asymptotic_state_2d( &
      mesh, state, [2.0_rk, 0.0_rk], endpoint, sampled, current)
    call require(abs(sampled(1, 1) - cmplx(0.38_rk, 0.0_rk, rk)) < &
                   2.0e-15_rk .and. &
                 abs(sampled(1, 3) - cmplx(0.12_rk, 0.0_rk, rk)) < &
                   2.0e-15_rk, &
      "asymptotic continuation is not continuous at the field boundary")
  end subroutine test_source_decay_powers


  subroutine test_axisymmetric_radial_reference()
    type(axial_radial_profile_t) :: profile
    type(cartesian_mesh_2d_t) :: mesh
    type(free_vortex_endpoint_2d_t) :: endpoint
    type(radial_mesh_t) :: radial_grid
    type(spinful_state_2d_t) :: state
    complex(rk) :: direct_gap(3, 3), endpoint_gap(3, 3)
    real(rk) :: direct_current(3), endpoint_current(3), radius
    integer :: point
    logical :: inside

    call make_test_state(mesh, state)
    call make_uniform_radial_mesh(4.0_rk, 8, radial_grid)
    call allocate_axial_radial_profile(radial_grid, profile)
    do point = 0, radial_grid%point_count() - 1
      radius = radial_grid%r(point)
      profile%harmonic(2, 2, point) = &
        cmplx(0.1_rk + 0.02_rk * radius, 0.01_rk * radius, rk)
      profile%azimuthal_mean_field(point) = 0.03_rk * radius
    end do
    endpoint = make_endpoint()
    call set_axisymmetric_radial_reference( &
      endpoint, radial_grid, profile, 3.0_rk)
    call require(endpoint%is_valid_for(mesh), &
      "axisymmetric radial endpoint was unexpectedly invalid")

    call sample_axial_radial_profile( &
      radial_grid, profile, 2.5_rk, 0.0_rk, endpoint%winding, &
      direct_gap, direct_current, inside)
    call require(inside, "direct radial endpoint reference sample failed")
    call sample_free_vortex_asymptotic_state_2d( &
      mesh, state, [2.5_rk, 0.0_rk], endpoint, endpoint_gap, endpoint_current)
    call require(maxval(abs(endpoint_gap - direct_gap)) < 2.0e-15_rk .and. &
                 maxval(abs(endpoint_current - direct_current)) < 2.0e-15_rk, &
      "radial-reference endpoint changed a profile sample outside the box")
  end subroutine test_axisymmetric_radial_reference


  subroutine test_trajectory_extension()
    type(cartesian_mesh_2d_t) :: mesh
    type(free_vortex_endpoint_2d_t) :: endpoint
    type(spinful_state_2d_t) :: state
    type(straight_trajectory_2d_t) :: trajectory
    type(trajectory_self_energy_2d_t) :: extended, interior
    real(rk), allocatable :: coordinate(:)
    integer :: entry_sample, exit_sample, sample, target

    call make_test_state(mesh, state)
    endpoint = make_endpoint()
    call build_straight_trajectory_2d( &
      mesh, [0.0_rk, 0.0_rk], [1.0_rk, 0.0_rk, 0.0_rk], &
      0.75_rk, trajectory)
    call sample_trajectory_self_energy_2d( &
      state, trajectory, 0.7_rk, interior)
    call extend_with_free_vortex_asymptotic_2d( &
      mesh, state, trajectory, interior, 0.7_rk, endpoint, coordinate, &
      target, extended)

    call require(size(coordinate) > trajectory%sample_count(), &
      "free-vortex path was not extended")
    call require(abs(coordinate(1) + endpoint%outer_radius) < 2.0e-14_rk .and. &
                 abs(coordinate(size(coordinate)) - endpoint%outer_radius) < &
                   2.0e-14_rk, &
      "free-vortex path does not reach the requested outer circle")
    call require(abs(coordinate(target)) < 2.0e-15_rk, &
      "free-vortex path lost its target point")
    call require(maxval(coordinate(2:) - &
                        coordinate(:size(coordinate) - 1)) <= &
                   0.75_rk * (1.0_rk + 32.0_rk * epsilon(1.0_rk)), &
      "free-vortex path contains an interval larger than its step bound")

    entry_sample = 0
    exit_sample = 0
    do sample = 1, size(coordinate)
      if (abs(coordinate(sample) - trajectory%entry_coordinate) < 1.0e-14_rk) &
        entry_sample = sample
      if (abs(coordinate(sample) - trajectory%exit_coordinate) < 1.0e-14_rk) &
        exit_sample = sample
    end do
    call require(entry_sample > 0 .and. exit_sample > 0, &
      "free-vortex path lost a field-boundary sample")
    call require(maxval(abs(extended%triplet_pair_potential(:, entry_sample) - &
                                interior%triplet_pair_potential(:, 1))) < &
                   2.0e-15_rk .and. &
                 maxval(abs(extended%triplet_pair_potential(:, exit_sample) - &
                                interior%triplet_pair_potential(:, &
                                  interior%sample_count()))) < 2.0e-15_rk, &
      "free-vortex extension changed an interior boundary value")
  end subroutine test_trajectory_extension


  subroutine make_test_state(mesh, state)
    type(cartesian_mesh_2d_t), intent(out) :: mesh
    type(spinful_state_2d_t), intent(out) :: state

    complex(rk) :: bulk_value
    real(rk) :: angle, radius, x, y
    integer :: component, point, x_node, y_node

    call make_uniform_cartesian_mesh( &
      -2.0_rk, 2.0_rk, 4, -2.0_rk, 2.0_rk, 4, mesh)
    call allocate_spinful_state_2d(mesh, state)
    do y_node = 0, mesh%y_cell_count()
      y = mesh%y_coordinate(y_node)
      do x_node = 0, mesh%x_cell_count()
        x = mesh%x_coordinate(x_node)
        point = mesh%point_index(x_node, y_node)
        radius = sqrt(x * x + y * y)
        if (radius > 0.0_rk) then
          angle = atan2(y, x)
        else
          angle = 0.0_rk
        end if
        bulk_value = cmplx(0.3_rk * cos(angle), &
                           0.3_rk * sin(angle), kind=rk)
        do component = 1, 3
          state%order_parameter(component, component, point) = bulk_value
        end do
        state%order_parameter(1, 1, point) = &
          state%order_parameter(1, 1, point) + cmplx(0.08_rk, 0.0_rk, rk)
        state%order_parameter(1, 2, point) = cmplx(0.03_rk, 0.0_rk, rk)
        state%order_parameter(1, 3, point) = cmplx(0.12_rk, 0.0_rk, rk)
        state%current_mean_field(:, point) = [0.04_rk, -0.02_rk, 0.01_rk]
      end do
    end do
  end subroutine make_test_state


  pure function make_endpoint() result(endpoint)
    type(free_vortex_endpoint_2d_t) :: endpoint

    endpoint%enabled = .true.
    endpoint%center = 0.0_rk
    endpoint%bulk_gap = 0.3_rk
    endpoint%winding = 1.0_rk
    endpoint%outer_radius = 5.0_rk
    endpoint%maximum_step = 0.4_rk
  end function make_endpoint


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_free_vortex_asymptotic_2d
