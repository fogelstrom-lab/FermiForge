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
    set_axisymmetric_radial_reference, set_state_asymptotic_fit, &
    set_state_asymptotic_elliptical_fit, &
    fit_state_bulk_reference, apply_free_vortex_asymptotic_halo_2d
  implicit none

  call test_source_decay_powers()
  call test_axisymmetric_radial_reference()
  call test_angle_dependent_two_term_fit()
  call test_eq20_block_powers_and_mean_field()
  call test_elliptical_fit_surfaces()
  call test_fitted_rotated_bulk_reference()
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


  subroutine test_angle_dependent_two_term_fit()
    type(cartesian_mesh_2d_t) :: mesh
    type(free_vortex_endpoint_2d_t) :: endpoint
    type(spinful_state_2d_t) :: completed, state
    complex(rk) :: expected_gap(3, 3), sampled_gap(3, 3)
    real(rk) :: expected_current(3), sampled_current(3), position(2)
    real(rk) :: radius
    logical, allocatable :: active(:)
    integer :: point, x_node, y_node

    call make_uniform_cartesian_mesh( &
      -3.0_rk, 3.0_rk, 12, -3.0_rk, 3.0_rk, 12, mesh)
    endpoint = make_endpoint()
    call set_state_asymptotic_fit(endpoint, 1.0_rk, 2.0_rk)
    call allocate_spinful_state_2d(mesh, state)
    do point = 1, mesh%point_count()
      position = mesh%point_coordinate(point)
      call evaluate_known_two_term_field( &
        position, endpoint, state%order_parameter(:, :, point), &
        state%current_mean_field(:, point))
    end do
    call require(endpoint%is_valid_for(mesh), &
      "two-term state endpoint was unexpectedly invalid")

    position = [4.0_rk, 0.0_rk]
    call sample_free_vortex_asymptotic_state_2d( &
      mesh, state, position, endpoint, sampled_gap, sampled_current)
    call evaluate_known_two_term_field( &
      position, endpoint, expected_gap, expected_current)
    call require(maxval(abs(sampled_gap - expected_gap)) < 3.0e-14_rk .and. &
                 maxval(abs(sampled_current - expected_current)) < &
                   3.0e-14_rk, &
      "two-term fit failed on the positive x axis")

    position = [-4.0_rk, 0.0_rk]
    call sample_free_vortex_asymptotic_state_2d( &
      mesh, state, position, endpoint, sampled_gap, sampled_current)
    call evaluate_known_two_term_field( &
      position, endpoint, expected_gap, expected_current)
    call require(maxval(abs(sampled_gap - expected_gap)) < 3.0e-14_rk .and. &
                 maxval(abs(sampled_current - expected_current)) < &
                   3.0e-14_rk, &
      "two-term fit lost the polar-angle dependence")

    allocate(active(mesh%point_count()), source=.false.)
    do point = 1, mesh%point_count()
      position = mesh%point_coordinate(point)
      radius = sqrt(dot_product(position, position))
      active(point) = radius <= 2.75_rk
    end do
    x_node = mesh%x_cell_count()
    y_node = mesh%y_cell_count() / 2
    point = mesh%point_index(x_node, y_node)
    state%order_parameter(:, :, point) = cmplx(99.0_rk, -17.0_rk, rk)
    state%current_mean_field(:, point) = 23.0_rk
    call apply_free_vortex_asymptotic_halo_2d( &
      mesh, state, active, endpoint, completed)
    position = mesh%point_coordinate(point)
    call evaluate_known_two_term_field( &
      position, endpoint, expected_gap, expected_current)
    call require(maxval(abs(completed%order_parameter(:, :, point) - &
                            expected_gap)) < 3.0e-14_rk .and. &
                 maxval(abs(completed%current_mean_field(:, point) - &
                            expected_current)) < 3.0e-14_rk, &
      "dependent asymptotic halo was not regenerated from the fit")
  end subroutine test_angle_dependent_two_term_fit


  subroutine test_eq20_block_powers_and_mean_field()
    type(cartesian_mesh_2d_t) :: mesh
    type(free_vortex_endpoint_2d_t) :: endpoint
    type(spinful_state_2d_t) :: state
    complex(rk) :: bulk_value, expected_diagonal, expected_mixed
    complex(rk) :: sampled_gap(3, 3)
    real(rk) :: angle, current_one(3), current_two(3)
    real(rk) :: expected_current(3), position(2), radius, ratio
    real(rk) :: sampled_current(3)
    integer :: point

    call make_uniform_cartesian_mesh( &
      -3.0_rk, 3.0_rk, 12, -3.0_rk, 3.0_rk, 12, mesh)
    endpoint = make_endpoint()
    call set_state_asymptotic_fit(endpoint, 1.0_rk, 2.0_rk)
    call allocate_spinful_state_2d(mesh, state)
    current_one = [0.07_rk, -0.04_rk, 0.03_rk]
    current_two = [0.02_rk, 0.05_rk, -0.01_rk]
    do point = 1, mesh%point_count()
      position = mesh%point_coordinate(point)
      radius = sqrt(dot_product(position, position))
      if (radius > 64.0_rk * epsilon(1.0_rk)) then
        angle = atan2(position(2), position(1))
        ratio = endpoint%matching_radius / radius
      else
        angle = 0.0_rk
        ratio = 0.0_rk
      end if
      bulk_value = cmplx(endpoint%bulk_gap, 0.0_rk, rk) * &
        exp(cmplx(0.0_rk, endpoint%winding * angle, rk))
      state%order_parameter(:, :, point) = cmplx(0.0_rk, 0.0_rk, rk)
      state%order_parameter(1, 1, point) = bulk_value + &
        ratio * cmplx(0.04_rk, 0.01_rk, rk) + &
        ratio**2 * cmplx(0.06_rk, -0.02_rk, rk)
      state%order_parameter(2, 2, point) = bulk_value
      state%order_parameter(3, 3, point) = bulk_value
      state%order_parameter(1, 3, point) = &
        ratio * cmplx(-0.03_rk, 0.02_rk, rk) + &
        ratio**2 * cmplx(0.05_rk, 0.01_rk, rk)
      state%current_mean_field(:, point) = &
        ratio * current_one + ratio**2 * current_two
    end do

    call sample_free_vortex_asymptotic_state_2d( &
      mesh, state, [4.0_rk, 0.0_rk], endpoint, sampled_gap, sampled_current)
    ratio = endpoint%matching_radius / 4.0_rk
    ! At Rc the deliberately contaminated diagonal and mixed departures are
    ! respectively 0.10-0.01i and 0.02+0.03i.  Equation (20) assigns all of
    ! the former to 1/r^2 and all of the latter to 1/r.
    expected_diagonal = cmplx(endpoint%bulk_gap, 0.0_rk, rk) + &
      ratio**2 * cmplx(0.10_rk, -0.01_rk, rk)
    expected_mixed = ratio * cmplx(0.02_rk, 0.03_rk, rk)
    expected_current = ratio * current_one + ratio**2 * current_two
    call require(abs(sampled_gap(1, 1) - expected_diagonal) < 3.0e-14_rk, &
      "Eq. (20) allowed a forbidden 1/r diagonal order-parameter tail")
    call require(abs(sampled_gap(1, 3) - expected_mixed) < 3.0e-14_rk, &
      "Eq. (20) allowed a forbidden 1/r^2 mixed order-parameter tail")
    call require(maxval(abs(sampled_current - expected_current)) < &
                   3.0e-14_rk, &
      "Fermi-liquid continuation lost one of its two radial powers")
  end subroutine test_eq20_block_powers_and_mean_field


  subroutine test_elliptical_fit_surfaces()
    type(cartesian_mesh_2d_t) :: mesh
    type(free_vortex_endpoint_2d_t) :: endpoint
    type(spinful_state_2d_t) :: completed, state
    complex(rk) :: expected_gap(3, 3), sampled_gap(3, 3)
    real(rk) :: expected_current(3), position(2), sampled_current(3)
    real(rk) :: scaled_radius_squared
    logical, allocatable :: active(:)
    integer :: point, x_node, y_node

    call make_uniform_cartesian_mesh( &
      -5.0_rk, 5.0_rk, 40, -5.0_rk, 5.0_rk, 40, mesh)
    endpoint = make_endpoint()
    endpoint%outer_radius = 9.0_rk
    call set_state_asymptotic_elliptical_fit( &
      endpoint, [1.0_rk, 1.5_rk], [2.0_rk, 3.0_rk])
    call allocate_spinful_state_2d(mesh, state)
    do point = 1, mesh%point_count()
      position = mesh%point_coordinate(point)
      call evaluate_known_elliptical_field( &
        position, endpoint, state%order_parameter(:, :, point), &
        state%current_mean_field(:, point))
    end do
    call require(endpoint%is_valid_for(mesh), &
      "elliptical state endpoint was unexpectedly invalid")

    position = [6.0_rk, 0.0_rk]
    call sample_free_vortex_asymptotic_state_2d( &
      mesh, state, position, endpoint, sampled_gap, sampled_current)
    call evaluate_known_elliptical_field( &
      position, endpoint, expected_gap, expected_current)
    call require(maxval(abs(sampled_gap - expected_gap)) < 3.0e-14_rk .and. &
                 maxval(abs(sampled_current - expected_current)) < &
                   3.0e-14_rk, &
      "elliptical inverse-radius fit failed on the x axis")

    position = [0.0_rk, 6.0_rk]
    call sample_free_vortex_asymptotic_state_2d( &
      mesh, state, position, endpoint, sampled_gap, sampled_current)
    call evaluate_known_elliptical_field( &
      position, endpoint, expected_gap, expected_current)
    call require(maxval(abs(sampled_gap - expected_gap)) < 3.0e-14_rk .and. &
                 maxval(abs(sampled_current - expected_current)) < &
                   3.0e-14_rk, &
      "elliptical inverse-radius fit failed on the y axis")

    allocate(active(mesh%point_count()), source=.false.)
    do point = 1, mesh%point_count()
      position = mesh%point_coordinate(point)
      scaled_radius_squared = (position(1) / 2.5_rk)**2 + &
        (position(2) / 3.5_rk)**2
      active(point) = scaled_radius_squared <= 1.0_rk
    end do
    x_node = mesh%x_cell_count()
    y_node = mesh%y_cell_count() / 2
    point = mesh%point_index(x_node, y_node)
    state%order_parameter(:, :, point) = cmplx(91.0_rk, -11.0_rk, rk)
    state%current_mean_field(:, point) = 19.0_rk
    call apply_free_vortex_asymptotic_halo_2d( &
      mesh, state, active, endpoint, completed)
    position = mesh%point_coordinate(point)
    call evaluate_known_elliptical_field( &
      position, endpoint, expected_gap, expected_current)
    call require(maxval(abs(completed%order_parameter(:, :, point) - &
                            expected_gap)) < 3.0e-14_rk .and. &
                 maxval(abs(completed%current_mean_field(:, point) - &
                            expected_current)) < 3.0e-14_rk, &
      "elliptical dependent halo was not regenerated from the fit")
  end subroutine test_elliptical_fit_surfaces


  subroutine test_fitted_rotated_bulk_reference()
    type(cartesian_mesh_2d_t) :: mesh
    type(free_vortex_endpoint_2d_t) :: endpoint, fixed_endpoint
    type(spinful_state_2d_t) :: state
    complex(rk) :: expected(3, 3), phase_factor, sampled(3, 3)
    real(rk) :: angle, current(3), fit_maximum, fit_rms, known_phase
    real(rk) :: known_rotation(3, 3), position(2)
    integer :: point

    call make_uniform_cartesian_mesh( &
      -4.0_rk, 4.0_rk, 160, -4.0_rk, 4.0_rk, 160, mesh)
    call allocate_spinful_state_2d(mesh, state)
    known_phase = 0.27_rk
    known_rotation = 0.0_rk
    known_rotation(1, 1) = cos(0.35_rk)
    known_rotation(1, 2) = -sin(0.35_rk)
    known_rotation(2, 1) = sin(0.35_rk)
    known_rotation(2, 2) = cos(0.35_rk)
    known_rotation(3, 3) = 1.0_rk
    endpoint = make_endpoint()
    endpoint%outer_radius = 7.0_rk
    do point = 1, mesh%point_count()
      position = mesh%point_coordinate(point)
      if (maxval(abs(position)) > 0.0_rk) then
        angle = atan2(position(2), position(1))
      else
        angle = 0.0_rk
      end if
      phase_factor = cmplx(endpoint%bulk_gap, 0.0_rk, rk) * &
        exp(cmplx(0.0_rk, endpoint%winding * angle + known_phase, rk))
      state%order_parameter(:, :, point) = phase_factor * &
        cmplx(known_rotation, 0.0_rk, rk)
    end do

    fixed_endpoint = make_endpoint()
    fixed_endpoint%outer_radius = 7.0_rk
    call fit_state_bulk_reference( &
      fixed_endpoint, mesh, state, 3.5_rk, 192, fit_rms, fit_maximum, &
      update_reference=.false.)
    call require(abs(fixed_endpoint%bulk_phase_offset) < 2.0e-15_rk .and. &
                 maxval(abs(fixed_endpoint%bulk_rotation - &
                   reshape([1.0_rk, 0.0_rk, 0.0_rk, &
                            0.0_rk, 1.0_rk, 0.0_rk, &
                            0.0_rk, 0.0_rk, 1.0_rk], [3, 3]))) < &
                   2.0e-15_rk, &
      "fixed bulk-reference diagnostic changed A0")

    call fit_state_bulk_reference( &
      endpoint, mesh, state, 3.5_rk, 192, fit_rms, fit_maximum)
    call require(abs(endpoint%bulk_phase_offset - known_phase) < 2.0e-5_rk, &
      "bulk-reference fit did not recover the common phase")
    call require(maxval(abs(endpoint%bulk_rotation - known_rotation)) < &
                   2.0e-5_rk, &
      "bulk-reference fit did not recover the spin-orbit rotation")
    call require(fit_rms < 2.0e-4_rk .and. fit_maximum < 5.0e-4_rk, &
      "bulk-reference fit has an unexpectedly large interpolation residual")

    call set_state_asymptotic_fit(endpoint, 2.5_rk, 3.5_rk)
    call require(endpoint%is_valid_for(mesh), &
      "fitted rotated endpoint was unexpectedly invalid")
    position = [6.0_rk, 0.0_rk]
    call sample_free_vortex_asymptotic_state_2d( &
      mesh, state, position, endpoint, sampled, current)
    expected = cmplx(endpoint%bulk_gap, 0.0_rk, rk) * &
      exp(cmplx(0.0_rk, known_phase, rk)) * &
      cmplx(known_rotation, 0.0_rk, rk)
    call require(maxval(abs(sampled - expected)) < 3.0e-5_rk, &
      "asymptotic continuation lost the fitted rotated bulk reference")
  end subroutine test_fitted_rotated_bulk_reference


  subroutine evaluate_known_two_term_field( &
      position, endpoint, order_parameter, current_mean_field)
    real(rk), intent(in) :: position(2)
    type(free_vortex_endpoint_2d_t), intent(in) :: endpoint
    complex(rk), intent(out) :: order_parameter(3, 3)
    real(rk), intent(out) :: current_mean_field(3)

    complex(rk) :: coefficient_one(3, 3), coefficient_two(3, 3)
    complex(rk) :: bulk_value
    real(rk) :: angle, current_one(3), current_two(3), radius, ratio

    radius = sqrt(dot_product(position, position))
    if (radius <= 64.0_rk * epsilon(1.0_rk)) then
      order_parameter = cmplx(0.0_rk, 0.0_rk, rk)
      current_mean_field = 0.0_rk
      return
    end if
    angle = atan2(position(2), position(1))
    bulk_value = cmplx(endpoint%bulk_gap, 0.0_rk, rk) * &
      exp(cmplx(0.0_rk, endpoint%winding * angle + &
        endpoint%bulk_phase_offset, rk))
    order_parameter = bulk_value * &
      cmplx(endpoint%bulk_rotation, 0.0_rk, rk)

    coefficient_one = cmplx(0.0_rk, 0.0_rk, rk)
    coefficient_two = cmplx(0.0_rk, 0.0_rk, rk)
    coefficient_one(1, 3) = cmplx( &
      0.08_rk * (1.0_rk + 0.25_rk * cos(angle)), &
      0.03_rk * sin(angle), rk)
    coefficient_one(3, 1) = cmplx( &
      -0.04_rk * sin(angle), 0.02_rk * cos(angle), rk)
    coefficient_two(1, 1) = cmplx( &
      0.05_rk * cos(2.0_rk * angle), 0.01_rk * sin(angle), rk)
    coefficient_two(2, 2) = cmplx(-0.03_rk, 0.02_rk * sin(angle), rk)
    ratio = endpoint%matching_radius / radius
    order_parameter = order_parameter + &
      cmplx(ratio, 0.0_rk, rk) * coefficient_one + &
      cmplx(ratio**2, 0.0_rk, rk) * coefficient_two

    current_one = [0.05_rk * cos(angle), 0.05_rk * sin(angle), &
                   0.01_rk * cos(2.0_rk * angle)]
    current_two = [0.02_rk * sin(angle), -0.03_rk * cos(angle), 0.015_rk]
    current_mean_field = ratio * current_one + ratio**2 * current_two
  end subroutine evaluate_known_two_term_field


  subroutine evaluate_known_elliptical_field( &
      position, endpoint, order_parameter, current_mean_field)
    real(rk), intent(in) :: position(2)
    type(free_vortex_endpoint_2d_t), intent(in) :: endpoint
    complex(rk), intent(out) :: order_parameter(3, 3)
    real(rk), intent(out) :: current_mean_field(3)

    complex(rk) :: bulk_value
    real(rk) :: angle, current_one(3), current_two(3), direction(2)
    real(rk) :: matching_radius, radius, ratio

    radius = sqrt(dot_product(position, position))
    if (radius <= 64.0_rk * epsilon(1.0_rk)) then
      order_parameter = cmplx(0.0_rk, 0.0_rk, rk)
      current_mean_field = 0.0_rk
      return
    end if
    direction = position / radius
    angle = atan2(direction(2), direction(1))
    matching_radius = 1.0_rk / sqrt( &
      (direction(1) / endpoint%matching_radii(1))**2 + &
      (direction(2) / endpoint%matching_radii(2))**2)
    ratio = matching_radius / radius
    bulk_value = cmplx(endpoint%bulk_gap, 0.0_rk, rk) * &
      exp(cmplx(0.0_rk, endpoint%winding * angle, rk))
    order_parameter = cmplx(0.0_rk, 0.0_rk, rk)
    order_parameter(1, 1) = bulk_value + &
      ratio**2 * cmplx(0.06_rk, -0.02_rk, rk)
    order_parameter(2, 2) = bulk_value
    order_parameter(3, 3) = bulk_value
    order_parameter(1, 3) = &
      ratio * cmplx(-0.03_rk, 0.02_rk, rk)
    current_one = [0.07_rk, -0.04_rk, 0.03_rk]
    current_two = [0.02_rk, 0.05_rk, -0.01_rk]
    current_mean_field = ratio * current_one + ratio**2 * current_two
  end subroutine evaluate_known_elliptical_field


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
