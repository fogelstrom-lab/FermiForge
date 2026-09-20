module free_vortex_asymptotic_2d
  use he3_kinds, only : rk
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t
  use spinful_state_2d, only : spinful_state_2d_t
  use radial_mesh, only : radial_mesh_t
  use axial_radial_embedding, only : axial_radial_profile_t, &
    sample_axial_radial_profile
  use bilinear_field_sampler_2d, only : sample_spinful_state_2d
  use straight_trajectory_2d, only : straight_trajectory_2d_t
  use trajectory_self_energy_2d, only : trajectory_self_energy_2d_t
  implicit none
  private

  complex(rk), parameter :: imaginary_unit = &
    cmplx(0.0_rk, 1.0_rk, kind=rk)
  real(rk), parameter :: identity_rotation(3, 3) = reshape([ &
    1.0_rk, 0.0_rk, 0.0_rk, &
    0.0_rk, 1.0_rk, 0.0_rk, &
    0.0_rk, 0.0_rk, 1.0_rk], [3, 3])

  type, public :: free_vortex_endpoint_2d_t
    logical :: enabled = .false.
    logical :: use_axisymmetric_radial_reference = .false.
    logical :: use_inverse_radius_fit = .false.
    logical :: use_elliptical_fit_surfaces = .false.
    real(rk) :: center(2) = 0.0_rk
    real(rk) :: bulk_gap = 0.0_rk
    real(rk) :: winding = 1.0_rk
    real(rk) :: bulk_phase_offset = 0.0_rk
    real(rk) :: bulk_rotation(3, 3) = identity_rotation
    real(rk) :: outer_radius = 0.0_rk
    real(rk) :: maximum_step = 0.0_rk
    real(rk) :: fit_inner_radius = 0.0_rk
    real(rk) :: matching_radius = 0.0_rk
    real(rk) :: fit_inner_radii(2) = 0.0_rk
    real(rk) :: matching_radii(2) = 0.0_rk
    type(radial_mesh_t) :: radial_grid
    type(axial_radial_profile_t) :: radial_profile
  contains
    procedure :: is_valid_for => free_vortex_endpoint_is_valid_for
  end type free_vortex_endpoint_2d_t

  public :: extend_with_free_vortex_asymptotic_2d
  public :: sample_free_vortex_asymptotic_state_2d
  public :: set_axisymmetric_radial_reference
  public :: set_axisymmetric_asymptotic_fit
  public :: set_state_asymptotic_fit
  public :: set_state_asymptotic_elliptical_fit
  public :: fit_state_bulk_reference
  public :: apply_free_vortex_asymptotic_halo_2d

contains

  subroutine set_axisymmetric_radial_reference( &
      endpoint, radial_grid, radial_profile, matching_radius)
    type(free_vortex_endpoint_2d_t), intent(inout) :: endpoint
    type(radial_mesh_t), intent(in) :: radial_grid
    type(axial_radial_profile_t), intent(in) :: radial_profile
    real(rk), intent(in) :: matching_radius

    if (.not. radial_profile%is_valid_for(radial_grid)) &
      error stop "cannot configure an invalid radial endpoint reference"
    if (matching_radius <= 0.0_rk .or. &
        matching_radius > radial_grid%outer_radius()) &
      error stop "radial endpoint matching radius is outside its profile"
    endpoint%use_axisymmetric_radial_reference = .true.
    endpoint%use_inverse_radius_fit = .false.
    endpoint%use_elliptical_fit_surfaces = .false.
    endpoint%fit_inner_radius = 0.0_rk
    endpoint%matching_radius = matching_radius
    endpoint%radial_grid = radial_grid
    endpoint%radial_profile = radial_profile
  end subroutine set_axisymmetric_radial_reference


  subroutine set_axisymmetric_asymptotic_fit( &
      endpoint, radial_grid, radial_profile, fit_inner_radius, matching_radius)
    type(free_vortex_endpoint_2d_t), intent(inout) :: endpoint
    type(radial_mesh_t), intent(in) :: radial_grid
    type(axial_radial_profile_t), intent(in) :: radial_profile
    real(rk), intent(in) :: fit_inner_radius, matching_radius

    if (.not. radial_profile%is_valid_for(radial_grid)) &
      error stop "cannot fit an invalid radial endpoint reference"
    call require_fit_radii( &
      fit_inner_radius, matching_radius, radial_grid%outer_radius())
    endpoint%use_axisymmetric_radial_reference = .true.
    endpoint%use_inverse_radius_fit = .true.
    endpoint%use_elliptical_fit_surfaces = .false.
    endpoint%fit_inner_radius = fit_inner_radius
    endpoint%matching_radius = matching_radius
    endpoint%radial_grid = radial_grid
    endpoint%radial_profile = radial_profile
  end subroutine set_axisymmetric_asymptotic_fit


  subroutine set_state_asymptotic_fit( &
      endpoint, fit_inner_radius, matching_radius)
    type(free_vortex_endpoint_2d_t), intent(inout) :: endpoint
    real(rk), intent(in) :: fit_inner_radius, matching_radius

    call require_fit_radii( &
      fit_inner_radius, matching_radius, endpoint%outer_radius)
    endpoint%use_axisymmetric_radial_reference = .false.
    endpoint%use_inverse_radius_fit = .true.
    endpoint%use_elliptical_fit_surfaces = .false.
    endpoint%fit_inner_radius = fit_inner_radius
    endpoint%matching_radius = matching_radius
  end subroutine set_state_asymptotic_fit


  subroutine set_state_asymptotic_elliptical_fit( &
      endpoint, fit_inner_radii, matching_radii)
    type(free_vortex_endpoint_2d_t), intent(inout) :: endpoint
    real(rk), intent(in) :: fit_inner_radii(2), matching_radii(2)

    if (any(fit_inner_radii <= 0.0_rk) .or. &
        any(matching_radii <= fit_inner_radii)) &
      error stop "elliptical inverse-radius fit radii are invalid"
    endpoint%use_axisymmetric_radial_reference = .false.
    endpoint%use_inverse_radius_fit = .true.
    endpoint%use_elliptical_fit_surfaces = .true.
    endpoint%fit_inner_radii = fit_inner_radii
    endpoint%matching_radii = matching_radii
    endpoint%fit_inner_radius = maxval(fit_inner_radii)
    endpoint%matching_radius = maxval(matching_radii)
  end subroutine set_state_asymptotic_elliptical_fit


  subroutine fit_state_bulk_reference( &
      endpoint, mesh, state, fitting_radius, sample_count, &
      rms_error, maximum_error, fitting_radii, update_reference)
    type(free_vortex_endpoint_2d_t), intent(inout) :: endpoint
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    real(rk), intent(in) :: fitting_radius
    integer, intent(in) :: sample_count
    real(rk), intent(out), optional :: rms_error, maximum_error
    real(rk), intent(in), optional :: fitting_radii(2)
    logical, intent(in), optional :: update_reference

    complex(rk) :: average(3, 3), normalized(3, 3), phase_factor
    complex(rk) :: phase_square, sampled_order_parameter(3, 3)
    real(rk) :: angle, candidate(3, 3), current_mean_field(3)
    real(rk) :: determinant, error_sum, local_maximum, position(2)
    real(rk) :: sample_radius
    real(rk) :: projected(3, 3), pi
    integer :: sample
    logical :: inside, should_update_reference

    if (.not. mesh%is_valid() .or. .not. state%is_valid_for(mesh)) &
      error stop "bulk-reference fit requires a valid field"
    if (.not. endpoint%enabled .or. endpoint%bulk_gap <= 0.0_rk) &
      error stop "bulk-reference fit requires an enabled endpoint and gap"
    if (sample_count < 8) &
      error stop "bulk-reference fit requires at least eight angular samples"
    if (present(fitting_radii)) then
      if (any(fitting_radii <= 0.0_rk) .or. &
          fitting_radii(1) >= min(endpoint%center(1) - mesh%x_minimum, &
                                  mesh%x_maximum - endpoint%center(1)) .or. &
          fitting_radii(2) >= min(endpoint%center(2) - mesh%y_minimum, &
                                  mesh%y_maximum - endpoint%center(2))) &
        error stop "bulk-reference fitting ellipse lies outside the field"
    else
      if (fitting_radius <= 0.0_rk .or. &
          fitting_radius >= minimum_edge_distance(mesh, endpoint%center)) &
        error stop "bulk-reference fitting circle lies outside the field"
    end if

    should_update_reference = .true.
    if (present(update_reference)) &
      should_update_reference = update_reference
    pi = acos(-1.0_rk)
    if (should_update_reference) then
      average = cmplx(0.0_rk, 0.0_rk, kind=rk)
      do sample = 0, sample_count - 1
        angle = 2.0_rk * pi * real(sample, rk) / real(sample_count, rk)
        if (present(fitting_radii)) then
          sample_radius = elliptical_radius_along_direction( &
            [cos(angle), sin(angle)], fitting_radii)
        else
          sample_radius = fitting_radius
        end if
        position = endpoint%center + sample_radius * [cos(angle), sin(angle)]
        call sample_spinful_state_2d( &
          mesh, state, position(1), position(2), sampled_order_parameter, &
          current_mean_field, inside)
        if (.not. inside) &
          error stop "bulk-reference fitting sample is outside the field"
        phase_factor = exp(-imaginary_unit * cmplx( &
          endpoint%winding * angle, 0.0_rk, kind=rk))
        average = average + phase_factor * sampled_order_parameter / &
          cmplx(endpoint%bulk_gap, 0.0_rk, kind=rk)
      end do
      average = average / cmplx(real(sample_count, rk), 0.0_rk, kind=rk)

      ! A B-phase bulk matrix has the form exp(i phase) R with real R in SO(3).
      ! Squaring all matrix elements removes the pi ambiguity in their common
      ! complex phase.  The determinant chooses between R and -R.
      phase_square = sum(average * average)
      if (abs(phase_square) <= 128.0_rk * epsilon(1.0_rk)) &
        error stop "bulk-reference phase cannot be inferred from the field"
      endpoint%bulk_phase_offset = 0.5_rk * atan2( &
        aimag(phase_square), real(phase_square, rk))
      phase_factor = exp(-imaginary_unit * &
        cmplx(endpoint%bulk_phase_offset, 0.0_rk, kind=rk))
      candidate = real(phase_factor * average, rk)
      determinant = determinant_3x3(candidate)
      if (determinant < 0.0_rk) then
        endpoint%bulk_phase_offset = endpoint%bulk_phase_offset + pi
        candidate = -candidate
      end if
      call project_to_proper_rotation(candidate, projected)
      endpoint%bulk_rotation = projected
    end if

    error_sum = 0.0_rk
    local_maximum = 0.0_rk
    phase_factor = exp(imaginary_unit * &
      cmplx(endpoint%bulk_phase_offset, 0.0_rk, kind=rk))
    do sample = 0, sample_count - 1
      angle = 2.0_rk * pi * real(sample, rk) / real(sample_count, rk)
      if (present(fitting_radii)) then
        sample_radius = elliptical_radius_along_direction( &
          [cos(angle), sin(angle)], fitting_radii)
      else
        sample_radius = fitting_radius
      end if
      position = endpoint%center + sample_radius * [cos(angle), sin(angle)]
      call sample_spinful_state_2d( &
        mesh, state, position(1), position(2), sampled_order_parameter, &
        current_mean_field, inside)
      normalized = exp(-imaginary_unit * cmplx( &
        endpoint%winding * angle, 0.0_rk, kind=rk)) * &
        sampled_order_parameter / cmplx(endpoint%bulk_gap, 0.0_rk, kind=rk)
      normalized = normalized - phase_factor * &
        cmplx(endpoint%bulk_rotation, 0.0_rk, kind=rk)
      error_sum = error_sum + sum(abs(normalized)**2)
      local_maximum = max(local_maximum, maxval(abs(normalized)))
    end do
    if (present(rms_error)) rms_error = sqrt( &
      error_sum / real(9 * sample_count, rk))
    if (present(maximum_error)) maximum_error = local_maximum
  end subroutine fit_state_bulk_reference

  subroutine sample_free_vortex_asymptotic_state_2d( &
      mesh, state, position, endpoint, order_parameter, current_mean_field)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    real(rk), intent(in) :: position(2)
    type(free_vortex_endpoint_2d_t), intent(in) :: endpoint
    complex(rk), intent(out) :: order_parameter(3, 3)
    real(rk), intent(out) :: current_mean_field(3)

    complex(rk) :: boundary_order_parameter(3, 3)
    complex(rk) :: bulk_order_parameter(3, 3)
    complex(rk) :: fit_outer_order_parameter(3, 3)
    real(rk) :: boundary_current(3), boundary_position(2)
    real(rk) :: fit_inner_current(3), fit_outer_current(3)
    real(rk) :: decay, direction(2), displacement(2), fit_inner_radius
    real(rk) :: matching_radius
    real(rk) :: phase_angle, radius, tolerance
    integer :: orbital, spin
    logical :: inside

    if (.not. endpoint%is_valid_for(mesh)) &
      error stop "invalid free-vortex endpoint condition"
    if (.not. state%is_valid_for(mesh)) &
      error stop "free-vortex endpoint state does not match its mesh"

    displacement = position - endpoint%center
    radius = sqrt(dot_product(displacement, displacement))
    tolerance = 128.0_rk * epsilon(1.0_rk) * &
      max(1.0_rk, endpoint%outer_radius)
    if (radius <= tolerance) &
      error stop "free-vortex asymptotic state is undefined at its center"

    direction = displacement / radius
    phase_angle = atan2(direction(2), direction(1))
    bulk_order_parameter = bulk_order_parameter_at_angle(endpoint, phase_angle)

    if (endpoint%use_inverse_radius_fit) then
      call fit_radii_along_direction( &
        endpoint, direction, fit_inner_radius, matching_radius)
      if (radius <= matching_radius + tolerance) then
        if (endpoint%use_axisymmetric_radial_reference) then
          call sample_axial_radial_profile( &
            endpoint%radial_grid, endpoint%radial_profile, &
            displacement(1), displacement(2), endpoint%winding, &
            order_parameter, current_mean_field, inside)
        else
          call sample_spinful_state_2d( &
            mesh, state, position(1), position(2), order_parameter, &
            current_mean_field, inside)
        end if
        if (.not. inside) &
          error stop "asymptotic fit received an unavailable interior sample"
        return
      end if

      ! The inner order-parameter sample is deliberately discarded.  Equation
      ! (20) of dcvlong assigns a single radial power to each Cartesian block,
      ! so its coefficient is fixed by continuity at the matching surface.
      ! The Fermi-liquid field can contain both powers and still needs the two
      ! radii.
      call sample_fit_source(fit_inner_radius, &
        fit_outer_order_parameter, fit_inner_current)
      call sample_fit_source(matching_radius, &
        fit_outer_order_parameter, fit_outer_current)
      call evaluate_inverse_radius_fit( &
        radius, phase_angle, fit_inner_radius, matching_radius, &
        fit_outer_order_parameter, &
        fit_inner_current, fit_outer_current, &
        order_parameter, current_mean_field)
      return
    end if

    if (endpoint%use_axisymmetric_radial_reference) then
      if (radius <= endpoint%matching_radius + tolerance) then
        call sample_axial_radial_profile( &
          endpoint%radial_grid, endpoint%radial_profile, &
          displacement(1), displacement(2), endpoint%winding, &
          order_parameter, current_mean_field, inside)
        if (.not. inside) &
          error stop "radial endpoint sample lies outside its source profile"
        return
      end if
      matching_radius = endpoint%matching_radius
      boundary_position = matching_radius * direction
      call sample_axial_radial_profile( &
        endpoint%radial_grid, endpoint%radial_profile, &
        boundary_position(1), boundary_position(2), endpoint%winding, &
        boundary_order_parameter, boundary_current, inside)
      if (.not. inside) &
        error stop "radial endpoint match lies outside its source profile"
    else
      matching_radius = radial_distance_to_rectangle( &
        mesh, endpoint%center, direction)
      if (radius < matching_radius - tolerance) &
        error stop "free-vortex asymptotic sample lies inside the field mesh"
      boundary_position = endpoint%center + matching_radius * direction
      call sample_spinful_state_2d( &
        mesh, state, boundary_position(1), boundary_position(2), &
        boundary_order_parameter, boundary_current, inside)
      if (.not. inside) &
        error stop "free-vortex radial match point lies outside the field mesh"
    end if

    order_parameter = cmplx(0.0_rk, 0.0_rk, kind=rk)
    do spin = 1, 3
      do orbital = 1, 3
        decay = matching_radius / radius
        ! This is the explicit Cartesian form of the free-vortex tail used
        ! by new_src/interpol.f90: mixed z/in-plane components decay as 1/r;
        ! all other departures from bulk B phase decay as 1/r^2.
        if ((spin == 3) .eqv. (orbital == 3)) decay = decay * decay
        order_parameter(spin, orbital) = &
          cmplx(decay, 0.0_rk, kind=rk) * &
          boundary_order_parameter(spin, orbital)
        order_parameter(spin, orbital) = order_parameter(spin, orbital) + &
          cmplx(1.0_rk - decay, 0.0_rk, kind=rk) * &
          bulk_order_parameter(spin, orbital)
      end do
    end do
    current_mean_field = (matching_radius / radius) * boundary_current

  contains

    subroutine sample_fit_source( &
        source_radius, source_order_parameter, source_current)
      real(rk), intent(in) :: source_radius
      complex(rk), intent(out) :: source_order_parameter(3, 3)
      real(rk), intent(out) :: source_current(3)

      real(rk) :: source_position(2)
      logical :: source_inside

      source_position = source_radius * direction
      if (endpoint%use_axisymmetric_radial_reference) then
        call sample_axial_radial_profile( &
          endpoint%radial_grid, endpoint%radial_profile, &
          source_position(1), source_position(2), endpoint%winding, &
          source_order_parameter, source_current, source_inside)
      else
        source_position = endpoint%center + source_position
        call sample_spinful_state_2d( &
          mesh, state, source_position(1), source_position(2), &
          source_order_parameter, source_current, source_inside)
      end if
      if (.not. source_inside) &
        error stop "inverse-radius fit sample lies outside its source field"
    end subroutine sample_fit_source


    subroutine evaluate_inverse_radius_fit( &
        evaluation_radius, angle, inner_radius, outer_radius, &
        outer_order_parameter, &
        inner_current, outer_current, &
        fitted_order_parameter, fitted_current)
      real(rk), intent(in) :: evaluation_radius, angle
      real(rk), intent(in) :: inner_radius, outer_radius
      complex(rk), intent(in) :: outer_order_parameter(3, 3)
      real(rk), intent(in) :: inner_current(3), outer_current(3)
      complex(rk), intent(out) :: fitted_order_parameter(3, 3)
      real(rk), intent(out) :: fitted_current(3)

      complex(rk) :: bulk_order_parameter(3, 3)
      complex(rk) :: coefficient_one(3, 3), coefficient_two(3, 3)
      complex(rk) :: outer_departure(3, 3)
      real(rk) :: current_one(3), current_two(3), denominator
      real(rk) :: fit_ratio, radial_ratio
      integer :: orbital, spin

      bulk_order_parameter = bulk_order_parameter_at_angle(endpoint, angle)

      fit_ratio = outer_radius / inner_radius
      denominator = fit_ratio * (fit_ratio - 1.0_rk)
      outer_departure = outer_order_parameter - bulk_order_parameter
      coefficient_one = cmplx(0.0_rk, 0.0_rk, kind=rk)
      coefficient_two = cmplx(0.0_rk, 0.0_rk, kind=rk)
      do spin = 1, 3
        do orbital = 1, 3
          ! dcvlong Eqs. (19)-(20): mixed z/in-plane elements contain the
          ! spin-rotation 1/r tail.  The in-plane block and zz contain only
          ! the 1/r^2 amplitude correction.  The coefficient functions remain
          ! angle dependent through outer_departure.
          if ((spin == 3) .neqv. (orbital == 3)) then
            coefficient_one(spin, orbital) = &
              outer_departure(spin, orbital)
          else
            coefficient_two(spin, orbital) = &
              outer_departure(spin, orbital)
          end if
        end do
      end do
      current_two = (inner_current - fit_ratio * outer_current) / denominator
      current_one = outer_current - current_two

      radial_ratio = outer_radius / evaluation_radius
      fitted_order_parameter = bulk_order_parameter + &
        cmplx(radial_ratio, 0.0_rk, kind=rk) * coefficient_one + &
        cmplx(radial_ratio**2, 0.0_rk, kind=rk) * coefficient_two
      fitted_current = radial_ratio * current_one + &
        radial_ratio**2 * current_two
    end subroutine evaluate_inverse_radius_fit
  end subroutine sample_free_vortex_asymptotic_state_2d


  subroutine apply_free_vortex_asymptotic_halo_2d( &
      mesh, source_state, active_point_mask, endpoint, completed_state)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: source_state
    logical, intent(in) :: active_point_mask(:)
    type(free_vortex_endpoint_2d_t), intent(in) :: endpoint
    type(spinful_state_2d_t), intent(out) :: completed_state

    real(rk) :: current_mean_field(3), direction(2), displacement(2)
    real(rk) :: fit_inner_radius, matching_radius, radius, tolerance
    complex(rk) :: order_parameter(3, 3)
    integer :: point

    if (.not. endpoint%use_inverse_radius_fit) &
      error stop "cannot construct a halo without an inverse-radius fit"
    if (.not. endpoint%is_valid_for(mesh)) &
      error stop "cannot construct a halo from an invalid endpoint"
    if (.not. source_state%is_valid_for(mesh)) &
      error stop "asymptotic halo source does not match its mesh"
    if (size(active_point_mask) /= mesh%point_count() .or. &
        .not. any(active_point_mask)) &
      error stop "asymptotic halo has an invalid active-point mask"

    completed_state = source_state
    do point = 1, mesh%point_count()
      if (active_point_mask(point)) cycle
      displacement = mesh%point_coordinate(point) - endpoint%center
      radius = sqrt(dot_product(displacement, displacement))
      tolerance = 128.0_rk * epsilon(1.0_rk) * &
        max(1.0_rk, endpoint%matching_radius)
      if (radius <= tolerance) &
        error stop "asymptotic halo cannot include the vortex center"
      direction = displacement / radius
      call fit_radii_along_direction( &
        endpoint, direction, fit_inner_radius, matching_radius)
      if (radius <= matching_radius + tolerance) &
        error stop "asymptotic halo overlaps the fit-source surface"
      call sample_free_vortex_asymptotic_state_2d( &
        mesh, source_state, mesh%point_coordinate(point), endpoint, &
        order_parameter, current_mean_field)
      completed_state%order_parameter(:, :, point) = order_parameter
      completed_state%current_mean_field(:, point) = current_mean_field
    end do
  end subroutine apply_free_vortex_asymptotic_halo_2d


  subroutine extend_with_free_vortex_asymptotic_2d( &
      mesh, state, trajectory, interior_self_energy, feedback_scale, endpoint, &
      path_coordinate, target_sample, extended_self_energy)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    type(straight_trajectory_2d_t), intent(in) :: trajectory
    type(trajectory_self_energy_2d_t), intent(in) :: interior_self_energy
    real(rk), intent(in) :: feedback_scale
    type(free_vortex_endpoint_2d_t), intent(in) :: endpoint
    real(rk), allocatable, intent(out) :: path_coordinate(:)
    integer, intent(out) :: target_sample
    type(trajectory_self_energy_2d_t), intent(out) :: extended_self_energy

    complex(rk) :: order_parameter(3, 3), pair_potential(3)
    real(rk) :: current_mean_field(3), discriminant, distance
    real(rk) :: offset(2), outer_entry, outer_exit, planar_norm_squared
    real(rk) :: position(2), projection
    integer :: entry_intervals, index, interior_count, orbital
    integer :: output_index, exit_intervals, total_count

    if (.not. endpoint%is_valid_for(mesh)) &
      error stop "cannot extend with an invalid free-vortex endpoint"
    if (.not. state%is_valid_for(mesh)) &
      error stop "free-vortex extension state does not match its mesh"
    if (.not. trajectory%is_valid()) &
      error stop "cannot extend an invalid straight trajectory"
    if (.not. interior_self_energy%is_valid() .or. &
        interior_self_energy%sample_count() /= trajectory%sample_count()) &
      error stop "free-vortex extension has incompatible self energies"

    interior_count = trajectory%sample_count()
    planar_norm_squared = dot_product( &
      trajectory%momentum(1:2), trajectory%momentum(1:2))
    if (planar_norm_squared <= 128.0_rk * epsilon(1.0_rk)) then
      allocate(path_coordinate(interior_count), &
        source=trajectory%path_coordinate)
      allocate(extended_self_energy%triplet_pair_potential(3, interior_count), &
        source=interior_self_energy%triplet_pair_potential)
      allocate(extended_self_energy%diagonal_shift(interior_count), &
        source=interior_self_energy%diagonal_shift)
      target_sample = trajectory%target_sample
      return
    end if

    offset = trajectory%origin - endpoint%center
    projection = dot_product(offset, trajectory%momentum(1:2))
    discriminant = projection * projection - planar_norm_squared * &
      (dot_product(offset, offset) - endpoint%outer_radius**2)
    if (discriminant <= 0.0_rk) &
      error stop "free-vortex outer circle does not cross the trajectory"
    outer_entry = (-projection - sqrt(discriminant)) / planar_norm_squared
    outer_exit = (-projection + sqrt(discriminant)) / planar_norm_squared
    if (outer_entry > trajectory%entry_coordinate .or. &
        outer_exit < trajectory%exit_coordinate) &
      error stop "free-vortex outer circle does not enclose the field mesh"

    distance = trajectory%entry_coordinate - outer_entry
    entry_intervals = max(1, ceiling(distance / endpoint%maximum_step))
    distance = outer_exit - trajectory%exit_coordinate
    exit_intervals = max(1, ceiling(distance / endpoint%maximum_step))
    total_count = entry_intervals + interior_count + exit_intervals
    allocate(path_coordinate(total_count))
    allocate(extended_self_energy%triplet_pair_potential(3, total_count), &
      source=cmplx(0.0_rk, 0.0_rk, kind=rk))
    allocate(extended_self_energy%diagonal_shift(total_count), source=0.0_rk)

    do index = 0, entry_intervals - 1
      output_index = index + 1
      path_coordinate(output_index) = outer_entry + &
        (trajectory%entry_coordinate - outer_entry) * &
        real(index, rk) / real(entry_intervals, rk)
      call sample_tail_self_energy(output_index)
    end do

    output_index = entry_intervals
    path_coordinate(output_index + 1:output_index + interior_count) = &
      trajectory%path_coordinate
    extended_self_energy%triplet_pair_potential( &
      :, output_index + 1:output_index + interior_count) = &
      interior_self_energy%triplet_pair_potential
    extended_self_energy%diagonal_shift( &
      output_index + 1:output_index + interior_count) = &
      interior_self_energy%diagonal_shift
    target_sample = output_index + trajectory%target_sample

    output_index = entry_intervals + interior_count
    do index = 1, exit_intervals
      path_coordinate(output_index + index) = trajectory%exit_coordinate + &
        (outer_exit - trajectory%exit_coordinate) * &
        real(index, rk) / real(exit_intervals, rk)
      call sample_tail_self_energy(output_index + index)
    end do

  contains

    subroutine sample_tail_self_energy(sample)
      integer, intent(in) :: sample

      position = trajectory%origin + path_coordinate(sample) * &
        trajectory%momentum(1:2)
      call sample_free_vortex_asymptotic_state_2d( &
        mesh, state, position, endpoint, order_parameter, current_mean_field)
      pair_potential = cmplx(0.0_rk, 0.0_rk, kind=rk)
      do orbital = 1, 3
        pair_potential = pair_potential + &
          cmplx(trajectory%momentum(orbital), 0.0_rk, kind=rk) * &
          order_parameter(:, orbital)
      end do
      extended_self_energy%triplet_pair_potential(:, sample) = pair_potential
      extended_self_energy%diagonal_shift(sample) = feedback_scale * &
        dot_product(trajectory%momentum, current_mean_field)
    end subroutine sample_tail_self_energy

  end subroutine extend_with_free_vortex_asymptotic_2d


  pure function bulk_order_parameter_at_angle(endpoint, angle) result(bulk)
    type(free_vortex_endpoint_2d_t), intent(in) :: endpoint
    real(rk), intent(in) :: angle
    complex(rk) :: bulk(3, 3)

    complex(rk) :: amplitude

    amplitude = cmplx(endpoint%bulk_gap, 0.0_rk, kind=rk) * &
      exp(imaginary_unit * cmplx( &
        endpoint%winding * angle + endpoint%bulk_phase_offset, &
        0.0_rk, kind=rk))
    bulk = amplitude * cmplx(endpoint%bulk_rotation, 0.0_rk, kind=rk)
  end function bulk_order_parameter_at_angle


  subroutine project_to_proper_rotation(matrix, rotation)
    real(rk), intent(in) :: matrix(3, 3)
    real(rk), intent(out) :: rotation(3, 3)

    real(rk) :: inverse(3, 3), next_rotation(3, 3)
    integer :: iteration
    logical :: succeeded

    rotation = matrix
    do iteration = 1, 32
      call inverse_3x3(rotation, inverse, succeeded)
      if (.not. succeeded) &
        error stop "bulk-reference matrix is singular"
      next_rotation = 0.5_rk * (rotation + transpose(inverse))
      if (maxval(abs(next_rotation - rotation)) <= &
          64.0_rk * epsilon(1.0_rk)) then
        rotation = next_rotation
        exit
      end if
      rotation = next_rotation
    end do
    if (determinant_3x3(rotation) <= 0.0_rk .or. &
        maxval(abs(matmul(transpose(rotation), rotation) - &
                   identity_rotation)) > 1.0e-10_rk) &
      error stop "bulk-reference matrix could not be projected onto SO(3)"
  end subroutine project_to_proper_rotation


  pure subroutine inverse_3x3(matrix, inverse, succeeded)
    real(rk), intent(in) :: matrix(3, 3)
    real(rk), intent(out) :: inverse(3, 3)
    logical, intent(out) :: succeeded

    real(rk) :: determinant, scale

    determinant = determinant_3x3(matrix)
    scale = max(1.0_rk, maxval(abs(matrix)))
    succeeded = abs(determinant) > 256.0_rk * epsilon(1.0_rk) * scale**3
    if (.not. succeeded) then
      inverse = 0.0_rk
      return
    end if
    inverse(1, 1) = matrix(2, 2) * matrix(3, 3) - &
                    matrix(2, 3) * matrix(3, 2)
    inverse(1, 2) = matrix(1, 3) * matrix(3, 2) - &
                    matrix(1, 2) * matrix(3, 3)
    inverse(1, 3) = matrix(1, 2) * matrix(2, 3) - &
                    matrix(1, 3) * matrix(2, 2)
    inverse(2, 1) = matrix(2, 3) * matrix(3, 1) - &
                    matrix(2, 1) * matrix(3, 3)
    inverse(2, 2) = matrix(1, 1) * matrix(3, 3) - &
                    matrix(1, 3) * matrix(3, 1)
    inverse(2, 3) = matrix(1, 3) * matrix(2, 1) - &
                    matrix(1, 1) * matrix(2, 3)
    inverse(3, 1) = matrix(2, 1) * matrix(3, 2) - &
                    matrix(2, 2) * matrix(3, 1)
    inverse(3, 2) = matrix(1, 2) * matrix(3, 1) - &
                    matrix(1, 1) * matrix(3, 2)
    inverse(3, 3) = matrix(1, 1) * matrix(2, 2) - &
                    matrix(1, 2) * matrix(2, 1)
    inverse = inverse / determinant
  end subroutine inverse_3x3


  pure real(rk) function determinant_3x3(matrix) result(determinant)
    real(rk), intent(in) :: matrix(3, 3)

    determinant = &
      matrix(1, 1) * (matrix(2, 2) * matrix(3, 3) - &
                      matrix(2, 3) * matrix(3, 2)) - &
      matrix(1, 2) * (matrix(2, 1) * matrix(3, 3) - &
                      matrix(2, 3) * matrix(3, 1)) + &
      matrix(1, 3) * (matrix(2, 1) * matrix(3, 2) - &
                      matrix(2, 2) * matrix(3, 1))
  end function determinant_3x3


  pure logical function free_vortex_endpoint_is_valid_for(self, mesh) &
      result(valid)
    class(free_vortex_endpoint_2d_t), intent(in) :: self
    type(cartesian_mesh_2d_t), intent(in) :: mesh

    real(rk) :: corner_radius, determinant, edge_radii(2)
    real(rk) :: orthogonality_error

    valid = mesh%is_valid()
    if (.not. valid .or. .not. self%enabled) return
    valid = self%bulk_gap > 0.0_rk .and. self%outer_radius > 0.0_rk .and. &
      self%maximum_step > 0.0_rk .and. &
      ieee_is_finite(self%bulk_gap) .and. &
      ieee_is_finite(self%winding) .and. &
      ieee_is_finite(self%bulk_phase_offset) .and. &
      all(ieee_is_finite(self%bulk_rotation)) .and. &
      self%center(1) > mesh%x_minimum .and. &
      self%center(1) < mesh%x_maximum .and. &
      self%center(2) > mesh%y_minimum .and. &
      self%center(2) < mesh%y_maximum
    if (.not. valid) return
    determinant = determinant_3x3(self%bulk_rotation)
    orthogonality_error = maxval(abs( &
      matmul(transpose(self%bulk_rotation), self%bulk_rotation) - &
      identity_rotation))
    valid = abs(determinant - 1.0_rk) <= 1.0e-10_rk .and. &
      orthogonality_error <= 1.0e-10_rk
    if (.not. valid) return
    if (self%use_axisymmetric_radial_reference) then
      valid = self%radial_profile%is_valid_for(self%radial_grid) .and. &
        self%matching_radius > 0.0_rk .and. &
        self%matching_radius <= self%radial_grid%outer_radius() .and. &
        self%outer_radius > self%matching_radius
      if (.not. valid) return
    end if
    if (self%use_inverse_radius_fit) then
      if (self%use_elliptical_fit_surfaces) then
        valid = .not. self%use_axisymmetric_radial_reference .and. &
          all(self%fit_inner_radii > 0.0_rk) .and. &
          all(self%fit_inner_radii < self%matching_radii)
        if (.not. valid) return
        edge_radii = [ &
          min(self%center(1) - mesh%x_minimum, &
              mesh%x_maximum - self%center(1)), &
          min(self%center(2) - mesh%y_minimum, &
              mesh%y_maximum - self%center(2))]
        valid = all(self%matching_radii < edge_radii)
      else
        valid = self%fit_inner_radius > 0.0_rk .and. &
          self%fit_inner_radius < self%matching_radius
        if (.not. valid) return
        if (.not. self%use_axisymmetric_radial_reference) &
          valid = self%matching_radius < &
            minimum_edge_distance(mesh, self%center)
      end if
      if (.not. valid) return
    end if
    corner_radius = maximum_corner_distance(mesh, self%center)
    valid = self%outer_radius > corner_radius * &
      (1.0_rk + 128.0_rk * epsilon(1.0_rk))
  end function free_vortex_endpoint_is_valid_for


  subroutine require_fit_radii(inner_radius, outer_radius, available_radius)
    real(rk), intent(in) :: inner_radius, outer_radius, available_radius

    if (inner_radius <= 0.0_rk .or. outer_radius <= inner_radius .or. &
        outer_radius >= available_radius) &
      error stop "inverse-radius fit radii are invalid"
  end subroutine require_fit_radii


  pure subroutine fit_radii_along_direction( &
      endpoint, direction, inner_radius, outer_radius)
    type(free_vortex_endpoint_2d_t), intent(in) :: endpoint
    real(rk), intent(in) :: direction(2)
    real(rk), intent(out) :: inner_radius, outer_radius

    if (endpoint%use_elliptical_fit_surfaces) then
      inner_radius = elliptical_radius_along_direction( &
        direction, endpoint%fit_inner_radii)
      outer_radius = elliptical_radius_along_direction( &
        direction, endpoint%matching_radii)
    else
      inner_radius = endpoint%fit_inner_radius
      outer_radius = endpoint%matching_radius
    end if
  end subroutine fit_radii_along_direction


  pure real(rk) function elliptical_radius_along_direction( &
      direction, radii) result(radius)
    real(rk), intent(in) :: direction(2), radii(2)

    radius = 1.0_rk / sqrt( &
      (direction(1) / radii(1))**2 + &
      (direction(2) / radii(2))**2)
  end function elliptical_radius_along_direction


  pure real(rk) function minimum_edge_distance(mesh, center) result(distance)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    real(rk), intent(in) :: center(2)

    distance = min(center(1) - mesh%x_minimum, &
                   mesh%x_maximum - center(1), &
                   center(2) - mesh%y_minimum, &
                   mesh%y_maximum - center(2))
  end function minimum_edge_distance


  pure real(rk) function radial_distance_to_rectangle(mesh, center, direction) &
      result(distance)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    real(rk), intent(in) :: center(2), direction(2)

    real(rk) :: candidate, tolerance

    distance = huge(1.0_rk)
    tolerance = 128.0_rk * epsilon(1.0_rk)
    if (direction(1) > tolerance) then
      candidate = (mesh%x_maximum - center(1)) / direction(1)
      distance = min(distance, candidate)
    else if (direction(1) < -tolerance) then
      candidate = (mesh%x_minimum - center(1)) / direction(1)
      distance = min(distance, candidate)
    end if
    if (direction(2) > tolerance) then
      candidate = (mesh%y_maximum - center(2)) / direction(2)
      distance = min(distance, candidate)
    else if (direction(2) < -tolerance) then
      candidate = (mesh%y_minimum - center(2)) / direction(2)
      distance = min(distance, candidate)
    end if
    if (distance >= 0.5_rk * huge(1.0_rk) .or. distance <= 0.0_rk) &
      error stop "could not intersect a radial line with the field mesh"
  end function radial_distance_to_rectangle


  pure real(rk) function maximum_corner_distance(mesh, center) result(distance)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    real(rk), intent(in) :: center(2)

    distance = max( &
      sqrt((mesh%x_minimum - center(1))**2 + &
           (mesh%y_minimum - center(2))**2), &
      sqrt((mesh%x_maximum - center(1))**2 + &
           (mesh%y_minimum - center(2))**2), &
      sqrt((mesh%x_minimum - center(1))**2 + &
           (mesh%y_maximum - center(2))**2), &
      sqrt((mesh%x_maximum - center(1))**2 + &
           (mesh%y_maximum - center(2))**2))
  end function maximum_corner_distance

end module free_vortex_asymptotic_2d
