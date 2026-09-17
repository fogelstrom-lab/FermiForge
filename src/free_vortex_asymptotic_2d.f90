module free_vortex_asymptotic_2d
  use he3_kinds, only : rk
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

  type, public :: free_vortex_endpoint_2d_t
    logical :: enabled = .false.
    logical :: use_axisymmetric_radial_reference = .false.
    real(rk) :: center(2) = 0.0_rk
    real(rk) :: bulk_gap = 0.0_rk
    real(rk) :: winding = 1.0_rk
    real(rk) :: outer_radius = 0.0_rk
    real(rk) :: maximum_step = 0.0_rk
    real(rk) :: matching_radius = 0.0_rk
    type(radial_mesh_t) :: radial_grid
    type(axial_radial_profile_t) :: radial_profile
  contains
    procedure :: is_valid_for => free_vortex_endpoint_is_valid_for
  end type free_vortex_endpoint_2d_t

  public :: extend_with_free_vortex_asymptotic_2d
  public :: sample_free_vortex_asymptotic_state_2d
  public :: set_axisymmetric_radial_reference

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
    endpoint%matching_radius = matching_radius
    endpoint%radial_grid = radial_grid
    endpoint%radial_profile = radial_profile
  end subroutine set_axisymmetric_radial_reference

  subroutine sample_free_vortex_asymptotic_state_2d( &
      mesh, state, position, endpoint, order_parameter, current_mean_field)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    real(rk), intent(in) :: position(2)
    type(free_vortex_endpoint_2d_t), intent(in) :: endpoint
    complex(rk), intent(out) :: order_parameter(3, 3)
    real(rk), intent(out) :: current_mean_field(3)

    complex(rk) :: boundary_order_parameter(3, 3), bulk_value
    real(rk) :: boundary_current(3), boundary_position(2)
    real(rk) :: decay, direction(2), displacement(2), matching_radius
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

    phase_angle = atan2(direction(2), direction(1))
    bulk_value = cmplx(endpoint%bulk_gap, 0.0_rk, kind=rk) * &
      exp(imaginary_unit * &
          cmplx(endpoint%winding * phase_angle, 0.0_rk, kind=rk))
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
        if (spin == orbital) order_parameter(spin, orbital) = &
          order_parameter(spin, orbital) + &
          cmplx(1.0_rk - decay, 0.0_rk, kind=rk) * bulk_value
      end do
    end do
    current_mean_field = (matching_radius / radius) * boundary_current
  end subroutine sample_free_vortex_asymptotic_state_2d


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


  pure logical function free_vortex_endpoint_is_valid_for(self, mesh) &
      result(valid)
    class(free_vortex_endpoint_2d_t), intent(in) :: self
    type(cartesian_mesh_2d_t), intent(in) :: mesh

    real(rk) :: corner_radius

    valid = mesh%is_valid()
    if (.not. valid .or. .not. self%enabled) return
    valid = self%bulk_gap > 0.0_rk .and. self%outer_radius > 0.0_rk .and. &
      self%maximum_step > 0.0_rk .and. &
      self%center(1) > mesh%x_minimum .and. &
      self%center(1) < mesh%x_maximum .and. &
      self%center(2) > mesh%y_minimum .and. &
      self%center(2) < mesh%y_maximum
    if (.not. valid) return
    if (self%use_axisymmetric_radial_reference) then
      valid = self%radial_profile%is_valid_for(self%radial_grid) .and. &
        self%matching_radius > 0.0_rk .and. &
        self%matching_radius <= self%radial_grid%outer_radius() .and. &
        self%outer_radius > self%matching_radius
      if (.not. valid) return
    end if
    corner_radius = maximum_corner_distance(mesh, self%center)
    valid = self%outer_radius > corner_radius * &
      (1.0_rk + 128.0_rk * epsilon(1.0_rk))
  end function free_vortex_endpoint_is_valid_for


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
