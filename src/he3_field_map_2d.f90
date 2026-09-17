module he3_field_map_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t
  use spinful_state_2d, only : spinful_state_2d_t, allocate_spinful_state_2d
  use he3_quadrature, only : angular_quadrature_3d_t, ozaki_quadrature_t
  use free_vortex_asymptotic_2d, only : free_vortex_endpoint_2d_t
  use he3_serial_point_map, only : he3_point_map_diagnostics_t, &
                                   evaluate_serial_he3_point_map
  implicit none
  private

  type, public :: he3_field_map_diagnostics_t
    integer :: requested_point_count = 0
    integer :: completed_point_count = 0
    integer :: failed_point_count = 0
    integer :: direction_count = 0
    integer :: pole_count = 0
    integer :: attempted_contributions = 0
    integer :: accepted_contributions = 0
    integer :: failed_propagators = 0
    real(rk) :: maximum_normalization_error = 0.0_rk
    real(rk) :: sampling_seconds = 0.0_rk
    real(rk) :: propagation_seconds = 0.0_rk
    real(rk) :: accumulation_seconds = 0.0_rk
  end type he3_field_map_diagnostics_t

  public :: evaluate_he3_field_map_subset
  public :: evaluate_serial_he3_field_map

contains

  subroutine evaluate_serial_he3_field_map( &
      mesh, state, angular_quadrature, energy_quadrature, feedback_scale, &
      trajectory_maximum_step, minimum_substep_count, &
      boundary_relaxation_distance, mapped_state, diagnostics, succeeded, &
      free_vortex_endpoint, active_point_mask)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    type(angular_quadrature_3d_t), intent(in) :: angular_quadrature
    type(ozaki_quadrature_t), intent(in) :: energy_quadrature
    real(rk), intent(in) :: feedback_scale, trajectory_maximum_step
    integer, intent(in) :: minimum_substep_count
    real(rk), intent(in) :: boundary_relaxation_distance
    type(spinful_state_2d_t), intent(out) :: mapped_state
    type(he3_field_map_diagnostics_t), intent(out) :: diagnostics
    logical, intent(out) :: succeeded
    type(free_vortex_endpoint_2d_t), intent(in), optional :: &
      free_vortex_endpoint
    logical, intent(in), optional :: active_point_mask(:)

    integer, allocatable :: point_indices(:)
    logical, allocatable :: point_completed(:)
    logical, allocatable :: active(:)
    integer :: list_index, point

    allocate(active(mesh%point_count()), source=.true.)
    if (present(active_point_mask)) then
      if (size(active_point_mask) /= mesh%point_count()) &
        error stop "serial field-map active mask has the wrong size"
      active = active_point_mask
      if (.not. any(active)) &
        error stop "serial field-map active mask contains no points"
    end if
    allocate(point_indices(count(active)))
    list_index = 0
    do point = 1, mesh%point_count()
      if (.not. active(point)) cycle
      list_index = list_index + 1
      point_indices(list_index) = point
    end do
    call evaluate_he3_field_map_subset( &
      mesh, state, point_indices, angular_quadrature, energy_quadrature, &
      feedback_scale, trajectory_maximum_step, minimum_substep_count, &
      boundary_relaxation_distance, mapped_state, point_completed, &
      diagnostics, succeeded, free_vortex_endpoint)
    do point = 1, mesh%point_count()
      if (active(point)) cycle
      mapped_state%order_parameter(:, :, point) = &
        state%order_parameter(:, :, point)
      mapped_state%current_mean_field(:, point) = &
        state%current_mean_field(:, point)
    end do
    succeeded = succeeded .and. all(point_completed .eqv. active)
  end subroutine evaluate_serial_he3_field_map


  subroutine evaluate_he3_field_map_subset( &
      mesh, state, point_indices, angular_quadrature, energy_quadrature, &
      feedback_scale, trajectory_maximum_step, minimum_substep_count, &
      boundary_relaxation_distance, mapped_state, point_completed, &
      diagnostics, succeeded, free_vortex_endpoint)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    integer, intent(in) :: point_indices(:)
    type(angular_quadrature_3d_t), intent(in) :: angular_quadrature
    type(ozaki_quadrature_t), intent(in) :: energy_quadrature
    real(rk), intent(in) :: feedback_scale, trajectory_maximum_step
    integer, intent(in) :: minimum_substep_count
    real(rk), intent(in) :: boundary_relaxation_distance
    type(spinful_state_2d_t), intent(out) :: mapped_state
    logical, allocatable, intent(out) :: point_completed(:)
    type(he3_field_map_diagnostics_t), intent(out) :: diagnostics
    logical, intent(out) :: succeeded
    type(free_vortex_endpoint_2d_t), intent(in), optional :: &
      free_vortex_endpoint

    type(he3_point_map_diagnostics_t) :: point_diagnostics
    complex(rk) :: mapped_gap(3, 3)
    real(rk) :: mapped_current(3), origin(2)
    logical, allocatable :: point_claimed(:)
    integer :: list_index, point
    logical :: point_succeeded

    if (.not. state%is_valid_for(mesh)) &
      error stop "field map state does not match its Cartesian mesh"
    allocate(point_claimed(mesh%point_count()), source=.false.)
    do list_index = 1, size(point_indices)
      point = point_indices(list_index)
      if (point < 1 .or. point > mesh%point_count()) &
        error stop "field-map target point is out of range"
      if (point_claimed(point)) &
        error stop "field-map target point was requested more than once"
      point_claimed(point) = .true.
    end do

    call allocate_spinful_state_2d(mesh, mapped_state)
    allocate(point_completed(mesh%point_count()), source=.false.)
    diagnostics = he3_field_map_diagnostics_t()
    diagnostics%requested_point_count = size(point_indices)
    diagnostics%direction_count = angular_quadrature%direction_count()
    diagnostics%pole_count = energy_quadrature%pole_count()

    do list_index = 1, size(point_indices)
      point = point_indices(list_index)
      origin = mesh%point_coordinate(point)
      call evaluate_serial_he3_point_map( &
        mesh, state, origin, angular_quadrature, energy_quadrature, &
        feedback_scale, trajectory_maximum_step, minimum_substep_count, &
        boundary_relaxation_distance, mapped_gap, mapped_current, &
        point_diagnostics, point_succeeded, free_vortex_endpoint)
      call add_point_diagnostics(diagnostics, point_diagnostics)
      if (point_succeeded) then
        mapped_state%order_parameter(:, :, point) = mapped_gap
        mapped_state%current_mean_field(:, point) = mapped_current
        point_completed(point) = .true.
        diagnostics%completed_point_count = &
          diagnostics%completed_point_count + 1
      else
        diagnostics%failed_point_count = diagnostics%failed_point_count + 1
      end if
    end do
    succeeded = diagnostics%failed_point_count == 0 .and. &
      diagnostics%completed_point_count == diagnostics%requested_point_count
  end subroutine evaluate_he3_field_map_subset


  subroutine add_point_diagnostics(field, point)
    type(he3_field_map_diagnostics_t), intent(inout) :: field
    type(he3_point_map_diagnostics_t), intent(in) :: point

    field%attempted_contributions = field%attempted_contributions + &
      point%attempted_contributions
    field%accepted_contributions = field%accepted_contributions + &
      point%accepted_contributions
    field%failed_propagators = field%failed_propagators + &
      point%failed_propagators
    field%maximum_normalization_error = max( &
      field%maximum_normalization_error, point%maximum_normalization_error)
    field%sampling_seconds = field%sampling_seconds + point%sampling_seconds
    field%propagation_seconds = field%propagation_seconds + &
      point%propagation_seconds
    field%accumulation_seconds = field%accumulation_seconds + &
      point%accumulation_seconds
  end subroutine add_point_diagnostics

end module he3_field_map_2d
