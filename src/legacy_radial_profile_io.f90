module legacy_radial_profile_io
  use he3_kinds, only : rk
  use radial_mesh, only : radial_mesh_t
  use axial_radial_embedding, only : axial_radial_profile_t, &
                                     set_axial_profile_from_axis_cartesian
  implicit none
  private

  public :: read_new_src_radial_profile

contains

  subroutine read_new_src_radial_profile( &
      order_parameter_path, current_path, mesh, profile, tangent_scale, &
      radial_grid_kind, detected_grid_kind)
    character(len=*), intent(in) :: order_parameter_path, current_path
    type(radial_mesh_t), intent(out) :: mesh
    type(axial_radial_profile_t), intent(out) :: profile
    real(rk), intent(in), optional :: tangent_scale
    character(len=*), intent(in), optional :: radial_grid_kind
    character(len=*), intent(out), optional :: detected_grid_kind

    complex(rk), allocatable :: axis_cartesian(:, :, :)
    real(rk), allocatable :: azimuthal_mean_field(:)
    real(rk) :: current_record(5), order_record(19)
    real(rk), allocatable :: raw_radius(:), tangent_radius(:), uniform_radius(:)
    real(rk) :: outer_radius, radial_scale, radius_tolerance, tangent_step
    real(rk) :: uniform_step
    character(len=16) :: requested_kind, selected_kind
    integer :: current_unit, input_status, orbital, order_unit, point
    integer :: point_count, spin
    logical :: matches_tangent, matches_uniform

    point_count = count_records(order_parameter_path)
    if (point_count < 3) &
      error stop "new_src radial order-parameter file is too short"
    if (count_records(current_path) /= point_count) &
      error stop "new_src radial state files have different lengths"

    allocate(mesh%r(0:point_count - 1))
    allocate(axis_cartesian(3, 3, 0:point_count - 1), &
             azimuthal_mean_field(0:point_count - 1), &
             raw_radius(0:point_count - 1), &
             tangent_radius(0:point_count - 1), &
             uniform_radius(0:point_count - 1))

    open(newunit=order_unit, file=trim(order_parameter_path), status="old", &
         action="read", iostat=input_status)
    if (input_status /= 0) &
      error stop "could not open new_src radial order-parameter file"
    open(newunit=current_unit, file=trim(current_path), status="old", &
         action="read", iostat=input_status)
    if (input_status /= 0) &
      error stop "could not open new_src radial current-field file"

    do point = 0, point_count - 1
      read(order_unit, *, iostat=input_status) order_record
      if (input_status /= 0) &
        error stop "could not read new_src radial order-parameter record"
      read(current_unit, *, iostat=input_status) current_record
      if (input_status /= 0) &
        error stop "could not read new_src radial current-field record"

      raw_radius(point) = order_record(1)
      radius_tolerance = 256.0_rk * epsilon(1.0_rk) * &
        max(1.0_rk, abs(raw_radius(point)))
      if (abs(current_record(1) - raw_radius(point)) > radius_tolerance) &
        error stop "new_src radial files use different coordinates"
      do spin = 1, 3
        do orbital = 1, 3
          axis_cartesian(spin, orbital, point) = cmplx( &
            order_record(2 * ((spin - 1) * 3 + orbital)), &
            order_record(2 * ((spin - 1) * 3 + orbital) + 1), kind=rk)
        end do
      end do
      ! On the positive x axis the source's azimuthal field is vy.  Its
      ! radial and axial columns are not part of the axial transport model.
      azimuthal_mean_field(point) = current_record(4)
      if (max(abs(current_record(3)), abs(current_record(5))) > &
          1.0e-9_rk * max(1.0_rk, abs(current_record(4)))) &
        error stop "new_src radial state is not purely azimuthal"
    end do
    close(order_unit)
    close(current_unit)

    ! The current new_src solver writes a tangent mesh, while the archived
    ! 50-point n-core files were produced by the earlier uniform-grid solver.
    ! Both write only three coordinate decimals, so identify the declared
    ! family and reconstruct its full-precision coordinates.
    radial_scale = 10.0_rk
    if (present(tangent_scale)) radial_scale = tangent_scale
    if (radial_scale <= 0.0_rk) &
      error stop "new_src tangent radial scale must be positive"
    outer_radius = raw_radius(point_count - 1)
    uniform_step = outer_radius / real(point_count - 1, rk)
    tangent_step = atan(outer_radius / radial_scale) / &
      real(point_count - 1, rk)
    do point = 0, point_count - 1
      uniform_radius(point) = uniform_step * real(point, rk)
      tangent_radius(point) = radial_scale * &
        tan(tangent_step * real(point, rk))
    end do
    matches_uniform = maxval(abs(raw_radius - uniform_radius)) <= 6.0e-4_rk
    matches_tangent = maxval(abs(raw_radius - tangent_radius)) <= 6.0e-4_rk

    requested_kind = "auto"
    if (present(radial_grid_kind)) &
      requested_kind = trim(adjustl(radial_grid_kind))
    select case (trim(requested_kind))
    case ("auto")
      if (matches_tangent .and. .not. matches_uniform) then
        selected_kind = "tangent"
      else if (matches_uniform .and. .not. matches_tangent) then
        selected_kind = "uniform"
      else if (matches_tangent .and. matches_uniform) then
        error stop "radial grid is ambiguous; select 'uniform' or 'tangent'"
      else
        error stop "radial coordinates match neither supported grid family"
      end if
    case ("tangent")
      if (.not. matches_tangent) &
        error stop "radial coordinates do not match the tangent grid"
      selected_kind = "tangent"
    case ("uniform")
      if (.not. matches_uniform) &
        error stop "radial coordinates do not match the uniform grid"
      selected_kind = "uniform"
    case default
      error stop "radial_grid_kind must be 'auto', 'uniform', or 'tangent'"
    end select

    select case (trim(selected_kind))
    case ("tangent")
      mesh%r = tangent_radius
    case ("uniform")
      mesh%r = uniform_radius
    end select
    if (present(detected_grid_kind)) detected_grid_kind = trim(selected_kind)

    if (.not. mesh%is_valid()) &
      error stop "new_src radial file does not define a valid mesh"
    call set_axial_profile_from_axis_cartesian( &
      mesh, axis_cartesian, azimuthal_mean_field, profile)
  end subroutine read_new_src_radial_profile


  integer function count_records(path) result(record_count)
    character(len=*), intent(in) :: path

    character(len=4096) :: line
    integer :: input_status, unit

    record_count = 0
    open(newunit=unit, file=trim(path), status="old", action="read", &
         iostat=input_status)
    if (input_status /= 0) error stop "could not open new_src radial data file"
    do
      read(unit, '(a)', iostat=input_status) line
      if (input_status < 0) exit
      if (input_status > 0) error stop "could not scan new_src radial data file"
      if (len_trim(line) > 0) record_count = record_count + 1
    end do
    close(unit)
  end function count_records

end module legacy_radial_profile_io
