module legacy_split_field_io_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, &
                                make_rectilinear_cartesian_mesh
  use spinful_state_2d, only : spinful_state_2d_t, &
                               allocate_spinful_state_2d
  implicit none
  private

  type, public :: legacy_split_field_report_2d_t
    integer :: number_of_points = 0
    integer :: number_of_x_points = 0
    integer :: number_of_y_points = 0
    real(rk) :: maximum_gap_norm_difference = 0.0_rk
    real(rk) :: maximum_in_plane_magnitude_difference = 0.0_rk
    real(rk) :: minimum_gap_norm = huge(1.0_rk)
    real(rk) :: minimum_gap_x = 0.0_rk
    real(rk) :: minimum_gap_y = 0.0_rk
  end type legacy_split_field_report_2d_t

  public :: read_legacy_split_field_directory_2d
  public :: read_legacy_split_field_map_2d

contains

  subroutine read_legacy_split_field_directory_2d( &
      directory, mesh, state, report)
    character(len=*), intent(in) :: directory
    type(cartesian_mesh_2d_t), intent(out) :: mesh
    type(spinful_state_2d_t), intent(out) :: state
    type(legacy_split_field_report_2d_t), intent(out) :: report

    call read_legacy_split_field_map_2d( &
      path_in_directory(directory, "op_x"), &
      path_in_directory(directory, "op_y"), &
      path_in_directory(directory, "op_z"), &
      path_in_directory(directory, "curr"), mesh, state, report)
  end subroutine read_legacy_split_field_directory_2d


  subroutine read_legacy_split_field_map_2d( &
      op_x_path, op_y_path, op_z_path, current_path, mesh, state, report)
    character(len=*), intent(in) :: op_x_path, op_y_path, op_z_path
    character(len=*), intent(in) :: current_path
    type(cartesian_mesh_2d_t), intent(out) :: mesh
    type(spinful_state_2d_t), intent(out) :: state
    type(legacy_split_field_report_2d_t), intent(out) :: report

    complex(rk), allocatable :: order_row(:, :)
    real(rk), allocatable :: coordinates(:, :), other_coordinates(:, :)
    real(rk), allocatable :: x_coordinates(:), y_coordinates(:)

    call read_order_parameter_row(op_x_path, coordinates, order_row)
    call infer_rectilinear_axes( &
      coordinates, x_coordinates, y_coordinates)
    call make_rectilinear_cartesian_mesh( &
      x_coordinates, y_coordinates, mesh)
    call allocate_spinful_state_2d(mesh, state)
    state%order_parameter(1, :, :) = order_row

    call read_order_parameter_row(op_y_path, other_coordinates, order_row)
    call require_matching_coordinates( &
      coordinates, other_coordinates, op_y_path)
    state%order_parameter(2, :, :) = order_row

    call read_order_parameter_row(op_z_path, other_coordinates, order_row)
    call require_matching_coordinates( &
      coordinates, other_coordinates, op_z_path)
    state%order_parameter(3, :, :) = order_row

    report = legacy_split_field_report_2d_t()
    report%number_of_points = mesh%point_count()
    report%number_of_x_points = mesh%x_point_count()
    report%number_of_y_points = mesh%y_point_count()
    call read_current_field( &
      current_path, coordinates, state, report)
  end subroutine read_legacy_split_field_map_2d


  subroutine read_order_parameter_row(path, coordinates, order_row)
    character(len=*), intent(in) :: path
    real(rk), allocatable, intent(out) :: coordinates(:, :)
    complex(rk), allocatable, intent(out) :: order_row(:, :)

    character(len=4096) :: line, message
    integer :: component, ios, point, record_count, unit
    real(rk) :: record(8)

    record_count = count_data_records(path)
    if (record_count < 4) &
      error stop "legacy order-parameter file has too few records: " // &
                 trim(path)
    allocate(coordinates(2, record_count))
    allocate(order_row(3, record_count))

    open(newunit=unit, file=trim(path), status="old", action="read", &
         iostat=ios, iomsg=message)
    if (ios /= 0) error stop "cannot open legacy field file: " // trim(message)
    point = 0
    do
      read(unit, '(a)', iostat=ios, iomsg=message) line
      if (ios < 0) exit
      if (ios > 0) error stop "cannot read legacy field file: " // trim(message)
      if (.not. is_data_record(line)) cycle
      point = point + 1
      read(line, *, iostat=ios, iomsg=message) record
      if (ios /= 0) error stop "invalid legacy order-parameter record: " // &
                               trim(message)
      coordinates(:, point) = record(1:2)
      do component = 1, 3
        order_row(component, point) = cmplx( &
          record(2 * component + 1), record(2 * component + 2), kind=rk)
      end do
    end do
    close(unit)
    if (point /= record_count) &
      error stop "legacy order-parameter record count changed while reading"
  end subroutine read_order_parameter_row


  subroutine infer_rectilinear_axes( &
      coordinates, x_coordinates, y_coordinates)
    real(rk), intent(in) :: coordinates(:, :)
    real(rk), allocatable, intent(out) :: x_coordinates(:), y_coordinates(:)

    integer :: number_of_points, number_of_x_points, number_of_y_points
    integer :: point, x_node, y_node

    if (size(coordinates, 1) /= 2) &
      error stop "legacy coordinate array must have two rows"
    number_of_points = size(coordinates, 2)
    number_of_x_points = 1
    do point = 2, number_of_points
      if (.not. same_coordinate( &
          coordinates(2, point), coordinates(2, 1))) exit
      number_of_x_points = number_of_x_points + 1
    end do
    if (number_of_x_points < 2 .or. &
        modulo(number_of_points, number_of_x_points) /= 0) &
      error stop "legacy field records do not form complete rectilinear rows"
    number_of_y_points = number_of_points / number_of_x_points
    if (number_of_y_points < 2) &
      error stop "legacy field records contain fewer than two y rows"

    allocate(x_coordinates(number_of_x_points))
    allocate(y_coordinates(number_of_y_points))
    x_coordinates = coordinates(1, 1:number_of_x_points)
    do y_node = 1, number_of_y_points
      point = 1 + (y_node - 1) * number_of_x_points
      y_coordinates(y_node) = coordinates(2, point)
    end do

    do y_node = 1, number_of_y_points
      do x_node = 1, number_of_x_points
        point = x_node + (y_node - 1) * number_of_x_points
        if (.not. same_coordinate( &
              coordinates(1, point), x_coordinates(x_node)) .or. &
            .not. same_coordinate( &
              coordinates(2, point), y_coordinates(y_node))) &
          error stop "legacy field coordinates are not an x-fast rectilinear grid"
      end do
    end do
  end subroutine infer_rectilinear_axes


  subroutine require_matching_coordinates(reference, actual, path)
    real(rk), intent(in) :: reference(:, :), actual(:, :)
    character(len=*), intent(in) :: path

    integer :: point

    if (any(shape(reference) /= shape(actual))) &
      error stop "legacy field point count differs in: " // trim(path)
    do point = 1, size(reference, 2)
      if (.not. same_coordinate(reference(1, point), actual(1, point)) .or. &
          .not. same_coordinate(reference(2, point), actual(2, point))) &
        error stop "legacy field coordinates differ in: " // trim(path)
    end do
  end subroutine require_matching_coordinates


  subroutine read_current_field(path, coordinates, state, report)
    character(len=*), intent(in) :: path
    real(rk), intent(in) :: coordinates(:, :)
    type(spinful_state_2d_t), intent(inout) :: state
    type(legacy_split_field_report_2d_t), intent(inout) :: report

    character(len=4096) :: line, message
    integer :: ios, point, record_count, unit
    real(rk) :: gap_norm, in_plane_magnitude, record(7)

    record_count = count_data_records(path)
    if (record_count /= size(coordinates, 2)) &
      error stop "legacy current file has the wrong number of records"
    open(newunit=unit, file=trim(path), status="old", action="read", &
         iostat=ios, iomsg=message)
    if (ios /= 0) error stop "cannot open legacy current file: " // trim(message)

    point = 0
    do
      read(unit, '(a)', iostat=ios, iomsg=message) line
      if (ios < 0) exit
      if (ios > 0) error stop "cannot read legacy current file: " // trim(message)
      if (.not. is_data_record(line)) cycle
      point = point + 1
      read(line, *, iostat=ios, iomsg=message) record
      if (ios /= 0) error stop "invalid legacy current record: " // trim(message)
      if (.not. same_coordinate(record(1), coordinates(1, point)) .or. &
          .not. same_coordinate(record(2), coordinates(2, point))) &
        error stop "legacy current coordinates differ from the order parameter"

      state%current_mean_field(:, point) = record(5:7)
      gap_norm = sqrt(sum(abs(state%order_parameter(:, :, point))**2) / &
                      3.0_rk)
      in_plane_magnitude = hypot(record(5), record(6))
      report%maximum_gap_norm_difference = max( &
        report%maximum_gap_norm_difference, abs(gap_norm - record(3)))
      report%maximum_in_plane_magnitude_difference = max( &
        report%maximum_in_plane_magnitude_difference, &
        abs(in_plane_magnitude - record(4)))
      if (gap_norm < report%minimum_gap_norm) then
        report%minimum_gap_norm = gap_norm
        report%minimum_gap_x = record(1)
        report%minimum_gap_y = record(2)
      end if
    end do
    close(unit)
  end subroutine read_current_field


  integer function count_data_records(path) result(record_count)
    character(len=*), intent(in) :: path

    character(len=4096) :: line, message
    integer :: ios, unit

    record_count = 0
    open(newunit=unit, file=trim(path), status="old", action="read", &
         iostat=ios, iomsg=message)
    if (ios /= 0) error stop "cannot open legacy field file: " // trim(message)
    do
      read(unit, '(a)', iostat=ios, iomsg=message) line
      if (ios < 0) exit
      if (ios > 0) error stop "cannot scan legacy field file: " // trim(message)
      if (is_data_record(line)) record_count = record_count + 1
    end do
    close(unit)
  end function count_data_records


  pure logical function is_data_record(line) result(is_data)
    character(len=*), intent(in) :: line

    character(len=len(line)) :: adjusted

    adjusted = adjustl(line)
    is_data = len_trim(adjusted) > 0
    if (is_data) is_data = adjusted(1:1) /= "#"
  end function is_data_record


  pure logical function same_coordinate(first, second) result(same)
    real(rk), intent(in) :: first, second

    real(rk) :: tolerance

    tolerance = 256.0_rk * epsilon(1.0_rk) * &
      max(1.0_rk, abs(first), abs(second))
    same = abs(first - second) <= tolerance
  end function same_coordinate


  pure function path_in_directory(directory, leaf) result(path)
    character(len=*), intent(in) :: directory, leaf
    character(len=:), allocatable :: path

    integer :: length

    length = len_trim(directory)
    if (length > 0) then
      if (directory(length:length) == "/") then
        path = trim(directory) // trim(leaf)
      else
        path = trim(directory) // "/" // trim(leaf)
      end if
    else
      path = trim(leaf)
    end if
  end function path_in_directory

end module legacy_split_field_io_2d
