module spinful_sparse_state_io_2d
  use, intrinsic :: iso_fortran_env, only : iostat_end
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t
  use spinful_state_2d, only : spinful_state_2d_t
  implicit none
  private

  integer, parameter :: record_value_count = 23

  public :: write_sparse_spinful_state_2d
  public :: read_sparse_spinful_state_2d

contains

  subroutine write_sparse_spinful_state_2d( &
      path, mesh, state, point_mask, label)
    character(len=*), intent(in) :: path
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    logical, intent(in) :: point_mask(:)
    character(len=*), intent(in), optional :: label

    character(len=2048) :: message
    integer :: point, status, unit
    real(rk) :: values(record_value_count)

    call require_valid_arguments(mesh, state, point_mask)
    open(newunit=unit, file=trim(path), status="replace", action="write", &
         iostat=status, iomsg=message)
    if (status /= 0) &
      error stop "cannot open sparse 2D state output: " // trim(message)
    write(unit, '(a)') "# FermiForge sparse spinful state version 1"
    if (present(label)) write(unit, '(a,a)') "# label=", trim(label)
    write(unit, '(a,i0)') "# point_count=", count(point_mask)
    write(unit, '(a)') &
      "# x y Re/Im[A(1,1)..A(3,3)] current_mean_field(1:3)"
    do point = 1, mesh%point_count()
      if (.not. point_mask(point)) cycle
      call pack_record(mesh, state, point, values)
      write(unit, '(23(es24.16e3,1x))') values
    end do
    close(unit)
  end subroutine write_sparse_spinful_state_2d


  subroutine read_sparse_spinful_state_2d( &
      path, mesh, state, applied_point_count)
    character(len=*), intent(in) :: path
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(inout) :: state
    integer, intent(out) :: applied_point_count

    character(len=2048) :: line, message
    integer :: point, status, unit, x_node, y_node
    logical, allocatable :: seen(:)
    real(rk) :: values(record_value_count)

    if (.not. state%is_valid_for(mesh)) &
      error stop "sparse 2D restart state does not match its mesh"
    allocate(seen(mesh%point_count()), source=.false.)
    applied_point_count = 0
    open(newunit=unit, file=trim(path), status="old", action="read", &
         iostat=status, iomsg=message)
    if (status /= 0) &
      error stop "cannot open sparse 2D restart: " // trim(message)
    do
      read(unit, '(a)', iostat=status, iomsg=message) line
      if (status == iostat_end) exit
      if (status /= 0) &
        error stop "cannot read sparse 2D restart: " // trim(message)
      if (len_trim(line) == 0 .or. line(1:1) == "#") cycle
      read(line, *, iostat=status, iomsg=message) values
      if (status /= 0) &
        error stop "invalid sparse 2D restart record: " // trim(message)
      call locate_mesh_node(mesh, values(1), values(2), x_node, y_node)
      point = mesh%point_index(x_node, y_node)
      if (seen(point)) error stop "duplicate sparse 2D restart coordinate"
      seen(point) = .true.
      call unpack_record(values, point, state)
      applied_point_count = applied_point_count + 1
    end do
    close(unit)
    if (applied_point_count < 1) &
      error stop "sparse 2D restart contains no field records"
  end subroutine read_sparse_spinful_state_2d


  subroutine pack_record(mesh, state, point, values)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    integer, intent(in) :: point
    real(rk), intent(out) :: values(record_value_count)

    integer :: component, orbital, position, spin

    values = 0.0_rk
    values(1:2) = mesh%point_coordinate(point)
    position = 2
    do spin = 1, 3
      do orbital = 1, 3
        position = position + 1
        values(position) = &
          real(state%order_parameter(spin, orbital, point), rk)
        position = position + 1
        values(position) = aimag(state%order_parameter(spin, orbital, point))
      end do
    end do
    do component = 1, 3
      position = position + 1
      values(position) = state%current_mean_field(component, point)
    end do
  end subroutine pack_record


  subroutine unpack_record(values, point, state)
    real(rk), intent(in) :: values(record_value_count)
    integer, intent(in) :: point
    type(spinful_state_2d_t), intent(inout) :: state

    integer :: component, orbital, position, spin

    position = 2
    do spin = 1, 3
      do orbital = 1, 3
        state%order_parameter(spin, orbital, point) = &
          cmplx(values(position + 1), values(position + 2), rk)
        position = position + 2
      end do
    end do
    do component = 1, 3
      position = position + 1
      state%current_mean_field(component, point) = values(position)
    end do
  end subroutine unpack_record


  subroutine locate_mesh_node(mesh, x, y, x_node, y_node)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    real(rk), intent(in) :: x, y
    integer, intent(out) :: x_node, y_node

    real(rk) :: tolerance

    x_node = nearest_axis_node(mesh, x, .true.)
    y_node = nearest_axis_node(mesh, y, .false.)
    tolerance = 512.0_rk * epsilon(1.0_rk) * max( &
      1.0_rk, abs(x), abs(y), abs(mesh%x_minimum), abs(mesh%x_maximum), &
      abs(mesh%y_minimum), abs(mesh%y_maximum))
    if (abs(mesh%x_coordinate(x_node) - x) > tolerance .or. &
        abs(mesh%y_coordinate(y_node) - y) > tolerance) &
      error stop "sparse 2D restart coordinate is not a mesh node"
  end subroutine locate_mesh_node


  integer function nearest_axis_node(mesh, coordinate, use_x) result(nearest)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    real(rk), intent(in) :: coordinate
    logical, intent(in) :: use_x

    integer :: node, node_count
    real(rk) :: difference, minimum_difference

    if (use_x) then
      node_count = mesh%x_cell_count()
    else
      node_count = mesh%y_cell_count()
    end if
    nearest = 0
    minimum_difference = huge(1.0_rk)
    do node = 0, node_count
      if (use_x) then
        difference = abs(mesh%x_coordinate(node) - coordinate)
      else
        difference = abs(mesh%y_coordinate(node) - coordinate)
      end if
      if (difference < minimum_difference) then
        nearest = node
        minimum_difference = difference
      end if
    end do
  end function nearest_axis_node


  subroutine require_valid_arguments(mesh, state, point_mask)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    logical, intent(in) :: point_mask(:)

    if (.not. state%is_valid_for(mesh)) &
      error stop "sparse 2D output state does not match its mesh"
    if (size(point_mask) /= mesh%point_count() .or. .not. any(point_mask)) &
      error stop "sparse 2D output mask is invalid"
  end subroutine require_valid_arguments

end module spinful_sparse_state_io_2d
