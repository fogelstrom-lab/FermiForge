program test_legacy_split_field_io_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t
  use spinful_state_2d, only : spinful_state_2d_t
  use legacy_split_field_io_2d, only : &
    legacy_split_field_report_2d_t, read_legacy_split_field_map_2d
  implicit none

  character(len=*), parameter :: op_x_path = "legacy_split_op_x.tmp"
  character(len=*), parameter :: op_y_path = "legacy_split_op_y.tmp"
  character(len=*), parameter :: op_z_path = "legacy_split_op_z.tmp"
  character(len=*), parameter :: current_path = "legacy_split_curr.tmp"
  real(rk), parameter :: x_coordinates(3) = [-2.0_rk, 0.0_rk, 3.0_rk]
  real(rk), parameter :: y_coordinates(2) = [-1.0_rk, 2.0_rk]

  complex(rk) :: expected_order_parameter(3, 3, 6)
  real(rk) :: expected_mean_field(3, 6)
  integer :: component, orbital, point, spin, x_node, y_node
  type(cartesian_mesh_2d_t) :: mesh
  type(spinful_state_2d_t) :: state
  type(legacy_split_field_report_2d_t) :: report

  point = 0
  do y_node = 1, size(y_coordinates)
    do x_node = 1, size(x_coordinates)
      point = point + 1
      do spin = 1, 3
        do orbital = 1, 3
          expected_order_parameter(spin, orbital, point) = cmplx( &
            0.01_rk * real(100 * spin + 10 * orbital + point, rk), &
            -0.02_rk * real(10 * spin + orbital + point, rk), kind=rk)
        end do
      end do
      do component = 1, 3
        expected_mean_field(component, point) = &
          0.003_rk * real(10 * component - point, rk)
      end do
    end do
  end do

  call write_order_row(op_x_path, 1, expected_order_parameter)
  call write_order_row(op_y_path, 2, expected_order_parameter)
  call write_order_row(op_z_path, 3, expected_order_parameter)
  call write_current(current_path, expected_order_parameter, expected_mean_field)

  call read_legacy_split_field_map_2d( &
    op_x_path, op_y_path, op_z_path, current_path, mesh, state, report)

  call require(mesh%x_point_count() == 3 .and. mesh%y_point_count() == 2, &
               "legacy split reader inferred the wrong grid shape")
  do x_node = 0, mesh%x_cell_count()
    call require(abs(mesh%x_coordinate(x_node) - &
                         x_coordinates(x_node + 1)) < 1.0e-15_rk, &
                 "legacy split reader inferred the wrong x coordinate")
  end do
  do y_node = 0, mesh%y_cell_count()
    call require(abs(mesh%y_coordinate(y_node) - &
                         y_coordinates(y_node + 1)) < 1.0e-15_rk, &
                 "legacy split reader inferred the wrong y coordinate")
  end do
  call require(maxval(abs(state%order_parameter - &
                              expected_order_parameter)) < 1.0e-14_rk, &
               "legacy split order parameter was not preserved")
  call require(maxval(abs(state%current_mean_field - &
                              expected_mean_field)) < 1.0e-14_rk, &
               "legacy split mean field was not preserved")
  call require(report%maximum_gap_norm_difference < 1.0e-14_rk, &
               "legacy split stored gap norm was interpreted incorrectly")
  call require(report%maximum_in_plane_magnitude_difference < 1.0e-14_rk, &
               "legacy split in-plane magnitude was interpreted incorrectly")

  call delete_file(op_x_path)
  call delete_file(op_y_path)
  call delete_file(op_z_path)
  call delete_file(current_path)
  print '(a)', "legacy split 2D field I/O tests passed"

contains

  subroutine write_order_row(path, spin, order_parameter)
    character(len=*), intent(in) :: path
    integer, intent(in) :: spin
    complex(rk), intent(in) :: order_parameter(:, :, :)

    integer :: orbital, output_unit, point, x_node, y_node
    real(rk) :: row(8)

    open(newunit=output_unit, file=path, status="replace", action="write")
    point = 0
    do y_node = 1, size(y_coordinates)
      do x_node = 1, size(x_coordinates)
        point = point + 1
        row = 0.0_rk
        row(1:2) = [x_coordinates(x_node), y_coordinates(y_node)]
        do orbital = 1, 3
          row(2 * orbital + 1) = &
            real(order_parameter(spin, orbital, point), rk)
          row(2 * orbital + 2) = aimag(order_parameter(spin, orbital, point))
        end do
        write(output_unit, '(*(es24.16e3,1x))') row
      end do
      write(output_unit, '(a)') ""
    end do
    close(output_unit)
  end subroutine write_order_row


  subroutine write_current(path, order_parameter, mean_field)
    character(len=*), intent(in) :: path
    complex(rk), intent(in) :: order_parameter(:, :, :)
    real(rk), intent(in) :: mean_field(:, :)

    integer :: output_unit, point, x_node, y_node
    real(rk) :: row(7)

    open(newunit=output_unit, file=path, status="replace", action="write")
    point = 0
    do y_node = 1, size(y_coordinates)
      do x_node = 1, size(x_coordinates)
        point = point + 1
        row(1:2) = [x_coordinates(x_node), y_coordinates(y_node)]
        row(3) = sqrt(sum(abs(order_parameter(:, :, point))**2) / 3.0_rk)
        row(4) = hypot(mean_field(1, point), mean_field(2, point))
        row(5:7) = mean_field(:, point)
        write(output_unit, '(*(es24.16e3,1x))') row
      end do
      write(output_unit, '(a)') ""
    end do
    close(output_unit)
  end subroutine write_current


  subroutine delete_file(path)
    character(len=*), intent(in) :: path

    integer :: unit

    open(newunit=unit, file=path, status="old")
    close(unit, status="delete")
  end subroutine delete_file


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_legacy_split_field_io_2d
