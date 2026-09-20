program test_spinful_field_io_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, &
                                make_symmetric_multiscale_cartesian_mesh
  use spinful_state_2d, only : spinful_state_2d_t, allocate_spinful_state_2d
  use spinful_field_io_2d, only : read_spinful_field_map_2d, &
                                  write_spinful_field_map_2d
  implicit none

  character(len=*), parameter :: path = "test_spinful_field_io_2d.tmp"
  type(cartesian_mesh_2d_t) :: mesh
  type(spinful_state_2d_t) :: actual, expected
  integer :: component, orbital, point, spin, unit

  call make_symmetric_multiscale_cartesian_mesh( &
    4.0_rk, 1.0_rk, 2.5_rk, 0.5_rk, 0.75_rk, 1.5_rk, mesh)
  call allocate_spinful_state_2d(mesh, expected)
  do point = 1, mesh%point_count()
    do spin = 1, 3
      do orbital = 1, 3
        expected%order_parameter(spin, orbital, point) = cmplx( &
          0.01_rk * real(100 * spin + 10 * orbital + point, rk), &
          -0.02_rk * real(10 * spin + orbital + point, rk), kind=rk)
      end do
    end do
    do component = 1, 3
      expected%current_mean_field(component, point) = &
        0.003_rk * real(10 * component - point, rk)
    end do
  end do

  call write_spinful_field_map_2d( &
    path, mesh, expected, "round_trip_test", "radial_reference")
  call allocate_spinful_state_2d(mesh, actual)
  call read_spinful_field_map_2d(path, mesh, actual)

  call require(maxval(abs(actual%order_parameter - &
                              expected%order_parameter)) < 1.0e-15_rk, &
               "2D order parameter did not survive a field-map round trip")
  call require(maxval(abs(actual%current_mean_field - &
                              expected%current_mean_field)) < 1.0e-15_rk, &
               "2D mean field did not survive a field-map round trip")

  open(newunit=unit, file=path, status="old")
  close(unit, status="delete")
  print '(a)', "spinful 2D field-map I/O tests passed"

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_spinful_field_io_2d
