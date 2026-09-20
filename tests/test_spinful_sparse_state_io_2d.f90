program test_spinful_sparse_state_io_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, &
    make_rectilinear_cartesian_mesh
  use spinful_state_2d, only : spinful_state_2d_t, allocate_spinful_state_2d
  use spinful_sparse_state_io_2d, only : &
    write_sparse_spinful_state_2d, read_sparse_spinful_state_2d
  implicit none

  character(len=*), parameter :: path = "test_spinful_sparse_state.tmp"
  type(cartesian_mesh_2d_t) :: mesh
  type(spinful_state_2d_t) :: restored, source
  logical, allocatable :: mask(:)
  integer :: applied, component, orbital, point, spin, unit

  call make_rectilinear_cartesian_mesh( &
    [-1.0_rk, -0.2_rk, 0.7_rk, 1.5_rk], &
    [-0.9_rk, 0.1_rk, 1.2_rk], mesh)
  call allocate_spinful_state_2d(mesh, source)
  do point = 1, mesh%point_count()
    do spin = 1, 3
      do orbital = 1, 3
        source%order_parameter(spin, orbital, point) = cmplx( &
          0.01_rk * real(100 * spin + 10 * orbital + point, rk), &
          -0.02_rk * real(10 * spin + orbital + point, rk), rk)
      end do
    end do
    do component = 1, 3
      source%current_mean_field(component, point) = &
        0.03_rk * real(10 * component + point, rk)
    end do
  end do
  allocate(mask(mesh%point_count()), source=.false.)
  mask(2) = .true.
  mask(7) = .true.
  mask(mesh%point_count()) = .true.
  call write_sparse_spinful_state_2d(path, mesh, source, mask, "round_trip")

  call allocate_spinful_state_2d(mesh, restored)
  restored%order_parameter = cmplx(-9.0_rk, 4.0_rk, rk)
  restored%current_mean_field = -7.0_rk
  call read_sparse_spinful_state_2d(path, mesh, restored, applied)
  call require(applied == count(mask), &
    "sparse state reader returned the wrong record count")
  do point = 1, mesh%point_count()
    if (mask(point)) then
      call require(maxval(abs(restored%order_parameter(:, :, point) - &
                              source%order_parameter(:, :, point))) < &
                     2.0e-15_rk, &
        "sparse state round trip changed an order-parameter value")
      call require(maxval(abs(restored%current_mean_field(:, point) - &
                              source%current_mean_field(:, point))) < &
                     2.0e-15_rk, &
        "sparse state round trip changed a mean-field value")
    else
      call require(maxval(abs(restored%order_parameter(:, :, point) - &
                              cmplx(-9.0_rk, 4.0_rk, rk))) < &
                     2.0e-15_rk .and. &
                   maxval(abs(restored%current_mean_field(:, point) + &
                              7.0_rk)) < 2.0e-15_rk, &
        "sparse state reader changed a point absent from the checkpoint")
    end if
  end do

  open(newunit=unit, file=path, status="old")
  close(unit, status="delete")
  print '(a)', "spinful sparse-state I/O tests passed"

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_spinful_sparse_state_io_2d
