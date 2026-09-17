program test_radial_state
  use he3_kinds, only : rk
  use radial_mesh, only : radial_mesh_t, make_uniform_radial_mesh
  use radial_state, only : radial_state_t, allocate_radial_state, &
                           legacy_iteration_vector_size, &
                           pack_legacy_iteration_vector, &
                           unpack_legacy_iteration_vector
  implicit none

  type(radial_mesh_t) :: mesh
  type(radial_state_t) :: original, restored
  real(rk), allocatable :: vector(:)
  integer :: i, orbital, spin

  call make_uniform_radial_mesh(20.0_rk, 49, mesh)
  call allocate_radial_state(mesh, original)
  call allocate_radial_state(mesh, restored)

  call require(original%is_valid_for(mesh), "new radial state is invalid")

  do i = 0, mesh%point_count() - 1
    do spin = 1, 3
      do orbital = 1, 3
        original%order_parameter(i, spin, orbital) = &
          cmplx(100.0_rk * real(i, rk) + 10.0_rk * real(spin, rk) + &
                real(orbital, rk), &
                -100.0_rk * real(i, rk) - 10.0_rk * real(spin, rk) - &
                real(orbital, rk), kind=rk)
      end do
      original%exchange_field(i, spin) = &
        cmplx(real(i + spin, rk) / 7.0_rk, real(i - spin, rk), kind=rk)
    end do
  end do

  call pack_legacy_iteration_vector(original, .true., vector)
  call require(size(vector) == 21 * mesh%point_count(), &
               "legacy vector size with exchange field is wrong")
  call require(abs(vector(1) - real(original%order_parameter(0, 1, 1))) < &
               epsilon(1.0_rk), "dxx real component is in the wrong position")
  call require(abs(vector(2) - aimag(original%order_parameter(0, 1, 1))) < &
               epsilon(1.0_rk), "dxx imaginary component is in the wrong position")
  call require(abs(vector(19) - 10.0_rk * real(original%exchange_field(0, 1))) < &
               epsilon(1.0_rk), "vx component is in the wrong position")

  call unpack_legacy_iteration_vector(vector, .true., restored)
  call require(maxval(abs(restored%order_parameter - original%order_parameter)) < &
               epsilon(1.0_rk), "order-parameter pack/unpack did not round-trip")
  call require(maxval(abs(real(restored%exchange_field) - &
                          real(original%exchange_field))) < &
               64.0_rk * epsilon(1.0_rk) * &
               max(1.0_rk, maxval(abs(real(original%exchange_field)))), &
               "exchange-field pack/unpack did not preserve its real part")
  call require(maxval(abs(aimag(restored%exchange_field))) < epsilon(1.0_rk), &
               "legacy unpack must produce a real exchange field")

  call pack_legacy_iteration_vector(original, .false., vector)
  call require(size(vector) == legacy_iteration_vector_size(original, .false.), &
               "legacy vector size without exchange field is wrong")
  call require(size(vector) == 18 * mesh%point_count(), &
               "legacy order-parameter-only vector must have 18 values per point")

  print '(a)', "radial state tests passed"

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_radial_state
