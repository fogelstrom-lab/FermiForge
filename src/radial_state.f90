module radial_state
  use he3_kinds, only : rk
  use radial_mesh, only : radial_mesh_t
  implicit none
  private

  integer, parameter, public :: spin_components = 3
  integer, parameter, public :: orbital_components = 3
  integer, parameter, public :: exchange_components = 3

  type, public :: radial_state_t
    ! The point index is the contiguous dimension. This favors batched field
    ! operations over many spatial points on CPUs and accelerators.
    complex(rk), allocatable :: order_parameter(:, :, :)
    complex(rk), allocatable :: exchange_field(:, :)
  contains
    procedure :: point_count => state_point_count
    procedure :: is_valid_for => state_is_valid_for
  end type radial_state_t

  public :: allocate_radial_state
  public :: legacy_iteration_vector_size
  public :: pack_legacy_iteration_vector
  public :: unpack_legacy_iteration_vector

contains

  subroutine allocate_radial_state(mesh, state)
    type(radial_mesh_t), intent(in) :: mesh
    type(radial_state_t), intent(out) :: state

    integer :: last_point

    if (.not. mesh%is_valid()) error stop "cannot allocate fields on an invalid mesh"

    last_point = mesh%point_count() - 1
    allocate(state%order_parameter(0:last_point, spin_components, orbital_components), &
             source=cmplx(0.0_rk, 0.0_rk, kind=rk))
    allocate(state%exchange_field(0:last_point, exchange_components), &
             source=cmplx(0.0_rk, 0.0_rk, kind=rk))
  end subroutine allocate_radial_state


  pure integer function legacy_iteration_vector_size(state, include_exchange) &
      result(vector_size)
    type(radial_state_t), intent(in) :: state
    logical, intent(in) :: include_exchange

    integer :: values_per_point

    values_per_point = 2 * spin_components * orbital_components
    if (include_exchange) values_per_point = values_per_point + exchange_components
    vector_size = values_per_point * state%point_count()
  end function legacy_iteration_vector_size


  subroutine pack_legacy_iteration_vector(state, include_exchange, vector)
    type(radial_state_t), intent(in) :: state
    logical, intent(in) :: include_exchange
    real(rk), allocatable, intent(out) :: vector(:)

    integer :: component, i, orbital, spin

    if (state%point_count() < 1) error stop "cannot pack an empty radial state"

    allocate(vector(legacy_iteration_vector_size(state, include_exchange)))
    component = 0
    do i = 0, state%point_count() - 1
      do spin = 1, spin_components
        do orbital = 1, orbital_components
          component = component + 1
          vector(component) = real(state%order_parameter(i, spin, orbital), kind=rk)
          component = component + 1
          vector(component) = aimag(state%order_parameter(i, spin, orbital))
        end do
      end do

      if (include_exchange) then
        do spin = 1, exchange_components
          component = component + 1
          ! The factor of ten and omission of the imaginary part exactly match
          ! legacy makearray. They are preserved here, not reinterpreted.
          vector(component) = 10.0_rk * real(state%exchange_field(i, spin), kind=rk)
        end do
      end if
    end do
  end subroutine pack_legacy_iteration_vector


  subroutine unpack_legacy_iteration_vector(vector, include_exchange, state)
    real(rk), intent(in) :: vector(:)
    logical, intent(in) :: include_exchange
    type(radial_state_t), intent(inout) :: state

    integer :: component, i, orbital, spin

    if (size(vector) /= legacy_iteration_vector_size(state, include_exchange)) &
      error stop "legacy iteration vector has the wrong size"

    component = 0
    do i = 0, state%point_count() - 1
      do spin = 1, spin_components
        do orbital = 1, orbital_components
          component = component + 2
          state%order_parameter(i, spin, orbital) = &
            cmplx(vector(component - 1), vector(component), kind=rk)
        end do
      end do

      if (include_exchange) then
        do spin = 1, exchange_components
          component = component + 1
          ! Legacy makeself reconstructs a real exchange field and divides by
          ! ten. Preserve that behavior until its physical convention is
          ! separately documented and tested.
          state%exchange_field(i, spin) = &
            cmplx(vector(component) / 10.0_rk, 0.0_rk, kind=rk)
        end do
      end if
    end do
  end subroutine unpack_legacy_iteration_vector


  pure integer function state_point_count(self) result(number_of_points)
    class(radial_state_t), intent(in) :: self

    if (allocated(self%order_parameter)) then
      number_of_points = size(self%order_parameter, 1)
    else
      number_of_points = 0
    end if
  end function state_point_count


  pure logical function state_is_valid_for(self, mesh) result(valid)
    class(radial_state_t), intent(in) :: self
    type(radial_mesh_t), intent(in) :: mesh

    valid = allocated(self%order_parameter) .and. allocated(self%exchange_field)
    if (.not. valid) return
    valid = lbound(self%order_parameter, 1) == 0 .and. &
            lbound(self%exchange_field, 1) == 0 .and. &
            size(self%order_parameter, 1) == mesh%point_count() .and. &
            size(self%order_parameter, 2) == spin_components .and. &
            size(self%order_parameter, 3) == orbital_components .and. &
            size(self%exchange_field, 1) == mesh%point_count() .and. &
            size(self%exchange_field, 2) == exchange_components
  end function state_is_valid_for

end module radial_state
