module spinful_state_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t
  implicit none
  private

  integer, parameter, public :: spin_components_2d = 3
  integer, parameter, public :: orbital_components_2d = 3
  integer, parameter, public :: mean_field_components_2d = 3

  type, public :: spinful_state_2d_t
    ! The logical field order follows A(spin, orbital, point). The flat point
    ! index belongs to the mesh and is independent of nonlinear-solver packing.
    complex(rk), allocatable :: order_parameter(:, :, :)
    real(rk), allocatable :: current_mean_field(:, :)
  contains
    procedure :: point_count => state_point_count
    procedure :: is_valid_for => state_is_valid_for
    procedure :: set_zero => zero_state
  end type spinful_state_2d_t

  public :: allocate_spinful_state_2d

contains

  subroutine allocate_spinful_state_2d(mesh, state)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(out) :: state

    integer :: number_of_points

    if (.not. mesh%is_valid()) &
      error stop "cannot allocate spinful fields on an invalid Cartesian mesh"

    number_of_points = mesh%point_count()
    allocate(state%order_parameter(spin_components_2d, orbital_components_2d, &
                                   number_of_points), &
             source=cmplx(0.0_rk, 0.0_rk, kind=rk))
    allocate(state%current_mean_field(mean_field_components_2d, number_of_points), &
             source=0.0_rk)
  end subroutine allocate_spinful_state_2d


  pure integer function state_point_count(self) result(number_of_points)
    class(spinful_state_2d_t), intent(in) :: self

    if (allocated(self%order_parameter)) then
      number_of_points = size(self%order_parameter, 3)
    else
      number_of_points = 0
    end if
  end function state_point_count


  pure logical function state_is_valid_for(self, mesh) result(valid)
    class(spinful_state_2d_t), intent(in) :: self
    type(cartesian_mesh_2d_t), intent(in) :: mesh

    valid = mesh%is_valid() .and. allocated(self%order_parameter) .and. &
            allocated(self%current_mean_field)
    if (.not. valid) return

    valid = size(self%order_parameter, 1) == spin_components_2d .and. &
            size(self%order_parameter, 2) == orbital_components_2d .and. &
            size(self%order_parameter, 3) == mesh%point_count() .and. &
            size(self%current_mean_field, 1) == mean_field_components_2d .and. &
            size(self%current_mean_field, 2) == mesh%point_count()
  end function state_is_valid_for


  subroutine zero_state(self)
    class(spinful_state_2d_t), intent(inout) :: self

    if (allocated(self%order_parameter)) &
      self%order_parameter = cmplx(0.0_rk, 0.0_rk, kind=rk)
    if (allocated(self%current_mean_field)) self%current_mean_field = 0.0_rk
  end subroutine zero_state

end module spinful_state_2d
