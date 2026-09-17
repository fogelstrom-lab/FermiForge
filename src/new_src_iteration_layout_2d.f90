module new_src_iteration_layout_2d
  use he3_kinds, only : rk
  use spinful_state_2d, only : spinful_state_2d_t, spin_components_2d, &
                               orbital_components_2d, mean_field_components_2d
  implicit none
  private

  integer, parameter, public :: new_src_values_per_2d_point = &
    2 * spin_components_2d * orbital_components_2d + mean_field_components_2d

  public :: new_src_iteration_vector_size_2d
  public :: pack_new_src_iteration_vector_2d
  public :: unpack_new_src_iteration_vector_2d

contains

  pure integer function new_src_iteration_vector_size_2d(state) result(vector_size)
    type(spinful_state_2d_t), intent(in) :: state

    vector_size = new_src_values_per_2d_point * state%point_count()
  end function new_src_iteration_vector_size_2d


  subroutine pack_new_src_iteration_vector_2d(state, vector)
    type(spinful_state_2d_t), intent(in) :: state
    real(rk), allocatable, intent(out) :: vector(:)

    integer :: block_end, block_start, component, number_of_points, orbital, spin

    call require_complete_state(state)
    number_of_points = state%point_count()
    allocate(vector(new_src_iteration_vector_size_2d(state)))

    block_end = 0
    do spin = 1, spin_components_2d
      do orbital = 1, orbital_components_2d
        block_start = block_end + 1
        block_end = block_end + number_of_points
        vector(block_start:block_end) = &
          real(state%order_parameter(spin, orbital, :), kind=rk)

        block_start = block_end + 1
        block_end = block_end + number_of_points
        vector(block_start:block_end) = &
          aimag(state%order_parameter(spin, orbital, :))
      end do
    end do

    do component = 1, mean_field_components_2d
      block_start = block_end + 1
      block_end = block_end + number_of_points
      vector(block_start:block_end) = state%current_mean_field(component, :)
    end do
  end subroutine pack_new_src_iteration_vector_2d


  subroutine unpack_new_src_iteration_vector_2d(vector, state)
    real(rk), intent(in) :: vector(:)
    type(spinful_state_2d_t), intent(inout) :: state

    integer :: block_end, block_start, component, number_of_points, orbital, spin
    real(rk), allocatable :: real_part(:)

    call require_complete_state(state)
    if (size(vector) /= new_src_iteration_vector_size_2d(state)) &
      error stop "new_src 2D iteration vector has the wrong size"

    number_of_points = state%point_count()
    allocate(real_part(number_of_points))

    block_end = 0
    do spin = 1, spin_components_2d
      do orbital = 1, orbital_components_2d
        block_start = block_end + 1
        block_end = block_end + number_of_points
        real_part = vector(block_start:block_end)

        block_start = block_end + 1
        block_end = block_end + number_of_points
        state%order_parameter(spin, orbital, :) = &
          cmplx(real_part, vector(block_start:block_end), kind=rk)
      end do
    end do

    do component = 1, mean_field_components_2d
      block_start = block_end + 1
      block_end = block_end + number_of_points
      state%current_mean_field(component, :) = vector(block_start:block_end)
    end do
  end subroutine unpack_new_src_iteration_vector_2d


  subroutine require_complete_state(state)
    type(spinful_state_2d_t), intent(in) :: state

    if (.not. allocated(state%order_parameter) .or. &
        .not. allocated(state%current_mean_field)) &
      error stop "cannot pack or unpack an empty 2D spinful state"
    if (state%point_count() < 1) &
      error stop "cannot pack or unpack a zero-length 2D spinful state"
    if (size(state%order_parameter, 1) /= spin_components_2d .or. &
        size(state%order_parameter, 2) /= orbital_components_2d .or. &
        size(state%current_mean_field, 1) /= mean_field_components_2d .or. &
        size(state%current_mean_field, 2) /= state%point_count()) &
      error stop "2D spinful state has inconsistent component dimensions"
  end subroutine require_complete_state

end module new_src_iteration_layout_2d
