module straight_trajectory_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t
  use spinful_state_2d, only : spinful_state_2d_t
  use bilinear_field_sampler_2d, only : bilinear_stencil_2d_t, &
                                        make_bilinear_stencil_2d, &
                                        sample_pair_potential_stencil_2d
  implicit none
  private

  type, public :: straight_trajectory_2d_t
    real(rk) :: origin(2) = 0.0_rk
    real(rk) :: momentum(3) = 0.0_rk
    real(rk) :: entry_coordinate = 0.0_rk
    real(rk) :: exit_coordinate = 0.0_rk
    integer :: target_sample = 0
    logical :: parallel_to_invariant_axis = .false.
    real(rk), allocatable :: path_coordinate(:)
    type(bilinear_stencil_2d_t), allocatable :: stencil(:)
  contains
    procedure :: sample_count => trajectory_sample_count
    procedure :: is_valid => trajectory_is_valid
  end type straight_trajectory_2d_t

  public :: build_straight_trajectory_2d
  public :: sample_pair_potential_along_trajectory_2d

contains

  subroutine build_straight_trajectory_2d(mesh, origin, momentum, maximum_step, &
                                          trajectory)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    real(rk), intent(in) :: origin(2), momentum(3), maximum_step
    type(straight_trajectory_2d_t), intent(out) :: trajectory

    type(bilinear_stencil_2d_t) :: origin_stencil
    real(rk) :: direction_norm, negative_length, planar_norm, positive_length
    real(rk) :: position(2), s_lower, s_upper, tolerance
    integer :: index, number_negative, number_positive, number_samples

    if (.not. mesh%is_valid()) &
      error stop "cannot build a trajectory on an invalid Cartesian mesh"
    if (maximum_step <= 0.0_rk) &
      error stop "trajectory maximum step must be positive"

    direction_norm = sqrt(dot_product(momentum, momentum))
    tolerance = 128.0_rk * epsilon(1.0_rk)
    if (abs(direction_norm - 1.0_rk) > tolerance) &
      error stop "trajectory momentum direction must have unit norm"

    call make_bilinear_stencil_2d(mesh, origin(1), origin(2), origin_stencil)
    if (.not. origin_stencil%inside) &
      error stop "trajectory origin lies outside the Cartesian mesh"

    trajectory%origin = origin
    trajectory%momentum = momentum
    planar_norm = sqrt(momentum(1) * momentum(1) + momentum(2) * momentum(2))
    if (planar_norm <= tolerance) then
      trajectory%parallel_to_invariant_axis = .true.
      trajectory%entry_coordinate = 0.0_rk
      trajectory%exit_coordinate = 0.0_rk
      trajectory%target_sample = 1
      allocate(trajectory%path_coordinate(1), source=0.0_rk)
      allocate(trajectory%stencil(1))
      trajectory%stencil(1) = origin_stencil
      return
    end if

    s_lower = -huge(1.0_rk)
    s_upper = huge(1.0_rk)
    call restrict_line_interval(origin(1), momentum(1), mesh%x_minimum, &
                                mesh%x_maximum, tolerance, s_lower, s_upper)
    call restrict_line_interval(origin(2), momentum(2), mesh%y_minimum, &
                                mesh%y_maximum, tolerance, s_lower, s_upper)
    if (s_lower > tolerance .or. s_upper < -tolerance .or. s_lower > s_upper) &
      error stop "trajectory does not cross its declared origin"

    if (abs(s_lower) <= tolerance) s_lower = 0.0_rk
    if (abs(s_upper) <= tolerance) s_upper = 0.0_rk
    trajectory%entry_coordinate = s_lower
    trajectory%exit_coordinate = s_upper

    negative_length = -s_lower
    positive_length = s_upper
    number_negative = 0
    number_positive = 0
    if (negative_length > tolerance) &
      number_negative = ceiling(negative_length / maximum_step)
    if (positive_length > tolerance) &
      number_positive = ceiling(positive_length / maximum_step)

    number_samples = number_negative + 1 + number_positive
    trajectory%target_sample = number_negative + 1
    allocate(trajectory%path_coordinate(number_samples))
    allocate(trajectory%stencil(number_samples))

    if (number_negative > 0) then
      do index = 0, number_negative
        trajectory%path_coordinate(index + 1) = s_lower + negative_length * &
          real(index, rk) / real(number_negative, rk)
      end do
    else
      trajectory%path_coordinate(1) = 0.0_rk
    end if
    do index = 1, number_positive
      trajectory%path_coordinate(trajectory%target_sample + index) = &
        positive_length * real(index, rk) / real(number_positive, rk)
    end do

    trajectory%path_coordinate(trajectory%target_sample) = 0.0_rk
    do index = 1, number_samples
      position = origin + trajectory%path_coordinate(index) * momentum(1:2)
      call make_bilinear_stencil_2d(mesh, position(1), position(2), &
                                    trajectory%stencil(index))
      if (.not. trajectory%stencil(index)%inside) &
        error stop "round-off moved a trajectory sample outside the mesh"
    end do
  end subroutine build_straight_trajectory_2d


  subroutine sample_pair_potential_along_trajectory_2d(state, trajectory, &
                                                       pair_potential, &
                                                       current_mean_field)
    type(spinful_state_2d_t), intent(in) :: state
    type(straight_trajectory_2d_t), intent(in) :: trajectory
    complex(rk), allocatable, intent(out) :: pair_potential(:, :)
    real(rk), allocatable, intent(out) :: current_mean_field(:, :)

    logical :: inside
    integer :: sample

    if (.not. trajectory%is_valid()) &
      error stop "cannot sample an invalid 2D trajectory"

    allocate(pair_potential(3, trajectory%sample_count()), &
             source=cmplx(0.0_rk, 0.0_rk, kind=rk))
    allocate(current_mean_field(3, trajectory%sample_count()), source=0.0_rk)
    do sample = 1, trajectory%sample_count()
      call sample_pair_potential_stencil_2d( &
        state, trajectory%stencil(sample), trajectory%momentum, &
        pair_potential(:, sample), current_mean_field(:, sample), inside)
      if (.not. inside) &
        error stop "cached trajectory stencil lies outside the 2D field"
    end do
  end subroutine sample_pair_potential_along_trajectory_2d


  pure subroutine restrict_line_interval(origin, direction, lower_bound, &
                                         upper_bound, tolerance, s_lower, s_upper)
    real(rk), intent(in) :: origin, direction, lower_bound, upper_bound, tolerance
    real(rk), intent(inout) :: s_lower, s_upper

    real(rk) :: first, second

    if (abs(direction) <= tolerance) then
      if (origin < lower_bound - tolerance .or. origin > upper_bound + tolerance) &
        error stop "line parallel to a domain side lies outside the domain"
      return
    end if

    first = (lower_bound - origin) / direction
    second = (upper_bound - origin) / direction
    s_lower = max(s_lower, min(first, second))
    s_upper = min(s_upper, max(first, second))
  end subroutine restrict_line_interval


  pure integer function trajectory_sample_count(self) result(number_of_samples)
    class(straight_trajectory_2d_t), intent(in) :: self

    if (allocated(self%path_coordinate)) then
      number_of_samples = size(self%path_coordinate)
    else
      number_of_samples = 0
    end if
  end function trajectory_sample_count


  pure logical function trajectory_is_valid(self) result(valid)
    class(straight_trajectory_2d_t), intent(in) :: self

    valid = allocated(self%path_coordinate) .and. allocated(self%stencil)
    if (.not. valid) return
    valid = size(self%path_coordinate) == size(self%stencil) .and. &
            size(self%path_coordinate) >= 1 .and. &
            self%target_sample >= 1 .and. &
            self%target_sample <= size(self%path_coordinate)
    if (.not. valid) return
    valid = abs(self%path_coordinate(self%target_sample)) <= &
            128.0_rk * epsilon(1.0_rk) .and. all(self%stencil%inside)
  end function trajectory_is_valid

end module straight_trajectory_2d
