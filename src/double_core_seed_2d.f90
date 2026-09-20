module double_core_seed_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t
  use spinful_state_2d, only : spinful_state_2d_t, &
    allocate_spinful_state_2d
  implicit none
  private

  type, public :: double_core_seed_report_t
    real(rk) :: center_pair_amplitude = 0.0_rk
    real(rk) :: positive_half_core_pair_amplitude = 0.0_rk
    real(rk) :: negative_half_core_pair_amplitude = 0.0_rk
    real(rk) :: outer_pair_amplitude = 0.0_rk
  end type double_core_seed_report_t

  public :: initialize_regularized_london_double_core_2d

contains

  subroutine initialize_regularized_london_double_core_2d( &
      mesh, bulk_gap, half_core_offset, core_width, domain_wall_width, &
      state, report)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    real(rk), intent(in) :: bulk_gap, half_core_offset
    real(rk), intent(in) :: core_width, domain_wall_width
    type(spinful_state_2d_t), intent(out) :: state
    type(double_core_seed_report_t), intent(out), optional :: report

    complex(rk) :: average_phase, relative_phase, vortex_phase
    real(rk) :: coordinate(2), distance_to_wall, radius
    real(rk) :: x, y, y_beyond_segment
    integer :: point

    if (.not. mesh%is_valid()) &
      error stop "cannot construct a double-core seed on an invalid mesh"
    if (bulk_gap <= 0.0_rk .or. half_core_offset <= 0.0_rk .or. &
        core_width <= 0.0_rk .or. domain_wall_width <= 0.0_rk) &
      error stop "regularized London seed scales must be positive"

    call allocate_spinful_state_2d(mesh, state)
    do point = 1, mesh%point_count()
      coordinate = mesh%point_coordinate(point)
      x = coordinate(1)
      y = coordinate(2)

      ! Equation (4) of dcvlong.pdf can be written without half-angle
      ! branches as
      !   exp(i phi) cos(delta) = (u_+ + u_-)/2,
      !   exp(i phi) sin(delta) = (u_+ - u_-)/(2 i),
      ! where u_+ and u_- are phases about the two half cores.  Replacing
      ! their unit denominators by sqrt(r^2+xi_c^2) regularizes each core.
      average_phase = cmplx(0.5_rk, 0.0_rk, kind=rk) * ( &
        regularized_phase(x, y - half_core_offset, core_width) + &
        regularized_phase(x, y + half_core_offset, core_width))
      relative_phase = cmplx(0.0_rk, -0.5_rk, kind=rk) * ( &
        regularized_phase(x, y + half_core_offset, core_width) - &
        regularized_phase(x, y - half_core_offset, core_width))

      ! The London A_yy component has a type-10 wall discontinuity between
      ! the half cores.  For a numerical seed, replace it by a continuous
      ! winding-one component that is suppressed over a wall of finite width.
      radius = sqrt(x * x + y * y)
      vortex_phase = cmplx(x, y, kind=rk) / cmplx( &
        sqrt(radius * radius + core_width * core_width), 0.0_rk, kind=rk)
      y_beyond_segment = max(abs(y) - half_core_offset, 0.0_rk)
      distance_to_wall = sqrt(x * x + y_beyond_segment * y_beyond_segment)

      state%order_parameter(1, 1, point) = &
        cmplx(bulk_gap, 0.0_rk, kind=rk) * average_phase
      state%order_parameter(3, 1, point) = &
        cmplx(bulk_gap, 0.0_rk, kind=rk) * relative_phase
      state%order_parameter(2, 2, point) = &
        cmplx(bulk_gap * tanh(distance_to_wall / domain_wall_width), &
              0.0_rk, kind=rk) * vortex_phase
      state%order_parameter(1, 3, point) = &
        -cmplx(bulk_gap, 0.0_rk, kind=rk) * relative_phase
      state%order_parameter(3, 3, point) = &
        cmplx(bulk_gap, 0.0_rk, kind=rk) * average_phase
    end do

    if (present(report)) call measure_seed(mesh, state, half_core_offset, report)
  end subroutine initialize_regularized_london_double_core_2d


  pure complex(rk) function regularized_phase(x, y, width) result(value)
    real(rk), intent(in) :: x, y, width

    value = cmplx(x, y, kind=rk) / cmplx( &
      sqrt(x * x + y * y + width * width), 0.0_rk, kind=rk)
  end function regularized_phase


  subroutine measure_seed(mesh, state, half_core_offset, report)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    real(rk), intent(in) :: half_core_offset
    type(double_core_seed_report_t), intent(out) :: report

    real(rk) :: outer_x

    outer_x = min(abs(mesh%x_minimum), abs(mesh%x_maximum))
    report%center_pair_amplitude = pair_amplitude_at(0.0_rk, 0.0_rk)
    report%positive_half_core_pair_amplitude = &
      pair_amplitude_at(0.0_rk, half_core_offset)
    report%negative_half_core_pair_amplitude = &
      pair_amplitude_at(0.0_rk, -half_core_offset)
    report%outer_pair_amplitude = pair_amplitude_at(outer_x, 0.0_rk)

  contains

    real(rk) function pair_amplitude_at(x, y) result(amplitude)
      real(rk), intent(in) :: x, y

      real(rk) :: coordinate(2), distance, minimum_distance
      integer :: point, selected

      selected = 1
      minimum_distance = huge(1.0_rk)
      do point = 1, mesh%point_count()
        coordinate = mesh%point_coordinate(point)
        distance = sum((coordinate - [x, y])**2)
        if (distance < minimum_distance) then
          selected = point
          minimum_distance = distance
        end if
      end do
      amplitude = sqrt(sum(abs(state%order_parameter(:, :, selected))**2) / &
        3.0_rk)
    end function pair_amplitude_at

  end subroutine measure_seed

end module double_core_seed_2d
