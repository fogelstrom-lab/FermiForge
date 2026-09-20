program test_historical_core_seed_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, &
    make_uniform_cartesian_mesh
  use historical_core_seed_2d, only : historical_core_seed_report_t, &
    initialize_historical_core_seed_2d
  use spinful_state_2d, only : spinful_state_2d_t
  implicit none

  type(cartesian_mesh_2d_t) :: mesh
  type(spinful_state_2d_t) :: nop_state, aop_state, dop_state
  type(historical_core_seed_report_t) :: nop_report, aop_report, dop_report
  integer :: center, x_axis, y_axis
  real(rk) :: aop_ar, aop_fill, aop_scale, dop_ar, dop_fill, dop_scale
  real(rk), parameter :: feedback = 1.9285714285714286_rk
  real(rk), parameter :: gap = 0.2799_rk
  real(rk), parameter :: temperature = 0.30_rk
  real(rk), parameter :: tolerance = 3.0e-13_rk

  call make_uniform_cartesian_mesh( &
    -2.0_rk, 2.0_rk, 4, -2.0_rk, 2.0_rk, 4, mesh)
  call initialize_historical_core_seed_2d( &
    mesh, "nop", gap, temperature, feedback, nop_state, nop_report)
  call initialize_historical_core_seed_2d( &
    mesh, "aop", gap, temperature, feedback, aop_state, aop_report)
  call initialize_historical_core_seed_2d( &
    mesh, "dop", gap, temperature, feedback, dop_state, dop_report)

  center = mesh%point_index(2, 2)
  x_axis = mesh%point_index(3, 2)
  y_axis = mesh%point_index(2, 3)
  aop_scale = 1.0_rk / sqrt(1.0_rk - temperature**2)
  dop_scale = (1.0_rk + feedback) * aop_scale
  aop_ar = 1.0_rk / (3.0_rk * aop_scale)
  dop_ar = 1.0_rk / (3.0_rk * dop_scale)
  aop_fill = tanh(aop_ar) / aop_ar
  dop_fill = tanh(dop_ar) / dop_ar

  call require(maxval(abs(nop_state%current_mean_field)) <= tiny(1.0_rk), &
    "historical nop mean field is not zero")
  call require(maxval(abs(aop_state%current_mean_field)) <= tiny(1.0_rk), &
    "historical aop mean field is not zero")
  call require(maxval(abs(dop_state%current_mean_field)) <= tiny(1.0_rk), &
    "historical dop mean field is not zero")

  call require(maxval(abs(nop_state%order_parameter(:, :, center))) <= &
               tiny(1.0_rk), "historical nop origin is not normal")
  call require_close(nop_state%order_parameter(1, 1, x_axis), &
    cmplx(gap * tanh(aop_ar), 0.0_rk, kind=rk), &
    "historical nop positive-x winding")
  call require_close(nop_state%order_parameter(1, 1, y_axis), &
    cmplx(0.0_rk, gap * tanh(aop_ar), kind=rk), &
    "historical nop positive-y winding")
  call require(maxval(abs(nop_state%order_parameter(:, :, x_axis) - &
    diagonal_matrix(nop_state%order_parameter(1, 1, x_axis)))) < tolerance, &
    "historical nop has a non-diagonal component")

  call require_close(aop_state%order_parameter(3, 1, center), &
    cmplx(gap, 0.0_rk, kind=rk), "historical aop origin A_zx")
  call require_close(aop_state%order_parameter(1, 3, center), &
    cmplx(-gap, 0.0_rk, kind=rk), "historical aop origin A_xz")
  call require_close(aop_state%order_parameter(3, 1, x_axis), &
    cmplx(gap * aop_fill, 0.0_rk, kind=rk), &
    "historical aop positive-x A_zx")
  call require_close(aop_state%order_parameter(3, 2, y_axis), &
    cmplx(0.0_rk, gap * aop_fill, kind=rk), &
    "historical aop positive-y A_zy")
  call require_close(aop_state%order_parameter(2, 3, y_axis), &
    -aop_state%order_parameter(3, 2, y_axis), &
    "historical aop lost yz/zy antisymmetry")

  call require_close(dop_state%order_parameter(3, 1, center), &
    cmplx(gap, 0.0_rk, kind=rk), "historical dop origin A_zx")
  call require_close(dop_state%order_parameter(3, 1, x_axis), &
    cmplx(gap * dop_fill, 0.0_rk, kind=rk), &
    "historical dop positive-x A_zx")
  call require(abs(dop_state%order_parameter(3, 1, y_axis)) < tolerance, &
    "historical dop cos(phi)^2 node is absent")
  call require_close(dop_state%order_parameter(1, 3, x_axis), &
    -dop_state%order_parameter(3, 1, x_axis), &
    "historical dop lost xz/zx antisymmetry")
  call require(abs(dop_report%core_length_scale - dop_scale) < tolerance, &
    "historical dop did not use the (1+aa0) radial scale")
  call require(abs(aop_report%core_length_scale - aop_scale) < tolerance, &
    "historical aop radial scale changed")
  call require(abs(nop_report%core_length_scale - aop_scale) < tolerance, &
    "historical nop radial scale changed")

  print '(a)', "historical core seed tests passed"

contains

  pure function diagonal_matrix(value) result(matrix)
    complex(rk), intent(in) :: value
    complex(rk) :: matrix(3, 3)
    integer :: component

    matrix = cmplx(0.0_rk, 0.0_rk, kind=rk)
    do component = 1, 3
      matrix(component, component) = value
    end do
  end function diagonal_matrix


  subroutine require_close(actual, expected, message)
    complex(rk), intent(in) :: actual, expected
    character(len=*), intent(in) :: message

    call require(abs(actual - expected) < tolerance, message)
  end subroutine require_close


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_historical_core_seed_2d
