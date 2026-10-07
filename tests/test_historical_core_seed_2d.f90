program test_historical_core_seed_2d
  use he3_kinds, only : rk
  use order_parameter_basis, only : cartesian_to_axial_harmonics
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, &
    make_uniform_cartesian_mesh
  use historical_core_seed_2d, only : historical_core_seed_report_t, &
    initialize_historical_core_seed_2d, initialize_localized_harmonic_seed_2d
  use spinful_state_2d, only : spinful_state_2d_t
  implicit none

  type(cartesian_mesh_2d_t) :: mesh
  type(spinful_state_2d_t) :: nop_state, aop_state, dop_state
  type(spinful_state_2d_t) :: qop_state
  type(spinful_state_2d_t) :: localized, background
  type(cartesian_mesh_2d_t) :: wide_mesh
  character(len=16), parameter :: modes(4)=[character(len=16) :: &
    'localized_0plus','localized_plus0','localized_0minus','localized_minus0']
  integer, parameter :: spins(4)=[2,1,2,3], orbitals(4)=[1,2,3,2]
  integer :: mode
  complex(rk) :: expected(3,3)
  complex(rk) :: harmonic(3,3)
  real(rk) :: xy(2), radius, fill
  integer :: point
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
  call initialize_historical_core_seed_2d( &
    mesh, "qop", gap, temperature, feedback, qop_state)

  center = mesh%point_index(2, 2)
  x_axis = mesh%point_index(3, 2)
  y_axis = mesh%point_index(2, 3)
  aop_scale = 1.0_rk / sqrt(1.0_rk - temperature**2)
  do point=1,mesh%point_count()
    xy=mesh%point_coordinate(point)
    radius=sqrt(sum(xy**2))/(3.0_rk*aop_scale)
    fill=1.0_rk
    if (radius > 0.0_rk) fill=tanh(radius)/radius
    call cartesian_to_axial_harmonics(qop_state%order_parameter(:,:,point),harmonic)
    call require_close(harmonic(2,3),cmplx(sqrt(2.0_rk)*gap*fill,0.0_rk,rk), &
      'quadrupole C_0- must be non-winding and finite at origin')
    call require_close(harmonic(3,2),-harmonic(2,3),'quadrupole C_-0 relative sign')
    call require(abs(harmonic(2,1))+abs(harmonic(1,2)) < tolerance, &
      'quadrupole seed must have zero C_0+ and C_+0')
    call require(maxval(abs([(qop_state%order_parameter(1,1,point)-nop_state%order_parameter(1,1,point)), &
      (qop_state%order_parameter(2,2,point)-nop_state%order_parameter(2,2,point)), &
      (qop_state%order_parameter(3,3,point)-nop_state%order_parameter(3,3,point))])) < tolerance, &
      'quadrupole seed changed the bulk winding')
  end do
  call require(maxval(abs(qop_state%current_mean_field)) <= tiny(1.0_rk), &
    'quadrupole initial mean field must vanish')
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

  call make_uniform_cartesian_mesh(-6.0_rk,6.0_rk,12,-6.0_rk,6.0_rk,12,wide_mesh)
  call initialize_historical_core_seed_2d(wide_mesh,'nop',gap,temperature,feedback,background)
  do mode=1,4
    call initialize_localized_harmonic_seed_2d( &
      wide_mesh,trim(modes(mode)),gap,temperature,5.0_rk,1.0_rk,localized)
    do point=1,wide_mesh%point_count()
      xy=wide_mesh%point_coordinate(point)
      radius=sqrt(sum(xy**2))
      expected=cmplx(0.0_rk,0.0_rk,rk)
      if(radius<5.0_rk) expected(spins(mode),orbitals(mode))= &
        cmplx(gap*(1.0_rk-(radius/5.0_rk)**2)**2,0.0_rk,rk)
      call cartesian_to_axial_harmonics( &
        localized%order_parameter(:,:,point)-background%order_parameter(:,:,point),harmonic)
      call require(maxval(abs(harmonic-expected))<tolerance, &
        'localized seed harmonic selection, phase or envelope incorrect')
      if(radius>=5.0_rk) call require( &
        maxval(abs(localized%order_parameter(:,:,point)-background%order_parameter(:,:,point))) &
        <=tiny(1.0_rk),'localized seed leaks beyond support')
    end do
    call require(maxval(abs(localized%current_mean_field))<=tiny(1.0_rk), &
      'localized initial mean field must be zero')
  end do
  print '(a)', "historical and localized core seed tests passed"

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
