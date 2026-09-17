program test_order_parameter_basis
  use he3_kinds, only : rk
  use order_parameter_basis, only : projection_plus, projection_zero, &
                                    projection_minus, &
                                    cartesian_to_axial_harmonics, &
                                    axial_harmonics_to_cartesian, &
                                    axial_harmonics_to_transport_cartesian, &
                                    phase_axial_harmonics, &
                                    reconstruct_axial_cartesian
  implicit none

  call test_cartesian_component_formulas()
  call test_complex_round_trip()
  call test_radial_transport_projection()
  call test_axial_phase_factors()
  call test_b_phase_far_field()

  print '(a)', "order-parameter basis tests passed"

contains

  subroutine test_cartesian_component_formulas()
    complex(rk) :: cartesian(3, 3), harmonic(3, 3)
    real(rk), parameter :: inverse_sqrt_two = 1.0_rk / sqrt(2.0_rk)

    cartesian = cmplx(0.0_rk, 0.0_rk, kind=rk)
    cartesian(1, 1) = cmplx(1.0_rk, 0.0_rk, kind=rk)
    call cartesian_to_axial_harmonics(cartesian, harmonic)
    call require_close(harmonic(projection_plus, projection_plus), &
                       cmplx(0.5_rk, 0.0_rk, kind=rk), &
                       "Axx to A++ formula differs from new_src")
    call require_close(harmonic(projection_plus, projection_minus), &
                       cmplx(0.5_rk, 0.0_rk, kind=rk), &
                       "Axx to A+- formula differs from new_src")
    call require_close(harmonic(projection_minus, projection_plus), &
                       cmplx(0.5_rk, 0.0_rk, kind=rk), &
                       "Axx to A-+ formula differs from new_src")
    call require_close(harmonic(projection_minus, projection_minus), &
                       cmplx(0.5_rk, 0.0_rk, kind=rk), &
                       "Axx to A-- formula differs from new_src")

    cartesian = cmplx(0.0_rk, 0.0_rk, kind=rk)
    cartesian(1, 2) = cmplx(1.0_rk, 0.0_rk, kind=rk)
    call cartesian_to_axial_harmonics(cartesian, harmonic)
    call require_close(harmonic(projection_plus, projection_plus), &
                       cmplx(0.0_rk, -0.5_rk, kind=rk), &
                       "Axy to A++ sign differs from new_src")
    call require_close(harmonic(projection_plus, projection_minus), &
                       cmplx(0.0_rk, 0.5_rk, kind=rk), &
                       "Axy to A+- sign differs from new_src")
    call require_close(harmonic(projection_minus, projection_plus), &
                       cmplx(0.0_rk, -0.5_rk, kind=rk), &
                       "Axy to A-+ sign differs from new_src")
    call require_close(harmonic(projection_minus, projection_minus), &
                       cmplx(0.0_rk, 0.5_rk, kind=rk), &
                       "Axy to A-- sign differs from new_src")

    cartesian = cmplx(0.0_rk, 0.0_rk, kind=rk)
    cartesian(1, 3) = cmplx(1.0_rk, 0.0_rk, kind=rk)
    call cartesian_to_axial_harmonics(cartesian, harmonic)
    call require_close(harmonic(projection_plus, projection_zero), &
                       cmplx(inverse_sqrt_two, 0.0_rk, kind=rk), &
                       "Axz to A+0 factor differs from new_src")
    call require_close(harmonic(projection_minus, projection_zero), &
                       cmplx(inverse_sqrt_two, 0.0_rk, kind=rk), &
                       "Axz to A-0 factor differs from new_src")

    cartesian = cmplx(0.0_rk, 0.0_rk, kind=rk)
    cartesian(3, 3) = cmplx(2.0_rk, -0.75_rk, kind=rk)
    call cartesian_to_axial_harmonics(cartesian, harmonic)
    call require_close(harmonic(projection_zero, projection_zero), &
                       cartesian(3, 3), "Azz to A00 formula differs from new_src")
  end subroutine test_cartesian_component_formulas


  subroutine test_complex_round_trip()
    complex(rk) :: cartesian(3, 3), harmonic(3, 3), restored(3, 3)
    integer :: orbital, spin

    do spin = 1, 3
      do orbital = 1, 3
        cartesian(spin, orbital) = cmplx( &
          0.31_rk * real(4 * spin - orbital, rk), &
          -0.17_rk * real(spin + 3 * orbital, rk), kind=rk)
      end do
    end do

    call cartesian_to_axial_harmonics(cartesian, harmonic)
    call axial_harmonics_to_cartesian(harmonic, restored)
    call require(maxval(abs(restored - cartesian)) < 32.0_rk * epsilon(1.0_rk), &
                 "Cartesian/harmonic complex round-trip failed")
  end subroutine test_complex_round_trip


  subroutine test_radial_transport_projection()
    complex(rk) :: cartesian(3, 3), harmonic(3, 3), transported(3, 3)
    integer :: orbital, spin

    do spin = 1, 3
      do orbital = 1, 3
        cartesian(spin, orbital) = cmplx( &
          0.21_rk * real(5 * spin + orbital, rk), &
          0.13_rk * real(2 * spin - 3 * orbital, rk), kind=rk)
      end do
    end do

    call cartesian_to_axial_harmonics(cartesian, harmonic)
    call axial_harmonics_to_transport_cartesian(harmonic, transported)

    call require_close(transported(1, 2), -cartesian(2, 1), &
                       "radial transport reconstruction changed its Axy=-Ayx rule")
    transported(1, 2) = cartesian(1, 2)
    call require(maxval(abs(transported - cartesian)) < &
                 64.0_rk * epsilon(1.0_rk), &
                 "radial transport reconstruction changed outside Axy")
  end subroutine test_radial_transport_projection


  subroutine test_axial_phase_factors()
    complex(rk) :: radial(3, 3), phased(3, 3)
    complex(rk), parameter :: amplitude = cmplx(0.7_rk, -0.2_rk, kind=rk)
    real(rk), parameter :: angle = 0.37_rk, winding = 1.0_rk

    radial = cmplx(0.0_rk, 0.0_rk, kind=rk)
    radial(projection_plus, projection_plus) = amplitude
    call phase_axial_harmonics(radial, angle, winding, phased)
    call require_close(phased(projection_plus, projection_plus), &
                       amplitude * exp(cmplx(0.0_rk, -angle, kind=rk)), &
                       "A++ axial phase exponent is wrong")

    radial = cmplx(0.0_rk, 0.0_rk, kind=rk)
    radial(projection_minus, projection_minus) = amplitude
    call phase_axial_harmonics(radial, angle, winding, phased)
    call require_close(phased(projection_minus, projection_minus), &
                       amplitude * exp(cmplx(0.0_rk, 3.0_rk * angle, kind=rk)), &
                       "A-- axial phase exponent is wrong")

    radial = cmplx(0.0_rk, 0.0_rk, kind=rk)
    radial(projection_plus, projection_minus) = amplitude
    call phase_axial_harmonics(radial, angle, winding, phased)
    call require_close(phased(projection_plus, projection_minus), &
                       amplitude * exp(cmplx(0.0_rk, angle, kind=rk)), &
                       "A+- axial phase exponent is wrong")
  end subroutine test_axial_phase_factors


  subroutine test_b_phase_far_field()
    complex(rk) :: cartesian_axis(3, 3), harmonic(3, 3), reconstructed(3, 3)
    real(rk), parameter :: angle = 0.63_rk, winding = 1.0_rk
    integer :: component

    cartesian_axis = cmplx(0.0_rk, 0.0_rk, kind=rk)
    do component = 1, 3
      cartesian_axis(component, component) = cmplx(1.0_rk, 0.0_rk, kind=rk)
    end do

    call cartesian_to_axial_harmonics(cartesian_axis, harmonic)
    call reconstruct_axial_cartesian(harmonic, angle, winding, reconstructed)
    call require(maxval(abs(reconstructed - &
                 exp(cmplx(0.0_rk, angle, kind=rk)) * cartesian_axis)) < &
                 64.0_rk * epsilon(1.0_rk), &
                 "axial reconstruction does not produce exp(i phi) times identity")
  end subroutine test_b_phase_far_field


  subroutine require_close(actual, expected, message)
    complex(rk), intent(in) :: actual, expected
    character(len=*), intent(in) :: message

    call require(abs(actual - expected) < 64.0_rk * epsilon(1.0_rk), message)
  end subroutine require_close


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_order_parameter_basis
