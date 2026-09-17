module order_parameter_basis
  use he3_kinds, only : rk
  implicit none
  private

  integer, parameter, public :: projection_plus = 1
  integer, parameter, public :: projection_zero = 2
  integer, parameter, public :: projection_minus = 3
  integer, parameter :: projection_value(3) = [1, 0, -1]
  real(rk), parameter :: inverse_sqrt_two = 1.0_rk / sqrt(2.0_rk)
  complex(rk), parameter :: complex_half = cmplx(0.5_rk, 0.0_rk, kind=rk)
  complex(rk), parameter :: complex_inverse_sqrt_two = &
    cmplx(inverse_sqrt_two, 0.0_rk, kind=rk)
  complex(rk), parameter :: imaginary_unit = cmplx(0.0_rk, 1.0_rk, kind=rk)

  public :: cartesian_to_axial_harmonics
  public :: axial_harmonics_to_cartesian
  public :: axial_harmonics_to_transport_cartesian
  public :: phase_axial_harmonics
  public :: reconstruct_axial_cartesian

contains

  pure subroutine cartesian_to_axial_harmonics(cartesian, harmonic)
    complex(rk), intent(in) :: cartesian(3, 3)
    complex(rk), intent(out) :: harmonic(3, 3)

    complex(rk) :: dxx, dxy, dxz, dyx, dyy, dyz, dzx, dzy, dzz

    dxx = cartesian(1, 1)
    dxy = cartesian(1, 2)
    dxz = cartesian(1, 3)
    dyx = cartesian(2, 1)
    dyy = cartesian(2, 2)
    dyz = cartesian(2, 3)
    dzx = cartesian(3, 1)
    dzy = cartesian(3, 2)
    dzz = cartesian(3, 3)

    ! These are the authoritative new_src/init_calc.f90 formulas. Harmonic
    ! indices are ordered (+,0,-) for both spin and orbital projection.
    harmonic(projection_minus, projection_minus) = &
      complex_half * (dxx - dyy + imaginary_unit * (dxy + dyx))
    harmonic(projection_minus, projection_zero) = &
      complex_inverse_sqrt_two * (dxz + imaginary_unit * dyz)
    harmonic(projection_minus, projection_plus) = &
      complex_half * (dxx + dyy - imaginary_unit * (dxy - dyx))

    harmonic(projection_zero, projection_minus) = &
      complex_inverse_sqrt_two * (dzx + imaginary_unit * dzy)
    harmonic(projection_zero, projection_zero) = dzz
    harmonic(projection_zero, projection_plus) = &
      complex_inverse_sqrt_two * (dzx - imaginary_unit * dzy)

    harmonic(projection_plus, projection_minus) = &
      complex_half * (dxx + dyy + imaginary_unit * (dxy - dyx))
    harmonic(projection_plus, projection_zero) = &
      complex_inverse_sqrt_two * (dxz - imaginary_unit * dyz)
    harmonic(projection_plus, projection_plus) = &
      complex_half * (dxx - dyy - imaginary_unit * (dxy + dyx))
  end subroutine cartesian_to_axial_harmonics


  pure subroutine axial_harmonics_to_cartesian(harmonic, cartesian)
    complex(rk), intent(in) :: harmonic(3, 3)
    complex(rk), intent(out) :: cartesian(3, 3)

    complex(rk) :: cmm, cmp, cpm, cpp

    cmm = harmonic(projection_minus, projection_minus)
    cmp = harmonic(projection_minus, projection_plus)
    cpm = harmonic(projection_plus, projection_minus)
    cpp = harmonic(projection_plus, projection_plus)

    ! Begin with the radial transport reconstruction, then restore the general
    ! Axy component that the axial source assumes to be -Ayx.
    call axial_harmonics_to_transport_cartesian(harmonic, cartesian)
    cartesian(1, 2) = imaginary_unit * complex_half * &
      (cpp - cmm + cmp - cpm)
  end subroutine axial_harmonics_to_cartesian


  pure subroutine axial_harmonics_to_transport_cartesian(harmonic, cartesian)
    complex(rk), intent(in) :: harmonic(3, 3)
    complex(rk), intent(out) :: cartesian(3, 3)

    complex(rk) :: cmm, cmo, cmp, com, coo, cop, cpm, cpo, cpp

    cmm = harmonic(projection_minus, projection_minus)
    cmo = harmonic(projection_minus, projection_zero)
    cmp = harmonic(projection_minus, projection_plus)
    com = harmonic(projection_zero, projection_minus)
    coo = harmonic(projection_zero, projection_zero)
    cop = harmonic(projection_zero, projection_plus)
    cpm = harmonic(projection_plus, projection_minus)
    cpo = harmonic(projection_plus, projection_zero)
    cpp = harmonic(projection_plus, projection_plus)

    ! These formulas are copied exactly from new_src/interpol.f90. In
    ! particular, the Axy formula enforces Axy=-Ayx. It is source-faithful on
    ! the axial radial subspace but is not a general inverse of the transform.
    cartesian(1, 1) = complex_half * (cmp + cpm + cmm + cpp)
    cartesian(1, 2) = &
      imaginary_unit * complex_half * (cmm - cpp + cmp - cpm)
    cartesian(1, 3) = complex_inverse_sqrt_two * (cmo + cpo)

    cartesian(2, 1) = &
      imaginary_unit * complex_half * (cpp - cmm - cmp + cpm)
    cartesian(2, 2) = complex_half * (cmp + cpm - cmm - cpp)
    cartesian(2, 3) = imaginary_unit * complex_inverse_sqrt_two * (cpo - cmo)

    cartesian(3, 1) = complex_inverse_sqrt_two * (com + cop)
    cartesian(3, 2) = imaginary_unit * complex_inverse_sqrt_two * (cop - com)
    cartesian(3, 3) = coo
  end subroutine axial_harmonics_to_transport_cartesian


  pure subroutine phase_axial_harmonics(radial_harmonic, angle, winding, &
                                        phased_harmonic)
    complex(rk), intent(in) :: radial_harmonic(3, 3)
    real(rk), intent(in) :: angle, winding
    complex(rk), intent(out) :: phased_harmonic(3, 3)

    real(rk) :: angular_momentum
    integer :: orbital_projection, spin_projection

    do spin_projection = projection_plus, projection_minus
      do orbital_projection = projection_plus, projection_minus
        angular_momentum = winding - &
          real(projection_value(spin_projection), rk) - &
          real(projection_value(orbital_projection), rk)
        phased_harmonic(spin_projection, orbital_projection) = &
          radial_harmonic(spin_projection, orbital_projection) * &
          exp(cmplx(0.0_rk, angular_momentum * angle, kind=rk))
      end do
    end do
  end subroutine phase_axial_harmonics


  pure subroutine reconstruct_axial_cartesian(radial_harmonic, angle, winding, &
                                              cartesian)
    complex(rk), intent(in) :: radial_harmonic(3, 3)
    real(rk), intent(in) :: angle, winding
    complex(rk), intent(out) :: cartesian(3, 3)

    complex(rk) :: phased_harmonic(3, 3)

    call phase_axial_harmonics(radial_harmonic, angle, winding, phased_harmonic)
    call axial_harmonics_to_transport_cartesian(phased_harmonic, cartesian)
  end subroutine reconstruct_axial_cartesian

end module order_parameter_basis
