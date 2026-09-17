module quasiclassical_propagator
  use he3_kinds, only : rk
  use spin_matrix_2x2, only : spin_identity_matrix, invert_spin_matrix, &
                              pair_coefficients_to_matrix, &
                              legacy_reverse_coefficients_to_matrix
  implicit none
  private

  complex(rk), parameter :: imaginary_unit = &
    cmplx(0.0_rk, 1.0_rk, kind=rk)
  complex(rk), parameter :: complex_two = &
    cmplx(2.0_rk, 0.0_rk, kind=rk)

  type, public :: quasiclassical_propagator_t
    complex(rk) :: particle(2, 2) = cmplx(0.0_rk, 0.0_rk, kind=rk)
    complex(rk) :: anomalous(2, 2) = cmplx(0.0_rk, 0.0_rk, kind=rk)
    complex(rk) :: tilde_anomalous(2, 2) = cmplx(0.0_rk, 0.0_rk, kind=rk)
    complex(rk) :: hole(2, 2) = cmplx(0.0_rk, 0.0_rk, kind=rk)
  contains
    procedure :: nambu_matrix => propagator_nambu_matrix
    procedure :: normalization_error => propagator_normalization_error
  end type quasiclassical_propagator_t

  public :: reconstruct_quasiclassical_propagator
  public :: reconstruct_legacy_quasiclassical_propagator

contains

  pure subroutine reconstruct_quasiclassical_propagator( &
      gamma, gamma_tilde, propagator, succeeded)
    complex(rk), intent(in) :: gamma(2, 2), gamma_tilde(2, 2)
    type(quasiclassical_propagator_t), intent(out) :: propagator
    logical, intent(out) :: succeeded

    complex(rk) :: determinant
    complex(rk) :: gamma_gamma_tilde(2, 2), gamma_tilde_gamma(2, 2)
    complex(rk) :: normalization_plus(2, 2), normalization_minus(2, 2)
    complex(rk) :: identity(2, 2)
    logical :: plus_succeeded, minus_succeeded

    identity = spin_identity_matrix()
    gamma_gamma_tilde = matmul(gamma, gamma_tilde)
    gamma_tilde_gamma = matmul(gamma_tilde, gamma)
    call invert_spin_matrix(identity - gamma_gamma_tilde, &
                            normalization_plus, determinant, plus_succeeded)
    call invert_spin_matrix(identity - gamma_tilde_gamma, &
                            normalization_minus, determinant, minus_succeeded)
    succeeded = plus_succeeded .and. minus_succeeded
    if (.not. succeeded) then
      propagator = quasiclassical_propagator_t()
      return
    end if

    ! Source normalization from new_src/riccati.f90:green.  The resulting
    ! dimensionless Nambu propagator obeys g_hat**2 = -I for nonsingular
    ! coherence amplitudes.
    propagator%particle = -imaginary_unit * &
      matmul(normalization_plus, identity + gamma_gamma_tilde)
    propagator%anomalous = -complex_two * imaginary_unit * &
      matmul(normalization_plus, gamma)
    propagator%tilde_anomalous = complex_two * imaginary_unit * &
      matmul(normalization_minus, gamma_tilde)
    propagator%hole = imaginary_unit * &
      matmul(normalization_minus, identity + gamma_tilde_gamma)
  end subroutine reconstruct_quasiclassical_propagator


  pure subroutine reconstruct_legacy_quasiclassical_propagator( &
      forward_coefficients, reverse_coefficients, propagator, succeeded)
    complex(rk), intent(in) :: forward_coefficients(0:3)
    complex(rk), intent(in) :: reverse_coefficients(0:3)
    type(quasiclassical_propagator_t), intent(out) :: propagator
    logical, intent(out) :: succeeded

    complex(rk) :: gamma(2, 2), gamma_tilde(2, 2)

    gamma = pair_coefficients_to_matrix(forward_coefficients)
    gamma_tilde = &
      legacy_reverse_coefficients_to_matrix(reverse_coefficients)
    call reconstruct_quasiclassical_propagator( &
      gamma, gamma_tilde, propagator, succeeded)
  end subroutine reconstruct_legacy_quasiclassical_propagator


  pure function propagator_nambu_matrix(self) result(matrix)
    class(quasiclassical_propagator_t), intent(in) :: self
    complex(rk) :: matrix(4, 4)

    matrix(1:2, 1:2) = self%particle
    matrix(1:2, 3:4) = self%anomalous
    matrix(3:4, 1:2) = self%tilde_anomalous
    matrix(3:4, 3:4) = self%hole
  end function propagator_nambu_matrix


  pure real(rk) function propagator_normalization_error(self) result(error)
    class(quasiclassical_propagator_t), intent(in) :: self

    complex(rk) :: identity(4, 4), matrix(4, 4)
    integer :: component

    identity = cmplx(0.0_rk, 0.0_rk, kind=rk)
    do component = 1, 4
      identity(component, component) = cmplx(1.0_rk, 0.0_rk, kind=rk)
    end do
    matrix = self%nambu_matrix()
    error = maxval(abs(matmul(matrix, matrix) + identity))
  end function propagator_normalization_error

end module quasiclassical_propagator
