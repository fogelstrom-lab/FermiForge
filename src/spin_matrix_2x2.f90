module spin_matrix_2x2
  use he3_kinds, only : rk
  implicit none
  private

  complex(rk), parameter :: imaginary_unit = &
    cmplx(0.0_rk, 1.0_rk, kind=rk)
  complex(rk), parameter :: complex_half = &
    cmplx(0.5_rk, 0.0_rk, kind=rk)

  public :: spin_identity_matrix
  public :: pauli_matrix
  public :: pauli_coefficients_to_matrix
  public :: matrix_to_pauli_coefficients
  public :: pair_coefficients_to_matrix
  public :: legacy_reverse_coefficients_to_matrix
  public :: invert_spin_matrix
  public :: make_spin_rotation
  public :: rotate_pair_matrix

contains

  pure function spin_identity_matrix() result(identity)
    complex(rk) :: identity(2, 2)

    identity = cmplx(0.0_rk, 0.0_rk, kind=rk)
    identity(1, 1) = cmplx(1.0_rk, 0.0_rk, kind=rk)
    identity(2, 2) = cmplx(1.0_rk, 0.0_rk, kind=rk)
  end function spin_identity_matrix


  pure function pauli_matrix(component) result(matrix)
    integer, intent(in) :: component
    complex(rk) :: matrix(2, 2)

    matrix = cmplx(0.0_rk, 0.0_rk, kind=rk)
    select case (component)
    case (1)
      matrix(1, 2) = cmplx(1.0_rk, 0.0_rk, kind=rk)
      matrix(2, 1) = cmplx(1.0_rk, 0.0_rk, kind=rk)
    case (2)
      matrix(1, 2) = -imaginary_unit
      matrix(2, 1) = imaginary_unit
    case (3)
      matrix(1, 1) = cmplx(1.0_rk, 0.0_rk, kind=rk)
      matrix(2, 2) = cmplx(-1.0_rk, 0.0_rk, kind=rk)
    case default
      error stop "Pauli-matrix component must be 1, 2, or 3"
    end select
  end function pauli_matrix


  pure function pauli_coefficients_to_matrix(coefficients) result(matrix)
    complex(rk), intent(in) :: coefficients(0:3)
    complex(rk) :: matrix(2, 2)

    matrix(1, 1) = coefficients(0) + coefficients(3)
    matrix(1, 2) = coefficients(1) - imaginary_unit * coefficients(2)
    matrix(2, 1) = coefficients(1) + imaginary_unit * coefficients(2)
    matrix(2, 2) = coefficients(0) - coefficients(3)
  end function pauli_coefficients_to_matrix


  pure function matrix_to_pauli_coefficients(matrix) result(coefficients)
    complex(rk), intent(in) :: matrix(2, 2)
    complex(rk) :: coefficients(0:3)

    coefficients(0) = complex_half * (matrix(1, 1) + matrix(2, 2))
    coefficients(1) = complex_half * (matrix(1, 2) + matrix(2, 1))
    coefficients(2) = -complex_half * imaginary_unit * &
      (matrix(2, 1) - matrix(1, 2))
    coefficients(3) = complex_half * (matrix(1, 1) - matrix(2, 2))
  end function matrix_to_pauli_coefficients


  pure function pair_coefficients_to_matrix(coefficients) result(matrix)
    complex(rk), intent(in) :: coefficients(0:3)
    complex(rk) :: matrix(2, 2)

    ! This is (coefficient_0 I + coefficient_alpha sigma_alpha) i sigma_y.
    ! coefficient_0 is the singlet amplitude and coefficients 1:3 are the
    ! triplet d vector.  It is also exactly the Gg construction in
    ! new_src/riccati.f90:green.
    matrix(1, 1) = -coefficients(1) + imaginary_unit * coefficients(2)
    matrix(1, 2) = coefficients(0) + coefficients(3)
    matrix(2, 1) = -coefficients(0) + coefficients(3)
    matrix(2, 2) = coefficients(1) + imaginary_unit * coefficients(2)
  end function pair_coefficients_to_matrix


  pure function legacy_reverse_coefficients_to_matrix(coefficients) result(matrix)
    complex(rk), intent(in) :: coefficients(0:3)
    complex(rk) :: matrix(2, 2)

    ! Exact transcription of the Gt construction in new_src/riccati.f90.
    ! This name deliberately does not assert that local conjugation is the
    ! project's eventual general equilibrium tilde operation.
    matrix(1, 1) = -coefficients(1) - imaginary_unit * coefficients(2)
    matrix(1, 2) = coefficients(0) + coefficients(3)
    matrix(2, 1) = -coefficients(0) + coefficients(3)
    matrix(2, 2) = coefficients(1) - imaginary_unit * coefficients(2)
  end function legacy_reverse_coefficients_to_matrix


  pure subroutine invert_spin_matrix(matrix, inverse, determinant, succeeded)
    complex(rk), intent(in) :: matrix(2, 2)
    complex(rk), intent(out) :: inverse(2, 2), determinant
    logical, intent(out) :: succeeded

    real(rk) :: scale

    determinant = matrix(1, 1) * matrix(2, 2) - &
                  matrix(1, 2) * matrix(2, 1)
    scale = max(1.0_rk, maxval(abs(matrix)))
    succeeded = abs(determinant) > &
      32.0_rk * epsilon(1.0_rk) * scale * scale
    if (.not. succeeded) then
      inverse = cmplx(0.0_rk, 0.0_rk, kind=rk)
      return
    end if

    inverse(1, 1) = matrix(2, 2) / determinant
    inverse(1, 2) = -matrix(1, 2) / determinant
    inverse(2, 1) = -matrix(2, 1) / determinant
    inverse(2, 2) = matrix(1, 1) / determinant
  end subroutine invert_spin_matrix


  subroutine make_spin_rotation(axis, angle, rotation)
    real(rk), intent(in) :: axis(3), angle
    complex(rk), intent(out) :: rotation(2, 2)

    complex(rk) :: generator(2, 2)
    real(rk) :: axis_norm
    integer :: component

    axis_norm = sqrt(dot_product(axis, axis))
    if (axis_norm <= 32.0_rk * epsilon(1.0_rk)) &
      error stop "spin-rotation axis cannot be zero"

    generator = cmplx(0.0_rk, 0.0_rk, kind=rk)
    do component = 1, 3
      generator = generator + &
        cmplx(axis(component) / axis_norm, 0.0_rk, kind=rk) * &
        pauli_matrix(component)
    end do
    rotation = cmplx(cos(0.5_rk * angle), 0.0_rk, kind=rk) * &
                 spin_identity_matrix() - &
               imaginary_unit * &
                 cmplx(sin(0.5_rk * angle), 0.0_rk, kind=rk) * generator
  end subroutine make_spin_rotation


  pure function rotate_pair_matrix(pair_matrix, rotation) result(rotated)
    complex(rk), intent(in) :: pair_matrix(2, 2), rotation(2, 2)
    complex(rk) :: rotated(2, 2)

    rotated = matmul(rotation, matmul(pair_matrix, transpose(rotation)))
  end function rotate_pair_matrix

end module spin_matrix_2x2
