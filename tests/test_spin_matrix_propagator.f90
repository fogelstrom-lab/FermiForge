program test_spin_matrix_propagator
  use he3_kinds, only : rk
  use spin_matrix_2x2, only : spin_identity_matrix, &
                              pauli_coefficients_to_matrix, &
                              matrix_to_pauli_coefficients, &
                              pair_coefficients_to_matrix, &
                              legacy_reverse_coefficients_to_matrix, &
                              invert_spin_matrix, make_spin_rotation, &
                              rotate_pair_matrix
  use quasiclassical_propagator, only : quasiclassical_propagator_t, &
                                        reconstruct_quasiclassical_propagator, &
                                        reconstruct_legacy_quasiclassical_propagator
  implicit none

  call test_pauli_round_trip()
  call test_pair_matrix_convention()
  call test_legacy_reverse_convention()
  call test_spin_matrix_inverse()
  call test_pair_spin_rotation()
  call test_spin_scalar_limit()
  call test_propagator_normalization()
  call test_conjugate_coherence_relation()
  call test_frozen_trajectory_coherences()

  print '(a)', "spin-matrix propagator tests passed"

contains

  subroutine test_pauli_round_trip()
    complex(rk) :: coefficients(0:3), recovered(0:3)

    coefficients = [ &
      cmplx(0.17_rk, -0.09_rk, kind=rk), &
      cmplx(-0.31_rk, 0.22_rk, kind=rk), &
      cmplx(0.07_rk, 0.41_rk, kind=rk), &
      cmplx(0.28_rk, -0.13_rk, kind=rk)]
    recovered = matrix_to_pauli_coefficients( &
      pauli_coefficients_to_matrix(coefficients))
    call require(maxval(abs(recovered - coefficients)) < 2.0e-15_rk, &
                 "Pauli coefficient round trip failed")
  end subroutine test_pauli_round_trip


  subroutine test_pair_matrix_convention()
    complex(rk) :: coefficients(0:3), expected(2, 2), matrix(2, 2)

    coefficients = [ &
      cmplx(0.13_rk, -0.04_rk, kind=rk), &
      cmplx(0.21_rk, 0.07_rk, kind=rk), &
      cmplx(-0.18_rk, 0.11_rk, kind=rk), &
      cmplx(0.06_rk, -0.23_rk, kind=rk)]
    expected(1, 1) = -coefficients(1) + &
      cmplx(0.0_rk, 1.0_rk, kind=rk) * coefficients(2)
    expected(1, 2) = coefficients(0) + coefficients(3)
    expected(2, 1) = -coefficients(0) + coefficients(3)
    expected(2, 2) = coefficients(1) + &
      cmplx(0.0_rk, 1.0_rk, kind=rk) * coefficients(2)
    matrix = pair_coefficients_to_matrix(coefficients)
    call require(maxval(abs(matrix - expected)) < tiny(1.0_rk), &
                 "triplet/singlet pair-matrix convention changed")
  end subroutine test_pair_matrix_convention


  subroutine test_legacy_reverse_convention()
    complex(rk) :: coefficients(0:3), expected(2, 2), matrix(2, 2)

    coefficients = [ &
      cmplx(-0.08_rk, 0.03_rk, kind=rk), &
      cmplx(0.19_rk, -0.12_rk, kind=rk), &
      cmplx(0.24_rk, 0.05_rk, kind=rk), &
      cmplx(-0.09_rk, 0.16_rk, kind=rk)]
    expected(1, 1) = -coefficients(1) - &
      cmplx(0.0_rk, 1.0_rk, kind=rk) * coefficients(2)
    expected(1, 2) = coefficients(0) + coefficients(3)
    expected(2, 1) = -coefficients(0) + coefficients(3)
    expected(2, 2) = coefficients(1) - &
      cmplx(0.0_rk, 1.0_rk, kind=rk) * coefficients(2)
    matrix = legacy_reverse_coefficients_to_matrix(coefficients)
    call require(maxval(abs(matrix - expected)) < tiny(1.0_rk), &
                 "legacy reverse coherence convention changed")
  end subroutine test_legacy_reverse_convention


  subroutine test_spin_matrix_inverse()
    complex(rk) :: determinant, inverse(2, 2), matrix(2, 2)
    logical :: succeeded

    matrix(1, :) = [cmplx(1.2_rk, 0.1_rk, kind=rk), &
                    cmplx(-0.2_rk, 0.3_rk, kind=rk)]
    matrix(2, :) = [cmplx(0.4_rk, -0.1_rk, kind=rk), &
                    cmplx(0.8_rk, 0.2_rk, kind=rk)]
    call invert_spin_matrix(matrix, inverse, determinant, succeeded)
    call require(succeeded, "well-conditioned spin matrix was called singular")
    call require(maxval(abs(matmul(matrix, inverse) - &
                 spin_identity_matrix())) < 2.0e-15_rk, &
                 "spin-matrix inverse is inaccurate")
  end subroutine test_spin_matrix_inverse


  subroutine test_pair_spin_rotation()
    real(rk), parameter :: angle = 0.63_rk
    complex(rk) :: coefficients(0:3), rotated_coefficients(0:3)
    complex(rk) :: rotation(2, 2), rotated_from_matrix(2, 2)

    coefficients = [ &
      cmplx(0.05_rk, -0.02_rk, kind=rk), &
      cmplx(0.27_rk, 0.08_rk, kind=rk), &
      cmplx(-0.14_rk, 0.19_rk, kind=rk), &
      cmplx(0.11_rk, -0.07_rk, kind=rk)]
    rotated_coefficients(0) = coefficients(0)
    rotated_coefficients(1) = &
      cmplx(cos(angle), 0.0_rk, kind=rk) * coefficients(1) - &
      cmplx(sin(angle), 0.0_rk, kind=rk) * coefficients(2)
    rotated_coefficients(2) = &
      cmplx(sin(angle), 0.0_rk, kind=rk) * coefficients(1) + &
      cmplx(cos(angle), 0.0_rk, kind=rk) * coefficients(2)
    rotated_coefficients(3) = coefficients(3)

    call make_spin_rotation([0.0_rk, 0.0_rk, 1.0_rk], angle, rotation)
    rotated_from_matrix = rotate_pair_matrix( &
      pair_coefficients_to_matrix(coefficients), rotation)
    call require(maxval(abs(rotated_from_matrix - &
                 pair_coefficients_to_matrix(rotated_coefficients))) < &
                 3.0e-15_rk, &
                 "triplet pair matrix is not covariant under spin rotation")
  end subroutine test_pair_spin_rotation


  subroutine test_spin_scalar_limit()
    type(quasiclassical_propagator_t) :: propagator
    complex(rk) :: denominator, expected_anomalous, expected_normal
    complex(rk) :: gamma(2, 2), gamma_tilde(2, 2)
    complex(rk), parameter :: amplitude = &
      cmplx(0.13_rk, -0.07_rk, kind=rk)
    complex(rk), parameter :: tilde_amplitude = &
      cmplx(-0.04_rk, 0.09_rk, kind=rk)
    complex(rk), parameter :: imaginary_unit = &
      cmplx(0.0_rk, 1.0_rk, kind=rk)
    logical :: succeeded

    gamma = amplitude * spin_identity_matrix()
    gamma_tilde = tilde_amplitude * spin_identity_matrix()
    call reconstruct_quasiclassical_propagator( &
      gamma, gamma_tilde, propagator, succeeded)
    call require(succeeded, "spin-scalar coherence limit is singular")

    denominator = cmplx(1.0_rk, 0.0_rk, kind=rk) - &
                  amplitude * tilde_amplitude
    expected_normal = -imaginary_unit * &
      (cmplx(1.0_rk, 0.0_rk, kind=rk) + amplitude * tilde_amplitude) / &
      denominator
    expected_anomalous = -cmplx(2.0_rk, 0.0_rk, kind=rk) * &
      imaginary_unit * amplitude / denominator
    call require(maxval(abs(propagator%particle - &
                 expected_normal * spin_identity_matrix())) < 2.0e-15_rk, &
                 "spin-scalar normal propagator has the wrong limit")
    call require(maxval(abs(propagator%anomalous - &
                 expected_anomalous * spin_identity_matrix())) < 2.0e-15_rk, &
                 "spin-scalar anomalous propagator has the wrong limit")
  end subroutine test_spin_scalar_limit


  subroutine test_propagator_normalization()
    type(quasiclassical_propagator_t) :: propagator
    complex(rk) :: forward(0:3), reverse(0:3)
    logical :: succeeded

    forward = [ &
      cmplx(0.04_rk, -0.03_rk, kind=rk), &
      cmplx(-0.12_rk, 0.09_rk, kind=rk), &
      cmplx(0.07_rk, 0.05_rk, kind=rk), &
      cmplx(0.11_rk, -0.02_rk, kind=rk)]
    reverse = [ &
      cmplx(-0.01_rk, 0.02_rk, kind=rk), &
      cmplx(0.08_rk, 0.06_rk, kind=rk), &
      cmplx(-0.09_rk, 0.03_rk, kind=rk), &
      cmplx(0.05_rk, -0.07_rk, kind=rk)]
    call reconstruct_legacy_quasiclassical_propagator( &
      forward, reverse, propagator, succeeded)
    call require(succeeded, "regular coherence matrices were called singular")
    call require(propagator%normalization_error() < 8.0e-15_rk, &
                 "source propagator does not satisfy g squared equals minus one")
  end subroutine test_propagator_normalization


  subroutine test_conjugate_coherence_relation()
    type(quasiclassical_propagator_t) :: propagator
    complex(rk) :: gamma(2, 2), gamma_tilde(2, 2)
    logical :: succeeded

    gamma(1, :) = [cmplx(0.07_rk, -0.02_rk, kind=rk), &
                    cmplx(-0.11_rk, 0.04_rk, kind=rk)]
    gamma(2, :) = [cmplx(0.03_rk, 0.08_rk, kind=rk), &
                    cmplx(0.09_rk, -0.05_rk, kind=rk)]
    gamma_tilde = transpose(conjg(gamma))
    call reconstruct_quasiclassical_propagator( &
      gamma, gamma_tilde, propagator, succeeded)
    call require(succeeded, "conjugate coherence pair is singular")
    call require(maxval(abs(propagator%tilde_anomalous - &
                 transpose(conjg(propagator%anomalous)))) < 4.0e-15_rk, &
                 "conjugate coherence relation is not preserved in f tilde")
  end subroutine test_conjugate_coherence_relation


  subroutine test_frozen_trajectory_coherences()
    type(quasiclassical_propagator_t) :: from_legacy, from_matrices
    complex(rk) :: forward(0:3), reverse(0:3)
    logical :: legacy_succeeded, matrix_succeeded

    forward = cmplx(0.0_rk, 0.0_rk, kind=rk)
    forward(1:3) = [ &
      cmplx(-3.45416635611957e-2_rk, 1.17489705712979e-1_rk, kind=rk), &
      cmplx(-9.01737612952456e-2_rk, -5.58797007439699e-2_rk, kind=rk), &
      cmplx(8.02906423983453e-2_rk, 7.10955985562741e-2_rk, kind=rk)]
    reverse = cmplx(0.0_rk, 0.0_rk, kind=rk)
    reverse(1:3) = [ &
      cmplx(1.62659859791399e-2_rk, 1.42678313138293e-1_rk, kind=rk), &
      cmplx(9.95709656830779e-2_rk, -4.16922195641589e-2_rk, kind=rk), &
      cmplx(-6.37043298757304e-2_rk, 5.90495834745839e-2_rk, kind=rk)]

    call reconstruct_legacy_quasiclassical_propagator( &
      forward, reverse, from_legacy, legacy_succeeded)
    call reconstruct_quasiclassical_propagator( &
      pair_coefficients_to_matrix(forward), &
      legacy_reverse_coefficients_to_matrix(reverse), &
      from_matrices, matrix_succeeded)
    call require(legacy_succeeded .and. matrix_succeeded, &
                 "frozen trajectory coherence matrices are singular")
    call require(maxval(abs(from_legacy%nambu_matrix() - &
                 from_matrices%nambu_matrix())) < tiny(1.0_rk), &
                 "legacy coefficient and explicit matrix paths differ")
    call require(from_legacy%normalization_error() < 8.0e-15_rk, &
                 "frozen trajectory propagator fails normalization")
  end subroutine test_frozen_trajectory_coherences


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_spin_matrix_propagator
