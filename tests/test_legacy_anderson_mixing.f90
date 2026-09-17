program test_legacy_anderson_mixing
  use he3_kinds, only : rk
  use legacy_anderson_mixing, only : legacy_anderson_t, anderson_report_t
  implicit none

  call test_first_relaxation_step()
  call test_two_point_subspace_step()
  call test_against_original_citerat_trace()
  call test_convergence_on_coupled_linear_map()
  call test_history_limit_discards_a_vector()
  call test_singular_subspace_fallback()
  call test_reset()

  print '(a)', "legacy Anderson-mixing tests passed"

contains

  subroutine test_first_relaxation_step()
    type(legacy_anderson_t) :: accelerator
    type(anderson_report_t) :: report
    real(rk) :: current(3), mapped(3), next(3)

    current = [1.0_rk, -2.0_rk, 0.5_rk]
    mapped = [2.0_rk, 1.0_rk, -0.5_rk]

    call accelerator%initialize(3, 5, 0.1_rk, 10.0_rk)
    call accelerator%update(current, mapped, 1.0e-12_rk, next, report)

    call require(maxval(abs(next - (current + 0.01_rk * (mapped - current)))) < &
                 16.0_rk * epsilon(1.0_rk), &
                 "first Anderson update does not match legacy relaxation")
    call require(report%history_size == 1, &
                 "first Anderson update retained the wrong history size")
    call require(abs(report%max_residual - 3.0_rk) < epsilon(1.0_rk), &
                 "legacy maximum residual was not reported")
    call require(.not. report%converged, &
                 "nonzero first residual was reported as converged")
  end subroutine test_first_relaxation_step


  subroutine test_two_point_subspace_step()
    type(legacy_anderson_t) :: accelerator
    type(anderson_report_t) :: report
    real(rk) :: current(2), mapped(2), next(2)
    real(rk) :: coefficient, displacement(2), expected(2), minimum_residual(2)
    real(rk) :: f0(2), f1(2), mixing, u1

    call accelerator%initialize(2, 5, 0.1_rk, 10.0_rk)

    current = [0.0_rk, 0.0_rk]
    call linear_map(current, mapped)
    f0 = mapped - current
    call accelerator%update(current, mapped, 1.0e-14_rk, next, report)

    current = next
    call linear_map(current, mapped)
    f1 = mapped - current

    u1 = dot_product(f1, f1)
    coefficient = -dot_product(f0 - f1, f1) / &
                  dot_product(f0 - f1, f0 - f1)
    displacement = coefficient * ([0.0_rk, 0.0_rk] - current)
    minimum_residual = f1 + coefficient * (f0 - f1)
    mixing = 0.5_rk * &
             (sqrt(dot_product(displacement, displacement) / u1) + 0.01_rk)
    expected = current + displacement + mixing * minimum_residual

    call accelerator%update(current, mapped, 1.0e-14_rk, next, report)

    call require(maxval(abs(next - expected)) < 2.0e-12_rk, &
                 "two-point Anderson step differs from the citerat equations")
    call require(abs(report%mixing - mixing) < 2.0e-13_rk, &
                 "adaptive mixing differs from the citerat p1 update")
    call require(.not. report%used_stochastic_step, &
                 "regular two-point subspace unexpectedly used fallback")
  end subroutine test_two_point_subspace_step


  subroutine test_against_original_citerat_trace()
    type(legacy_anderson_t) :: accelerator
    type(anderson_report_t) :: report
    real(rk) :: current(3), mapped(3), next(3)
    real(rk), parameter :: expected_after_four(3) = [ &
      1.0006160362342167_rk, -1.9993673406499604_rk, &
      0.5002819786829410_rk]
    real(rk), parameter :: expected_fourth_residual = &
      6.111195348963e-4_rk
    integer :: iteration

    ! These values were produced by compiling and running the unchanged
    ! fixed-form citerat, excluse, rearrange, deps1, p1 and LINPACK routines
    ! on this same map with mmax=5, prog=0.1 and pmax=2.0.
    current = [3.0_rk, -4.0_rk, 2.0_rk]
    call accelerator%initialize(3, 5, 0.1_rk, 2.0_rk)
    do iteration = 1, 4
      call coupled_map(current, mapped)
      call accelerator%update(current, mapped, 0.0_rk, next, report)
      current = next
    end do

    call require(maxval(abs(current - expected_after_four)) < 2.0e-9_rk, &
                 "modern Anderson trace differs from original citerat")
    call require(abs(report%max_residual - expected_fourth_residual) < &
                 2.0e-9_rk, &
                 "modern residual trace differs from original citerat")
  end subroutine test_against_original_citerat_trace


  subroutine test_convergence_on_coupled_linear_map()
    type(legacy_anderson_t) :: accelerator
    type(anderson_report_t) :: report
    real(rk) :: current(3), mapped(3), next(3), solution(3)
    integer :: iteration

    current = [3.0_rk, -4.0_rk, 2.0_rk]
    solution = [1.0_rk, -2.0_rk, 0.5_rk]
    call accelerator%initialize(3, 5, 0.1_rk, 2.0_rk)

    do iteration = 1, 100
      call coupled_map(current, mapped)
      call accelerator%update(current, mapped, 1.0e-11_rk, next, report)
      current = next
      if (report%converged) exit
    end do

    call require(report%converged, &
                 "Anderson accelerator did not converge on a contractive map")
    call require(maxval(abs(current - solution)) < 1.0e-9_rk, &
                 "Anderson accelerator converged to the wrong fixed point")
    call require(iteration < 100, &
                 "Anderson accelerator exhausted the test iteration limit")
  end subroutine test_convergence_on_coupled_linear_map


  subroutine test_history_limit_discards_a_vector()
    type(legacy_anderson_t) :: accelerator
    type(anderson_report_t) :: report
    real(rk) :: current(3), mapped(3), next(3)

    call accelerator%initialize(3, 2, 0.1_rk, 2.0_rk)

    current = 0.0_rk
    mapped = current + [1.0_rk, 0.0_rk, 0.0_rk]
    call accelerator%update(current, mapped, 0.0_rk, next, report)

    current = next
    mapped = current + [0.0_rk, 1.0_rk, 0.0_rk]
    call accelerator%update(current, mapped, 0.0_rk, next, report)

    current = next
    mapped = current + [0.0_rk, 0.0_rk, 1.0_rk]
    call accelerator%update(current, mapped, 0.0_rk, next, report)

    call require(report%discarded_vectors >= 1, &
                 "full Anderson history did not discard a vector")
    call require(report%history_size <= 2, &
                 "Anderson history exceeded its configured limit")
    call require(all(abs(next) < huge(1.0_rk)), &
                 "history-limited Anderson update produced a non-finite value")
  end subroutine test_history_limit_discards_a_vector


  subroutine test_singular_subspace_fallback()
    type(legacy_anderson_t) :: accelerator
    type(anderson_report_t) :: report
    real(rk) :: current(2), mapped(2), next(2)

    call accelerator%initialize(2, 3, 0.1_rk, 2.0_rk)

    current = 0.0_rk
    mapped = current + [1.0_rk, -1.0_rk]
    call accelerator%update(current, mapped, 0.0_rk, next, report)

    current = next
    mapped = current + [1.0_rk, -1.0_rk]
    call accelerator%update(current, mapped, 0.0_rk, next, report)

    call require(report%used_stochastic_step, &
                 "singular Anderson subspace did not use the legacy fallback")
    call require(abs(report%mixing - 1.0_rk) < epsilon(1.0_rk), &
                 "legacy singular-subspace fallback did not reset mixing")
    call require(all(abs(next) < huge(1.0_rk)), &
                 "singular-subspace fallback produced a non-finite value")
  end subroutine test_singular_subspace_fallback


  subroutine test_reset()
    type(legacy_anderson_t) :: accelerator
    type(anderson_report_t) :: report
    real(rk) :: current(2), mapped(2), next(2), expected(2)

    current = [0.25_rk, -0.75_rk]
    mapped = [1.0_rk, 0.5_rk]
    expected = current + 0.01_rk * (mapped - current)

    call accelerator%initialize(2, 3)
    call accelerator%update(current, mapped, 0.0_rk, next, report)
    call accelerator%reset()
    call accelerator%update(current, mapped, 0.0_rk, next, report)

    call require(maxval(abs(next - expected)) < 16.0_rk * epsilon(1.0_rk), &
                 "reset did not restore the initial legacy mixing step")
    call require(report%iteration == 1 .and. report%history_size == 1, &
                 "reset did not clear Anderson iteration history")
  end subroutine test_reset


  pure subroutine linear_map(x, mapped)
    real(rk), intent(in) :: x(2)
    real(rk), intent(out) :: mapped(2)

    mapped(1) = 0.5_rk * x(1) + 1.0_rk
    mapped(2) = 0.25_rk * x(2) - 2.0_rk
  end subroutine linear_map


  pure subroutine coupled_map(x, mapped)
    real(rk), intent(in) :: x(3)
    real(rk), intent(out) :: mapped(3)
    real(rk), parameter :: fixed_point(3) = [1.0_rk, -2.0_rk, 0.5_rk]
    real(rk) :: difference(3)

    difference = x - fixed_point
    mapped(1) = fixed_point(1) + 0.35_rk * difference(1) + &
                0.10_rk * difference(2)
    mapped(2) = fixed_point(2) + 0.05_rk * difference(1) + &
                0.30_rk * difference(2) + 0.04_rk * difference(3)
    mapped(3) = fixed_point(3) + 0.08_rk * difference(2) + &
                0.25_rk * difference(3)
  end subroutine coupled_map


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_legacy_anderson_mixing
