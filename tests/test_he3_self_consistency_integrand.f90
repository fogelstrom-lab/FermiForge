program test_he3_self_consistency_integrand
  use he3_kinds, only : rk
  use quasiclassical_propagator, only : quasiclassical_propagator_t, &
                                        reconstruct_legacy_quasiclassical_propagator
  use he3_self_consistency_integrand, only : he3_point_map_accumulator_t, &
                                             accumulate_he3_point_map_contribution, &
                                             finalize_he3_point_map
  implicit none

  call test_frozen_legacy_green_contribution()
  call test_accumulation_and_finalization()

  print '(a)', "3He self-consistency integrand tests passed"

contains

  subroutine test_frozen_legacy_green_contribution()
    type(he3_point_map_accumulator_t) :: accumulator
    type(quasiclassical_propagator_t) :: propagator
    complex(rk) :: expected_gap(3, 3), forward(0:3), reverse(0:3)
    real(rk) :: expected_current(3), momentum(3), transverse_momentum
    logical :: succeeded

    call frozen_coherences(forward, reverse)
    call reconstruct_legacy_quasiclassical_propagator( &
      forward, reverse, propagator, succeeded)
    call require(succeeded, "frozen legacy coherence matrices are singular")

    transverse_momentum = sqrt(1.0_rk - 0.35_rk**2)
    momentum = [cos(0.41_rk) * transverse_momentum, &
                sin(0.41_rk) * transverse_momentum, 0.35_rk]
    call accumulate_he3_point_map_contribution( &
      accumulator, propagator, momentum, 0.11_rk * 0.17_rk, 0.42_rk, 0.37_rk)

    expected_gap(1, :) = [ &
      cmplx(4.92682998500639e-3_rk, 8.89835106135951e-4_rk, kind=rk), &
      cmplx(2.14135405124218e-3_rk, 3.86750104075139e-4_rk, kind=rk), &
      cmplx(2.00717622290815e-3_rk, 3.62516237170849e-4_rk, kind=rk)]
    expected_gap(2, :) = [ &
      cmplx(-1.89752698425036e-3_rk, 3.56259996796113e-3_rk, kind=rk), &
      cmplx(-8.24724438925529e-4_rk, 1.54841711558247e-3_rk, kind=rk), &
      cmplx(-7.73046980858824e-4_rk, 1.45139287720228e-3_rk, kind=rk)]
    expected_gap(3, :) = [ &
      cmplx(2.50020067590406e-3_rk, -2.71553398666873e-3_rk, kind=rk), &
      cmplx(1.08666523150964e-3_rk, -1.18025580775777e-3_rk, kind=rk), &
      cmplx(1.01857449200513e-3_rk, -1.10630065724368e-3_rk, kind=rk)]
    expected_current = [ &
      -1.95017605137306e-5_rk, &
      -8.47607366390128e-6_rk, &
      -7.94496057853277e-6_rk]

    call require(maxval(abs(accumulator%gap - expected_gap)) < 5.0e-17_rk, &
                 "modern gap integrand differs from new_src:green")
    call require(maxval(abs(accumulator%current_mean_field - &
                 expected_current)) < 5.0e-19_rk, &
                 "modern current integrand differs from new_src:green")
    call require(accumulator%contribution_count == 1, &
                 "one integrand evaluation has the wrong contribution count")
  end subroutine test_frozen_legacy_green_contribution


  subroutine test_accumulation_and_finalization()
    type(he3_point_map_accumulator_t) :: accumulator
    type(quasiclassical_propagator_t) :: propagator
    complex(rk) :: first_gap(3, 3), forward(0:3), mapped_gap(3, 3), reverse(0:3)
    real(rk) :: first_current(3), mapped_current(3), momentum(3)
    logical :: succeeded

    call frozen_coherences(forward, reverse)
    call reconstruct_legacy_quasiclassical_propagator( &
      forward, reverse, propagator, succeeded)
    call require(succeeded, "accumulation coherence matrices are singular")
    momentum = [0.36_rk, 0.48_rk, 0.8_rk]

    call accumulate_he3_point_map_contribution( &
      accumulator, propagator, momentum, 0.07_rk, 0.31_rk, 0.42_rk)
    first_gap = accumulator%gap
    first_current = accumulator%current_mean_field
    call accumulate_he3_point_map_contribution( &
      accumulator, propagator, momentum, 0.07_rk, 0.31_rk, 0.42_rk)
    call require(maxval(abs(accumulator%gap - &
                 cmplx(2.0_rk, 0.0_rk, kind=rk) * first_gap)) < 2.0e-17_rk, &
                 "gap contributions do not accumulate")
    call require(maxval(abs(accumulator%current_mean_field - &
                 2.0_rk * first_current)) < 2.0e-19_rk, &
                 "current contributions do not accumulate")

    call finalize_he3_point_map(accumulator, 0.73_rk, mapped_gap, mapped_current)
    call require(maxval(abs(mapped_gap - &
                 cmplx(0.73_rk, 0.0_rk, kind=rk) * accumulator%gap)) < &
                 tiny(1.0_rk), "gap prefactor was not applied at finalization")
    call require(maxval(abs(mapped_current - &
                 accumulator%current_mean_field)) < tiny(1.0_rk), &
                 "gap prefactor incorrectly changed the current field")

    call accumulator%reset()
    call require(maxval(abs(accumulator%gap)) < tiny(1.0_rk) .and. &
                 maxval(abs(accumulator%current_mean_field)) < tiny(1.0_rk) .and. &
                 accumulator%contribution_count == 0, &
                 "point-map accumulator reset failed")
  end subroutine test_accumulation_and_finalization


  subroutine frozen_coherences(forward, reverse)
    complex(rk), intent(out) :: forward(0:3), reverse(0:3)

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
  end subroutine frozen_coherences


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_he3_self_consistency_integrand
