program test_riccati_segment
  use he3_kinds, only : rk
  use riccati_segment, only : advance_legacy_riccati_segment
  implicit none

  call test_frozen_legacy_checkpoint()
  call test_zero_distance()

  print '(a)', "Riccati segment tests passed"

contains

  subroutine test_frozen_legacy_checkpoint()
    complex(rk) :: coherence_in(0:3), coherence_out(0:3)
    complex(rk) :: pair_a_left(0:3), pair_a_right(0:3)
    complex(rk) :: pair_b_left(0:3), pair_b_right(0:3)
    complex(rk) :: exchange_left(3), exchange_right(3)
    complex(rk) :: expected(0:3)

    coherence_in = [cmplx(0.013_rk, -0.021_rk, kind=rk), &
                    cmplx(-0.14_rk, 0.032_rk, kind=rk), &
                    cmplx(0.071_rk, 0.055_rk, kind=rk), &
                    cmplx(-0.023_rk, -0.091_rk, kind=rk)]
    pair_a_left = [cmplx(0.017_rk, -0.011_rk, kind=rk), &
                   cmplx(0.31_rk, 0.07_rk, kind=rk), &
                   cmplx(-0.08_rk, 0.19_rk, kind=rk), &
                   cmplx(0.12_rk, -0.16_rk, kind=rk)]
    pair_a_right = [cmplx(-0.009_rk, 0.015_rk, kind=rk), &
                    cmplx(0.27_rk, 0.11_rk, kind=rk), &
                    cmplx(-0.03_rk, 0.21_rk, kind=rk), &
                    cmplx(0.16_rk, -0.10_rk, kind=rk)]
    pair_b_left = [cmplx(-0.013_rk, 0.005_rk, kind=rk), &
                   cmplx(0.28_rk, -0.09_rk, kind=rk), &
                   cmplx(-0.06_rk, -0.17_rk, kind=rk), &
                   cmplx(0.10_rk, 0.13_rk, kind=rk)]
    pair_b_right = [cmplx(0.006_rk, -0.012_rk, kind=rk), &
                    cmplx(0.25_rk, -0.13_rk, kind=rk), &
                    cmplx(-0.01_rk, -0.20_rk, kind=rk), &
                    cmplx(0.14_rk, 0.08_rk, kind=rk)]
    exchange_left = [cmplx(0.012_rk, -0.003_rk, kind=rk), &
                     cmplx(-0.007_rk, 0.004_rk, kind=rk), &
                     cmplx(0.005_rk, 0.002_rk, kind=rk)]
    exchange_right = [cmplx(0.009_rk, -0.001_rk, kind=rk), &
                      cmplx(-0.011_rk, 0.006_rk, kind=rk), &
                      cmplx(0.003_rk, -0.004_rk, kind=rk)]
    expected = [ &
      cmplx(5.85570667597123e-3_rk, -9.72258425874292e-3_rk, kind=rk), &
      cmplx(-9.73495235104594e-2_rk, 9.65147393449474e-2_rk, kind=rk), &
      cmplx(-1.66791927689114e-2_rk, 1.47821668870606e-2_rk, kind=rk), &
      cmplx(2.08188996921606e-2_rk, -9.57975429655730e-3_rk, kind=rk)]

    call advance_legacy_riccati_segment( &
      coherence_in, coherence_out, pair_a_left, pair_a_right, &
      pair_b_left, pair_b_right, &
      cmplx(-0.035_rk, 0.82_rk, kind=rk), &
      cmplx(-0.019_rk, 0.91_rk, kind=rk), &
      exchange_left, exchange_right, 0.37_rk, 7)

    call require(maxval(abs(coherence_out - expected)) < 5.0e-15_rk, &
                 "modern Riccati segment differs from frozen new_src output")
  end subroutine test_frozen_legacy_checkpoint


  subroutine test_zero_distance()
    complex(rk) :: coherence_in(0:3), coherence_out(0:3)
    complex(rk) :: pair_a(0:3), pair_b(0:3), exchange(3)

    coherence_in = [cmplx(0.2_rk, -0.1_rk, kind=rk), &
                    cmplx(-0.3_rk, 0.4_rk, kind=rk), &
                    cmplx(0.5_rk, 0.1_rk, kind=rk), &
                    cmplx(-0.2_rk, -0.6_rk, kind=rk)]
    pair_a = cmplx(0.1_rk, 0.2_rk, kind=rk)
    pair_b = cmplx(-0.2_rk, 0.3_rk, kind=rk)
    exchange = cmplx(0.05_rk, -0.02_rk, kind=rk)

    call advance_legacy_riccati_segment( &
      coherence_in, coherence_out, pair_a, pair_a, pair_b, pair_b, &
      cmplx(0.0_rk, 1.0_rk, kind=rk), &
      cmplx(0.0_rk, 1.0_rk, kind=rk), exchange, exchange, 0.0_rk, 3)
    call require(maxval(abs(coherence_out - coherence_in)) < tiny(1.0_rk), &
                 "zero-length Riccati segment changed its input")
  end subroutine test_zero_distance


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_riccati_segment
