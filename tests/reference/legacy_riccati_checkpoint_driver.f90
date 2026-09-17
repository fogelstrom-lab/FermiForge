program legacy_riccati_checkpoint_driver
  use riccati, only : ricc
  implicit none

  complex :: coherence_in(0:3), coherence_out(0:3)
  complex :: pair_a_left(0:3), pair_a_right(0:3)
  complex :: pair_b_left(0:3), pair_b_right(0:3)
  complex :: exchange_left(3), exchange_right(3)
  complex :: energy_left, energy_right
  real :: distance
  integer :: component

  coherence_in = [cmplx(0.013, -0.021), cmplx(-0.14, 0.032), &
                  cmplx(0.071, 0.055), cmplx(-0.023, -0.091)]
  pair_a_left = [cmplx(0.017, -0.011), cmplx(0.31, 0.07), &
                 cmplx(-0.08, 0.19), cmplx(0.12, -0.16)]
  pair_a_right = [cmplx(-0.009, 0.015), cmplx(0.27, 0.11), &
                  cmplx(-0.03, 0.21), cmplx(0.16, -0.10)]
  pair_b_left = [cmplx(-0.013, 0.005), cmplx(0.28, -0.09), &
                 cmplx(-0.06, -0.17), cmplx(0.10, 0.13)]
  pair_b_right = [cmplx(0.006, -0.012), cmplx(0.25, -0.13), &
                  cmplx(-0.01, -0.20), cmplx(0.14, 0.08)]
  exchange_left = [cmplx(0.012, -0.003), cmplx(-0.007, 0.004), &
                   cmplx(0.005, 0.002)]
  exchange_right = [cmplx(0.009, -0.001), cmplx(-0.011, 0.006), &
                    cmplx(0.003, -0.004)]
  energy_left = cmplx(-0.035, 0.82)
  energy_right = cmplx(-0.019, 0.91)
  distance = 0.37

  call ricc(coherence_in, coherence_out, pair_a_left, pair_a_right, &
            pair_b_left, pair_b_right, energy_left, energy_right, &
            exchange_left, exchange_right, distance, 7)

  do component = 0, 3
    write (*, '(i1,2(1x,es25.17))') component, &
      real(coherence_out(component)), aimag(coherence_out(component))
  end do
end program legacy_riccati_checkpoint_driver
