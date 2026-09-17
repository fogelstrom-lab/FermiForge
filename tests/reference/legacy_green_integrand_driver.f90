program legacy_green_integrand_driver
  use global_variables
  use riccati, only : green
  implicit none

  complex :: accumulator(12)
  complex :: gamma_forward_0, gamma_forward_1
  complex :: gamma_forward_2, gamma_forward_3
  complex :: gamma_reverse_0, gamma_reverse_1
  complex :: gamma_reverse_2, gamma_reverse_3
  integer :: component

  ! A deterministic one-quadrature-point fixture for the self-consistency
  ! contribution in new_src/riccati.f90:green. The coherence values are the
  ! target values frozen by legacy_riccati_trajectory_driver.f90.
  t = 0.37
  kx(1) = cos(0.41)
  ky(1) = sin(0.41)
  awei(1) = 0.11
  kz(1) = 0.35
  pwei(1) = 0.17
  Rp(1) = 0.42

  gamma_forward_0 = cmplx(0.0, 0.0)
  gamma_forward_1 = cmplx(-3.45416635611957e-2, 1.17489705712979e-1)
  gamma_forward_2 = cmplx(-9.01737612952456e-2, -5.58797007439699e-2)
  gamma_forward_3 = cmplx(8.02906423983453e-2, 7.10955985562741e-2)

  gamma_reverse_0 = cmplx(0.0, 0.0)
  gamma_reverse_1 = cmplx(1.62659859791399e-2, 1.42678313138293e-1)
  gamma_reverse_2 = cmplx(9.95709656830779e-2, -4.16922195641589e-2)
  gamma_reverse_3 = cmplx(-6.37043298757304e-2, 5.90495834745839e-2)

  accumulator = cmplx(0.0, 0.0)
  call green(1, 1, 1, accumulator, &
             gamma_forward_0, gamma_forward_1, gamma_forward_2, &
             gamma_forward_3, gamma_reverse_0, gamma_reverse_1, &
             gamma_reverse_2, gamma_reverse_3)

  do component = 1, 12
    write(*, '(i3,2(1x,es25.17))') component, accumulator(component)
  end do
end program legacy_green_integrand_driver
