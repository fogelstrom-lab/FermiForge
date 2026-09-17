program legacy_riccati_trajectory_driver
  use Global_variables, only : mx
  use riccati, only : bulk_gammas, ricc, setenergy1, setenergy2
  implicit none

  complex :: sem(4, -mx:mx), coherence(0:3), updated(0:3)
  complex :: pair_a_left(0:3), pair_a_right(0:3)
  complex :: pair_b_left(0:3), pair_b_right(0:3)
  complex :: exchange_left(3), exchange_right(3)
  complex :: energy_left, energy_right, spectral_energy
  integer :: component, point

  sem = cmplx(0.0, 0.0)
  do point = -2, 2
    sem(1, point) = cmplx(0.23 + 0.015 * point, 0.04 - 0.012 * point)
    sem(2, point) = cmplx(-0.08 + 0.009 * point, 0.17 + 0.006 * point)
    sem(3, point) = cmplx(0.11 - 0.007 * point, -0.13 + 0.01 * point)
    sem(4, point) = cmplx(0.018 * point, 0.0)
  end do
  spectral_energy = cmplx(0.0, 0.83)

  call bulk_gammas(coherence, sem, spectral_energy, -2, 1)
  call write_checkpoint('F', -2, coherence)
  call setenergy1(-2, spectral_energy, pair_a_left, pair_b_left, &
                  energy_left, exchange_left, sem)
  do point = -1, 0
    call setenergy1(point, spectral_energy, pair_a_right, pair_b_right, &
                    energy_right, exchange_right, sem)
    call ricc(coherence, updated, pair_a_left, pair_a_right, &
              pair_b_left, pair_b_right, energy_left, energy_right, &
              exchange_left, exchange_right, 0.2, 3)
    coherence = updated
    pair_a_left = pair_a_right
    pair_b_left = pair_b_right
    energy_left = energy_right
    exchange_left = exchange_right
    call write_checkpoint('F', point, coherence)
  end do

  call bulk_gammas(coherence, sem, spectral_energy, 2, 2)
  call write_checkpoint('R', 2, coherence)
  call setenergy2(2, spectral_energy, pair_a_left, pair_b_left, &
                  energy_left, exchange_left, sem)
  do point = 1, 0, -1
    call setenergy2(point, spectral_energy, pair_a_right, pair_b_right, &
                    energy_right, exchange_right, sem)
    call ricc(coherence, updated, pair_a_left, pair_a_right, &
              pair_b_left, pair_b_right, energy_left, energy_right, &
              exchange_left, exchange_right, 0.2, 3)
    coherence = updated
    pair_a_left = pair_a_right
    pair_b_left = pair_b_right
    energy_left = energy_right
    exchange_left = exchange_right
    call write_checkpoint('R', point, coherence)
  end do

contains

  subroutine write_checkpoint(label, point_index, values)
    character(len=1), intent(in) :: label
    integer, intent(in) :: point_index
    complex, intent(in) :: values(0:3)

    do component = 0, 3
      write (*, '(a1,1x,i3,1x,i1,2(1x,es25.17))') label, point_index, &
        component, real(values(component)), aimag(values(component))
    end do
  end subroutine write_checkpoint

end program legacy_riccati_trajectory_driver
