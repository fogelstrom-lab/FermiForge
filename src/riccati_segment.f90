module riccati_segment
  use he3_kinds, only : rk
  implicit none
  private

  public :: advance_legacy_riccati_segment
  public :: legacy_bulk_coherence

contains

  pure subroutine legacy_bulk_coherence(pair_a, pair_b, effective_energy, &
                                         exchange_field, reverse_branch, &
                                         coherence)
    complex(rk), intent(in) :: pair_a(0:3), pair_b(0:3)
    complex(rk), intent(in) :: effective_energy, exchange_field(3)
    logical, intent(in) :: reverse_branch
    complex(rk), intent(out) :: coherence(0:3)

    complex(rk), parameter :: imaginary_unit = &
      cmplx(0.0_rk, 1.0_rk, kind=rk)
    complex(rk), parameter :: complex_one = cmplx(1.0_rk, 0.0_rk, kind=rk)
    complex(rk) :: denominator, inverse_denominator, lambda

    lambda = pair_a(1) * pair_b(1) + pair_a(2) * pair_b(2) + &
             pair_a(3) * pair_b(3)
    lambda = imaginary_unit * sqrt(lambda - effective_energy * effective_energy)
    inverse_denominator = complex_one / (effective_energy + lambda)
    coherence(1:3) = -pair_a(1:3) * inverse_denominator
    if (reverse_branch) then
      denominator = effective_energy + coherence(1) * pair_b(1) + &
                    coherence(2) * pair_b(2) + coherence(3) * pair_b(3)
    else
      denominator = effective_energy - coherence(1) * pair_b(1) - &
                    coherence(2) * pair_b(2) - coherence(3) * pair_b(3)
    end if
    coherence(0) = (coherence(1) * exchange_field(1) + &
                    coherence(2) * exchange_field(2) + &
                    coherence(3) * exchange_field(3)) / denominator
  end subroutine legacy_bulk_coherence

  pure subroutine advance_legacy_riccati_segment( &
      coherence_in, coherence_out, pair_a_left, pair_a_right, &
      pair_b_left, pair_b_right, energy_left, energy_right, &
      exchange_left, exchange_right, distance, substep_count)
    complex(rk), intent(in) :: coherence_in(0:3)
    complex(rk), intent(out) :: coherence_out(0:3)
    complex(rk), intent(in) :: pair_a_left(0:3), pair_a_right(0:3)
    complex(rk), intent(in) :: pair_b_left(0:3), pair_b_right(0:3)
    complex(rk), intent(in) :: energy_left, energy_right
    complex(rk), intent(in) :: exchange_left(3), exchange_right(3)
    real(rk), intent(in) :: distance
    integer, intent(in) :: substep_count

    complex(rk), parameter :: imaginary_unit = &
      cmplx(0.0_rk, 1.0_rk, kind=rk)
    complex(rk), parameter :: complex_two = cmplx(2.0_rk, 0.0_rk, kind=rk)
    complex(rk), parameter :: complex_half = cmplx(0.5_rk, 0.0_rk, kind=rk)
    complex(rk), parameter :: complex_one_sixth = &
      cmplx(1.0_rk / 6.0_rk, 0.0_rk, kind=rk)
    complex(rk) :: energy_1, energy_2, energy_3, energy_increment
    complex(rk) :: exchange_1(3), exchange_2(3), exchange_3(3)
    complex(rk) :: exchange_increment(3)
    complex(rk) :: pair_a_1(0:3), pair_a_2(0:3), pair_a_3(0:3)
    complex(rk) :: pair_a_increment(0:3)
    complex(rk) :: pair_b_1(0:3), pair_b_2(0:3), pair_b_3(0:3)
    complex(rk) :: pair_b_increment(0:3)
    complex(rk) :: stage_base(0:3), stage_1(0:3), stage_2(0:3)
    complex(rk) :: stage_3(0:3), stage_4(0:3)
    complex(rk) :: slope_1(0:3), slope_2(0:3)
    complex(rk) :: slope_3(0:3), slope_4(0:3)
    complex(rk) :: scalar_product, vector_product, sum_1, sum_2
    complex(rk) :: exchange_dot, exchange_x, exchange_y, exchange_z
    complex(rk) :: complex_half_fraction, complex_step, complex_twice_step
    real(rk) :: half_fraction, fraction, step, twice_step
    integer :: substep

    if (substep_count < 1) &
      error stop "Riccati segment needs at least one RK4 substep"

    ! This is a direct explicit-kind transcription of new_src/riccati.f90:
    ! ricc. The seemingly unusual half-step endpoint increments are retained
    ! exactly so that the modern kernel remains a numerical oracle match.
    half_fraction = 1.0_rk / real(2 * substep_count, rk)
    fraction = half_fraction + half_fraction
    step = distance * fraction
    twice_step = step + step
    complex_half_fraction = cmplx(half_fraction, 0.0_rk, kind=rk)
    complex_step = cmplx(step, 0.0_rk, kind=rk)
    complex_twice_step = cmplx(twice_step, 0.0_rk, kind=rk)

    energy_1 = imaginary_unit * energy_left * complex_twice_step
    energy_2 = imaginary_unit * energy_right * complex_twice_step
    energy_3 = energy_1
    energy_increment = (energy_2 - energy_1) * complex_half_fraction

    exchange_1 = imaginary_unit * exchange_left * complex_twice_step
    exchange_2 = imaginary_unit * exchange_right * complex_twice_step
    exchange_3 = exchange_1
    exchange_increment = (exchange_2 - exchange_1) * complex_half_fraction

    pair_a_1 = imaginary_unit * pair_a_left * complex_step
    pair_a_2 = imaginary_unit * pair_a_right * complex_step
    pair_a_3 = pair_a_1
    pair_a_increment = (pair_a_2 - pair_a_1) * complex_half_fraction

    pair_b_1 = imaginary_unit * pair_b_left * complex_step
    pair_b_2 = imaginary_unit * pair_b_right * complex_step
    pair_b_3 = pair_b_1
    pair_b_increment = (pair_b_2 - pair_b_1) * complex_half_fraction

    stage_3 = coherence_in
    do substep = 1, substep_count
      energy_1 = energy_3
      energy_2 = energy_1 + energy_increment
      energy_3 = energy_2 + energy_increment

      exchange_1 = exchange_3
      exchange_2 = exchange_1 + exchange_increment
      exchange_3 = exchange_2 + exchange_increment

      pair_a_1 = pair_a_3
      pair_a_2 = pair_a_1 + pair_a_increment
      pair_a_3 = pair_a_2 + pair_a_increment

      pair_b_1 = pair_b_3
      pair_b_2 = pair_b_1 + pair_b_increment
      pair_b_3 = pair_b_2 + pair_b_increment

      stage_base = stage_3
      stage_1 = stage_base

      scalar_product = pair_b_1(0) * stage_1(0)
      vector_product = pair_b_1(1) * stage_1(1) + &
                       pair_b_1(2) * stage_1(2) + &
                       pair_b_1(3) * stage_1(3)
      sum_1 = energy_1 + complex_two * vector_product
      sum_2 = sum_1 + complex_two * scalar_product
      scalar_product = stage_1(0) * stage_1(0)
      vector_product = stage_1(1) * stage_1(1) + &
                       stage_1(2) * stage_1(2) + &
                       stage_1(3) * stage_1(3)
      exchange_dot = exchange_1(1) * stage_1(1) + &
                     exchange_1(2) * stage_1(2) + &
                     exchange_1(3) * stage_1(3)
      exchange_x = exchange_1(1) * stage_1(0)
      exchange_y = exchange_1(2) * stage_1(0)
      exchange_z = exchange_1(3) * stage_1(0)
      slope_1(0) = (scalar_product + vector_product) * pair_b_1(0) + &
                   stage_1(0) * sum_1 + pair_a_1(0) - exchange_dot
      slope_1(1) = (scalar_product - vector_product) * pair_b_1(1) + &
                   stage_1(1) * sum_2 + pair_a_1(1) - exchange_x
      slope_1(2) = (scalar_product - vector_product) * pair_b_1(2) + &
                   stage_1(2) * sum_2 + pair_a_1(2) - exchange_y
      slope_1(3) = (scalar_product - vector_product) * pair_b_1(3) + &
                   stage_1(3) * sum_2 + pair_a_1(3) - exchange_z

      stage_2 = stage_base + complex_half * slope_1
      scalar_product = pair_b_2(0) * stage_2(0)
      vector_product = pair_b_2(1) * stage_2(1) + &
                       pair_b_2(2) * stage_2(2) + &
                       pair_b_2(3) * stage_2(3)
      sum_1 = energy_2 + complex_two * vector_product
      sum_2 = sum_1 + complex_two * scalar_product
      scalar_product = stage_2(0) * stage_2(0)
      vector_product = stage_2(1) * stage_2(1) + &
                       stage_2(2) * stage_2(2) + &
                       stage_2(3) * stage_2(3)
      exchange_dot = exchange_2(1) * stage_2(1) + &
                     exchange_2(2) * stage_2(2) + &
                     exchange_2(3) * stage_2(3)
      exchange_x = exchange_2(1) * stage_2(0)
      exchange_y = exchange_2(2) * stage_2(0)
      exchange_z = exchange_2(3) * stage_2(0)
      slope_2(0) = (scalar_product + vector_product) * pair_b_2(0) + &
                   stage_2(0) * sum_1 + pair_a_2(0) - exchange_dot
      slope_2(1) = (scalar_product - vector_product) * pair_b_2(1) + &
                   stage_2(1) * sum_2 + pair_a_2(1) - exchange_x
      slope_2(2) = (scalar_product - vector_product) * pair_b_2(2) + &
                   stage_2(2) * sum_2 + pair_a_2(2) - exchange_y
      slope_2(3) = (scalar_product - vector_product) * pair_b_2(3) + &
                   stage_2(3) * sum_2 + pair_a_2(3) - exchange_z

      stage_3 = stage_base + complex_half * slope_2
      scalar_product = pair_b_2(0) * stage_3(0)
      vector_product = pair_b_2(1) * stage_3(1) + &
                       pair_b_2(2) * stage_3(2) + &
                       pair_b_2(3) * stage_3(3)
      sum_1 = energy_2 + complex_two * vector_product
      sum_2 = sum_1 + complex_two * scalar_product
      scalar_product = stage_3(0) * stage_3(0)
      vector_product = stage_3(1) * stage_3(1) + &
                       stage_3(2) * stage_3(2) + &
                       stage_3(3) * stage_3(3)
      exchange_dot = exchange_2(1) * stage_3(1) + &
                     exchange_2(2) * stage_3(2) + &
                     exchange_2(3) * stage_3(3)
      exchange_x = exchange_2(1) * stage_3(0)
      exchange_y = exchange_2(2) * stage_3(0)
      exchange_z = exchange_2(3) * stage_3(0)
      slope_3(0) = (scalar_product + vector_product) * pair_b_2(0) + &
                   stage_3(0) * sum_1 + pair_a_2(0) - exchange_dot
      slope_3(1) = (scalar_product - vector_product) * pair_b_2(1) + &
                   stage_3(1) * sum_2 + pair_a_2(1) - exchange_x
      slope_3(2) = (scalar_product - vector_product) * pair_b_2(2) + &
                   stage_3(2) * sum_2 + pair_a_2(2) - exchange_y
      slope_3(3) = (scalar_product - vector_product) * pair_b_2(3) + &
                   stage_3(3) * sum_2 + pair_a_2(3) - exchange_z

      stage_4 = stage_base + slope_3
      scalar_product = pair_b_3(0) * stage_4(0)
      vector_product = pair_b_3(1) * stage_4(1) + &
                       pair_b_3(2) * stage_4(2) + &
                       pair_b_3(3) * stage_4(3)
      sum_1 = energy_3 + complex_two * vector_product
      sum_2 = sum_1 + complex_two * scalar_product
      scalar_product = stage_4(0) * stage_4(0)
      vector_product = stage_4(1) * stage_4(1) + &
                       stage_4(2) * stage_4(2) + &
                       stage_4(3) * stage_4(3)
      exchange_dot = exchange_3(1) * stage_4(1) + &
                     exchange_3(2) * stage_4(2) + &
                     exchange_3(3) * stage_4(3)
      exchange_x = exchange_3(1) * stage_4(0)
      exchange_y = exchange_3(2) * stage_4(0)
      exchange_z = exchange_3(3) * stage_4(0)
      slope_4(0) = (scalar_product + vector_product) * pair_b_3(0) + &
                   stage_4(0) * sum_1 + pair_a_3(0) - exchange_dot
      slope_4(1) = (scalar_product - vector_product) * pair_b_3(1) + &
                   stage_4(1) * sum_2 + pair_a_3(1) - exchange_x
      slope_4(2) = (scalar_product - vector_product) * pair_b_3(2) + &
                   stage_4(2) * sum_2 + pair_a_3(2) - exchange_y
      slope_4(3) = (scalar_product - vector_product) * pair_b_3(3) + &
                   stage_4(3) * sum_2 + pair_a_3(3) - exchange_z

      stage_3 = stage_base + &
        (slope_1 + slope_2 + slope_2 + slope_3 + slope_3 + slope_4) * &
        complex_one_sixth
    end do

    coherence_out = stage_3
  end subroutine advance_legacy_riccati_segment

end module riccati_segment
