module he3_self_consistency_integrand
  use he3_kinds, only : rk
  use quasiclassical_propagator, only : quasiclassical_propagator_t
  implicit none
  private

  complex(rk), parameter :: complex_quarter = &
    cmplx(0.25_rk, 0.0_rk, kind=rk)
  complex(rk), parameter :: negative_imaginary_quarter = &
    cmplx(0.0_rk, -0.25_rk, kind=rk)

  type, public :: he3_point_map_accumulator_t
    complex(rk) :: gap(3, 3) = cmplx(0.0_rk, 0.0_rk, kind=rk)
    real(rk) :: current_mean_field(3) = 0.0_rk
    integer :: contribution_count = 0
  contains
    procedure :: reset => reset_he3_point_map_accumulator
  end type he3_point_map_accumulator_t

  public :: accumulate_he3_point_map_contribution
  public :: finalize_he3_point_map

contains

  pure subroutine reset_he3_point_map_accumulator(self)
    class(he3_point_map_accumulator_t), intent(inout) :: self

    self%gap = cmplx(0.0_rk, 0.0_rk, kind=rk)
    self%current_mean_field = 0.0_rk
    self%contribution_count = 0
  end subroutine reset_he3_point_map_accumulator


  pure subroutine accumulate_he3_point_map_contribution( &
      accumulator, propagator, momentum, direction_weight, pole_weight, &
      temperature)
    type(he3_point_map_accumulator_t), intent(inout) :: accumulator
    type(quasiclassical_propagator_t), intent(in) :: propagator
    real(rk), intent(in) :: momentum(3), direction_weight
    real(rk), intent(in) :: pole_weight, temperature

    complex(rk) :: spin_integrand(3)
    complex(rk) :: weighted_spin_integrand(3)
    real(rk) :: current_weight, gap_weight, particle_trace
    integer :: orbital, spin

    ! Exact component extraction from new_src/riccati.f90:green.  The
    ! anomalous blocks use the source's dimensionless g_hat**2=-1
    ! normalization.  direction_weight is awei*pwei and pole_weight is Rp.
    spin_integrand(1) = ( &
      propagator%anomalous(2, 2) - propagator%anomalous(1, 1) - &
      conjg(propagator%tilde_anomalous(2, 2) - &
            propagator%tilde_anomalous(1, 1))) * complex_quarter
    spin_integrand(2) = ( &
      propagator%anomalous(2, 2) + propagator%anomalous(1, 1) - &
      conjg(propagator%tilde_anomalous(2, 2) + &
            propagator%tilde_anomalous(1, 1))) * &
      negative_imaginary_quarter
    spin_integrand(3) = ( &
      propagator%anomalous(1, 2) + propagator%anomalous(2, 1) - &
      conjg(propagator%tilde_anomalous(2, 1) + &
            propagator%tilde_anomalous(1, 2))) * complex_quarter

    gap_weight = 3.0_rk * direction_weight * pole_weight
    weighted_spin_integrand = &
      cmplx(gap_weight, 0.0_rk, kind=rk) * spin_integrand
    do spin = 1, 3
      do orbital = 1, 3
        accumulator%gap(spin, orbital) = &
          accumulator%gap(spin, orbital) + &
          weighted_spin_integrand(spin) * &
          cmplx(momentum(orbital), 0.0_rk, kind=rk)
      end do
    end do

    particle_trace = real( &
      propagator%particle(1, 1) + propagator%particle(2, 2), kind=rk)
    current_weight = 0.5_rk * temperature * direction_weight * &
                     pole_weight * particle_trace
    accumulator%current_mean_field = accumulator%current_mean_field + &
      current_weight * momentum
    accumulator%contribution_count = accumulator%contribution_count + 1
  end subroutine accumulate_he3_point_map_contribution


  pure subroutine finalize_he3_point_map( &
      accumulator, gap_prefactor, mapped_gap, mapped_current_mean_field)
    type(he3_point_map_accumulator_t), intent(in) :: accumulator
    real(rk), intent(in) :: gap_prefactor
    complex(rk), intent(out) :: mapped_gap(3, 3)
    real(rk), intent(out) :: mapped_current_mean_field(3)

    mapped_gap = cmplx(gap_prefactor, 0.0_rk, kind=rk) * accumulator%gap
    mapped_current_mean_field = accumulator%current_mean_field
  end subroutine finalize_he3_point_map

end module he3_self_consistency_integrand
