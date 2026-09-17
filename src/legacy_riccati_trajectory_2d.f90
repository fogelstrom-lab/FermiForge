module legacy_riccati_trajectory_2d
  use he3_kinds, only : rk
  use straight_trajectory_2d, only : straight_trajectory_2d_t
  use trajectory_self_energy_2d, only : trajectory_self_energy_2d_t, &
                                        legacy_forward_branch, &
                                        legacy_reverse_branch, &
                                        make_legacy_riccati_coefficients
  use riccati_segment, only : advance_legacy_riccati_segment, &
                              legacy_bulk_coherence
  implicit none
  private

  public :: propagate_legacy_riccati_to_target
  public :: propagate_legacy_riccati_path_to_target

contains

  subroutine propagate_legacy_riccati_to_target( &
      trajectory, self_energy, spectral_energy, substep_count, &
      boundary_relaxation_distance, coherence_from_entry, coherence_from_exit)
    type(straight_trajectory_2d_t), intent(in) :: trajectory
    type(trajectory_self_energy_2d_t), intent(in) :: self_energy
    complex(rk), intent(in) :: spectral_energy
    integer, intent(in) :: substep_count
    real(rk), intent(in) :: boundary_relaxation_distance
    complex(rk), allocatable, intent(out) :: coherence_from_entry(:, :)
    complex(rk), allocatable, intent(out) :: coherence_from_exit(:, :)

    if (.not. trajectory%is_valid()) &
      error stop "cannot propagate on an invalid 2D trajectory"
    call propagate_legacy_riccati_path_to_target( &
      trajectory%path_coordinate, trajectory%target_sample, self_energy, &
      spectral_energy, substep_count, boundary_relaxation_distance, &
      coherence_from_entry, coherence_from_exit)
  end subroutine propagate_legacy_riccati_to_target


  subroutine propagate_legacy_riccati_path_to_target( &
      path_coordinate, target_sample, self_energy, spectral_energy, &
      substep_count, boundary_relaxation_distance, coherence_from_entry, &
      coherence_from_exit)
    real(rk), intent(in) :: path_coordinate(:)
    integer, intent(in) :: target_sample
    type(trajectory_self_energy_2d_t), intent(in) :: self_energy
    complex(rk), intent(in) :: spectral_energy
    integer, intent(in) :: substep_count
    real(rk), intent(in) :: boundary_relaxation_distance
    complex(rk), allocatable, intent(out) :: coherence_from_entry(:, :)
    complex(rk), allocatable, intent(out) :: coherence_from_exit(:, :)

    complex(rk) :: coherence(0:3), updated(0:3)
    complex(rk) :: pair_a_left(0:3), pair_a_right(0:3)
    complex(rk) :: pair_b_left(0:3), pair_b_right(0:3)
    complex(rk) :: energy_left, energy_right
    complex(rk) :: exchange_left(3), exchange_right(3)
    real(rk) :: distance
    integer :: number_of_samples, sample

    if (.not. self_energy%is_valid()) &
      error stop "cannot propagate invalid trajectory self energies"
    if (size(path_coordinate) < 1 .or. &
        self_energy%sample_count() /= size(path_coordinate)) &
      error stop "trajectory and self-energy sample counts differ"
    if (target_sample < 1 .or. target_sample > size(path_coordinate)) &
      error stop "Riccati target sample lies outside the path"
    if (substep_count < 1) &
      error stop "Riccati trajectory needs at least one substep per interval"
    if (boundary_relaxation_distance < 0.0_rk) &
      error stop "Riccati boundary relaxation distance cannot be negative"
    if (size(path_coordinate) > 1) then
      if (any(path_coordinate(2:) <= &
              path_coordinate(:size(path_coordinate) - 1))) &
        error stop "Riccati trajectory path coordinates are not increasing"
    end if

    number_of_samples = size(path_coordinate)
    allocate(coherence_from_entry(0:3, number_of_samples), &
             source=cmplx(0.0_rk, 0.0_rk, kind=rk))
    allocate(coherence_from_exit(0:3, number_of_samples), &
             source=cmplx(0.0_rk, 0.0_rk, kind=rk))

    call make_legacy_riccati_coefficients( &
      self_energy, 1, legacy_forward_branch, spectral_energy, &
      pair_a_left, pair_b_left, energy_left, exchange_left)
    call legacy_bulk_coherence(pair_a_left, pair_b_left, energy_left, &
                               exchange_left, .false., coherence)
    if (boundary_relaxation_distance > 0.0_rk) then
      call advance_legacy_riccati_segment( &
        coherence, updated, pair_a_left, pair_a_left, pair_b_left, pair_b_left, &
        energy_left, energy_left, exchange_left, exchange_left, &
        boundary_relaxation_distance, substep_count)
      coherence = updated
    end if
    coherence_from_entry(:, 1) = coherence

    do sample = 2, target_sample
      call make_legacy_riccati_coefficients( &
        self_energy, sample, legacy_forward_branch, spectral_energy, &
        pair_a_right, pair_b_right, energy_right, exchange_right)
      distance = path_coordinate(sample) - path_coordinate(sample - 1)
      call advance_legacy_riccati_segment( &
        coherence, updated, pair_a_left, pair_a_right, pair_b_left, pair_b_right, &
        energy_left, energy_right, exchange_left, exchange_right, distance, &
        substep_count)
      coherence = updated
      coherence_from_entry(:, sample) = coherence
      pair_a_left = pair_a_right
      pair_b_left = pair_b_right
      energy_left = energy_right
      exchange_left = exchange_right
    end do

    call make_legacy_riccati_coefficients( &
      self_energy, number_of_samples, legacy_reverse_branch, spectral_energy, &
      pair_a_left, pair_b_left, energy_left, exchange_left)
    call legacy_bulk_coherence(pair_a_left, pair_b_left, energy_left, &
                               exchange_left, .true., coherence)
    if (boundary_relaxation_distance > 0.0_rk) then
      call advance_legacy_riccati_segment( &
        coherence, updated, pair_a_left, pair_a_left, pair_b_left, pair_b_left, &
        energy_left, energy_left, exchange_left, exchange_left, &
        boundary_relaxation_distance, substep_count)
      coherence = updated
    end if
    coherence_from_exit(:, number_of_samples) = coherence

    do sample = number_of_samples - 1, target_sample, -1
      call make_legacy_riccati_coefficients( &
        self_energy, sample, legacy_reverse_branch, spectral_energy, &
        pair_a_right, pair_b_right, energy_right, exchange_right)
      distance = path_coordinate(sample + 1) - path_coordinate(sample)
      call advance_legacy_riccati_segment( &
        coherence, updated, pair_a_left, pair_a_right, pair_b_left, pair_b_right, &
        energy_left, energy_right, exchange_left, exchange_right, distance, &
        substep_count)
      coherence = updated
      coherence_from_exit(:, sample) = coherence
      pair_a_left = pair_a_right
      pair_b_left = pair_b_right
      energy_left = energy_right
      exchange_left = exchange_right
    end do
  end subroutine propagate_legacy_riccati_path_to_target

end module legacy_riccati_trajectory_2d
