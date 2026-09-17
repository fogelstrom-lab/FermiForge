module trajectory_self_energy_2d
  use he3_kinds, only : rk
  use spinful_state_2d, only : spinful_state_2d_t
  use straight_trajectory_2d, only : straight_trajectory_2d_t, &
                                     sample_pair_potential_along_trajectory_2d
  implicit none
  private

  integer, parameter, public :: legacy_forward_branch = 1
  integer, parameter, public :: legacy_reverse_branch = 2

  type, public :: trajectory_self_energy_2d_t
    complex(rk), allocatable :: triplet_pair_potential(:, :)
    real(rk), allocatable :: diagonal_shift(:)
  contains
    procedure :: sample_count => self_energy_sample_count
    procedure :: is_valid => self_energy_is_valid
  end type trajectory_self_energy_2d_t

  public :: sample_trajectory_self_energy_2d
  public :: make_legacy_riccati_coefficients

contains

  subroutine sample_trajectory_self_energy_2d(state, trajectory, &
                                              feedback_scale, self_energy)
    type(spinful_state_2d_t), intent(in) :: state
    type(straight_trajectory_2d_t), intent(in) :: trajectory
    real(rk), intent(in) :: feedback_scale
    type(trajectory_self_energy_2d_t), intent(out) :: self_energy

    real(rk), allocatable :: sampled_mean_field(:, :)
    integer :: sample

    call sample_pair_potential_along_trajectory_2d( &
      state, trajectory, self_energy%triplet_pair_potential, sampled_mean_field)
    allocate(self_energy%diagonal_shift(trajectory%sample_count()), &
             source=0.0_rk)
    do sample = 1, trajectory%sample_count()
      self_energy%diagonal_shift(sample) = feedback_scale * &
        dot_product(trajectory%momentum, sampled_mean_field(:, sample))
    end do
  end subroutine sample_trajectory_self_energy_2d


  pure subroutine make_legacy_riccati_coefficients( &
      self_energy, sample, branch, spectral_energy, pair_a, pair_b, &
      effective_energy, exchange_field)
    type(trajectory_self_energy_2d_t), intent(in) :: self_energy
    integer, intent(in) :: sample, branch
    complex(rk), intent(in) :: spectral_energy
    complex(rk), intent(out) :: pair_a(0:3), pair_b(0:3)
    complex(rk), intent(out) :: effective_energy
    complex(rk), intent(out) :: exchange_field(3)

    if (.not. self_energy%is_valid()) &
      error stop "cannot form Riccati coefficients from invalid self-energy data"
    if (sample < 1 .or. sample > self_energy%sample_count()) &
      error stop "Riccati self-energy sample index is outside the trajectory"

    pair_a = cmplx(0.0_rk, 0.0_rk, kind=rk)
    pair_b = cmplx(0.0_rk, 0.0_rk, kind=rk)
    select case (branch)
    case (legacy_forward_branch)
      pair_a(1:3) = self_energy%triplet_pair_potential(:, sample)
      pair_b = conjg(pair_a)
    case (legacy_reverse_branch)
      pair_b(1:3) = self_energy%triplet_pair_potential(:, sample)
      pair_a = conjg(pair_b)
    case default
      error stop "unknown legacy Riccati propagation branch"
    end select

    effective_energy = spectral_energy - &
      cmplx(self_energy%diagonal_shift(sample), 0.0_rk, kind=rk)
    ! The active new_src calculation zeros the spin-vector diagonal
    ! self-energy in setenergy1/2. Keep that source behavior explicit here.
    exchange_field = cmplx(0.0_rk, 0.0_rk, kind=rk)
  end subroutine make_legacy_riccati_coefficients


  pure integer function self_energy_sample_count(self) result(number_of_samples)
    class(trajectory_self_energy_2d_t), intent(in) :: self

    if (allocated(self%diagonal_shift)) then
      number_of_samples = size(self%diagonal_shift)
    else
      number_of_samples = 0
    end if
  end function self_energy_sample_count


  pure logical function self_energy_is_valid(self) result(valid)
    class(trajectory_self_energy_2d_t), intent(in) :: self

    valid = allocated(self%triplet_pair_potential) .and. &
            allocated(self%diagonal_shift)
    if (.not. valid) return
    valid = size(self%triplet_pair_potential, 1) == 3 .and. &
            size(self%triplet_pair_potential, 2) == &
              size(self%diagonal_shift) .and. &
            size(self%diagonal_shift) >= 1
  end function self_energy_is_valid

end module trajectory_self_energy_2d
