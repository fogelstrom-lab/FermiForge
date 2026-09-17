module legacy_anderson_mixing
  use he3_kinds, only : rk
  implicit none
  private

  real(rk), parameter :: initial_mixing = 0.01_rk
  real(rk), parameter :: stochastic_residual_marker = 1.0e29_rk

  type, public :: anderson_report_t
    integer :: iteration = 0
    integer :: history_size = 0
    integer :: discarded_vectors = 0
    real(rk) :: max_residual = huge(1.0_rk)
    real(rk) :: residual_norm = huge(1.0_rk)
    real(rk) :: legacy_relative_residual = huge(1.0_rk)
    real(rk) :: mixing = initial_mixing
    logical :: converged = .false.
    logical :: used_stochastic_step = .false.
  end type anderson_report_t

  type, public :: legacy_anderson_t
    private
    integer :: number_of_values = 0
    integer :: maximum_history = 0
    integer :: history_dimension = 0
    integer :: iteration = 0
    real(rk) :: progress_threshold = 0.1_rk
    real(rk) :: maximum_mixing = 1.0_rk
    real(rk) :: mixing = initial_mixing
    real(rk) :: minimum_residual_norm_squared = 0.0_rk
    real(rk) :: stochastic_state = 2.0_rk * acos(-1.0_rk)
    real(rk), allocatable :: iterates(:, :)
    real(rk), allocatable :: residuals(:, :)
    real(rk), allocatable :: residual_norm_squared(:)
  contains
    procedure :: initialize => initialize_anderson
    procedure :: reset => reset_anderson
    procedure :: update => update_anderson
    procedure :: vector_size => anderson_vector_size
    procedure :: history_size => anderson_history_size
  end type legacy_anderson_t

contains

  subroutine initialize_anderson(self, number_of_values, maximum_history, &
                                 progress_threshold, maximum_mixing)
    class(legacy_anderson_t), intent(inout) :: self
    integer, intent(in) :: number_of_values
    integer, intent(in) :: maximum_history
    real(rk), intent(in), optional :: progress_threshold
    real(rk), intent(in), optional :: maximum_mixing

    real(rk) :: requested_maximum_mixing, requested_progress_threshold

    requested_progress_threshold = 0.1_rk
    if (present(progress_threshold)) &
      requested_progress_threshold = progress_threshold

    requested_maximum_mixing = 1.0_rk
    if (present(maximum_mixing)) requested_maximum_mixing = maximum_mixing

    if (number_of_values < 1) &
      error stop "Anderson vector size must be positive"
    if (maximum_history < 1) &
      error stop "Anderson history limit must be positive"
    if (requested_progress_threshold <= 0.0_rk .or. &
        requested_progress_threshold > 1.0_rk) &
      error stop "Anderson progress threshold must be in (0,1]"
    if (requested_maximum_mixing < initial_mixing) &
      error stop "Anderson maximum mixing cannot be smaller than 0.01"

    self%number_of_values = number_of_values
    self%maximum_history = maximum_history
    self%progress_threshold = requested_progress_threshold
    self%maximum_mixing = requested_maximum_mixing

    if (allocated(self%iterates)) deallocate(self%iterates)
    if (allocated(self%residuals)) deallocate(self%residuals)
    if (allocated(self%residual_norm_squared)) &
      deallocate(self%residual_norm_squared)

    ! The original citerat uses slots 0:mmax.  The extra slot holds the newly
    ! evaluated point before the least useful history vector is discarded.
    allocate(self%iterates(number_of_values, maximum_history + 1), source=0.0_rk)
    allocate(self%residuals(number_of_values, maximum_history + 1), source=0.0_rk)
    allocate(self%residual_norm_squared(maximum_history + 1), source=0.0_rk)

    call self%reset()
  end subroutine initialize_anderson


  subroutine reset_anderson(self)
    class(legacy_anderson_t), intent(inout) :: self

    if (.not. allocated(self%iterates)) &
      error stop "Anderson accelerator must be initialized before reset"

    self%history_dimension = 0
    self%iteration = 0
    self%mixing = initial_mixing
    self%minimum_residual_norm_squared = 0.0_rk
    self%stochastic_state = 2.0_rk * acos(-1.0_rk)
    self%iterates = 0.0_rk
    self%residuals = 0.0_rk
    self%residual_norm_squared = 0.0_rk
  end subroutine reset_anderson


  subroutine update_anderson(self, current, mapped, tolerance, next, report)
    class(legacy_anderson_t), intent(inout) :: self
    real(rk), intent(in) :: current(:)
    real(rk), intent(in) :: mapped(:)
    real(rk), intent(in) :: tolerance
    real(rk), intent(out) :: next(:)
    type(anderson_report_t), intent(out) :: report

    real(rk), allocatable :: coefficients(:), displacement(:), &
                             minimum_residual(:)
    real(rk) :: mapped_norm_squared, mixing_ratio, residual_norm_squared
    integer :: discarded, m, slot
    logical :: solved

    call require_compatible_vectors(self, current, mapped, next)
    if (tolerance < 0.0_rk) &
      error stop "Anderson convergence tolerance cannot be negative"

    self%iteration = self%iteration + 1
    m = self%history_dimension
    slot = m + 1

    self%iterates(:, slot) = current
    self%residuals(:, slot) = mapped - current
    residual_norm_squared = dot_product(self%residuals(:, slot), &
                                        self%residuals(:, slot))
    self%residual_norm_squared(slot) = residual_norm_squared
    mapped_norm_squared = dot_product(mapped, mapped)

    report = anderson_report_t()
    report%iteration = self%iteration
    report%max_residual = maxval(abs(self%residuals(:, slot)))
    report%residual_norm = sqrt(max(residual_norm_squared, 0.0_rk))
    if (mapped_norm_squared < tolerance) mapped_norm_squared = 1.0_rk
    report%legacy_relative_residual = &
      sqrt(max(residual_norm_squared / mapped_norm_squared, 0.0_rk))
    report%mixing = self%mixing

    ! citerat restores the current iterate, rather than G(current), when its
    ! unnormalised maximum-component residual reaches the requested tolerance.
    if (report%max_residual <= tolerance) then
      next = current
      report%converged = .true.
      report%history_size = m
      return
    end if

    discarded = 0
    if (m == self%maximum_history) discarded = 1
    if (m > 1) then
      if (self%minimum_residual_norm_squared / residual_norm_squared < &
            self%progress_threshold .or. &
          self%residual_norm_squared(1) / residual_norm_squared < 1.0_rk) &
        discarded = m / 2
    end if

    mixing_ratio = min(1.0_rk, &
      sqrt(max(self%minimum_residual_norm_squared / residual_norm_squared, &
               0.0_rk)))

    if (discarded > 0) then
      call discard_largest_residuals(self, m, discarded)
      m = m - discarded
    end if
    report%discarded_vectors = discarded

    if (m == 0) then
      next = self%iterates(:, 1) + self%mixing * self%residuals(:, 1)
      self%minimum_residual_norm_squared = residual_norm_squared
      self%history_dimension = 1
      report%history_size = self%history_dimension
      return
    end if

    ! Previous calls leave the retained vectors in descending residual norm.
    ! Insert the newly evaluated point into that ordering, as rearrange did.
    call insert_current_by_descending_residual(self, m)
    residual_norm_squared = self%residual_norm_squared(m + 1)

    allocate(coefficients(m), displacement(self%number_of_values), &
             minimum_residual(self%number_of_values))
    call residual_minimizing_coefficients(self, m, residual_norm_squared, &
                                          coefficients, solved)

    if (.not. solved) then
      call stochastic_step(self, m, next)
      self%minimum_residual_norm_squared = stochastic_residual_marker
      self%history_dimension = m + 1
      report%history_size = self%history_dimension
      report%mixing = self%mixing
      report%used_stochastic_step = .true.
      return
    end if

    displacement = 0.0_rk
    do slot = 1, m
      displacement = displacement + &
        (self%iterates(:, slot) - self%iterates(:, m + 1)) * &
        coefficients(slot)
    end do

    self%mixing = mixing_ratio * &
      (sqrt(max(dot_product(displacement, displacement) / &
                residual_norm_squared, 0.0_rk)) + self%mixing) * 0.5_rk
    self%mixing = min(self%maximum_mixing, self%mixing)
    self%mixing = max(initial_mixing, self%mixing)

    minimum_residual = self%residuals(:, m + 1)
    do slot = 1, m
      minimum_residual = minimum_residual + &
        (self%residuals(:, slot) - self%residuals(:, m + 1)) * &
        coefficients(slot)
    end do
    self%minimum_residual_norm_squared = &
      dot_product(minimum_residual, minimum_residual)

    next = self%iterates(:, m + 1) + displacement + &
           self%mixing * minimum_residual

    discarded = 0
    if (self%minimum_residual_norm_squared / residual_norm_squared < &
        self%progress_threshold) then
      discarded = m / 2
      if (discarded > 0) then
        call discard_largest_residuals(self, m, discarded)
        m = m - discarded
      end if
    end if

    self%history_dimension = m + 1
    report%discarded_vectors = report%discarded_vectors + discarded
    report%history_size = self%history_dimension
    report%mixing = self%mixing
  end subroutine update_anderson


  subroutine residual_minimizing_coefficients(self, m, reference_norm_squared, &
                                               coefficients, solved)
    class(legacy_anderson_t), intent(in) :: self
    integer, intent(in) :: m
    real(rk), intent(in) :: reference_norm_squared
    real(rk), intent(out) :: coefficients(m)
    logical, intent(out) :: solved

    real(rk), allocatable :: differences(:, :), matrix(:, :), right_hand_side(:)
    real(rk) :: difference_norm_ratio
    integer :: i

    allocate(differences(self%number_of_values, m), matrix(m, m), &
             right_hand_side(m))

    do i = 1, m
      differences(:, i) = self%residuals(:, i) - self%residuals(:, m + 1)
    end do

    matrix = matmul(transpose(differences), differences) / reference_norm_squared
    right_hand_side = -matmul(transpose(differences), &
                              self%residuals(:, m + 1)) / reference_norm_squared

    if (m == 1) then
      difference_norm_ratio = sqrt(max(matrix(1, 1), 0.0_rk))
      if (difference_norm_ratio < 1.0e-15_rk) then
        solved = .false.
        coefficients = 0.0_rk
      else
        solved = .true.
        coefficients(1) = right_hand_side(1) / matrix(1, 1)
      end if
      return
    end if

    call solve_dense_system(matrix, right_hand_side, coefficients, solved)
  end subroutine residual_minimizing_coefficients


  subroutine solve_dense_system(matrix, right_hand_side, solution, solved)
    real(rk), intent(in) :: matrix(:, :)
    real(rk), intent(in) :: right_hand_side(:)
    real(rk), intent(out) :: solution(:)
    logical, intent(out) :: solved

    real(rk), allocatable :: a(:, :), b(:), row(:)
    real(rk) :: factor, matrix_scale, pivot_floor, value
    integer :: i, k, n, pivot

    n = size(right_hand_side)
    if (size(matrix, 1) /= n .or. size(matrix, 2) /= n .or. &
        size(solution) /= n) error stop "invalid Anderson subspace system"

    allocate(a(n, n), b(n), row(n))
    a = matrix
    b = right_hand_side
    matrix_scale = maxval(abs(a))
    pivot_floor = epsilon(1.0_rk) * max(1.0_rk, matrix_scale) * real(n, rk)

    solved = .false.
    solution = 0.0_rk
    do k = 1, n
      pivot = k - 1 + maxloc(abs(a(k:n, k)), dim=1)
      if (abs(a(pivot, k)) <= pivot_floor) return

      if (pivot /= k) then
        row = a(k, :)
        a(k, :) = a(pivot, :)
        a(pivot, :) = row
        value = b(k)
        b(k) = b(pivot)
        b(pivot) = value
      end if

      do i = k + 1, n
        factor = a(i, k) / a(k, k)
        a(i, k:n) = a(i, k:n) - factor * a(k, k:n)
        b(i) = b(i) - factor * b(k)
      end do
    end do

    do i = n, 1, -1
      if (abs(a(i, i)) <= pivot_floor) return
      solution(i) = (b(i) - dot_product(a(i, i + 1:n), &
                                        solution(i + 1:n))) / a(i, i)
    end do
    solved = .true.
  end subroutine solve_dense_system


  subroutine insert_current_by_descending_residual(self, m)
    class(legacy_anderson_t), intent(inout) :: self
    integer, intent(in) :: m

    real(rk), allocatable :: vector(:)
    real(rk) :: norm_squared
    integer :: j

    allocate(vector(self%number_of_values))
    do j = m, 1, -1
      if (self%residual_norm_squared(j + 1) <= &
          self%residual_norm_squared(j)) exit

      vector = self%iterates(:, j + 1)
      self%iterates(:, j + 1) = self%iterates(:, j)
      self%iterates(:, j) = vector

      vector = self%residuals(:, j + 1)
      self%residuals(:, j + 1) = self%residuals(:, j)
      self%residuals(:, j) = vector

      norm_squared = self%residual_norm_squared(j + 1)
      self%residual_norm_squared(j + 1) = self%residual_norm_squared(j)
      self%residual_norm_squared(j) = norm_squared
    end do
  end subroutine insert_current_by_descending_residual


  subroutine discard_largest_residuals(self, m, number_to_discard)
    class(legacy_anderson_t), intent(inout) :: self
    integer, intent(in) :: m
    integer, intent(in) :: number_to_discard

    integer :: retained_points

    if (number_to_discard < 1 .or. number_to_discard > m) &
      error stop "invalid number of Anderson history vectors to discard"

    retained_points = m + 1 - number_to_discard
    self%iterates(:, 1:retained_points) = &
      self%iterates(:, number_to_discard + 1:m + 1)
    self%residuals(:, 1:retained_points) = &
      self%residuals(:, number_to_discard + 1:m + 1)
    self%residual_norm_squared(1:retained_points) = &
      self%residual_norm_squared(number_to_discard + 1:m + 1)
  end subroutine discard_largest_residuals


  subroutine stochastic_step(self, m, next)
    class(legacy_anderson_t), intent(inout) :: self
    integer, intent(in) :: m
    real(rk), intent(out) :: next(:)

    real(rk) :: argument, pi2
    integer :: i

    pi2 = 2.0_rk * acos(-1.0_rk)
    self%stochastic_state = &
      (self%stochastic_state - real(nint(self%stochastic_state), rk)) * &
      10.0_rk + pi2 * 1.0e-10_rk
    argument = pi2 * 1000.0_rk * self%stochastic_state / &
               real(self%number_of_values, rk)

    do i = 1, self%number_of_values
      next(i) = self%iterates(i, m + 1) + sin(argument * real(i, rk)) * &
                self%residuals(i, m + 1)
    end do
    self%mixing = 1.0_rk
  end subroutine stochastic_step


  subroutine require_compatible_vectors(self, current, mapped, next)
    class(legacy_anderson_t), intent(in) :: self
    real(rk), intent(in) :: current(:), mapped(:), next(:)

    if (.not. allocated(self%iterates)) &
      error stop "Anderson accelerator has not been initialized"
    if (size(current) /= self%number_of_values .or. &
        size(mapped) /= self%number_of_values .or. &
        size(next) /= self%number_of_values) &
      error stop "Anderson vector has the wrong size"
  end subroutine require_compatible_vectors


  pure integer function anderson_vector_size(self) result(number_of_values)
    class(legacy_anderson_t), intent(in) :: self

    number_of_values = self%number_of_values
  end function anderson_vector_size


  pure integer function anderson_history_size(self) result(number_of_vectors)
    class(legacy_anderson_t), intent(in) :: self

    number_of_vectors = self%history_dimension
  end function anderson_history_size

end module legacy_anderson_mixing
