module he3_quadrature
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use he3_kinds, only : rk
  implicit none
  private

  real(rk), parameter :: pi = acos(-1.0_rk)

  type, public :: angular_quadrature_3d_t
    integer :: azimuth_count = 0
    integer :: polar_count = 0
    real(rk), allocatable :: momentum(:, :)
    real(rk), allocatable :: weight(:)
  contains
    procedure :: direction_count => angular_direction_count
    procedure :: is_valid => angular_quadrature_is_valid
    procedure :: second_moment => angular_second_moment
  end type angular_quadrature_3d_t

  type, public :: ozaki_quadrature_t
    real(rk) :: temperature = 0.0_rk
    integer :: cutoff_index = 0
    real(rk) :: gap_prefactor = 0.0_rk
    real(rk), allocatable :: pole(:)
    real(rk), allocatable :: residue(:)
  contains
    procedure :: pole_count => ozaki_pole_count
    procedure :: is_valid => ozaki_quadrature_is_valid
  end type ozaki_quadrature_t

  public :: make_legacy_angular_quadrature
  public :: read_legacy_gauss_table
  public :: make_ozaki_quadrature
  public :: read_legacy_ozaki_table

contains

  subroutine make_legacy_angular_quadrature( &
      azimuth_count, polar_node, raw_polar_weight, quadrature)
    integer, intent(in) :: azimuth_count
    real(rk), intent(in) :: polar_node(:), raw_polar_weight(:)
    type(angular_quadrature_3d_t), intent(out) :: quadrature

    real(rk) :: azimuth, planar_magnitude
    integer :: azimuth_index, direction, polar_index

    if (azimuth_count < 1) &
      error stop "angular quadrature needs at least one azimuth"
    if (size(polar_node) < 1 .or. &
        size(raw_polar_weight) /= size(polar_node)) &
      error stop "angular quadrature polar arrays are inconsistent"
    if (any(abs(polar_node) > 1.0_rk + 64.0_rk * epsilon(1.0_rk))) &
      error stop "angular quadrature polar node lies outside [-1,1]"
    if (any(raw_polar_weight <= 0.0_rk)) &
      error stop "angular quadrature weights must be positive"

    quadrature%azimuth_count = azimuth_count
    quadrature%polar_count = size(polar_node)
    allocate(quadrature%momentum(3, azimuth_count * size(polar_node)))
    allocate(quadrature%weight(azimuth_count * size(polar_node)))

    direction = 0
    do azimuth_index = 1, azimuth_count
      azimuth = (real(azimuth_index - 1, rk) + 0.5_rk) * &
                pi / real(azimuth_count, rk)
      do polar_index = 1, size(polar_node)
        direction = direction + 1
        planar_magnitude = sqrt(max(0.0_rk, &
          1.0_rk - polar_node(polar_index)**2))
        quadrature%momentum(:, direction) = [ &
          cos(azimuth) * planar_magnitude, &
          sin(azimuth) * planar_magnitude, polar_node(polar_index)]
        ! new_src divides the Gauss-Legendre weights by two and uses a
        ! uniform azimuthal weight 1/tmax.
        quadrature%weight(direction) = &
          0.5_rk * raw_polar_weight(polar_index) / real(azimuth_count, rk)
      end do
    end do
  end subroutine make_legacy_angular_quadrature


  subroutine read_legacy_gauss_table( &
      path, azimuth_count, polar_count, quadrature)
    character(len=*), intent(in) :: path
    integer, intent(in) :: azimuth_count, polar_count
    type(angular_quadrature_3d_t), intent(out) :: quadrature

    real(rk), allocatable :: node(:), raw_weight(:)
    integer :: input_status, polar_index, unit

    if (polar_count < 1) error stop "Gauss table polar count must be positive"
    allocate(node(polar_count), raw_weight(polar_count))
    open(newunit=unit, file=trim(path), status="old", action="read", &
         iostat=input_status)
    if (input_status /= 0) error stop "could not open legacy Gauss table"
    do polar_index = 1, polar_count
      read(unit, *, iostat=input_status) node(polar_index), raw_weight(polar_index)
      if (input_status /= 0) error stop "could not read legacy Gauss table"
    end do
    close(unit)
    call make_legacy_angular_quadrature( &
      azimuth_count, node, raw_weight, quadrature)
  end subroutine read_legacy_gauss_table


  subroutine make_ozaki_quadrature( &
      temperature, cutoff_index, pole, residue, quadrature)
    real(rk), intent(in) :: temperature
    integer, intent(in) :: cutoff_index
    real(rk), intent(in) :: pole(:), residue(:)
    type(ozaki_quadrature_t), intent(out) :: quadrature

    real(rk) :: denominator

    if (temperature <= 0.0_rk) &
      error stop "Ozaki temperature must be positive"
    if (cutoff_index < 0) error stop "Ozaki cutoff index cannot be negative"
    if (size(pole) < 1 .or. size(residue) /= size(pole)) &
      error stop "Ozaki pole and residue arrays are inconsistent"
    if (any(pole <= 0.0_rk) .or. any(residue <= 0.0_rk)) &
      error stop "Ozaki poles and residues must be positive"

    quadrature%temperature = temperature
    quadrature%cutoff_index = cutoff_index
    allocate(quadrature%pole, source=pole)
    allocate(quadrature%residue, source=residue)
    denominator = log(temperature) + &
      sum(temperature * residue / pole)
    if (abs(denominator) <= 64.0_rk * epsilon(1.0_rk)) &
      error stop "Ozaki gap-prefactor denominator is singular"
    quadrature%gap_prefactor = temperature / denominator
  end subroutine make_ozaki_quadrature


  subroutine read_legacy_ozaki_table(path, quadrature)
    character(len=*), intent(in) :: path
    type(ozaki_quadrature_t), intent(out) :: quadrature

    real(rk), allocatable :: pole(:), residue(:)
    real(rk) :: cutoff, temperature
    integer :: count, input_status, pole_index, unit

    open(newunit=unit, file=trim(path), status="old", action="read", &
         iostat=input_status)
    if (input_status /= 0) error stop "could not open legacy Ozaki table"
    read(unit, *, iostat=input_status) count, temperature, cutoff
    if (input_status /= 0 .or. count < 1) &
      error stop "could not read legacy Ozaki header"
    allocate(pole(count), residue(count))
    do pole_index = 1, count
      read(unit, *, iostat=input_status) pole(pole_index), residue(pole_index)
      if (input_status /= 0) error stop "could not read legacy Ozaki table"
    end do
    close(unit)
    call make_ozaki_quadrature( &
      temperature, int(cutoff), pole, residue, quadrature)
  end subroutine read_legacy_ozaki_table


  pure integer function angular_direction_count(self) result(count)
    class(angular_quadrature_3d_t), intent(in) :: self

    if (allocated(self%weight)) then
      count = size(self%weight)
    else
      count = 0
    end if
  end function angular_direction_count


  pure logical function angular_quadrature_is_valid(self) result(valid)
    class(angular_quadrature_3d_t), intent(in) :: self

    real(rk) :: norm_error
    integer :: direction

    valid = self%azimuth_count > 0 .and. self%polar_count > 0 .and. &
            allocated(self%momentum) .and. allocated(self%weight)
    if (.not. valid) return
    valid = size(self%momentum, 1) == 3 .and. &
            size(self%momentum, 2) == size(self%weight) .and. &
            size(self%weight) == self%azimuth_count * self%polar_count .and. &
            all(self%weight > 0.0_rk) .and. &
            all(ieee_is_finite(self%momentum)) .and. &
            all(ieee_is_finite(self%weight))
    if (.not. valid) return
    norm_error = 0.0_rk
    do direction = 1, self%direction_count()
      norm_error = max(norm_error, &
        abs(dot_product(self%momentum(:, direction), &
                        self%momentum(:, direction)) - 1.0_rk))
    end do
    valid = norm_error <= 256.0_rk * epsilon(1.0_rk)
  end function angular_quadrature_is_valid


  pure function angular_second_moment(self) result(moment)
    class(angular_quadrature_3d_t), intent(in) :: self
    real(rk) :: moment(3)
    integer :: component, direction

    if (.not. self%is_valid()) &
      error stop "cannot take moments of invalid angular quadrature"
    moment = 0.0_rk
    do direction = 1, self%direction_count()
      do component = 1, 3
        moment(component) = moment(component) + self%weight(direction) * &
          self%momentum(component, direction)**2
      end do
    end do
  end function angular_second_moment


  pure integer function ozaki_pole_count(self) result(count)
    class(ozaki_quadrature_t), intent(in) :: self

    if (allocated(self%pole)) then
      count = size(self%pole)
    else
      count = 0
    end if
  end function ozaki_pole_count


  pure logical function ozaki_quadrature_is_valid(self) result(valid)
    class(ozaki_quadrature_t), intent(in) :: self

    valid = self%temperature > 0.0_rk .and. self%cutoff_index >= 0 .and. &
            allocated(self%pole) .and. allocated(self%residue)
    if (.not. valid) return
    valid = size(self%pole) == size(self%residue) .and. &
            size(self%pole) >= 1 .and. all(self%pole > 0.0_rk) .and. &
            all(self%residue > 0.0_rk) .and. &
            all(ieee_is_finite(self%pole)) .and. &
            all(ieee_is_finite(self%residue)) .and. &
            ieee_is_finite(self%gap_prefactor)
  end function ozaki_quadrature_is_valid

end module he3_quadrature
