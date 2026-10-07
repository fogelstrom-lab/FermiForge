module historical_core_seed_2d
  use he3_kinds, only : rk
  use order_parameter_basis, only : axial_harmonics_to_cartesian
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t
  use spinful_state_2d, only : spinful_state_2d_t, &
    allocate_spinful_state_2d
  implicit none
  private

  type, public :: historical_core_seed_report_t
    character(len=3) :: kind = ""
    real(rk) :: core_length_scale = 0.0_rk
    real(rk) :: center_pair_amplitude = 0.0_rk
    real(rk) :: outer_pair_amplitude = 0.0_rk
  end type historical_core_seed_report_t

  public :: initialize_historical_core_seed_2d
  public :: initialize_localized_harmonic_seed_2d

contains

  subroutine initialize_localized_harmonic_seed_2d( &
      mesh, mode, bulk_gap, reduced_temperature, core_radius, amplitude, state)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    character(len=*), intent(in) :: mode
    real(rk), intent(in) :: bulk_gap, reduced_temperature, core_radius, amplitude
    type(spinful_state_2d_t), intent(out) :: state
    complex(rk) :: harmonic(3,3), addition(3,3)
    real(rk) :: coordinate(2), radius_squared, envelope
    integer :: spin, orbital, point

    if (.not. ieee_is_finite(core_radius) .or. .not. ieee_is_finite(amplitude)) &
      error stop 'localized seed controls must be finite'
    if (core_radius<=0.0_rk.or.amplitude<=0.0_rk) &
      error stop 'localized seed radius and amplitude must be positive'
    select case(trim(mode))
    case('localized_0plus')
      spin=2; orbital=1
    case('localized_plus0')
      spin=1; orbital=2
    case('localized_0minus')
      spin=2; orbital=3
    case('localized_minus0')
      spin=3; orbital=2
    case default
      error stop 'unknown localized harmonic seed'
    end select
    ! Preserve the phase-wound diagonal B background, but no historical tail.
    call initialize_historical_core_seed_2d( &
      mesh,'nop',bulk_gap,reduced_temperature,0.0_rk,state)
    do point=1,mesh%point_count()
      coordinate=mesh%point_coordinate(point)
      radius_squared=sum(coordinate**2)
      if(radius_squared>=core_radius**2) cycle
      ! Compact C1 envelope: value and radial derivative vanish at the edge.
      envelope=(1.0_rk-radius_squared/core_radius**2)**2
      harmonic=cmplx(0.0_rk,0.0_rk,rk)
      harmonic(spin,orbital)=cmplx(amplitude*bulk_gap*envelope,0.0_rk,rk)
      call axial_harmonics_to_cartesian(harmonic,addition)
      state%order_parameter(:,:,point)=state%order_parameter(:,:,point)+addition
    end do
  end subroutine initialize_localized_harmonic_seed_2d

  subroutine initialize_historical_core_seed_2d( &
      mesh, kind, bulk_gap, reduced_temperature, feedback_parameter, &
      state, report)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    character(len=*), intent(in) :: kind
    real(rk), intent(in) :: bulk_gap, reduced_temperature
    real(rk), intent(in) :: feedback_parameter
    type(spinful_state_2d_t), intent(out) :: state
    type(historical_core_seed_report_t), intent(out), optional :: report

    complex(rk) :: diagonal, vortex_phase, core_harmonic(3,3), core_cartesian(3,3)
    real(rk) :: angular_cosine, angular_sine, ar, core_fill
    real(rk) :: coordinate(2), core_length_scale, phase, radius, x, y
    integer :: point

    if (.not. mesh%is_valid()) &
      error stop "cannot construct a historical seed on an invalid mesh"
    if (bulk_gap <= 0.0_rk .or. reduced_temperature < 0.0_rk .or. &
        reduced_temperature >= 1.0_rk) &
      error stop "invalid historical core seed gap or temperature"
    if (trim(kind) /= "nop" .and. trim(kind) /= "aop" .and. &
        trim(kind) /= "dop" .and. trim(kind) /= "qop") &
      error stop "core seed kind must be nop, aop, dop, or qop"
    if (trim(kind) == "dop" .and. 1.0_rk + feedback_parameter <= 0.0_rk) &
      error stop "historical dop seed has a nonpositive feedback scale"

    ! This is a literal modern transcription of nop/aop/dop in
    ! new_src/how_I_initialise_the_different_cores.f.  The symbol aa0 in
    ! that source is the already transformed feedback parameter.  Only dop
    ! multiplies its radial scale by (1+aa0).
    core_length_scale = 1.0_rk / sqrt(1.0_rk - reduced_temperature**2)
    if (trim(kind) == "dop") &
      core_length_scale = core_length_scale * (1.0_rk + feedback_parameter)

    call allocate_spinful_state_2d(mesh, state)
    do point = 1, mesh%point_count()
      coordinate = mesh%point_coordinate(point)
      x = coordinate(1)
      y = coordinate(2)
      radius = sqrt(x * x + y * y)
      ar = radius / (3.0_rk * core_length_scale)
      phase = 0.0_rk
      diagonal = cmplx(0.0_rk, 0.0_rk, kind=rk)
      if (radius > 0.0_rk) then
        phase = atan2(y, x)
        vortex_phase = exp(cmplx(0.0_rk, phase, kind=rk))
        diagonal = cmplx(bulk_gap * tanh(ar), 0.0_rk, kind=rk) * &
          vortex_phase
      end if
      state%order_parameter(1, 1, point) = diagonal
      state%order_parameter(2, 2, point) = diagonal
      state%order_parameter(3, 3, point) = diagonal

      if (trim(kind) == "nop") cycle
      core_fill = historical_core_fill(ar)
      select case (trim(kind))
      case ("qop")
        ! JLTP 116 (1999), Fig. 1: non-winding C_0- and C_-0.
        ! Indices are (+,0,-). Do NOT apply the axisymmetric phase law:
        ! it would attach exp(2 i phi) and remove the intended finite core.
        ! This is a new seed, not a transcription of the historical aop.
        core_harmonic = cmplx(0.0_rk,0.0_rk,rk)
        core_harmonic(2,3) = sqrt(2.0_rk)*bulk_gap*core_fill
        core_harmonic(3,2) = -core_harmonic(2,3)
        call axial_harmonics_to_cartesian(core_harmonic,core_cartesian)
        state%order_parameter(:,:,point) = state%order_parameter(:,:,point) + core_cartesian
      case ("aop")
        angular_cosine = cos(2.0_rk * phase)
        angular_sine = sin(2.0_rk * phase)
        state%order_parameter(3, 1, point) = cmplx( &
          0.5_rk * bulk_gap * core_fill * (angular_cosine + 1.0_rk), &
          0.5_rk * bulk_gap * core_fill * angular_sine, kind=rk)
        state%order_parameter(1, 3, point) = &
          -state%order_parameter(3, 1, point)
        state%order_parameter(3, 2, point) = cmplx( &
          0.5_rk * bulk_gap * core_fill * angular_sine, &
          0.5_rk * bulk_gap * core_fill * (1.0_rk - angular_cosine), &
          kind=rk)
        state%order_parameter(2, 3, point) = &
          -state%order_parameter(3, 2, point)
      case ("dop")
        state%order_parameter(3, 1, point) = cmplx( &
          bulk_gap * core_fill * cos(phase)**2, 0.0_rk, kind=rk)
        state%order_parameter(1, 3, point) = &
          -state%order_parameter(3, 1, point)
      end select
    end do

    if (present(report)) call measure_historical_seed( &
      mesh, state, trim(kind), core_length_scale, report)
  end subroutine initialize_historical_core_seed_2d


  pure real(rk) function historical_core_fill(ar) result(value)
    real(rk), intent(in) :: ar

    if (abs(ar) <= sqrt(epsilon(1.0_rk))) then
      value = 1.0_rk
    else
      value = tanh(ar) / ar
    end if
  end function historical_core_fill


  subroutine measure_historical_seed( &
      mesh, state, kind, core_length_scale, report)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    character(len=*), intent(in) :: kind
    real(rk), intent(in) :: core_length_scale
    type(historical_core_seed_report_t), intent(out) :: report

    report%kind = kind
    report%core_length_scale = core_length_scale
    report%center_pair_amplitude = pair_amplitude_at(0.0_rk, 0.0_rk)
    report%outer_pair_amplitude = pair_amplitude_at( &
      min(abs(mesh%x_minimum), abs(mesh%x_maximum)), 0.0_rk)

  contains

    real(rk) function pair_amplitude_at(x, y) result(amplitude)
      real(rk), intent(in) :: x, y

      real(rk) :: coordinate(2), distance, minimum_distance
      integer :: point, selected

      selected = 1
      minimum_distance = huge(1.0_rk)
      do point = 1, mesh%point_count()
        coordinate = mesh%point_coordinate(point)
        distance = sum((coordinate - [x, y])**2)
        if (distance < minimum_distance) then
          selected = point
          minimum_distance = distance
        end if
      end do
      amplitude = sqrt(sum(abs(state%order_parameter(:, :, selected))**2) / &
        3.0_rk)
    end function pair_amplitude_at

  end subroutine measure_historical_seed

end module historical_core_seed_2d
