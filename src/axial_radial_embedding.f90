module axial_radial_embedding
  use he3_kinds, only : rk
  use radial_mesh, only : radial_mesh_t, locate_radial_cell
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t
  use spinful_state_2d, only : spinful_state_2d_t
  use order_parameter_basis, only : projection_plus, projection_zero, &
                                    projection_minus, &
                                    cartesian_to_axial_harmonics, &
                                    axial_harmonics_to_transport_cartesian
  implicit none
  private

  integer, parameter :: projection_value(3) = [1, 0, -1]
  complex(rk), parameter :: complex_zero = cmplx(0.0_rk, 0.0_rk, kind=rk)

  type, public :: axial_radial_profile_t
    complex(rk), allocatable :: harmonic(:, :, :)
    real(rk), allocatable :: azimuthal_mean_field(:)
  contains
    procedure :: point_count => axial_profile_point_count
    procedure :: is_valid_for => axial_profile_is_valid_for
  end type axial_radial_profile_t

  public :: allocate_axial_radial_profile
  public :: set_axial_profile_from_axis_cartesian
  public :: sample_axial_radial_profile
  public :: embed_axial_radial_profile_2d

contains

  subroutine allocate_axial_radial_profile(mesh, profile)
    type(radial_mesh_t), intent(in) :: mesh
    type(axial_radial_profile_t), intent(out) :: profile

    integer :: last_point

    if (.not. mesh%is_valid()) &
      error stop "cannot allocate an axial profile on an invalid radial mesh"
    if (mesh%cell_count() < 2) &
      error stop "source-faithful axial interpolation needs at least two cells"

    last_point = mesh%point_count() - 1
    allocate(profile%harmonic(3, 3, 0:last_point), source=complex_zero)
    allocate(profile%azimuthal_mean_field(0:last_point), source=0.0_rk)
  end subroutine allocate_axial_radial_profile


  subroutine set_axial_profile_from_axis_cartesian(mesh, axis_cartesian, &
                                                   azimuthal_mean_field, profile)
    type(radial_mesh_t), intent(in) :: mesh
    complex(rk), intent(in) :: axis_cartesian(:, :, 0:)
    real(rk), intent(in) :: azimuthal_mean_field(0:)
    type(axial_radial_profile_t), intent(out) :: profile

    integer :: point

    if (size(axis_cartesian, 1) /= 3 .or. size(axis_cartesian, 2) /= 3) &
      error stop "axis Cartesian order parameter must be 3 by 3 at every point"
    if (size(axis_cartesian, 3) /= mesh%point_count() .or. &
        size(azimuthal_mean_field) /= mesh%point_count()) &
      error stop "axis Cartesian profile and radial mesh sizes differ"

    call allocate_axial_radial_profile(mesh, profile)
    do point = 0, mesh%point_count() - 1
      call cartesian_to_axial_harmonics(axis_cartesian(:, :, point), &
                                        profile%harmonic(:, :, point))
    end do
    profile%azimuthal_mean_field = azimuthal_mean_field
  end subroutine set_axial_profile_from_axis_cartesian


  pure subroutine sample_axial_radial_profile(mesh, profile, x, y, winding, &
                                              order_parameter, &
                                              current_mean_field, inside)
    type(radial_mesh_t), intent(in) :: mesh
    type(axial_radial_profile_t), intent(in) :: profile
    real(rk), intent(in) :: x, y, winding
    complex(rk), intent(out) :: order_parameter(3, 3)
    real(rk), intent(out) :: current_mean_field(3)
    logical, intent(out) :: inside

    complex(rk) :: phased_harmonic(3, 3)
    real(rk) :: angle, azimuthal_value, radius, tolerance

    if (.not. profile%is_valid_for(mesh)) &
      error stop "cannot sample an axial profile that does not match its mesh"

    order_parameter = complex_zero
    current_mean_field = 0.0_rk
    radius = sqrt(x * x + y * y)
    tolerance = 64.0_rk * epsilon(1.0_rk) * &
                max(1.0_rk, mesh%outer_radius())
    inside = radius <= mesh%outer_radius() + tolerance
    if (.not. inside) return
    radius = min(radius, mesh%outer_radius())

    if (radius > tolerance) then
      angle = atan2(y, x)
    else
      angle = 0.0_rk
    end if

    call interpolate_phased_harmonics(mesh, profile, radius, angle, winding, &
                                      phased_harmonic)
    call axial_harmonics_to_transport_cartesian(phased_harmonic, order_parameter)
    call interpolate_azimuthal_mean_field(mesh, profile, radius, azimuthal_value)

    if (radius > tolerance) then
      current_mean_field(1) = -sin(angle) * azimuthal_value
      current_mean_field(2) = cos(angle) * azimuthal_value
    end if
  end subroutine sample_axial_radial_profile


  subroutine embed_axial_radial_profile_2d(radial_grid, profile, winding, &
                                           cartesian_grid, state)
    type(radial_mesh_t), intent(in) :: radial_grid
    type(axial_radial_profile_t), intent(in) :: profile
    real(rk), intent(in) :: winding
    type(cartesian_mesh_2d_t), intent(in) :: cartesian_grid
    type(spinful_state_2d_t), intent(inout) :: state

    complex(rk) :: order_parameter(3, 3)
    real(rk) :: current_mean_field(3), x, y
    logical :: inside
    integer :: point, x_node, y_node

    if (.not. profile%is_valid_for(radial_grid)) &
      error stop "cannot embed an invalid axial radial profile"
    if (.not. state%is_valid_for(cartesian_grid)) &
      error stop "2D state does not match the Cartesian embedding mesh"

    do y_node = 0, cartesian_grid%y_cell_count()
      y = cartesian_grid%y_coordinate(y_node)
      do x_node = 0, cartesian_grid%x_cell_count()
        x = cartesian_grid%x_coordinate(x_node)
        point = cartesian_grid%point_index(x_node, y_node)
        call sample_axial_radial_profile(radial_grid, profile, x, y, winding, &
                                         order_parameter, current_mean_field, &
                                         inside)
        if (.not. inside) &
          error stop "Cartesian embedding mesh extends beyond radial profile"
        state%order_parameter(:, :, point) = order_parameter
        state%current_mean_field(:, point) = current_mean_field
      end do
    end do
  end subroutine embed_axial_radial_profile_2d


  pure subroutine interpolate_phased_harmonics(mesh, profile, radius, angle, &
                                               winding, phased_harmonic)
    type(radial_mesh_t), intent(in) :: mesh
    type(axial_radial_profile_t), intent(in) :: profile
    real(rk), intent(in) :: radius, angle, winding
    complex(rk), intent(out) :: phased_harmonic(3, 3)

    complex(rk) :: phase, value(3)
    real(rk) :: angular_momentum, coordinate(3), weight(3)
    integer :: cell, index(3), orbital_projection, sample, spin_projection
    logical :: reflected(3)

    call radial_source_stencil(mesh, radius, cell, index, coordinate, &
                               reflected, weight)
    phased_harmonic = complex_zero

    do spin_projection = projection_plus, projection_minus
      do orbital_projection = projection_plus, projection_minus
        angular_momentum = winding - &
          real(projection_value(spin_projection), rk) - &
          real(projection_value(orbital_projection), rk)
        do sample = 1, 3
          phase = exp(cmplx(0.0_rk, angular_momentum * angle, kind=rk))
          if (reflected(sample)) phase = phase * &
            exp(cmplx(0.0_rk, angular_momentum * acos(-1.0_rk), kind=rk))
          value(sample) = profile%harmonic(spin_projection, orbital_projection, &
                                           index(sample)) * phase
        end do
        phased_harmonic(spin_projection, orbital_projection) = &
          sum(cmplx(weight, 0.0_rk, kind=rk) * value)
      end do
    end do
  end subroutine interpolate_phased_harmonics


  pure subroutine interpolate_azimuthal_mean_field(mesh, profile, radius, value)
    type(radial_mesh_t), intent(in) :: mesh
    type(axial_radial_profile_t), intent(in) :: profile
    real(rk), intent(in) :: radius
    real(rk), intent(out) :: value

    real(rk) :: coordinate(3), sample_value(3), weight(3)
    integer :: cell, index(3), sample
    logical :: reflected(3)

    call radial_source_stencil(mesh, radius, cell, index, coordinate, &
                               reflected, weight)
    do sample = 1, 3
      sample_value(sample) = profile%azimuthal_mean_field(index(sample))
      if (reflected(sample)) sample_value(sample) = -sample_value(sample)
    end do
    value = dot_product(weight, sample_value)
  end subroutine interpolate_azimuthal_mean_field


  pure subroutine radial_source_stencil(mesh, radius, cell, index, coordinate, &
                                        reflected, weight)
    type(radial_mesh_t), intent(in) :: mesh
    real(rk), intent(in) :: radius
    integer, intent(out) :: cell, index(3)
    real(rk), intent(out) :: coordinate(3), weight(3)
    logical, intent(out) :: reflected(3)

    cell = locate_radial_cell(mesh, radius)
    reflected = .false.
    if (cell == 0) then
      index = [1, 0, 1]
      coordinate = [-mesh%r(1), mesh%r(0), mesh%r(1)]
      reflected(1) = .true.
    else
      index = [cell - 1, cell, cell + 1]
      coordinate = mesh%r(index)
    end if

    weight(1) = (radius - coordinate(2)) * (radius - coordinate(3)) / &
                ((coordinate(1) - coordinate(2)) * &
                 (coordinate(1) - coordinate(3)))
    weight(2) = (radius - coordinate(1)) * (radius - coordinate(3)) / &
                ((coordinate(2) - coordinate(1)) * &
                 (coordinate(2) - coordinate(3)))
    weight(3) = (radius - coordinate(1)) * (radius - coordinate(2)) / &
                ((coordinate(3) - coordinate(1)) * &
                 (coordinate(3) - coordinate(2)))
  end subroutine radial_source_stencil


  pure integer function axial_profile_point_count(self) result(number_of_points)
    class(axial_radial_profile_t), intent(in) :: self

    if (allocated(self%harmonic)) then
      number_of_points = size(self%harmonic, 3)
    else
      number_of_points = 0
    end if
  end function axial_profile_point_count


  pure logical function axial_profile_is_valid_for(self, mesh) result(valid)
    class(axial_radial_profile_t), intent(in) :: self
    type(radial_mesh_t), intent(in) :: mesh

    valid = mesh%is_valid() .and. mesh%cell_count() >= 2 .and. &
            allocated(self%harmonic) .and. &
            allocated(self%azimuthal_mean_field)
    if (.not. valid) return
    valid = lbound(self%harmonic, 3) == 0 .and. &
            lbound(self%azimuthal_mean_field, 1) == 0 .and. &
            size(self%harmonic, 1) == 3 .and. &
            size(self%harmonic, 2) == 3 .and. &
            size(self%harmonic, 3) == mesh%point_count() .and. &
            size(self%azimuthal_mean_field) == mesh%point_count()
  end function axial_profile_is_valid_for

end module axial_radial_embedding
