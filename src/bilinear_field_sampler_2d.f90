module bilinear_field_sampler_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t
  use spinful_state_2d, only : spinful_state_2d_t, spin_components_2d, &
                               orbital_components_2d, mean_field_components_2d
  implicit none
  private

  integer, parameter :: bilinear_corners = 4

  type, public :: bilinear_stencil_2d_t
    logical :: inside = .false.
    integer :: x_cell = -1
    integer :: y_cell = -1
    integer :: count = bilinear_corners
    integer :: point(16) = 0
    real(rk) :: weight(16) = 0.0_rk
  end type bilinear_stencil_2d_t

  public :: make_bilinear_stencil_2d
  public :: make_quadratic_stencil_2d
  public :: sample_spinful_state_2d
  public :: sample_spinful_stencil_2d
  public :: sample_pair_potential_2d
  public :: sample_pair_potential_stencil_2d

contains

  pure subroutine make_bilinear_stencil_2d(mesh, x, y, stencil)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    real(rk), intent(in) :: x, y
    type(bilinear_stencil_2d_t), intent(out) :: stencil

    real(rk) :: local_x, local_y, sample_x, sample_y
    real(rk) :: tolerance_x, tolerance_y
    integer :: cell_x, cell_y

    if (.not. mesh%is_valid()) error stop "cannot sample an invalid Cartesian mesh"

    stencil = bilinear_stencil_2d_t()
    tolerance_x = 64.0_rk * epsilon(1.0_rk) * &
                  max(1.0_rk, abs(mesh%x_minimum), abs(mesh%x_maximum))
    tolerance_y = 64.0_rk * epsilon(1.0_rk) * &
                  max(1.0_rk, abs(mesh%y_minimum), abs(mesh%y_maximum))

    if (x < mesh%x_minimum - tolerance_x .or. &
        x > mesh%x_maximum + tolerance_x .or. &
        y < mesh%y_minimum - tolerance_y .or. &
        y > mesh%y_maximum + tolerance_y) return

    ! Clamp only round-off-sized excursions at an explicitly included boundary.
    sample_x = min(mesh%x_maximum, max(mesh%x_minimum, x))
    sample_y = min(mesh%y_maximum, max(mesh%y_minimum, y))

    cell_x = mesh%x_cell_index(sample_x)
    cell_y = mesh%y_cell_index(sample_y)

    local_x = (sample_x - mesh%x_coordinate(cell_x)) / &
              mesh%x_cell_spacing(cell_x)
    local_y = (sample_y - mesh%y_coordinate(cell_y)) / &
              mesh%y_cell_spacing(cell_y)
    local_x = min(1.0_rk, max(0.0_rk, local_x))
    local_y = min(1.0_rk, max(0.0_rk, local_y))

    stencil%inside = .true.
    stencil%x_cell = cell_x
    stencil%y_cell = cell_y
    stencil%point(1:4) = [mesh%point_index(cell_x, cell_y), &
                     mesh%point_index(cell_x + 1, cell_y), &
                     mesh%point_index(cell_x, cell_y + 1), &
                     mesh%point_index(cell_x + 1, cell_y + 1)]
    stencil%weight(1:4) = [(1.0_rk - local_x) * (1.0_rk - local_y), &
                       local_x * (1.0_rk - local_y), &
                      (1.0_rk - local_x) * local_y, &
                       local_x * local_y]
  end subroutine make_bilinear_stencil_2d

  pure subroutine make_quadratic_stencil_2d(mesh, x, y, stencil)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    real(rk), intent(in) :: x, y
    type(bilinear_stencil_2d_t), intent(out) :: stencil
    integer :: ix(4), iy(4), nx, ny, i, j, k
    real(rk) :: wx(4), wy(4)
    call make_bilinear_stencil_2d(mesh, x, y, stencil)
    if (.not. stencil%inside) return
    if (mesh%x_point_count() < 3 .or. mesh%y_point_count() < 3) &
      error stop 'quadratic interpolation needs three nodes on each axis'
    call axis_weights(.true., min(mesh%x_maximum,max(mesh%x_minimum,x)), &
                      stencil%x_cell, ix, wx, nx)
    call axis_weights(.false., min(mesh%y_maximum,max(mesh%y_minimum,y)), &
                      stencil%y_cell, iy, wy, ny)
    stencil%point = 0
    stencil%weight = 0.0_rk
    stencil%count = nx*ny
    k = 0
    do j = 1, ny
      do i = 1, nx
        k = k+1
        stencil%point(k) = mesh%point_index(ix(i),iy(j))
        stencil%weight(k) = wx(i)*wy(j)
      end do
    end do
  contains
    pure subroutine axis_weights(along_x, value, cell, nodes, weights, n)
      logical, intent(in) :: along_x
      real(rk), intent(in) :: value
      integer, intent(in) :: cell
      integer, intent(out) :: nodes(4), n
      real(rk), intent(out) :: weights(4)
      integer :: last, first, left, right, start, a, b, index
      real(rk) :: coordinates(4), w
      last = mesh%y_cell_count()
      if (along_x) last = mesh%x_cell_count()
      ! Average the two bracketing quadratics in an interior cell.
      ! Both interpolate its endpoints: continuous at grid nodes, no
      ! midpoint stencil-switch discontinuity. Use one-sided at boundaries.
      left = max(0,cell-1)
      right = min(cell,last-2)
      first = left
      n = right+3-first
      weights = 0.0_rk
      nodes = 0
      do a = 1,n
        nodes(a) = first+a-1
        if (along_x) then
          coordinates(a) = mesh%x_coordinate(nodes(a))
        else
          coordinates(a) = mesh%y_coordinate(nodes(a))
        end if
      end do
      do start = left,right
        do a = start,start+2
          index = a-first+1
          w = 1.0_rk / real(right-left+1,rk)
          do b = start,start+2
            if (a == b) cycle
            w = w*(value-coordinates(b-first+1)) / &
                  (coordinates(index)-coordinates(b-first+1))
          end do
          weights(index) = weights(index)+w
        end do
      end do
    end subroutine axis_weights
  end subroutine make_quadratic_stencil_2d


  pure subroutine sample_spinful_state_2d(mesh, state, x, y, order_parameter, &
                                          current_mean_field, inside)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    real(rk), intent(in) :: x, y
    complex(rk), intent(out) :: order_parameter(spin_components_2d, &
                                                 orbital_components_2d)
    real(rk), intent(out) :: current_mean_field(mean_field_components_2d)
    logical, intent(out) :: inside

    type(bilinear_stencil_2d_t) :: local_stencil

    if (allocated(state%radial_point)) then
      call state%sample_radial(x,y,order_parameter,current_mean_field,inside)
      return
    end if

    if (.not. state%is_valid_for(mesh)) &
      error stop "cannot sample a 2D state that does not match its mesh"

    call make_bilinear_stencil_2d(mesh, x, y, local_stencil)
    call sample_spinful_stencil_2d(state, local_stencil, order_parameter, &
                                   current_mean_field, inside)
  end subroutine sample_spinful_state_2d


  pure subroutine sample_spinful_stencil_2d(state, stencil, order_parameter, &
                                            current_mean_field, inside)
    type(spinful_state_2d_t), intent(in) :: state
    type(bilinear_stencil_2d_t), intent(in) :: stencil
    complex(rk), intent(out) :: order_parameter(spin_components_2d, &
                                                 orbital_components_2d)
    real(rk), intent(out) :: current_mean_field(mean_field_components_2d)
    logical, intent(out) :: inside

    integer :: corner, point

    order_parameter = cmplx(0.0_rk, 0.0_rk, kind=rk)
    current_mean_field = 0.0_rk
    inside = stencil%inside
    if (.not. inside) return

    if (.not. allocated(state%order_parameter) .or. &
        .not. allocated(state%current_mean_field)) &
      error stop "cannot apply a bilinear stencil to an empty 2D state"
    if (any(stencil%point(1:stencil%count) < 1) .or. &
        any(stencil%point(1:stencil%count) > state%point_count())) &
      error stop "bilinear stencil point index is outside the 2D state"

    do corner = 1, stencil%count
      point = stencil%point(corner)
      order_parameter = order_parameter + &
        cmplx(stencil%weight(corner), 0.0_rk, kind=rk) * &
        state%order_parameter(:, :, point)
      current_mean_field = current_mean_field + &
        stencil%weight(corner) * state%current_mean_field(:, point)
    end do
  end subroutine sample_spinful_stencil_2d


  pure subroutine sample_pair_potential_2d(mesh, state, x, y, momentum, &
                                           pair_potential, current_mean_field, &
                                           inside)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    real(rk), intent(in) :: x, y
    real(rk), intent(in) :: momentum(orbital_components_2d)
    complex(rk), intent(out) :: pair_potential(spin_components_2d)
    real(rk), intent(out) :: current_mean_field(mean_field_components_2d)
    logical, intent(out) :: inside

    complex(rk) :: order_parameter(spin_components_2d, orbital_components_2d)
    integer :: orbital

    call sample_spinful_state_2d(mesh, state, x, y, order_parameter, &
                                 current_mean_field, inside)

    pair_potential = cmplx(0.0_rk, 0.0_rk, kind=rk)
    if (.not. inside) return

    do orbital = 1, orbital_components_2d
      pair_potential = pair_potential + &
        cmplx(momentum(orbital), 0.0_rk, kind=rk) * &
        order_parameter(:, orbital)
    end do
  end subroutine sample_pair_potential_2d


  pure subroutine sample_pair_potential_stencil_2d(state, stencil, momentum, &
                                                   pair_potential, &
                                                   current_mean_field, inside)
    type(spinful_state_2d_t), intent(in) :: state
    type(bilinear_stencil_2d_t), intent(in) :: stencil
    real(rk), intent(in) :: momentum(orbital_components_2d)
    complex(rk), intent(out) :: pair_potential(spin_components_2d)
    real(rk), intent(out) :: current_mean_field(mean_field_components_2d)
    logical, intent(out) :: inside

    complex(rk) :: order_parameter(spin_components_2d, orbital_components_2d)
    integer :: orbital

    call sample_spinful_stencil_2d(state, stencil, order_parameter, &
                                   current_mean_field, inside)

    pair_potential = cmplx(0.0_rk, 0.0_rk, kind=rk)
    if (.not. inside) return
    do orbital = 1, orbital_components_2d
      pair_potential = pair_potential + &
        cmplx(momentum(orbital), 0.0_rk, kind=rk) * &
        order_parameter(:, orbital)
    end do
  end subroutine sample_pair_potential_stencil_2d

end module bilinear_field_sampler_2d
