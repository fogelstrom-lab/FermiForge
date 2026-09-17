module radial_mesh
  use he3_kinds, only : rk
  implicit none
  private

  type, public :: radial_mesh_t
    real(rk), allocatable :: r(:)
  contains
    procedure :: cell_count => radial_cell_count
    procedure :: point_count => radial_point_count
    procedure :: outer_radius => radial_outer_radius
    procedure :: minimum_spacing => radial_minimum_spacing
    procedure :: maximum_spacing => radial_maximum_spacing
    procedure :: is_valid => radial_is_valid
  end type radial_mesh_t

  public :: make_uniform_radial_mesh
  public :: make_boundary_refined_radial_mesh
  public :: refine_radial_mesh
  public :: locate_radial_cell

contains

  subroutine make_uniform_radial_mesh(radius, number_of_cells, mesh)
    real(rk), intent(in) :: radius
    integer, intent(in) :: number_of_cells
    type(radial_mesh_t), intent(out) :: mesh

    integer :: i

    if (radius <= 0.0_rk) error stop "radial mesh radius must be positive"
    if (number_of_cells < 1) error stop "radial mesh needs at least one cell"

    allocate(mesh%r(0:number_of_cells))
    do i = 0, number_of_cells
      mesh%r(i) = radius * real(i, rk) / real(number_of_cells, rk)
    end do
  end subroutine make_uniform_radial_mesh


  subroutine make_boundary_refined_radial_mesh(radius, core_width, &
                                                surface_width, fine_spacing, &
                                                bulk_spacing, mesh)
    real(rk), intent(in) :: radius
    real(rk), intent(in) :: core_width
    real(rk), intent(in) :: surface_width
    real(rk), intent(in) :: fine_spacing
    real(rk), intent(in) :: bulk_spacing
    type(radial_mesh_t), intent(out) :: mesh

    logical, allocatable :: marked(:)
    real(rk) :: left, right, width
    integer :: i, initial_cells

    if (radius <= 0.0_rk) error stop "radial mesh radius must be positive"
    if (core_width < 0.0_rk .or. surface_width < 0.0_rk) &
      error stop "refined-region widths cannot be negative"
    if (fine_spacing <= 0.0_rk .or. bulk_spacing <= 0.0_rk) &
      error stop "radial mesh spacings must be positive"
    if (fine_spacing > bulk_spacing) &
      error stop "fine spacing cannot exceed bulk spacing"

    initial_cells = ceiling(radius / bulk_spacing)
    call make_uniform_radial_mesh(radius, initial_cells, mesh)

    do
      allocate(marked(0:mesh%cell_count() - 1), source=.false.)
      do i = 0, mesh%cell_count() - 1
        left = mesh%r(i)
        right = mesh%r(i + 1)
        width = right - left
        if (left < core_width .or. right > radius - surface_width) then
          marked(i) = width > fine_spacing * (1.0_rk + 32.0_rk * epsilon(1.0_rk))
        end if
      end do

      if (.not. any(marked)) then
        deallocate(marked)
        exit
      end if

      call refine_radial_mesh(mesh, marked)
      deallocate(marked)
    end do

    ! Enforce a two-to-one limit between adjacent cell widths. This is the
    ! one-dimensional analogue of the balance rule used by block AMR meshes.
    do
      allocate(marked(0:mesh%cell_count() - 1), source=.false.)
      do i = 0, mesh%cell_count() - 2
        left = mesh%r(i + 1) - mesh%r(i)
        right = mesh%r(i + 2) - mesh%r(i + 1)
        if (left > 2.0_rk * right * (1.0_rk + 32.0_rk * epsilon(1.0_rk))) &
          marked(i) = .true.
        if (right > 2.0_rk * left * (1.0_rk + 32.0_rk * epsilon(1.0_rk))) &
          marked(i + 1) = .true.
      end do

      if (.not. any(marked)) then
        deallocate(marked)
        exit
      end if

      call refine_radial_mesh(mesh, marked)
      deallocate(marked)
    end do
  end subroutine make_boundary_refined_radial_mesh


  subroutine refine_radial_mesh(mesh, marked)
    type(radial_mesh_t), intent(inout) :: mesh
    logical, intent(in) :: marked(0:)

    real(rk), allocatable :: refined_nodes(:)
    integer :: i, j, old_cells, new_cells

    if (.not. mesh%is_valid()) error stop "cannot refine an invalid radial mesh"

    old_cells = mesh%cell_count()
    if (size(marked) /= old_cells) &
      error stop "one radial refinement flag is required per cell"

    new_cells = old_cells + count(marked)
    allocate(refined_nodes(0:new_cells))

    j = 0
    refined_nodes(j) = mesh%r(0)
    do i = 0, old_cells - 1
      if (marked(i)) then
        j = j + 1
        refined_nodes(j) = 0.5_rk * (mesh%r(i) + mesh%r(i + 1))
      end if
      j = j + 1
      refined_nodes(j) = mesh%r(i + 1)
    end do

    call move_alloc(refined_nodes, mesh%r)
  end subroutine refine_radial_mesh


  pure integer function locate_radial_cell(mesh, radius) result(cell)
    type(radial_mesh_t), intent(in) :: mesh
    real(rk), intent(in) :: radius

    integer :: lower, middle, upper

    upper = mesh%cell_count()
    if (radius <= mesh%r(0)) then
      cell = 0
      return
    end if
    if (radius >= mesh%r(upper)) then
      cell = upper - 1
      return
    end if

    lower = 0
    do while (upper - lower > 1)
      middle = (lower + upper) / 2
      if (radius < mesh%r(middle)) then
        upper = middle
      else
        lower = middle
      end if
    end do
    cell = lower
  end function locate_radial_cell


  pure integer function radial_cell_count(self) result(number_of_cells)
    class(radial_mesh_t), intent(in) :: self

    if (allocated(self%r)) then
      number_of_cells = size(self%r) - 1
    else
      number_of_cells = 0
    end if
  end function radial_cell_count


  pure integer function radial_point_count(self) result(number_of_points)
    class(radial_mesh_t), intent(in) :: self

    if (allocated(self%r)) then
      number_of_points = size(self%r)
    else
      number_of_points = 0
    end if
  end function radial_point_count


  pure real(rk) function radial_outer_radius(self) result(radius)
    class(radial_mesh_t), intent(in) :: self

    if (allocated(self%r)) then
      radius = self%r(ubound(self%r, 1))
    else
      radius = 0.0_rk
    end if
  end function radial_outer_radius


  pure real(rk) function radial_minimum_spacing(self) result(spacing)
    class(radial_mesh_t), intent(in) :: self
    integer :: n

    n = self%cell_count()
    if (n > 0) then
      spacing = minval(self%r(1:n) - self%r(0:n - 1))
    else
      spacing = 0.0_rk
    end if
  end function radial_minimum_spacing


  pure real(rk) function radial_maximum_spacing(self) result(spacing)
    class(radial_mesh_t), intent(in) :: self
    integer :: n

    n = self%cell_count()
    if (n > 0) then
      spacing = maxval(self%r(1:n) - self%r(0:n - 1))
    else
      spacing = 0.0_rk
    end if
  end function radial_maximum_spacing


  pure logical function radial_is_valid(self) result(valid)
    class(radial_mesh_t), intent(in) :: self
    integer :: n

    valid = .false.
    if (.not. allocated(self%r)) return
    if (lbound(self%r, 1) /= 0) return

    n = self%cell_count()
    if (n < 1) return
    if (abs(self%r(0)) > 32.0_rk * epsilon(1.0_rk)) return
    if (any(self%r(1:n) <= self%r(0:n - 1))) return

    valid = .true.
  end function radial_is_valid

end module radial_mesh
