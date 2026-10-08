module cartesian_mesh_2d
  use he3_kinds, only : rk
  use specular_cylinder_2d, only: cylinder_geometry_t
  implicit none
  private

  type, public :: cartesian_mesh_2d_t
    type(cylinder_geometry_t) :: cylinder
    integer :: trajectory_interpolation_order = 1
    real(rk) :: x_minimum = 0.0_rk
    real(rk) :: x_maximum = 0.0_rk
    real(rk) :: y_minimum = 0.0_rk
    real(rk) :: y_maximum = 0.0_rk
    integer :: number_of_x_cells = 0
    integer :: number_of_y_cells = 0
    real(rk), allocatable :: x_node_coordinates(:)
    real(rk), allocatable :: y_node_coordinates(:)
    ! Optional non-tensor-product disk nodes. Transport and nonlinear packing
    ! use point_count/point_coordinate; Cartesian-only algorithms must not
    ! interpret these as an x-by-y product.
    real(rk), allocatable :: disk_xy(:,:), ring_radius(:)
    integer, allocatable :: ring_start(:), ring_count(:)
  contains
    procedure :: x_cell_count => cartesian_x_cell_count
    procedure :: y_cell_count => cartesian_y_cell_count
    procedure :: x_point_count => cartesian_x_point_count
    procedure :: y_point_count => cartesian_y_point_count
    procedure :: point_count => cartesian_point_count
    procedure :: x_spacing => cartesian_x_spacing
    procedure :: y_spacing => cartesian_y_spacing
    procedure :: x_cell_spacing => cartesian_x_cell_spacing
    procedure :: y_cell_spacing => cartesian_y_cell_spacing
    procedure :: x_cell_index => cartesian_x_cell_index
    procedure :: y_cell_index => cartesian_y_cell_index
    procedure :: x_coordinate => cartesian_x_coordinate
    procedure :: y_coordinate => cartesian_y_coordinate
    procedure :: point_index => cartesian_point_index
    procedure :: point_coordinate => cartesian_point_coordinate
    procedure :: minimum_spacing => cartesian_minimum_spacing
    procedure :: maximum_spacing => cartesian_maximum_spacing
    procedure :: is_uniform => cartesian_is_uniform
    procedure :: is_valid => cartesian_is_valid
  end type cartesian_mesh_2d_t

  public :: make_uniform_cartesian_mesh
  public :: make_rectilinear_cartesian_mesh
  public :: make_symmetric_multiscale_cartesian_mesh
  public :: make_smooth_cartesian_mesh
  public :: make_extended_smooth_cartesian_mesh
  public :: make_circular_active_mask
  public :: make_annular_disk_mesh

contains

  subroutine make_annular_disk_mesh(radius,rings,stretch,tangent_spacing,mesh,radial_layout)
    real(rk), intent(in) :: radius,stretch,tangent_spacing
    integer, intent(in) :: rings
    type(cartesian_mesh_2d_t), intent(out) :: mesh
    character(len=*), optional, intent(in) :: radial_layout
    logical :: core_wall
    real(rk) :: u,theta,pi
    integer :: j,k,n,p
    core_wall=.false.
    if(present(radial_layout)) then
      if(radial_layout/='wall'.and.radial_layout/='core_wall') error stop 'unknown annular radial layout'
      core_wall=radial_layout=='core_wall'
    end if
    if(radius<=0.or.rings<4.or.stretch<0.or.stretch>4.or.tangent_spacing<=0) &
      error stop 'invalid annular disk controls'
    pi=acos(-1._rk)
    allocate(mesh%ring_radius(0:rings),mesh%ring_start(0:rings),mesh%ring_count(0:rings))
    mesh%ring_radius(0)=0; mesh%ring_count(0)=1; mesh%ring_start(0)=1
    n=1
    do j=1,rings
      u=real(j,rk)/rings
      if(stretch<1.e-8_rk) then
        mesh%ring_radius(j)=radius*u
      else if(core_wall) then
        mesh%ring_radius(j)=radius*(1+tanh(stretch*(2*u-1))/tanh(stretch))/2
      else
        mesh%ring_radius(j)=radius*(1-(exp(stretch*(1-u))-1)/(exp(stretch)-1))
      end if
      mesh%ring_count(j)=max(8,4*ceiling(2*pi*mesh%ring_radius(j)/tangent_spacing/4))
      ! Retain angular resolution of the winding in small inner rings.
      if(core_wall) mesh%ring_count(j)=max(24,mesh%ring_count(j))
      mesh%ring_start(j)=n+1; n=n+mesh%ring_count(j)
    end do
    mesh%ring_radius(rings)=radius
    allocate(mesh%disk_xy(2,n)); mesh%disk_xy(:,1)=0
    do j=1,rings
      do k=0,mesh%ring_count(j)-1
        theta=2*pi*k/mesh%ring_count(j); p=mesh%ring_start(j)+k
        mesh%disk_xy(:,p)=mesh%ring_radius(j)*[cos(theta),sin(theta)]
      end do
    end do
    mesh%x_minimum=-radius; mesh%x_maximum=radius
    mesh%y_minimum=-radius; mesh%y_maximum=radius
    mesh%cylinder%enabled=.true.; mesh%cylinder%radius=radius
  end subroutine

  subroutine make_extended_smooth_cartesian_mesh(core_width,cells,stretch,outer_width,growth,mesh)
    real(rk), intent(in) :: core_width,stretch,outer_width,growth
    integer, intent(in) :: cells
    type(cartesian_mesh_2d_t), intent(out) :: mesh
    type(cartesian_mesh_2d_t) :: core
    real(rk), allocatable :: positive(:),coordinates(:)
    real(rk) :: step, distance, factor
    integer :: n,i,j
    if (outer_width<=core_width.or.growth<1.0_rk.or.growth>2.0_rk) &
      error stop 'invalid extended smooth grid controls'
    call make_smooth_cartesian_mesh(core_width,cells,stretch,core)
    step=core%x_cell_spacing(cells-1)
    distance=outer_width-core_width
    n=1
    do
      if (growth==1.0_rk) then
        factor=real(n,rk)
      else
        factor=(growth**n-1.0_rk)/(growth-1.0_rk)
      end if
      if(step*factor>=distance) exit
      n=n+1
      if(n>10000) error stop 'too many outer grid bands'
    end do
    ! Rescale outer steps together to land exactly at the boundary, avoiding
    ! a tiny final cell. Preserve all old coordinates, including +/-core_width.
    step=distance/factor
    allocate(positive(n),coordinates(cells+1+2*n))
    positive(1)=core_width+step
    do i=2,n
      positive(i)=positive(i-1)+step*growth**(i-1)
    end do
    positive(n)=outer_width
    do i=1,n
      coordinates(i)=-positive(n-i+1)
    end do
    do j=0,cells
      coordinates(n+j+1)=core%x_coordinate(j)
    end do
    coordinates(n+cells+2:)=positive
    call make_rectilinear_cartesian_mesh(coordinates,coordinates,mesh)
  end subroutine make_extended_smooth_cartesian_mesh

  subroutine make_smooth_cartesian_mesh(half_width, cells, stretch, mesh)
    real(rk), intent(in) :: half_width, stretch
    integer, intent(in) :: cells
    type(cartesian_mesh_2d_t), intent(out) :: mesh
    real(rk), allocatable :: coordinates(:)
    real(rk) :: u
    integer :: i
    if (half_width <= 0.0_rk .or. cells < 2 .or. modulo(cells,2) /= 0 .or. &
        stretch < 0.0_rk .or. stretch > 10.0_rk) &
      error stop 'smooth grid requires positive width, even cells, stretch in [0,10]'
    allocate(coordinates(cells+1))
    do i=0,cells
      u=real(2*i-cells,rk)/real(cells,rk)
      if (stretch < 1.0e-8_rk) then
        coordinates(i+1)=half_width*u
      else
        coordinates(i+1)=half_width*sinh(stretch*u)/sinh(stretch)
      end if
    end do
    coordinates(1)=-half_width
    coordinates(cells+1)=half_width
    coordinates(cells/2+1)=0.0_rk
    call make_rectilinear_cartesian_mesh(coordinates,coordinates,mesh)
  end subroutine make_smooth_cartesian_mesh

  subroutine make_uniform_cartesian_mesh(x_minimum, x_maximum, number_of_x_cells, &
                                         y_minimum, y_maximum, number_of_y_cells, mesh)
    real(rk), intent(in) :: x_minimum, x_maximum
    real(rk), intent(in) :: y_minimum, y_maximum
    integer, intent(in) :: number_of_x_cells, number_of_y_cells
    type(cartesian_mesh_2d_t), intent(out) :: mesh

    real(rk), allocatable :: x_coordinates(:), y_coordinates(:)
    integer :: node

    if (x_maximum <= x_minimum) &
      error stop "Cartesian mesh x bounds must be strictly increasing"
    if (y_maximum <= y_minimum) &
      error stop "Cartesian mesh y bounds must be strictly increasing"
    if (number_of_x_cells < 1 .or. number_of_y_cells < 1) &
      error stop "Cartesian mesh needs at least one cell in each direction"

    allocate(x_coordinates(number_of_x_cells + 1))
    allocate(y_coordinates(number_of_y_cells + 1))
    do node = 0, number_of_x_cells
      x_coordinates(node + 1) = x_minimum + &
        (x_maximum - x_minimum) * real(node, rk) / &
        real(number_of_x_cells, rk)
    end do
    do node = 0, number_of_y_cells
      y_coordinates(node + 1) = y_minimum + &
        (y_maximum - y_minimum) * real(node, rk) / &
        real(number_of_y_cells, rk)
    end do
    call make_rectilinear_cartesian_mesh(x_coordinates, y_coordinates, mesh)
  end subroutine make_uniform_cartesian_mesh


  subroutine make_rectilinear_cartesian_mesh( &
      x_coordinates, y_coordinates, mesh)
    real(rk), intent(in) :: x_coordinates(:), y_coordinates(:)
    type(cartesian_mesh_2d_t), intent(out) :: mesh

    integer :: node

    if (size(x_coordinates) < 2 .or. size(y_coordinates) < 2) &
      error stop "rectilinear Cartesian mesh needs two nodes per axis"
    do node = 2, size(x_coordinates)
      if (x_coordinates(node) <= x_coordinates(node - 1)) &
        error stop "rectilinear x coordinates must be strictly increasing"
    end do
    do node = 2, size(y_coordinates)
      if (y_coordinates(node) <= y_coordinates(node - 1)) &
        error stop "rectilinear y coordinates must be strictly increasing"
    end do

    mesh%number_of_x_cells = size(x_coordinates) - 1
    mesh%number_of_y_cells = size(y_coordinates) - 1
    mesh%x_minimum = x_coordinates(1)
    mesh%x_maximum = x_coordinates(size(x_coordinates))
    mesh%y_minimum = y_coordinates(1)
    mesh%y_maximum = y_coordinates(size(y_coordinates))
    allocate(mesh%x_node_coordinates(0:mesh%number_of_x_cells))
    allocate(mesh%y_node_coordinates(0:mesh%number_of_y_cells))
    mesh%x_node_coordinates = x_coordinates
    mesh%y_node_coordinates = y_coordinates
  end subroutine make_rectilinear_cartesian_mesh


  subroutine make_symmetric_multiscale_cartesian_mesh( &
      half_width, fine_region_half_width, medium_region_half_width, &
      fine_spacing, medium_spacing, coarse_spacing, mesh)
    real(rk), intent(in) :: half_width, fine_region_half_width
    real(rk), intent(in) :: medium_region_half_width
    real(rk), intent(in) :: fine_spacing, medium_spacing, coarse_spacing
    type(cartesian_mesh_2d_t), intent(out) :: mesh

    real(rk), allocatable :: coordinates(:)

    if (half_width <= 0.0_rk .or. fine_region_half_width <= 0.0_rk .or. &
        medium_region_half_width <= fine_region_half_width .or. &
        half_width <= medium_region_half_width) &
      error stop "invalid symmetric multiscale region extents"
    if (fine_spacing <= 0.0_rk .or. medium_spacing < fine_spacing .or. &
        coarse_spacing < medium_spacing) &
      error stop "multiscale target spacings must be positive and ordered"

    call make_symmetric_axis_coordinates( &
      half_width, fine_region_half_width, medium_region_half_width, &
      fine_spacing, medium_spacing, coarse_spacing, coordinates)
    call make_rectilinear_cartesian_mesh(coordinates, coordinates, mesh)
  end subroutine make_symmetric_multiscale_cartesian_mesh


  subroutine make_symmetric_axis_coordinates( &
      half_width, fine_region_half_width, medium_region_half_width, &
      fine_spacing, medium_spacing, coarse_spacing, coordinates)
    real(rk), intent(in) :: half_width, fine_region_half_width
    real(rk), intent(in) :: medium_region_half_width
    real(rk), intent(in) :: fine_spacing, medium_spacing, coarse_spacing
    real(rk), allocatable, intent(out) :: coordinates(:)

    real(rk) :: lower, upper
    real(rk) :: region_end(3), target_spacing(3)
    integer :: center, count(3), index, offset, region

    region_end = [fine_region_half_width, medium_region_half_width, half_width]
    target_spacing = [fine_spacing, medium_spacing, coarse_spacing]
    lower = 0.0_rk
    do region = 1, 3
      count(region) = ceiling((region_end(region) - lower) / &
                              target_spacing(region))
      lower = region_end(region)
    end do

    center = sum(count) + 1
    allocate(coordinates(2 * sum(count) + 1))
    coordinates(center) = 0.0_rk
    lower = 0.0_rk
    offset = 0
    do region = 1, 3
      upper = region_end(region)
      do index = 1, count(region)
        offset = offset + 1
        coordinates(center + offset) = lower + (upper - lower) * &
          real(index, rk) / real(count(region), rk)
      end do
      lower = upper
    end do
    do index = 1, center - 1
      coordinates(center - index) = -coordinates(center + index)
    end do
  end subroutine make_symmetric_axis_coordinates


  subroutine make_circular_active_mask(mesh, center, radius, active)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    real(rk), intent(in) :: center(2), radius
    logical, allocatable, intent(out) :: active(:)

    real(rk) :: coordinate(2), tolerance
    integer :: point

    if (.not. mesh%is_valid()) &
      error stop "cannot mask an invalid Cartesian mesh"
    if (radius <= 0.0_rk) &
      error stop "circular active radius must be positive"
    tolerance = 64.0_rk * epsilon(1.0_rk) * &
      max(1.0_rk, radius, maxval(abs(center)))
    allocate(active(mesh%point_count()))
    do point = 1, mesh%point_count()
      coordinate = mesh%point_coordinate(point)
      active(point) = sum((coordinate - center)**2) <= &
        (radius + tolerance)**2
    end do
    if (.not. any(active)) error stop "circular active mask contains no points"
  end subroutine make_circular_active_mask


  pure integer function cartesian_x_cell_count(self) result(number_of_cells)
    class(cartesian_mesh_2d_t), intent(in) :: self

    number_of_cells = self%number_of_x_cells
  end function cartesian_x_cell_count


  pure integer function cartesian_y_cell_count(self) result(number_of_cells)
    class(cartesian_mesh_2d_t), intent(in) :: self

    number_of_cells = self%number_of_y_cells
  end function cartesian_y_cell_count


  pure integer function cartesian_x_point_count(self) result(number_of_points)
    class(cartesian_mesh_2d_t), intent(in) :: self

    if (self%number_of_x_cells > 0) then
      number_of_points = self%number_of_x_cells + 1
    else
      number_of_points = 0
    end if
  end function cartesian_x_point_count


  pure integer function cartesian_y_point_count(self) result(number_of_points)
    class(cartesian_mesh_2d_t), intent(in) :: self

    if (self%number_of_y_cells > 0) then
      number_of_points = self%number_of_y_cells + 1
    else
      number_of_points = 0
    end if
  end function cartesian_y_point_count


  pure integer function cartesian_point_count(self) result(number_of_points)
    class(cartesian_mesh_2d_t), intent(in) :: self

    number_of_points = self%x_point_count() * self%y_point_count()
    if(allocated(self%disk_xy)) number_of_points=size(self%disk_xy,2)
  end function cartesian_point_count


  pure real(rk) function cartesian_x_spacing(self) result(spacing)
    class(cartesian_mesh_2d_t), intent(in) :: self

    if (self%number_of_x_cells > 0) then
      spacing = (self%x_maximum - self%x_minimum) / &
                real(self%number_of_x_cells, rk)
    else
      spacing = 0.0_rk
    end if
  end function cartesian_x_spacing


  pure real(rk) function cartesian_y_spacing(self) result(spacing)
    class(cartesian_mesh_2d_t), intent(in) :: self

    if (self%number_of_y_cells > 0) then
      spacing = (self%y_maximum - self%y_minimum) / &
                real(self%number_of_y_cells, rk)
    else
      spacing = 0.0_rk
    end if
  end function cartesian_y_spacing


  pure real(rk) function cartesian_x_cell_spacing(self, cell) result(spacing)
    class(cartesian_mesh_2d_t), intent(in) :: self
    integer, intent(in) :: cell

    if (cell < 0 .or. cell >= self%number_of_x_cells) &
      error stop "Cartesian x-cell index is out of range"
    spacing = self%x_node_coordinates(cell + 1) - &
              self%x_node_coordinates(cell)
  end function cartesian_x_cell_spacing


  pure real(rk) function cartesian_y_cell_spacing(self, cell) result(spacing)
    class(cartesian_mesh_2d_t), intent(in) :: self
    integer, intent(in) :: cell

    if (cell < 0 .or. cell >= self%number_of_y_cells) &
      error stop "Cartesian y-cell index is out of range"
    spacing = self%y_node_coordinates(cell + 1) - &
              self%y_node_coordinates(cell)
  end function cartesian_y_cell_spacing


  pure integer function cartesian_x_cell_index(self, coordinate) result(cell)
    class(cartesian_mesh_2d_t), intent(in) :: self
    real(rk), intent(in) :: coordinate

    cell = axis_cell_index(self%x_node_coordinates, coordinate)
  end function cartesian_x_cell_index


  pure integer function cartesian_y_cell_index(self, coordinate) result(cell)
    class(cartesian_mesh_2d_t), intent(in) :: self
    real(rk), intent(in) :: coordinate

    cell = axis_cell_index(self%y_node_coordinates, coordinate)
  end function cartesian_y_cell_index


  pure integer function axis_cell_index(coordinates, coordinate) result(cell)
    real(rk), intent(in) :: coordinates(0:)
    real(rk), intent(in) :: coordinate

    integer :: left, middle, right

    if (coordinate < coordinates(0) .or. &
        coordinate > coordinates(ubound(coordinates, 1))) &
      error stop "Cartesian coordinate lies outside its mesh axis"
    if (coordinate >= coordinates(ubound(coordinates, 1))) then
      cell = ubound(coordinates, 1) - 1
      return
    end if
    left = 0
    right = ubound(coordinates, 1)
    do while (right - left > 1)
      middle = (left + right) / 2
      if (coordinate >= coordinates(middle)) then
        left = middle
      else
        right = middle
      end if
    end do
    cell = left
  end function axis_cell_index


  pure real(rk) function cartesian_x_coordinate(self, node) result(coordinate)
    class(cartesian_mesh_2d_t), intent(in) :: self
    integer, intent(in) :: node

    if (node < 0 .or. node > self%number_of_x_cells) &
      error stop "Cartesian x-node index is out of range"
    coordinate = self%x_node_coordinates(node)
  end function cartesian_x_coordinate


  pure real(rk) function cartesian_y_coordinate(self, node) result(coordinate)
    class(cartesian_mesh_2d_t), intent(in) :: self
    integer, intent(in) :: node

    if (node < 0 .or. node > self%number_of_y_cells) &
      error stop "Cartesian y-node index is out of range"
    coordinate = self%y_node_coordinates(node)
  end function cartesian_y_coordinate


  pure integer function cartesian_point_index(self, x_node, y_node) result(point)
    class(cartesian_mesh_2d_t), intent(in) :: self
    integer, intent(in) :: x_node, y_node

    if (x_node < 0 .or. x_node > self%number_of_x_cells) &
      error stop "Cartesian x-node index is out of range"
    if (y_node < 0 .or. y_node > self%number_of_y_cells) &
      error stop "Cartesian y-node index is out of range"

    ! The x index varies fastest. This gives one deterministic flat point
    ! numbering for MPI work assignment, interpolation, and solver packing.
    point = 1 + x_node + self%x_point_count() * y_node
  end function cartesian_point_index


  pure function cartesian_point_coordinate(self, point) result(coordinate)
    class(cartesian_mesh_2d_t), intent(in) :: self
    integer, intent(in) :: point
    real(rk) :: coordinate(2)

    integer :: point_zero_based, x_node, y_node

    if (point < 1 .or. point > self%point_count()) &
      error stop "Cartesian flat point index is out of range"
    if(allocated(self%disk_xy)) then
      coordinate=self%disk_xy(:,point)
      return
    end if
    point_zero_based = point - 1
    x_node = modulo(point_zero_based, self%x_point_count())
    y_node = point_zero_based / self%x_point_count()
    coordinate = [self%x_coordinate(x_node), self%y_coordinate(y_node)]
  end function cartesian_point_coordinate


  pure real(rk) function cartesian_minimum_spacing(self) result(spacing)
    class(cartesian_mesh_2d_t), intent(in) :: self

    if(allocated(self%disk_xy)) then
      spacing=minval(self%ring_radius(1:)-self%ring_radius(:size(self%ring_radius)-2))
      spacing=min(spacing,minval(2*self%ring_radius(1:)*sin(acos(-1._rk)/self%ring_count(1:))))
      return
    end if

    if (.not. self%is_valid()) then
      spacing = 0.0_rk
      return
    end if
    spacing = min( &
      minval(self%x_node_coordinates(1:) - self%x_node_coordinates(:self%number_of_x_cells - 1)), &
      minval(self%y_node_coordinates(1:) - self%y_node_coordinates(:self%number_of_y_cells - 1)))
  end function cartesian_minimum_spacing


  pure real(rk) function cartesian_maximum_spacing(self) result(spacing)
    class(cartesian_mesh_2d_t), intent(in) :: self

    if(allocated(self%disk_xy)) then
      spacing=maxval(self%ring_radius(1:)-self%ring_radius(:size(self%ring_radius)-2))
      spacing=max(spacing,maxval(2*self%ring_radius(1:)*sin(acos(-1._rk)/self%ring_count(1:))))
      return
    end if

    if (.not. self%is_valid()) then
      spacing = 0.0_rk
      return
    end if
    spacing = max( &
      maxval(self%x_node_coordinates(1:) - self%x_node_coordinates(:self%number_of_x_cells - 1)), &
      maxval(self%y_node_coordinates(1:) - self%y_node_coordinates(:self%number_of_y_cells - 1)))
  end function cartesian_maximum_spacing


  pure logical function cartesian_is_uniform(self) result(uniform)
    class(cartesian_mesh_2d_t), intent(in) :: self

    real(rk) :: tolerance

    if(allocated(self%disk_xy)) then
      uniform=.false.; return
    end if
    uniform = self%is_valid()
    if (.not. uniform) return
    tolerance = 64.0_rk * epsilon(1.0_rk) * &
      max(1.0_rk, abs(self%x_minimum), abs(self%x_maximum), &
          abs(self%y_minimum), abs(self%y_maximum))
    uniform = maxval(abs( &
      self%x_node_coordinates(1:) - &
      self%x_node_coordinates(:self%number_of_x_cells - 1) - &
      self%x_spacing())) <= tolerance .and. maxval(abs( &
      self%y_node_coordinates(1:) - &
      self%y_node_coordinates(:self%number_of_y_cells - 1) - &
      self%y_spacing())) <= tolerance
  end function cartesian_is_uniform


  pure logical function cartesian_is_valid(self) result(valid)
    class(cartesian_mesh_2d_t), intent(in) :: self

    if(allocated(self%disk_xy)) then
      valid=size(self%disk_xy,1)==2.and.size(self%disk_xy,2)>6.and.allocated(self%ring_radius)
      return
    end if

    valid = self%number_of_x_cells > 0 .and. &
            self%number_of_y_cells > 0 .and. &
            self%x_maximum > self%x_minimum .and. &
            self%y_maximum > self%y_minimum .and. &
            allocated(self%x_node_coordinates) .and. &
            allocated(self%y_node_coordinates)
    if (.not. valid) return
    valid = lbound(self%x_node_coordinates, 1) == 0 .and. &
            ubound(self%x_node_coordinates, 1) == self%number_of_x_cells .and. &
            lbound(self%y_node_coordinates, 1) == 0 .and. &
            ubound(self%y_node_coordinates, 1) == self%number_of_y_cells .and. &
            all(self%x_node_coordinates(1:) > &
                self%x_node_coordinates(:self%number_of_x_cells - 1)) .and. &
            all(self%y_node_coordinates(1:) > &
                self%y_node_coordinates(:self%number_of_y_cells - 1))
  end function cartesian_is_valid

end module cartesian_mesh_2d
