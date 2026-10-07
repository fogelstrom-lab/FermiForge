module spinful_state_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t
  use order_parameter_basis, only : cartesian_to_axial_harmonics, &
    axial_harmonics_to_cartesian, reconstruct_axial_cartesian
  implicit none
  private

  integer, parameter, public :: spin_components_2d = 3
  integer, parameter, public :: orbital_components_2d = 3
  integer, parameter, public :: mean_field_components_2d = 3

  type, public :: spinful_state_2d_t
    ! The logical field order follows A(spin, orbital, point). The flat point
    ! index belongs to the mesh and is independent of nonlinear-solver packing.
    complex(rk), allocatable :: order_parameter(:, :, :)
    real(rk), allocatable :: current_mean_field(:, :)
    ! Optional independent +x ray. Off-ray storage is output/cache only.
    integer, allocatable :: radial_point(:)
    real(rk), allocatable :: radial_coordinate(:)
    real(rk) :: radial_winding = 1.0_rk
    integer :: radial_interpolation_order = 2
  contains
    procedure :: point_count => state_point_count
    procedure :: is_valid_for => state_is_valid_for
    procedure :: set_zero => zero_state
    procedure :: sample_radial => sample_radial_state
  end type spinful_state_2d_t

  public :: allocate_spinful_state_2d
  public :: configure_radial_symmetry, project_radial_origin

contains

  subroutine configure_radial_symmetry(mesh, state, radius, winding, active)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(inout) :: state
    real(rk), intent(in) :: radius, winding
    logical, intent(out) :: active(:)
    real(rk) :: xy(2), tol
    integer :: p, n, i
    tol=128.0_rk*epsilon(1.0_rk)*max(1.0_rk,radius)
    active=.false.
    do p=1,mesh%point_count()
      xy=mesh%point_coordinate(p)
      active(p)=abs(xy(2))<=tol .and. xy(1)>=-tol .and. xy(1)<=radius+tol
    end do
    n=count(active)
    if(n<3) error stop 'radial symmetry requires at least three +x nodes'
    if(allocated(state%radial_point)) deallocate(state%radial_point,state%radial_coordinate)
    allocate(state%radial_point(n),state%radial_coordinate(n))
    i=0
    do p=1,mesh%point_count()
      if(.not.active(p)) cycle
      i=i+1
      state%radial_point(i)=p
      xy=mesh%point_coordinate(p)
      state%radial_coordinate(i)=max(0.0_rk,xy(1))
    end do
    if(abs(state%radial_coordinate(1))>tol) error stop 'radial ray must include origin'
    if(any(state%radial_coordinate(2:)<=state%radial_coordinate(:n-1))) &
      error stop 'radial nodes must increase'
    state%radial_winding=winding
    state%radial_interpolation_order=mesh%trajectory_interpolation_order
  end subroutine configure_radial_symmetry

  subroutine project_radial_origin(state)
    type(spinful_state_2d_t), intent(inout) :: state
    complex(rk) :: h(3,3)
    integer, parameter :: m(3)=[1,0,-1]
    integer :: i,j,p
    if(.not.allocated(state%radial_point)) return
    p=state%radial_point(1)
    call cartesian_to_axial_harmonics(state%order_parameter(:,:,p),h)
    do j=1,3
      do i=1,3
        if(abs(state%radial_winding-real(m(i)+m(j),rk))>1.e-12_rk) h(i,j)=0
      end do
    end do
    call axial_harmonics_to_cartesian(h,state%order_parameter(:,:,p))
    state%current_mean_field(1:2,p)=0
  end subroutine project_radial_origin

  pure subroutine sample_radial_state(self,x,y,gap,current,inside)
    class(spinful_state_2d_t), intent(in) :: self
    real(rk), intent(in) :: x,y
    complex(rk), intent(out) :: gap(3,3)
    real(rk), intent(out) :: current(3)
    logical, intent(out) :: inside
    complex(rk) :: axis_gap(3,3), h(3,3)
    real(rk) :: r,theta,w,axis_current(3),tol
    integer :: n,cell,lo,hi,mid,left,right,start,i,j,p
    gap=0; current=0; inside=.false.
    if(.not.allocated(self%radial_point)) return
    n=size(self%radial_point)
    r=sqrt(x*x+y*y)
    tol=128*epsilon(1.0_rk)*max(1.0_rk,self%radial_coordinate(n))
    if(r>self%radial_coordinate(n)+tol) return
    r=min(r,self%radial_coordinate(n)); theta=0
    if(r>tol) theta=atan2(y,x)
    lo=1; hi=n
    do while(hi-lo>1)
      mid=(lo+hi)/2
      if(self%radial_coordinate(mid)<=r) then
        lo=mid
      else
        hi=mid
      end if
    end do
    cell=min(lo,n-1)
    axis_gap=0; axis_current=0
    if(self%radial_interpolation_order==1) then
      left=cell; right=cell
    else
      left=max(1,cell-1); right=min(cell,n-2)
    end if
    ! Average bracketing quadratics: continuous across radial cell boundaries.
    do start=left,right
      do i=start,start+self%radial_interpolation_order
        w=1.0_rk/real(right-left+1,rk)
        do j=start,start+self%radial_interpolation_order
          if(i==j) cycle
          w=w*(r-self%radial_coordinate(j))/(self%radial_coordinate(i)-self%radial_coordinate(j))
        end do
        p=self%radial_point(i)
        axis_gap=axis_gap+w*self%order_parameter(:,:,p)
        axis_current=axis_current+w*self%current_mean_field(:,p)
      end do
    end do
    call cartesian_to_axial_harmonics(axis_gap,h)
    call reconstruct_axial_cartesian(h,theta,self%radial_winding,gap)
    current(1)=cos(theta)*axis_current(1)-sin(theta)*axis_current(2)
    current(2)=sin(theta)*axis_current(1)+cos(theta)*axis_current(2)
    current(3)=axis_current(3)
    inside=.true.
  end subroutine sample_radial_state

  subroutine allocate_spinful_state_2d(mesh, state)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(out) :: state

    integer :: number_of_points

    if (.not. mesh%is_valid()) &
      error stop "cannot allocate spinful fields on an invalid Cartesian mesh"

    number_of_points = mesh%point_count()
    allocate(state%order_parameter(spin_components_2d, orbital_components_2d, &
                                   number_of_points), &
             source=cmplx(0.0_rk, 0.0_rk, kind=rk))
    allocate(state%current_mean_field(mean_field_components_2d, number_of_points), &
             source=0.0_rk)
  end subroutine allocate_spinful_state_2d


  pure integer function state_point_count(self) result(number_of_points)
    class(spinful_state_2d_t), intent(in) :: self

    if (allocated(self%order_parameter)) then
      number_of_points = size(self%order_parameter, 3)
    else
      number_of_points = 0
    end if
  end function state_point_count


  pure logical function state_is_valid_for(self, mesh) result(valid)
    class(spinful_state_2d_t), intent(in) :: self
    type(cartesian_mesh_2d_t), intent(in) :: mesh

    valid = mesh%is_valid() .and. allocated(self%order_parameter) .and. &
            allocated(self%current_mean_field)
    if (.not. valid) return

    valid = size(self%order_parameter, 1) == spin_components_2d .and. &
            size(self%order_parameter, 2) == orbital_components_2d .and. &
            size(self%order_parameter, 3) == mesh%point_count() .and. &
            size(self%current_mean_field, 1) == mean_field_components_2d .and. &
            size(self%current_mean_field, 2) == mesh%point_count()
  end function state_is_valid_for


  subroutine zero_state(self)
    class(spinful_state_2d_t), intent(inout) :: self

    if (allocated(self%order_parameter)) &
      self%order_parameter = cmplx(0.0_rk, 0.0_rk, kind=rk)
    if (allocated(self%current_mean_field)) self%current_mean_field = 0.0_rk
  end subroutine zero_state

end module spinful_state_2d
