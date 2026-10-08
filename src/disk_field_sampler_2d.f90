! Interior-only interpolating moving least squares on Cartesian/annular disks.
! All six polynomials of total degree <=2 are reproduced. This is not a wall
! boundary condition: the reflected transport imposes specular scattering.
module disk_field_sampler_2d
  use he3_kinds, only: rk
  use cartesian_mesh_2d, only: cartesian_mesh_2d_t
  use spinful_state_2d, only: spinful_state_2d_t
  implicit none
  private
  public :: sample_disk_state, make_disk_stencil, disk_stencil_t
  type :: disk_stencil_t
    integer :: count=0, point(81)=0
    real(rk) :: weight(81)=0
    logical :: valid=.false.
  end type
contains
  subroutine make_disk_stencil(mesh,x,y,stencil)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    real(rk), intent(in) :: x,y
    type(disk_stencil_t), intent(out) :: stencil
    real(rk) :: h,radius,xy(2),dx,dy,d2,w(81),q(81,6),rr(6,6),v(81),z(6),b(6),norm
    integer :: ix,iy,i,j,k,n,cx,cy,pass
    stencil=disk_stencil_t()
    if(.not.mesh%cylinder%enabled) error stop 'disk sampler requires cylinder'
    if(.not.allocated(mesh%disk_xy).and..not.mesh%is_uniform()) &
      error stop 'disk sampler needs uniform Cartesian cylinder mesh'
    radius=mesh%cylinder%radius
    if(hypot(x,y)>radius+128*epsilon(1._rk)*max(1._rk,radius)) return
    n=0; q=0; rr=0; w=0
    if(allocated(mesh%disk_xy)) then
      call annular_candidates(mesh,x,y,stencil,n,q,w)
      if(stencil%valid) return
    else
    h=mesh%x_spacing()
    if(abs(h-mesh%y_spacing())>128*epsilon(h)*max(1._rk,h)) &
      error stop 'disk sampler needs equal x/y spacing'
    cx=nint((x-mesh%x_minimum)/h); cy=nint((y-mesh%y_minimum)/h)
    n=0; q=0; rr=0; w=0
    do iy=max(0,cy-4),min(mesh%y_cell_count(),cy+4)
      do ix=max(0,cx-4),min(mesh%x_cell_count(),cx+4)
        xy=[mesh%x_coordinate(ix),mesh%y_coordinate(iy)]
        if(norm2(xy)>radius+128*epsilon(radius)*max(1._rk,radius)) cycle
        dx=(xy(1)-x)/h; dy=(xy(2)-y)/h; d2=dx*dx+dy*dy
        if(d2>=3.5_rk**2) cycle
        if(d2<1.e-24_rk) then
          stencil%count=1; stencil%point(1)=mesh%point_index(ix,iy)
          stencil%weight(1)=1; stencil%valid=.true.; return
        end if
        n=n+1; stencil%point(n)=mesh%point_index(ix,iy)
        ! Square root of the compact, singular MLS weight: nodal interpolation
        ! in the limit, taper smoothly to zero at the support boundary.
        ! Bound the weight ratio near a node (notably the displaced wall
        ! target); an unbounded ratio destroys QR rank numerically.
        w(n)=(1-d2/3.5_rk**2)**2/max(d2,1.e-8_rk)
        q(n,:)=[1._rk,dx,dy,dx*dx,dx*dy,dy*dy]
      end do
    end do
    end if
    if(n<6) return
    w(:n)=w(:n)/maxval(w(:n))
    do j=1,6
      q(:n,j)=w(:n)*q(:n,j)
    end do
    ! Reorthogonalized modified Gram-Schmidt avoids normal equations.
    do j=1,6
      v=q(:,j)
      do pass=1,2
        do k=1,j-1
          norm=dot_product(q(:n,k),v(:n)); rr(k,j)=rr(k,j)+norm
          v(:n)=v(:n)-norm*q(:n,k)
        end do
      end do
      rr(j,j)=norm2(v(:n))
      if(rr(j,j)<1.e-13_rk) return
      q(:n,j)=v(:n)/rr(j,j)
    end do
    b=0; b(1)=1; z=0
    ! R^T z=e_1; interpolation weights = W^(1/2) Q z.
    do j=1,6
      z(j)=(b(j)-dot_product(rr(:j-1,j),z(:j-1)))/rr(j,j)
    end do
    stencil%count=n
    stencil%weight(:n)=w(:n)*matmul(q(:n,:),z)
    stencil%valid=.true.
  end subroutine

  subroutine annular_candidates(mesh,x,y,stencil,n,q,w)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    real(rk), intent(in) :: x,y
    type(disk_stencil_t), intent(inout) :: stencil
    integer, intent(out) :: n
    real(rk), intent(out) :: q(81,6),w(81)
    real(rk) :: r,theta,pi,hr,ht,er(2),et(2),d(2),dx,dy,d2
    real(rk) :: fraction,h0,h1,t,dr,da,taper,angle
    integer :: lo,hi,mid,j,k,m,c,p,first,last,nr
    n=0; q=0; w=0
    r=hypot(x,y); pi=acos(-1._rk); theta=atan2(y,x)
    nr=ubound(mesh%ring_radius,1); lo=0; hi=nr
    do while(hi-lo>1)
      mid=(hi+lo)/2
      if(mesh%ring_radius(mid)<r) then
        lo=mid
      else
        hi=mid
      end if
    end do
    fraction=max(0._rk,min(1._rk,(r-mesh%ring_radius(lo))/(mesh%ring_radius(hi)-mesh%ring_radius(lo))))
    t=lo+fraction
    h0=(mesh%ring_radius(lo+1)-mesh%ring_radius(max(0,lo-1)))/real(lo+1-max(0,lo-1),rk)
    h1=(mesh%ring_radius(min(nr,hi+1))-mesh%ring_radius(hi-1))/real(min(nr,hi+1)-hi+1,rk)
    hr=(1-fraction)*h0+fraction*h1
    ht=max(hr,2*pi*((1-fraction)*mesh%ring_radius(lo)/mesh%ring_count(lo)+ &
      fraction*mesh%ring_radius(hi)/mesh%ring_count(hi)))
    er=[cos(theta),sin(theta)]; et=[-er(2),er(1)]
    ! Seven neighbouring rings, at most nine points per ring. The radial and
    ! tangential scales are independent; work does not grow with total nodes.
    first=max(0,lo-3); last=min(nr,lo+3)
    do j=first,last
      dr=(j-t)/3
      if(abs(dr)>=1) cycle
      m=mesh%ring_count(j); c=floor(theta*m/(2*pi))
      do k=-min(4,(m-1)/2),min(4,m/2)
        p=mesh%ring_start(j)+modulo(c+k,m)
        angle=2*pi*modulo(c+k,m)/m-theta
        da=atan2(sin(angle),cos(angle))*m/(8*pi)
        if(j==0.or.m<=8) da=0 ! Full small rings: no angular stencil switching.
        if(abs(da)>=1) cycle
        d=mesh%disk_xy(:,p)-[x,y]
        dx=dot_product(d,er)/hr; dy=dot_product(d,et)/ht; d2=dx*dx+dy*dy
        if(d2<1.e-24_rk) then
          stencil%count=1; stencil%point(1)=p; stencil%weight(1)=1
          stencil%valid=.true.; return
        end if
        n=n+1; stencil%point(n)=p
        ! Vanishing support weights avoid jumps when a ring or angular
        ! neighbour enters/leaves the stencil. Metric scales are continuous.
        taper=(1-dr*dr)**2*(1-da*da)**2
        w(n)=taper*exp(-d2/16)/max(d2,1.e-8_rk)
        q(n,:)=[1._rk,dx,dy,dx*dx,dx*dy,dy*dy]
      end do
    end do
  end subroutine

  subroutine sample_disk_state(mesh,state,x,y,gap,current,inside)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    real(rk), intent(in) :: x,y
    complex(rk), intent(out) :: gap(3,3)
    real(rk), intent(out) :: current(3)
    logical, intent(out) :: inside
    type(disk_stencil_t) :: stencil
    integer :: i,p
    call make_disk_stencil(mesh,x,y,stencil)
    inside=stencil%valid; gap=0; current=0
    if(.not.inside) return
    do i=1,stencil%count
      p=stencil%point(i)
      gap=gap+stencil%weight(i)*state%order_parameter(:,:,p)
      current=current+stencil%weight(i)*state%current_mean_field(:,p)
    end do
  end subroutine
end module
