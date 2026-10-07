! Infinite cylinder: lengths in xi_0, no end caps or exterior reservoir.
module specular_cylinder_2d
  use he3_kinds, only: rk
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  type, public :: cylinder_geometry_t
    logical :: enabled=.false.
    real(rk) :: radius=0.0_rk, half_length=0.0_rk
    ! False reproduces uniform-step sampling for comparison to new_src.
    logical :: collision_aligned=.true.
  end type
  type, public :: cylinder_path_t
    real(rk), allocatable :: s(:), position(:,:), momentum(:,:)
    integer :: target=0
  end type
  public :: build_cylinder_path
contains
  subroutine build_cylinder_path(geometry,origin,momentum,maximum_step,path)
    type(cylinder_geometry_t), intent(in) :: geometry
    real(rk), intent(in) :: origin(2),momentum(3),maximum_step
    type(cylinder_path_t), intent(out) :: path
    type(cylinder_path_t) :: forward,backward
    real(rk) :: x(2),r,tol
    integer :: n,m,i
    if(.not.all(ieee_is_finite([geometry%radius,geometry%half_length,maximum_step,origin,momentum]))) &
      error stop 'nonfinite cylinder geometry'
    if(geometry%radius<=0.or.geometry%half_length<=0.or.maximum_step<=0) &
      error stop 'cylinder radius, half length and step must be positive'
    if(abs(norm2(momentum)-1._rk)>1.e-12_rk) error stop 'cylinder momentum must be unit length'
    r=norm2(origin); tol=128*epsilon(1._rk)*max(1._rk,geometry%radius)
    if(r>geometry%radius+tol) error stop 'cylinder target outside wall'
    x=origin
    ! Same interior-limit wall convention as new_src; configurable convergence
    ! is by step refinement. The sampling field still includes the wall node.
    if(r>=geometry%radius-tol) x=origin*((geometry%radius-1.e-7_rk* &
      min(maximum_step,geometry%radius))/r)
    call half_path(x,momentum,geometry,maximum_step,forward)
    call half_path(x,-momentum,geometry,maximum_step,backward)
    n=size(backward%s); m=size(forward%s)
    allocate(path%s(n+m-1),path%position(2,n+m-1),path%momentum(3,n+m-1))
    do i=1,n
      path%s(i)=-backward%s(n+1-i)
      path%position(:,i)=backward%position(:,n+1-i)
      ! Backward tracing velocity is opposite to physical momentum.
      path%momentum(:,i)=-backward%momentum(:,n+1-i)
    end do
    path%s(n:)=forward%s
    path%position(:,n:)=forward%position
    path%momentum(:,n:)=forward%momentum
    path%target=n
  end subroutine

  subroutine half_path(origin,direction,geometry,step,path)
    real(rk), intent(in) :: origin(2),direction(3),step
    type(cylinder_geometry_t), intent(in) :: geometry
    type(cylinder_path_t), intent(out) :: path
    real(rk), allocatable :: buffer(:,:),grown(:,:)
    real(rk) :: x(2),q(3),normal(2),s,next_sample,hit,a,b,c,disc,ds,tol
    integer :: used,capacity,events,sample_index
    capacity=256; allocate(buffer(6,capacity)); used=0
    x=origin; q=direction; s=0; events=0; sample_index=1
    tol=64*epsilon(1._rk)*max(1._rk,geometry%half_length)
    call append()
    next_sample=min(step,geometry%half_length)
    do while(s<geometry%half_length)
      a=dot_product(q(1:2),q(1:2))
      hit=huge(1._rk)
      if(a>tiny(1._rk)) then
        b=dot_product(x,q(1:2)); c=dot_product(x,x)-geometry%radius**2
        disc=sqrt(max(0._rk,b*b-a*c))
        if(b>0) then
          hit=max(0._rk,-c/(b+disc))
        else
          hit=max(0._rk,(-b+disc)/a)
        end if
      end if
      ds=next_sample-s
      if(hit<=ds) then
        x=x+hit*q(1:2); s=s+hit
        normal=x/norm2(x); x=geometry%radius*normal
        if(geometry%collision_aligned) call append() ! incoming wall limit
        ! Legacy zigzag stores its first exactly coincident wall sample on
        ! the incoming segment. Preserve that convention in comparison mode.
        if(.not.geometry%collision_aligned.and.events==0.and.abs(next_sample-s)<=tol) call append()
        q(1:2)=q(1:2)-2*dot_product(q(1:2),normal)*normal
        if(geometry%collision_aligned) call append() ! outgoing, same s
        events=events+1
        if(events>100000) error stop 'cylinder path has excessive grazing collisions'
        if(abs(next_sample-s)<=tol) then
          if(.not.geometry%collision_aligned.and.events>1) call append()
          sample_index=sample_index+1
          next_sample=min(step*real(sample_index,rk),geometry%half_length)
        end if
      else
        x=x+ds*q(1:2); s=next_sample
        call append()
        sample_index=sample_index+1
        next_sample=min(step*real(sample_index,rk),geometry%half_length)
      end if
    end do
    path%s=buffer(1,:used)
    path%position=buffer(2:3,:used)
    path%momentum=buffer(4:6,:used)
  contains
    subroutine append()
      if(used==capacity) then
        allocate(grown(6,2*capacity)); grown(:,:capacity)=buffer
        call move_alloc(grown,buffer); capacity=2*capacity
      end if
      used=used+1
      buffer(:,used)=[s,x,q]
    end subroutine
  end subroutine
end module
