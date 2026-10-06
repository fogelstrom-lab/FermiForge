program test_radial_cylinder_sampling
  use Global_variables
  use Interpolations
  use cylindrical_trajectories, only : trajectories
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  real :: bx(0:mx),by(0:mx),bz(0:mx),fx(0:mx),fy(0:mx),fz(0:mx)
  real :: bpx(0:mx),bpy(0:mx),bpz(0:mx),fpx(0:mx),fpy(0:mx),fpz(0:mx)
  real :: origin, sx,sy,sz,sd,factor,r,ph,tolerance
  complex :: sem(4,-mx:mx),expected(4),ses(10),dummy
  integer :: i,iy,direction,case_index
  integer, parameter :: targets(4) = [0, nx/2, nx-1, nx]

  cyl=.true.; icyl=1; Rx=4.0; dx=Rx/real(nx); tx=0.0
  vort=0.0; aa0=0.7; deltat=1.0
  tolerance=2.0e-10
  do i=0,nx
    xgrid(i)=Rx*real(i)/real(nx)
    ! A(r)=(1+0.02*r^2)*I and an azimuthal mean field vy=0.3*r.
    ! Their contractions have analytic answers on every reflected segment.
    factor=1.0+0.02*xgrid(i)**2
    cmp(i)=factor; cpm(i)=factor; coo(i)=factor
    cpp(i)=0.0; cpo(i)=0.0; cop(i)=0.0
    com(i)=0.0; cmo(i)=0.0; cmm(i)=0.0
    vy(i)=0.3*xgrid(i)
  end do
  kx(1)=0.6; ky(1)=0.8; kz(1)=0.6
  kx(2)=1.0; ky(2)=0.0; kz(2)=0.0
  kx(3)=0.0; ky(3)=1.0; kz(3)=1.0
  do direction=1,3
    do case_index=1,size(targets)
      iy=targets(case_index)
      sz=kz(direction); sd=sqrt(1.0-sz*sz)
      sx=kx(direction)*sd; sy=ky(direction)*sd
      origin=xgrid(iy)
      if (iy==nx) origin=Rx-1.0e-7*dx
      call trajectories(origin,0.0,0.0,sx,sy,sz,Rx,Rx*Rx,dx,mx, &
        bx,by,bz,fx,fy,fz,bpx,bpy,bpz,fpx,fpy,fpz)
      call intord_c(direction,direction,iy,sem)
      if (.not. all(ieee_is_finite(real(sem)))) error stop 'Nonfinite real sample'
      if (.not. all(ieee_is_finite(aimag(sem)))) error stop 'Nonfinite imaginary sample'
      do i=0,mx
        factor=1.0+0.02*(fx(i)**2+fy(i)**2)
        expected(1:3)=factor*[fpx(i),fpy(i),fpz(i)]
        expected(4)=0.3*aa0*(-fy(i)*fpx(i)+fx(i)*fpy(i))
        if (maxval(abs(sem(:,i)-expected))>tolerance) error stop 'Forward analytic sample'
        if(i==0) cycle
        factor=1.0+0.02*(bx(i)**2+by(i)**2)
        expected(1:3)=factor*[bpx(i),bpy(i),bpz(i)]
        expected(4)=0.3*aa0*(-by(i)*bpx(i)+bx(i)*bpy(i))
        if (maxval(abs(sem(:,-i)-expected))>tolerance) error stop 'Reverse analytic sample'
      end do
    end do
  end do
  ! Exact wall interpolation must use the last three nodes, not node nx+1.
  r=Rx; ph=0.7
  dummy=intpol(r,ph,ses,1)
  expected(1)=1.0+0.02*r*r
  if(abs(ses(1)-expected(1))>tolerance) error stop 'Wall interpolation'
  if(abs(ses(5)-expected(1))>tolerance) error stop 'Wall interpolation yy'
  if(abs(ses(9)-expected(1))>tolerance) error stop 'Wall interpolation zz'
  print *, 'PASS: cylinder samples match analytic fields at core, interior, and wall'
end program test_radial_cylinder_sampling
