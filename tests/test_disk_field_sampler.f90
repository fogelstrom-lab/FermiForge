program test_disk_field_sampler
  use he3_kinds, only: rk
  use cartesian_mesh_2d
  use spinful_state_2d
  use disk_field_sampler_2d
  use, intrinsic :: ieee_arithmetic
  implicit none
  type(cartesian_mesh_2d_t) :: mesh
  type(spinful_state_2d_t) :: state
  type(disk_stencil_t) :: stencil
  real(rk) :: xy(2),r,theta,cur(3),worst,expected(6),actual(6),d(2),w
  complex(rk) :: gap(3,3)
  integer :: p,i,j,k
  logical :: ok
  if(command_argument_count()>0) then
    call make_annular_disk_mesh(1._rk,20,2._rk,0.1_rk,mesh)
  else
  call make_uniform_cartesian_mesh(-1._rk,1._rk,32,-1._rk,1._rk,32,mesh)
  end if
  mesh%cylinder%enabled=.true.; mesh%cylinder%radius=1
  call allocate_spinful_state_2d(mesh,state)
  do p=1,mesh%point_count()
    xy=mesh%point_coordinate(p)
    if(norm2(xy)>1+1.e-14_rk) then
      state%order_parameter(:,:,p)=cmplx(ieee_value(0._rk,ieee_quiet_nan),0._rk,rk)
      state%current_mean_field(:,p)=ieee_value(0._rk,ieee_quiet_nan)
    else
      state%order_parameter(:,:,p)=cmplx(xy(1)**2+xy(2),xy(1)*xy(2),rk)
      state%current_mean_field(:,p)=xy(1)-2*xy(2)**2
    end if
  end do
  worst=0
  do i=0,100
    theta=2*acos(-1._rk)*i/101
    do j=0,10
      r=real(j,rk)/10; xy=r*[cos(theta),sin(theta)]
      call make_disk_stencil(mesh,xy(1),xy(2),stencil)
      if(.not.stencil%valid) error stop 'invalid disk stencil'
      actual=0
      do k=1,stencil%count
        d=mesh%point_coordinate(stencil%point(k)); w=stencil%weight(k)
        if(norm2(d)>1+1.e-14_rk) error stop 'exterior node used'
        actual=actual+w*[1._rk,d(1),d(2),d(1)**2,d(1)*d(2),d(2)**2]
      end do
      expected=[1._rk,xy(1),xy(2),xy(1)**2,xy(1)*xy(2),xy(2)**2]
      worst=max(worst,maxval(abs(actual-expected)))
      call sample_disk_state(mesh,state,xy(1),xy(2),gap,cur,ok)
      if(.not.ok.or..not.all(ieee_is_finite(real(gap))).or..not.all(ieee_is_finite(cur))) &
        error stop 'exterior poison reached sample'
      if(maxval(abs(gap-cmplx(xy(1)**2+xy(2),xy(1)*xy(2),rk)))>1.e-10_rk) &
        error stop 'complex polynomial not reproduced'
    end do
  end do
  if(worst>1.e-10_rk) error stop 'quadratic reproduction failed'
  do i=1,12
    xy=[1._rk-10._rk**(-i),0._rk]
    call sample_disk_state(mesh,state,xy(1),xy(2),gap,cur,ok)
    if(.not.ok) error stop 'near-node stencil lost rank'
    if(maxval(abs(gap-cmplx(xy(1)**2,0._rk,rk)))>1.e-9_rk) error stop 'near-node error'
  end do
  call make_disk_stencil(mesh,1.001_rk,0._rk,stencil)
  if(stencil%valid) error stop 'outside query accepted'
  print *, 'PASS interior-only disk quadratics, maximum error:',worst
end program
