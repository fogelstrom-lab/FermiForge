program test_specular_cylinder_2d
  use he3_kinds, only: rk
  use specular_cylinder_2d
  implicit none
  type(cylinder_geometry_t) :: geometry
  type(cylinder_path_t) :: path
  real(rk) :: p(3),x(2),normal(2),expected(3),ds
  integer :: i,j,collisions
  geometry%enabled=.true.; geometry%radius=2; geometry%half_length=12
  do j=1,4
    select case(j)
    case(1)
      x=[0._rk,0._rk]; p=[1._rk,0._rk,0._rk]
    case(2)
      x=[1.3_rk,0.4_rk]; p=[0.6_rk,0._rk,0.8_rk]
    case(3)
      x=[2._rk,0._rk]; p=[0.8_rk,0.6_rk,0._rk]
    case(4)
      x=[1._rk,0._rk]; p=[0._rk,0._rk,1._rk]
    end select
    call build_cylinder_path(geometry,x,p,0.17_rk,path)
    if(abs(path%s(1)+12)>1.e-12_rk.or.abs(path%s(size(path%s))-12)>1.e-12_rk) &
      error stop 'wrong endpoints'
    if(path%s(path%target)/=0.or.maxval(abs(path%momentum(:,path%target)-p))>1.e-12_rk) &
      error stop 'wrong target'
    collisions=0
    do i=1,size(path%s)
      if(norm2(path%position(:,i))>2+1.e-12_rk) error stop 'escaped wall'
      if(abs(norm2(path%momentum(:,i))-1)>1.e-12_rk) error stop 'lost unit momentum'
      if(abs(path%momentum(3,i)-p(3))>1.e-12_rk) error stop 'changed axial momentum'
      if(i==1) cycle
      ds=path%s(i)-path%s(i-1)
      if(ds<0.or.ds>0.17_rk+1.e-12_rk) error stop 'invalid step'
      if(ds==0) then
        collisions=collisions+1
        normal=path%position(:,i)/2
        expected=path%momentum(:,i-1)
        expected(1:2)=expected(1:2)-2*dot_product(expected(1:2),normal)*normal
        if(maxval(abs(expected-path%momentum(:,i)))>1.e-12_rk) error stop 'reflection law'
      else
        if(maxval(abs(path%position(:,i)-path%position(:,i-1)-ds*path%momentum(1:2,i))) &
          >1.e-11_rk) error stop 'interval crosses collision'
      end if
    end do
    if(j<4.and.collisions<2) error stop 'missing reflections'
    if(j==4.and.collisions/=0) error stop 'axial ray reflected'
  end do
  geometry%collision_aligned=.false.
  call build_cylinder_path(geometry,[0._rk,0._rk],[1._rk,0._rk,0._rk],0.5_rk,path)
  if(size(path%s)/=49.or.any(path%s(2:)<=path%s(:48))) error stop 'uniform compatibility grid'
  print *, 'PASS cylinder geometry, reflection, axial and wall limits'
end program
