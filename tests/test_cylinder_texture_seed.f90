program test_cylinder_texture_seed
  use he3_kinds, only: rk
  use cylinder_bulk_seed
  implicit none
  complex(rk) :: a(3,3),b(3,3),v(3),cross(3)
  real(rk) :: r,phi,pi,ell(3),beta,normal(3),eps
  integer :: i,j
  pi=acos(-1._rk); eps=1.e-6_rk
  do i=0,20
    r=i/20._rk
    do j=0,31
      phi=2*pi*j/32
      call make_cylinder_texture_seed('a_mermin_ho',r*cos(phi),r*sin(phi),1._rk,1._rk,a)
      v=a(3,:); beta=pi*r/2
      if(maxval(abs(a(:2,:)))>1.e-14_rk) error stop 'wrong spin direction'
      if(abs(sum(v*v))>1.e-13_rk.or.abs(sum(abs(v)**2)-2)>1.e-13_rk) error stop 'not pure A'
      cross=[v(2)*conjg(v(3))-v(3)*conjg(v(2)),v(3)*conjg(v(1))-v(1)*conjg(v(3)), &
        v(1)*conjg(v(2))-v(2)*conjg(v(1))]
      ell=real(cmplx(0._rk,1._rk,rk)*cross/2,rk)
      if(maxval(abs(ell-[sin(beta)*cos(phi),sin(beta)*sin(phi),cos(beta)]))>1.e-13_rk) &
        error stop 'wrong l texture'
      if(i==20) then
        normal=[cos(phi),sin(phi),0._rk]
        if(abs(sum(v*normal))>1.e-13_rk) error stop 'nonzero wall normal pairing'
        call make_cylinder_texture_seed('a_mermin_ho',cos(phi+eps),sin(phi+eps),1._rk,1._rk,b)
        if(abs(aimag(sum(conjg(v)*(b(3,:)-v)))/(2*eps)-1)>1.e-5_rk) error stop 'wall circulation'
      end if
      call make_cylinder_texture_seed('a_panam',r*cos(phi),r*sin(phi),1._rk,1._rk,b)
      if(abs(sum(abs(b)**2)-2)>1.e-13_rk.or.maxval(abs(b(3,:)))>1.e-13_rk) &
        error stop 'PanAm spin norm/plane'
      if(maxval(abs(b(1,:)-cos(-pi*r*r*cos(phi)*sin(phi)/2)*v))>1.e-13_rk) &
        error stop 'PanAm spin angle'
      if(maxval(abs(b(2,:)-sin(-pi*r*r*cos(phi)*sin(phi)/2)*v))>1.e-13_rk) &
        error stop 'PanAm spin angle y'
    end do
  end do
  call make_cylinder_texture_seed('a_mermin_ho',0._rk,0._rk,1._rk,1._rk,a)
  if(abs(a(3,1)-1)>1.e-14_rk.or.abs(a(3,2)-cmplx(0._rk,1._rk,rk))>1.e-14_rk) &
    error stop 'origin limit'
  call make_cylinder_texture_seed('a_planar',0._rk,0.5_rk,1._rk,1._rk,b)
  if(maxval(abs(b(3,:)-[cmplx(0._rk,0._rk,rk),cmplx(0._rk,1/sqrt(2._rk),rk), &
    cmplx(0._rk,1/sqrt(2._rk),rk)]))>1.e-14_rk) error stop 'original polar point changed'
  call make_cylinder_texture_seed('a_planar',0._rk,0._rk,1._rk,1._rk,b)
  if(maxval(abs(a-b))>1.e-14_rk) error stop 'seed central scales differ'
  print *, 'PASS A texture geometry, null condition, origin, wall, circulation and original polar point'
end program
