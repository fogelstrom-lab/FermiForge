program test_cylinder_bulk_seed
  use he3_kinds, only: rk
  use cylinder_bulk_seed
  use historical_core_seed_2d, only: initialize_historical_core_seed_2d
  use cartesian_mesh_2d, only: cartesian_mesh_2d_t, make_uniform_cartesian_mesh
  use spinful_state_2d, only: spinful_state_2d_t, allocate_spinful_state_2d, &
    configure_radial_symmetry,project_radial_origin
  implicit none
  type(cartesian_mesh_2d_t) :: mesh
  type(spinful_state_2d_t) :: state
  complex(rk) :: seed(3,3),actual(3,3)
  real(rk) :: m,current(3),theta
  logical, allocatable :: active(:)
  logical :: inside
  integer :: init,i
  call make_uniform_cartesian_mesh(0._rk,10._rk,10,-10._rk,10._rk,2,mesh)
  allocate(active(mesh%point_count()))
  do init=-2,-1
    call allocate_spinful_state_2d(mesh,state)
    call make_cylinder_bulk_seed(init,0._rk,0.28_rk,seed,m)
    call configure_radial_symmetry(mesh,state,10._rk,m,active)
    do i=1,size(state%radial_point)
      state%order_parameter(:,:,state%radial_point(i))=seed
    end do
    call project_radial_origin(state)
    if(maxval(abs(state%order_parameter(:,:,state%radial_point(1))-seed))>1.e-14_rk) &
      error stop 'bulk seed removed at origin'
    do i=0,32
      theta=2*acos(-1._rk)*i/32
      call state%sample_radial(5*cos(theta),5*sin(theta),actual,current,inside)
      if(.not.inside.or.maxval(abs(actual-seed))>1.e-14_rk) error stop 'unwanted seed phase winding'
    end do
    if(init==-2) then
      if(abs(seed(3,2)-cmplx(0._rk,1._rk,rk)*seed(3,1))>1.e-14_rk) error stop 'wrong chirality'
    end if
  end do
  call initialize_historical_core_seed_2d(mesh,'aop',0.28_rk,0.3_rk,5.4_rk/2.8_rk,state)
  call configure_radial_symmetry(mesh,state,10._rk,1._rk,active)
  call project_radial_origin(state)
  seed=state%order_parameter(:,:,state%radial_point(1))
  if(maxval(abs(seed))<0.01_rk) error stop 'A core removed at origin'
  do i=0,16
    theta=2*acos(-1._rk)*i/16
    call state%sample_radial(0._rk*cos(theta),0._rk*sin(theta),actual,current,inside)
    if(.not.inside.or.maxval(abs(actual-seed))>1.e-14_rk) error stop 'invalid vortex origin'
  end do
  if(maxval(abs(state%order_parameter(:,:,state%radial_point(11))))<0.2_rk) &
    error stop 'missing B background'
  call initialize_historical_core_seed_2d(mesh,'nop',0.28_rk,0.3_rk,5.4_rk/2.8_rk,state)
  call configure_radial_symmetry(mesh,state,10._rk,1._rk,active)
  call project_radial_origin(state)
  if(maxval(abs(state%order_parameter(:,:,state%radial_point(1))))>1.e-14_rk) &
    error stop 'normal core must start with zero gap at origin'
  call state%sample_radial(5._rk,0._rk,seed,current,inside)
  call state%sample_radial(0._rk,5._rk,actual,current,inside)
  if(.not.inside.or.maxval(abs(actual-cmplx(0._rk,1._rk,rk)*seed))>1.e-14_rk) &
    error stop 'normal seed must have unit phase winding'
  print *, 'PASS bulk A/B and normal/A-core cylinder seeds'
end program
