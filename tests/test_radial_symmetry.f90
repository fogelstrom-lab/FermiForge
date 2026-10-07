program test_radial_symmetry
  use he3_kinds, only: rk
  use cartesian_mesh_2d, only: cartesian_mesh_2d_t, make_uniform_cartesian_mesh
  use spinful_state_2d, only: spinful_state_2d_t, allocate_spinful_state_2d, &
    configure_radial_symmetry, project_radial_origin
  use order_parameter_basis, only: axial_harmonics_to_cartesian, cartesian_to_axial_harmonics, &
    reconstruct_axial_cartesian
  use radial_mesh, only: radial_mesh_t, make_uniform_radial_mesh
  use axial_radial_embedding, only: axial_radial_profile_t, allocate_axial_radial_profile, &
    sample_axial_radial_profile
  use bilinear_field_sampler_2d, only: sample_spinful_state_2d
  implicit none
  type(cartesian_mesh_2d_t) :: mesh
  type(spinful_state_2d_t) :: state
  type(radial_mesh_t) :: radial
  type(axial_radial_profile_t) :: profile
  logical, allocatable :: active(:)
  complex(rk) :: h(3,3), expected(3,3), actual(3,3), reference(3,3)
  real(rk) :: r,theta,v(3),vr(3)
  logical :: inside
  integer :: i,p,j
  call make_uniform_cartesian_mesh(-4._rk,4._rk,16,-4._rk,4._rk,16,mesh)
  mesh%trajectory_interpolation_order=2
  call allocate_spinful_state_2d(mesh,state)
  allocate(active(mesh%point_count()))
  call configure_radial_symmetry(mesh,state,4._rk,1._rk,active)
  if(count(active)/=9) error stop 'ray count'
  call make_uniform_radial_mesh(4._rk,8,radial)
  call allocate_axial_radial_profile(radial,profile)
  do i=1,size(state%radial_point)
    r=state%radial_coordinate(i); p=state%radial_point(i)
    h=0
    h(1,2)=1+0.1_rk*r*r
    h(2,2)=cmplx(r,0.2_rk*r,rk)
    h(3,2)=0.3_rk*r*r
    call axial_harmonics_to_cartesian(h,state%order_parameter(:,:,p))
    state%current_mean_field(:,p)=[0._rk,0.2_rk*r,0._rk]
    profile%harmonic(:,:,i-1)=h
    profile%azimuthal_mean_field(i-1)=0.2_rk*r
  end do
  ! All off-ray data is zero: sampling must not use the Cartesian cache.
  do j=1,17
    r=0.17_rk*real(j,rk); theta=0.37_rk*real(j,rk)
    h=0; h(1,2)=1+0.1_rk*r*r; h(2,2)=cmplx(r,0.2_rk*r,rk); h(3,2)=0.3_rk*r*r
    call reconstruct_axial_cartesian(h,theta,1._rk,expected)
    call sample_spinful_state_2d(mesh,state,r*cos(theta),r*sin(theta),actual,v,inside)
    if(.not.inside.or.maxval(abs(actual-expected))>1.e-12_rk) error stop 'analytic radial rotation'
    call sample_axial_radial_profile(radial,profile,r*cos(theta),r*sin(theta),1._rk,reference,vr,inside)
    if(maxval(abs(actual-reference))>1.e-12_rk.or.maxval(abs(v-vr))>1.e-12_rk) &
      error stop 'legacy radial embedding agreement'
  end do
  call state%sample_radial(5._rk,0._rk,actual,v,inside)
  if(inside) error stop 'outside ray support'
  p=state%radial_point(1)
  state%order_parameter(:,:,p)=1
  state%current_mean_field(:,p)=1
  call project_radial_origin(state)
  call cartesian_to_axial_harmonics(state%order_parameter(:,:,p),h)
  h(1,2)=0; h(2,1)=0
  if(maxval(abs(h))>1.e-12_rk.or.any(state%current_mean_field(1:2,p)/=0)) &
    error stop 'origin regularity'
  print *, 'radial symmetry sampling tests passed'
end program test_radial_symmetry
