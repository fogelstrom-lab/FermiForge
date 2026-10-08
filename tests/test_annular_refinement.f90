program test_annular_refinement
  use he3_kinds, only: rk
  use cartesian_mesh_2d
  use spinful_state_2d
  use spinful_field_io_2d
  use disk_field_sampler_2d
  implicit none
  type(cartesian_mesh_2d_t) :: mesh
  type(spinful_state_2d_t) :: state,reloaded
  complex(rk) :: gap(3,3),left(3,3),right(3,3)
  real(rk) :: xy(2),theta,r,current(3),errors(3),expected
  integer :: level,p,i,j
  logical :: ok
  character(len=16) :: layout='wall'
  if(command_argument_count()>0) layout='core_wall'
  do level=1,3
    call make_annular_disk_mesh(1._rk,8*2**(level-1),2._rk,0.25_rk/2**(level-1),mesh,layout)
    if(layout=='core_wall') then
      if(abs(mesh%ring_radius(1)-(1-mesh%ring_radius(ubound(mesh%ring_radius,1)-1)))>1.e-14_rk) &
        error stop 'core/wall spacing not symmetric'
      if(minval(mesh%ring_count(1:))<24) error stop 'inner angular coverage'
    end if
    if(.not.mesh%is_valid().or.mesh%is_uniform()) error stop 'annular mesh identity'
    if(mesh%point_count()/=sum(mesh%ring_count)) error stop 'annular count'
    if(mesh%ring_radius(ubound(mesh%ring_radius,1))/=1) error stop 'wall missing'
    call allocate_spinful_state_2d(mesh,state)
    do p=1,mesh%point_count()
      xy=mesh%point_coordinate(p)
      state%order_parameter(:,:,p)=profile(xy)
      state%current_mean_field(:,p)=profile(xy)
    end do
    call write_spinful_field_map_2d('annular_roundtrip.dat',mesh,state,'test','disk')
    call allocate_spinful_state_2d(mesh,reloaded)
    call read_spinful_field_map_2d('annular_roundtrip.dat',mesh,reloaded)
    if(maxval(abs(state%order_parameter-reloaded%order_parameter))>1.e-15_rk) error stop 'restart gap'
    if(maxval(abs(state%current_mean_field-reloaded%current_mean_field))>1.e-15_rk) error stop 'restart field'
    errors(level)=0
    do i=0,60
      theta=2*acos(-1._rk)*(i+0.31_rk)/61
      do j=0,40
        r=real(j,rk)/40; xy=r*[cos(theta),sin(theta)]
        call sample_disk_state(mesh,state,xy(1),xy(2),gap,current,ok)
        if(.not.ok) error stop 'annular refinement stencil'
        expected=profile(xy)
        errors(level)=max(errors(level),maxval(abs(gap-expected)))
      end do
    end do
    ! Cross radial and angular stencil switches away from exact nodes.
    do j=1,ubound(mesh%ring_radius,1)-1
      r=mesh%ring_radius(j); theta=0.137_rk
      xy=(r-1.e-9_rk)*[cos(theta),sin(theta)]
      call sample_disk_state(mesh,state,xy(1),xy(2),left,current,ok)
      if(.not.ok) error stop 'left ring stencil'
      xy=(r+1.e-9_rk)*[cos(theta),sin(theta)]
      call sample_disk_state(mesh,state,xy(1),xy(2),right,current,ok)
      if(.not.ok.or.maxval(abs(left-right))>1.e-6_rk) error stop 'ring stencil jump'
      r=(mesh%ring_radius(j)+mesh%ring_radius(j+1))/2
      theta=2*acos(-1._rk)/mesh%ring_count(j)
      xy=r*[cos(theta-1.e-9_rk),sin(theta-1.e-9_rk)]
      call sample_disk_state(mesh,state,xy(1),xy(2),left,current,ok)
      if(.not.ok) error stop 'left angular stencil'
      xy=r*[cos(theta+1.e-9_rk),sin(theta+1.e-9_rk)]
      call sample_disk_state(mesh,state,xy(1),xy(2),right,current,ok)
      if(.not.ok.or.maxval(abs(left-right))>1.e-6_rk) error stop 'angular stencil jump'
    end do
  end do
  print *, 'annular nonpolynomial max errors:',errors
  if(any(errors(2:)>0.7_rk*errors(:2))) error stop 'annular refinement did not reduce error'
contains
  pure real(rk) function profile(xy) result(value)
    real(rk),intent(in)::xy(2)
    value=exp(-(1-sum(xy**2))/0.2_rk)*(1+0.1_rk*xy(1))
    if(layout=='core_wall') value=value+0.2_rk*exp(-sum(xy**2)/0.02_rk)
  end function
end program
