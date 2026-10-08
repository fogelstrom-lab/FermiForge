program test_disk_point_map
  use he3_kinds, only: rk
  use cartesian_mesh_2d
  use spinful_state_2d
  use he3_quadrature
  use he3_serial_point_map
  implicit none
  type(cartesian_mesh_2d_t) :: disk,ray
  type(spinful_state_2d_t) :: planar,radial
  type(angular_quadrature_3d_t) :: angular
  type(ozaki_quadrature_t) :: energy
  type(he3_point_map_diagnostics_t) :: diag
  complex(rk) :: a(3,3),b(3,3)
  real(rk) :: u(3),v(3),xy(2),r2,err
  integer :: p,j,k,mode
  logical, allocatable :: active(:)
  logical :: ok
  call make_uniform_cartesian_mesh(-4._rk,4._rk,16,-4._rk,4._rk,16,disk)
  if(command_argument_count()>0) call make_annular_disk_mesh(4._rk,12,2._rk,0.5_rk,disk)
  disk%cylinder%enabled=.true.; disk%cylinder%radius=4; disk%cylinder%half_length=20
  call make_uniform_cartesian_mesh(0._rk,4._rk,40,-4._rk,4._rk,2,ray)
  ray%cylinder=disk%cylinder; ray%trajectory_interpolation_order=2
  call allocate_spinful_state_2d(disk,planar)
  call allocate_spinful_state_2d(ray,radial)
  allocate(active(ray%point_count()))
  call configure_radial_symmetry(ray,radial,4._rk,0._rk,active)
  do p=1,disk%point_count()
    xy=disk%point_coordinate(p); r2=sum(xy**2)
    do j=1,3
      planar%order_parameter(j,j,p)=0.2_rk+0.001_rk*r2
    end do
    planar%current_mean_field(:,p)=[-0.001_rk*xy(2),0.001_rk*xy(1),0._rk]
  end do
  do p=1,size(radial%radial_point)
    k=radial%radial_point(p); r2=radial%radial_coordinate(p)**2
    do j=1,3
      radial%order_parameter(j,j,k)=0.2_rk+0.001_rk*r2
    end do
    radial%current_mean_field(:,k)=[0._rk,0.001_rk*radial%radial_coordinate(p),0._rk]
  end do
  call make_legacy_angular_quadrature(4,[-0.35_rk,0.35_rk],[1._rk,1._rk],angular)
  call make_ozaki_quadrature(0.3_rk,12,[0.55_rk],[0.8_rk],energy)
  err=0
  do mode=0,1
    disk%cylinder%collision_aligned=mode==1; ray%cylinder=disk%cylinder
    do p=1,3
      select case(p)
      case(1); xy=[0._rk,0._rk]
      case(2); xy=[0.7_rk,-1.1_rk]
      case(3); xy=[4._rk,0._rk]
      end select
      call evaluate_serial_he3_point_map(disk,planar,xy,angular,energy,1.9_rk, &
        0.1_rk,1,1._rk,a,u,diag,ok)
      if(.not.ok) error stop '2D disk map failed'
      call evaluate_serial_he3_point_map(ray,radial,xy,angular,energy,1.9_rk, &
        0.1_rk,1,1._rk,b,v,diag,ok)
      if(.not.ok) error stop 'radial map failed'
      err=max(err,maxval(abs(a-b)),maxval(abs(u-v)))
    end do
  end do
  if(err>1.e-10_rk) error stop 'radial and 2D polynomial map disagree'
  print *, 'PASS radial/2D reflected polynomial map, error:',err
end program
