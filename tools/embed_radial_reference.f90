! Solver-free radial embedding onto an existing field-map grid.
program embed_radial_reference
  use he3_kinds, only: rk
  use radial_mesh, only: radial_mesh_t
  use cartesian_mesh_2d, only: cartesian_mesh_2d_t, make_rectilinear_cartesian_mesh
  use axial_radial_embedding, only: axial_radial_profile_t, embed_axial_radial_profile_2d
  use legacy_radial_profile_io, only: read_new_src_radial_profile
  use spinful_state_2d, only: spinful_state_2d_t, allocate_spinful_state_2d
  use spinful_field_io_2d, only: read_spinful_field_map_2d, write_spinful_field_map_2d
  implicit none
  type(radial_mesh_t) :: radial
  type(axial_radial_profile_t) :: profile
  type(cartesian_mesh_2d_t) :: mesh
  type(spinful_state_2d_t) :: state
  character(len=2048) :: op, curr, template, output, line
  integer :: u, ios, nx, ny, n, ix, iy
  real(rk) :: row(24)
  real(rk), allocatable :: x(:), y(:)
  if (command_argument_count()/=4) error stop 'usage: embed_radial_reference op curr template output (unit winding)'
  call get_command_argument(1,op)
  call get_command_argument(2,curr)
  call get_command_argument(3,template)
  call get_command_argument(4,output)
  nx=0; ny=0
  open(newunit=u,file=trim(template),status='old',action='read')
  do
    read(u,'(a)',iostat=ios) line
    if (ios/=0) exit
    if (index(line,'# nx_points=')==1) read(line(13:),*) nx
    if (index(line,'# ny_points=')==1) read(line(13:),*) ny
  end do
  if (nx<2.or.ny<2) error stop 'missing grid dimensions'
  allocate(x(nx),y(ny))
  rewind(u)
  n=0
  do
    read(u,'(a)',iostat=ios) line
    if (ios/=0) exit
    if (len_trim(line)==0.or.line(1:1)=='#') cycle
    n=n+1
    if (n>nx*ny) error stop 'too many grid rows'
    read(line,*) row
    ix=modulo(n-1,nx)+1; iy=(n-1)/nx+1
    if (iy==1) x(ix)=row(1)
    if (ix==1) y(iy)=row(2)
  end do
  close(u)
  if(n/=nx*ny) error stop 'incomplete grid'
  call make_rectilinear_cartesian_mesh(x,y,mesh)
  call allocate_spinful_state_2d(mesh,state)
  ! Validate every template coordinate using the production reader.
  call read_spinful_field_map_2d(trim(template),mesh,state)
  call read_new_src_radial_profile(trim(op),trim(curr),radial,profile)
  call embed_axial_radial_profile_2d(radial,profile,1.0_rk,mesh,state)
  call write_spinful_field_map_2d(trim(output),mesh,state,'updated_radial_reference','reference_only')
end program embed_radial_reference
