! Reduced-quadrature execution smoke, not a converged physical benchmark.
program test_radial_cylinder_map
  use Global_variables
  use MPI_variables
  use Initialisation
  use mpicalls
  use NewSES
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  integer :: ierr, i
  complex :: fields(12,0:nx)
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  if(myid==0) then
    call init_calc
    ! Keep this test bounded: two azimuths, eleven polar nodes, one pole.
    Ncmax=1
  end if
  call input_bcast
  ! Full-precision initial fields for an independent modern-stack map check.
  if(myid==0) then
    open(81,file='initial_op_xyz',status='replace')
    open(82,file='initial_curr',status='replace')
    do i=0,nx
      write(81,'(19es26.17e3)') xgrid(i),dxx(i),dxy(i),dxz(i),dyx(i),dyy(i),dyz(i),dzx(i),dzy(i),dzz(i)
      write(82,'(5es26.17e3)') xgrid(i),0.0,real(vx(i)),real(vy(i)),real(vz(i))
    end do
    close(81); close(82)
  end if
  if (.not. cyl) error stop 'Cylinder switch was not enabled'
  if (abs(xgrid(nx)-Rx)>1.0e-12) error stop 'Grid does not end at cylinder wall'
  do i=0,nx
    if(abs(xgrid(i)-Rx*real(i)/real(nx))>1.0e-12) error stop 'Cylinder grid mismatch'
  end do
  call getnewop
  if(myid==0) then
    fields(1,:)=dxx; fields(2,:)=dxy; fields(3,:)=dxz
    fields(4,:)=dyx; fields(5,:)=dyy; fields(6,:)=dyz
    fields(7,:)=dzx; fields(8,:)=dzy; fields(9,:)=dzz
    fields(10,:)=vx; fields(11,:)=vy; fields(12,:)=vz
    if(.not. all(ieee_is_finite(real(fields)))) error stop 'Nonfinite cylinder map'
    if(.not. all(ieee_is_finite(aimag(fields)))) error stop 'Nonfinite cylinder map'
    if(maxval(abs(fields(1:9,:)))==0.0) error stop 'Cylinder map is identically zero'
    open(80,file='cylinder_map.dat',status='replace')
    do i=0,nx
      write(80,'(25(es25.16,1x))') xgrid(i),fields(:,i)
    end do
    close(80)
    print *, 'PASS: reduced-quadrature cylinder map is finite on every radial node'
  end if
  call MPI_FINALIZE(ierr)
end program test_radial_cylinder_map
