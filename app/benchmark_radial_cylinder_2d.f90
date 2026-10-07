! Radially constrained entry point to the modern point-map/MPI/AA stack.
program benchmark_radial_cylinder_2d
  use mpi_f08
  use he3_kinds, only: rk
  use cylinder_bulk_seed, only: make_cylinder_bulk_seed
  use historical_core_seed_2d, only: initialize_historical_core_seed_2d
  use cartesian_mesh_2d, only: cartesian_mesh_2d_t, make_uniform_cartesian_mesh
  use spinful_state_2d, only: spinful_state_2d_t, allocate_spinful_state_2d, &
    configure_radial_symmetry, project_radial_origin
  use he3_quadrature, only: angular_quadrature_3d_t, ozaki_quadrature_t, &
    read_legacy_gauss_table, read_legacy_ozaki_table
  use he3_bulk_gap, only: bulk_gap_legacy
  use he3_mpi_field_map_2d, only: evaluate_mpi_he3_field_map, &
    update_mpi_he3_state_with_anderson, he3_mpi_field_map_diagnostics_t
  use he3_nonlinear_iteration_2d, only: compute_he3_field_residual, he3_field_residual_report_t
  use legacy_anderson_mixing, only: legacy_anderson_t, anderson_report_t
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(cartesian_mesh_2d_t) :: mesh
  type(spinful_state_2d_t) :: state,mapped,next
  type(angular_quadrature_3d_t) :: angular
  type(ozaki_quadrature_t) :: energy
  type(he3_mpi_field_map_diagnostics_t) :: diagnostics
  type(he3_field_residual_report_t) :: residual
  type(legacy_anderson_t) :: accelerator
  type(anderson_report_t) :: report
  real(rk) :: radius=10,half_length=8000._rk/99,maximum_step=10._rk/99
  real(rk) :: symmetry_winding
  complex(rk) :: seed(3,3)
  integer :: initialization=-1
  real(rk) :: temperature=0.30_rk,fs1=0,winding=0,tolerance=2.e-6_rk,p_max=3
  real(rk) :: boundary_relaxation_distance=100._rk/99,feedback,gap,r,values(19),curr(5)
  integer :: radial_points=100,azimuths=32,max_iterations=0,pole_limit=0
  logical :: collision_aligned=.true.,success
  character(len=1024) :: gauss_file='new_src/gauss11.dat',ozaki_file='new_src/ozaki_T=0.3.dat'
  character(len=1024) :: initial_gap_file='',initial_current_file='',input,output_prefix='cylinder'
  character(len=32) :: status
  integer :: ierr,rank,unit,cu,history,i,j,k,p,iteration,ios
  logical, allocatable :: active(:)
  namelist /radial_cylinder/ radius,half_length,maximum_step,temperature,fs1,winding, &
    tolerance,p_max,boundary_relaxation_distance,radial_points,azimuths,max_iterations, &
    pole_limit,collision_aligned,gauss_file,ozaki_file,initial_gap_file,initial_current_file,output_prefix,initialization
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  if(command_argument_count()/=1) error stop 'usage: benchmark_radial_cylinder_2d input.nml'
  call get_command_argument(1,input)
  open(newunit=unit,file=trim(input),status='old',action='read')
  read(unit,nml=radial_cylinder); close(unit)
  if(.not.all(ieee_is_finite([radius,half_length,maximum_step,temperature,fs1,winding, &
    tolerance,p_max,boundary_relaxation_distance]))) error stop 'nonfinite cylinder input'
  if(radius<=0.or.half_length<=0.or.maximum_step<=0.or.radial_points<3.or.azimuths<1.or. &
    fs1<0.or.tolerance<=0.or.max_iterations<0.or.pole_limit<0) error stop 'invalid cylinder input'
  if(initialization==0.or.initialization==1) then
    if(winding/=1) error stop 'vortex cylinder seeds require unit winding'
  else
    if(winding/=0) error stop 'bulk cylinder driver requires zero physical winding'
  end if
  if((len_trim(initial_gap_file)==0).neqv.(len_trim(initial_current_file)==0)) &
    error stop 'provide both initial radial files'
  call read_legacy_gauss_table(trim(gauss_file),azimuths,11,angular)
  call read_legacy_ozaki_table(trim(ozaki_file),energy)
  if(abs(energy%temperature-temperature)>1.e-10_rk) error stop 'Ozaki temperature mismatch'
  ! A pole limit is solely for reduced-work tests; keep the full-table prefactor.
  if(pole_limit>energy%pole_count()) error stop 'pole limit exceeds table'
  if(pole_limit>0) then
    energy%pole=energy%pole(:pole_limit); energy%residue=energy%residue(:pole_limit)
  end if
  feedback=fs1/(1+fs1/3)
  call make_uniform_cartesian_mesh(0._rk,radius,radial_points-1,-radius,radius,2,mesh)
  mesh%trajectory_interpolation_order=2
  mesh%cylinder%enabled=.true.; mesh%cylinder%radius=radius
  mesh%cylinder%half_length=half_length; mesh%cylinder%collision_aligned=collision_aligned
  call allocate_spinful_state_2d(mesh,state)
  allocate(active(mesh%point_count()))
  gap=bulk_gap_legacy(temperature)
  if(initialization==0.or.initialization==1) then
    if(initialization==0) then
      call initialize_historical_core_seed_2d(mesh,'nop',gap,temperature,feedback,state)
    else
      call initialize_historical_core_seed_2d(mesh,'aop',gap,temperature,feedback,state)
    end if
    symmetry_winding=1
  else
    call make_cylinder_bulk_seed(initialization,winding,gap,seed,symmetry_winding)
  end if
  call configure_radial_symmetry(mesh,state,radius,symmetry_winding,active)
  do i=1,radial_points
    p=state%radial_point(i)
    if(initialization<0) state%order_parameter(:,:,p)=seed
  end do
  if(len_trim(initial_gap_file)>0) then
    open(newunit=unit,file=trim(initial_gap_file),status='old',action='read')
    open(newunit=cu,file=trim(initial_current_file),status='old',action='read')
    do i=1,radial_points
      read(unit,*,iostat=ios) values
      if(ios/=0) error stop 'invalid initial gap profile'
      read(cu,*,iostat=ios) curr
      if(ios/=0) error stop 'invalid initial current profile'
      r=state%radial_coordinate(i); p=state%radial_point(i)
      if(abs(values(1)-r)>5.1e-4_rk.or.abs(curr(1)-r)>5.1e-4_rk) &
        error stop 'initial profile grid does not match cylinder'
      if(.not.all(ieee_is_finite(values)).or..not.all(ieee_is_finite(curr))) &
        error stop 'nonfinite initial field'
      do j=1,3
        do k=1,3
          state%order_parameter(j,k,p)=cmplx(values(2*((j-1)*3+k)),values(2*((j-1)*3+k)+1),rk)
        end do
      end do
      state%current_mean_field(:,p)=curr(3:5)
    end do
    read(unit,*,iostat=ios) values
    if(ios==0) error stop 'extra initial gap rows'
    read(cu,*,iostat=ios) curr
    if(ios==0) error stop 'extra initial current rows'
    close(unit); close(cu)
  end if
  call project_radial_origin(state)
  call accelerator%initialize(21*count(active),10,maximum_mixing=p_max)
  if(rank==0) then
    print *, 'initialization (-1 bulk B, -2 bulk A, 0 normal-core, 1 A-core): ',initialization
    print *, 'physical phase winding: ',winding,' radial symmetry label: ',symmetry_winding
    call write_profile(state,'initial')
    open(newunit=history,file=trim(output_prefix)//'_history.dat',status='replace')
    write(history,'(a)') '# evaluated_state rms max relative_l2 normalization map_seconds'
  end if
  status='iteration_limit'
  ! max_iterations counts updates, not maps. Evaluate the final saved state too.
  do iteration=0,max_iterations
    call evaluate_mpi_he3_field_map(MPI_COMM_WORLD,mesh,state,angular,energy,feedback, &
      maximum_step,1,boundary_relaxation_distance,mapped,diagnostics,success,active_point_mask=active)
    if(.not.success) call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    call compute_he3_field_residual(state,mapped,residual,active)
    if(rank==0) then
      write(history,'(i6,5es22.12)') iteration,residual%rms_residual, &
        residual%maximum_absolute_residual,residual%relative_l2_residual, &
        diagnostics%global%maximum_normalization_error,diagnostics%maximum_rank_elapsed_seconds
      flush(history)
      print '(a,i5,a,2es12.4,a,f10.2)', 'cylinder state ',iteration,' residual(rms,max)=', &
        residual%rms_residual,residual%maximum_absolute_residual,' map[s]=', &
        diagnostics%maximum_rank_elapsed_seconds
      flush(6)
      call write_profile(state,'final')
      call write_profile(mapped,'mapped')
    end if
    if(max_iterations==0) then
      status='single_map'; exit
    end if
    if(residual%maximum_absolute_residual<=tolerance) then
      status='converged'; exit
    end if
    if(iteration==max_iterations) exit
    call update_mpi_he3_state_with_anderson(MPI_COMM_WORLD,mesh,accelerator,state,mapped, &
      tolerance,next,report,residual,active)
    state=next
    call configure_radial_symmetry(mesh,state,radius,symmetry_winding,active)
    call project_radial_origin(state)
  end do
  if(rank==0) then
    close(history)
    print '(a,a)', 'terminal status: ',trim(status)
  end if
  call MPI_Finalize(ierr)
contains
  subroutine write_profile(field,label)
    type(spinful_state_2d_t), intent(in) :: field
    character(len=*), intent(in) :: label
    integer :: u,v,idx,point,spin,orb
    open(newunit=u,file=trim(output_prefix)//'_'//label//'_op_xyz',status='replace')
    open(newunit=v,file=trim(output_prefix)//'_'//label//'_curr',status='replace')
    do idx=1,size(state%radial_point)
      point=state%radial_point(idx)
      write(u,'(19es26.17e3)') state%radial_coordinate(idx), &
        ((field%order_parameter(spin,orb,point),orb=1,3),spin=1,3)
      write(v,'(5es26.17e3)') state%radial_coordinate(idx), &
        sqrt(sum(abs(field%order_parameter(:,:,point))**2)/3),field%current_mean_field(:,point)
    end do
    close(u); close(v)
  end subroutine
end program
