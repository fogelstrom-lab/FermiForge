program fermiforge_2d_mpi_smoke
  use mpi_f08
  use he3_kinds, only : rk
  use radial_mesh, only : radial_mesh_t
  use axial_radial_embedding, only : axial_radial_profile_t, &
                                     embed_axial_radial_profile_2d
  use legacy_radial_profile_io, only : read_new_src_radial_profile
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, &
                                make_uniform_cartesian_mesh
  use spinful_state_2d, only : spinful_state_2d_t, allocate_spinful_state_2d
  use spinful_field_io_2d, only : write_spinful_field_map_2d
  use he3_quadrature, only : angular_quadrature_3d_t, ozaki_quadrature_t, &
                             read_legacy_gauss_table, &
                             read_legacy_ozaki_table
  use he3_mpi_field_map_2d, only : he3_mpi_field_map_diagnostics_t, &
                                   evaluate_mpi_he3_field_map, &
                                   update_mpi_he3_state_with_anderson
  use new_src_iteration_layout_2d, only : new_src_iteration_vector_size_2d
  use legacy_anderson_mixing, only : legacy_anderson_t, anderson_report_t
  use he3_nonlinear_iteration_2d, only : he3_field_residual_report_t
  use free_vortex_asymptotic_2d, only : free_vortex_endpoint_2d_t
  implicit none

  type(angular_quadrature_3d_t) :: angular
  type(anderson_report_t) :: anderson_report
  type(cartesian_mesh_2d_t) :: mesh
  type(he3_mpi_field_map_diagnostics_t) :: map_diagnostics
  type(he3_field_residual_report_t) :: residual_report
  type(legacy_anderson_t) :: accelerator
  type(ozaki_quadrature_t) :: energy
  type(axial_radial_profile_t) :: profile
  type(radial_mesh_t) :: radial_grid
  type(spinful_state_2d_t) :: mapped, next, state
  type(free_vortex_endpoint_2d_t) :: free_vortex_endpoint

  character(len=32) :: endpoint_policy
  character(len=512) :: current_field_file, gauss_file, input_file
  character(len=512) :: input_message, order_parameter_file, output_file
  character(len=512) :: output_override, ozaki_file
  integer :: anderson_history_limit, ierr, input_status, input_unit
  integer :: iteration, iteration_count, minimum_substep_count
  integer :: number_of_cells, rank, rank_count
  real(rk) :: anderson_maximum_mixing, anderson_progress_threshold
  real(rk) :: asymptotic_bulk_gap, asymptotic_center_x, asymptotic_center_y
  real(rk) :: asymptotic_maximum_step, asymptotic_outer_radius
  real(rk) :: asymptotic_winding
  real(rk) :: boundary_relaxation_distance, convergence_tolerance
  real(rk) :: feedback_scale, half_width, trajectory_maximum_step
  logical :: map_succeeded

  namelist /field_map_2d/ half_width, number_of_cells, feedback_scale, &
    trajectory_maximum_step, minimum_substep_count, &
    boundary_relaxation_distance, iteration_count, convergence_tolerance, &
    anderson_history_limit, anderson_progress_threshold, &
    anderson_maximum_mixing, order_parameter_file, current_field_file, &
    gauss_file, ozaki_file, output_file, endpoint_policy, &
    asymptotic_center_x, asymptotic_center_y, asymptotic_bulk_gap, &
    asymptotic_winding, asymptotic_outer_radius, asymptotic_maximum_step

  half_width = 8.0_rk
  number_of_cells = 16
  feedback_scale = 1.9285714285714286_rk
  trajectory_maximum_step = 0.25_rk
  minimum_substep_count = 1
  boundary_relaxation_distance = 1.75_rk
  iteration_count = 1
  convergence_tolerance = 1.0e-6_rk
  anderson_history_limit = 10
  anderson_progress_threshold = 0.1_rk
  anderson_maximum_mixing = 5.0_rk
  endpoint_policy = "local"
  asymptotic_center_x = 0.0_rk
  asymptotic_center_y = 0.0_rk
  asymptotic_bulk_gap = 0.2799_rk
  asymptotic_winding = 1.0_rk
  asymptotic_outer_radius = 24.0_rk
  asymptotic_maximum_step = 0.25_rk
  order_parameter_file = "new_src/op_xyz"
  current_field_file = "new_src/curr"
  gauss_file = "new_src/gauss11.dat"
  ozaki_file = "new_src/ozaki.dat"
  output_file = "work/fermiforge_2d_mpi_smoke.dat"

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI initialization failed"
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI rank query failed"
  call MPI_Comm_size(MPI_COMM_WORLD, rank_count, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI size query failed"

  call get_command_argument(1, input_file)
  if (len_trim(input_file) == 0) input_file = "examples/2d_mpi_smoke.nml"
  open(newunit=input_unit, file=trim(input_file), status="old", action="read", &
       iostat=input_status, iomsg=input_message)
  if (input_status /= 0) &
    error stop "cannot open MPI field-map input: " // trim(input_message)
  read(input_unit, nml=field_map_2d, iostat=input_status, iomsg=input_message)
  close(input_unit)
  if (input_status /= 0) &
    error stop "cannot read MPI field-map input: " // trim(input_message)
  call get_command_argument(2, output_override)
  if (len_trim(output_override) > 0) output_file = output_override
  call configure_endpoint()
  call validate_input()

  call read_new_src_radial_profile(trim(order_parameter_file), &
    trim(current_field_file), radial_grid, profile)
  if (sqrt(2.0_rk) * half_width > radial_grid%outer_radius()) &
    error stop "Cartesian smoke mesh extends beyond the radial profile"
  call read_legacy_gauss_table(trim(gauss_file), 48, 11, angular)
  call read_legacy_ozaki_table(trim(ozaki_file), energy)
  call make_uniform_cartesian_mesh( &
    -half_width, half_width, number_of_cells, &
    -half_width, half_width, number_of_cells, mesh)
  call allocate_spinful_state_2d(mesh, state)
  call embed_axial_radial_profile_2d( &
    radial_grid, profile, 1.0_rk, mesh, state)
  if (rank == 0) call accelerator%initialize( &
    new_src_iteration_vector_size_2d(state), anderson_history_limit, &
    anderson_progress_threshold, anderson_maximum_mixing)

  do iteration = 1, iteration_count
    call evaluate_mpi_he3_field_map( &
      MPI_COMM_WORLD, mesh, state, angular, energy, feedback_scale, &
      trajectory_maximum_step, minimum_substep_count, &
      boundary_relaxation_distance, mapped, map_diagnostics, map_succeeded, &
      free_vortex_endpoint)
    if (.not. map_succeeded) error stop "MPI field map failed"
    call update_mpi_he3_state_with_anderson( &
      MPI_COMM_WORLD, mesh, accelerator, state, mapped, &
      convergence_tolerance, next, anderson_report, residual_report)
    if (rank == 0) then
      write(*, '(a,i4,a,2(es12.4,1x),a,f7.3,a,i3,a,f8.3)') &
        "iteration ", iteration, " residual(rms,max)= ", &
        residual_report%rms_residual, &
        residual_report%maximum_absolute_residual, " p=", &
        anderson_report%mixing, " md=", anderson_report%history_size, &
        " map-time=", map_diagnostics%maximum_rank_elapsed_seconds
      write(*, '(a,i0,a,es12.4)') "  worst point ", &
        residual_report%maximum_point, " point-rms=", &
        residual_report%maximum_point_rms_residual
    end if
    state = next
    if (anderson_report%converged) exit
  end do

  if (rank == 0) then
    call write_spinful_field_map_2d( &
      trim(output_file), mesh, state, &
      "mpi_self_consistency_iteration_smoke", trim(endpoint_policy))
    write(*, '(a,i0)') "MPI ranks: ", rank_count
    write(*, '(a,a)') "output: ", trim(output_file)
    write(*, '(a,a)') "trajectory endpoint: ", trim(endpoint_policy)
    write(*, '(a)') "status: iteration smoke test; not a converged solution"
  end if

  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI finalization failed"

contains

  subroutine validate_input()
    if (half_width <= 0.0_rk .or. number_of_cells < 2) &
      error stop "invalid MPI field-map mesh"
    if (trajectory_maximum_step <= 0.0_rk .or. &
        minimum_substep_count < 1 .or. boundary_relaxation_distance < 0.0_rk) &
      error stop "invalid MPI field-map trajectory controls"
    if (iteration_count < 1 .or. anderson_history_limit < 1 .or. &
        convergence_tolerance < 0.0_rk) &
      error stop "invalid MPI field-map iteration controls"
    if (free_vortex_endpoint%enabled .and. &
        (asymptotic_bulk_gap <= 0.0_rk .or. &
         asymptotic_outer_radius <= 0.0_rk .or. &
         asymptotic_maximum_step <= 0.0_rk)) &
      error stop "invalid free-vortex asymptotic controls"
  end subroutine validate_input


  subroutine configure_endpoint()
    select case (trim(adjustl(endpoint_policy)))
    case ("local")
      free_vortex_endpoint%enabled = .false.
      endpoint_policy = "local"
    case ("free_vortex")
      free_vortex_endpoint%enabled = .true.
      endpoint_policy = "free_vortex"
    case default
      error stop "endpoint_policy must be 'local' or 'free_vortex'"
    end select
    free_vortex_endpoint%center = &
      [asymptotic_center_x, asymptotic_center_y]
    free_vortex_endpoint%bulk_gap = asymptotic_bulk_gap
    free_vortex_endpoint%winding = asymptotic_winding
    free_vortex_endpoint%outer_radius = asymptotic_outer_radius
    free_vortex_endpoint%maximum_step = asymptotic_maximum_step
  end subroutine configure_endpoint


end program fermiforge_2d_mpi_smoke
