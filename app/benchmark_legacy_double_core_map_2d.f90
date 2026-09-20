program benchmark_legacy_double_core_map_2d
  use mpi_f08
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, &
    make_rectilinear_cartesian_mesh, make_uniform_cartesian_mesh, &
    make_symmetric_multiscale_cartesian_mesh
  use spinful_state_2d, only : spinful_state_2d_t, &
    allocate_spinful_state_2d
  use spinful_field_io_2d, only : read_spinful_field_map_2d, &
    write_spinful_field_map_2d
  use spinful_sparse_state_io_2d, only : &
    read_sparse_spinful_state_2d, write_sparse_spinful_state_2d
  use spinful_tangent_projection_2d, only : &
    zero_mode_tangent_basis_2d_t, zero_mode_projection_report_t, &
    build_zero_mode_tangent_basis_2d, &
    project_spinful_residual_away_from_tangents_2d
  use legacy_split_field_io_2d, only : &
    legacy_split_field_report_2d_t, &
    read_legacy_split_field_directory_2d
  use double_core_seed_2d, only : double_core_seed_report_t, &
    initialize_regularized_london_double_core_2d
  use historical_core_seed_2d, only : historical_core_seed_report_t, &
    initialize_historical_core_seed_2d
  use he3_quadrature, only : angular_quadrature_3d_t, ozaki_quadrature_t, &
    read_legacy_gauss_table, read_legacy_ozaki_table
  use he3_mpi_field_map_2d, only : he3_mpi_field_map_diagnostics_t, &
    evaluate_mpi_he3_field_map, update_mpi_he3_state_with_anderson
  use he3_nonlinear_iteration_2d, only : he3_field_residual_report_t, &
    compute_he3_field_residual
  use free_vortex_asymptotic_2d, only : free_vortex_endpoint_2d_t, &
    fit_state_bulk_reference, set_state_asymptotic_fit, &
    set_state_asymptotic_elliptical_fit, &
    apply_free_vortex_asymptotic_halo_2d
  use new_src_iteration_layout_2d, only : &
    new_src_masked_iteration_vector_size_2d
  use legacy_anderson_mixing, only : legacy_anderson_t, anderson_report_t
  implicit none

  type(angular_quadrature_3d_t) :: angular
  type(cartesian_mesh_2d_t) :: mesh
  type(he3_field_residual_report_t) :: core_residual, outer_residual
  type(he3_field_residual_report_t) :: residual, update_residual
  type(he3_field_residual_report_t) :: shadow_probe_residual
  type(he3_mpi_field_map_diagnostics_t) :: diagnostics
  type(he3_mpi_field_map_diagnostics_t) :: passive_probe_diagnostics
  type(he3_mpi_field_map_diagnostics_t) :: shadow_probe_diagnostics
  type(anderson_report_t) :: anderson_report
  type(anderson_report_t) :: shadow_probe_anderson_report
  type(legacy_anderson_t) :: accelerator, shadow_probe_accelerator
  type(legacy_split_field_report_2d_t) :: import_report
  type(double_core_seed_report_t) :: seed_report
  type(historical_core_seed_report_t) :: historical_seed_report
  type(ozaki_quadrature_t) :: energy
  type(free_vortex_endpoint_2d_t) :: endpoint
  type(spinful_state_2d_t) :: evaluated_state, mapped, next
  type(spinful_state_2d_t) :: passive_probe_baseline, passive_probe_mapped
  type(spinful_state_2d_t) :: shadow_probe_baseline
  type(spinful_state_2d_t) :: shadow_probe_first_mapped
  type(spinful_state_2d_t) :: shadow_probe_mapped
  type(spinful_state_2d_t) :: shadow_probe_next, shadow_probe_state
  type(spinful_state_2d_t) :: projected_mapped, state
  type(zero_mode_tangent_basis_2d_t) :: zero_mode_basis
  type(zero_mode_projection_report_t) :: zero_mode_report

  character(len=1024) :: benchmark_directory, gauss_file, input_file
  character(len=1024) :: input_message, metrics_file, ozaki_file
  character(len=1024) :: asymptotic_probe_file
  character(len=1024) :: asymptotic_probe_history_file
  character(len=1024) :: asymptotic_shadow_file
  character(len=1024) :: output_override, sample_file
  character(len=1024) :: final_field_file, history_file, restart_field_file
  character(len=1024) :: sparse_final_file, sparse_restart_file
  character(len=1024) :: checkpoint_file, initial_field_file
  character(len=32) :: bulk_reference_policy, checkpoint_scope
  character(len=32) :: endpoint_policy, inactive_halo_policy
  character(len=32) :: initialization_mode, terminal_status
  character(len=32) :: passive_probe_source
  character(len=32) :: shadow_probe_terminal_status
  character(len=32) :: scratch_mesh_kind
  character(len=32) :: update_region, zero_mode_projection
  integer :: azimuth_count, ierr, input_status, input_unit
  integer :: anderson_history_limit, history_unit
  integer :: asymptotic_probe_iteration_stencil_radius
  integer :: asymptotic_probe_iterations
  integer :: asymptotic_probe_anderson_history_limit
  integer :: asymptotic_ray_probe_count
  integer :: bulk_reference_sample_count, center_probe_stencil_radius
  integer :: checkpoint_interval, iteration, iterations_completed, map_evaluations
  integer :: maximum_iterations, polar_count, probe_stencil_radius
  integer :: passive_probe_map_evaluations
  integer :: shadow_probe_history_unit, shadow_probe_iterations_completed
  integer :: shadow_probe_map_evaluations
  integer :: scratch_number_of_cells
  integer :: rank, rank_count, sparse_restart_applied_points
  integer :: result_integer
  real(rk) :: asymptotic_bulk_gap, asymptotic_fit_inner_radius
  real(rk) :: asymptotic_matching_radius, asymptotic_outer_radius
  real(rk) :: asymptotic_fit_inner_radius_x
  real(rk) :: asymptotic_fit_inner_radius_y
  real(rk) :: asymptotic_matching_radius_x
  real(rk) :: asymptotic_matching_radius_y
  real(rk) :: asymptotic_probe_convergence_tolerance
  real(rk) :: asymptotic_probe_anderson_maximum_mixing
  real(rk) :: asymptotic_fit_inner_radii(2), asymptotic_matching_radii(2)
  real(rk) :: anderson_maximum_mixing, anderson_progress_threshold
  real(rk) :: boundary_relaxation_distance, feedback_scale
  real(rk) :: bulk_reference_fit_maximum, bulk_reference_fit_radius
  real(rk) :: bulk_reference_fit_rms
  real(rk) :: convergence_tolerance
  real(rk) :: cumulative_map_seconds, maximum_normalization_error_over_solve
  real(rk) :: passive_probe_map_seconds
  real(rk) :: passive_probe_maximum_normalization_error
  real(rk) :: shadow_probe_cumulative_map_seconds
  real(rk) :: shadow_probe_maximum_normalization_error
  real(rk) :: global_probe_radius, half_core_offset_x, half_core_offset_y
  real(rk) :: measured_center_pair_amplitude
  real(rk) :: measured_half_core_separation
  real(rk) :: measured_negative_half_core_amplitude
  real(rk) :: measured_negative_half_core_y
  real(rk) :: measured_positive_half_core_amplitude
  real(rk) :: measured_positive_half_core_y
  real(rk) :: initial_halo_maximum_change, last_halo_maximum_change
  real(rk) :: update_radius, update_radius_x, update_radius_y
  real(rk) :: scratch_half_width, scratch_fine_region_half_width
  real(rk) :: scratch_medium_region_half_width, scratch_fine_spacing
  real(rk) :: scratch_medium_spacing, scratch_coarse_spacing
  real(rk) :: seed_core_width, seed_domain_wall_width
  real(rk) :: seed_reduced_temperature
  real(rk) :: trajectory_maximum_step
  integer :: minimum_substep_count
  integer, allocatable :: asymptotic_ray_id(:)
  logical, allocatable :: active_point(:), asymptotic_ray_point(:)
  logical, allocatable :: asymptotic_probe_iteration_point(:)
  logical, allocatable :: asymptotic_source_point(:)
  logical, allocatable :: checkpoint_point(:), core_point(:)
  logical, allocatable :: dependent_halo_point(:)
  logical, allocatable :: outer_point(:), update_point(:)
  logical :: benchmark_passed, map_succeeded, pin_center, require_convergence
  logical :: passive_probe_map_succeeded, passive_probe_output_requested
  logical :: shadow_probe_map_succeeded
  logical :: solve_requested, update_outer_probes, use_zero_mode_projection

  namelist /legacy_double_core_map_benchmark/ benchmark_directory, &
    gauss_file, ozaki_file, azimuth_count, polar_count, feedback_scale, &
    trajectory_maximum_step, minimum_substep_count, &
    boundary_relaxation_distance, global_probe_radius, &
    asymptotic_ray_probe_count, asymptotic_probe_file, &
    asymptotic_probe_iterations, asymptotic_probe_iteration_stencil_radius, &
    asymptotic_probe_convergence_tolerance, &
    asymptotic_probe_anderson_history_limit, &
    asymptotic_probe_anderson_maximum_mixing, &
    asymptotic_probe_history_file, &
    asymptotic_shadow_file, &
    half_core_offset_x, half_core_offset_y, center_probe_stencil_radius, &
    probe_stencil_radius, endpoint_policy, asymptotic_bulk_gap, &
    asymptotic_fit_inner_radius, asymptotic_matching_radius, &
    asymptotic_fit_inner_radius_x, asymptotic_fit_inner_radius_y, &
    asymptotic_matching_radius_x, asymptotic_matching_radius_y, &
    asymptotic_outer_radius, bulk_reference_fit_radius, &
    bulk_reference_sample_count, bulk_reference_policy, &
    initialization_mode, restart_field_file, &
    sparse_restart_file, sparse_final_file, &
    scratch_mesh_kind, scratch_half_width, scratch_number_of_cells, &
    scratch_fine_region_half_width, scratch_medium_region_half_width, &
    scratch_fine_spacing, scratch_medium_spacing, scratch_coarse_spacing, &
    seed_core_width, seed_domain_wall_width, seed_reduced_temperature, &
    inactive_halo_policy, &
    maximum_iterations, convergence_tolerance, anderson_history_limit, &
    anderson_progress_threshold, anderson_maximum_mixing, &
    update_region, update_radius, update_radius_x, update_radius_y, &
    update_outer_probes, pin_center, &
    zero_mode_projection, &
    require_convergence, &
    checkpoint_interval, checkpoint_file, checkpoint_scope, &
    initial_field_file, final_field_file, history_file, &
    metrics_file, sample_file

  benchmark_directory = &
    "2D_benchmarks/Double_core_vortex_T=0.30_Fs1=5.4"
  gauss_file = trim(benchmark_directory) // "/gauss11.dat"
  ozaki_file = trim(benchmark_directory) // "/ozaki.dat"
  azimuth_count = 64
  polar_count = 11
  feedback_scale = 5.4_rk / (1.0_rk + 5.4_rk / 3.0_rk)
  trajectory_maximum_step = 0.25_rk
  minimum_substep_count = 1
  boundary_relaxation_distance = 1.75_rk
  global_probe_radius = 40.0_rk
  asymptotic_ray_probe_count = 0
  asymptotic_probe_file = ""
  asymptotic_probe_iterations = 0
  asymptotic_probe_iteration_stencil_radius = 0
  asymptotic_probe_convergence_tolerance = 1.0e-7_rk
  asymptotic_probe_anderson_history_limit = 5
  asymptotic_probe_anderson_maximum_mixing = 1.0_rk
  asymptotic_probe_history_file = ""
  asymptotic_shadow_file = ""
  half_core_offset_x = 0.0_rk
  half_core_offset_y = 11.8_rk
  center_probe_stencil_radius = 0
  probe_stencil_radius = 1
  endpoint_policy = "state_asymptotic"
  asymptotic_bulk_gap = 0.2799_rk
  asymptotic_fit_inner_radius = 50.0_rk
  asymptotic_matching_radius = 58.0_rk
  asymptotic_fit_inner_radius_x = 0.0_rk
  asymptotic_fit_inner_radius_y = 0.0_rk
  asymptotic_matching_radius_x = 0.0_rk
  asymptotic_matching_radius_y = 0.0_rk
  asymptotic_outer_radius = 100.0_rk
  bulk_reference_fit_radius = 58.0_rk
  bulk_reference_sample_count = 128
  bulk_reference_policy = "fit_initial"
  bulk_reference_fit_rms = 0.0_rk
  bulk_reference_fit_maximum = 0.0_rk
  initialization_mode = "legacy_split"
  restart_field_file = ""
  sparse_restart_file = ""
  sparse_final_file = ""
  scratch_mesh_kind = "multiscale"
  scratch_half_width = 12.0_rk
  scratch_number_of_cells = 48
  scratch_fine_region_half_width = 4.0_rk
  scratch_medium_region_half_width = 8.0_rk
  scratch_fine_spacing = 0.5_rk
  scratch_medium_spacing = 1.0_rk
  scratch_coarse_spacing = 2.0_rk
  seed_core_width = 0.8_rk
  seed_domain_wall_width = 0.6_rk
  seed_reduced_temperature = 0.30_rk
  inactive_halo_policy = "frozen"
  maximum_iterations = 0
  convergence_tolerance = 1.0e-6_rk
  anderson_history_limit = 10
  anderson_progress_threshold = 0.1_rk
  anderson_maximum_mixing = 5.0_rk
  update_region = "half_core_stencils"
  update_radius = 0.4_rk
  update_radius_x = 0.4_rk
  update_radius_y = 0.4_rk
  update_outer_probes = .false.
  pin_center = .true.
  zero_mode_projection = "none"
  require_convergence = .false.
  checkpoint_interval = 0
  checkpoint_file = ""
  checkpoint_scope = "update_region"
  initial_field_file = ""
  final_field_file = ""
  history_file = "work/double-core-one-map/iteration_history.dat"
  metrics_file = "work/double-core-one-map/metrics.txt"
  sample_file = "work/double-core-one-map/sampled_residuals.dat"

  call MPI_Init(ierr)
  call require_mpi(ierr == MPI_SUCCESS, "MPI initialization failed")
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call require_mpi(ierr == MPI_SUCCESS, "MPI rank query failed")
  call MPI_Comm_size(MPI_COMM_WORLD, rank_count, ierr)
  call require_mpi(ierr == MPI_SUCCESS, "MPI size query failed")

  call get_command_argument(1, input_file)
  if (len_trim(input_file) == 0) &
    input_file = "examples/2d_double_core_one_map.nml"
  open(newunit=input_unit, file=trim(input_file), status="old", &
       action="read", iostat=input_status, iomsg=input_message)
  if (input_status /= 0) &
    error stop "cannot open double-core benchmark input: " // &
      trim(input_message)
  read(input_unit, nml=legacy_double_core_map_benchmark, &
       iostat=input_status, iomsg=input_message)
  close(input_unit)
  if (input_status /= 0) &
    error stop "cannot read double-core benchmark input: " // &
      trim(input_message)
  call get_command_argument(2, output_override)
  if (len_trim(output_override) > 0) metrics_file = output_override
  call get_command_argument(3, output_override)
  if (len_trim(output_override) > 0) sample_file = output_override
  call validate_scalar_input()

  call import_and_broadcast_state()
  call validate_mesh_controls()
  call configure_endpoint()
  call read_legacy_gauss_table( &
    trim(gauss_file), azimuth_count, polar_count, angular)
  call read_legacy_ozaki_table(trim(ozaki_file), energy)
  call make_probe_masks()
  call make_asymptotic_halo_masks()
  call make_asymptotic_probe_iteration_mask()
  call refresh_dependent_halo(state, initial_halo_maximum_change)
  last_halo_maximum_change = initial_halo_maximum_change
  if (rank == 0) call measure_current_core_shape()
  if (rank == 0 .and. len_trim(initial_field_file) > 0) &
    call write_spinful_field_map_2d( &
      trim(initial_field_file), mesh, state, &
      "double_core_initial_state", trim(endpoint_policy))
  call configure_zero_mode_basis()

  solve_requested = maximum_iterations > 0
  terminal_status = "single_map"
  iterations_completed = 0
  map_evaluations = 0
  cumulative_map_seconds = 0.0_rk
  maximum_normalization_error_over_solve = 0.0_rk
  zero_mode_report = zero_mode_projection_report_t()
  shadow_probe_terminal_status = "disabled"
  shadow_probe_iterations_completed = 0
  shadow_probe_map_evaluations = 0
  shadow_probe_cumulative_map_seconds = 0.0_rk
  shadow_probe_maximum_normalization_error = 0.0_rk
  shadow_probe_map_succeeded = .true.
  passive_probe_output_requested = asymptotic_ray_probe_count > 0 .and. &
    len_trim(asymptotic_probe_file) > 0
  passive_probe_source = "disabled"
  passive_probe_map_evaluations = 0
  passive_probe_map_seconds = 0.0_rk
  passive_probe_maximum_normalization_error = 0.0_rk
  passive_probe_map_succeeded = .true.
  if (rank == 0) then
    call open_iteration_history()
    if (solve_requested) call accelerator%initialize( &
      new_src_masked_iteration_vector_size_2d(state, update_point), &
      anderson_history_limit, anderson_progress_threshold, &
      anderson_maximum_mixing)
  end if

  do iteration = 1, max(1, maximum_iterations)
    evaluated_state = state
    if (trim(endpoint_policy) == "state_asymptotic") then
      call evaluate_mpi_he3_field_map( &
        MPI_COMM_WORLD, mesh, state, angular, energy, feedback_scale, &
        trajectory_maximum_step, minimum_substep_count, &
        boundary_relaxation_distance, mapped, diagnostics, map_succeeded, &
        free_vortex_endpoint=endpoint, active_point_mask=active_point)
    else
      call evaluate_mpi_he3_field_map( &
        MPI_COMM_WORLD, mesh, state, angular, energy, feedback_scale, &
        trajectory_maximum_step, minimum_substep_count, &
        boundary_relaxation_distance, mapped, diagnostics, map_succeeded, &
        active_point_mask=active_point)
    end if
    map_evaluations = map_evaluations + 1
    cumulative_map_seconds = cumulative_map_seconds + &
      diagnostics%maximum_rank_elapsed_seconds
    maximum_normalization_error_over_solve = max( &
      maximum_normalization_error_over_solve, &
      diagnostics%global%maximum_normalization_error)
    if (.not. map_succeeded) then
      terminal_status = "map_failure"
      exit
    end if
    if (rank == 0) call compute_current_residuals()
    if (.not. solve_requested) then
      if (rank == 0) call write_single_map_history(iteration)
      exit
    end if

    if (use_zero_mode_projection) then
      if (rank == 0) then
        call project_spinful_residual_away_from_tangents_2d( &
          state, mapped, update_point, zero_mode_basis, &
          projected_mapped, zero_mode_report)
        call update_mpi_he3_state_with_anderson( &
          MPI_COMM_WORLD, mesh, accelerator, state, projected_mapped, &
          convergence_tolerance, next, anderson_report, update_residual, &
          update_point)
      else
        call update_mpi_he3_state_with_anderson( &
          MPI_COMM_WORLD, mesh, accelerator, state, mapped, &
          convergence_tolerance, next, anderson_report, update_residual, &
          update_point)
      end if
    else
      call update_mpi_he3_state_with_anderson( &
        MPI_COMM_WORLD, mesh, accelerator, state, mapped, &
        convergence_tolerance, next, anderson_report, update_residual, &
        update_point)
    end if
    iterations_completed = iteration
    state = next
    call refresh_endpoint_reference()
    call refresh_dependent_halo(state, last_halo_maximum_change)
    if (rank == 0) then
      call measure_current_core_shape()
      call write_iteration_history(iteration)
      call print_iteration(iteration)
    end if
    if (rank == 0 .and. checkpoint_interval > 0) then
      if (mod(iteration, checkpoint_interval) == 0) &
        call write_sparse_spinful_state_2d( &
          trim(checkpoint_file), mesh, state, checkpoint_point, &
          "double_core_iteration_checkpoint")
    end if
    if (anderson_report%converged) then
      terminal_status = "converged"
      exit
    end if
    if (iteration == maximum_iterations) terminal_status = "iteration_limit"
  end do

  if (map_succeeded .and. asymptotic_probe_iterations > 0) &
    call iterate_asymptotic_probe_shadow()
  if (map_succeeded .and. passive_probe_output_requested .and. &
      asymptotic_probe_iterations == 0) &
    call evaluate_final_passive_probe_map()

  benchmark_passed = map_succeeded .and. shadow_probe_map_succeeded .and. &
    passive_probe_map_succeeded .and. &
    maximum_normalization_error_over_solve <= 1.0e-10_rk
  if (asymptotic_probe_iterations > 0) benchmark_passed = &
    benchmark_passed .and. &
    shadow_probe_maximum_normalization_error <= 1.0e-10_rk
  if (passive_probe_output_requested) benchmark_passed = benchmark_passed .and. &
    passive_probe_maximum_normalization_error <= 1.0e-10_rk
  if (require_convergence .and. solve_requested) benchmark_passed = &
    benchmark_passed .and. trim(terminal_status) == "converged"
  if (rank == 0 .and. map_succeeded) then
    close(history_unit)
    call write_metrics()
    call write_sampled_residuals()
    if (passive_probe_output_requested .and. passive_probe_map_succeeded) &
      call write_asymptotic_probe_states()
    if (asymptotic_probe_iterations > 0 .and. &
        shadow_probe_map_succeeded) &
      call write_asymptotic_shadow_states()
    if (len_trim(final_field_file) > 0) call write_spinful_field_map_2d( &
      trim(final_field_file), mesh, state, &
      "double_core_final_state", trim(endpoint_policy))
    if (len_trim(sparse_final_file) > 0) &
      call write_sparse_spinful_state_2d( &
        trim(sparse_final_file), mesh, state, checkpoint_point, &
        "double_core_sparse_final_state")
    if (checkpoint_interval > 0 .and. len_trim(checkpoint_file) > 0) &
      call write_sparse_spinful_state_2d( &
        trim(checkpoint_file), mesh, state, checkpoint_point, &
        "double_core_iteration_checkpoint")
    call print_summary()
  end if

  result_integer = merge(1, 0, benchmark_passed)
  call MPI_Bcast(result_integer, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call require_mpi(ierr == MPI_SUCCESS, &
    "double-core benchmark result broadcast failed")
  benchmark_passed = result_integer == 1
  call MPI_Finalize(ierr)
  call require_mpi(ierr == MPI_SUCCESS, "MPI finalization failed")
  if (.not. benchmark_passed) &
    error stop "double-core transport/iteration evaluation failed"

contains

  subroutine validate_scalar_input()
    logical :: any_elliptical_radius, all_elliptical_radii

    asymptotic_fit_inner_radii = &
      [asymptotic_fit_inner_radius_x, asymptotic_fit_inner_radius_y]
    asymptotic_matching_radii = &
      [asymptotic_matching_radius_x, asymptotic_matching_radius_y]
    any_elliptical_radius = any(asymptotic_fit_inner_radii /= 0.0_rk) .or. &
      any(asymptotic_matching_radii /= 0.0_rk)
    all_elliptical_radii = all(asymptotic_fit_inner_radii > 0.0_rk) .and. &
      all(asymptotic_matching_radii > 0.0_rk)
    if (len_trim(gauss_file) == 0 .or. len_trim(ozaki_file) == 0) &
      error stop "double-core quadrature paths cannot be blank"
    if (azimuth_count < 1 .or. polar_count < 1) &
      error stop "double-core quadrature counts must be positive"
    if (feedback_scale < 0.0_rk .or. trajectory_maximum_step <= 0.0_rk .or. &
        minimum_substep_count < 1 .or. &
        boundary_relaxation_distance < 0.0_rk) &
      error stop "invalid double-core transport controls"
    if (global_probe_radius <= 0.0_rk .or. &
        asymptotic_ray_probe_count < 0 .or. &
        asymptotic_probe_iterations < 0 .or. &
        asymptotic_probe_iteration_stencil_radius < 0 .or. &
        asymptotic_probe_convergence_tolerance <= 0.0_rk .or. &
        asymptotic_probe_anderson_history_limit < 1 .or. &
        asymptotic_probe_anderson_maximum_mixing < 0.01_rk .or. &
        center_probe_stencil_radius < 0 .or. probe_stencil_radius < 0) &
      error stop "invalid double-core probe controls"
    if (asymptotic_probe_iterations > 0 .and. &
        (asymptotic_ray_probe_count == 0 .or. &
         trim(inactive_halo_policy) /= "state_asymptotic" .or. &
         trim(endpoint_policy) /= "state_asymptotic" .or. &
         len_trim(asymptotic_probe_history_file) == 0 .or. &
         len_trim(asymptotic_shadow_file) == 0)) &
      error stop "shadow probe iteration needs ray probes and an asymptotic halo"
    if (trim(endpoint_policy) /= "local_supplied_field" .and. &
        trim(endpoint_policy) /= "state_asymptotic") &
      error stop "unknown double-core endpoint policy"
    if (trim(inactive_halo_policy) /= "frozen" .and. &
        trim(inactive_halo_policy) /= "state_asymptotic") &
      error stop "unknown double-core inactive-halo policy"
    if (trim(inactive_halo_policy) == "state_asymptotic" .and. &
        trim(endpoint_policy) /= "state_asymptotic") &
      error stop "dependent halo requires the state-asymptotic endpoint"
    if (trim(endpoint_policy) == "state_asymptotic" .and. &
        (asymptotic_bulk_gap <= 0.0_rk .or. &
         asymptotic_fit_inner_radius <= 0.0_rk .or. &
         asymptotic_matching_radius <= asymptotic_fit_inner_radius .or. &
         asymptotic_outer_radius <= asymptotic_matching_radius .or. &
         bulk_reference_fit_radius <= 0.0_rk .or. &
         bulk_reference_sample_count < 8)) &
      error stop "invalid double-core asymptotic controls"
    if (trim(bulk_reference_policy) /= "theoretical" .and. &
        trim(bulk_reference_policy) /= "fit_initial") &
      error stop "bulk_reference_policy must be theoretical or fit_initial"
    if (any_elliptical_radius .neqv. all_elliptical_radii) &
      error stop "all four elliptical asymptotic radii must be supplied"
    if (all_elliptical_radii .and. &
        trim(update_region) /= "centered_ellipse") &
      error stop "elliptical asymptotics require a centered-ellipse update"
    if (all_elliptical_radii .and. &
        (any(asymptotic_matching_radii <= asymptotic_fit_inner_radii) .or. &
         asymptotic_outer_radius <= maxval(asymptotic_matching_radii))) &
      error stop "invalid nested elliptical asymptotic surfaces"
    if (trim(initialization_mode) /= "legacy_split" .and. &
        trim(initialization_mode) /= "restart" .and. &
        trim(initialization_mode) /= "regularized_london" .and. &
        trim(initialization_mode) /= "historical_nop" .and. &
        trim(initialization_mode) /= "historical_aop" .and. &
        trim(initialization_mode) /= "historical_dop") &
      error stop "unknown double-core initialization mode"
    if ((trim(initialization_mode) == "legacy_split" .or. &
         trim(initialization_mode) == "restart") .and. &
        len_trim(benchmark_directory) == 0) &
      error stop "legacy double-core initialization needs a benchmark directory"
    if (trim(initialization_mode) == "restart" .and. &
        len_trim(restart_field_file) == 0) &
      error stop "double-core restart needs a field file"
    if (trim(initialization_mode) == "regularized_london" .or. &
        trim(initialization_mode) == "historical_nop" .or. &
        trim(initialization_mode) == "historical_aop" .or. &
        trim(initialization_mode) == "historical_dop") then
      select case (trim(scratch_mesh_kind))
      case ("uniform")
        if (scratch_half_width <= 0.0_rk .or. &
            scratch_number_of_cells < 2 .or. &
            modulo(scratch_number_of_cells, 2) /= 0) &
          error stop "scratch uniform mesh needs a positive width and even cells"
      case ("multiscale")
        if (scratch_half_width <= 0.0_rk .or. &
            scratch_fine_region_half_width <= 0.0_rk .or. &
            scratch_medium_region_half_width <= &
              scratch_fine_region_half_width .or. &
            scratch_half_width <= scratch_medium_region_half_width .or. &
            scratch_fine_spacing <= 0.0_rk .or. &
            scratch_medium_spacing < scratch_fine_spacing .or. &
            scratch_coarse_spacing < scratch_medium_spacing) &
          error stop "invalid scratch multiscale mesh controls"
      case default
        error stop "scratch_mesh_kind must be uniform or multiscale"
      end select
      if (abs(half_core_offset_y) >= scratch_half_width) &
        error stop "scratch core probes lie outside the generated mesh"
    end if
    if (trim(initialization_mode) == "regularized_london") then
      if (abs(half_core_offset_x) > 128.0_rk * epsilon(1.0_rk) * &
          max(1.0_rk, abs(half_core_offset_y))) &
        error stop "regularized London seed currently places half cores on y"
      if (abs(half_core_offset_y) <= 0.0_rk .or. &
          seed_core_width <= 0.0_rk .or. seed_domain_wall_width <= 0.0_rk) &
        error stop "invalid regularized London seed controls"
    end if
    if ((trim(initialization_mode) == "historical_nop" .or. &
         trim(initialization_mode) == "historical_aop" .or. &
         trim(initialization_mode) == "historical_dop") .and. &
        (seed_reduced_temperature < 0.0_rk .or. &
         seed_reduced_temperature >= 1.0_rk)) &
      error stop "historical seed reduced temperature must be in [0,1)"
    if (maximum_iterations < 0 .or. convergence_tolerance <= 0.0_rk) &
      error stop "invalid double-core iteration controls"
    if (checkpoint_interval < 0 .or. &
        (checkpoint_interval > 0 .and. len_trim(checkpoint_file) == 0)) &
      error stop "invalid double-core checkpoint controls"
    if (trim(checkpoint_scope) /= "update_region" .and. &
        trim(checkpoint_scope) /= "full_state") &
      error stop "checkpoint_scope must be update_region or full_state"
    if (anderson_history_limit < 1 .or. &
        anderson_progress_threshold <= 0.0_rk .or. &
        anderson_progress_threshold > 1.0_rk .or. &
        anderson_maximum_mixing < 0.01_rk) &
      error stop "invalid double-core Anderson controls"
    if (trim(update_region) /= "half_core_stencils" .and. &
        trim(update_region) /= "half_core_disks" .and. &
        trim(update_region) /= "centered_disk" .and. &
        trim(update_region) /= "centered_ellipse") &
      error stop "unknown double-core update region"
    if ((trim(update_region) == "half_core_disks" .or. &
         trim(update_region) == "centered_disk") .and. &
        update_radius <= 0.0_rk) &
      error stop "double-core disk update radius must be positive"
    if (trim(update_region) == "centered_ellipse" .and. &
        (update_radius_x <= 0.0_rk .or. update_radius_y <= 0.0_rk)) &
      error stop "double-core ellipse update radii must be positive"
    if (all_elliptical_radii .and. &
        any(asymptotic_matching_radii >= &
          [update_radius_x, update_radius_y])) &
      error stop "asymptotic matching ellipse must lie inside update ellipse"
    if (trim(zero_mode_projection) /= "none" .and. &
        trim(zero_mode_projection) /= "translations" .and. &
        trim(zero_mode_projection) /= "translation_orientation" .and. &
        trim(zero_mode_projection) /= &
          "gauge_translation_orientation") &
      error stop "unknown double-core zero-mode projection"
    if (len_trim(metrics_file) == 0 .or. len_trim(sample_file) == 0 .or. &
        len_trim(history_file) == 0) &
      error stop "double-core output paths cannot be blank"
  end subroutine validate_scalar_input


  subroutine validate_mesh_controls()
    real(rk) :: edge_distance

    edge_distance = min( &
      abs(mesh%x_minimum), abs(mesh%x_maximum), &
      abs(mesh%y_minimum), abs(mesh%y_maximum))
    if (global_probe_radius >= edge_distance) &
      error stop "global double-core probes must lie inside the field boundary"
    if (abs(half_core_offset_x) >= edge_distance .or. &
        abs(half_core_offset_y) >= edge_distance) &
      error stop "double-core half-core probes lie outside the field mesh"
    if (asymptotic_ray_probe_count > 0) then
      if (trim(update_region) == "centered_ellipse" .and. &
          (update_radius_x >= min(abs(mesh%x_minimum), abs(mesh%x_maximum)) .or. &
           update_radius_y >= min(abs(mesh%y_minimum), abs(mesh%y_maximum)))) &
        error stop "asymptotic ray probes need mesh space outside the active ellipse"
      if (trim(update_region) == "centered_disk" .and. &
          update_radius >= edge_distance) &
        error stop "asymptotic ray probes need mesh space outside the active disk"
    end if
  end subroutine validate_mesh_controls


  subroutine import_and_broadcast_state()
    integer :: dimensions(2), node
    real(rk), allocatable :: x_coordinates(:), y_coordinates(:)

    if (rank == 0) then
      select case (trim(initialization_mode))
      case ("regularized_london")
        call make_scratch_mesh()
        call initialize_regularized_london_double_core_2d( &
          mesh, asymptotic_bulk_gap, abs(half_core_offset_y), &
          seed_core_width, seed_domain_wall_width, state, seed_report)
        import_report = legacy_split_field_report_2d_t()
        historical_seed_report = historical_core_seed_report_t()
      case ("historical_nop", "historical_aop", "historical_dop")
        call make_scratch_mesh()
        call initialize_historical_core_seed_2d( &
          mesh, initialization_mode(12:14), asymptotic_bulk_gap, &
          seed_reduced_temperature, feedback_scale, state, &
          historical_seed_report)
        seed_report = double_core_seed_report_t()
        import_report = legacy_split_field_report_2d_t()
      case default
        call read_legacy_split_field_directory_2d( &
          trim(benchmark_directory), mesh, state, import_report)
        seed_report = double_core_seed_report_t()
        historical_seed_report = historical_core_seed_report_t()
      end select
      sparse_restart_applied_points = 0
      if (trim(initialization_mode) == "restart") &
        call read_spinful_field_map_2d(trim(restart_field_file), mesh, state)
      if (len_trim(sparse_restart_file) > 0) &
        call read_sparse_spinful_state_2d( &
          trim(sparse_restart_file), mesh, state, &
          sparse_restart_applied_points)
      dimensions = [mesh%x_point_count(), mesh%y_point_count()]
    else
      dimensions = 0
    end if
    call MPI_Bcast(dimensions, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call require_mpi(ierr == MPI_SUCCESS, &
      "could not broadcast double-core mesh dimensions")
    allocate(x_coordinates(dimensions(1)), y_coordinates(dimensions(2)))
    if (rank == 0) then
      do node = 0, mesh%x_cell_count()
        x_coordinates(node + 1) = mesh%x_coordinate(node)
      end do
      do node = 0, mesh%y_cell_count()
        y_coordinates(node + 1) = mesh%y_coordinate(node)
      end do
    end if
    call MPI_Bcast(x_coordinates, dimensions(1), MPI_DOUBLE_PRECISION, &
      0, MPI_COMM_WORLD, ierr)
    call require_mpi(ierr == MPI_SUCCESS, &
      "could not broadcast double-core x coordinates")
    call MPI_Bcast(y_coordinates, dimensions(2), MPI_DOUBLE_PRECISION, &
      0, MPI_COMM_WORLD, ierr)
    call require_mpi(ierr == MPI_SUCCESS, &
      "could not broadcast double-core y coordinates")
    if (rank /= 0) then
      call make_rectilinear_cartesian_mesh( &
        x_coordinates, y_coordinates, mesh)
      call allocate_spinful_state_2d(mesh, state)
    end if
    call MPI_Bcast(state%order_parameter, size(state%order_parameter), &
      MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
    call require_mpi(ierr == MPI_SUCCESS, &
      "could not broadcast double-core order parameter")
    call MPI_Bcast(state%current_mean_field, &
      size(state%current_mean_field), MPI_DOUBLE_PRECISION, 0, &
      MPI_COMM_WORLD, ierr)
    call require_mpi(ierr == MPI_SUCCESS, &
      "could not broadcast double-core mean field")
  end subroutine import_and_broadcast_state


  subroutine make_scratch_mesh()
    select case (trim(scratch_mesh_kind))
    case ("uniform")
      call make_uniform_cartesian_mesh( &
        -scratch_half_width, scratch_half_width, scratch_number_of_cells, &
        -scratch_half_width, scratch_half_width, scratch_number_of_cells, &
        mesh)
    case ("multiscale")
      call make_symmetric_multiscale_cartesian_mesh( &
        scratch_half_width, scratch_fine_region_half_width, &
        scratch_medium_region_half_width, scratch_fine_spacing, &
        scratch_medium_spacing, scratch_coarse_spacing, mesh)
    end select
  end subroutine make_scratch_mesh


  subroutine configure_endpoint()
    if (trim(endpoint_policy) == "local_supplied_field") return

    endpoint%enabled = .true.
    endpoint%center = 0.0_rk
    endpoint%bulk_gap = asymptotic_bulk_gap
    endpoint%winding = 1.0_rk
    endpoint%outer_radius = asymptotic_outer_radius
    endpoint%maximum_step = trajectory_maximum_step
    if (trim(bulk_reference_policy) == "fit_initial") then
      if (use_elliptical_asymptotic_fit()) then
        call fit_state_bulk_reference( &
          endpoint, mesh, state, bulk_reference_fit_radius, &
          bulk_reference_sample_count, bulk_reference_fit_rms, &
          bulk_reference_fit_maximum, fitting_radii=asymptotic_matching_radii)
      else
        call fit_state_bulk_reference( &
          endpoint, mesh, state, bulk_reference_fit_radius, &
          bulk_reference_sample_count, bulk_reference_fit_rms, &
          bulk_reference_fit_maximum)
      end if
    else if (use_elliptical_asymptotic_fit()) then
      call fit_state_bulk_reference( &
        endpoint, mesh, state, bulk_reference_fit_radius, &
        bulk_reference_sample_count, bulk_reference_fit_rms, &
        bulk_reference_fit_maximum, fitting_radii=asymptotic_matching_radii, &
        update_reference=.false.)
    else
      call fit_state_bulk_reference( &
        endpoint, mesh, state, bulk_reference_fit_radius, &
        bulk_reference_sample_count, bulk_reference_fit_rms, &
        bulk_reference_fit_maximum, update_reference=.false.)
    end if
    if (use_elliptical_asymptotic_fit()) then
      call set_state_asymptotic_elliptical_fit( &
        endpoint, asymptotic_fit_inner_radii, asymptotic_matching_radii)
    else
      call set_state_asymptotic_fit( &
        endpoint, asymptotic_fit_inner_radius, asymptotic_matching_radius)
    end if
    if (.not. endpoint%is_valid_for(mesh)) &
      error stop "configured double-core endpoint is invalid for its mesh"
  end subroutine configure_endpoint


  subroutine make_probe_masks()
    integer :: x_sign, y_sign

    allocate(active_point(mesh%point_count()), source=.false.)
    allocate(core_point(mesh%point_count()), source=.false.)
    do y_sign = -1, 1
      do x_sign = -1, 1
        call mark_probe( &
          real(x_sign, rk) * global_probe_radius, &
          real(y_sign, rk) * global_probe_radius, active_point, 0)
      end do
    end do
    call mark_probe( &
      0.0_rk, 0.0_rk, core_point, center_probe_stencil_radius)
    call mark_probe(half_core_offset_x, half_core_offset_y, &
      core_point, probe_stencil_radius)
    call mark_probe(-half_core_offset_x, -half_core_offset_y, &
      core_point, probe_stencil_radius)
    allocate(update_point(mesh%point_count()), source=.false.)
    select case (trim(update_region))
    case ("half_core_stencils")
      update_point = core_point
    case ("half_core_disks")
      call mark_disk( &
        half_core_offset_x, half_core_offset_y, update_radius, update_point)
      call mark_disk( &
        -half_core_offset_x, -half_core_offset_y, update_radius, update_point)
    case ("centered_disk")
      call mark_disk(0.0_rk, 0.0_rk, update_radius, update_point)
    case ("centered_ellipse")
      call mark_ellipse( &
        0.0_rk, 0.0_rk, update_radius_x, update_radius_y, update_point)
    end select
    core_point = core_point .or. update_point
    allocate(asymptotic_ray_point(mesh%point_count()), source=.false.)
    allocate(asymptotic_ray_id(mesh%point_count()), source=0)
    call mark_asymptotic_ray_probes()
    active_point = active_point .or. core_point .or. asymptotic_ray_point
    allocate(outer_point, source=active_point .and. .not. core_point)
    if (update_outer_probes) update_point = update_point .or. &
      (outer_point .and. .not. asymptotic_ray_point)
    if (pin_center) update_point( &
      mesh%point_index(nearest_x_node(0.0_rk), nearest_y_node(0.0_rk))) = &
      .false.
    if (.not. any(core_point) .or. .not. any(outer_point)) &
      error stop "double-core benchmark needs both core and outer probes"
    if (.not. any(update_point)) &
      error stop "double-core iteration mask contains no points"
    allocate(checkpoint_point(mesh%point_count()))
    checkpoint_point = update_point
    if (trim(checkpoint_scope) == "full_state") checkpoint_point = .true.
  end subroutine make_probe_masks


  subroutine mark_asymptotic_ray_probes()
    real(rk) :: active_x, active_y

    if (asymptotic_ray_probe_count == 0) return
    select case (trim(update_region))
    case ("centered_ellipse")
      active_x = update_radius_x
      active_y = update_radius_y
    case ("centered_disk")
      active_x = update_radius
      active_y = update_radius
    case default
      error stop "asymptotic ray probes require a centered active domain"
    end select
    if (trim(inactive_halo_policy) == "state_asymptotic" .and. &
        .not. use_elliptical_asymptotic_fit()) then
      active_x = max(active_x, asymptotic_matching_radius)
      active_y = max(active_y, asymptotic_matching_radius)
    end if

    call mark_asymptotic_axis_ray_probes(.true., 1, active_x, 1)
    call mark_asymptotic_axis_ray_probes(.true., -1, active_x, -1)
    call mark_asymptotic_axis_ray_probes(.false., 1, active_y, 2)
    call mark_asymptotic_axis_ray_probes(.false., -1, active_y, -2)
    if (.not. any(asymptotic_ray_point)) &
      error stop "no asymptotic ray probes fit outside the active domain"
  end subroutine mark_asymptotic_ray_probes


  subroutine mark_asymptotic_axis_ray_probes( &
      along_x, sign, active_radius, ray_id)
    logical, intent(in) :: along_x
    integer, intent(in) :: sign, ray_id
    real(rk), intent(in) :: active_radius

    integer, allocatable :: candidate(:)
    integer :: candidate_count, node, point, probe, selected_count
    real(rk) :: coordinate, tolerance

    tolerance = 64.0_rk * epsilon(1.0_rk) * max(1.0_rk, active_radius)
    candidate_count = 0
    if (along_x) then
      do node = 1, mesh%x_cell_count() - 1
        coordinate = mesh%x_coordinate(node)
        if (real(sign, rk) * coordinate > active_radius + tolerance) &
          candidate_count = candidate_count + 1
      end do
    else
      do node = 1, mesh%y_cell_count() - 1
        coordinate = mesh%y_coordinate(node)
        if (real(sign, rk) * coordinate > active_radius + tolerance) &
          candidate_count = candidate_count + 1
      end do
    end if
    if (candidate_count == 0) return

    allocate(candidate(candidate_count))
    candidate_count = 0
    if (along_x) then
      do node = 1, mesh%x_cell_count() - 1
        coordinate = mesh%x_coordinate(node)
        if (real(sign, rk) * coordinate > active_radius + tolerance) then
          candidate_count = candidate_count + 1
          candidate(candidate_count) = node
        end if
      end do
    else
      do node = 1, mesh%y_cell_count() - 1
        coordinate = mesh%y_coordinate(node)
        if (real(sign, rk) * coordinate > active_radius + tolerance) then
          candidate_count = candidate_count + 1
          candidate(candidate_count) = node
        end if
      end do
    end if

    selected_count = min(asymptotic_ray_probe_count, candidate_count)
    do probe = 1, selected_count
      node = candidate(ceiling(real(probe * candidate_count, rk) / &
        real(selected_count + 1, rk)))
      if (along_x) then
        point = mesh%point_index(node, nearest_y_node(0.0_rk))
      else
        point = mesh%point_index(nearest_x_node(0.0_rk), node)
      end if
      asymptotic_ray_point(point) = .true.
      asymptotic_ray_id(point) = ray_id
    end do
  end subroutine mark_asymptotic_axis_ray_probes


  subroutine make_asymptotic_halo_masks()
    real(rk) :: displacement(2), radius, scaled_radius_squared, tolerance
    integer :: point

    allocate(dependent_halo_point(mesh%point_count()), source=.false.)
    allocate(asymptotic_source_point(mesh%point_count()), source=.true.)
    if (trim(inactive_halo_policy) == "frozen") return
    tolerance = 128.0_rk * epsilon(1.0_rk)
    do point = 1, mesh%point_count()
      displacement = mesh%point_coordinate(point) - endpoint%center
      if (use_elliptical_asymptotic_fit()) then
        scaled_radius_squared = &
          (displacement(1) / update_radius_x)**2 + &
          (displacement(2) / update_radius_y)**2
        dependent_halo_point(point) = &
          scaled_radius_squared > 1.0_rk + tolerance
      else
        radius = sqrt(dot_product(displacement, displacement))
        dependent_halo_point(point) = radius > &
          asymptotic_matching_radius + tolerance * &
            max(1.0_rk, asymptotic_matching_radius)
      end if
    end do
    if (any(dependent_halo_point .and. update_point)) &
      error stop "double-core update region crosses the asymptotic halo"
    asymptotic_source_point = .not. dependent_halo_point
    if (.not. any(dependent_halo_point)) &
      error stop "state-asymptotic halo policy selected no dependent points"
  end subroutine make_asymptotic_halo_masks


  subroutine make_asymptotic_probe_iteration_mask()
    real(rk) :: coordinate(2)
    integer :: point

    allocate(asymptotic_probe_iteration_point( &
      mesh%point_count()), source=.false.)
    if (asymptotic_probe_iterations == 0) return
    if (any(asymptotic_ray_point .and. .not. dependent_halo_point)) &
      error stop "asymptotic ray center is not in the extrapolated halo"
    do point = 1, mesh%point_count()
      if (.not. asymptotic_ray_point(point)) cycle
      coordinate = mesh%point_coordinate(point)
      call mark_probe( &
        coordinate(1), coordinate(2), asymptotic_probe_iteration_point, &
        asymptotic_probe_iteration_stencil_radius)
    end do
    asymptotic_probe_iteration_point = &
      asymptotic_probe_iteration_point .and. dependent_halo_point
    if (.not. any(asymptotic_probe_iteration_point)) &
      error stop "shadow asymptotic probe iteration mask is empty"
  end subroutine make_asymptotic_probe_iteration_mask


  subroutine refresh_dependent_halo(field, maximum_change)
    type(spinful_state_2d_t), intent(inout) :: field
    real(rk), intent(out) :: maximum_change

    type(spinful_state_2d_t) :: completed
    integer :: point

    maximum_change = 0.0_rk
    if (trim(inactive_halo_policy) == "state_asymptotic") then
      if (rank == 0) then
        call apply_free_vortex_asymptotic_halo_2d( &
          mesh, field, asymptotic_source_point, endpoint, completed)
        do point = 1, mesh%point_count()
          if (.not. dependent_halo_point(point)) cycle
          maximum_change = max(maximum_change, maxval(abs( &
            completed%order_parameter(:, :, point) - &
            field%order_parameter(:, :, point))))
          maximum_change = max(maximum_change, maxval(abs( &
            completed%current_mean_field(:, point) - &
            field%current_mean_field(:, point))))
        end do
        field = completed
      end if
      call MPI_Bcast(field%order_parameter, size(field%order_parameter), &
        MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
      call require_mpi(ierr == MPI_SUCCESS, &
        "could not broadcast regenerated asymptotic order parameter")
      call MPI_Bcast(field%current_mean_field, &
        size(field%current_mean_field), MPI_DOUBLE_PRECISION, 0, &
        MPI_COMM_WORLD, ierr)
      call require_mpi(ierr == MPI_SUCCESS, &
        "could not broadcast regenerated asymptotic mean field")
      call MPI_Bcast(maximum_change, 1, MPI_DOUBLE_PRECISION, 0, &
        MPI_COMM_WORLD, ierr)
      call require_mpi(ierr == MPI_SUCCESS, &
        "could not broadcast asymptotic halo diagnostic")
    end if
  end subroutine refresh_dependent_halo


  subroutine configure_zero_mode_basis()
    logical :: include_gauge, include_orientation, include_translations

    include_gauge = .false.
    include_translations = .false.
    include_orientation = .false.
    select case (trim(zero_mode_projection))
    case ("none")
      use_zero_mode_projection = .false.
      return
    case ("translations")
      include_translations = .true.
    case ("translation_orientation")
      include_translations = .true.
      include_orientation = .true.
    case ("gauge_translation_orientation")
      include_gauge = .true.
      include_translations = .true.
      include_orientation = .true.
    end select
    use_zero_mode_projection = .true.
    if (rank == 0) call build_zero_mode_tangent_basis_2d( &
      mesh, state, update_point, include_gauge, include_translations, &
      include_orientation, zero_mode_basis)
  end subroutine configure_zero_mode_basis


  subroutine mark_disk(center_x, center_y, radius, mask)
    real(rk), intent(in) :: center_x, center_y, radius
    logical, intent(inout) :: mask(:)

    real(rk) :: coordinate(2), displacement(2)
    integer :: point

    do point = 1, mesh%point_count()
      coordinate = mesh%point_coordinate(point)
      displacement = coordinate - [center_x, center_y]
      if (dot_product(displacement, displacement) <= &
          radius**2 * (1.0_rk + 64.0_rk * epsilon(1.0_rk))) &
        mask(point) = .true.
    end do
  end subroutine mark_disk


  subroutine mark_ellipse(center_x, center_y, radius_x, radius_y, mask)
    real(rk), intent(in) :: center_x, center_y, radius_x, radius_y
    logical, intent(inout) :: mask(:)

    real(rk) :: coordinate(2), scaled_radius_squared
    integer :: point

    do point = 1, mesh%point_count()
      coordinate = mesh%point_coordinate(point)
      scaled_radius_squared = &
        ((coordinate(1) - center_x) / radius_x)**2 + &
        ((coordinate(2) - center_y) / radius_y)**2
      if (scaled_radius_squared <= &
          1.0_rk + 64.0_rk * epsilon(1.0_rk)) &
        mask(point) = .true.
    end do
  end subroutine mark_ellipse


  subroutine mark_probe(x, y, mask, stencil_radius)
    real(rk), intent(in) :: x, y
    logical, intent(inout) :: mask(:)
    integer, intent(in) :: stencil_radius

    integer :: center_x, center_y, point, x_node, y_node

    center_x = nearest_x_node(x)
    center_y = nearest_y_node(y)
    do y_node = max(0, center_y - stencil_radius), &
                min(mesh%y_cell_count(), center_y + stencil_radius)
      do x_node = max(0, center_x - stencil_radius), &
                  min(mesh%x_cell_count(), center_x + stencil_radius)
        point = mesh%point_index(x_node, y_node)
        mask(point) = .true.
      end do
    end do
  end subroutine mark_probe


  integer function nearest_x_node(x) result(nearest)
    real(rk), intent(in) :: x

    integer :: node
    real(rk) :: difference, minimum_difference

    nearest = 0
    minimum_difference = huge(1.0_rk)
    do node = 0, mesh%x_cell_count()
      difference = abs(mesh%x_coordinate(node) - x)
      if (difference < minimum_difference) then
        nearest = node
        minimum_difference = difference
      end if
    end do
  end function nearest_x_node


  integer function nearest_y_node(y) result(nearest)
    real(rk), intent(in) :: y

    integer :: node
    real(rk) :: difference, minimum_difference

    nearest = 0
    minimum_difference = huge(1.0_rk)
    do node = 0, mesh%y_cell_count()
      difference = abs(mesh%y_coordinate(node) - y)
      if (difference < minimum_difference) then
        nearest = node
        minimum_difference = difference
      end if
    end do
  end function nearest_y_node


  subroutine refresh_endpoint_reference()
    if (trim(endpoint_policy) /= "state_asymptotic") return
    if (use_elliptical_asymptotic_fit()) then
      call fit_state_bulk_reference( &
        endpoint, mesh, state, bulk_reference_fit_radius, &
        bulk_reference_sample_count, bulk_reference_fit_rms, &
        bulk_reference_fit_maximum, fitting_radii=asymptotic_matching_radii, &
        update_reference=.false.)
    else
      call fit_state_bulk_reference( &
        endpoint, mesh, state, bulk_reference_fit_radius, &
        bulk_reference_sample_count, bulk_reference_fit_rms, &
        bulk_reference_fit_maximum, update_reference=.false.)
    end if
  end subroutine refresh_endpoint_reference


  pure logical function use_elliptical_asymptotic_fit() result(enabled)
    enabled = trim(update_region) == "centered_ellipse" .and. &
      all(asymptotic_fit_inner_radii > 0.0_rk) .and. &
      all(asymptotic_matching_radii > 0.0_rk)
  end function use_elliptical_asymptotic_fit


  subroutine compute_current_residuals()
    call compute_he3_field_residual( &
      evaluated_state, mapped, residual, active_point_mask=active_point)
    call compute_he3_field_residual( &
      evaluated_state, mapped, core_residual, active_point_mask=core_point)
    call compute_he3_field_residual( &
      evaluated_state, mapped, outer_residual, active_point_mask=outer_point)
  end subroutine compute_current_residuals


  subroutine iterate_asymptotic_probe_shadow()
    character(len=2048) :: message
    integer :: point, shadow_iteration, status

    shadow_probe_baseline = state
    shadow_probe_state = state
    shadow_probe_terminal_status = "iteration_limit"
    shadow_probe_anderson_report = anderson_report_t()
    shadow_probe_residual = he3_field_residual_report_t()
    if (passive_probe_output_requested) then
      passive_probe_source = "shadow_first_map"
      passive_probe_map_succeeded = .false.
    end if
    if (rank == 0) then
      call shadow_probe_accelerator%initialize( &
        new_src_masked_iteration_vector_size_2d( &
          shadow_probe_state, asymptotic_probe_iteration_point), &
        asymptotic_probe_anderson_history_limit, &
        anderson_progress_threshold, &
        asymptotic_probe_anderson_maximum_mixing)
      open(newunit=shadow_probe_history_unit, &
           file=trim(asymptotic_probe_history_file), status="replace", &
           action="write", iostat=status, iomsg=message)
      if (status /= 0) error stop &
        "cannot open shadow probe history: " // trim(message)
      write(shadow_probe_history_unit, '(a)') &
        "# map update rms maximum relative_l2 mixing history_size " // &
        "map_seconds normalization_error"
    end if

    do shadow_iteration = 1, asymptotic_probe_iterations
      call evaluate_shadow_probe_map()
      if (.not. shadow_probe_map_succeeded) then
        shadow_probe_terminal_status = "map_failure"
        exit
      end if
      if (shadow_iteration == 1) then
        shadow_probe_first_mapped = shadow_probe_mapped
        if (passive_probe_output_requested) then
          passive_probe_baseline = shadow_probe_baseline
          passive_probe_mapped = shadow_probe_first_mapped
          passive_probe_diagnostics = shadow_probe_diagnostics
          passive_probe_map_succeeded = .true.
          passive_probe_map_evaluations = 1
          passive_probe_map_seconds = &
            passive_probe_diagnostics%maximum_rank_elapsed_seconds
          passive_probe_maximum_normalization_error = &
            passive_probe_diagnostics%global%maximum_normalization_error
        end if
      end if
      call update_mpi_he3_state_with_anderson( &
        MPI_COMM_WORLD, mesh, shadow_probe_accelerator, &
        shadow_probe_state, shadow_probe_mapped, &
        asymptotic_probe_convergence_tolerance, shadow_probe_next, &
        shadow_probe_anderson_report, shadow_probe_residual, &
        asymptotic_probe_iteration_point)
      shadow_probe_iterations_completed = shadow_iteration
      if (rank == 0) call write_shadow_probe_history_row(shadow_iteration)
      shadow_probe_state = shadow_probe_next
      if (shadow_probe_anderson_report%converged) then
        shadow_probe_terminal_status = "converged"
        exit
      end if
    end do

    if (shadow_probe_map_succeeded) then
      call evaluate_shadow_probe_map()
      if (shadow_probe_map_succeeded) then
        call compute_he3_field_residual( &
          shadow_probe_state, shadow_probe_mapped, shadow_probe_residual, &
          active_point_mask=asymptotic_probe_iteration_point)
        if (shadow_probe_residual%maximum_absolute_residual <= &
            asymptotic_probe_convergence_tolerance) &
          shadow_probe_terminal_status = "converged"
        if (rank == 0) call write_shadow_probe_history_row( &
          shadow_probe_iterations_completed)
      else
        shadow_probe_terminal_status = "map_failure"
      end if
    end if

    do point = 1, mesh%point_count()
      if (asymptotic_probe_iteration_point(point)) cycle
      if (maxval(abs(shadow_probe_state%order_parameter(:, :, point) - &
                     shadow_probe_baseline%order_parameter(:, :, point))) > &
            0.0_rk .or. &
          maxval(abs(shadow_probe_state%current_mean_field(:, point) - &
                     shadow_probe_baseline%current_mean_field(:, point))) > &
            0.0_rk) &
        error stop "shadow probe iteration modified its frozen background"
    end do
    if (rank == 0) close(shadow_probe_history_unit)
  end subroutine iterate_asymptotic_probe_shadow


  subroutine evaluate_final_passive_probe_map()
    passive_probe_source = "fresh_final_map"
    passive_probe_baseline = state
    passive_probe_map_succeeded = .false.
    if (trim(endpoint_policy) == "state_asymptotic") then
      call evaluate_mpi_he3_field_map( &
        MPI_COMM_WORLD, mesh, passive_probe_baseline, angular, energy, &
        feedback_scale, trajectory_maximum_step, minimum_substep_count, &
        boundary_relaxation_distance, passive_probe_mapped, &
        passive_probe_diagnostics, passive_probe_map_succeeded, &
        free_vortex_endpoint=endpoint, &
        active_point_mask=asymptotic_ray_point)
    else
      call evaluate_mpi_he3_field_map( &
        MPI_COMM_WORLD, mesh, passive_probe_baseline, angular, energy, &
        feedback_scale, trajectory_maximum_step, minimum_substep_count, &
        boundary_relaxation_distance, passive_probe_mapped, &
        passive_probe_diagnostics, passive_probe_map_succeeded, &
        active_point_mask=asymptotic_ray_point)
    end if
    passive_probe_map_evaluations = 1
    passive_probe_map_seconds = &
      passive_probe_diagnostics%maximum_rank_elapsed_seconds
    passive_probe_maximum_normalization_error = &
      passive_probe_diagnostics%global%maximum_normalization_error
  end subroutine evaluate_final_passive_probe_map


  subroutine evaluate_shadow_probe_map()
    call evaluate_mpi_he3_field_map( &
      MPI_COMM_WORLD, mesh, shadow_probe_state, angular, energy, &
      feedback_scale, trajectory_maximum_step, minimum_substep_count, &
      boundary_relaxation_distance, shadow_probe_mapped, &
      shadow_probe_diagnostics, shadow_probe_map_succeeded, &
      free_vortex_endpoint=endpoint, &
      active_point_mask=asymptotic_probe_iteration_point)
    shadow_probe_map_evaluations = shadow_probe_map_evaluations + 1
    shadow_probe_cumulative_map_seconds = &
      shadow_probe_cumulative_map_seconds + &
      shadow_probe_diagnostics%maximum_rank_elapsed_seconds
    shadow_probe_maximum_normalization_error = max( &
      shadow_probe_maximum_normalization_error, &
      shadow_probe_diagnostics%global%maximum_normalization_error)
  end subroutine evaluate_shadow_probe_map


  subroutine write_shadow_probe_history_row(update)
    integer, intent(in) :: update

    write(shadow_probe_history_unit, &
      '(2(i8,1x),4(es24.16e3,1x),i8,1x,2(es24.16e3,1x))') &
      shadow_probe_map_evaluations, update, &
      shadow_probe_residual%rms_residual, &
      shadow_probe_residual%maximum_absolute_residual, &
      shadow_probe_residual%relative_l2_residual, &
      shadow_probe_anderson_report%mixing, &
      shadow_probe_anderson_report%history_size, &
      shadow_probe_diagnostics%maximum_rank_elapsed_seconds, &
      shadow_probe_diagnostics%global%maximum_normalization_error
    flush(shadow_probe_history_unit)
  end subroutine write_shadow_probe_history_row


  subroutine measure_current_core_shape()
    real(rk) :: amplitude, coordinate(2), measurement_radius, tolerance
    integer :: center_point, point
    logical :: found_negative, found_positive

    tolerance = 128.0_rk * epsilon(1.0_rk) * &
      max(1.0_rk, mesh%maximum_spacing())
    measured_positive_half_core_amplitude = huge(1.0_rk)
    measured_negative_half_core_amplitude = huge(1.0_rk)
    measured_positive_half_core_y = 0.0_rk
    measured_negative_half_core_y = 0.0_rk
    found_positive = .false.
    found_negative = .false.
    measurement_radius = update_radius
    if (trim(update_region) == "centered_ellipse") &
      measurement_radius = update_radius_y
    do point = 1, mesh%point_count()
      coordinate = mesh%point_coordinate(point)
      if (abs(coordinate(1)) > tolerance) cycle
      if (abs(coordinate(2)) > measurement_radius + tolerance) cycle
      amplitude = sqrt(sum(abs(state%order_parameter(:, :, point))**2) / &
        3.0_rk)
      if (coordinate(2) > tolerance .and. &
          amplitude < measured_positive_half_core_amplitude) then
        measured_positive_half_core_amplitude = amplitude
        measured_positive_half_core_y = coordinate(2)
        found_positive = .true.
      else if (coordinate(2) < -tolerance .and. &
               amplitude < measured_negative_half_core_amplitude) then
        measured_negative_half_core_amplitude = amplitude
        measured_negative_half_core_y = coordinate(2)
        found_negative = .true.
      end if
    end do
    if (.not. found_positive .or. .not. found_negative) &
      error stop "cannot locate both double-core minima on the y axis"
    measured_half_core_separation = measured_positive_half_core_y - &
      measured_negative_half_core_y
    center_point = mesh%point_index(nearest_x_node(0.0_rk), &
                                    nearest_y_node(0.0_rk))
    measured_center_pair_amplitude = sqrt(sum(abs( &
      state%order_parameter(:, :, center_point))**2) / 3.0_rk)
  end subroutine measure_current_core_shape


  subroutine open_iteration_history()
    character(len=2048) :: message
    integer :: status

    open(newunit=history_unit, file=trim(history_file), status="replace", &
         action="write", iostat=status, iomsg=message)
    if (status /= 0) &
      error stop "cannot open double-core iteration history: " // &
        trim(message)
    write(history_unit, '(a)') &
      "# iteration sample_rms sample_max sample_relative_l2 " // &
      "update_rms update_max update_relative_l2 anderson_norm p history " // &
      "normalization_error map_seconds zero_mode_removed " // &
      "zero_mode_fraction zero_mode_gauge zero_mode_translation_x " // &
      "zero_mode_translation_y zero_mode_orientation " // &
      "half_core_y_minus half_core_y_plus half_core_separation " // &
      "center_pair_amplitude half_core_amplitude_minus " // &
      "half_core_amplitude_plus"
  end subroutine open_iteration_history


  subroutine write_single_map_history(iteration_number)
    integer, intent(in) :: iteration_number

    write(history_unit, &
      '(i8,1x,8(es24.16e3,1x),i8,1x,14(es24.16e3,1x))') &
      iteration_number, residual%rms_residual, &
      residual%maximum_absolute_residual, residual%relative_l2_residual, &
      0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0, &
      diagnostics%global%maximum_normalization_error, &
      diagnostics%maximum_rank_elapsed_seconds, &
      0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, &
      measured_negative_half_core_y, measured_positive_half_core_y, &
      measured_half_core_separation, measured_center_pair_amplitude, &
      measured_negative_half_core_amplitude, &
      measured_positive_half_core_amplitude
    flush(history_unit)
  end subroutine write_single_map_history


  subroutine write_iteration_history(iteration_number)
    integer, intent(in) :: iteration_number

    write(history_unit, &
      '(i8,1x,8(es24.16e3,1x),i8,1x,14(es24.16e3,1x))') &
      iteration_number, residual%rms_residual, &
      residual%maximum_absolute_residual, residual%relative_l2_residual, &
      update_residual%rms_residual, &
      update_residual%maximum_absolute_residual, &
      update_residual%relative_l2_residual, &
      anderson_report%residual_norm, anderson_report%mixing, &
      anderson_report%history_size, &
      diagnostics%global%maximum_normalization_error, &
      diagnostics%maximum_rank_elapsed_seconds, &
      zero_mode_report%removed_residual_norm, &
      zero_mode_report%removed_fraction, zero_mode_report%coefficient, &
      measured_negative_half_core_y, measured_positive_half_core_y, &
      measured_half_core_separation, measured_center_pair_amplitude, &
      measured_negative_half_core_amplitude, &
      measured_positive_half_core_amplitude
    flush(history_unit)
  end subroutine write_iteration_history


  subroutine print_iteration(iteration_number)
    integer, intent(in) :: iteration_number

    print '(a,i4,a,2(es11.3,1x),a,f6.2,a,i3,a,f8.2,a,f7.2)', &
      "solved (double-core AA): ", iteration_number, &
      " update error(avg,max)= (", &
      update_residual%mean_absolute_residual, &
      update_residual%maximum_absolute_residual, ") p=", &
      anderson_report%mixing, " md=", anderson_report%history_size, &
      " map[s]=", diagnostics%maximum_rank_elapsed_seconds, &
      " a=", measured_half_core_separation
    if (use_zero_mode_projection) &
      print '(a,es11.3,a,es11.3)', "  zero-mode removed norm=", &
        zero_mode_report%removed_residual_norm, " fraction=", &
        zero_mode_report%removed_fraction
  end subroutine print_iteration


  subroutine write_metrics()
    character(len=2048) :: message
    integer :: status, unit
    real(rk) :: maximum_coordinate(2)

    maximum_coordinate = mesh%point_coordinate(residual%maximum_point)
    open(newunit=unit, file=trim(metrics_file), status="replace", &
         action="write", iostat=status, iomsg=message)
    if (status /= 0) &
      error stop "cannot open double-core metrics: " // trim(message)
    if (trim(initialization_mode) == "regularized_london" .or. &
        trim(initialization_mode) == "historical_nop" .or. &
        trim(initialization_mode) == "historical_aop" .or. &
        trim(initialization_mode) == "historical_dop") then
      if (solve_requested) then
        write(unit, '(a)') "kind=double_core_from_scratch_iteration"
      else
        write(unit, '(a)') "kind=double_core_from_scratch_one_map"
      end if
    else if (solve_requested) then
      write(unit, '(a)') "kind=legacy_double_core_sparse_iteration"
    else
      write(unit, '(a)') "kind=legacy_double_core_sparse_one_map"
    end if
    write(unit, '(a,a)') "terminal_status=", trim(terminal_status)
    write(unit, '(a,a)') "initialization_mode=", trim(initialization_mode)
    write(unit, '(a,a)') "scratch_mesh_kind=", trim(scratch_mesh_kind)
    write(unit, '(a,es24.16e3)') "scratch_half_width=", &
      scratch_half_width
    write(unit, '(a,es24.16e3)') "seed_half_core_offset=", &
      abs(half_core_offset_y)
    write(unit, '(a,es24.16e3)') "seed_core_width=", seed_core_width
    write(unit, '(a,es24.16e3)') "seed_domain_wall_width=", &
      seed_domain_wall_width
    write(unit, '(a,es24.16e3)') "seed_reduced_temperature=", &
      seed_reduced_temperature
    if (trim(initialization_mode) == "regularized_london") then
      write(unit, '(a,es24.16e3)') "seed_center_pair_amplitude=", &
        seed_report%center_pair_amplitude
      write(unit, '(a,es24.16e3)') &
        "seed_positive_half_core_pair_amplitude=", &
        seed_report%positive_half_core_pair_amplitude
      write(unit, '(a,es24.16e3)') &
        "seed_negative_half_core_pair_amplitude=", &
        seed_report%negative_half_core_pair_amplitude
      write(unit, '(a,es24.16e3)') "seed_outer_pair_amplitude=", &
        seed_report%outer_pair_amplitude
    else if (trim(initialization_mode) == "historical_nop" .or. &
             trim(initialization_mode) == "historical_aop" .or. &
             trim(initialization_mode) == "historical_dop") then
      write(unit, '(a,a)') "historical_seed_kind=", &
        historical_seed_report%kind
      write(unit, '(a,es24.16e3)') "historical_seed_core_length_scale=", &
        historical_seed_report%core_length_scale
      write(unit, '(a,es24.16e3)') "seed_center_pair_amplitude=", &
        historical_seed_report%center_pair_amplitude
      write(unit, '(a,es24.16e3)') "seed_outer_pair_amplitude=", &
        historical_seed_report%outer_pair_amplitude
    end if
    write(unit, '(a,i0)') "sparse_restart_applied_points=", &
      sparse_restart_applied_points
    write(unit, '(a,a)') "endpoint_policy=", trim(endpoint_policy)
    write(unit, '(a,a)') "inactive_halo_policy=", &
      trim(inactive_halo_policy)
    write(unit, '(a,i0)') "dependent_asymptotic_halo_points=", &
      count(dependent_halo_point)
    write(unit, '(a,i0)') "frozen_nonupdated_points=", count( &
      .not. update_point .and. .not. dependent_halo_point)
    write(unit, '(a,es24.16e3)') "initial_halo_maximum_change=", &
      initial_halo_maximum_change
    write(unit, '(a,es24.16e3)') "last_halo_maximum_change=", &
      last_halo_maximum_change
    write(unit, '(a,i0)') "mpi_ranks=", rank_count
    write(unit, '(a,i0)') "mesh_x_points=", mesh%x_point_count()
    write(unit, '(a,i0)') "mesh_y_points=", mesh%y_point_count()
    write(unit, '(a,i0)') "sampled_points=", count(active_point)
    write(unit, '(a,i0)') "core_points=", count(core_point)
    write(unit, '(a,i0)') "outer_points=", count(outer_point)
    write(unit, '(a,i0)') "asymptotic_ray_probe_count_requested=", &
      asymptotic_ray_probe_count
    write(unit, '(a,i0)') "asymptotic_ray_probe_points=", &
      count(asymptotic_ray_point)
    write(unit, '(a,a)') "asymptotic_probe_file=", &
      trim(asymptotic_probe_file)
    write(unit, '(a,a)') "asymptotic_probe_output_source=", &
      trim(passive_probe_source)
    write(unit, '(a,l1)') "asymptotic_probe_output_map_succeeded=", &
      passive_probe_map_succeeded
    write(unit, '(a,i0)') "asymptotic_probe_output_map_evaluations=", &
      passive_probe_map_evaluations
    write(unit, '(a,es24.16e3)') &
      "asymptotic_probe_output_map_seconds=", passive_probe_map_seconds
    write(unit, '(a,es24.16e3)') &
      "asymptotic_probe_output_maximum_normalization_error=", &
      passive_probe_maximum_normalization_error
    write(unit, '(a,i0)') "asymptotic_probe_iterations_requested=", &
      asymptotic_probe_iterations
    write(unit, '(a,i0)') "asymptotic_probe_iteration_stencil_radius=", &
      asymptotic_probe_iteration_stencil_radius
    write(unit, '(a,i0)') "asymptotic_probe_iteration_points=", &
      count(asymptotic_probe_iteration_point)
    write(unit, '(a,es24.16e3)') &
      "asymptotic_probe_convergence_tolerance=", &
      asymptotic_probe_convergence_tolerance
    write(unit, '(a,a)') "asymptotic_probe_history_file=", &
      trim(asymptotic_probe_history_file)
    write(unit, '(a,a)') "asymptotic_shadow_file=", &
      trim(asymptotic_shadow_file)
    write(unit, '(a,a)') "asymptotic_probe_terminal_status=", &
      trim(shadow_probe_terminal_status)
    write(unit, '(a,i0)') "asymptotic_probe_iterations_completed=", &
      shadow_probe_iterations_completed
    write(unit, '(a,i0)') "asymptotic_probe_map_evaluations=", &
      shadow_probe_map_evaluations
    write(unit, '(a,es24.16e3)') &
      "asymptotic_probe_terminal_rms=", shadow_probe_residual%rms_residual
    write(unit, '(a,es24.16e3)') &
      "asymptotic_probe_terminal_maximum=", &
      shadow_probe_residual%maximum_absolute_residual
    write(unit, '(a,es24.16e3)') &
      "asymptotic_probe_terminal_relative_l2=", &
      shadow_probe_residual%relative_l2_residual
    write(unit, '(a,es24.16e3)') &
      "asymptotic_probe_maximum_normalization_error=", &
      shadow_probe_maximum_normalization_error
    write(unit, '(a,es24.16e3)') &
      "asymptotic_probe_cumulative_map_seconds=", &
      shadow_probe_cumulative_map_seconds
    write(unit, '(a,i0)') "updated_points=", count(update_point)
    write(unit, '(a,a)') "checkpoint_scope=", trim(checkpoint_scope)
    write(unit, '(a,i0)') "checkpoint_points=", count(checkpoint_point)
    write(unit, '(a,a)') "update_region=", trim(update_region)
    write(unit, '(a,es24.16e3)') "update_radius=", update_radius
    write(unit, '(a,es24.16e3)') "update_radius_x=", update_radius_x
    write(unit, '(a,es24.16e3)') "update_radius_y=", update_radius_y
    write(unit, '(a,l1)') "elliptical_asymptotic_fit=", &
      use_elliptical_asymptotic_fit()
    write(unit, '(a,es24.16e3)') "asymptotic_fit_inner_radius=", &
      asymptotic_fit_inner_radius
    write(unit, '(a,es24.16e3)') "asymptotic_matching_radius=", &
      asymptotic_matching_radius
    write(unit, '(a,2(es24.16e3,1x))') &
      "asymptotic_fit_inner_radii=", asymptotic_fit_inner_radii
    write(unit, '(a,2(es24.16e3,1x))') &
      "asymptotic_matching_radii=", asymptotic_matching_radii
    write(unit, '(a,es24.16e3)') "measured_half_core_y_minus=", &
      measured_negative_half_core_y
    write(unit, '(a,es24.16e3)') "measured_half_core_y_plus=", &
      measured_positive_half_core_y
    write(unit, '(a,es24.16e3)') "measured_half_core_separation=", &
      measured_half_core_separation
    write(unit, '(a,es24.16e3)') "measured_center_pair_amplitude=", &
      measured_center_pair_amplitude
    write(unit, '(a,es24.16e3)') &
      "measured_negative_half_core_amplitude=", &
      measured_negative_half_core_amplitude
    write(unit, '(a,es24.16e3)') &
      "measured_positive_half_core_amplitude=", &
      measured_positive_half_core_amplitude
    write(unit, '(a,l1)') "center_pinned=", pin_center
    write(unit, '(a,a)') "zero_mode_projection=", &
      trim(zero_mode_projection)
    if (use_zero_mode_projection) then
      write(unit, '(a)') &
        "orientation_constraint=inactive_interior_plus_template_tangents"
      write(unit, '(a)') &
        "zero_mode_metric=unweighted_active_Anderson_vector_Euclidean"
      write(unit, '(a)') &
        "zero_mode_orientation_action=gauge_compensated_coordinate_rotation"
      write(unit, '(a,i0)') "zero_mode_requested_count=", &
        zero_mode_basis%requested_count
      write(unit, '(a,i0)') "zero_mode_retained_count=", &
        zero_mode_basis%retained_count
      write(unit, '(a,4(es24.16e3,1x))') "zero_mode_raw_norms=", &
        zero_mode_basis%raw_norm
      write(unit, '(a,es24.16e3)') "zero_mode_original_residual_norm=", &
        zero_mode_report%original_residual_norm
      write(unit, '(a,es24.16e3)') "zero_mode_projected_residual_norm=", &
        zero_mode_report%projected_residual_norm
      write(unit, '(a,es24.16e3)') "zero_mode_removed_residual_norm=", &
        zero_mode_report%removed_residual_norm
      write(unit, '(a,es24.16e3)') "zero_mode_removed_fraction=", &
        zero_mode_report%removed_fraction
      write(unit, '(a,4(es24.16e3,1x))') "zero_mode_coefficients=", &
        zero_mode_report%coefficient
    else
      write(unit, '(a)') "orientation_constraint=inactive_interior"
    end if
    write(unit, '(a,i0)') "anderson_vector_size=", &
      new_src_masked_iteration_vector_size_2d(state, update_point)
    write(unit, '(a,i0)') "full_vector_size=", 21 * mesh%point_count()
    write(unit, '(a,i0)') "maximum_iterations=", maximum_iterations
    write(unit, '(a,i0)') "iterations_completed=", iterations_completed
    write(unit, '(a,i0)') "map_evaluations=", map_evaluations
    write(unit, '(a,i0)') "checkpoint_interval=", checkpoint_interval
    write(unit, '(a,a)') "checkpoint_file=", trim(checkpoint_file)
    write(unit, '(a,i0)') "directions=", angular%direction_count()
    write(unit, '(a,i0)') "poles=", energy%pole_count()
    write(unit, '(a,es24.16e3)') "import_gap_norm_difference=", &
      import_report%maximum_gap_norm_difference
    write(unit, '(a,es24.16e3)') "import_in_plane_difference=", &
      import_report%maximum_in_plane_magnitude_difference
    if (trim(endpoint_policy) == "state_asymptotic") then
      write(unit, '(a,a)') "bulk_reference_policy=", &
        trim(bulk_reference_policy)
      write(unit, '(a,es24.16e3)') "bulk_reference_phase_offset=", &
        endpoint%bulk_phase_offset
      write(unit, '(a,es24.16e3)') &
        "bulk_reference_annular_rms_deviation=", &
        bulk_reference_fit_rms
      write(unit, '(a,es24.16e3)') &
        "bulk_reference_annular_maximum_deviation=", &
        bulk_reference_fit_maximum
      write(unit, '(a,3(es24.16e3,1x))') "bulk_rotation_row_1=", &
        endpoint%bulk_rotation(1, :)
      write(unit, '(a,3(es24.16e3,1x))') "bulk_rotation_row_2=", &
        endpoint%bulk_rotation(2, :)
      write(unit, '(a,3(es24.16e3,1x))') "bulk_rotation_row_3=", &
        endpoint%bulk_rotation(3, :)
    end if
    call write_residual_metrics(unit, "sample", residual)
    call write_residual_metrics(unit, "core", core_residual)
    call write_residual_metrics(unit, "outer", outer_residual)
    if (solve_requested) then
      write(unit, '(a,es24.16e3)') "final_anderson_mixing=", &
        anderson_report%mixing
      write(unit, '(a,i0)') "final_anderson_history=", &
        anderson_report%history_size
    end if
    write(unit, '(a,es24.16e3)') "maximum_normalization_error=", &
      maximum_normalization_error_over_solve
    write(unit, '(a,es24.16e3)') "maximum_rank_seconds=", &
      diagnostics%maximum_rank_elapsed_seconds
    write(unit, '(a,es24.16e3)') "cumulative_map_seconds=", &
      cumulative_map_seconds
    write(unit, '(a,2(es24.16e3,1x))') &
      "maximum_residual_coordinate=", maximum_coordinate
    close(unit)
  end subroutine write_metrics


  subroutine write_residual_metrics(unit, prefix, report)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: prefix
    type(he3_field_residual_report_t), intent(in) :: report

    write(unit, '(a,es24.16e3)') trim(prefix) // "_mean_absolute=", &
      report%mean_absolute_residual
    write(unit, '(a,es24.16e3)') trim(prefix) // "_rms=", &
      report%rms_residual
    write(unit, '(a,es24.16e3)') trim(prefix) // "_maximum_absolute=", &
      report%maximum_absolute_residual
    write(unit, '(a,es24.16e3)') trim(prefix) // "_relative_l2=", &
      report%relative_l2_residual
    write(unit, '(a,es24.16e3)') trim(prefix) // "_maximum_point_rms=", &
      report%maximum_point_rms_residual
  end subroutine write_residual_metrics


  subroutine write_sampled_residuals()
    character(len=2048) :: message
    complex(rk) :: difference(3, 3)
    integer :: orbital, point, spin, status, unit
    real(rk) :: coordinate(2), current_gap, current_mean_norm
    real(rk) :: mapped_gap, mapped_mean_norm, maximum_absolute, point_rms
    real(rk) :: sum_squared, value

    open(newunit=unit, file=trim(sample_file), status="replace", &
         action="write", iostat=status, iomsg=message)
    if (status /= 0) &
      error stop "cannot open double-core sampled residuals: " // trim(message)
    write(unit, '(a)') &
      "# x y region input_gap mapped_gap gap_difference " // &
      "input_mean_field_norm mapped_mean_field_norm point_rms maximum_absolute"
    do point = 1, mesh%point_count()
      if (.not. active_point(point)) cycle
      coordinate = mesh%point_coordinate(point)
      current_gap = sqrt(sum(abs( &
        evaluated_state%order_parameter(:, :, point))**2) / &
                         3.0_rk)
      mapped_gap = sqrt(sum(abs(mapped%order_parameter(:, :, point))**2) / &
                        3.0_rk)
      current_mean_norm = sqrt(sum( &
        evaluated_state%current_mean_field(:, point)**2))
      mapped_mean_norm = sqrt(sum(mapped%current_mean_field(:, point)**2))
      difference = mapped%order_parameter(:, :, point) - &
        evaluated_state%order_parameter(:, :, point)
      sum_squared = sum(abs(difference)**2) + sum(( &
        mapped%current_mean_field(:, point) - &
        evaluated_state%current_mean_field(:, point))**2)
      point_rms = sqrt(sum_squared / 21.0_rk)
      maximum_absolute = maxval(abs( &
        mapped%current_mean_field(:, point) - &
        evaluated_state%current_mean_field(:, point)))
      do spin = 1, 3
        do orbital = 1, 3
          value = abs(real(difference(spin, orbital), rk))
          maximum_absolute = max(maximum_absolute, value)
          value = abs(aimag(difference(spin, orbital)))
          maximum_absolute = max(maximum_absolute, value)
        end do
      end do
      write(unit, '(2(es24.16e3,1x),i1,1x,7(es24.16e3,1x))') &
        coordinate, merge(1, 0, core_point(point)), current_gap, mapped_gap, &
        mapped_gap - current_gap, current_mean_norm, mapped_mean_norm, &
        point_rms, maximum_absolute
    end do
    close(unit)
  end subroutine write_sampled_residuals


  subroutine write_asymptotic_probe_states()
    character(len=2048) :: message
    integer :: component, column, orbital, point, spin, status, unit
    real(rk) :: coordinate(2), row(45)

    open(newunit=unit, file=trim(asymptotic_probe_file), status="replace", &
         action="write", iostat=status, iomsg=message)
    if (status /= 0) &
      error stop "cannot open asymptotic ray probe output: " // trim(message)
    write(unit, '(a)', advance='no') "# x y ray_id"
    do spin = 1, 3
      do orbital = 1, 3
        write(unit, '(a,i0,i0,a)', advance='no') &
          " input_A", spin, orbital, "_re"
        write(unit, '(a,i0,i0,a)', advance='no') &
          " input_A", spin, orbital, "_im"
      end do
    end do
    do spin = 1, 3
      do orbital = 1, 3
        write(unit, '(a,i0,i0,a)', advance='no') &
          " mapped_A", spin, orbital, "_re"
        write(unit, '(a,i0,i0,a)', advance='no') &
          " mapped_A", spin, orbital, "_im"
      end do
    end do
    do component = 1, 3
      write(unit, '(a,i0)', advance='no') " input_mf", component
    end do
    do component = 1, 3
      write(unit, '(a,i0)', advance='no') " mapped_mf", component
    end do
    write(unit, *)

    do point = 1, mesh%point_count()
      if (.not. asymptotic_ray_point(point)) cycle
      coordinate = mesh%point_coordinate(point)
      row = 0.0_rk
      row(1:3) = [coordinate, real(asymptotic_ray_id(point), rk)]
      column = 4
      do spin = 1, 3
        do orbital = 1, 3
          row(column) = real( &
            passive_probe_baseline%order_parameter(spin, orbital, point), rk)
          row(column + 1) = aimag( &
            passive_probe_baseline%order_parameter(spin, orbital, point))
          column = column + 2
        end do
      end do
      do spin = 1, 3
        do orbital = 1, 3
          row(column) = real( &
            passive_probe_mapped%order_parameter(spin, orbital, point), rk)
          row(column + 1) = aimag( &
            passive_probe_mapped%order_parameter(spin, orbital, point))
          column = column + 2
        end do
      end do
      row(column:column + 2) = &
        passive_probe_baseline%current_mean_field(:, point)
      column = column + 3
      row(column:column + 2) = &
        passive_probe_mapped%current_mean_field(:, point)
      write(unit, '(*(es24.16e3,1x))') row
    end do
    close(unit)
  end subroutine write_asymptotic_probe_states


  subroutine write_asymptotic_shadow_states()
    character(len=2048) :: message
    character(len=16), parameter :: label(4) = [ &
      character(len=16) :: "baseline", "first_map", "shadow", "final_map"]
    integer :: column, component, field_index, orbital, point
    integer :: spin, status, unit
    real(rk) :: coordinate(2), row(87)

    open(newunit=unit, file=trim(asymptotic_shadow_file), status="replace", &
         action="write", iostat=status, iomsg=message)
    if (status /= 0) &
      error stop "cannot open asymptotic shadow output: " // trim(message)
    write(unit, '(a)', advance='no') "# x y ray_id"
    do field_index = 1, 4
      do spin = 1, 3
        do orbital = 1, 3
          write(unit, '(3a,i0,i0,a)', advance='no') &
            " ", trim(label(field_index)), "_A", spin, orbital, "_re"
          write(unit, '(3a,i0,i0,a)', advance='no') &
            " ", trim(label(field_index)), "_A", spin, orbital, "_im"
        end do
      end do
      do component = 1, 3
        write(unit, '(3a,i0)', advance='no') &
          " ", trim(label(field_index)), "_mf", component
      end do
    end do
    write(unit, *)

    do point = 1, mesh%point_count()
      if (.not. asymptotic_ray_point(point)) cycle
      coordinate = mesh%point_coordinate(point)
      row = 0.0_rk
      row(1:3) = [coordinate, real(asymptotic_ray_id(point), rk)]
      column = 4
      call append_shadow_state_values( &
        row, column, shadow_probe_baseline, point)
      call append_shadow_state_values( &
        row, column, shadow_probe_first_mapped, point)
      call append_shadow_state_values( &
        row, column, shadow_probe_state, point)
      call append_shadow_state_values( &
        row, column, shadow_probe_mapped, point)
      write(unit, '(*(es24.16e3,1x))') row
    end do
    close(unit)
  end subroutine write_asymptotic_shadow_states


  subroutine append_shadow_state_values(row, column, field, point)
    real(rk), intent(inout) :: row(:)
    integer, intent(inout) :: column
    type(spinful_state_2d_t), intent(in) :: field
    integer, intent(in) :: point

    integer :: component, orbital, spin

    do spin = 1, 3
      do orbital = 1, 3
        row(column) = real(field%order_parameter(spin, orbital, point), rk)
        row(column + 1) = aimag(field%order_parameter(spin, orbital, point))
        column = column + 2
      end do
    end do
    do component = 1, 3
      row(column) = field%current_mean_field(component, point)
      column = column + 1
    end do
  end subroutine append_shadow_state_values


  subroutine print_summary()
    real(rk) :: maximum_coordinate(2)

    maximum_coordinate = mesh%point_coordinate(residual%maximum_point)
    print '(a,i0,a,i0)', "field grid: ", mesh%x_point_count(), " x ", &
      mesh%y_point_count()
    print '(a,i0,a,i0)', "sampled points: ", count(active_point), &
      " on MPI ranks: ", rank_count
    if (asymptotic_ray_probe_count > 0) &
      print '(a,i0)', "passive asymptotic ray probes: ", &
        count(asymptotic_ray_point)
    if (passive_probe_output_requested) &
      print '(a,a,a,l1,a,es11.3)', "passive probe output: ", &
        trim(passive_probe_source), ", map succeeded: ", &
        passive_probe_map_succeeded, ", normalization: ", &
        passive_probe_maximum_normalization_error
    if (asymptotic_probe_iterations > 0) then
      print '(a,a,a,i0,a,i0)', "shadow probe status: ", &
        trim(shadow_probe_terminal_status), ", updates: ", &
        shadow_probe_iterations_completed, ", maps: ", &
        shadow_probe_map_evaluations
      print '(a,2(es11.3,1x))', &
        "shadow terminal residual [rms,max]: ", &
        shadow_probe_residual%rms_residual, &
        shadow_probe_residual%maximum_absolute_residual
    end if
    print '(a,a,a,i0)', "terminal status: ", trim(terminal_status), &
      ", iterations: ", iterations_completed
    print '(a,i0,a,i0)', "updated points: ", count(update_point), &
      ", Anderson values: ", &
      new_src_masked_iteration_vector_size_2d(state, update_point)
    print '(a,f8.3,a,2(f8.3,1x))', &
      "measured hard-core separation a: ", &
      measured_half_core_separation, " at y-/y+: ", &
      measured_negative_half_core_y, measured_positive_half_core_y
    print '(a,a)', "endpoint policy: ", trim(endpoint_policy)
    print '(a,a,a,i0)', "inactive halo: ", trim(inactive_halo_policy), &
      ", dependent points: ", count(dependent_halo_point)
    if (any(dependent_halo_point)) &
      print '(a,2(es12.4,1x))', &
        "halo maximum change [initial,last]: ", &
        initial_halo_maximum_change, last_halo_maximum_change
    print '(a,a)', "zero-mode projection: ", trim(zero_mode_projection)
    if (use_zero_mode_projection) then
      print '(a,i0,a,i0)', "zero-mode tangents retained: ", &
        zero_mode_basis%retained_count, " of ", &
        zero_mode_basis%requested_count
      print '(a,es12.4,a,es11.3)', "final removed tangent norm: ", &
        zero_mode_report%removed_residual_norm, " fraction: ", &
        zero_mode_report%removed_fraction
    end if
    if (trim(endpoint_policy) == "state_asymptotic") then
      print '(a,a,a,es12.4)', "bulk reference: ", &
        trim(bulk_reference_policy), " phase offset: ", &
        endpoint%bulk_phase_offset
      print '(a,es12.4)', "bulk-reference annular RMS deviation: ", &
        bulk_reference_fit_rms
      print '(a,es12.4)', "bulk-reference annular maximum deviation: ", &
        bulk_reference_fit_maximum
    end if
    print '(a,i0,a,i0)', "directions: ", angular%direction_count(), &
      " poles: ", energy%pole_count()
    print '(a,3(es12.4,1x))', "residual [mean,rms,max]: ", &
      residual%mean_absolute_residual, residual%rms_residual, &
      residual%maximum_absolute_residual
    print '(a,es12.4)', "relative L2 residual: ", &
      residual%relative_l2_residual
    print '(a,es12.4,a,2(f8.2,1x))', "maximum point RMS: ", &
      residual%maximum_point_rms_residual, " at ", maximum_coordinate
    print '(a,es12.4)', "core relative L2 residual: ", &
      core_residual%relative_l2_residual
    print '(a,es12.4)', "outer relative L2 residual: ", &
      outer_residual%relative_l2_residual
    print '(a,es12.4)', "maximum normalization error: ", &
      maximum_normalization_error_over_solve
    print '(a,f10.2)', "maximum rank elapsed time [s]: ", &
      diagnostics%maximum_rank_elapsed_seconds
    print '(a,a)', "metrics: ", trim(metrics_file)
    print '(a,a)', "sampled residuals: ", trim(sample_file)
    print '(a,a)', "iteration history: ", trim(history_file)
    if (len_trim(initial_field_file) > 0) &
      print '(a,a)', "initial field: ", trim(initial_field_file)
    if (len_trim(final_field_file) > 0) &
      print '(a,a)', "final field: ", trim(final_field_file)
    if (len_trim(sparse_final_file) > 0) &
      print '(a,a)', "sparse final state: ", trim(sparse_final_file)
  end subroutine print_summary


  subroutine require_mpi(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require_mpi

end program benchmark_legacy_double_core_map_2d
