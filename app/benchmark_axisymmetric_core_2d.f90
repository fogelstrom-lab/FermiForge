program benchmark_axisymmetric_core_2d
  use mpi_f08
  use he3_kinds, only : rk
  use radial_mesh, only : radial_mesh_t
  use axial_radial_embedding, only : axial_radial_profile_t, &
                                     embed_axial_radial_profile_2d, &
                                     sample_axial_radial_profile
  use legacy_radial_profile_io, only : read_new_src_radial_profile
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, &
                                make_uniform_cartesian_mesh, &
                                make_symmetric_multiscale_cartesian_mesh, &
                                make_circular_active_mask
  use spinful_state_2d, only : spinful_state_2d_t, allocate_spinful_state_2d
  use spinful_field_io_2d, only : read_spinful_field_map_2d, &
                                  write_spinful_field_map_2d
  use he3_quadrature, only : angular_quadrature_3d_t, ozaki_quadrature_t, &
                             read_legacy_gauss_table, &
                             read_legacy_ozaki_table
  use he3_mpi_field_map_2d, only : he3_mpi_field_map_diagnostics_t, &
                                   evaluate_mpi_he3_field_map, &
                                   update_mpi_he3_state_with_anderson
  use he3_nonlinear_iteration_2d, only : he3_field_residual_report_t, &
                                         compute_he3_field_residual
  use new_src_iteration_layout_2d, only : &
    new_src_masked_iteration_vector_size_2d
  use legacy_anderson_mixing, only : legacy_anderson_t, anderson_report_t
  use free_vortex_asymptotic_2d, only : free_vortex_endpoint_2d_t, &
    set_axisymmetric_radial_reference, set_axisymmetric_asymptotic_fit, &
    apply_free_vortex_asymptotic_halo_2d
  implicit none

  type :: area_weighted_residual_t
    real(rk) :: rms = 0.0_rk
    real(rk) :: relative_l2 = 0.0_rk
    real(rk) :: sampled_area = 0.0_rk
  end type area_weighted_residual_t

  type(angular_quadrature_3d_t) :: angular
  type(cartesian_mesh_2d_t) :: mesh
  type(free_vortex_endpoint_2d_t) :: endpoint
  type(he3_field_residual_report_t) :: fine_residual, medium_residual
  type(he3_field_residual_report_t) :: outer_residual, residual
  type(area_weighted_residual_t) :: weighted_fine, weighted_medium
  type(area_weighted_residual_t) :: weighted_outer, weighted_residual
  type(he3_mpi_field_map_diagnostics_t) :: diagnostics
  type(anderson_report_t) :: anderson_report
  type(legacy_anderson_t) :: accelerator
  type(ozaki_quadrature_t) :: energy
  type(axial_radial_profile_t) :: profile
  type(radial_mesh_t) :: radial_grid
  type(spinful_state_2d_t) :: evaluated_state, mapped, next
  type(spinful_state_2d_t) :: reference_state, state

  character(len=32) :: detected_grid_kind, endpoint_policy, mesh_kind
  character(len=16) :: radial_grid_kind
  character(len=32) :: initialization_mode, terminal_status
  character(len=64) :: core_label
  character(len=512) :: current_field_file, gauss_file, input_file
  character(len=512) :: checkpoint_file, final_output_file, history_file
  character(len=512) :: input_message, input_output_file, mapped_output_file
  character(len=512) :: metrics_file, order_parameter_file, output_override
  character(len=512) :: ozaki_file, restart_field_file
  integer :: anderson_history_limit, checkpoint_interval, history_unit
  integer :: ierr, input_status, input_unit, iteration, iterations_completed
  integer :: map_evaluations, maximum_iterations, minimum_substep_count
  integer :: number_of_cells, rank, rank_count, result_integer
  real(rk) :: asymptotic_bulk_gap, asymptotic_center_x, asymptotic_center_y
  real(rk) :: asymptotic_fit_inner_radius, asymptotic_matching_radius
  real(rk) :: asymptotic_maximum_step
  real(rk) :: asymptotic_outer_radius
  real(rk) :: asymptotic_winding, boundary_relaxation_distance
  real(rk) :: feedback_scale, half_width, input_current_embedding_error
  real(rk) :: input_gap_embedding_error, mapped_current_radial_error
  real(rk) :: mapped_gap_radial_error, maximum_embedding_error
  real(rk) :: initial_current_radial_error, initial_gap_radial_error
  real(rk) :: final_current_radial_error, final_gap_radial_error
  real(rk) :: input_current_embedding_relative_l2
  real(rk) :: input_gap_embedding_relative_l2
  real(rk) :: initial_current_radial_relative_l2
  real(rk) :: initial_gap_radial_relative_l2
  real(rk) :: final_current_radial_relative_l2
  real(rk) :: final_gap_radial_relative_l2
  real(rk) :: final_halo_current_relative_error
  real(rk) :: final_halo_current_relative_l2
  real(rk) :: final_halo_gap_relative_error
  real(rk) :: final_halo_gap_relative_l2
  real(rk) :: mapped_current_radial_relative_l2
  real(rk) :: mapped_gap_radial_relative_l2
  real(rk) :: active_radius, coarse_spacing, fine_region_half_width
  real(rk) :: fine_spacing, medium_region_half_width, medium_spacing
  real(rk) :: anderson_maximum_mixing, anderson_progress_threshold
  real(rk) :: convergence_tolerance, cumulative_map_seconds
  real(rk) :: initial_perturbation_amplitude, initial_perturbation_radius
  real(rk) :: maximum_normalization_error_over_solve
  real(rk) :: trajectory_maximum_step
  logical, allocatable :: active_point(:), fine_point(:)
  logical, allocatable :: medium_point(:), outer_point(:)
  logical :: benchmark_passed, history_is_open, map_succeeded, solve_requested

  namelist /axisymmetric_core_benchmark/ core_label, half_width, &
    number_of_cells, mesh_kind, fine_region_half_width, &
    medium_region_half_width, fine_spacing, medium_spacing, coarse_spacing, &
    active_radius, feedback_scale, trajectory_maximum_step, &
    minimum_substep_count, boundary_relaxation_distance, &
    maximum_embedding_error, order_parameter_file, current_field_file, &
    radial_grid_kind, gauss_file, ozaki_file, endpoint_policy, &
    asymptotic_center_x, asymptotic_center_y, asymptotic_bulk_gap, &
    asymptotic_winding, asymptotic_fit_inner_radius, &
    asymptotic_matching_radius, &
    asymptotic_outer_radius, asymptotic_maximum_step, &
    initialization_mode, restart_field_file, initial_perturbation_amplitude, &
    initial_perturbation_radius, maximum_iterations, convergence_tolerance, &
    anderson_history_limit, anderson_progress_threshold, &
    anderson_maximum_mixing, checkpoint_interval, input_output_file, &
    mapped_output_file, final_output_file, metrics_file, history_file, &
    checkpoint_file

  core_label = "normal-core"
  half_width = 8.0_rk
  number_of_cells = 16
  mesh_kind = "uniform"
  fine_region_half_width = 2.0_rk
  medium_region_half_width = 5.0_rk
  fine_spacing = 0.5_rk
  medium_spacing = 1.0_rk
  coarse_spacing = 2.0_rk
  active_radius = 0.0_rk
  feedback_scale = 5.4_rk / (1.0_rk + 5.4_rk / 3.0_rk)
  trajectory_maximum_step = 0.25_rk
  minimum_substep_count = 1
  boundary_relaxation_distance = 1.75_rk
  maximum_embedding_error = 1.0e-12_rk
  radial_grid_kind = "auto"
  endpoint_policy = "radial_reference"
  asymptotic_center_x = 0.0_rk
  asymptotic_center_y = 0.0_rk
  asymptotic_bulk_gap = 0.2799_rk
  asymptotic_winding = 1.0_rk
  asymptotic_fit_inner_radius = 0.0_rk
  asymptotic_matching_radius = 0.0_rk
  asymptotic_outer_radius = 70.0_rk
  asymptotic_maximum_step = 0.25_rk
  order_parameter_file = "benchmarks/normal_core_2d/reference/op_xyz"
  current_field_file = "benchmarks/normal_core_2d/reference/curr"
  gauss_file = "new_src/gauss11.dat"
  ozaki_file = "new_src/ozaki.dat"
  input_output_file = "work/normal-core-2d/input_fields_2d.dat"
  mapped_output_file = "work/normal-core-2d/mapped_fields_2d.dat"
  final_output_file = "work/normal-core-2d/final_fields_2d.dat"
  metrics_file = "work/normal-core-2d/metrics.txt"
  history_file = "work/normal-core-2d/iteration_history.dat"
  checkpoint_file = "work/normal-core-2d/checkpoint_fields_2d.dat"
  initialization_mode = "radial_reference"
  restart_field_file = ""
  initial_perturbation_amplitude = 0.0_rk
  initial_perturbation_radius = 3.0_rk
  maximum_iterations = 0
  convergence_tolerance = 2.0e-6_rk
  anderson_history_limit = 10
  anderson_progress_threshold = 0.1_rk
  anderson_maximum_mixing = 5.0_rk
  checkpoint_interval = 5

  call MPI_Init(ierr)
  call require_mpi(ierr == MPI_SUCCESS, "MPI initialization failed")
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call require_mpi(ierr == MPI_SUCCESS, "MPI rank query failed")
  call MPI_Comm_size(MPI_COMM_WORLD, rank_count, ierr)
  call require_mpi(ierr == MPI_SUCCESS, "MPI size query failed")

  call get_command_argument(1, input_file)
  if (len_trim(input_file) == 0) &
    input_file = "examples/2d_normal_core_benchmark.nml"
  open(newunit=input_unit, file=trim(input_file), status="old", action="read", &
       iostat=input_status, iomsg=input_message)
  if (input_status /= 0) &
    error stop "cannot open axisymmetric benchmark input: " // &
      trim(input_message)
  read(input_unit, nml=axisymmetric_core_benchmark, &
       iostat=input_status, iomsg=input_message)
  close(input_unit)
  if (input_status /= 0) &
    error stop "cannot read axisymmetric benchmark input: " // &
      trim(input_message)
  call get_command_argument(2, output_override)
  if (len_trim(output_override) > 0) input_output_file = output_override
  call get_command_argument(3, output_override)
  if (len_trim(output_override) > 0) mapped_output_file = output_override
  call get_command_argument(4, output_override)
  if (len_trim(output_override) > 0) metrics_file = output_override
  call get_command_argument(5, output_override)
  if (len_trim(output_override) > 0) final_output_file = output_override
  call get_command_argument(6, output_override)
  if (len_trim(output_override) > 0) history_file = output_override
  call get_command_argument(7, output_override)
  if (len_trim(output_override) > 0) checkpoint_file = output_override
  call configure_endpoint()
  call validate_input()

  call read_new_src_radial_profile( &
    trim(order_parameter_file), trim(current_field_file), radial_grid, &
    profile, radial_grid_kind=trim(radial_grid_kind), &
    detected_grid_kind=detected_grid_kind)
  if (trim(endpoint_policy) == "radial_reference") then
    if (asymptotic_matching_radius <= 0.0_rk) &
      asymptotic_matching_radius = &
        radial_grid%r(radial_grid%point_count() - 6)
    call set_axisymmetric_radial_reference( &
      endpoint, radial_grid, profile, asymptotic_matching_radius)
  else if (trim(endpoint_policy) == "radial_asymptotic") then
    call set_axisymmetric_asymptotic_fit( &
      endpoint, radial_grid, profile, asymptotic_fit_inner_radius, &
      asymptotic_matching_radius)
  end if
  if (sqrt(2.0_rk) * half_width > radial_grid%outer_radius()) &
    error stop "benchmark Cartesian mesh extends beyond the radial profile"
  call read_legacy_gauss_table(trim(gauss_file), 48, 11, angular)
  call read_legacy_ozaki_table(trim(ozaki_file), energy)
  call configure_mesh()
  if (active_radius > 0.0_rk) then
    call make_circular_active_mask( &
      mesh, [asymptotic_center_x, asymptotic_center_y], active_radius, &
      active_point)
  else
    allocate(active_point(mesh%point_count()), source=.true.)
  end if
  call make_residual_zone_masks()
  call allocate_spinful_state_2d(mesh, reference_state)
  call embed_axial_radial_profile_2d( &
    radial_grid, profile, asymptotic_winding, mesh, reference_state)
  call radial_reference_error( &
    mesh, radial_grid, profile, reference_state, asymptotic_winding, &
    active_point, input_gap_embedding_error, input_current_embedding_error, &
    input_gap_embedding_relative_l2, input_current_embedding_relative_l2)

  call initialize_state()
  call refresh_dependent_halo(state)
  call radial_reference_error( &
    mesh, radial_grid, profile, state, asymptotic_winding, active_point, &
    initial_gap_radial_error, initial_current_radial_error, &
    initial_gap_radial_relative_l2, initial_current_radial_relative_l2)

  solve_requested = maximum_iterations > 0
  terminal_status = "single_map"
  iterations_completed = 0
  map_evaluations = 0
  cumulative_map_seconds = 0.0_rk
  maximum_normalization_error_over_solve = 0.0_rk
  history_unit = -1
  history_is_open = .false.
  if (rank == 0) then
    call write_spinful_field_map_2d( &
      trim(input_output_file), mesh, state, &
      trim(core_label) // "_initial_state", trim(endpoint_policy))
    call open_iteration_history()
    if (solve_requested) call accelerator%initialize( &
      new_src_masked_iteration_vector_size_2d(state, active_point), &
      anderson_history_limit, &
      anderson_progress_threshold, anderson_maximum_mixing)
  end if

  do iteration = 1, max(1, maximum_iterations)
    evaluated_state = state
    call evaluate_mpi_he3_field_map( &
      MPI_COMM_WORLD, mesh, state, angular, energy, feedback_scale, &
      trajectory_maximum_step, minimum_substep_count, &
      boundary_relaxation_distance, mapped, diagnostics, map_succeeded, &
      free_vortex_endpoint=endpoint, active_point_mask=active_point)
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

    if (.not. solve_requested) then
      if (rank == 0) then
        call compute_he3_field_residual( &
          state, mapped, residual, active_point_mask=active_point)
        call write_single_map_history(iteration)
      end if
      exit
    end if

    call update_mpi_he3_state_with_anderson( &
      MPI_COMM_WORLD, mesh, accelerator, state, mapped, &
      convergence_tolerance, next, anderson_report, residual, active_point)
    iterations_completed = iteration
    if (rank == 0) then
      call compute_he3_field_residual( &
        state, mapped, residual, active_point_mask=active_point)
      call write_iteration_history(iteration)
      call print_iteration(iteration)
    end if

    if (anderson_report%converged) then
      terminal_status = "converged"
      state = next
      call refresh_dependent_halo(state)
      if (rank == 0) call write_checkpoint(state, iteration)
      exit
    end if
    if (iteration == maximum_iterations) then
      terminal_status = "iteration_limit"
      state = next
      call refresh_dependent_halo(state)
      if (rank == 0) call write_checkpoint(next, iteration)
      exit
    end if
    state = next
    call refresh_dependent_halo(state)
    if (rank == 0 .and. checkpoint_interval > 0) then
      if (modulo(iteration, checkpoint_interval) == 0) &
        call write_checkpoint(state, iteration)
    end if
  end do

  benchmark_passed = map_succeeded
  if (solve_requested) benchmark_passed = benchmark_passed .and. &
    trim(terminal_status) == "converged"
  if (rank == 0) then
    if (history_is_open) close(history_unit)
    call radial_reference_error( &
      mesh, radial_grid, profile, state, asymptotic_winding, active_point, &
      final_gap_radial_error, final_current_radial_error, &
      final_gap_radial_relative_l2, final_current_radial_relative_l2)
    call radial_reference_error( &
      mesh, radial_grid, profile, mapped, asymptotic_winding, active_point, &
      mapped_gap_radial_error, mapped_current_radial_error, &
      mapped_gap_radial_relative_l2, mapped_current_radial_relative_l2)
    final_halo_gap_relative_error = 0.0_rk
    final_halo_current_relative_error = 0.0_rk
    final_halo_gap_relative_l2 = 0.0_rk
    final_halo_current_relative_l2 = 0.0_rk
    if (any(.not. active_point)) call radial_reference_error( &
      mesh, radial_grid, profile, state, asymptotic_winding, &
      .not. active_point, final_halo_gap_relative_error, &
      final_halo_current_relative_error, final_halo_gap_relative_l2, &
      final_halo_current_relative_l2)
    call compute_final_residuals()
    benchmark_passed = benchmark_passed .and. &
      input_gap_embedding_error <= maximum_embedding_error .and. &
      input_current_embedding_error <= maximum_embedding_error .and. &
      maximum_normalization_error_over_solve <= 1.0e-10_rk
    call write_spinful_field_map_2d( &
      trim(mapped_output_file), mesh, mapped, &
      trim(core_label) // "_last_quasiclassical_map", trim(endpoint_policy))
    call write_spinful_field_map_2d( &
      trim(final_output_file), mesh, state, &
      trim(core_label) // "_final_evaluated_state", trim(endpoint_policy))
    call write_metrics()
    call print_metrics()
  end if

  result_integer = merge(1, 0, benchmark_passed)
  call MPI_Bcast(result_integer, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call require_mpi(ierr == MPI_SUCCESS, "benchmark result broadcast failed")
  benchmark_passed = result_integer == 1
  call MPI_Finalize(ierr)
  call require_mpi(ierr == MPI_SUCCESS, "MPI finalization failed")
  if (.not. benchmark_passed) &
    error stop "axisymmetric core benchmark acceptance check failed"

contains

  subroutine initialize_state()
    real(rk) :: coordinate(2), radius, scale
    integer :: point

    call allocate_spinful_state_2d(mesh, state)
    select case (trim(adjustl(initialization_mode)))
    case ("radial_reference")
      initialization_mode = "radial_reference"
      state = reference_state
      if (abs(initial_perturbation_amplitude) <= tiny(1.0_rk)) return
      do point = 1, mesh%point_count()
        if (.not. active_point(point)) cycle
        coordinate = mesh%point_coordinate(point) - &
          [asymptotic_center_x, asymptotic_center_y]
        radius = sqrt(sum(coordinate**2))
        scale = 1.0_rk - initial_perturbation_amplitude * &
          exp(-(radius / initial_perturbation_radius)**2)
        state%order_parameter(:, :, point) = &
          cmplx(scale, 0.0_rk, kind=rk) * &
          state%order_parameter(:, :, point)
        state%current_mean_field(:, point) = &
          scale * state%current_mean_field(:, point)
      end do
    case ("restart")
      initialization_mode = "restart"
      call read_spinful_field_map_2d(trim(restart_field_file), mesh, state)
    case default
      error stop "initialization_mode must be 'radial_reference' or 'restart'"
    end select
  end subroutine initialize_state


  subroutine open_iteration_history()
    character(len=512) :: message
    integer :: status

    open(newunit=history_unit, file=trim(history_file), status="replace", &
         action="write", iostat=status, iomsg=message)
    if (status /= 0) &
      error stop "cannot open iteration history: " // trim(message)
    history_is_open = .true.
    write(history_unit, '(a)') &
      "# iteration map_rms map_max map_relative_l2 legacy_scaled_max " // &
      "anderson_vector_norm anderson_relative p history discarded " // &
      "stochastic_fallback normalization_error map_seconds"
  end subroutine open_iteration_history


  subroutine write_single_map_history(iteration_number)
    integer, intent(in) :: iteration_number

    write(history_unit, '(i8,1x,7(es24.16e3,1x),3(i8,1x),2(es24.16e3,1x))') &
      iteration_number, residual%rms_residual, &
      residual%maximum_absolute_residual, residual%relative_l2_residual, &
      residual%legacy_scaled_maximum_residual, 0.0_rk, 0.0_rk, 0.0_rk, &
      0, 0, 0, diagnostics%global%maximum_normalization_error, &
      diagnostics%maximum_rank_elapsed_seconds
    flush(history_unit)
  end subroutine write_single_map_history


  subroutine write_iteration_history(iteration_number)
    integer, intent(in) :: iteration_number

    write(history_unit, '(i8,1x,7(es24.16e3,1x),3(i8,1x),2(es24.16e3,1x))') &
      iteration_number, residual%rms_residual, &
      residual%maximum_absolute_residual, residual%relative_l2_residual, &
      residual%legacy_scaled_maximum_residual, &
      anderson_report%residual_norm, &
      anderson_report%legacy_relative_residual, anderson_report%mixing, &
      anderson_report%history_size, anderson_report%discarded_vectors, &
      merge(1, 0, anderson_report%used_stochastic_step), &
      diagnostics%global%maximum_normalization_error, &
      diagnostics%maximum_rank_elapsed_seconds
    flush(history_unit)
  end subroutine write_iteration_history


  subroutine print_iteration(iteration_number)
    integer, intent(in) :: iteration_number

    print '(a,i5,a,2(es11.3,1x),a,f6.2,a,i3,a,f9.2)', &
      "solved (2D AA): ", iteration_number, " error(avg,max)= (", &
      residual%mean_absolute_residual, residual%maximum_absolute_residual, &
      ") p=", anderson_report%mixing, &
      " md=", anderson_report%history_size, &
      " map[s]=", diagnostics%maximum_rank_elapsed_seconds
  end subroutine print_iteration


  subroutine write_checkpoint(field, iteration_number)
    type(spinful_state_2d_t), intent(in) :: field
    integer, intent(in) :: iteration_number

    character(len=64) :: iteration_label

    write(iteration_label, '(a,i0)') "normal_core_checkpoint_iteration_", &
      iteration_number
    call write_spinful_field_map_2d( &
      trim(checkpoint_file), mesh, field, trim(iteration_label), &
      trim(endpoint_policy))
  end subroutine write_checkpoint


  subroutine compute_final_residuals()
    call compute_he3_field_residual( &
      evaluated_state, mapped, residual, active_point_mask=active_point)
    call compute_he3_field_residual( &
      evaluated_state, mapped, fine_residual, active_point_mask=fine_point)
    call compute_he3_field_residual( &
      evaluated_state, mapped, medium_residual, &
      active_point_mask=medium_point)
    call compute_he3_field_residual( &
      evaluated_state, mapped, outer_residual, active_point_mask=outer_point)
    call compute_area_weighted_residual( &
      mesh, evaluated_state, mapped, active_point, weighted_residual)
    call compute_area_weighted_residual( &
      mesh, evaluated_state, mapped, fine_point, weighted_fine)
    call compute_area_weighted_residual( &
      mesh, evaluated_state, mapped, medium_point, weighted_medium)
    call compute_area_weighted_residual( &
      mesh, evaluated_state, mapped, outer_point, weighted_outer)
  end subroutine compute_final_residuals

  subroutine make_residual_zone_masks()
    real(rk) :: coordinate(2), radius
    integer :: point

    allocate(fine_point(mesh%point_count()), source=.false.)
    allocate(medium_point(mesh%point_count()), source=.false.)
    allocate(outer_point(mesh%point_count()), source=.false.)
    do point = 1, mesh%point_count()
      if (.not. active_point(point)) cycle
      coordinate = mesh%point_coordinate(point) - &
        [asymptotic_center_x, asymptotic_center_y]
      radius = sqrt(sum(coordinate**2))
      if (radius <= fine_region_half_width) then
        fine_point(point) = .true.
      else if (radius <= medium_region_half_width) then
        medium_point(point) = .true.
      else
        outer_point(point) = .true.
      end if
    end do
    if (.not. any(fine_point) .or. .not. any(medium_point) .or. &
        .not. any(outer_point)) &
      error stop "benchmark residual zones must each contain active points"
  end subroutine make_residual_zone_masks


  subroutine compute_area_weighted_residual( &
      field_mesh, current, mapped_field, mask, report)
    type(cartesian_mesh_2d_t), intent(in) :: field_mesh
    type(spinful_state_2d_t), intent(in) :: current, mapped_field
    logical, intent(in) :: mask(:)
    type(area_weighted_residual_t), intent(out) :: report

    real(rk) :: current_squared, difference_squared, weight, weight_x, weight_y
    integer :: component, orbital, point, spin, x_node, y_node

    if (size(mask) /= field_mesh%point_count() .or. .not. any(mask)) &
      error stop "invalid area-weighted residual mask"
    report = area_weighted_residual_t()
    difference_squared = 0.0_rk
    current_squared = 0.0_rk
    do y_node = 0, field_mesh%y_cell_count()
      weight_y = nodal_axis_weight(field_mesh, y_node, .false.)
      do x_node = 0, field_mesh%x_cell_count()
        point = field_mesh%point_index(x_node, y_node)
        if (.not. mask(point)) cycle
        weight_x = nodal_axis_weight(field_mesh, x_node, .true.)
        weight = weight_x * weight_y
        report%sampled_area = report%sampled_area + weight
        do spin = 1, 3
          do orbital = 1, 3
            difference_squared = difference_squared + weight * ( &
              real(mapped_field%order_parameter(spin, orbital, point) - &
                   current%order_parameter(spin, orbital, point), rk)**2 + &
              aimag(mapped_field%order_parameter(spin, orbital, point) - &
                    current%order_parameter(spin, orbital, point))**2)
            current_squared = current_squared + weight * &
              abs(current%order_parameter(spin, orbital, point))**2
          end do
        end do
        do component = 1, 3
          difference_squared = difference_squared + weight * ( &
            mapped_field%current_mean_field(component, point) - &
            current%current_mean_field(component, point))**2
          current_squared = current_squared + weight * &
            current%current_mean_field(component, point)**2
        end do
      end do
    end do
    report%rms = sqrt(difference_squared / &
      (21.0_rk * report%sampled_area))
    report%relative_l2 = sqrt(difference_squared / &
      max(current_squared, tiny(1.0_rk)))
  end subroutine compute_area_weighted_residual


  pure real(rk) function nodal_axis_weight( &
      field_mesh, node, use_x_axis) result(weight)
    type(cartesian_mesh_2d_t), intent(in) :: field_mesh
    integer, intent(in) :: node
    logical, intent(in) :: use_x_axis

    integer :: cell_count

    if (use_x_axis) then
      cell_count = field_mesh%x_cell_count()
      if (node == 0) then
        weight = 0.5_rk * field_mesh%x_cell_spacing(0)
      else if (node == cell_count) then
        weight = 0.5_rk * field_mesh%x_cell_spacing(cell_count - 1)
      else
        weight = 0.5_rk * (field_mesh%x_cell_spacing(node - 1) + &
                           field_mesh%x_cell_spacing(node))
      end if
    else
      cell_count = field_mesh%y_cell_count()
      if (node == 0) then
        weight = 0.5_rk * field_mesh%y_cell_spacing(0)
      else if (node == cell_count) then
        weight = 0.5_rk * field_mesh%y_cell_spacing(cell_count - 1)
      else
        weight = 0.5_rk * (field_mesh%y_cell_spacing(node - 1) + &
                           field_mesh%y_cell_spacing(node))
      end if
    end if
  end function nodal_axis_weight

  subroutine configure_mesh()
    select case (trim(adjustl(mesh_kind)))
    case ("uniform")
      mesh_kind = "uniform"
      call make_uniform_cartesian_mesh( &
        -half_width, half_width, number_of_cells, &
        -half_width, half_width, number_of_cells, mesh)
    case ("multiscale")
      mesh_kind = "multiscale"
      call make_symmetric_multiscale_cartesian_mesh( &
        half_width, fine_region_half_width, medium_region_half_width, &
        fine_spacing, medium_spacing, coarse_spacing, mesh)
    case default
      error stop "mesh_kind must be 'uniform' or 'multiscale'"
    end select
  end subroutine configure_mesh

  subroutine configure_endpoint()
    select case (trim(adjustl(endpoint_policy)))
    case ("local")
      endpoint%enabled = .false.
      endpoint_policy = "local"
    case ("free_vortex")
      endpoint%enabled = .true.
      endpoint_policy = "free_vortex"
    case ("radial_reference")
      endpoint%enabled = .true.
      endpoint_policy = "radial_reference"
    case ("radial_asymptotic")
      endpoint%enabled = .true.
      endpoint_policy = "radial_asymptotic"
      if (asymptotic_matching_radius <= 0.0_rk) then
        if (active_radius > 0.0_rk) then
          asymptotic_matching_radius = active_radius
        else
          asymptotic_matching_radius = 0.8_rk * half_width
        end if
      end if
      if (asymptotic_fit_inner_radius <= 0.0_rk) &
        asymptotic_fit_inner_radius = 0.75_rk * &
          asymptotic_matching_radius
    case default
      error stop &
        "endpoint_policy must be local, free_vortex, radial_reference, or radial_asymptotic"
    end select
    endpoint%center = [asymptotic_center_x, asymptotic_center_y]
    endpoint%bulk_gap = asymptotic_bulk_gap
    endpoint%winding = asymptotic_winding
    endpoint%outer_radius = asymptotic_outer_radius
    endpoint%maximum_step = asymptotic_maximum_step
  end subroutine configure_endpoint


  subroutine validate_input()
    if (half_width <= 0.0_rk) &
      error stop "benchmark mesh half width must be positive"
    select case (trim(adjustl(mesh_kind)))
    case ("uniform")
      if (number_of_cells < 2 .or. modulo(number_of_cells, 2) /= 0) &
        error stop "uniform benchmark mesh needs an even cell count"
    case ("multiscale")
      if (fine_region_half_width <= 0.0_rk .or. &
          medium_region_half_width <= fine_region_half_width .or. &
          half_width <= medium_region_half_width .or. &
          fine_spacing <= 0.0_rk .or. medium_spacing < fine_spacing .or. &
          coarse_spacing < medium_spacing) &
        error stop "invalid multiscale benchmark mesh controls"
    case default
      error stop "mesh_kind must be 'uniform' or 'multiscale'"
    end select
    if (active_radius < 0.0_rk .or. active_radius > half_width) &
      error stop "active radius must be zero or lie inside the mesh halo"
    if (trajectory_maximum_step <= 0.0_rk .or. &
        minimum_substep_count < 1 .or. boundary_relaxation_distance < 0.0_rk) &
      error stop "invalid benchmark trajectory controls"
    if (maximum_embedding_error <= 0.0_rk) &
      error stop "benchmark embedding tolerance must be positive"
    if (maximum_iterations < 0) &
      error stop "maximum_iterations cannot be negative"
    if (convergence_tolerance < 0.0_rk) &
      error stop "convergence_tolerance cannot be negative"
    if (maximum_iterations > 0 .and. convergence_tolerance <= 0.0_rk) &
      error stop "a self-consistent solve needs a positive tolerance"
    if (anderson_history_limit < 1) &
      error stop "anderson_history_limit must be positive"
    if (anderson_progress_threshold <= 0.0_rk .or. &
        anderson_progress_threshold > 1.0_rk) &
      error stop "anderson_progress_threshold must be in (0,1]"
    if (anderson_maximum_mixing < 0.01_rk) &
      error stop "anderson_maximum_mixing cannot be smaller than 0.01"
    if (checkpoint_interval < 0) &
      error stop "checkpoint_interval cannot be negative"
    if (trim(endpoint_policy) == "radial_asymptotic") then
      if (asymptotic_fit_inner_radius <= 0.0_rk .or. &
          asymptotic_matching_radius <= asymptotic_fit_inner_radius .or. &
          asymptotic_matching_radius > half_width) &
        error stop "invalid radial-asymptotic fit annulus"
      if (active_radius > 0.0_rk .and. &
          active_radius < asymptotic_matching_radius) &
        error stop "active radius must enclose the asymptotic matching circle"
    end if
    select case (trim(adjustl(initialization_mode)))
    case ("radial_reference")
      if (abs(initial_perturbation_amplitude) >= 1.0_rk) &
        error stop "initial perturbation amplitude must have magnitude below one"
      if (abs(initial_perturbation_amplitude) > 0.0_rk .and. &
          initial_perturbation_radius <= 0.0_rk) &
        error stop "a nonzero initial perturbation needs a positive radius"
    case ("restart")
      if (len_trim(restart_field_file) == 0) &
        error stop "restart initialization needs restart_field_file"
      if (abs(initial_perturbation_amplitude) > 0.0_rk) &
        error stop "restart initialization cannot also apply a perturbation"
    case default
      error stop "initialization_mode must be 'radial_reference' or 'restart'"
    end select
    if (len_trim(input_output_file) == 0 .or. &
        len_trim(mapped_output_file) == 0 .or. &
        len_trim(final_output_file) == 0 .or. &
        len_trim(metrics_file) == 0 .or. len_trim(history_file) == 0 .or. &
        len_trim(checkpoint_file) == 0) &
      error stop "benchmark output paths cannot be blank"
    if (endpoint%enabled .and. &
        (asymptotic_bulk_gap <= 0.0_rk .or. &
         asymptotic_outer_radius <= sqrt(2.0_rk) * half_width .or. &
         asymptotic_maximum_step <= 0.0_rk)) &
      error stop "invalid free-vortex benchmark endpoint controls"
  end subroutine validate_input


  subroutine refresh_dependent_halo(field)
    type(spinful_state_2d_t), intent(inout) :: field

    type(spinful_state_2d_t) :: completed

    if (trim(endpoint_policy) /= "radial_asymptotic") return
    if (all(active_point)) return
    call apply_free_vortex_asymptotic_halo_2d( &
      mesh, field, active_point, endpoint, completed)
    field = completed
  end subroutine refresh_dependent_halo


  subroutine radial_reference_error( &
      field_mesh, source_mesh, source_profile, field, winding, mask, &
      gap_error, current_error, gap_relative_l2, current_relative_l2)
    type(cartesian_mesh_2d_t), intent(in) :: field_mesh
    type(radial_mesh_t), intent(in) :: source_mesh
    type(axial_radial_profile_t), intent(in) :: source_profile
    type(spinful_state_2d_t), intent(in) :: field
    real(rk), intent(in) :: winding
    logical, intent(in) :: mask(:)
    real(rk), intent(out) :: gap_error, current_error
    real(rk), intent(out) :: gap_relative_l2, current_relative_l2

    complex(rk) :: expected_gap(3, 3)
    real(rk) :: current_difference_squared, current_reference_squared
    real(rk) :: expected_current(3), gap_difference_squared
    real(rk) :: gap_reference_squared, gap_scale, current_scale
    real(rk) :: weight, weight_x, weight_y, x, y
    integer :: point, x_node, y_node
    logical :: inside

    if (size(mask) /= field_mesh%point_count() .or. .not. any(mask)) &
      error stop "invalid radial-reference error mask"
    gap_error = 0.0_rk
    current_error = 0.0_rk
    gap_scale = 0.0_rk
    current_scale = 0.0_rk
    gap_difference_squared = 0.0_rk
    gap_reference_squared = 0.0_rk
    current_difference_squared = 0.0_rk
    current_reference_squared = 0.0_rk
    do y_node = 0, field_mesh%y_cell_count()
      y = field_mesh%y_coordinate(y_node)
      weight_y = nodal_axis_weight(field_mesh, y_node, .false.)
      do x_node = 0, field_mesh%x_cell_count()
        x = field_mesh%x_coordinate(x_node)
        point = field_mesh%point_index(x_node, y_node)
        if (.not. mask(point)) cycle
        call sample_axial_radial_profile( &
          source_mesh, source_profile, x, y, winding, expected_gap, &
          expected_current, inside)
        if (.not. inside) &
          error stop "benchmark reference sample left the radial domain"
        gap_error = max(gap_error, maxval(abs( &
          field%order_parameter(:, :, point) - expected_gap)))
        current_error = max(current_error, maxval(abs( &
          field%current_mean_field(:, point) - expected_current)))
        gap_scale = max(gap_scale, maxval(abs(expected_gap)))
        current_scale = max(current_scale, maxval(abs(expected_current)))
        weight_x = nodal_axis_weight(field_mesh, x_node, .true.)
        weight = weight_x * weight_y
        gap_difference_squared = gap_difference_squared + weight * &
          sum(abs(field%order_parameter(:, :, point) - expected_gap)**2)
        gap_reference_squared = gap_reference_squared + weight * &
          sum(abs(expected_gap)**2)
        current_difference_squared = current_difference_squared + weight * &
          sum((field%current_mean_field(:, point) - expected_current)**2)
        current_reference_squared = current_reference_squared + weight * &
          sum(expected_current**2)
      end do
    end do
    gap_error = gap_error / max(gap_scale, tiny(1.0_rk))
    current_error = current_error / max(current_scale, tiny(1.0_rk))
    gap_relative_l2 = sqrt(gap_difference_squared / &
      max(gap_reference_squared, tiny(1.0_rk)))
    current_relative_l2 = sqrt(current_difference_squared / &
      max(current_reference_squared, tiny(1.0_rk)))
  end subroutine radial_reference_error


  subroutine write_metrics()
    character(len=512) :: message
    integer :: status, unit

    open(newunit=unit, file=trim(metrics_file), status="replace", &
         action="write", iostat=status, iomsg=message)
    if (status /= 0) &
      error stop "cannot open benchmark metrics: " // trim(message)
    write(unit, '(a)') "benchmark=axisymmetric_core_2d"
    write(unit, '(a,a)') "core_label=", trim(core_label)
    write(unit, '(a,a)') "radial_grid_kind=", trim(detected_grid_kind)
    write(unit, '(a,a)') "mesh_kind=", trim(mesh_kind)
    write(unit, '(a,a)') "endpoint_policy=", trim(endpoint_policy)
    write(unit, '(a,a)') "initialization_mode=", trim(initialization_mode)
    write(unit, '(a,a)') "terminal_status=", trim(terminal_status)
    write(unit, '(a,i0)') "mpi_ranks=", rank_count
    write(unit, '(a,i0)') "maximum_iterations=", maximum_iterations
    write(unit, '(a,i0)') "iterations_completed=", iterations_completed
    write(unit, '(a,i0)') "map_evaluations=", map_evaluations
    write(unit, '(a,es24.16e3)') "convergence_tolerance=", &
      convergence_tolerance
    write(unit, '(a,i0)') "anderson_history_limit=", &
      anderson_history_limit
    write(unit, '(a,es24.16e3)') "anderson_progress_threshold=", &
      anderson_progress_threshold
    write(unit, '(a,es24.16e3)') "anderson_maximum_mixing=", &
      anderson_maximum_mixing
    write(unit, '(a,es24.16e3)') "initial_perturbation_amplitude=", &
      initial_perturbation_amplitude
    write(unit, '(a,es24.16e3)') "initial_perturbation_radius=", &
      initial_perturbation_radius
    write(unit, '(a,i0)') "cartesian_points=", mesh%point_count()
    write(unit, '(a,i0)') "active_points=", count(active_point)
    if (trim(endpoint_policy) == "radial_asymptotic") then
      write(unit, '(a,i0)') "frozen_halo_points=", 0
      write(unit, '(a,i0)') "dependent_asymptotic_halo_points=", &
        mesh%point_count() - count(active_point)
    else
      write(unit, '(a,i0)') "frozen_halo_points=", &
        mesh%point_count() - count(active_point)
      write(unit, '(a,i0)') "dependent_asymptotic_halo_points=", 0
    end if
    write(unit, '(a,i0)') "fine_zone_points=", count(fine_point)
    write(unit, '(a,i0)') "medium_zone_points=", count(medium_point)
    write(unit, '(a,i0)') "outer_zone_points=", count(outer_point)
    write(unit, '(a,i0)') "x_points=", mesh%x_point_count()
    write(unit, '(a,i0)') "y_points=", mesh%y_point_count()
    write(unit, '(a,es24.16e3)') "minimum_spacing=", &
      mesh%minimum_spacing()
    write(unit, '(a,es24.16e3)') "maximum_spacing=", &
      mesh%maximum_spacing()
    write(unit, '(a,es24.16e3)') "active_radius=", active_radius
    write(unit, '(a,es24.16e3)') "asymptotic_fit_inner_radius=", &
      asymptotic_fit_inner_radius
    write(unit, '(a,es24.16e3)') "asymptotic_matching_radius=", &
      asymptotic_matching_radius
    write(unit, '(a,es24.16e3)') "input_gap_embedding_relative_error=", &
      input_gap_embedding_error
    write(unit, '(a,es24.16e3)') "input_current_embedding_relative_error=", &
      input_current_embedding_error
    write(unit, '(a,es24.16e3)') "input_gap_embedding_relative_l2=", &
      input_gap_embedding_relative_l2
    write(unit, '(a,es24.16e3)') "input_current_embedding_relative_l2=", &
      input_current_embedding_relative_l2
    write(unit, '(a,es24.16e3)') "initial_gap_radial_relative_error=", &
      initial_gap_radial_error
    write(unit, '(a,es24.16e3)') "initial_current_radial_relative_error=", &
      initial_current_radial_error
    write(unit, '(a,es24.16e3)') "initial_gap_radial_relative_l2=", &
      initial_gap_radial_relative_l2
    write(unit, '(a,es24.16e3)') "initial_current_radial_relative_l2=", &
      initial_current_radial_relative_l2
    write(unit, '(a,es24.16e3)') "final_gap_radial_relative_error=", &
      final_gap_radial_error
    write(unit, '(a,es24.16e3)') "final_current_radial_relative_error=", &
      final_current_radial_error
    write(unit, '(a,es24.16e3)') "final_gap_radial_relative_l2=", &
      final_gap_radial_relative_l2
    write(unit, '(a,es24.16e3)') "final_current_radial_relative_l2=", &
      final_current_radial_relative_l2
    write(unit, '(a,es24.16e3)') &
      "final_asymptotic_halo_gap_relative_error=", &
      final_halo_gap_relative_error
    write(unit, '(a,es24.16e3)') &
      "final_asymptotic_halo_current_relative_error=", &
      final_halo_current_relative_error
    write(unit, '(a,es24.16e3)') &
      "final_asymptotic_halo_gap_relative_l2=", final_halo_gap_relative_l2
    write(unit, '(a,es24.16e3)') &
      "final_asymptotic_halo_current_relative_l2=", &
      final_halo_current_relative_l2
    write(unit, '(a,es24.16e3)') "mapped_gap_radial_relative_error=", &
      mapped_gap_radial_error
    write(unit, '(a,es24.16e3)') "mapped_current_radial_relative_error=", &
      mapped_current_radial_error
    write(unit, '(a,es24.16e3)') "mapped_gap_radial_relative_l2=", &
      mapped_gap_radial_relative_l2
    write(unit, '(a,es24.16e3)') "mapped_current_radial_relative_l2=", &
      mapped_current_radial_relative_l2
    write(unit, '(a,es24.16e3)') "map_residual_rms=", residual%rms_residual
    write(unit, '(a,es24.16e3)') "map_residual_max=", &
      residual%maximum_absolute_residual
    write(unit, '(a,es24.16e3)') "map_residual_relative_l2=", &
      residual%relative_l2_residual
    write(unit, '(a,i0)') "map_residual_worst_point=", residual%maximum_point
    call write_weighted_metrics(unit, "active", weighted_residual)
    call write_zone_metrics(unit, "fine", fine_residual)
    call write_zone_metrics(unit, "medium", medium_residual)
    call write_zone_metrics(unit, "outer", outer_residual)
    call write_weighted_metrics(unit, "fine", weighted_fine)
    call write_weighted_metrics(unit, "medium", weighted_medium)
    call write_weighted_metrics(unit, "outer", weighted_outer)
    write(unit, '(a,es24.16e3)') "maximum_normalization_error=", &
      diagnostics%global%maximum_normalization_error
    write(unit, '(a,es24.16e3)') &
      "maximum_normalization_error_over_solve=", &
      maximum_normalization_error_over_solve
    write(unit, '(a,es24.16e3)') "maximum_rank_elapsed_seconds=", &
      diagnostics%maximum_rank_elapsed_seconds
    write(unit, '(a,es24.16e3)') "cumulative_map_seconds=", &
      cumulative_map_seconds
    write(unit, '(a,l1)') "passed=", benchmark_passed
    close(unit)
  end subroutine write_metrics


  subroutine write_zone_metrics(unit, label, zone)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: label
    type(he3_field_residual_report_t), intent(in) :: zone

    write(unit, '(a,es24.16e3)') trim(label) // "_residual_rms=", &
      zone%rms_residual
    write(unit, '(a,es24.16e3)') trim(label) // "_residual_max=", &
      zone%maximum_absolute_residual
    write(unit, '(a,es24.16e3)') trim(label) // "_residual_relative_l2=", &
      zone%relative_l2_residual
  end subroutine write_zone_metrics


  subroutine write_weighted_metrics(unit, label, zone)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: label
    type(area_weighted_residual_t), intent(in) :: zone

    write(unit, '(a,es24.16e3)') &
      trim(label) // "_area_weighted_residual_rms=", zone%rms
    write(unit, '(a,es24.16e3)') &
      trim(label) // "_area_weighted_residual_relative_l2=", zone%relative_l2
    write(unit, '(a,es24.16e3)') &
      trim(label) // "_sampled_nodal_area=", zone%sampled_area
  end subroutine write_weighted_metrics


  subroutine print_metrics()
    print '(a,a)', "axisymmetric core: ", trim(core_label)
    print '(a,a)', "terminal status: ", trim(terminal_status)
    print '(a,i0,a,i0)', "iterations/maps: ", iterations_completed, "/", &
      map_evaluations
    print '(a,a)', "radial grid: ", trim(detected_grid_kind)
    print '(a,a)', "field mesh: ", trim(mesh_kind)
    print '(a,i0,a,i0,a,2(es10.3,1x))', "points (active/total): ", &
      count(active_point), "/", mesh%point_count(), " spacing (min,max): ", &
      mesh%minimum_spacing(), mesh%maximum_spacing()
    print '(a,2(es12.4,1x))', "input embedding errors (gap,current): ", &
      input_gap_embedding_error, input_current_embedding_error
    print '(a,2(es12.4,1x))', "mapped radial errors (gap,current): ", &
      mapped_gap_radial_error, mapped_current_radial_error
    print '(a,2(es12.4,1x))', "final radial errors (gap,current): ", &
      final_gap_radial_error, final_current_radial_error
    print '(a,2(es12.4,1x))', "final radial relative-L2 (gap,current): ", &
      final_gap_radial_relative_l2, final_current_radial_relative_l2
    print '(a,3(es12.4,1x))', "map residual (rms,max,relative-L2): ", &
      residual%rms_residual, residual%maximum_absolute_residual, &
      residual%relative_l2_residual
    print '(a,3(es12.4,1x))', "zone relative-L2 (fine,medium,outer): ", &
      fine_residual%relative_l2_residual, &
      medium_residual%relative_l2_residual, &
      outer_residual%relative_l2_residual
    print '(a,4(es12.4,1x))', &
      "area-weighted relative-L2 (all,fine,medium,outer): ", &
      weighted_residual%relative_l2, weighted_fine%relative_l2, &
      weighted_medium%relative_l2, weighted_outer%relative_l2
    print '(a,es12.4)', "maximum normalization error over solve: ", &
      maximum_normalization_error_over_solve
    print '(a,f11.3)', "cumulative map time [s]: ", &
      cumulative_map_seconds
    print '(a,l1)', "benchmark passed: ", benchmark_passed
  end subroutine print_metrics


  subroutine require_mpi(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require_mpi

end program benchmark_axisymmetric_core_2d
