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
  use spinful_field_io_2d, only : write_spinful_field_map_2d
  use he3_quadrature, only : angular_quadrature_3d_t, ozaki_quadrature_t, &
                             read_legacy_gauss_table, &
                             read_legacy_ozaki_table
  use he3_mpi_field_map_2d, only : he3_mpi_field_map_diagnostics_t, &
                                   evaluate_mpi_he3_field_map
  use he3_nonlinear_iteration_2d, only : he3_field_residual_report_t, &
                                         compute_he3_field_residual
  use free_vortex_asymptotic_2d, only : free_vortex_endpoint_2d_t, &
    set_axisymmetric_radial_reference
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
  type(ozaki_quadrature_t) :: energy
  type(axial_radial_profile_t) :: profile
  type(radial_mesh_t) :: radial_grid
  type(spinful_state_2d_t) :: mapped, state

  character(len=16) :: detected_grid_kind, endpoint_policy, mesh_kind
  character(len=16) :: radial_grid_kind
  character(len=64) :: core_label
  character(len=512) :: current_field_file, gauss_file, input_file
  character(len=512) :: input_message, input_output_file, mapped_output_file
  character(len=512) :: metrics_file, order_parameter_file, output_override
  character(len=512) :: ozaki_file
  integer :: ierr, input_status, input_unit, minimum_substep_count
  integer :: number_of_cells, rank, rank_count, result_integer
  real(rk) :: asymptotic_bulk_gap, asymptotic_center_x, asymptotic_center_y
  real(rk) :: asymptotic_matching_radius, asymptotic_maximum_step
  real(rk) :: asymptotic_outer_radius
  real(rk) :: asymptotic_winding, boundary_relaxation_distance
  real(rk) :: feedback_scale, half_width, input_current_embedding_error
  real(rk) :: input_gap_embedding_error, mapped_current_radial_error
  real(rk) :: mapped_gap_radial_error, maximum_embedding_error
  real(rk) :: active_radius, coarse_spacing, fine_region_half_width
  real(rk) :: fine_spacing, medium_region_half_width, medium_spacing
  real(rk) :: trajectory_maximum_step
  logical, allocatable :: active_point(:), fine_point(:)
  logical, allocatable :: medium_point(:), outer_point(:)
  logical :: benchmark_passed, map_succeeded

  namelist /axisymmetric_core_benchmark/ core_label, half_width, &
    number_of_cells, mesh_kind, fine_region_half_width, &
    medium_region_half_width, fine_spacing, medium_spacing, coarse_spacing, &
    active_radius, feedback_scale, trajectory_maximum_step, &
    minimum_substep_count, boundary_relaxation_distance, &
    maximum_embedding_error, order_parameter_file, current_field_file, &
    radial_grid_kind, gauss_file, ozaki_file, endpoint_policy, &
    asymptotic_center_x, asymptotic_center_y, asymptotic_bulk_gap, &
    asymptotic_winding, asymptotic_matching_radius, &
    asymptotic_outer_radius, asymptotic_maximum_step, &
    input_output_file, mapped_output_file, metrics_file

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
  asymptotic_matching_radius = 0.0_rk
  asymptotic_outer_radius = 70.0_rk
  asymptotic_maximum_step = 0.25_rk
  order_parameter_file = "benchmarks/normal_core_2d/reference/op_xyz"
  current_field_file = "benchmarks/normal_core_2d/reference/curr"
  gauss_file = "new_src/gauss11.dat"
  ozaki_file = "new_src/ozaki.dat"
  input_output_file = "work/normal-core-2d/input_fields_2d.dat"
  mapped_output_file = "work/normal-core-2d/mapped_fields_2d.dat"
  metrics_file = "work/normal-core-2d/metrics.txt"

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
  call allocate_spinful_state_2d(mesh, state)
  call embed_axial_radial_profile_2d( &
    radial_grid, profile, asymptotic_winding, mesh, state)

  call evaluate_mpi_he3_field_map( &
    MPI_COMM_WORLD, mesh, state, angular, energy, feedback_scale, &
    trajectory_maximum_step, minimum_substep_count, &
    boundary_relaxation_distance, mapped, diagnostics, map_succeeded, &
    free_vortex_endpoint=endpoint, active_point_mask=active_point)

  benchmark_passed = map_succeeded
  if (rank == 0) then
    call radial_reference_error( &
      mesh, radial_grid, profile, state, asymptotic_winding, &
      input_gap_embedding_error, input_current_embedding_error)
    call radial_reference_error( &
      mesh, radial_grid, profile, mapped, asymptotic_winding, &
      mapped_gap_radial_error, mapped_current_radial_error)
    call compute_he3_field_residual( &
      state, mapped, residual, active_point_mask=active_point)
    call compute_he3_field_residual( &
      state, mapped, fine_residual, active_point_mask=fine_point)
    call compute_he3_field_residual( &
      state, mapped, medium_residual, active_point_mask=medium_point)
    call compute_he3_field_residual( &
      state, mapped, outer_residual, active_point_mask=outer_point)
    call compute_area_weighted_residual( &
      mesh, state, mapped, active_point, weighted_residual)
    call compute_area_weighted_residual( &
      mesh, state, mapped, fine_point, weighted_fine)
    call compute_area_weighted_residual( &
      mesh, state, mapped, medium_point, weighted_medium)
    call compute_area_weighted_residual( &
      mesh, state, mapped, outer_point, weighted_outer)
    benchmark_passed = benchmark_passed .and. &
      input_gap_embedding_error <= maximum_embedding_error .and. &
      input_current_embedding_error <= maximum_embedding_error .and. &
      diagnostics%global%maximum_normalization_error <= 1.0e-10_rk
    call write_spinful_field_map_2d( &
      trim(input_output_file), mesh, state, &
      trim(core_label) // "_embedded_radial_reference", trim(endpoint_policy))
    call write_spinful_field_map_2d( &
      trim(mapped_output_file), mesh, mapped, &
      trim(core_label) // "_quasiclassical_point_map", trim(endpoint_policy))
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
    case default
      error stop &
        "endpoint_policy must be 'local', 'free_vortex', or 'radial_reference'"
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
    if (endpoint%enabled .and. &
        (asymptotic_bulk_gap <= 0.0_rk .or. &
         asymptotic_outer_radius <= sqrt(2.0_rk) * half_width .or. &
         asymptotic_maximum_step <= 0.0_rk)) &
      error stop "invalid free-vortex benchmark endpoint controls"
  end subroutine validate_input


  subroutine radial_reference_error( &
      field_mesh, source_mesh, source_profile, field, winding, gap_error, &
      current_error)
    type(cartesian_mesh_2d_t), intent(in) :: field_mesh
    type(radial_mesh_t), intent(in) :: source_mesh
    type(axial_radial_profile_t), intent(in) :: source_profile
    type(spinful_state_2d_t), intent(in) :: field
    real(rk), intent(in) :: winding
    real(rk), intent(out) :: gap_error, current_error

    complex(rk) :: expected_gap(3, 3)
    real(rk) :: expected_current(3), gap_scale, current_scale, x, y
    integer :: point, x_node, y_node
    logical :: inside

    gap_error = 0.0_rk
    current_error = 0.0_rk
    gap_scale = 0.0_rk
    current_scale = 0.0_rk
    do y_node = 0, field_mesh%y_cell_count()
      y = field_mesh%y_coordinate(y_node)
      do x_node = 0, field_mesh%x_cell_count()
        x = field_mesh%x_coordinate(x_node)
        point = field_mesh%point_index(x_node, y_node)
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
      end do
    end do
    gap_error = gap_error / max(gap_scale, tiny(1.0_rk))
    current_error = current_error / max(current_scale, tiny(1.0_rk))
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
    write(unit, '(a,i0)') "mpi_ranks=", rank_count
    write(unit, '(a,i0)') "cartesian_points=", mesh%point_count()
    write(unit, '(a,i0)') "active_points=", count(active_point)
    write(unit, '(a,i0)') "frozen_halo_points=", &
      mesh%point_count() - count(active_point)
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
    write(unit, '(a,es24.16e3)') "input_gap_embedding_relative_error=", &
      input_gap_embedding_error
    write(unit, '(a,es24.16e3)') "input_current_embedding_relative_error=", &
      input_current_embedding_error
    write(unit, '(a,es24.16e3)') "mapped_gap_radial_relative_error=", &
      mapped_gap_radial_error
    write(unit, '(a,es24.16e3)') "mapped_current_radial_relative_error=", &
      mapped_current_radial_error
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
    write(unit, '(a,es24.16e3)') "maximum_rank_elapsed_seconds=", &
      diagnostics%maximum_rank_elapsed_seconds
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
    print '(a,a)', "radial grid: ", trim(detected_grid_kind)
    print '(a,a)', "field mesh: ", trim(mesh_kind)
    print '(a,i0,a,i0,a,2(es10.3,1x))', "points (active/total): ", &
      count(active_point), "/", mesh%point_count(), " spacing (min,max): ", &
      mesh%minimum_spacing(), mesh%maximum_spacing()
    print '(a,2(es12.4,1x))', "input embedding errors (gap,current): ", &
      input_gap_embedding_error, input_current_embedding_error
    print '(a,2(es12.4,1x))', "mapped radial errors (gap,current): ", &
      mapped_gap_radial_error, mapped_current_radial_error
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
    print '(a,es12.4)', "maximum normalization error: ", &
      diagnostics%global%maximum_normalization_error
    print '(a,f9.3)', "slowest rank [s]: ", &
      diagnostics%maximum_rank_elapsed_seconds
    print '(a,l1)', "benchmark passed: ", benchmark_passed
  end subroutine print_metrics


  subroutine require_mpi(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require_mpi

end program benchmark_axisymmetric_core_2d
