program test_he3_mpi_field_map
  use mpi_f08
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, &
    make_uniform_cartesian_mesh, make_circular_active_mask
  use spinful_state_2d, only : spinful_state_2d_t, allocate_spinful_state_2d
  use he3_quadrature, only : angular_quadrature_3d_t, ozaki_quadrature_t, &
                             read_legacy_gauss_table, &
                             read_legacy_ozaki_table
  use he3_field_map_2d, only : he3_field_map_diagnostics_t, &
                               evaluate_serial_he3_field_map
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
  type(he3_field_map_diagnostics_t) :: serial_diagnostics
  type(he3_mpi_field_map_diagnostics_t) :: mpi_diagnostics
  type(he3_mpi_field_map_diagnostics_t) :: masked_diagnostics
  type(he3_field_residual_report_t) :: residual_report
  type(free_vortex_endpoint_2d_t) :: endpoint
  type(legacy_anderson_t) :: accelerator
  type(ozaki_quadrature_t) :: energy
  type(spinful_state_2d_t) :: masked_mapped, mpi_mapped
  type(spinful_state_2d_t) :: next, serial_mapped, state
  character(len=1024) :: gauss_path, ozaki_path
  real(rk) :: checksum, x, y
  integer :: ierr, passed_integer, point, rank, rank_count, spin, orbital
  integer :: x_node, y_node
  logical, allocatable :: active_point(:)
  logical :: masked_succeeded, mpi_succeeded, passed, serial_succeeded

  call MPI_Init(ierr)
  call require_mpi(ierr == MPI_SUCCESS, "MPI initialization failed")
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call require_mpi(ierr == MPI_SUCCESS, "MPI rank query failed")
  call MPI_Comm_size(MPI_COMM_WORLD, rank_count, ierr)
  call require_mpi(ierr == MPI_SUCCESS, "MPI size query failed")

  if (command_argument_count() /= 2) &
    error stop "MPI field-map test needs Gauss and Ozaki table paths"
  call get_command_argument(1, gauss_path)
  call get_command_argument(2, ozaki_path)
  call read_legacy_gauss_table(trim(gauss_path), 48, 11, angular)
  call read_legacy_ozaki_table(trim(ozaki_path), energy)
  call make_uniform_cartesian_mesh( &
    -0.4_rk, 0.4_rk, 3, -0.4_rk, 0.4_rk, 3, mesh)
  call allocate_spinful_state_2d(mesh, state)
  do y_node = 0, mesh%y_cell_count()
    y = mesh%y_coordinate(y_node)
    do x_node = 0, mesh%x_cell_count()
      x = mesh%x_coordinate(x_node)
      point = mesh%point_index(x_node, y_node)
      do spin = 1, 3
        do orbital = 1, 3
          state%order_parameter(spin, orbital, point) = cmplx( &
            0.08_rk * merge(1.0_rk, 0.2_rk, spin == orbital) + &
              0.006_rk * x - 0.004_rk * y, &
            0.003_rk * real(spin - orbital, rk) + 0.002_rk * (x + y), rk)
        end do
      end do
      state%current_mean_field(:, point) = &
        [0.01_rk * x, -0.012_rk * y, 0.004_rk * (x - y)]
    end do
  end do

  endpoint%enabled = .true.
  endpoint%bulk_gap = 0.1_rk
  endpoint%winding = 1.0_rk
  endpoint%outer_radius = 1.0_rk
  endpoint%maximum_step = 0.2_rk

  call evaluate_mpi_he3_field_map( &
    MPI_COMM_WORLD, mesh, state, angular, energy, 0.7_rk, 0.2_rk, 1, &
    0.2_rk, mpi_mapped, mpi_diagnostics, mpi_succeeded, endpoint)
  call make_circular_active_mask( &
    mesh, [0.0_rk, 0.0_rk], 0.2_rk, active_point)
  call evaluate_mpi_he3_field_map( &
    MPI_COMM_WORLD, mesh, state, angular, energy, 0.7_rk, 0.2_rk, 1, &
    0.2_rk, masked_mapped, masked_diagnostics, masked_succeeded, &
    free_vortex_endpoint=endpoint, active_point_mask=active_point)
  if (rank == 0) call accelerator%initialize( &
    new_src_iteration_vector_size_2d(state), 5, 0.1_rk, 5.0_rk)
  call update_mpi_he3_state_with_anderson( &
    MPI_COMM_WORLD, mesh, accelerator, state, mpi_mapped, 0.0_rk, next, &
    anderson_report, residual_report)
  passed = mpi_succeeded
  passed = passed .and. mpi_diagnostics%rank_count == rank_count
  passed = passed .and. &
    mpi_diagnostics%global%requested_point_count == mesh%point_count()
  passed = passed .and. &
    mpi_diagnostics%global%completed_point_count == mesh%point_count()
  passed = passed .and. mpi_diagnostics%global%attempted_contributions == &
    mesh%point_count() * angular%direction_count() * energy%pole_count()
  passed = passed .and. masked_succeeded
  passed = passed .and. masked_diagnostics%global%requested_point_count == &
    count(active_point)
  passed = passed .and. masked_diagnostics%global%attempted_contributions == &
    count(active_point) * angular%direction_count() * energy%pole_count()
  do point = 1, mesh%point_count()
    if (active_point(point)) cycle
    passed = passed .and. maxval(abs( &
      masked_mapped%order_parameter(:, :, point) - &
      state%order_parameter(:, :, point))) < epsilon(1.0_rk)
    passed = passed .and. maxval(abs( &
      masked_mapped%current_mean_field(:, point) - &
      state%current_mean_field(:, point))) < epsilon(1.0_rk)
  end do
  passed = passed .and. &
    mpi_diagnostics%maximum_points_per_rank - &
    mpi_diagnostics%minimum_points_per_rank <= 1
  passed = passed .and. anderson_report%iteration == 1 .and. &
    anderson_report%history_size == 1
  passed = passed .and. residual_report%value_count == &
    21 * mesh%point_count()
  passed = passed .and. maxval(abs(next%order_parameter - &
    (state%order_parameter + cmplx(0.01_rk, 0.0_rk, rk) * &
     (mpi_mapped%order_parameter - state%order_parameter)))) < 3.0e-15_rk
  passed = passed .and. maxval(abs(next%current_mean_field - &
    (state%current_mean_field + 0.01_rk * &
     (mpi_mapped%current_mean_field - state%current_mean_field)))) < 3.0e-15_rk

  if (rank == 0) then
    call evaluate_serial_he3_field_map( &
      mesh, state, angular, energy, 0.7_rk, 0.2_rk, 1, 0.2_rk, &
      serial_mapped, serial_diagnostics, serial_succeeded, endpoint)
    passed = passed .and. serial_succeeded
    passed = passed .and. &
      maxval(abs(mpi_mapped%order_parameter - &
                 serial_mapped%order_parameter)) < 3.0e-14_rk
    passed = passed .and. &
      maxval(abs(mpi_mapped%current_mean_field - &
                 serial_mapped%current_mean_field)) < 3.0e-15_rk
    checksum = sum(real(mpi_mapped%order_parameter, rk)) + &
      sum(aimag(mpi_mapped%order_parameter)) + &
      sum(mpi_mapped%current_mean_field)
    print '(a,i0)', "MPI ranks: ", rank_count
    print '(a,i0,a,i0)', "points per rank: ", &
      mpi_diagnostics%minimum_points_per_rank, " to ", &
      mpi_diagnostics%maximum_points_per_rank
    print '(a,es25.17)', "field checksum: ", checksum
    print '(a,es12.4)', "maximum normalization error: ", &
      mpi_diagnostics%global%maximum_normalization_error
    print '(a,f8.3)', "maximum rank elapsed time [s]: ", &
      mpi_diagnostics%maximum_rank_elapsed_seconds
  end if
  passed_integer = merge(1, 0, passed)
  call MPI_Bcast(passed_integer, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call require_mpi(ierr == MPI_SUCCESS, "MPI result broadcast failed")
  passed = passed_integer == 1

  call MPI_Finalize(ierr)
  call require_mpi(ierr == MPI_SUCCESS, "MPI finalization failed")
  call require_mpi(passed, "MPI field map differs from serial field map")

contains

  subroutine require_mpi(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require_mpi

end program test_he3_mpi_field_map
