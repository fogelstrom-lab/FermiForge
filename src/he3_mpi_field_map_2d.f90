module he3_mpi_field_map_2d
  use mpi_f08
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t
  use spinful_state_2d, only : spinful_state_2d_t, allocate_spinful_state_2d
  use he3_quadrature, only : angular_quadrature_3d_t, ozaki_quadrature_t
  use free_vortex_asymptotic_2d, only : free_vortex_endpoint_2d_t
  use he3_field_map_2d, only : he3_field_map_diagnostics_t, &
                               evaluate_he3_field_map_subset
  use legacy_anderson_mixing, only : legacy_anderson_t, anderson_report_t
  use he3_nonlinear_iteration_2d, only : he3_field_residual_report_t, &
                                         update_he3_state_with_anderson
  implicit none
  private

  type, public :: he3_mpi_field_map_diagnostics_t
    integer :: rank_count = 0
    integer :: minimum_points_per_rank = 0
    integer :: maximum_points_per_rank = 0
    real(rk) :: maximum_rank_elapsed_seconds = 0.0_rk
    type(he3_field_map_diagnostics_t) :: global
  end type he3_mpi_field_map_diagnostics_t

  public :: evaluate_mpi_he3_field_map
  public :: update_mpi_he3_state_with_anderson

contains

  subroutine evaluate_mpi_he3_field_map( &
      communicator, mesh, state, angular_quadrature, energy_quadrature, &
      feedback_scale, trajectory_maximum_step, minimum_substep_count, &
      boundary_relaxation_distance, mapped_state, diagnostics, succeeded, &
      free_vortex_endpoint, active_point_mask)
    type(MPI_Comm), intent(in) :: communicator
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    type(angular_quadrature_3d_t), intent(in) :: angular_quadrature
    type(ozaki_quadrature_t), intent(in) :: energy_quadrature
    real(rk), intent(in) :: feedback_scale, trajectory_maximum_step
    integer, intent(in) :: minimum_substep_count
    real(rk), intent(in) :: boundary_relaxation_distance
    type(spinful_state_2d_t), intent(out) :: mapped_state
    type(he3_mpi_field_map_diagnostics_t), intent(out) :: diagnostics
    logical, intent(out) :: succeeded
    type(free_vortex_endpoint_2d_t), intent(in), optional :: &
      free_vortex_endpoint
    logical, intent(in), optional :: active_point_mask(:)

    type(he3_field_map_diagnostics_t) :: local_diagnostics
    type(spinful_state_2d_t) :: local_state
    integer, allocatable :: global_completion(:), local_completion(:)
    integer, allocatable :: local_points(:)
    logical, allocatable :: active(:)
    integer :: global_counts(6), ierr, local_count, local_counts(6)
    integer :: active_ordinal, local_index, local_point_minimum
    integer :: point, rank, rank_count
    real(rk) :: global_times(3), local_elapsed, local_times(3)
    real(rk) :: start_time
    logical, allocatable :: local_point_completed(:)
    logical :: local_succeeded

    call MPI_Comm_rank(communicator, rank, ierr)
    call require_mpi_success(ierr, "could not obtain MPI rank")
    call MPI_Comm_size(communicator, rank_count, ierr)
    call require_mpi_success(ierr, "could not obtain MPI size")

    allocate(active(mesh%point_count()), source=.true.)
    if (present(active_point_mask)) then
      if (size(active_point_mask) /= mesh%point_count()) &
        error stop "MPI field-map active mask has the wrong size"
      active = active_point_mask
      if (.not. any(active)) &
        error stop "MPI field-map active mask contains no points"
    end if

    local_count = 0
    active_ordinal = 0
    do point = 1, mesh%point_count()
      if (.not. active(point)) cycle
      active_ordinal = active_ordinal + 1
      if (modulo(active_ordinal - 1, rank_count) == rank) &
        local_count = local_count + 1
    end do
    allocate(local_points(local_count))
    active_ordinal = 0
    local_index = 0
    do point = 1, mesh%point_count()
      if (.not. active(point)) cycle
      active_ordinal = active_ordinal + 1
      if (modulo(active_ordinal - 1, rank_count) /= rank) cycle
      local_index = local_index + 1
      local_points(local_index) = point
    end do

    start_time = real(MPI_Wtime(), rk)
    call evaluate_he3_field_map_subset( &
      mesh, state, local_points, angular_quadrature, energy_quadrature, &
      feedback_scale, trajectory_maximum_step, minimum_substep_count, &
      boundary_relaxation_distance, local_state, local_point_completed, &
      local_diagnostics, local_succeeded, free_vortex_endpoint)
    local_elapsed = real(MPI_Wtime(), rk) - start_time

    call allocate_spinful_state_2d(mesh, mapped_state)
    call MPI_Allreduce( &
      local_state%order_parameter, mapped_state%order_parameter, &
      size(mapped_state%order_parameter), MPI_DOUBLE_COMPLEX, MPI_SUM, &
      communicator, ierr)
    call require_mpi_success(ierr, "could not reduce mapped order parameter")
    call MPI_Allreduce( &
      local_state%current_mean_field, mapped_state%current_mean_field, &
      size(mapped_state%current_mean_field), MPI_DOUBLE_PRECISION, MPI_SUM, &
      communicator, ierr)
    call require_mpi_success(ierr, "could not reduce mapped mean field")
    do point = 1, mesh%point_count()
      if (active(point)) cycle
      mapped_state%order_parameter(:, :, point) = &
        state%order_parameter(:, :, point)
      mapped_state%current_mean_field(:, point) = &
        state%current_mean_field(:, point)
    end do

    allocate(local_completion(mesh%point_count()), &
             global_completion(mesh%point_count()))
    local_completion = merge(1, 0, local_point_completed)
    call MPI_Allreduce(local_completion, global_completion, &
      mesh%point_count(), MPI_INTEGER, MPI_SUM, communicator, ierr)
    call require_mpi_success(ierr, "could not reduce point completion mask")

    local_counts = [local_diagnostics%requested_point_count, &
      local_diagnostics%completed_point_count, &
      local_diagnostics%failed_point_count, &
      local_diagnostics%attempted_contributions, &
      local_diagnostics%accepted_contributions, &
      local_diagnostics%failed_propagators]
    call MPI_Allreduce( &
      local_counts, global_counts, size(local_counts), MPI_INTEGER, MPI_SUM, &
      communicator, ierr)
    call require_mpi_success(ierr, "could not reduce field-map counts")
    local_times = [local_diagnostics%sampling_seconds, &
      local_diagnostics%propagation_seconds, &
      local_diagnostics%accumulation_seconds]
    call MPI_Allreduce( &
      local_times, global_times, size(local_times), MPI_DOUBLE_PRECISION, &
      MPI_SUM, communicator, ierr)
    call require_mpi_success(ierr, "could not reduce field-map timings")

    diagnostics = he3_mpi_field_map_diagnostics_t()
    diagnostics%rank_count = rank_count
    call MPI_Allreduce(local_count, diagnostics%maximum_points_per_rank, 1, &
      MPI_INTEGER, MPI_MAX, communicator, ierr)
    call require_mpi_success(ierr, "could not reduce maximum rank work")
    local_point_minimum = local_count
    call MPI_Allreduce(local_point_minimum, &
      diagnostics%minimum_points_per_rank, 1, MPI_INTEGER, MPI_MIN, &
      communicator, ierr)
    call require_mpi_success(ierr, "could not reduce minimum rank work")
    call MPI_Allreduce(local_elapsed, &
      diagnostics%maximum_rank_elapsed_seconds, 1, MPI_DOUBLE_PRECISION, &
      MPI_MAX, communicator, ierr)
    call require_mpi_success(ierr, "could not reduce maximum rank time")

    diagnostics%global%requested_point_count = global_counts(1)
    diagnostics%global%completed_point_count = global_counts(2)
    diagnostics%global%failed_point_count = global_counts(3)
    diagnostics%global%direction_count = angular_quadrature%direction_count()
    diagnostics%global%pole_count = energy_quadrature%pole_count()
    diagnostics%global%attempted_contributions = global_counts(4)
    diagnostics%global%accepted_contributions = global_counts(5)
    diagnostics%global%failed_propagators = global_counts(6)
    diagnostics%global%sampling_seconds = global_times(1)
    diagnostics%global%propagation_seconds = global_times(2)
    diagnostics%global%accumulation_seconds = global_times(3)
    call MPI_Allreduce(local_diagnostics%maximum_normalization_error, &
      diagnostics%global%maximum_normalization_error, 1, &
      MPI_DOUBLE_PRECISION, MPI_MAX, communicator, ierr)
    call require_mpi_success(ierr, "could not reduce normalization error")

    succeeded = local_succeeded .and. &
      diagnostics%global%failed_point_count == 0 .and. &
      diagnostics%global%failed_propagators == 0 .and. &
      all((global_completion == 1) .eqv. active)
  end subroutine evaluate_mpi_he3_field_map


  subroutine update_mpi_he3_state_with_anderson( &
      communicator, mesh, accelerator, current, mapped, tolerance, next, &
      anderson_report, residual_report)
    type(MPI_Comm), intent(in) :: communicator
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(legacy_anderson_t), intent(inout) :: accelerator
    type(spinful_state_2d_t), intent(in) :: current, mapped
    real(rk), intent(in) :: tolerance
    type(spinful_state_2d_t), intent(out) :: next
    type(anderson_report_t), intent(out) :: anderson_report
    type(he3_field_residual_report_t), intent(out) :: residual_report

    integer :: anderson_integers(4), ierr, rank, residual_integers(2)
    real(rk) :: anderson_reals(4), residual_reals(6)

    call MPI_Comm_rank(communicator, rank, ierr)
    call require_mpi_success(ierr, "could not obtain Anderson MPI rank")
    if (rank == 0) then
      call update_he3_state_with_anderson( &
        mesh, accelerator, current, mapped, tolerance, next, &
        anderson_report, residual_report)
      anderson_integers = [anderson_report%iteration, &
        anderson_report%history_size, anderson_report%discarded_vectors, &
        merge(1, 0, anderson_report%converged) + &
        2 * merge(1, 0, anderson_report%used_stochastic_step)]
      anderson_reals = [anderson_report%max_residual, &
        anderson_report%residual_norm, &
        anderson_report%legacy_relative_residual, anderson_report%mixing]
      residual_integers = [residual_report%value_count, &
                           residual_report%maximum_point]
      residual_reals = [residual_report%mean_absolute_residual, &
        residual_report%rms_residual, &
        residual_report%maximum_absolute_residual, &
        residual_report%relative_l2_residual, &
        residual_report%legacy_scaled_maximum_residual, &
        residual_report%maximum_point_rms_residual]
    else
      call allocate_spinful_state_2d(mesh, next)
      anderson_integers = 0
      anderson_reals = 0.0_rk
      residual_integers = 0
      residual_reals = 0.0_rk
    end if

    call MPI_Bcast(next%order_parameter, size(next%order_parameter), &
      MPI_DOUBLE_COMPLEX, 0, communicator, ierr)
    call require_mpi_success(ierr, "could not broadcast Anderson gap state")
    call MPI_Bcast(next%current_mean_field, size(next%current_mean_field), &
      MPI_DOUBLE_PRECISION, 0, communicator, ierr)
    call require_mpi_success(ierr, &
      "could not broadcast Anderson current-field state")
    call MPI_Bcast(anderson_integers, size(anderson_integers), MPI_INTEGER, &
      0, communicator, ierr)
    call require_mpi_success(ierr, "could not broadcast Anderson report")
    call MPI_Bcast(anderson_reals, size(anderson_reals), MPI_DOUBLE_PRECISION, &
      0, communicator, ierr)
    call require_mpi_success(ierr, "could not broadcast Anderson metrics")
    call MPI_Bcast(residual_integers, size(residual_integers), MPI_INTEGER, &
      0, communicator, ierr)
    call require_mpi_success(ierr, "could not broadcast residual report")
    call MPI_Bcast(residual_reals, size(residual_reals), &
      MPI_DOUBLE_PRECISION, 0, communicator, ierr)
    call require_mpi_success(ierr, "could not broadcast residual metrics")

    if (rank /= 0) then
      anderson_report = anderson_report_t()
      anderson_report%iteration = anderson_integers(1)
      anderson_report%history_size = anderson_integers(2)
      anderson_report%discarded_vectors = anderson_integers(3)
      anderson_report%converged = modulo(anderson_integers(4), 2) == 1
      anderson_report%used_stochastic_step = anderson_integers(4) >= 2
      anderson_report%max_residual = anderson_reals(1)
      anderson_report%residual_norm = anderson_reals(2)
      anderson_report%legacy_relative_residual = anderson_reals(3)
      anderson_report%mixing = anderson_reals(4)
      residual_report = he3_field_residual_report_t()
      residual_report%value_count = residual_integers(1)
      residual_report%maximum_point = residual_integers(2)
      residual_report%mean_absolute_residual = residual_reals(1)
      residual_report%rms_residual = residual_reals(2)
      residual_report%maximum_absolute_residual = residual_reals(3)
      residual_report%relative_l2_residual = residual_reals(4)
      residual_report%legacy_scaled_maximum_residual = residual_reals(5)
      residual_report%maximum_point_rms_residual = residual_reals(6)
    end if
  end subroutine update_mpi_he3_state_with_anderson


  subroutine require_mpi_success(error_code, message)
    integer, intent(in) :: error_code
    character(len=*), intent(in) :: message

    if (error_code /= MPI_SUCCESS) error stop message
  end subroutine require_mpi_success

end module he3_mpi_field_map_2d
