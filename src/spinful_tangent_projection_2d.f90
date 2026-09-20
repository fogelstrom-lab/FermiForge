module spinful_tangent_projection_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t
  use spinful_state_2d, only : spinful_state_2d_t
  use new_src_iteration_layout_2d, only : &
    new_src_masked_iteration_vector_size_2d, &
    pack_new_src_masked_iteration_vector_2d, &
    unpack_new_src_masked_iteration_vector_2d
  implicit none
  private

  integer, parameter, public :: zero_mode_gauge = 1
  integer, parameter, public :: zero_mode_translation_x = 2
  integer, parameter, public :: zero_mode_translation_y = 3
  integer, parameter, public :: zero_mode_orientation = 4
  integer, parameter, public :: zero_mode_count = 4

  type, public :: zero_mode_tangent_basis_2d_t
    integer :: requested_count = 0
    integer :: retained_count = 0
    integer, allocatable :: mode_id(:)
    real(rk), allocatable :: vector(:, :)
    real(rk) :: raw_norm(zero_mode_count) = 0.0_rk
  end type zero_mode_tangent_basis_2d_t

  type, public :: zero_mode_projection_report_t
    integer :: mode_count = 0
    real(rk) :: original_residual_norm = 0.0_rk
    real(rk) :: projected_residual_norm = 0.0_rk
    real(rk) :: removed_residual_norm = 0.0_rk
    real(rk) :: removed_fraction = 0.0_rk
    real(rk) :: coefficient(zero_mode_count) = 0.0_rk
  end type zero_mode_projection_report_t

  public :: build_zero_mode_tangent_basis_2d
  public :: project_spinful_residual_away_from_tangents_2d
  public :: zero_mode_name

contains

  subroutine build_zero_mode_tangent_basis_2d( &
      mesh, reference, active_point_mask, include_gauge, &
      include_translations, include_orientation, basis)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: reference
    logical, intent(in) :: active_point_mask(:)
    logical, intent(in) :: include_gauge, include_translations
    logical, intent(in) :: include_orientation
    type(zero_mode_tangent_basis_2d_t), intent(out) :: basis

    integer :: candidate_count, candidate_index, mode, pass, retained
    integer, allocatable :: candidate_mode(:), retained_mode(:)
    real(rk) :: candidate_norm, rejection_threshold
    real(rk), allocatable :: candidate(:, :), orthonormal(:, :), work(:)

    call require_valid_arguments(mesh, reference, active_point_mask)
    candidate_count = merge(1, 0, include_gauge) + &
      2 * merge(1, 0, include_translations) + &
      merge(1, 0, include_orientation)
    basis = zero_mode_tangent_basis_2d_t()
    basis%requested_count = candidate_count
    allocate(basis%mode_id(0), basis%vector( &
      new_src_masked_iteration_vector_size_2d( &
        reference, active_point_mask), 0))
    if (candidate_count == 0) return

    allocate(candidate_mode(candidate_count))
    candidate_index = 0
    if (include_gauge) then
      candidate_index = candidate_index + 1
      candidate_mode(candidate_index) = zero_mode_gauge
    end if
    if (include_translations) then
      candidate_index = candidate_index + 1
      candidate_mode(candidate_index) = zero_mode_translation_x
      candidate_index = candidate_index + 1
      candidate_mode(candidate_index) = zero_mode_translation_y
    end if
    if (include_orientation) then
      candidate_index = candidate_index + 1
      candidate_mode(candidate_index) = zero_mode_orientation
    end if

    allocate(candidate(size(basis%vector, 1), candidate_count))
    allocate(orthonormal(size(basis%vector, 1), candidate_count), &
             source=0.0_rk)
    allocate(retained_mode(candidate_count), source=0)
    allocate(work(size(basis%vector, 1)))
    do candidate_index = 1, candidate_count
      mode = candidate_mode(candidate_index)
      call build_mode_vector( &
        mesh, reference, active_point_mask, mode, &
        candidate(:, candidate_index))
      basis%raw_norm(mode) = sqrt(max(dot_product( &
        candidate(:, candidate_index), candidate(:, candidate_index)), &
        0.0_rk))
    end do
    rejection_threshold = 4096.0_rk * epsilon(1.0_rk) * &
      max(1.0_rk, maxval(basis%raw_norm))

    retained = 0
    do candidate_index = 1, candidate_count
      work = candidate(:, candidate_index)
      ! Reorthogonalize once. This keeps the small collective-coordinate
      ! basis stable when orientation contains a substantial gauge part.
      do pass = 1, 2
        do mode = 1, retained
          work = work - dot_product(orthonormal(:, mode), work) * &
            orthonormal(:, mode)
        end do
      end do
      candidate_norm = sqrt(max(dot_product(work, work), 0.0_rk))
      if (candidate_norm <= rejection_threshold) cycle
      retained = retained + 1
      orthonormal(:, retained) = work / candidate_norm
      retained_mode(retained) = candidate_mode(candidate_index)
    end do
    basis%retained_count = retained
    deallocate(basis%mode_id, basis%vector)
    allocate(basis%mode_id(retained))
    allocate(basis%vector(size(candidate, 1), retained))
    if (retained > 0) then
      basis%mode_id = retained_mode(:retained)
      basis%vector = orthonormal(:, :retained)
    end if
  end subroutine build_zero_mode_tangent_basis_2d


  subroutine project_spinful_residual_away_from_tangents_2d( &
      current, mapped, active_point_mask, basis, projected_mapped, report)
    type(spinful_state_2d_t), intent(in) :: current, mapped
    logical, intent(in) :: active_point_mask(:)
    type(zero_mode_tangent_basis_2d_t), intent(in) :: basis
    type(spinful_state_2d_t), intent(out) :: projected_mapped
    type(zero_mode_projection_report_t), intent(out) :: report

    integer :: mode
    real(rk) :: original_squared, projected_squared
    real(rk), allocatable :: current_vector(:), mapped_vector(:), residual(:)

    if (current%point_count() /= mapped%point_count() .or. &
        current%point_count() < 1) &
      error stop "tangent projection states have incompatible point counts"
    if (size(active_point_mask) /= current%point_count() .or. &
        .not. any(active_point_mask)) &
      error stop "tangent projection mask is invalid"
    call pack_new_src_masked_iteration_vector_2d( &
      current, active_point_mask, current_vector)
    call pack_new_src_masked_iteration_vector_2d( &
      mapped, active_point_mask, mapped_vector)
    if (.not. allocated(basis%vector) .or. &
        size(basis%vector, 1) /= size(current_vector) .or. &
        size(basis%vector, 2) /= basis%retained_count .or. &
        .not. allocated(basis%mode_id) .or. &
        size(basis%mode_id) /= basis%retained_count) &
      error stop "zero-mode tangent basis has incompatible dimensions"

    residual = mapped_vector - current_vector
    original_squared = dot_product(residual, residual)
    report = zero_mode_projection_report_t()
    report%mode_count = basis%retained_count
    do mode = 1, basis%retained_count
      report%coefficient(basis%mode_id(mode)) = &
        dot_product(basis%vector(:, mode), residual)
      residual = residual - &
        report%coefficient(basis%mode_id(mode)) * basis%vector(:, mode)
    end do
    projected_squared = dot_product(residual, residual)
    report%original_residual_norm = sqrt(max(original_squared, 0.0_rk))
    report%projected_residual_norm = sqrt(max(projected_squared, 0.0_rk))
    report%removed_residual_norm = sqrt(max( &
      original_squared - projected_squared, 0.0_rk))
    if (report%original_residual_norm > tiny(1.0_rk)) &
      report%removed_fraction = report%removed_residual_norm / &
        report%original_residual_norm

    projected_mapped = current
    call unpack_new_src_masked_iteration_vector_2d( &
      current_vector + residual, active_point_mask, projected_mapped)
  end subroutine project_spinful_residual_away_from_tangents_2d


  subroutine build_mode_vector( &
      mesh, state, active_point_mask, mode, vector)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    logical, intent(in) :: active_point_mask(:)
    integer, intent(in) :: mode
    real(rk), intent(out) :: vector(:)

    complex(rk) :: tangent
    integer :: component, orbital, point, position, spin
    real(rk) :: coordinate(2), real_tangent

    if (size(vector) /= new_src_masked_iteration_vector_size_2d( &
        state, active_point_mask)) &
      error stop "zero-mode tangent vector has the wrong size"
    position = 0
    do spin = 1, 3
      do orbital = 1, 3
        do point = 1, state%point_count()
          if (.not. active_point_mask(point)) cycle
          tangent = order_parameter_tangent( &
            mesh, state, spin, orbital, point, mode)
          position = position + 1
          vector(position) = real(tangent, rk)
        end do
        do point = 1, state%point_count()
          if (.not. active_point_mask(point)) cycle
          tangent = order_parameter_tangent( &
            mesh, state, spin, orbital, point, mode)
          position = position + 1
          vector(position) = aimag(tangent)
        end do
      end do
    end do
    do component = 1, 3
      do point = 1, state%point_count()
        if (.not. active_point_mask(point)) cycle
        select case (mode)
        case (zero_mode_gauge)
          real_tangent = 0.0_rk
        case (zero_mode_translation_x)
          real_tangent = mean_field_derivative( &
            mesh, state, component, point, 1)
        case (zero_mode_translation_y)
          real_tangent = mean_field_derivative( &
            mesh, state, component, point, 2)
        case (zero_mode_orientation)
          coordinate = mesh%point_coordinate(point)
          real_tangent = coordinate(2) * mean_field_derivative( &
            mesh, state, component, point, 1) - &
            coordinate(1) * mean_field_derivative( &
            mesh, state, component, point, 2)
        case default
          error stop "unknown zero-mode tangent identifier"
        end select
        position = position + 1
        vector(position) = real_tangent
      end do
    end do
  end subroutine build_mode_vector


  function order_parameter_tangent( &
      mesh, state, spin, orbital, point, mode) result(tangent)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    integer, intent(in) :: spin, orbital, point, mode
    complex(rk) :: tangent

    complex(rk) :: value
    real(rk) :: coordinate(2)

    value = state%order_parameter(spin, orbital, point)
    select case (mode)
    case (zero_mode_gauge)
      tangent = cmplx(-aimag(value), real(value, rk), rk)
    case (zero_mode_translation_x)
      tangent = order_parameter_derivative( &
        mesh, state, spin, orbital, point, 1)
    case (zero_mode_translation_y)
      tangent = order_parameter_derivative( &
        mesh, state, spin, orbital, point, 2)
    case (zero_mode_orientation)
      coordinate = mesh%point_coordinate(point)
      ! This is the infinitesimal form of
      ! exp(i theta) A(R(-theta) r). The gauge factor keeps the winding-one
      ! bulk field exp(i phi) I fixed while the nonaxisymmetric core rotates.
      tangent = cmplx(-aimag(value), real(value, rk), rk) + &
        cmplx(coordinate(2), 0.0_rk, rk) * order_parameter_derivative( &
          mesh, state, spin, orbital, point, 1) - &
        cmplx(coordinate(1), 0.0_rk, rk) * order_parameter_derivative( &
          mesh, state, spin, orbital, point, 2)
    case default
      error stop "unknown zero-mode tangent identifier"
    end select
  end function order_parameter_tangent


  function order_parameter_derivative( &
      mesh, state, spin, orbital, point, axis) result(derivative)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    integer, intent(in) :: spin, orbital, point, axis
    complex(rk) :: derivative

    integer :: nodes(3), sample_point, sample_count, sample, x_node, y_node
    real(rk) :: weights(3)

    call point_nodes(mesh, point, x_node, y_node)
    call derivative_stencil(mesh, axis, merge(x_node, y_node, axis == 1), &
      nodes, weights, sample_count)
    derivative = cmplx(0.0_rk, 0.0_rk, rk)
    do sample = 1, sample_count
      if (axis == 1) then
        sample_point = mesh%point_index(nodes(sample), y_node)
      else
        sample_point = mesh%point_index(x_node, nodes(sample))
      end if
      derivative = derivative + cmplx(weights(sample), 0.0_rk, rk) * &
        state%order_parameter(spin, orbital, sample_point)
    end do
  end function order_parameter_derivative


  real(rk) function mean_field_derivative( &
      mesh, state, component, point, axis) result(derivative)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    integer, intent(in) :: component, point, axis

    integer :: nodes(3), sample_point, sample_count, sample, x_node, y_node
    real(rk) :: weights(3)

    call point_nodes(mesh, point, x_node, y_node)
    call derivative_stencil(mesh, axis, merge(x_node, y_node, axis == 1), &
      nodes, weights, sample_count)
    derivative = 0.0_rk
    do sample = 1, sample_count
      if (axis == 1) then
        sample_point = mesh%point_index(nodes(sample), y_node)
      else
        sample_point = mesh%point_index(x_node, nodes(sample))
      end if
      derivative = derivative + weights(sample) * &
        state%current_mean_field(component, sample_point)
    end do
  end function mean_field_derivative


  subroutine derivative_stencil( &
      mesh, axis, center_node, nodes, weights, sample_count)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    integer, intent(in) :: axis, center_node
    integer, intent(out) :: nodes(3), sample_count
    real(rk), intent(out) :: weights(3)

    integer :: last_node, sample
    real(rk) :: coordinate(3), evaluation

    if (axis == 1) then
      last_node = mesh%x_cell_count()
    else if (axis == 2) then
      last_node = mesh%y_cell_count()
    else
      error stop "Cartesian derivative axis must be one or two"
    end if
    nodes = 0
    weights = 0.0_rk
    if (last_node == 1) then
      sample_count = 2
      nodes(1:2) = [0, 1]
    else
      sample_count = 3
      if (center_node == 0) then
        nodes = [0, 1, 2]
      else if (center_node == last_node) then
        nodes = [last_node - 2, last_node - 1, last_node]
      else
        nodes = [center_node - 1, center_node, center_node + 1]
      end if
    end if
    do sample = 1, sample_count
      if (axis == 1) then
        coordinate(sample) = mesh%x_coordinate(nodes(sample))
      else
        coordinate(sample) = mesh%y_coordinate(nodes(sample))
      end if
    end do
    if (axis == 1) then
      evaluation = mesh%x_coordinate(center_node)
    else
      evaluation = mesh%y_coordinate(center_node)
    end if
    if (sample_count == 2) then
      weights(1) = -1.0_rk / (coordinate(2) - coordinate(1))
      weights(2) = -weights(1)
    else
      weights(1) = (2.0_rk * evaluation - coordinate(2) - coordinate(3)) / &
        ((coordinate(1) - coordinate(2)) * &
         (coordinate(1) - coordinate(3)))
      weights(2) = (2.0_rk * evaluation - coordinate(1) - coordinate(3)) / &
        ((coordinate(2) - coordinate(1)) * &
         (coordinate(2) - coordinate(3)))
      weights(3) = (2.0_rk * evaluation - coordinate(1) - coordinate(2)) / &
        ((coordinate(3) - coordinate(1)) * &
         (coordinate(3) - coordinate(2)))
    end if
  end subroutine derivative_stencil


  pure subroutine point_nodes(mesh, point, x_node, y_node)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    integer, intent(in) :: point
    integer, intent(out) :: x_node, y_node

    x_node = modulo(point - 1, mesh%x_point_count())
    y_node = (point - 1) / mesh%x_point_count()
  end subroutine point_nodes


  function zero_mode_name(mode) result(name)
    integer, intent(in) :: mode
    character(len=32) :: name

    select case (mode)
    case (zero_mode_gauge)
      name = "gauge"
    case (zero_mode_translation_x)
      name = "translation_x"
    case (zero_mode_translation_y)
      name = "translation_y"
    case (zero_mode_orientation)
      name = "gauge_compensated_orientation"
    case default
      name = "unknown"
    end select
  end function zero_mode_name


  subroutine require_valid_arguments(mesh, state, active_point_mask)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(in) :: state
    logical, intent(in) :: active_point_mask(:)

    if (.not. state%is_valid_for(mesh)) &
      error stop "zero-mode reference state does not match its mesh"
    if (size(active_point_mask) /= mesh%point_count() .or. &
        .not. any(active_point_mask)) &
      error stop "zero-mode active-point mask is invalid"
  end subroutine require_valid_arguments

end module spinful_tangent_projection_2d
