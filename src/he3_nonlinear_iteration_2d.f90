module he3_nonlinear_iteration_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t
  use spinful_state_2d, only : spinful_state_2d_t, allocate_spinful_state_2d
  use new_src_iteration_layout_2d, only : &
    new_src_iteration_vector_size_2d, pack_new_src_iteration_vector_2d, &
    unpack_new_src_iteration_vector_2d
  use legacy_anderson_mixing, only : legacy_anderson_t, anderson_report_t
  implicit none
  private

  type, public :: he3_field_residual_report_t
    integer :: value_count = 0
    integer :: maximum_point = 0
    real(rk) :: mean_absolute_residual = 0.0_rk
    real(rk) :: rms_residual = 0.0_rk
    real(rk) :: maximum_absolute_residual = 0.0_rk
    real(rk) :: relative_l2_residual = 0.0_rk
    real(rk) :: legacy_scaled_maximum_residual = 0.0_rk
    real(rk) :: maximum_point_rms_residual = 0.0_rk
  end type he3_field_residual_report_t

  public :: compute_he3_field_residual
  public :: update_he3_state_with_anderson

contains

  subroutine compute_he3_field_residual( &
      current, mapped, report, active_point_mask)
    type(spinful_state_2d_t), intent(in) :: current, mapped
    type(he3_field_residual_report_t), intent(out) :: report
    logical, intent(in), optional :: active_point_mask(:)

    real(rk), allocatable :: current_vector(:), mapped_vector(:), residual(:)
    real(rk) :: current_norm, point_sum_squared, value
    logical, allocatable :: active(:)
    integer :: component, orbital, point, position, spin

    if (current%point_count() /= mapped%point_count() .or. &
        current%point_count() < 1) &
      error stop "residual states have incompatible point counts"
    allocate(active(current%point_count()), source=.true.)
    if (present(active_point_mask)) then
      if (size(active_point_mask) /= current%point_count()) &
        error stop "residual active-point mask has the wrong size"
      active = active_point_mask
      if (.not. any(active)) &
        error stop "residual active-point mask contains no points"
      allocate(current_vector(21 * count(active)))
      allocate(mapped_vector(21 * count(active)))
      position = 0
      do point = 1, current%point_count()
        if (.not. active(point)) cycle
        do spin = 1, 3
          do orbital = 1, 3
            position = position + 1
            current_vector(position) = &
              real(current%order_parameter(spin, orbital, point), rk)
            mapped_vector(position) = &
              real(mapped%order_parameter(spin, orbital, point), rk)
            position = position + 1
            current_vector(position) = &
              aimag(current%order_parameter(spin, orbital, point))
            mapped_vector(position) = &
              aimag(mapped%order_parameter(spin, orbital, point))
          end do
        end do
        do component = 1, 3
          position = position + 1
          current_vector(position) = &
            current%current_mean_field(component, point)
          mapped_vector(position) = &
            mapped%current_mean_field(component, point)
        end do
      end do
    else
      call pack_new_src_iteration_vector_2d(current, current_vector)
      call pack_new_src_iteration_vector_2d(mapped, mapped_vector)
    end if
    residual = mapped_vector - current_vector

    report = he3_field_residual_report_t()
    report%value_count = size(residual)
    report%mean_absolute_residual = sum(abs(residual)) / &
      real(report%value_count, rk)
    report%rms_residual = sqrt(dot_product(residual, residual) / &
      real(report%value_count, rk))
    report%maximum_absolute_residual = maxval(abs(residual))
    current_norm = sqrt(max(dot_product(current_vector, current_vector), &
                            0.0_rk))
    if (current_norm <= tiny(1.0_rk)) current_norm = 1.0_rk
    report%relative_l2_residual = &
      sqrt(max(dot_product(residual, residual), 0.0_rk)) / current_norm
    report%legacy_scaled_maximum_residual = &
      real(report%value_count, rk) * report%maximum_absolute_residual / &
      current_norm

    do point = 1, current%point_count()
      if (.not. active(point)) cycle
      point_sum_squared = 0.0_rk
      do spin = 1, 3
        do orbital = 1, 3
          value = real(mapped%order_parameter(spin, orbital, point) - &
                       current%order_parameter(spin, orbital, point), kind=rk)
          point_sum_squared = point_sum_squared + value * value
          value = aimag(mapped%order_parameter(spin, orbital, point) - &
                        current%order_parameter(spin, orbital, point))
          point_sum_squared = point_sum_squared + value * value
        end do
      end do
      do component = 1, 3
        value = mapped%current_mean_field(component, point) - &
                current%current_mean_field(component, point)
        point_sum_squared = point_sum_squared + value * value
      end do
      value = sqrt(point_sum_squared / 21.0_rk)
      if (value > report%maximum_point_rms_residual) then
        report%maximum_point_rms_residual = value
        report%maximum_point = point
      end if
    end do
  end subroutine compute_he3_field_residual


  subroutine update_he3_state_with_anderson( &
      mesh, accelerator, current, mapped, tolerance, next, &
      anderson_report, residual_report)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(legacy_anderson_t), intent(inout) :: accelerator
    type(spinful_state_2d_t), intent(in) :: current, mapped
    real(rk), intent(in) :: tolerance
    type(spinful_state_2d_t), intent(out) :: next
    type(anderson_report_t), intent(out) :: anderson_report
    type(he3_field_residual_report_t), intent(out) :: residual_report

    real(rk), allocatable :: current_vector(:), mapped_vector(:), next_vector(:)

    if (.not. current%is_valid_for(mesh) .or. &
        .not. mapped%is_valid_for(mesh)) &
      error stop "Anderson state does not match the Cartesian mesh"
    if (accelerator%vector_size() /= &
        new_src_iteration_vector_size_2d(current)) &
      error stop "Anderson accelerator has the wrong 2D vector size"

    call compute_he3_field_residual(current, mapped, residual_report)
    call pack_new_src_iteration_vector_2d(current, current_vector)
    call pack_new_src_iteration_vector_2d(mapped, mapped_vector)
    allocate(next_vector(size(current_vector)))
    call accelerator%update( &
      current_vector, mapped_vector, tolerance, next_vector, anderson_report)
    call allocate_spinful_state_2d(mesh, next)
    call unpack_new_src_iteration_vector_2d(next_vector, next)
  end subroutine update_he3_state_with_anderson

end module he3_nonlinear_iteration_2d
