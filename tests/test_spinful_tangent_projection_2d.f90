program test_spinful_tangent_projection_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, &
    make_rectilinear_cartesian_mesh
  use spinful_state_2d, only : spinful_state_2d_t, allocate_spinful_state_2d
  use new_src_iteration_layout_2d, only : &
    pack_new_src_masked_iteration_vector_2d, &
    unpack_new_src_masked_iteration_vector_2d
  use spinful_tangent_projection_2d, only : &
    zero_mode_count, zero_mode_tangent_basis_2d_t, &
    zero_mode_projection_report_t, build_zero_mode_tangent_basis_2d, &
    project_spinful_residual_away_from_tangents_2d
  implicit none

  type(cartesian_mesh_2d_t) :: mesh
  type(spinful_state_2d_t) :: current, mapped, projected
  type(zero_mode_tangent_basis_2d_t) :: basis
  type(zero_mode_projection_report_t) :: report
  logical, allocatable :: active(:)
  real(rk), allocatable :: current_vector(:), mapped_vector(:), &
                           physical_vector(:), projected_vector(:)
  real(rk) :: gram(zero_mode_count, zero_mode_count)
  integer :: component, mode, orbital, point, spin, x_node, y_node
  real(rk) :: coefficient(4), component_value, coordinate(2)
  real(rk) :: orbital_value, spin_value, x, y
  complex(rk) :: derivative_x, derivative_y, tangent
  real(rk) :: derivative_mean_x, derivative_mean_y, mean_tangent

  call make_rectilinear_cartesian_mesh( &
    [-2.0_rk, -0.7_rk, 0.2_rk, 1.4_rk, 2.8_rk], &
    [-1.9_rk, -0.4_rk, 0.6_rk, 1.7_rk, 3.1_rk], mesh)
  call allocate_spinful_state_2d(mesh, current)
  do point = 1, mesh%point_count()
    coordinate = mesh%point_coordinate(point)
    x = coordinate(1)
    y = coordinate(2)
    do spin = 1, 3
      spin_value = real(spin, rk)
      do orbital = 1, 3
        orbital_value = real(orbital, rk)
        current%order_parameter(spin, orbital, point) = cmplx( &
          0.12_rk * spin_value - 0.07_rk * orbital_value + &
            0.03_rk * spin_value * x + 0.02_rk * orbital_value * y + &
            0.011_rk * x * x - 0.009_rk * y * y + &
            0.005_rk * x * y, &
          -0.04_rk * spin_value + 0.08_rk * orbital_value - &
            0.017_rk * orbital_value * x + 0.013_rk * spin_value * y - &
            0.006_rk * x * x + 0.007_rk * y * y - &
            0.004_rk * x * y, rk)
      end do
    end do
    do component = 1, 3
      component_value = real(component, rk)
      current%current_mean_field(component, point) = &
        0.02_rk * component_value + 0.014_rk * component_value * x - &
        0.009_rk * y + 0.003_rk * x * x + 0.002_rk * x * y
    end do
  end do
  allocate(active(mesh%point_count()), source=.false.)
  do y_node = 1, mesh%y_cell_count() - 1
    do x_node = 1, mesh%x_cell_count() - 1
      active(mesh%point_index(x_node, y_node)) = .true.
    end do
  end do

  call build_zero_mode_tangent_basis_2d( &
    mesh, current, active, .true., .true., .true., basis)
  call require(basis%requested_count == zero_mode_count .and. &
               basis%retained_count == zero_mode_count, &
    "analytic reference did not retain all four zero modes")
  gram = matmul(transpose(basis%vector), basis%vector)
  do mode = 1, zero_mode_count
    gram(mode, mode) = gram(mode, mode) - 1.0_rk
  end do
  call require(maxval(abs(gram)) < 2.0e-14_rk, &
    "zero-mode basis is not orthonormal")

  coefficient = [0.31_rk, -0.23_rk, 0.19_rk, 0.27_rk]
  mapped = current
  do point = 1, mesh%point_count()
    if (.not. active(point)) cycle
    coordinate = mesh%point_coordinate(point)
    x = coordinate(1)
    y = coordinate(2)
    do spin = 1, 3
      spin_value = real(spin, rk)
      do orbital = 1, 3
        orbital_value = real(orbital, rk)
        derivative_x = cmplx( &
          0.03_rk * spin_value + 0.022_rk * x + 0.005_rk * y, &
          -0.017_rk * orbital_value - 0.012_rk * x - 0.004_rk * y, rk)
        derivative_y = cmplx( &
          0.02_rk * orbital_value - 0.018_rk * y + 0.005_rk * x, &
          0.013_rk * spin_value + 0.014_rk * y - 0.004_rk * x, rk)
        tangent = cmplx(coefficient(1), 0.0_rk, rk) * cmplx( &
          -aimag(current%order_parameter(spin, orbital, point)), &
          real(current%order_parameter(spin, orbital, point), rk), rk) + &
          cmplx(coefficient(2), 0.0_rk, rk) * derivative_x + &
          cmplx(coefficient(3), 0.0_rk, rk) * derivative_y + &
          cmplx(coefficient(4), 0.0_rk, rk) * ( &
            cmplx( &
              -aimag(current%order_parameter(spin, orbital, point)), &
              real(current%order_parameter(spin, orbital, point), rk), rk) + &
            cmplx(y, 0.0_rk, rk) * derivative_x - &
            cmplx(x, 0.0_rk, rk) * derivative_y)
        mapped%order_parameter(spin, orbital, point) = &
          current%order_parameter(spin, orbital, point) + tangent
      end do
    end do
    do component = 1, 3
      component_value = real(component, rk)
      derivative_mean_x = 0.014_rk * component_value + &
        0.006_rk * x + 0.002_rk * y
      derivative_mean_y = -0.009_rk + 0.002_rk * x
      mean_tangent = coefficient(2) * derivative_mean_x + &
        coefficient(3) * derivative_mean_y + coefficient(4) * ( &
          y * derivative_mean_x - x * derivative_mean_y)
      mapped%current_mean_field(component, point) = &
        current%current_mean_field(component, point) + mean_tangent
    end do
  end do
  call project_spinful_residual_away_from_tangents_2d( &
    current, mapped, active, basis, projected, report)
  call require(report%removed_fraction > 1.0_rk - 2.0e-13_rk, &
    "collective-coordinate residual was not removed")
  call require(maxval(abs(projected%order_parameter - &
                          current%order_parameter)) < 2.0e-13_rk .and. &
               maxval(abs(projected%current_mean_field - &
                          current%current_mean_field)) < 2.0e-13_rk, &
    "projected collective-coordinate update changed the reference state")

  call pack_new_src_masked_iteration_vector_2d( &
    current, active, current_vector)
  allocate(physical_vector(size(current_vector)))
  do point = 1, size(physical_vector)
    physical_vector(point) = sin(0.173_rk * real(point, rk)) + &
      0.37_rk * cos(0.071_rk * real(point, rk))
  end do
  do mode = 1, basis%retained_count
    physical_vector = physical_vector - &
      dot_product(basis%vector(:, mode), physical_vector) * &
      basis%vector(:, mode)
  end do
  mapped = current
  call unpack_new_src_masked_iteration_vector_2d( &
    current_vector + physical_vector, active, mapped)
  call project_spinful_residual_away_from_tangents_2d( &
    current, mapped, active, basis, projected, report)
  call pack_new_src_masked_iteration_vector_2d( &
    projected, active, projected_vector)
  mapped_vector = current_vector + physical_vector
  call require(maxval(abs(projected_vector - mapped_vector)) < 3.0e-14_rk, &
    "projection changed a residual orthogonal to the zero modes")
  call require(report%removed_fraction < 2.0e-14_rk, &
    "orthogonal residual has a spurious zero-mode component")

  print '(a)', "spinful tangent-projection tests passed"

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_spinful_tangent_projection_2d
