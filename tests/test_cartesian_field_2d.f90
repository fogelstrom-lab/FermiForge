program test_cartesian_field_2d
  use he3_kinds, only : rk
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, &
    make_uniform_cartesian_mesh, &
    make_symmetric_multiscale_cartesian_mesh, make_circular_active_mask
  use spinful_state_2d, only : spinful_state_2d_t, allocate_spinful_state_2d
  use bilinear_field_sampler_2d, only : bilinear_stencil_2d_t, &
                                        make_bilinear_stencil_2d, &
                                        sample_spinful_state_2d, &
                                        sample_spinful_stencil_2d, &
                                        sample_pair_potential_2d
  use new_src_iteration_layout_2d, only : new_src_values_per_2d_point, &
                                          new_src_iteration_vector_size_2d, &
                                          pack_new_src_iteration_vector_2d, &
                                          unpack_new_src_iteration_vector_2d
  implicit none

  call test_mesh_geometry_and_numbering()
  call test_multiscale_mesh_and_active_circle()
  call test_boundary_safe_stencils()
  call test_constant_and_affine_fields()
  call test_bilinear_field_exactness()
  call test_quadratic_grid_convergence()
  call test_axial_vortex_ray_convergence()
  call test_precomputed_stencil()
  call test_full_momentum_contraction()
  call test_new_src_iteration_layout()

  print '(a)', "Cartesian 2D field-sampling tests passed"

contains

  subroutine test_mesh_geometry_and_numbering()
    type(cartesian_mesh_2d_t) :: mesh
    real(rk), parameter :: tolerance = 64.0_rk * epsilon(1.0_rk)

    call make_uniform_cartesian_mesh(-2.0_rk, 3.0_rk, 5, &
                                     -4.0_rk, 2.0_rk, 3, mesh)

    call require(mesh%is_valid(), "uniform Cartesian mesh is invalid")
    call require(mesh%x_point_count() == 6, "wrong number of x points")
    call require(mesh%y_point_count() == 4, "wrong number of y points")
    call require(mesh%point_count() == 24, "wrong total Cartesian point count")
    call require(abs(mesh%x_spacing() - 1.0_rk) < tolerance, &
                 "wrong Cartesian x spacing")
    call require(abs(mesh%y_spacing() - 2.0_rk) < tolerance, &
                 "wrong Cartesian y spacing")
    call require(mesh%point_index(0, 0) == 1, &
                 "first Cartesian point index is wrong")
    call require(mesh%point_index(5, 0) == 6, &
                 "x must be the fastest Cartesian point index")
    call require(mesh%point_index(0, 1) == 7, &
                 "second Cartesian row starts at the wrong index")
    call require(mesh%point_index(5, 3) == 24, &
                 "last Cartesian point index is wrong")
  end subroutine test_mesh_geometry_and_numbering


  subroutine test_multiscale_mesh_and_active_circle()
    type(cartesian_mesh_2d_t) :: mesh
    type(spinful_state_2d_t) :: state
    complex(rk) :: expected_order(3, 3), sampled_order(3, 3)
    real(rk) :: expected_mean(3), sampled_mean(3)
    logical, allocatable :: active(:)
    logical :: inside
    integer :: center_point

    call make_symmetric_multiscale_cartesian_mesh( &
      8.0_rk, 2.0_rk, 5.0_rk, 0.5_rk, 1.0_rk, 2.0_rk, mesh)
    call require(mesh%is_valid(), "multiscale Cartesian mesh is invalid")
    call require(.not. mesh%is_uniform(), &
                 "multiscale Cartesian mesh was reported as uniform")
    call require(mesh%x_point_count() == 19 .and. &
                 mesh%y_point_count() == 19, &
                 "multiscale Cartesian point count is wrong")
    call require(abs(mesh%x_coordinate(9)) < epsilon(1.0_rk), &
                 "multiscale Cartesian origin is not an exact node")
    call require(abs(mesh%x_coordinate(13) - 2.0_rk) < epsilon(1.0_rk), &
                 "fine-region boundary is not an exact node")
    call require(abs(mesh%x_coordinate(16) - 5.0_rk) < epsilon(1.0_rk), &
                 "medium-region boundary is not an exact node")
    call require(abs(mesh%minimum_spacing() - 0.5_rk) < epsilon(1.0_rk), &
                 "multiscale minimum spacing is wrong")
    call require(abs(mesh%maximum_spacing() - 1.5_rk) < epsilon(1.0_rk), &
                 "multiscale maximum spacing is wrong")
    call require(mesh%x_cell_index(-8.0_rk) == 0 .and. &
                 mesh%x_cell_index(8.0_rk) == 17, &
                 "multiscale boundary cell lookup is wrong")

    call allocate_spinful_state_2d(mesh, state)
    call fill_state(mesh, state, bilinear_order, bilinear_mean)
    call sample_spinful_state_2d( &
      mesh, state, 3.2_rk, -6.1_rk, sampled_order, sampled_mean, inside)
    expected_order = bilinear_order(3.2_rk, -6.1_rk)
    expected_mean = bilinear_mean(3.2_rk, -6.1_rk)
    call require(inside, "multiscale in-domain sample was rejected")
    call require(maxval(abs(sampled_order - expected_order)) < 2.0e-12_rk, &
                 "multiscale Q1 sampler changed a bilinear complex field")
    call require(maxval(abs(sampled_mean - expected_mean)) < 2.0e-12_rk, &
                 "multiscale Q1 sampler changed a bilinear real field")

    call make_circular_active_mask(mesh, [0.0_rk, 0.0_rk], 5.0_rk, active)
    center_point = mesh%point_index(9, 9)
    call require(active(center_point), &
                 "circular active mask rejected the origin")
    call require(.not. active(mesh%point_index(0, 0)), &
                 "circular active mask accepted a corner halo point")
    call require(count(active) < mesh%point_count(), &
                 "circular active mask did not create a frozen halo")
  end subroutine test_multiscale_mesh_and_active_circle


  subroutine test_boundary_safe_stencils()
    type(cartesian_mesh_2d_t) :: mesh
    type(bilinear_stencil_2d_t) :: stencil
    real(rk), parameter :: tolerance = 128.0_rk * epsilon(1.0_rk)

    call make_uniform_cartesian_mesh(-1.0_rk, 1.0_rk, 4, &
                                     -2.0_rk, 2.0_rk, 8, mesh)

    call make_bilinear_stencil_2d(mesh, -1.0_rk, -2.0_rk, stencil)
    call require(stencil%inside, "lower-left boundary was rejected")
    call require(stencil%point(1) == mesh%point_index(0, 0), &
                 "lower-left boundary uses the wrong node")
    call require(abs(stencil%weight(1) - 1.0_rk) < tolerance, &
                 "lower-left boundary is not nodally exact")

    call make_bilinear_stencil_2d(mesh, 1.0_rk, 2.0_rk, stencil)
    call require(stencil%inside, "upper-right boundary was rejected")
    call require(stencil%point(4) == mesh%point_index(4, 8), &
                 "upper-right boundary uses the wrong node")
    call require(abs(stencil%weight(4) - 1.0_rk) < tolerance, &
                 "upper-right boundary is not nodally exact")
    call require(abs(sum(stencil%weight) - 1.0_rk) < tolerance, &
                 "bilinear weights do not form a partition of unity")

    call make_bilinear_stencil_2d(mesh, 1.01_rk, 0.0_rk, stencil)
    call require(.not. stencil%inside, "point beyond x boundary was accepted")
    call make_bilinear_stencil_2d(mesh, 0.0_rk, -2.01_rk, stencil)
    call require(.not. stencil%inside, "point beyond y boundary was accepted")
  end subroutine test_boundary_safe_stencils


  subroutine test_constant_and_affine_fields()
    type(cartesian_mesh_2d_t) :: mesh
    type(spinful_state_2d_t) :: state
    complex(rk) :: sampled_order(3, 3), expected_order(3, 3)
    real(rk) :: sampled_mean(3), expected_mean(3)
    real(rk), parameter :: points(2, 7) = reshape([ &
      -2.0_rk, -1.5_rk, &
       3.0_rk,  2.5_rk, &
      -2.0_rk,  2.5_rk, &
       3.0_rk, -1.5_rk, &
       0.0_rk,  0.0_rk, &
       1.125_rk, -0.625_rk, &
      -0.75_rk,  1.875_rk], [2, 7])
    real(rk) :: x, y
    logical :: inside
    integer :: sample

    call make_uniform_cartesian_mesh(-2.0_rk, 3.0_rk, 10, &
                                     -1.5_rk, 2.5_rk, 8, mesh)
    call allocate_spinful_state_2d(mesh, state)
    call fill_state(mesh, state, affine_order, affine_mean)

    call require(state%is_valid_for(mesh), "new 2D spinful state is invalid")
    do sample = 1, size(points, 2)
      x = points(1, sample)
      y = points(2, sample)
      call sample_spinful_state_2d(mesh, state, x, y, sampled_order, &
                                   sampled_mean, inside)
      expected_order = affine_order(x, y)
      expected_mean = affine_mean(x, y)
      call require(inside, "in-domain affine sample was rejected")
      call require(maxval(abs(sampled_order - expected_order)) < 3.0e-13_rk, &
                   "bilinear sampler does not reproduce an affine complex field")
      call require(maxval(abs(sampled_mean - expected_mean)) < 3.0e-13_rk, &
                   "bilinear sampler does not reproduce an affine real field")
    end do
  end subroutine test_constant_and_affine_fields


  subroutine test_bilinear_field_exactness()
    type(cartesian_mesh_2d_t) :: mesh
    type(spinful_state_2d_t) :: state
    complex(rk) :: sampled_order(3, 3), expected_order(3, 3)
    real(rk) :: sampled_mean(3), expected_mean(3)
    real(rk) :: x, y
    logical :: inside
    integer :: i

    call make_uniform_cartesian_mesh(-1.0_rk, 2.0_rk, 7, &
                                     -3.0_rk, 1.0_rk, 9, mesh)
    call allocate_spinful_state_2d(mesh, state)
    call fill_state(mesh, state, bilinear_order, bilinear_mean)

    do i = 0, 40
      x = -1.0_rk + 3.0_rk * real(mod(17 * i, 41), rk) / 40.0_rk
      y = -3.0_rk + 4.0_rk * real(mod(29 * i, 41), rk) / 40.0_rk
      call sample_spinful_state_2d(mesh, state, x, y, sampled_order, &
                                   sampled_mean, inside)
      expected_order = bilinear_order(x, y)
      expected_mean = bilinear_mean(x, y)
      call require(inside, "in-domain bilinear sample was rejected")
      call require(maxval(abs(sampled_order - expected_order)) < 2.0e-12_rk, &
                   "Q1 sampler is not exact for a bilinear complex field")
      call require(maxval(abs(sampled_mean - expected_mean)) < 2.0e-12_rk, &
                   "Q1 sampler is not exact for a bilinear real field")
    end do
  end subroutine test_bilinear_field_exactness


  subroutine test_quadratic_grid_convergence()
    real(rk) :: coarse_error, fine_error

    coarse_error = quadratic_error(8, 6)
    fine_error = quadratic_error(16, 12)

    call require(coarse_error > 0.0_rk, &
                 "quadratic interpolation test produced no coarse-grid error")
    call require(fine_error < coarse_error / 3.9_rk, &
                 "bilinear interpolation did not show second-order grid convergence")
  end subroutine test_quadratic_grid_convergence


  real(rk) function quadratic_error(number_of_x_cells, number_of_y_cells) &
      result(maximum_error)
    integer, intent(in) :: number_of_x_cells, number_of_y_cells

    type(cartesian_mesh_2d_t) :: mesh
    type(spinful_state_2d_t) :: state
    complex(rk) :: sampled_order(3, 3), expected_order(3, 3)
    real(rk) :: sampled_mean(3)
    real(rk) :: x, y
    logical :: inside
    integer :: cell_x, cell_y

    call make_uniform_cartesian_mesh(-1.0_rk, 1.0_rk, number_of_x_cells, &
                                     -2.0_rk, 2.0_rk, number_of_y_cells, mesh)
    call allocate_spinful_state_2d(mesh, state)
    call fill_state(mesh, state, quadratic_order, quadratic_mean)

    maximum_error = 0.0_rk
    do cell_y = 0, mesh%y_cell_count() - 1
      y = mesh%y_coordinate(cell_y) + 0.5_rk * mesh%y_spacing()
      do cell_x = 0, mesh%x_cell_count() - 1
        x = mesh%x_coordinate(cell_x) + 0.5_rk * mesh%x_spacing()
        call sample_spinful_state_2d(mesh, state, x, y, sampled_order, &
                                     sampled_mean, inside)
        expected_order = quadratic_order(x, y)
        call require(inside, "quadratic cell-center sample was rejected")
        maximum_error = max(maximum_error, &
                            maxval(abs(sampled_order - expected_order)))
      end do
    end do
  end function quadratic_error


  subroutine test_axial_vortex_ray_convergence()
    real(rk) :: coarse_error, fine_error

    coarse_error = axial_vortex_ray_error(24)
    fine_error = axial_vortex_ray_error(48)

    call require(coarse_error > 0.0_rk, &
                 "axial-vortex ray test produced no coarse-grid error")
    call require(fine_error < 0.45_rk * coarse_error, &
                 "axial-vortex ray sampling did not converge under grid refinement")
  end subroutine test_axial_vortex_ray_convergence


  real(rk) function axial_vortex_ray_error(number_of_cells) result(maximum_error)
    integer, intent(in) :: number_of_cells

    type(cartesian_mesh_2d_t) :: mesh
    type(spinful_state_2d_t) :: state
    complex(rk) :: sampled_order(3, 3), expected_order(3, 3)
    real(rk) :: sampled_mean(3)
    real(rk), parameter :: origin(2) = [0.137_rk, -0.219_rk]
    real(rk), parameter :: direction(2) = [0.8_rk, 0.6_rk]
    real(rk) :: path_coordinate, x, y
    logical :: inside
    integer :: sample

    call make_uniform_cartesian_mesh(-4.0_rk, 4.0_rk, number_of_cells, &
                                     -4.0_rk, 4.0_rk, number_of_cells, mesh)
    call allocate_spinful_state_2d(mesh, state)
    call fill_state(mesh, state, axial_vortex_order, axial_vortex_mean)

    maximum_error = 0.0_rk
    do sample = 0, 120
      path_coordinate = -3.0_rk + 6.0_rk * real(sample, rk) / 120.0_rk
      x = origin(1) + direction(1) * path_coordinate
      y = origin(2) + direction(2) * path_coordinate
      call sample_spinful_state_2d(mesh, state, x, y, sampled_order, &
                                   sampled_mean, inside)
      expected_order = axial_vortex_order(x, y)
      call require(inside, "axial-vortex ray left its declared test domain")
      maximum_error = max(maximum_error, &
                          maxval(abs(sampled_order - expected_order)))
    end do
  end function axial_vortex_ray_error


  subroutine test_precomputed_stencil()
    type(cartesian_mesh_2d_t) :: mesh
    type(spinful_state_2d_t) :: state
    type(bilinear_stencil_2d_t) :: stencil
    complex(rk) :: direct_order(3, 3), cached_order(3, 3)
    real(rk) :: direct_mean(3), cached_mean(3)
    logical :: direct_inside, cached_inside

    call make_uniform_cartesian_mesh(-2.0_rk, 2.0_rk, 8, &
                                     -2.0_rk, 2.0_rk, 8, mesh)
    call allocate_spinful_state_2d(mesh, state)
    call fill_state(mesh, state, bilinear_order, bilinear_mean)
    call make_bilinear_stencil_2d(mesh, 0.37_rk, -1.14_rk, stencil)

    call sample_spinful_state_2d(mesh, state, 0.37_rk, -1.14_rk, direct_order, &
                                 direct_mean, direct_inside)
    call sample_spinful_stencil_2d(state, stencil, cached_order, cached_mean, &
                                   cached_inside)

    call require(direct_inside .and. cached_inside, &
                 "precomputed in-domain stencil was rejected")
    call require(maxval(abs(cached_order - direct_order)) < epsilon(1.0_rk), &
                 "precomputed stencil changed the sampled order parameter")
    call require(maxval(abs(cached_mean - direct_mean)) < epsilon(1.0_rk), &
                 "precomputed stencil changed the sampled mean field")
  end subroutine test_precomputed_stencil


  subroutine test_full_momentum_contraction()
    type(cartesian_mesh_2d_t) :: mesh
    type(spinful_state_2d_t) :: state
    complex(rk) :: pair_potential(3), expected(3)
    real(rk) :: current_mean_field(3)
    real(rk), parameter :: momentum(3) = [0.0_rk, 0.0_rk, 0.75_rk]
    logical :: inside
    integer :: point, spin

    call make_uniform_cartesian_mesh(-1.0_rk, 1.0_rk, 2, &
                                     -1.0_rk, 1.0_rk, 2, mesh)
    call allocate_spinful_state_2d(mesh, state)
    do point = 1, mesh%point_count()
      do spin = 1, 3
        state%order_parameter(spin, 3, point) = &
          cmplx(real(spin, rk), -0.5_rk * real(spin, rk), kind=rk)
      end do
    end do

    call sample_pair_potential_2d(mesh, state, 0.2_rk, -0.3_rk, momentum, &
                                  pair_potential, current_mean_field, inside)
    do spin = 1, 3
      expected(spin) = cmplx(momentum(3), 0.0_rk, kind=rk) * &
        cmplx(real(spin, rk), -0.5_rk * real(spin, rk), kind=rk)
    end do

    call require(inside, "pair-potential sample was rejected")
    call require(maxval(abs(pair_potential - expected)) < 64.0_rk * epsilon(1.0_rk), &
                 "p_z was not retained in the triplet momentum contraction")
  end subroutine test_full_momentum_contraction


  subroutine test_new_src_iteration_layout()
    type(cartesian_mesh_2d_t) :: mesh
    type(spinful_state_2d_t) :: original, restored
    real(rk), allocatable :: vector(:)
    integer :: number_of_points, orbital, point, spin

    call make_uniform_cartesian_mesh(-1.0_rk, 1.0_rk, 3, &
                                     -2.0_rk, 2.0_rk, 2, mesh)
    call allocate_spinful_state_2d(mesh, original)
    call allocate_spinful_state_2d(mesh, restored)
    number_of_points = mesh%point_count()

    do point = 1, number_of_points
      do spin = 1, 3
        do orbital = 1, 3
          original%order_parameter(spin, orbital, point) = cmplx( &
            1000.0_rk * real(spin, rk) + 100.0_rk * real(orbital, rk) + &
            real(point, rk), &
            -1000.0_rk * real(spin, rk) - 100.0_rk * real(orbital, rk) - &
            real(point, rk), kind=rk)
        end do
      end do
      original%current_mean_field(:, point) = &
        [real(point, rk), -real(point, rk), 0.25_rk * real(point, rk)]
    end do

    call pack_new_src_iteration_vector_2d(original, vector)
    call require(new_src_values_per_2d_point == 21, &
                 "new_src layout must retain 21 real values per 2D point")
    call require(size(vector) == 21 * number_of_points, &
                 "new_src 2D iteration-vector length is wrong")
    call require(size(vector) == new_src_iteration_vector_size_2d(original), &
                 "reported new_src 2D iteration-vector size is inconsistent")
    call require(abs(vector(1) - real(original%order_parameter(1, 1, 1))) < &
                 epsilon(1.0_rk), "dxx real block starts at the wrong position")
    call require(abs(vector(number_of_points + 1) - &
                     aimag(original%order_parameter(1, 1, 1))) < &
                 epsilon(1.0_rk), "dxx imaginary block starts at the wrong position")
    call require(abs(vector(18 * number_of_points + 1) - &
                     original%current_mean_field(1, 1)) < epsilon(1.0_rk), &
                 "vx block starts at the wrong position")

    call unpack_new_src_iteration_vector_2d(vector, restored)
    call require(maxval(abs(restored%order_parameter - original%order_parameter)) < &
                 epsilon(1.0_rk), &
                 "new_src 2D order-parameter packing did not round-trip")
    call require(maxval(abs(restored%current_mean_field - &
                            original%current_mean_field)) < epsilon(1.0_rk), &
                 "new_src 2D mean-field packing did not round-trip")
  end subroutine test_new_src_iteration_layout


  subroutine fill_state(mesh, state, order_function, mean_function)
    type(cartesian_mesh_2d_t), intent(in) :: mesh
    type(spinful_state_2d_t), intent(inout) :: state
    interface
      pure function order_function(x, y) result(value)
        import rk
        real(rk), intent(in) :: x, y
        complex(rk) :: value(3, 3)
      end function order_function
      pure function mean_function(x, y) result(value)
        import rk
        real(rk), intent(in) :: x, y
        real(rk) :: value(3)
      end function mean_function
    end interface

    real(rk) :: x, y
    integer :: point, x_node, y_node

    do y_node = 0, mesh%y_cell_count()
      y = mesh%y_coordinate(y_node)
      do x_node = 0, mesh%x_cell_count()
        x = mesh%x_coordinate(x_node)
        point = mesh%point_index(x_node, y_node)
        state%order_parameter(:, :, point) = order_function(x, y)
        state%current_mean_field(:, point) = mean_function(x, y)
      end do
    end do
  end subroutine fill_state


  pure function affine_order(x, y) result(value)
    real(rk), intent(in) :: x, y
    complex(rk) :: value(3, 3)
    integer :: orbital, spin

    do orbital = 1, 3
      do spin = 1, 3
        value(spin, orbital) = cmplx( &
          real(10 * spin + orbital, rk) + 0.7_rk * x - 0.3_rk * y, &
          real(spin - 2 * orbital, rk) - 0.2_rk * x + 0.9_rk * y, kind=rk)
      end do
    end do
  end function affine_order


  pure function affine_mean(x, y) result(value)
    real(rk), intent(in) :: x, y
    real(rk) :: value(3)
    integer :: component

    do component = 1, 3
      value(component) = real(component, rk) + 0.4_rk * x - 0.8_rk * y
    end do
  end function affine_mean


  pure function bilinear_order(x, y) result(value)
    real(rk), intent(in) :: x, y
    complex(rk) :: value(3, 3)
    integer :: orbital, spin

    do orbital = 1, 3
      do spin = 1, 3
        value(spin, orbital) = cmplx( &
          real(spin + orbital, rk) + 0.3_rk * x - 0.6_rk * y + &
          0.11_rk * real(spin, rk) * x * y, &
          real(spin - orbital, rk) - 0.4_rk * x + 0.2_rk * y - &
          0.07_rk * real(orbital, rk) * x * y, kind=rk)
      end do
    end do
  end function bilinear_order


  pure function bilinear_mean(x, y) result(value)
    real(rk), intent(in) :: x, y
    real(rk) :: value(3)
    integer :: component

    do component = 1, 3
      value(component) = real(component, rk) - 0.2_rk * x + 0.5_rk * y + &
                         0.13_rk * real(component, rk) * x * y
    end do
  end function bilinear_mean


  pure function quadratic_order(x, y) result(value)
    real(rk), intent(in) :: x, y
    complex(rk) :: value(3, 3)
    integer :: orbital, spin

    do orbital = 1, 3
      do spin = 1, 3
        value(spin, orbital) = cmplx( &
          real(spin + orbital, rk) + 0.25_rk * x * x + 0.4_rk * y * y, &
          real(spin - orbital, rk) - 0.15_rk * x * x + 0.3_rk * y * y, kind=rk)
      end do
    end do
  end function quadratic_order


  pure function quadratic_mean(x, y) result(value)
    real(rk), intent(in) :: x, y
    real(rk) :: value(3)
    integer :: component

    do component = 1, 3
      value(component) = real(component, rk) + 0.2_rk * x * x - 0.3_rk * y * y
    end do
  end function quadratic_mean


  pure function axial_vortex_order(x, y) result(value)
    real(rk), intent(in) :: x, y
    complex(rk) :: value(3, 3)
    complex(rk) :: vortex_amplitude
    real(rk), parameter :: core_width = 0.7_rk
    integer :: component

    value = cmplx(0.0_rk, 0.0_rk, kind=rk)
    ! This smooth manufactured B-like field has unit phase winding and a core
    ! amplitude that vanishes without evaluating atan2 at the origin.
    vortex_amplitude = cmplx(x, y, kind=rk) / cmplx( &
      sqrt(x * x + y * y + core_width * core_width), 0.0_rk, kind=rk)
    do component = 1, 3
      value(component, component) = vortex_amplitude
    end do
  end function axial_vortex_order


  pure function axial_vortex_mean(x, y) result(value)
    real(rk), intent(in) :: x, y
    real(rk) :: value(3)
    real(rk), parameter :: core_width = 0.7_rk
    real(rk) :: denominator

    denominator = x * x + y * y + core_width * core_width
    value = [-y / denominator, x / denominator, 0.0_rk]
  end function axial_vortex_mean


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_cartesian_field_2d
