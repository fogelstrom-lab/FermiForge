program test_axial_radial_embedding
  use he3_kinds, only : rk
  use radial_mesh, only : radial_mesh_t, make_uniform_radial_mesh
  use cartesian_mesh_2d, only : cartesian_mesh_2d_t, make_uniform_cartesian_mesh
  use spinful_state_2d, only : spinful_state_2d_t, allocate_spinful_state_2d
  use bilinear_field_sampler_2d, only : sample_spinful_state_2d
  use order_parameter_basis, only : projection_plus, projection_zero, &
                                    projection_minus, &
                                    axial_harmonics_to_transport_cartesian
  use axial_radial_embedding, only : axial_radial_profile_t, &
                                     allocate_axial_radial_profile, &
                                     sample_axial_radial_profile, &
                                     embed_axial_radial_profile_2d
  implicit none

  call test_source_parity_interpolation()
  call test_embedding_ray_convergence()
  call test_outside_detection()

  print '(a)', "axial radial embedding tests passed"

contains

  subroutine test_source_parity_interpolation()
    type(radial_mesh_t) :: mesh
    type(axial_radial_profile_t) :: profile
    complex(rk) :: actual(3, 3), expected(3, 3), harmonic(3, 3)
    real(rk) :: actual_mean(3), angle, radius, x, y
    logical :: inside
    integer :: point

    call make_uniform_radial_mesh(4.0_rk, 16, mesh)
    call allocate_axial_radial_profile(mesh, profile)
    do point = 0, mesh%point_count() - 1
      radius = mesh%r(point)
      profile%harmonic(projection_plus, projection_zero, point) = &
        cmplx(1.0_rk + 0.2_rk * radius * radius, 0.0_rk, kind=rk)
      profile%harmonic(projection_zero, projection_zero, point) = &
        cmplx(radius, 0.0_rk, kind=rk)
      profile%harmonic(projection_plus, projection_plus, point) = &
        cmplx(radius, 0.0_rk, kind=rk)
      profile%harmonic(projection_minus, projection_plus, point) = &
        cmplx(radius, 0.0_rk, kind=rk)
      profile%azimuthal_mean_field(point) = radius
    end do

    x = -0.11_rk
    y = 0.07_rk
    radius = sqrt(x * x + y * y)
    angle = atan2(y, x)
    harmonic = cmplx(0.0_rk, 0.0_rk, kind=rk)
    harmonic(projection_plus, projection_zero) = &
      cmplx(1.0_rk + 0.2_rk * radius * radius, 0.0_rk, kind=rk)
    harmonic(projection_zero, projection_zero) = &
      cmplx(radius, 0.0_rk, kind=rk) * &
      exp(cmplx(0.0_rk, angle, kind=rk))
    harmonic(projection_plus, projection_plus) = &
      cmplx(radius, 0.0_rk, kind=rk) * &
      exp(cmplx(0.0_rk, -angle, kind=rk))
    harmonic(projection_minus, projection_plus) = &
      cmplx(radius, 0.0_rk, kind=rk) * &
      exp(cmplx(0.0_rk, angle, kind=rk))
    call axial_harmonics_to_transport_cartesian(harmonic, expected)

    call sample_axial_radial_profile(mesh, profile, x, y, 1.0_rk, actual, &
                                     actual_mean, inside)
    call require(inside, "near-origin radial sample was rejected")
    call require(maxval(abs(actual - expected)) < 2.0e-13_rk, &
                 "source parity extension changed an analytic axial profile")
    call require(maxval(abs(actual_mean - [-y, x, 0.0_rk])) < 2.0e-13_rk, &
                 "azimuthal radial mean field was not rotated into Cartesian form")
  end subroutine test_source_parity_interpolation


  subroutine test_embedding_ray_convergence()
    real(rk) :: coarse_error, fine_error

    coarse_error = embedding_ray_error(20)
    fine_error = embedding_ray_error(40)
    call require(coarse_error > 0.0_rk, &
                 "embedded axial ray test produced no coarse-grid error")
    call require(fine_error < 0.45_rk * coarse_error, &
                 "embedded axial ray sampling did not converge under refinement")
  end subroutine test_embedding_ray_convergence


  real(rk) function embedding_ray_error(number_of_cells) result(maximum_error)
    integer, intent(in) :: number_of_cells

    type(radial_mesh_t) :: radial_grid
    type(axial_radial_profile_t) :: profile
    type(cartesian_mesh_2d_t) :: cartesian_grid
    type(spinful_state_2d_t) :: state
    complex(rk) :: direct_order(3, 3), sampled_order(3, 3)
    real(rk) :: direct_mean(3), sampled_mean(3), s, x, y
    logical :: direct_inside, sampled_inside
    integer :: point, sample

    call make_uniform_radial_mesh(4.0_rk, 64, radial_grid)
    call allocate_axial_radial_profile(radial_grid, profile)
    do point = 0, radial_grid%point_count() - 1
      s = radial_grid%r(point)
      profile%harmonic(projection_plus, projection_zero, point) = &
        cmplx(1.0_rk + 0.2_rk * s * s, 0.1_rk * s * s, kind=rk)
      profile%harmonic(projection_zero, projection_zero, point) = &
        cmplx(s, 0.0_rk, kind=rk)
      profile%azimuthal_mean_field(point) = s
    end do

    call make_uniform_cartesian_mesh(-2.0_rk, 2.0_rk, number_of_cells, &
                                     -2.0_rk, 2.0_rk, number_of_cells, &
                                     cartesian_grid)
    call allocate_spinful_state_2d(cartesian_grid, state)
    call embed_axial_radial_profile_2d(radial_grid, profile, 1.0_rk, &
                                       cartesian_grid, state)

    maximum_error = 0.0_rk
    do sample = 0, 100
      s = -1.7_rk + 3.4_rk * real(sample, rk) / 100.0_rk
      x = 0.137_rk + 0.8_rk * s
      y = -0.219_rk + 0.6_rk * s
      call sample_axial_radial_profile(radial_grid, profile, x, y, 1.0_rk, &
                                       direct_order, direct_mean, direct_inside)
      call sample_spinful_state_2d(cartesian_grid, state, x, y, sampled_order, &
                                   sampled_mean, sampled_inside)
      call require(direct_inside .and. sampled_inside, &
                   "embedded ray sample left its declared domain")
      maximum_error = max(maximum_error, &
                          maxval(abs(sampled_order - direct_order)))
    end do
  end function embedding_ray_error


  subroutine test_outside_detection()
    type(radial_mesh_t) :: mesh
    type(axial_radial_profile_t) :: profile
    complex(rk) :: order_parameter(3, 3)
    real(rk) :: current_mean_field(3)
    logical :: inside

    call make_uniform_radial_mesh(2.0_rk, 8, mesh)
    call allocate_axial_radial_profile(mesh, profile)
    call sample_axial_radial_profile(mesh, profile, 2.1_rk, 0.0_rk, 1.0_rk, &
                                     order_parameter, current_mean_field, inside)
    call require(.not. inside, "point outside radial support was accepted")
    call require(maxval(abs(order_parameter)) < tiny(1.0_rk), &
                 "outside radial order-parameter sample is not zero")
    call require(maxval(abs(current_mean_field)) < tiny(1.0_rk), &
                 "outside radial mean-field sample is not zero")
  end subroutine test_outside_detection


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_axial_radial_embedding
