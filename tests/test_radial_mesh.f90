program test_radial_mesh
  use he3_kinds, only : rk
  use radial_interpolation, only : interpolate_radial_complex
  use radial_mesh, only : radial_mesh_t, make_uniform_radial_mesh, &
                          make_boundary_refined_radial_mesh, &
                          locate_radial_cell
  implicit none

  call test_legacy_uniform_grid()
  call test_boundary_refinement()
  call test_nonuniform_quadratic_interpolation()

  print '(a)', "radial mesh tests passed"

contains

  subroutine test_legacy_uniform_grid()
    type(radial_mesh_t) :: mesh
    real(rk), parameter :: radius = 20.0_rk
    real(rk) :: expected
    integer :: i

    call make_uniform_radial_mesh(radius, 49, mesh)
    call require(mesh%is_valid(), "uniform mesh is invalid")
    call require(mesh%point_count() == 50, "legacy grid must have 50 points")
    call require(mesh%cell_count() == 49, "legacy grid must have 49 cells")

    do i = 0, 49
      expected = radius * real(i, rk) / 49.0_rk
      call require(abs(mesh%r(i) - expected) < 32.0_rk * epsilon(radius), &
                   "uniform node differs from the legacy coordinate")
    end do

    call require(locate_radial_cell(mesh, 0.0_rk) == 0, &
                 "origin cell lookup failed")
    call require(locate_radial_cell(mesh, radius) == 48, &
                 "surface cell lookup failed")
  end subroutine test_legacy_uniform_grid


  subroutine test_boundary_refinement()
    type(radial_mesh_t) :: mesh
    real(rk), allocatable :: spacing(:)
    real(rk), parameter :: radius = 80.0_rk
    integer :: i, n

    call make_boundary_refined_radial_mesh(radius=radius, &
                                            core_width=8.0_rk, &
                                            surface_width=8.0_rk, &
                                            fine_spacing=0.125_rk, &
                                            bulk_spacing=4.0_rk, &
                                            mesh=mesh)

    call require(mesh%is_valid(), "adaptive mesh is invalid")
    call require(abs(mesh%r(0)) < epsilon(radius), &
                 "adaptive mesh does not start at the origin")
    call require(abs(mesh%outer_radius() - radius) < epsilon(radius), &
                 "adaptive mesh does not end at the surface")
    call require(mesh%minimum_spacing() <= 0.125_rk, &
                 "core and surface refinement is too coarse")
    call require(mesh%maximum_spacing() <= 4.0_rk, &
                 "bulk spacing exceeds its requested maximum")
    call require(mesh%cell_count() < ceiling(radius / 0.125_rk), &
                 "adaptive mesh did not reduce the uniform fine-grid size")

    n = mesh%cell_count()
    allocate(spacing(0:n - 1))
    spacing = mesh%r(1:n) - mesh%r(0:n - 1)
    do i = 0, n - 2
      call require(max(spacing(i), spacing(i + 1)) <= &
                   2.0_rk * min(spacing(i), spacing(i + 1)) * &
                   (1.0_rk + 64.0_rk * epsilon(1.0_rk)), &
                   "adaptive mesh violates two-to-one balance")
    end do
  end subroutine test_boundary_refinement


  subroutine test_nonuniform_quadratic_interpolation()
    type(radial_mesh_t) :: mesh
    complex(rk), allocatable :: field(:)
    complex(rk) :: expected, interpolated
    real(rk) :: x
    integer :: i

    call make_boundary_refined_radial_mesh(radius=80.0_rk, &
                                            core_width=8.0_rk, &
                                            surface_width=8.0_rk, &
                                            fine_spacing=0.125_rk, &
                                            bulk_spacing=4.0_rk, &
                                            mesh=mesh)

    allocate(field(0:mesh%cell_count()))
    do i = 0, mesh%cell_count()
      field(i) = quadratic(mesh%r(i))
    end do

    do i = 0, 200
      x = 80.0_rk * real(i, rk) / 200.0_rk
      expected = quadratic(x)
      interpolated = interpolate_radial_complex(mesh, field, x)
      call require(abs(interpolated - expected) < 2.0e-11_rk, &
                   "quadratic interpolation is not exact on the adaptive mesh")
    end do
  end subroutine test_nonuniform_quadratic_interpolation


  pure complex(rk) function quadratic(x) result(value)
    real(rk), intent(in) :: x

    value = cmplx(x * x + 2.0_rk * x + 1.0_rk, &
                  -0.5_rk * x * x + 3.0_rk, kind=rk)
  end function quadratic


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_radial_mesh
