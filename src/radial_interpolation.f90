module radial_interpolation
  use he3_kinds, only : rk
  use radial_mesh, only : radial_mesh_t, locate_radial_cell
  implicit none
  private

  public :: interpolate_radial_complex

contains

  function interpolate_radial_complex(mesh, values, radius) result(value)
    type(radial_mesh_t), intent(in) :: mesh
    complex(rk), intent(in) :: values(0:)
    real(rk), intent(in) :: radius
    complex(rk) :: value

    real(rk) :: x0, x1, x2, w0, w1, w2
    integer :: cell, i0, i1, i2, n

    if (.not. mesh%is_valid()) error stop "cannot interpolate an invalid mesh"
    if (size(values) /= mesh%point_count()) &
      error stop "radial field and mesh sizes differ"

    n = mesh%cell_count()
    cell = locate_radial_cell(mesh, radius)

    if (n == 1) then
      w1 = (radius - mesh%r(0)) / (mesh%r(1) - mesh%r(0))
      value = cmplx(1.0_rk - w1, 0.0_rk, kind=rk) * values(0) + &
              cmplx(w1, 0.0_rk, kind=rk) * values(1)
      return
    end if

    if (cell == 0) then
      i0 = 0
      i1 = 1
      i2 = 2
    else if (cell == n - 1) then
      i0 = n - 2
      i1 = n - 1
      i2 = n
    else
      i0 = cell - 1
      i1 = cell
      i2 = cell + 1
    end if

    x0 = mesh%r(i0)
    x1 = mesh%r(i1)
    x2 = mesh%r(i2)
    w0 = (radius - x1) * (radius - x2) / ((x0 - x1) * (x0 - x2))
    w1 = (radius - x0) * (radius - x2) / ((x1 - x0) * (x1 - x2))
    w2 = (radius - x0) * (radius - x1) / ((x2 - x0) * (x2 - x1))

    value = cmplx(w0, 0.0_rk, kind=rk) * values(i0) + &
            cmplx(w1, 0.0_rk, kind=rk) * values(i1) + &
            cmplx(w2, 0.0_rk, kind=rk) * values(i2)
  end function interpolate_radial_complex

end module radial_interpolation
