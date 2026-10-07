module cylinder_bulk_seed
  use he3_kinds, only: rk
  implicit none
  private
  public :: make_cylinder_bulk_seed
contains
  subroutine make_cylinder_bulk_seed(initialization,physical_winding,gap,seed,symmetry_winding)
    integer, intent(in) :: initialization
    real(rk), intent(in) :: physical_winding,gap
    complex(rk), intent(out) :: seed(3,3)
    real(rk), intent(out) :: symmetry_winding
    integer :: j
    if(physical_winding/=0) error stop 'bulk cylinder seeds require zero physical phase winding'
    seed=cmplx(0._rk,0._rk,rk)
    select case(initialization)
    case(-1)
      do j=1,3
        seed(j,j)=gap
      end do
      symmetry_winding=0
    case(-2)
      ! new_src bulkA: d_z (p_x+i p_y), with its original seed amplitude.
      seed(3,1)=sqrt(2._rk)*gap
      seed(3,2)=cmplx(0._rk,sqrt(2._rk)*gap,rk)
      ! Only C_(0,+) is nonzero. In exp[i(m-s-k)phi], m=1
      ! makes its Cartesian seed constant: intrinsic chirality, NOT a vortex.
      symmetry_winding=1
    case default
      error stop 'cylinder initialization must be -1 (B) or -2 (A)'
    end select
  end subroutine
end module
