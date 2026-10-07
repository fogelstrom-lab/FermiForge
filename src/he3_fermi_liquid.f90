module he3_fermi_liquid
  use he3_kinds, only: rk
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: resolve_fs1
contains
  subroutine resolve_fs1(fs1,feedback_scale)
    real(rk), intent(in) :: fs1
    real(rk), intent(inout) :: feedback_scale
    ! Omitted Fs1 preserves historical feedback_scale inputs.
    if(fs1==huge(1.0_rk)) return
    if(.not.ieee_is_finite(fs1).or.fs1<0.0_rk) &
      error stop 'Fs1 must be finite and nonnegative for the current vortex drivers'
    feedback_scale=fs1/(1.0_rk+fs1/3.0_rk)
  end subroutine
end module
