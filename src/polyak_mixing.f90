! Polyak heavy-ball fixed-point iteration, SuperConga Appendix E, E4/E6/E7.
! Residual r=F(x)-x; v_next=(1-drag)*v+step_size*r; x_next=x+v_next.
! This is NOT the objective-value-based Polyak step-size method.
module polyak_mixing
  use he3_kinds, only: rk
  use legacy_anderson_mixing, only: anderson_report_t
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: polyak_mixing_t
  type :: polyak_mixing_t
    private
    real(rk), allocatable :: velocity(:)
    real(rk) :: step_size=2.0_rk, drag=0.5_rk
    integer :: iteration=0
  contains
    procedure :: initialize
    procedure :: reset
    procedure :: update
  end type
contains
  subroutine initialize(self,n,step_size,drag)
    class(polyak_mixing_t), intent(inout) :: self
    integer, intent(in) :: n
    real(rk), intent(in) :: step_size,drag
    if(n<1) error stop 'Polyak vector must be nonempty'
    if(.not.all(ieee_is_finite([step_size,drag]))) error stop 'nonfinite Polyak parameters'
    ! drag=1 is also permitted as the exactly memoryless Picard limit.
    if(step_size<=0 .or. drag<=0 .or. drag>1) &
      error stop 'Polyak requires step_size>0 and 0<drag<=1'
    if(allocated(self%velocity)) deallocate(self%velocity)
    allocate(self%velocity(n))
    self%step_size=step_size
    self%drag=drag
    call self%reset()
  end subroutine

  subroutine reset(self)
    class(polyak_mixing_t), intent(inout) :: self
    if(allocated(self%velocity)) self%velocity=0
    self%iteration=0
  end subroutine

  subroutine update(self,x,mapped,tolerance,next,report)
    class(polyak_mixing_t), intent(inout) :: self
    real(rk), intent(in) :: x(:),mapped(:),tolerance
    real(rk), intent(out) :: next(:)
    type(anderson_report_t), intent(out) :: report
    real(rk), allocatable :: residual(:),move(:)
    real(rk) :: denominator
    if(.not.allocated(self%velocity)) error stop 'Polyak not initialized'
    if(size(x)/=size(self%velocity).or.size(mapped)/=size(x).or.size(next)/=size(x)) &
      error stop 'Polyak vector size mismatch'
    if(.not.ieee_is_finite(tolerance).or.tolerance<0) error stop 'invalid Polyak tolerance'
    if(.not.all(ieee_is_finite(x)).or..not.all(ieee_is_finite(mapped))) &
      error stop 'nonfinite Polyak input'
    residual=mapped-x
    if(.not.all(ieee_is_finite(residual))) error stop 'nonfinite Polyak residual'
    self%iteration=self%iteration+1
    report=anderson_report_t()
    report%iteration=self%iteration
    report%max_residual=maxval(abs(residual))
    report%residual_norm=norm2(residual)
    if(.not.ieee_is_finite(report%residual_norm)) error stop 'nonfinite Polyak norm'
    denominator=norm2(mapped)
    if(.not.ieee_is_finite(denominator)) error stop 'nonfinite Polyak mapped norm'
    if(denominator<=sqrt(max(tolerance,tiny(1.0_rk)))) denominator=1
    report%legacy_relative_residual=report%residual_norm/denominator
    report%converged=report%max_residual<=tolerance
    next=x
    if(report%converged) then
      ! A converged fixed point must not be displaced by stale momentum.
      self%velocity=0
      return
    end if
    move=(1-self%drag)*self%velocity+self%step_size*residual
    next=x+move
    if(.not.all(ieee_is_finite(move)).or..not.all(ieee_is_finite(next))) &
      error stop 'nonfinite Polyak update; reduce step size or increase drag'
    self%velocity=move
    report%mixing=self%step_size
  end subroutine
end module
