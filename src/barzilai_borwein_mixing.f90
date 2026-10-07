module barzilai_borwein_mixing
  use he3_kinds, only: rk
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use legacy_anderson_mixing, only: anderson_report_t
  implicit none
  private
  public :: bb_mixing_t

  type :: bb_mixing_t
    private
    real(rk), allocatable :: previous_x(:), previous_r(:)
    real(rk) :: initial=0.1_rk, lower=0.001_rk, upper=5.0_rk, growth=2.0_rk
    real(rk) :: previous_mixing=0.1_rk
    logical :: absolute_curvature=.false., have_previous=.false.
    integer :: iteration=0
    integer, public :: status=0, formula=0
    logical, public :: limited=.false.
    logical, public :: diagnostics_enabled=.false., secant_valid=.false., alignment_valid=.false.
    real(rk), public :: signed_cosine=0, residual_alignment=0, step_norm=0, secant_norm=0
    real(rk), public :: raw_bb1=0, raw_bb2=0
  contains
    procedure :: initialize
    procedure :: reset
    procedure :: update
  end type
contains
  subroutine initialize(self, n, initial, lower, upper, growth, absolute_curvature)
    class(bb_mixing_t), intent(inout) :: self
    integer, intent(in) :: n
    real(rk), intent(in) :: initial, lower, upper, growth
    logical, intent(in) :: absolute_curvature
    if (n < 1) error stop "BB vector must be nonempty"
    if (.not. all(ieee_is_finite([initial,lower,upper,growth]))) &
      error stop "BB parameters must be finite"
    if (lower <= 0 .or. initial < lower .or. upper < initial .or. &
        (growth /= 0 .and. growth <= 1)) &
      error stop "BB requires 0 < min <= initial <= max; growth=0 disables or must exceed 1"
    if (allocated(self%previous_x)) deallocate(self%previous_x,self%previous_r)
    allocate(self%previous_x(n),self%previous_r(n))
    self%initial=initial
    self%lower=lower
    self%upper=upper
    self%growth=growth
    self%absolute_curvature=absolute_curvature
    call self%reset()
  end subroutine

  subroutine reset(self)
    class(bb_mixing_t), intent(inout) :: self
    self%have_previous=.false.
    self%iteration=0
    self%status=0
    self%formula=0
    self%limited=.false.
    self%previous_mixing=self%initial
  end subroutine

  subroutine update(self, x, mapped, tolerance, next, report)
    class(bb_mixing_t), intent(inout) :: self
    real(rk), intent(in) :: x(:), mapped(:), tolerance
    real(rk), intent(out) :: next(:)
    type(anderson_report_t), intent(out) :: report
    real(rk) :: r(size(x)), s(size(x)), y(size(x))
    real(rk) :: nr, old_nr, ns, ny, cosine, alpha, proposal, denominator
    logical :: excessive_growth
    if (.not. allocated(self%previous_x)) error stop "BB not initialized"
    if (size(x) /= size(self%previous_x) .or. size(mapped) /= size(x) .or. &
        size(next) /= size(x)) error stop "BB vector size mismatch"
    if (.not. ieee_is_finite(tolerance) .or. tolerance < 0) error stop "invalid BB tolerance"
    if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(mapped))) &
      error stop "nonfinite BB input"
    r=mapped-x
    if (.not. all(ieee_is_finite(r))) error stop "nonfinite BB residual"
    nr=norm2(r)
    self%secant_valid=.false.
    self%alignment_valid=.false.
    self%signed_cosine=0; self%residual_alignment=0
    self%step_norm=0; self%secant_norm=0
    self%raw_bb1=0; self%raw_bb2=0
    if(self%diagnostics_enabled.and.self%have_previous) then
      s=x-self%previous_x
      y=self%previous_r-r
      ns=norm2(s); ny=norm2(y); old_nr=norm2(self%previous_r)
      self%step_norm=ns; self%secant_norm=ny
      if(nr>tiny(1.0_rk).and.old_nr>tiny(1.0_rk)) then
        self%residual_alignment=dot_product(r/nr,self%previous_r/old_nr)
        self%alignment_valid=ieee_is_finite(self%residual_alignment)
      end if
      if(ieee_is_finite(ns).and.ieee_is_finite(ny).and. &
         ns>epsilon(1.0_rk)*max(1.0_rk,norm2(x),norm2(self%previous_x)).and. &
         ny>epsilon(1.0_rk)*max(nr,old_nr,tiny(1.0_rk))) then
        self%signed_cosine=dot_product(s/ns,y/ny)
        if(abs(self%signed_cosine)>sqrt(epsilon(1.0_rk))) then
          self%raw_bb1=(ns/ny)/self%signed_cosine
          self%raw_bb2=(ns/ny)*self%signed_cosine
          self%secant_valid=all(ieee_is_finite([self%raw_bb1,self%raw_bb2]))
        end if
      end if
    end if
    if (.not. ieee_is_finite(nr)) error stop "nonfinite BB residual norm"
    self%iteration=self%iteration+1
    self%status=0
    self%formula=0
    self%limited=.false.
    report=anderson_report_t()
    report%iteration=self%iteration
    report%max_residual=maxval(abs(r))
    report%residual_norm=nr
    denominator=norm2(mapped)
    if (denominator <= sqrt(max(tolerance,tiny(1.0_rk)))) denominator=1
    report%legacy_relative_residual=nr/denominator
    report%converged=report%max_residual <= tolerance
    alpha=self%initial
    next=x
    if (report%converged) then
      report%mixing=0
      return
    end if
    if (self%have_previous) then
      s=x-self%previous_x
      ! Root residual is X-F(X), hence y = r_previous-r_current.
      y=self%previous_r-r
      ns=norm2(s)
      ny=norm2(y)
      old_nr=norm2(self%previous_r)
      excessive_growth=.false.
      if (self%growth > 0) excessive_growth=nr/self%growth > old_nr
      if (excessive_growth) then
        self%status=5 ! damp the NEXT move; this is not a rejection/line search
        alpha=max(self%lower,min(self%initial,0.5_rk*self%previous_mixing))
      else if (.not. ieee_is_finite(ns) .or. .not. ieee_is_finite(ny)) then
        self%status=6
      else if (ns <= epsilon(1.0_rk)*max(1.0_rk,norm2(x),norm2(self%previous_x)) .or. &
               ny <= epsilon(1.0_rk)*max(nr,old_nr,tiny(1.0_rk))) then
        self%status=3
      else
        cosine=dot_product(s/ns,y/ny)
        if (self%absolute_curvature) cosine=abs(cosine)
        if (cosine <= sqrt(epsilon(1.0_rk))) then
          self%status=4
        else
          self%formula=merge(1,2,mod(self%iteration,2)==1)
          if (self%formula==1) then
            proposal=(ns/ny)/cosine
          else
            proposal=(ns/ny)*cosine
          end if
          if (.not. ieee_is_finite(proposal)) then
            self%status=6
          else
            alpha=max(self%lower,min(self%upper,proposal))
            self%limited=alpha /= proposal
          end if
        end if
      end if
    end if
    next=x+alpha*r
    if (.not. all(ieee_is_finite(next))) then
      self%status=6
      alpha=self%lower
      next=x+alpha*r
      if (.not. all(ieee_is_finite(next))) error stop "nonfinite BB fallback update"
    end if
    report%mixing=alpha
    self%previous_x=x
    self%previous_r=r
    self%previous_mixing=alpha
    self%have_previous=.true.
  end subroutine
end module
