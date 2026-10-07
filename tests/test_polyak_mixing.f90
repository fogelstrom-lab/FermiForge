program test_polyak_mixing
  use he3_kinds, only: rk
  use polyak_mixing, only: polyak_mixing_t
  use legacy_anderson_mixing, only: anderson_report_t
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  type(polyak_mixing_t) :: p
  type(anderson_report_t) :: report
  real(rk) :: x(2), mapped(2), next(2), old(2), v(2), expected(2)
  integer :: k
  character(len=32) :: invalid
  call get_command_argument(1,invalid)
  if(len_trim(invalid)>0) then
    select case(trim(invalid))
    case('zero_drag')
      call p%initialize(2,2.0_rk,0.0_rk)
    case('negative_step')
      call p%initialize(2,-1.0_rk,0.5_rk)
    case('large_drag')
      call p%initialize(2,2.0_rk,1.1_rk)
    case('nonfinite_step')
      call p%initialize(2,ieee_value(0.0_rk,ieee_quiet_nan),0.5_rk)
    case('nonfinite_input')
      call p%initialize(2,2.0_rk,0.5_rk)
      x=ieee_value(0.0_rk,ieee_quiet_nan)
      call p%update(x,x,0.0_rk,next,report)
    end select
    stop 0 ! The WILL_FAIL regression must fail if invalid input was accepted.
  end if
  call p%initialize(2,2.0_rk,0.5_rk)
  x=[1.0_rk,2.0_rk]; mapped=0.75_rk*x
  call p%update(x,mapped,0.0_rk,next,report)
  call check(maxval(abs(next-0.5_rk*x))<1e-14_rk,'zero initial velocity')
  old=x; x=next
  call p%update(x,0.75_rk*x,0.0_rk,next,report)
  call check(maxval(abs(next))<1e-14_rk,'second step includes momentum')
  call check(report%history_size==0.and.report%iteration==2,'report has no Anderson history')
  call check(abs(report%mixing-2.0_rk)<1e-14_rk,'report alpha')
  ! Stale velocity must not move an already converged state.
  call p%update(x,x,1e-12_rk,next,report)
  call check(report%converged.and.all(next==x),'converged point does not drift')
  call p%update(x,0.75_rk*x,0.0_rk,next,report)
  call check(maxval(abs(next-0.5_rk*x))<1e-14_rk,'convergence clears velocity')
  call p%reset()
  call p%update(x,0.75_rk*x,0.0_rk,next,report)
  call check(report%iteration==1.and.maxval(abs(next-0.5_rk*x))<1e-14_rk,'reset')
  call p%initialize(2,0.4_rk,1.0_rk)
  do k=1,4
    mapped=0.75_rk*x
    call p%update(x,mapped,0.0_rk,next,report)
    call check(maxval(abs(next-(x+0.4_rk*(mapped-x))))<1e-14_rk,'Picard limit')
    x=next
  end do
  call p%initialize(2,2.0_rk,0.5_rk)
  x=old; v=0
  do k=1,100
    mapped=[0.7_rk,0.9_rk]*x
    v=0.5_rk*v+2.0_rk*(mapped-x)
    expected=x+v
    call p%update(x,mapped,0.0_rk,next,report)
    call check(maxval(abs(next-expected))<1e-14_rk,'exact E6 E7 recurrence')
    x=next
  end do
  call check(norm2(x)<1e-12_rk,'linear contraction converges')
  print *, 'Polyak tests passed'
contains
  subroutine check(ok,message)
    logical,intent(in)::ok
    character(*),intent(in)::message
    if(.not.ok) error stop message
  end subroutine
end program
