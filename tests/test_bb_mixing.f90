program test_bb_mixing
  use he3_kinds, only: rk
  use barzilai_borwein_mixing, only: bb_mixing_t
  use legacy_anderson_mixing, only: anderson_report_t
  implicit none
  type(bb_mixing_t) :: bb
  type(bb_mixing_t) :: quiet, recorded
  type(anderson_report_t) :: quiet_report, recorded_report
  real(rk) :: quiet_next(2),recorded_next(2)
  integer :: k
  type(anderson_report_t) :: report
  real(rk) :: x(2), mapped(2), next(2)
  call bb%initialize(2,0.1_rk,0.001_rk,5.0_rk,2.0_rk,.false.)
  x=[1.0_rk,2.0_rk]
  mapped=0.5_rk*x
  call bb%update(x,mapped,0.0_rk,next,report)
  call check(maxval(abs(next-0.95_rk*x))<1e-14_rk,"initial damping")
  x=next
  call bb%update(x,0.5_rk*x,0.0_rk,next,report)
  call check(abs(report%mixing-2.0_rk)<1e-13_rk .and. bb%formula==2,"BB2 secant")
  call bb%update(next,next,1e-12_rk,x,report)
  call check(report%converged,"convergence")
  bb%diagnostics_enabled=.true.
  call bb%reset()
  x=1
  call bb%update(x,x+1,0.0_rk,next,report)
  x=next
  call bb%update(x,x+1,0.0_rk,next,report)
  call check(bb%status==3 .and. abs(report%mixing-0.1_rk)<1e-14_rk,"constant residual fallback")
  call bb%reset()
  x=1
  call bb%update(x,2*x,0.0_rk,next,report)
  x=next
  call bb%update(x,2*x,0.0_rk,next,report)
  call check(bb%status==4,"negative curvature fallback")
  call bb%initialize(2,0.1_rk,0.001_rk,5.0_rk,2.0_rk,.true.)
  x=1
  call bb%update(x,2*x,0.0_rk,next,report)
  x=next
  call bb%update(x,2*x,0.0_rk,next,report)
  call check(abs(report%mixing-1.0_rk)<1e-13_rk,"absolute curvature mode")
  call check(bb%secant_valid.and.bb%raw_bb1<0.and.bb%raw_bb2<0, "preserve signed curvature")
  call check(bb%signed_cosine<0.and.bb%residual_alignment>0, "signed alignment")
  call bb%reset()
  x=1
  call bb%update(x,x+1,0.0_rk,next,report)
  x=next
  call bb%update(x,x+3,0.0_rk,next,report)
  call check(bb%status==5,"growth damping")
  call bb%initialize(2,0.1_rk,0.001_rk,1.0_rk,2.0_rk,.false.)
  x=1
  call bb%update(x,0.5_rk*x,0.0_rk,next,report)
  x=next
  call bb%update(x,0.5_rk*x,0.0_rk,next,report)
  call check(bb%limited .and. abs(report%mixing-1.0_rk)<1e-14_rk,"upper bound")
  x=next
  call bb%update(x,0.5_rk*x,0.0_rk,next,report)
  call check(bb%formula==1 .and. bb%limited,"BB1 alternation")
  call bb%initialize(2,0.5_rk,0.5_rk,5.0_rk,100.0_rk,.false.)
  x=1
  call bb%update(x,-9*x,0.0_rk,next,report)
  x=next
  call bb%update(x,-9*x,0.0_rk,next,report)
  call check(bb%limited .and. abs(report%mixing-0.5_rk)<1e-14_rk,"lower bound")
  call bb%reset()
  x=0
  call bb%update(x,x,0.0_rk,next,report)
  call check(report%converged .and. report%legacy_relative_residual==0,"zero-vector convergence")
  call bb%initialize(2,0.1_rk,0.001_rk,5.0_rk,2.0_rk,.false.)
  x=[1.0_rk,2.0_rk]
  mapped=x-[0.5_rk,1.0_rk]*x
  call bb%update(x,mapped,0.0_rk,next,report)
  ! Deliberately alter the displacement: use actual s, not alpha_previous*r_previous.
  x=[0.9_rk,1.7_rk]
  mapped=x-[0.5_rk,1.0_rk]*x
  call bb%update(x,mapped,0.0_rk,next,report)
  call check(abs(report%mixing-0.095_rk/0.0925_rk)<1e-13_rk,"actual displacement secant")
  call bb%initialize(2,1.0_rk,0.001_rk,100.0_rk,0.0_rk,.true.)
  x=1
  call bb%update(x,x+1,0.0_rk,next,report)
  call check(abs(report%mixing-1.0_rk)<1e-14_rk,"radial-style startup")
  x=next
  call bb%update(x,x+3,0.0_rk,next,report)
  call check(bb%status==0 .and. abs(report%mixing-0.5_rk)<1e-14_rk, &
    "disabled growth guard must allow BB on increased residual")
  call bb%reset()
  x=1
  call bb%update(x,0.999_rk*x,0.0_rk,next,report)
  x=next
  call bb%update(x,0.999_rk*x,0.0_rk,next,report)
  call check(bb%limited .and. abs(report%mixing-100.0_rk)<1e-12_rk,"radial-style cap")
  call quiet%initialize(2,1.0_rk,0.001_rk,100.0_rk,0.0_rk,.true.)
  call recorded%initialize(2,1.0_rk,0.001_rk,100.0_rk,0.0_rk,.true.)
  recorded%diagnostics_enabled=.true.
  x=[1.0_rk,2.0_rk]
  do k=1,20
    mapped=[0.7_rk,0.95_rk]*x+0.03_rk*sin(x)
    call quiet%update(x,mapped,0.0_rk,quiet_next,quiet_report)
    call recorded%update(x,mapped,0.0_rk,recorded_next,recorded_report)
    call check(all(quiet_next==recorded_next),"diagnostics must leave updates bitwise unchanged")
    call check(quiet_report%mixing==recorded_report%mixing,"diagnostics changed alpha")
    x=quiet_next
  end do
  print *, "BB mixing tests passed"
contains
  subroutine check(ok,message)
    logical, intent(in) :: ok
    character(*), intent(in) :: message
    if (.not. ok) error stop message
  end subroutine
end program
