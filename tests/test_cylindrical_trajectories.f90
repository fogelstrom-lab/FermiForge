! Standalone regression against the frozen fixed-form cylinder routine.
! Link with trajectories.f after renaming trajectories/zigzag to legacy_*.
program test_cylindrical_trajectories
  use cylindrical_trajectories, only : trajectories
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  integer, parameter :: nsn = 120, legacy_size = 20000
  real :: modern(0:nsn,12), legacy(0:legacy_size,12)
  real :: xo, yo, zo, px, py, pz, old_pz, r, step, tolerance
  integer :: sample, case_index
  interface
    subroutine legacy_trajectories(xo,yo,zo,px,py,pz,r,rsq,step,nsn, &
        rx,ry,rz,fx,fy,fz,rpx,rpy,rpz,fpx,fpy,fpz)
      integer :: nsn
      real :: xo,yo,zo,px,py,pz,r,rsq,step
      real :: rx(0:20000),ry(0:20000),rz(0:20000)
      real :: fx(0:20000),fy(0:20000),fz(0:20000)
      real :: rpx(0:20000),rpy(0:20000),rpz(0:20000)
      real :: fpx(0:20000),fpy(0:20000),fpz(0:20000)
    end subroutine legacy_trajectories
  end interface

  r = 2.0
  tolerance = 5000.0 * epsilon(r)
  do case_index = 1, 9
    xo = 0.3; yo = -0.2; zo = 0.7
    px = 0.48; py = 0.64; pz = 0.6
    step = 0.13
    select case (case_index)
    case (2)
      pz = -pz
    case (3)
      px = 1.0; py = 0.0; pz = 0.0
    case (4)
      px = 0.0; py = -1.0; pz = 0.0
    case (5)
      px = 0.0; py = 0.0; pz = 1.0
    case (6)
      px = 0.0; py = 0.0; pz = -1.0
    case (7)
      ! Exact wall hits on a diameter; equality branch matters.
      xo = 0.0; yo = 0.0; zo = 0.0
      px = 1.0; py = 0.0; pz = 0.0; step = 0.5
    case (8)
      ! Multiple reflections within one integration step.
      step = 11.0
    case (9)
      xo = 1.999; yo = 0.0
      px = 0.6; py = 0.0; pz = 0.8
    end select
    old_pz = pz
    call legacy_trajectories(xo,yo,zo,px,py,old_pz,r,r*r,step,nsn, &
      legacy(:,1),legacy(:,2),legacy(:,3),legacy(:,4),legacy(:,5),legacy(:,6), &
      legacy(:,7),legacy(:,8),legacy(:,9),legacy(:,10),legacy(:,11),legacy(:,12))
    call trajectories(xo,yo,zo,px,py,pz,r,r*r,step,nsn, &
      modern(:,1),modern(:,2),modern(:,3),modern(:,4),modern(:,5),modern(:,6), &
      modern(:,7),modern(:,8),modern(:,9),modern(:,10),modern(:,11),modern(:,12))
    if (.not. all(ieee_is_finite(modern))) error stop 'Nonfinite trajectory'
    if (any(modern /= legacy(0:nsn,:))) error stop 'Legacy trajectory mismatch'
    if (pz /= old_pz) error stop 'Legacy pz side effect mismatch'
    if (maxval(modern(:,1)**2 + modern(:,2)**2) > r*r + tolerance) &
      error stop 'Reverse ray escaped cylinder'
    if (maxval(modern(:,4)**2 + modern(:,5)**2) > r*r + tolerance) &
      error stop 'Forward ray escaped cylinder'
    do sample = 0, nsn
      if (abs(sum(modern(sample,7:9)**2) - 1.0) > tolerance) &
        error stop 'Reverse momentum norm'
      if (abs(sum(modern(sample,10:12)**2) - 1.0) > tolerance) &
        error stop 'Forward momentum norm'
    end do
  end do
  print *, 'PASS: nine cylinder cases match the legacy routine exactly'
end program test_cylindrical_trajectories
