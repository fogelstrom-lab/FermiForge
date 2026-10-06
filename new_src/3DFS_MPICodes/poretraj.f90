!     Up to the subroutine trajectories, the beginning
!     is a test program for the subroutine.  I have used
!     it indebugging.

program poretraj
      use cylindrical_trajectories, only : trajectories
      implicit none
      real :: r, rsq, step, xo, yo, zo, px, py, pz, tp, tn, fp, fn
      integer :: nsn, ii
      real, allocatable :: rx(:), ry(:), rz(:), fx(:), fy(:), fz(:)
      real, allocatable :: rpx(:), rpy(:), rpz(:), fpx(:), fpy(:), fpz(:)

      write(*,*) ' give the radius '
      read(*,*) r
      write(*,*) ' give the step '
      read(*,*) step
      write(*,*) ' give the number of steps '
      read(*,*)  nsn
      if (r <= 0.0 .or. step <= 0.0 .or. nsn < 0) &
        error stop 'Require positive radius/step and nonnegative step count'
      allocate(rx(0:nsn), ry(0:nsn), rz(0:nsn), &
        fx(0:nsn), fy(0:nsn), fz(0:nsn), &
        rpx(0:nsn), rpy(0:nsn), rpz(0:nsn), &
        fpx(0:nsn), fpy(0:nsn), fpz(0:nsn))
      rsq=r*r

      write(6,*) 'initial position'
      read(5,*) xo,yo,zo
      write(6,*) 'direction cosines'
      read(5,*) px,py,pz

      if (xo*xo + yo*yo >= rsq) &
        error stop 'Initial point must lie strictly inside the cylinder'
      if (abs(px*px + py*py + pz*pz - 1.0) > 1.0e-5) &
        error stop 'Direction cosines must form a unit vector'
      call trajectories(xo,yo,zo,px,py,pz,r,rsq,step,nsn, &
        rx,ry,rz,fx,fy,fz,rpx,rpy,rpz,fpx,fpy,fpz)

      open(10,status='unknown',file='traj')
      write(10,200)xo,yo,zo
      write(10,200)px,py,pz
      write(10,*)
      do ii=1,nsn
       write(10,100) rx(ii),ry(ii),rz(ii),fx(ii),fy(ii),fz(ii)
      enddo
      write(10,*)
      do ii=1,nsn
       write(10,100) rpx(ii),rpy(ii),rpz(ii), &
      & -fpx(ii),-fpy(ii),-fpz(ii)
      enddo
      100 format(3(1x,f6.3),2x,3(1x,f6.3))
      200 format(3(1x,f6.3))
      tp=atan2(ry(nsn),rx(nsn))
      tn=atan2(fy(nsn),fx(nsn))
      fp=atan2((-rpx(nsn)*sin(tp)+rpy(nsn)*cos(tp)), &
      & (rpx(nsn)*cos(tp)+rpy(nsn)*sin(tp)))
      fn=atan2((fpx(nsn)*sin(tn)-fpy(nsn)*cos(tn)), &
      & (-fpx(nsn)*cos(tn)-fpy(nsn)*sin(tn)))
      write(*,300) fp,tp,fn,tn
      write(10,300) fp,tp,fn,tn
      300 format(2(1x,f6.3),2x,2(1x,f6.3))
      close(10)
end program poretraj
