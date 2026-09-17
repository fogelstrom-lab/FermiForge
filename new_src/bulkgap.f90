MODULE Bulkop
   implicit none

   integer, parameter ::  stksz = 512
   real, parameter :: tol = 1.0e-12
   complex, parameter :: czero = cmplx(0.0,0.0), cone = cmplx(1.0,0.0), ww = cmplx(0.0,1.0)

   integer :: ienmax, inttyp, irep
   real :: temp, pi2, s2, dv, norm
   real :: rop1, theta, ps, js
   complex :: ag2, ag3, en

contains
!
!------------------------------------------------------------------
!
   subroutine bulkgap(tt,deltat)
!
!------------------------------------------------------------------
!
      real, intent(in) :: tt
      real, intent(inout) :: deltat

      integer :: ien, iemax
      real :: psO, jsO, thetaO, eta0
      

!-- Initialize 

      temp = tt
      thetaO = 0.0
      eta0 = 0.0
      psO = 0.0
      jsO = 0.0
      iemax = int(15.0/temp+0.00001)

      rop1=0.20  
      theta=thetaO
      ps=psO
      ienmax=iemax

      pi2 = 2. * acos(-1.0)
      s2=sqrt(2.0)

      dv=0.0
      do ien=0,ienmax,1
         dv=dv+1.0/(ien+0.5)
      enddo

      call getaverages(10)

      eta0=norm
      call getop
      call getcurr
      deltat=rop1
      jsO=js
!
!---------------------------------------------------------------------
!
   end subroutine bulkgap
!
!---------------------------------------------------------------------
!
   subroutine getop
!
!---------------------------------------------------------------------
!
      integer :: icnt
      real, parameter :: delta=1.e-10
      real :: err, gap1, F, dF
!
!-- Initial guesses for the gaps
!
      gap1=rop1
      icnt=0

!-- Set up the matrices for the gap equation Newton step

 1    icnt=icnt+1

      if(icnt.gt.50) then
         print *, ' The Sigmas do not converge'
         return
      endif
      rop1=gap1 
      call getF(F)
      rop1=gap1+delta
      call getF(dF)

!-- Numerical derivatives

      dF=(dF-F)/delta
      gap1=gap1-F/dF
      err=abs(F/dF)

      if(err.gt.tol.and.abs(gap1).gt.tol) goto 1

!---------------------------------------------------------------------
!
   end subroutine getop
!
!---------------------------------------------------------------------
!
   subroutine getF(F)
!
!---------------------------------------------------------------------
!
      real, intent(inout) :: F
      integer, parameter :: avtyp = 0
      integer :: ien
      real :: nrop1, vc1

      vc1=rop1*(log(temp)+dv)
      nrop1=0.0
      do ien=0,ienmax,1
         en=temp*(float(ien)+0.5)
         call getaverages(avtyp)
         nrop1=nrop1+aimag(ag2)
      enddo
      F=vc1-temp*nrop1
!
!---------------------------------------------------------------------
!
   end subroutine getF
!
!---------------------------------------------------------------------
!
   subroutine getcurr
!
!---------------------------------------------------------------------
!
      integer, parameter :: avtyp = 1 
      integer :: ien
      real :: njs

      njs=0.0
      do ien=0,ienmax,1
         en=temp*(float(ien)+0.5)
         call getaverages(avtyp)
         njs=njs+real(ag3)
      enddo
      js=2.0*temp*njs
!
!---------------------------------------------------------------------
!
   end subroutine getcurr
!
!---------------------------------------------------------------------
!
   subroutine getaverages(avtyp)
!                          
!---------------------------------------------------------------------
!
      complex :: Atol, Rtol
      real :: dlim, ulim, llim
      integer :: ierr, avtyp, i
      
      Atol=(1.e-8,1.e-8)
      Rtol=(1.e-4,1.e-4)
      dlim=pi2/7.0
      ag2=0.0
      ag3=0.0
      if(avtyp.eq.10) then
         norm=0.0
         do i=1,7
            llim=dlim*(i-1)
            ulim=dlim*i
            inttyp=0
            norm=norm+real(gint(llim,ulim,Atol,Rtol,ierr)/pi2)
         enddo
         norm=1./sqrt(norm)
      elseif(avtyp.eq.0) then
         do i=1,7
            llim=dlim*(i-1)
            ulim=dlim*i
            inttyp=11
            ag2=ag2+gint(llim,ulim,Atol,Rtol,ierr)/pi2
         enddo
      else
         do i=1,7
            llim=dlim*(i-1)
            ulim=dlim*i
            inttyp=17
            ag3=ag3+gint(llim,ulim,Atol,Rtol,ierr)/pi2
         enddo
      endif
!
!-----------------------------------------------------------------------
!
   end subroutine getaverages
!
!-----------------------------------------------------------------------
!
   complex function integrand(ang)
!
!---------------------------------------------------------------------
!
      real, intent(in) :: ang
      real :: u, x, v
      complex :: g2, g3, v2, lam, hr, zen

      integrand = czero

      u = 0.0
      v = 0.0
      x=ang-theta
      if(irep.eq.0) then
         u=cos(2.*x)                     ! Dominant B1g
         v=0.0      
      endif
      if(irep.eq.1)then
         u=cos(x)       ! Dominant E2u: 
         v=sin(x)       ! Dominant  
      endif

      if(inttyp.eq.0) then
        integrand=(u+ww*v)*(u-ww*v)
        goto 100
      endif

      u=u*norm 
      v=v*norm
      v2=(u+ww*v)*(u-ww*v)
      hr=rop1*rop1*v2    ! Set the Greens function
      zen=en-ww*ps*sin(ang)
      lam=sqrt(zen*zen+hr)

      zen=-ww*zen
      if(abs(lam) > 1.0e-18) then
         g3=gap_cdiv(zen,lam)
      else
         g3=-ww
      endif
!
!-- For the Order parameter
!
      if(inttyp.eq.11) then
         hr=ww*rop1*(u-ww*v)
         g2=gap_cdiv(hr,lam)
         integrand=(u+ww*v)*g2
         goto 100
      endif

      if(inttyp.eq.17) then
         if(aimag(g3).gt.0.0)  g3= conjg(g3)
         integrand=g3*sin(ang)
         goto 100
      endif

 100  return
!
!---------------------------------------------------------------------
!
   end function integrand
!
!---------------------------------------------------------------------
!
   complex function gap_cdiv(d,n)
!
!---------------------------------------------------------------------
!
      complex, intent(in) :: d ,n
      complex ::  cn, cd
      real :: c

      cn=conjg(n)
      cd=d*cn
      c=real(n*cn)
      c=1.0/c
      gap_cdiv=cd*c
!
!---------------------------------------------------------------------
!
   end function gap_cdiv
!
!---------------------------------------------------------------------
!
   complex function gint(llim,ulim,Atol,Rtol,ierr)
!
!---------------------------------------------------------------------
!
      integer, intent(out) :: ierr
      complex, intent(in) :: Atol, Rtol
      real, intent(in) :: llim, ulim
 
      real, dimension(stksz) :: xstk 
      complex, dimension(stksz) :: zstk 
 
      integer :: itask, istk
      complex :: lz, mz, uz, zres
      real :: lx, mx, ux
      real, parameter :: two=2.0
!
!------------------------------------------------------------------------
!
!   **Set some things to zero.**
      istk = 0
      ierr = 0
      gint = (0.0,0.0)
!
!   **Get the three points for the first Simpson integral.**
      lx = llim
      mx = (llim + ulim) / two
      ux = ulim
      lz = integrand(lx)
      mz = integrand(mx)
      uz = integrand(ux)
!
!   **Initialize the stacks.**
      xstk(1) = lx
      xstk(2) = mx
      xstk(3) = ux
      zstk(1) = lz
      zstk(2) = mz
      zstk(3) = uz
!
 10   continue
!
!   **Call the recursive part.**
      call gintrp(zres,lx,lz,mx,mz,ux,uz,istk,xstk,zstk,Atol,Rtol,itask)
!
!   **Check if we've exceeded the stack size.**
      if (3*(istk+1) .ge. stksz) then
        write (*,*) 'gint_: Stack size exceeded!'
        stop
      end if
!
!   **ierr keeps track of how big the stack gets.**
      if (istk .gt. ierr)  ierr=istk
!
!   **See if the current branch has converged.**
      if (itask .eq. 0) then
        gint = gint + zres
      else
        goto 10
      end if
!
!   **If the stack size is not zero, we have more integrating to do.**
!   **Pop the top three guys off the stack and integrate over them. **
      if (istk .ne. 0) then
        istk = istk - 1
        lx = xstk(3*istk+1)
        mx = xstk(3*istk+2)
        ux = xstk(3*istk+3)
        lz = zstk(3*istk+1)
        mz = zstk(3*istk+2)
        uz = zstk(3*istk+3)
        goto 10
      end if
!
!---------------------------------------------------------------------
!
   end function gint
!
!---------------------------------------------------------------------
!
   subroutine gintrp(zres,lx,lz,mx,mz,ux,uz,istk,xstk,zstk,Atol,Rtol,itask)
!
!---------------------------------------------------------------------
!
      real, dimension(stksz), intent(inout) :: xstk 
      complex, dimension(stksz), intent(inout) :: zstk 
 
      integer, intent(inout) :: istk, itask
      real, intent(inout) :: lx, mx, ux 
      complex, intent(in) :: Atol, Rtol
      complex, intent(out) :: zres
      complex, intent(inout) :: lz, mz, uz
 
      complex :: ctest, lmz, umz
      real :: dx, lmx, umx, tmpr, tmpi
      real, parameter :: two=2.0, three=3.0, four=4.0
 
 
!   **Do the initial three-point Simpson's rule for reference.
      dx = (ux-lx) / two
      ctest = dx*(lz + four*mz + uz) / three
 
!   **Do the five-point Simpson's rule as a check.
      dx = dx / two
      lmx = (lx + mx) / two
      lmz = integrand(lmx)
      umx = (ux + mx) / two
      umz = integrand(umx)
      zres = dx*(lz + two*mz + four*(lmz+umz) + uz) / three
 
!   **Check the absolute tolerance and then the relative tolerance.**
 
      if ( abs(real(zres-ctest)) .lt. real(Atol) ) then
        if ( abs(aimag(zres-ctest)) .lt. aimag(Atol) ) then
          itask = 0
        else
          tmpi = two*abs(aimag(zres-ctest)) / abs(aimag(zres+ctest))
          if ( tmpi .lt. aimag(Rtol) ) then
            itask = 0
          else
            itask = 1
          end if
        end if
      else
        tmpr = two*abs(real(zres-ctest))  / abs(real(zres+ctest))
        if ( tmpr .lt. real(Rtol) ) then
          if ( abs(aimag(zres-ctest)) .lt. aimag(Atol) ) then
            itask = 0
          else
            tmpi = two*abs(aimag(zres-ctest)) / abs(aimag(zres+ctest))
            if ( tmpi .lt. aimag(Rtol) ) then
              itask = 0
            else
              itask = 1
            end if
          end if
        else
          itask = 1
        end if
      end if
 
!   **If we must go on, push the upper three guys onto the stack and **
!   **integrate over the lower three.                                **
      if (itask .eq. 1) then
        xstk(3*istk+1) = mx
        xstk(3*istk+2) = umx
        xstk(3*istk+3) = ux
        zstk(3*istk+1) = mz
        zstk(3*istk+2) = umz
        zstk(3*istk+3) = uz
        istk = istk + 1
        ux = mx
        mx = lmx
        uz = mz
        mz = lmz
      end if
!
!---------------------------------------------------------------------
!
   end subroutine gintrp
!
!---------------------------------------------------------------------
!
end MODULE Bulkop
!
!---------------------------------------------------------------------

