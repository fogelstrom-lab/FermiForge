c------------------------------------------------------------------
c
c     subroutine gap(tt,iirep,thetaO,deltat,eta0,psO,jsO,iemax)
      subroutine gapmag(tt,deltat,iirep)
c
c------------------------------------------------------------------
c
      include 'gap.dat'

      integer ien,iirep,iemax
      real deltat,tt,psO,jsO,thetaO,eta0

c-- Initialize the MPI world
      
      thetaO=0.0
      eta0=0.0
      psO=0.0
      jsO=0.0
      iemax=int(15.0/tt+0.00001)

c---------------------------
      
      irep=iirep
      t=tt
      rop1=0.20  
      theta=thetaO
      ps=psO
      ienmax=iemax

      pi=acos(-1.0)
      pi2=2.0*pi
      pih=0.5*pi
      s2=sqrt(2.0)
      eps=1.0e-8

      dv=0.0
      do ien=0,ienmax,1
         dv=dv+1.0/(ien+0.5)
      enddo
c     write(*,1000) t,psO,jsO,deltat,t/(log(t)+dv)

      call getaverages(10)
c     write(*,*) ' norm =',norm
      eta0=norm
      call getop
      call getcurr
      deltat=rop1
      jsO=js
c     write(*,1000) t,deltat,norm
 1000 format(4(x,f13.7))
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine getop
c
c---------------------------------------------------------------------
c
      include 'gap.dat'

      integer icnt
      real err,delta,gap1,F,dF
      parameter (delta=1.e-8)
c
c-- Initial guesses for the gaps
c
      gap1=rop1
      icnt=0

c-- Set up the matrices for the gap equation Newton step

 1    icnt=icnt+1

      if(icnt.gt.50) then
         write(*,*) ' The Sigmas do not converge'
         return
      endif
      rop1=gap1 
      call getF(F)
      rop1=gap1+delta
      call getF(dF)

c-- Numerical derivatives

      dF=(dF-F)/delta
      gap1=gap1-F/dF
      err=abs(F/dF)
c     write(*,1000) t,gap1,eps,F,dF

      if(err.gt.eps.and.abs(gap1).gt.eps) goto 1

 1000 format(9(1x,e12.5))
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine getF(F)
c
c------------------------------------------------------------------
c
      include 'gap.dat'

      integer ien,avtyp
      parameter (avtyp=0)
      real nrop1,vc1
      real F

      vc1=rop1*(log(t)+dv)
      nrop1=0.0
      do ien=0,ienmax,1
         en=t*(float(ien)+0.5)
         call getaverages(avtyp)
         nrop1=nrop1+aimag(ag2)
      enddo
      F=vc1-nrop1*t

 1000 format(9(1x,e12.5))
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine getcurr
c
c------------------------------------------------------------------
c
      include 'gap.dat'

      integer ien,avtyp
      parameter (avtyp=1)
      real njs

      njs=0.0
      do ien=0,ienmax,1
         en=t*(float(ien)+0.5)
         call getaverages(avtyp)
         njs=njs+real(ag3)
      enddo
      js=2.0*t*njs
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine getaverages(avtyp)
c                          
c---------------------------------------------------------------------
c
      include 'gap.dat'
      complex Atol,Rtol,gint
      real dlim,ulim,llim
      integer ierr,avtyp,i
      
      Atol=(1.e-6,1.e-6)
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
      return
      end
c
c-----------------------------------------------------------------------
c
      complex function gintegrand(ang)
c
c---------------------------------------------------------------------
c
      include 'gap.dat'
      real ang,u,x,v
      complex g2,g3,gap_cdiv,v2
      complex lam,hr,zen

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
        gintegrand=(u+w*v)*(u-w*v)
        goto 100
      endif

      u=u*norm 
      v=v*norm
      v2=(u+w*v)*(u-w*v)
      hr=rop1*rop1*v2    ! Set the Greens function
      zen=en-w*ps*sin(ang)
      lam=sqrt(zen*zen+hr)

      zen=-w*zen
      if(abs(lam).ne.0.0) then
         g3=gap_cdiv(zen,lam)
      else
	 g3=-w
      endif
c
c-- For the Order parameter
c
      if(inttyp.eq.11) then
         hr=w*rop1*(u-w*v)
         g2=gap_cdiv(hr,lam)
         gintegrand=(u+w*v)*g2
         goto 100
      endif

      if(inttyp.eq.17) then
         if(aimag(g3).gt.0.0)  g3= conjg(g3)
         gintegrand=g3*sin(ang)
         goto 100
      endif

 100  return
      end
c
c---------------------------------------------------------------------
c
      complex function gap_cdiv(d,n)
c
c------------------------------------------------------------------------
c
      implicit none
      complex  d,n,cn,cd
      real c
      cn=conjg(n)
      cd=d*cn
      c=real(n*cn)
      c=1.0/c
      gap_cdiv=cd*c
      return
      end
c
c---------------------------------------------------------------------
c
      complex function gint(llim,ulim,Atol,Rtol,ierr)
c
c------------------------------------------------------------------------
c
      implicit none
      integer ierr
      complex Atol, Rtol
      real llim, ulim
c
      integer stksz
      parameter (stksz = 512)
      dimension xstk(stksz), zstk(stksz)
      real xstk
      complex zstk
c
      integer itask, istk
      complex lz, mz, uz, zres, gintegrand
      real lx, mx, ux
      real  zero, two
      parameter (zero=0.0, two=2.0)
c
c------------------------------------------------------------------------
c
c   **Set some things to zero.**
      istk = 0
      ierr = 0
      gint = (0.0,0.0)
c
c   **Get the three points for the first Simpson integral.**
      lx = llim
      mx = (llim + ulim) / two
      ux = ulim
      lz = gintegrand(lx)
      mz = gintegrand(mx)
      uz = gintegrand(ux)
c
c   **Initialize the stacks.**
      xstk(1) = lx
      xstk(2) = mx
      xstk(3) = ux
      zstk(1) = lz
      zstk(2) = mz
      zstk(3) = uz
c
 10   continue
c
c   **Call the recursive part.**
      call gintrp(zres,lx,lz,mx,mz,ux,uz,istk,xstk,zstk,Atol,Rtol,itask)
c
c   **Check if we've exceeded the stack size.**
      if (3*(istk+1) .ge. stksz) then
        write (*,*) 'gint_: Stack size exceeded!'
        stop
      end if
c
c   **ierr keeps track of how big the stack gets.**
      if (istk .gt. ierr)  ierr=istk
c
c   **See if the current branch has converged.**
      if (itask .eq. 0) then
        gint = gint + zres
      else
        goto 10
      end if
c
c   **If the stack size is not zero, we have more integrating to do.**
c   **Pop the top three guys off the stack and integrate over them. **
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
c
      return
      end
c
c**********************************************************************
c
      subroutine gintrp(zres,lx,lz,mx,mz,ux,uz,istk,xstk,zstk,
     &                  Atol,Rtol,itask)
c
      implicit none
      integer stksz
      parameter (stksz = 512)

      integer istk, itask
      complex zres, lz, mz, uz, Atol, Rtol, zstk(stksz)
      real lx, mx, ux, xstk(stksz)
c
      complex ctest, lmz, umz, gintegrand
      real dx, lmx, umx, tmpr, tmpi
      real two, three, four
      parameter (two=2.0, three=3.0, four=4.0)
c
c
c   **Do the initial three-point Simpson's rule for reference.
      dx = (ux-lx) / two
      ctest = dx*(lz + four*mz + uz) / three
c
c   **Do the five-point Simpson's rule as a check.
      dx = dx / two
      lmx = (lx + mx) / two
      lmz = gintegrand(lmx)
      umx = (ux + mx) / two
      umz = gintegrand(umx)
      zres = dx*(lz + two*mz + four*(lmz+umz) + uz) / three
c
c   **Check the absolute tolerance and then the relative tolerance.**
c
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
c
c   **If we must go on, push the upper three guys onto the stack and **
c   **integrate over the lower three.                                **
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
c
      return
      end
c
c------------------------------------------------------------------------

