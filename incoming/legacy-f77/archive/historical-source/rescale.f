c-------------------------------------------------------------------
c
      program rescale
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'
      integer ix,iy
      real x,y,ndx,Rxnew,ax,bx,ay,by,az,bz,d,c
      real avx(-nx:nx,-nx:nx),avy(-nx:nx,-nx:nx),avz(-nx:nx,-nx:nx)
      complex axx(-nx:nx,-nx:nx),axy(-nx:nx,-nx:nx),axz(-nx:nx,-nx:nx)
      complex ayx(-nx:nx,-nx:nx),ayy(-nx:nx,-nx:nx),ayz(-nx:nx,-nx:nx)
      complex azx(-nx:nx,-nx:nx),azy(-nx:nx,-nx:nx),azz(-nx:nx,-nx:nx)
      complex intpol

      pi=acos(-1.)
      sqrt2=sqrt(2.0)
      myid=0
      nproc=1
c
c-- get input ops
c
      open(1,file='op_x',status='unknown')
      open(2,file='op_y',status='unknown')
      open(3,file='op_z',status='unknown')
      open(4,file='curr',status='unknown')

      do iy=-nx,nx,1
         do ix=-nx,nx,1
            d=0.0
            read(1,*) x,y,ax,bx,ay,by,az,bz
            dxx(ix,iy)=cmplx(ax,bx)
            dxy(ix,iy)=cmplx(ay,by)
            dxz(ix,iy)=cmplx(az,bz)

            read(2,*) x,y,ax,bx,ay,by,az,bz
            dyx(ix,iy)=cmplx(ax,bx)
            dyy(ix,iy)=cmplx(ay,by)
            dyz(ix,iy)=cmplx(az,bz)

            read(3,*) x,y,ax,bx,ay,by,az,bz
            dzx(ix,iy)=cmplx(ax,bx)
            dzy(ix,iy)=cmplx(ay,by)
            dzz(ix,iy)=cmplx(az,bz)

            read(4,*) x,y,bx,bx,ax,ay,az
            vx(ix,iy)=ax
            vy(ix,iy)=ay
            vz(ix,iy)=az
         enddo
         read(1,*)
         read(2,*)
         read(3,*)
         read(4,*)
      enddo
      close(1)
      close(2)
      close(3)
      close(4)

      write(*,*) ' Give the old Rx'
      read(*,*) rx
      dx=Rx/float(nx) 
      write(*,*) ' Give the new Rx'
      read(*,*) Rxnew
      if(Rxnew.gt.Rx) stop ' can only shrink the lattice '
      ndx=Rxnew/float(nx) 
      do iy=-nx,nx,1
         y=ndx*float(iy)
         do ix=-nx,nx,1
            x=ndx*float(ix)
            axx(ix,iy)=intpol(x,y,1)
            axy(ix,iy)=intpol(x,y,2)
            axz(ix,iy)=intpol(x,y,3)

            ayx(ix,iy)=intpol(x,y,4)
            ayy(ix,iy)=intpol(x,y,5)
            ayz(ix,iy)=intpol(x,y,6)

            azx(ix,iy)=intpol(x,y,7)
            azy(ix,iy)=intpol(x,y,8)
            azz(ix,iy)=intpol(x,y,9)

            avx(ix,iy)=real(intpol(x,y,10))
            avy(ix,iy)=real(intpol(x,y,11))
            avz(ix,iy)=real(intpol(x,y,12))
         enddo
      enddo

      open(1,file='op_x',status='unknown')
      open(2,file='op_y',status='unknown')
      open(3,file='op_z',status='unknown')
      open(4,file='curr',status='unknown')

      do iy=-nx,nx,1
         y=ndx*float(iy)
         do ix=-nx,nx,1
            x=ndx*float(ix)

            write(1,1000) x,y,axx(ix,iy),axy(ix,iy),axz(ix,iy)
            write(2,1000) x,y,ayx(ix,iy),ayy(ix,iy),ayz(ix,iy)
            write(3,1000) x,y,azx(ix,iy),azy(ix,iy),azz(ix,iy)

            d=axx(ix,iy)*conjg(axx(ix,iy))
     +       +axy(ix,iy)*conjg(axy(ix,iy))
     +       +axz(ix,iy)*conjg(axz(ix,iy))
     +       +ayx(ix,iy)*conjg(ayx(ix,iy))
     +       +ayy(ix,iy)*conjg(ayy(ix,iy))
     +       +ayz(ix,iy)*conjg(ayz(ix,iy))
     +       +azx(ix,iy)*conjg(azx(ix,iy))
     +       +azy(ix,iy)*conjg(azy(ix,iy))
     +       +azz(ix,iy)*conjg(azz(ix,iy))
            d=sqrt(d/3.0)
            c=sqrt(avx(ix,iy)**2+avy(ix,iy)**2)

            write(4,1000) x,y,d,c,avx(ix,iy)
     +                           ,avy(ix,iy),avz(ix,iy)
         enddo
         write(1,1100)
         write(2,1100)
         write(3,1100)
         write(4,1100)
      enddo

 1000 format(2(1x,f8.3),30(1x,e14.6))
 2000 format(1x,f8.3,30(1x,e14.6))
 1100 format()

      end        
c
c---------------------------------------------------------------------
c
      complex function intpol(x,y,icomp)
c
c--------------------------------------------------------------------
c
      include 'qcv.dat'

      integer ix,iy,idx,idy,icomp
      real x,y,p,q
      complex axx(-nx:nx,-nx:nx),ayy(-nx:nx,-nx:nx),azz(-nx:nx,-nx:nx)
      common /as/axx,ayy,azz

      intpol=0.0
      if(x.ge.0) then
         ix=int(x/dx)
         idx=1
         p=(x-dx*ix)/dx
      else
         ix=-int(abs(x)/dx)
         idx=-1
         p=(dx*ix-x)/dx
      endif
      if(y.ge.0) then
         iy=int(y/dx)
         idy=1
         q=(y-dx*iy)/dx
      else
         iy=-int(abs(y)/dx)
         idy=-1
         q=(dx*iy-y)/dx
      endif
      if(abs(p).gt.1) write(*,*) ' intpol: p ',p,ix,x
      if(abs(q).gt.1) write(*,*) ' intpol: q ',q,iy,y

c     write(*,1000) ' intpol ',icomp,ix,iy,x,y
      if(icomp.eq. 1) intpol=(1.-p)*(1-q)*dxx(ix,iy)+
     +                       p*(1-q)*dxx(ix+idx,iy)+
     +                       q*(1-p)*dxx(ix,iy+idy)+
     +                       p*q*dxx(ix+idx,iy+idy)
      if(icomp.eq. 2) intpol=(1.-p)*(1-q)*dxy(ix,iy)+
     +                       p*(1-q)*dxy(ix+idx,iy)+
     +                       q*(1-p)*dxy(ix,iy+idy)+
     +                       p*q*dxy(ix+idx,iy+idy)
      if(icomp.eq. 3) intpol=(1.-p)*(1-q)*dxz(ix,iy)+
     +                       p*(1-q)*dxz(ix+idx,iy)+
     +                       q*(1-p)*dxz(ix,iy+idy)+
     +                       p*q*dxz(ix+idx,iy+idy)
      if(icomp.eq. 4) intpol=(1.-p)*(1-q)*dyx(ix,iy)+
     +                       p*(1-q)*dyx(ix+idx,iy)+
     +                       q*(1-p)*dyx(ix,iy+idy)+
     +                       p*q*dyx(ix+idx,iy+idy)
      if(icomp.eq. 5) intpol=(1.-p)*(1-q)*dyy(ix,iy)+
     +                       p*(1-q)*dyy(ix+idx,iy)+
     +                       q*(1-p)*dyy(ix,iy+idy)+
     +                       p*q*dyy(ix+idx,iy+idy)
      if(icomp.eq. 6) intpol=(1.-p)*(1-q)*dyz(ix,iy)+
     +                       p*(1-q)*dyz(ix+idx,iy)+
     +                       q*(1-p)*dyz(ix,iy+idy)+
     +                       p*q*dyz(ix+idx,iy+idy)
      if(icomp.eq. 7) intpol=(1.-p)*(1-q)*dzx(ix,iy)+
     +                       p*(1-q)*dzx(ix+idx,iy)+
     +                       q*(1-p)*dzx(ix,iy+idy)+
     +                       p*q*dzx(ix+idx,iy+idy)
      if(icomp.eq. 8) intpol=(1.-p)*(1-q)*dzy(ix,iy)+
     +                       p*(1-q)*dzy(ix+idx,iy)+
     +                       q*(1-p)*dzy(ix,iy+idy)+
     +                       p*q*dzy(ix+idx,iy+idy)
      if(icomp.eq. 9) intpol=(1.-p)*(1-q)*dzz(ix,iy)+
     +                       p*(1-q)*dzz(ix+idx,iy)+
     +                       q*(1-p)*dzz(ix,iy+idy)+
     +                       p*q*dzz(ix+idx,iy+idy)
      if(icomp.eq.10) intpol=(1.-p)*(1-q)*vx(ix,iy)+
     +                       p*(1-q)*vx(ix+idx,iy)+
     +                       q*(1-p)*vx(ix,iy+idy)+
     +                       p*q*vx(ix+idx,iy+idy)
      if(icomp.eq.11) intpol=(1.-p)*(1-q)*vy(ix,iy)+
     +                       p*(1-q)*vy(ix+idx,iy)+
     +                       q*(1-p)*vy(ix,iy+idy)+
     +                       p*q*vy(ix+idx,iy+idy)
      if(icomp.eq.12) intpol=(1.-p)*(1-q)*vz(ix,iy)+
     +                       p*(1-q)*vz(ix+idx,iy)+
     +                       q*(1-p)*vz(ix,iy+idy)+
     +                       p*q*vz(ix+idx,iy+idy)

      if(icomp.eq. -1) intpol=(1.-p)*(1-q)*axx(ix,iy)+
     +                        p*(1-q)*axx(ix+idx,iy)+
     +                        q*(1-p)*axx(ix,iy+idy)+
     +                        p*q*axx(ix+idx,iy+idy)
      if(icomp.eq. -5) intpol=(1.-p)*(1-q)*ayy(ix,iy)+
     +                        p*(1-q)*ayy(ix+idx,iy)+
     +                        q*(1-p)*ayy(ix,iy+idy)+
     +                        p*q*ayy(ix+idx,iy+idy)
      if(icomp.eq. -9) intpol=(1.-p)*(1-q)*azz(ix,iy)+
     +                        p*(1-q)*azz(ix+idx,iy)+
     +                        q*(1-p)*azz(ix,iy+idy)+
     +                        p*q*azz(ix+idx,iy+idy)

 1000 format(a,3(x,i5),2(x,f7.3))
      return
      end
c
c---------------------------------------------------------------------
