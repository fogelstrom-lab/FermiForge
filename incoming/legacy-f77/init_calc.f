c==========================================================================
c
      subroutine init_calc
c
c==========================================================================
c
      include 'qcv.dat' 
c     
c---  Give the input and rig the calculation
c     
      integer ii,i,j,ix,ien,inittyp
      real r,x,y,z,phi,dv,d,temp
      real a1x,b1x,a1y,b1y,a1z,b1z
      real a2x,b2x,a2y,b2y,a2z,b2z
      real a3x,b3x,a3y,b3y,a3z,b3z
      complex xx,xy,xz,yx,yy,yz,zx,zy,zz
c     
      read(*,*) t
      read(*,*) vort
      read(*,*) icyl
      read(*,*) aa0
      read(*,*) rx
      read(*,*) tmax
      read(*,*) eps 
      read(*,*) pmax(1),pmax(2),pmax(3),pmax(4),pmax(5)
      read(*,*) inmax(1),inmax(2),inmax(3),inmax(4),inmax(5)
      read(*,*) istart
      read(*,*) itmax
c      
c--- Some constants
c      
      pwave=.false.
      dwave=.false.
      cyl=.false.
      if(icyl.eq.1) cyl=.true.
      irep=1
      pwave=.true.
c        
c--- Get the Ozaki parameters
c           
      open(10,file='ozaki.dat',status='old')
      read(10,*) Ncmax,temp,x
      write(*,*) Ncmax,temp,x
      NcCutof=int(x)
      if(abs(temp-t).gt.1e-5) write(*,*) ' WRONG TEMP!!!!'
      do i=1,Ncmax
         read(10,*) zp(i),Rp(i)
      enddo       
      close(10)
c           
c--- get the bulk gap
c        
      maxfrec=int(NcCutof/t+.00001)
      aa0=aa0/(1.0+aa0/3.0)
      dx=rx/float(nx)
      call gapmag(t,deltat,irep)
c
c---  Set the coupling constant
c
      dv=0.0
      do ien=1,Ncmax,1
         dv=dv+t*Rp(ien)/zp(ien)   ! must be the same as in gap equation
      enddo
      esum=t/(log(t)+dv)
c        
c---  Initialize the order parameter
c
      xx=0.0
      xy=0.0
      xz=0.0
      yx=0.0
      yy=0.0
      yz=0.0
      zx=0.0
      zy=0.0
      zz=0.0

      if(istart.le.2) then
         do ix=0,nx,1
            x= float(ix)*dx

            if(istart.eq.-2) call bulkA(x,xx,xy,xz,yx,yy,yz,zx,zy,zz)
            if(istart.eq.-1) call bulk(x,xx,xy,xz,yx,yy,yz,zx,zy,zz)
            if(istart.eq.0) call nop(x,xx,xy,xz,yx,yy,yz,zx,zy,zz)
            if(istart.eq.1) call aop(x,xx,xy,xz,yx,yy,yz,zx,zy,zz)
            if(istart.eq.2) call dop(x,xx,xy,xz,yx,yy,yz,zx,zy,zz)

            dxx(ix)=xx
            dxy(ix)=xy
            dxz(ix)=xz

            dyx(ix)=yx
            dyy(ix)=yy
            dyz(ix)=yz

            dzx(ix)=zx
            dzy(ix)=zy
            dzz(ix)=zz

            vx(ix)=0.0
            vy(ix)=0.0
            vz(ix)=0.0
         enddo
      else
         open(1,file='op_xyz',status='unknown')
         open(3,file='curr',status='unknown')
         do ix=0,nx,1
            read(1,*) x,a1x,b1x,a1y,b1y,a1z,b1z
     +                 ,a2x,b2x,a2y,b2y,a2z,b2z
     +                 ,a3x,b3x,a3y,b3y,a3z,b3z
            dxx(ix)=cmplx(a1x,b1x)
            dxy(ix)=cmplx(a1y,b1y)
            dxz(ix)=cmplx(a1z,b1z)
            dyx(ix)=cmplx(a2x,b2x)
            dyy(ix)=cmplx(a2y,b2y)
            dyz(ix)=cmplx(a2z,b2z)
            dzx(ix)=cmplx(a3x,b3x)
            dzy(ix)=cmplx(a3y,b3y)
            dzz(ix)=cmplx(a3z,b3z)
            read(3,*) x,x,x,a1x,a1y,a1z
            vx(ix)=cmplx(a1x,0.0)
            vy(ix)=cmplx(a1y,0.0)
            vz(ix)=cmplx(a1z,0.0)
         enddo
         close(1)
         close(3)
      endif
c
c---  Pick the trajectories
c
      x=pi/tmax
      do i=1,tmax
         phi=float(i-1)*x+0.5*x
         kx(i)=cos(phi)
         ky(i)=sin(phi)
         awei(i)=1.0/tmax
      enddo

      open(10,file='gauss11.dat',status='old')
      x=0.0
      do i=1,11
         read(10,*) kz(i),pwei(i)
         pwei(i)=pwei(i)/2.0
         x=x+pwei(i)
         write(*,1300) i,kz(i),pwei(i),x
      enddo
      close(10)
      x=0.0
      y=0.0
      z=0.0
      do i=1,tmax
         do j=1,11
            x=x+pwei(j)*awei(i)
            y=y+pwei(j)*awei(i)*kx(i)*kx(i)*(1.0-kz(j)*kz(j))
            z=z+pwei(j)*awei(i)*kz(j)*kz(j)
         enddo
      enddo
      write(*,*) 'quads =',x,y,z
c
c---  Go with this setup
c
      eps=eps*deltat
      write(*,5200) ' n-c(0)/a-c(1)/d-c(2)  :',istart
      write(*,5000) ' Temperature           :',t
      write(*,5000) ' Vorticity             :',vort
      write(*,5050) '     Vcs0              :',esum
      write(*,5050) '     As1               :',aa0
      write(*,5050) ' Delta(T)              :',deltat
      write(*,5000) ' Grid radius           :',rx
      write(*,5100) ' Accuracy              :',eps
      write(*,5200) ' Nr. of trajs          :',tmax
      write(*,5200) ' Nr. of freqs          :',Ncmax
      write(*,5200) ' Max iterations        :',itmax
c
c--- Some format statements
c
 1200 format(20(1x,e14.6))
 1300 format(x,i3,20(1x,e14.6))
 2000 format(1x,f8.3,6(1x,e14.6))
 3000 format(2(1x,i4),6(1x,e14.6))
 5000 format(a,2x,f8.4)
 5050 format(a,2x,f8.4)
 5100 format(a,4x,e10.2)
 5200 format(a,2x,i5)
      return
      end
c
c---------------------------------------------------------------------
