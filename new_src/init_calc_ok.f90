Module Initialisation
   use Global_variables
   use Bulkop
contains
!--------------------------------------------------------------------------
!
   subroutine init_calc
!
!--------------------------------------------------------------------------
!
! --  Get the input and rig the calculation
!     
      integer :: i,j,ix,ien
      real ::  x,y,z,phi,dv,temp, tol
      real ::  a1x,b1x,a1y,b1y,a1z,b1z
      real ::  a2x,b2x,a2y,b2y,a2z,b2z
      real ::  a3x,b3x,a3y,b3y,a3z,b3z
      complex :: xx,xy,xz,yx,yy,yz,zx,zy,zz
      complex :: cpp,cpo,cpm,cop,coo,com,cmp,cmo,cmm
!
! -- input data given in qcv.inp
!
      read *, temp
      read *, vort
      read *, icyl
      read *, aa0
      read *, rx
      read *, tmax
      read *, tol
      read *, istart
      read *, ittyp
      read *, itmax
!      
!--- Some constants
!      
      t = temp 
      errtol = tol
      irep = 1

      pwave=.false.
      dwave=.false.
      cyl=.false.
      if(icyl.eq.1) cyl=.true.
      pwave=.true.
!        
!--- Get the Ozaki parameters
!           
      open(10,file='ozaki.dat',status='old')
      read(10,*) Ncmax,temp,x
!     print *, Ncmax,temp,x
      NcCutof=int(x)
      if(abs(temp-t).gt.1e-8) print *, ' WRONG TEMP!!!!'

      do i=1,Ncmax
         read(10,*) zp(i),Rp(i)
      end do       
      close(10)
!           
!--- get the bulk gap
!        
      maxfrec=int(NcCutof/t+.00001)
      aa0=aa0/(1.0+aa0/3.0)
      dx=rx/float(nx)
      call bulkgap(t,deltat)
!
!---  Set the coupling constant
!
      dv=0.0
      do ien=1,Ncmax,1
         dv=dv+t*Rp(ien)/zp(ien)   ! must be the same as in gap equation
      end do
      esum=t/(log(t)+dv)
!        
!---  Initialise the order parameter
!
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

            if(istart.eq.-2) call bulkA(xx,xy,xz,yx,yy,yz,zx,zy,zz)
            if(istart.eq.-1) call bulkB(xx,xy,xz,yx,yy,yz,zx,zy,zz)
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
         end do
      else
         open(1,file='op_xyz',status='unknown')
         open(3,file='curr',status='unknown')
         do ix=0,nx,1
            read(1,*) x,a1x,b1x,a1y,b1y,a1z,b1z,a2x,b2x,a2y,b2y,a2z,b2z,a3x,b3x,a3y,b3y,a3z,b3z
            dxx(ix)=cmplx(a1x,b1x)
            dxy(ix)=cmplx(a1y,b1y)
            dxz(ix)=cmplx(a1z,b1z)
            dyx(ix)=cmplx(a2x,b2x)
            dyy(ix)=cmplx(a2y,b2y)
            dyz(ix)=cmplx(a2z,b2z)
            dzx(ix)=cmplx(a3x,b3x)
            dzy(ix)=cmplx(a3y,b3y)
            dzz(ix)=cmplx(a3z,b3z)
            read(3,*) x,x,a1x,a1y,a1z
            vx(ix)=cmplx(a1x,0.0)
            vy(ix)=cmplx(a1y,0.0)
            vz(ix)=cmplx(a1z,0.0)
         end do
         close(1)
         close(3)
      endif
!
!---  Pick the trajectories
!
      x = pi / tmax
      do i=1,tmax
         phi=float(i-1)*x+0.5*x
         kx(i)=cos(phi)
         ky(i)=sin(phi)
         awei(i)=1.0/tmax
      end do

      open(10,file='gauss11.dat',status='old')
      x=0.0
      do i=1,11
         read(10,*) kz(i),pwei(i)
         pwei(i)=pwei(i)/2.0
         x=x+pwei(i)
!        print 1300, i,kz(i),pwei(i),x
      end do
      close(10)
      x=0.0
      y=0.0
      z=0.0
      do i=1,tmax
         do j=1,11
            x=x+pwei(j)*awei(i)
            y=y+pwei(j)*awei(i)*kx(i)*kx(i)*(1.0-kz(j)*kz(j))
            z=z+pwei(j)*awei(i)*kz(j)*kz(j)
         end do
      end do
!      print *, 'quads =',x,y,z
!
!---  Go with this setup
!
      errtol=errtol*deltat
      print 5200, ' n-c(0)/a-c(1)/d-c(2)/old(3)    : ',istart
      print 5000, ' Temperature                    : ',t
      print 5000, ' Vorticity                      : ',vort
      print 5050, '     Vcs0                       : ',esum
      print 5050, '     As1                        : ',aa0
      print 5050, ' Delta(T)                       : ',deltat
      print 5000, ' Grid radius                    : ',rx
      print 5100, ' Accuracy                       : ',errtol
      print 5200, ' Nr. of trajs                   : ',tmax
      print 5200, ' Nr. of freqs                   : ',Ncmax
      print 5200, ' Ittyp (0=NN, 1=BR, 2=BB, AM=3) : ',ittyp
      print 5200, ' Max iterations                 : ',itmax

      open(1,file='op_xyz',status='unknown')
      open(2,file='op_harm',status='unknown')
      open(3,file='curr',status='unknown')

      do ix=0,nx,1
         x=float(ix)*dx

         cmm=0.5*(dxx(ix)-dyy(ix)+w*(dxy(ix)+dyx(ix)))
         cmo=sqrth*(dxz(ix)+w*dyz(ix))
         cmp=0.5*(dxx(ix)+dyy(ix)-w*(dxy(ix)-dyx(ix)))

         com=sqrth*(dzx(ix)+w*dzy(ix))
         coo=dzz(ix)
         cop=sqrth*(dzx(ix)-w*dzy(ix))

         cpm=0.5*(dxx(ix)+dyy(ix)+w*(dxy(ix)-dyx(ix)))
         cpo=sqrth*(dxz(ix)-w*dyz(ix))
         cpp=0.5*(dxx(ix)-dyy(ix)-w*(dxy(ix)+dyx(ix)))

         write(1,1000) x,dxx(ix),dxy(ix),dxz(ix),dyx(ix),dyy(ix),dyz(ix),dzx(ix),dzy(ix),dzz(ix)
         write(2,1000) x,cpp,cpo,cpm,cop,coo,com,cmp,cmo,cmm
         d=  abs(dxx(ix))**2+abs(dxy(ix))**2+abs(dxz(ix))**2
         d=d+abs(dyx(ix))**2+abs(dyy(ix))**2+abs(dyz(ix))**2
         d=d+abs(dzx(ix))**2+abs(dzy(ix))**2+abs(dzz(ix))**2
         d=sqrt(d/3.0)

         write(3,1000) x,d,real(vx(ix)),real(vy(ix)),real(vz(ix))
      end do
      close(1)
      close(2)
      close(3)
      close(10)
!
!--- Some format statements
!
 1000 format((1x,f8.3),30(1x,e14.6))
!1300 format(1x,i3,20(1x,e14.6))
 5000 format(a,f8.4)
 5050 format(a,f8.4)
 5100 format(a,e10.2)
 5200 format(a,i5)
!
!---------------------------------------------------------------------
!
   end subroutine init_calc
!
!---------------------------------------------------------------------
!
   subroutine dop(x,xx,xy,xz,yx,yy,yz,zx,zy,zz)
!
!---------------------------------------------------------------------
!
      real, intent(in) :: x
      complex, intent(out) :: xx, xy, xz, yx, yy, yz, zx, zy, zz

      real :: r, dc4, xgl, ar

      xgl=1.0/sqrt(1.0-t**2)
      r=sqrt(x*x)
      ar=r/xgl/3.0

      dc4=1.0
      if(r > 1.0e-12) dc4=tanh(ar)/ar

      xx=tanh(ar)
      yy=xx
      zz=xx
      yx=0.0
      xy=0.0
      zx=cmplx( dc4, 0.0)
      zy=cmplx( 0.0, 0.07*dc4)
      xz=cmplx(-dc4, 0.0)
      yz=cmplx( 0.0, 0.05*dc4)

      xx=deltat*xx
      xy=deltat*xy
      xz=deltat*xz
      yx=deltat*yx
      yy=deltat*yy
      yz=deltat*yz
      zx=deltat*zx
      zy=deltat*zy
      zz=deltat*zz
!
!---------------------------------------------------------------------
!
   end subroutine dop
!
!---------------------------------------------------------------------
!
   subroutine aop(x,xx,xy,xz,yx,yy,yz,zx,zy,zz)
!
!---------------------------------------------------------------------
!
      real, intent(in) :: x
      complex, intent(out) :: xx, xy, xz, yx, yy, yz, zx, zy, zz

      real :: r, ph, dc4, xgl, ar, c, s

      xgl=1.0/sqrt(1.0-t**2)
      r=sqrt(x**2)
      ar=r/xgl/3.0

      dc4=1.0
      if(r > 1.0e-12) dc4=tanh(ar)/ar

      xx=tanh(ar)
      ph=0.0

      yy=xx
      zz=xx
      yx=0.0
      xy=0.0

      s=sin(2.0*ph)
      c=cos(2.0*ph)
      rz=dc4*(c+1.0)/2.0
      az=dc4*s/2.0
      zx= cmplx(rz,az)
      xz=-zx

      rz=dc4*s/2.0
      az=dc4*(1.0-c)/2.0
      zy=cmplx(rz,az)
      yz=-zy

      xx=deltat*xx
      xy=deltat*xy
      xz=deltat*xz
      yx=deltat*yx
      yy=deltat*yy
      yz=deltat*yz
      zx=deltat*zx
      zy=deltat*zy
      zz=deltat*zz
!
!---------------------------------------------------------------------
!
   end subroutine aop
!
!---------------------------------------------------------------------
!
   subroutine nop(x,xx,xy,xz,yx,yy,yz,zx,zy,zz)
!
!---------------------------------------------------------------------
!
      real, intent(in) :: x
      complex, intent(out) :: xx, xy, xz, yx, yy, yz, zx, zy, zz

      real :: r, xgl, ar, cc

      xgl=1.0/sqrt(1.0-t**2)
      r=sqrt(x**2)
      ar=r/(3.*xgl)

      cc=1.0
      if(vort.gt.0.0) cc=tanh(ar)

      xx=deltat*cc
      xy=0.0
      xz=0.0
      yx=0.0
      yy=deltat*cc
      yz=0.0
      zx=0.0
      zy=0.0
      zz=deltat*cc
!
!---------------------------------------------------------------------
!
   end subroutine nop
!
!---------------------------------------------------------------------
!
   subroutine bulkB(xx,xy,xz,yx,yy,yz,zx,zy,zz)
!
!---------------------------------------------------------------------
!
      complex, intent(out) :: xx, xy, xz, yx, yy, yz, zx, zy, zz

      xx = cmplx(deltat,0.0)
      xy = 0.0
      xz = 0.0

      yx = 0.0
      yy = xx
      yz = 0.0

      zx = 0.0
      zy = 0.0
      zz = xx
!
!---------------------------------------------------------------------
!
   end subroutine bulkB
!
!---------------------------------------------------------------------
!
   subroutine bulkA(xx,xy,xz,yx,yy,yz,zx,zy,zz)
!
!---------------------------------------------------------------------
!
      complex, intent(out) :: xx, xy, xz, yx, yy, yz, zx, zy, zz
      real :: ar

      ar= sqrt2*deltat

      xx = 0.0
      xy = 0.0
      xz = 0.0

      yx = 0.0
      yy = 0.0
      yz = 0.0

      zx = cmplx(ar,0.0)
      zy = cmplx(0.0,ar)
      zz = 0.0
!
!--------------------------------------------------------------------
!
   end subroutine bulkA
!
!--------------------------------------------------------------------
!
end MODULE Initialisation
!
!---------------------------------------------------------------------

