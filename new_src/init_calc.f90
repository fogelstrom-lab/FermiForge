Module Initialisation
   use Global_variables
   use Bulkop
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
contains
   subroutine read_radial_input(temp,tol,gridM)
      ! One numeric value per record, with optional trailing comments.
      ! Legacy: 9/10 free records, 11 cylinder records (explicit radius).
      ! Extended: gridM,Rx after icyl; 11/12 free or 13 cylinder records.
      use, intrinsic :: iso_fortran_env, only : input_unit
      real, intent(out) :: temp,tol,gridM
      real :: values(13)
      integer :: count, status, offset, pos, k
      character(len=1024) :: line

      count=0
      do
         read(input_unit,'(A)',iostat=status) line
         if (status < 0) exit
         if (status /= 0) error stop 'Could not read radial input'
         pos=scan(line,'!#:')
         if (pos > 0) line=line(:pos-1)
         if (len_trim(line)==0) cycle
         count=count+1
         if (count > size(values)) error stop 'Too many radial input records'
         read(line,*,iostat=status) values(count)
         if (status /= 0) error stop 'Invalid numeric radial input record'
      end do
      if (count < 9) error stop 'Radial input requires at least nine records'
      if (.not. ieee_is_finite(values(3))) error stop 'Invalid icyl'
      if (values(3)/=0.0 .and. values(3)/=1.0) error stop 'icyl must be 0 or 1'
      icyl=int(values(3))
      cyl=icyl==1
      offset=0
      if (cyl) then
         if (count/=11 .and. count/=13) &
            error stop 'icyl=1 requires a cylinder radius after AA p_max (11 or 13 records)'
         if (count==13) offset=2
         if (.not. ieee_is_finite(values(count))) error stop 'Cylinder radius must be finite'
         if (values(count)<=0.0) error stop 'Cylinder radius must be positive'
      else
         if (count>12) error stop 'Free radial input requires 9,10,11 or 12 records'
         if (count>=11) offset=2
      end if
      if (.not. all(ieee_is_finite(values(:count)))) error stop 'Nonfinite radial input'
      gridM=70.0
      Rx=10.0
      if (offset==2) then
         gridM=values(4)
         Rx=values(5)
      end if
      if (gridM<=0.0 .or. Rx<=0.0) error stop 'Radial grid scales must be positive'
      temp=values(1)
      vort=values(2)
      aa0=values(4+offset)
      do k=1,5
         ! Check integer records before conversion, excluding tolerance.
         if (k==2) cycle
         pos=4+offset+k
         if (values(pos)/=real(int(values(pos)))) error stop 'Expected integer radial input'
      end do
      tmax=int(values(5+offset))
      tol=values(6+offset)
      istart=int(values(7+offset))
      ittyp=int(values(8+offset))
      itmax=int(values(9+offset))
      aa_pmax=1.0
      if (count>=10+offset) aa_pmax=values(10+offset)
      if (aa_pmax<0.01) error stop 'AA p_max must be at least 0.01'
      if (cyl) Rx=values(count)
   end subroutine read_radial_input

!--------------------------------------------------------------------------
!
   subroutine init_calc
!
!--------------------------------------------------------------------------
!
! --  Get the input and rig the calculation
!     
      integer :: i,j,ix,ien
      real ::  x,y,z,d,phi,dv,temp, tol, gridM, restart_radius, current_radius
      real ::  a1x,b1x,a1y,b1y,a1z,b1z
      real ::  a2x,b2x,a2y,b2y,a2z,b2z
      real ::  a3x,b3x,a3y,b3y,a3z,b3z
      complex :: xx,xy,xz,yx,yy,yz,zx,zy,zz
!
! -- input data given in qcv.inp
!
      call read_radial_input(temp,tol,gridM)
!      
!--- Some constants
!      
      tx = atan(gridM/Rx)/float(nx)
      if (cyl) then
         ! read_radial_input sets Rx to the explicit physical wall radius.
         gridM = Rx
      end if

      open(1,file='xgrid.dat',status='unknown')
      do i = 0, nx, 1
         if (cyl) then
            xgrid(i)=Rx*real(i)/real(nx)
         else
            xgrid(i)=Rx*tan(tx*i)
         end if
         write(1,*) i,xgrid(i)
      end do
      close(1)
      dx = 2.0*gridM/float(mx)
      if (cyl) dx = Rx/real(nx)

      t = temp 
      errtol = tol
      irep = 1

      pwave=.false.
      dwave=.false.
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
            x= xgrid(ix)

            if(istart.eq.-2) call bulkA(xx,xy,xz,yx,yy,yz,zx,zy,zz)
            if(istart.eq.-1) call bulkB(xx,xy,xz,yx,yy,yz,zx,zy,zz)
            if(istart.eq.0) call nc_op(x,xx,xy,xz,yx,yy,yz,zx,zy,zz)
            if(istart.eq.1) call ac_op(x,xx,xy,xz,yx,yy,yz,zx,zy,zz)
            if(istart.eq.2) call dc_op(x,xx,xy,xz,yx,yy,yz,zx,zy,zz)

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
            restart_radius = x
            if (cyl) then
               ! Existing text checkpoints round radii to three decimals.
               if (abs(restart_radius-xgrid(ix)) > 5.1e-4) &
                  error stop 'Cylinder restart order-parameter grid does not match radius'
            end if
            dxx(ix)=cmplx(a1x,b1x)
            dxy(ix)=cmplx(a1y,b1y)
            dxz(ix)=cmplx(a1z,b1z)
            dyx(ix)=cmplx(a2x,b2x)
            dyy(ix)=cmplx(a2y,b2y)
            dyz(ix)=cmplx(a2z,b2z)
            dzx(ix)=cmplx(a3x,b3x)
            dzy(ix)=cmplx(a3y,b3y)
            dzz(ix)=cmplx(a3z,b3z)
            read(3,*) current_radius,x,a1x,a1y,a1z
            if (cyl) then
               if (abs(current_radius-xgrid(ix)) > 5.1e-4) &
                  error stop 'Cylinder restart mean-field grid does not match radius'
            end if
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
!     errtol=errtol*deltat  ! we use relative errors ||dg||/||g||
!
      print 5200, ' n-c(0)/a-c(1)/d-c(2)/old(3)    : ',istart
      print 5000, ' Temperature                    : ',t
      print 5000, ' Vorticity                      : ',vort
      print 5050, '     Vcs0                       : ',esum
      print 5051, '     As1(Fs1)                   : ',aa0,'   (',aa0/(1.-aa0/3.0),')'
      print 5050, ' Delta(T)                       : ',deltat
      print 5000, ' Grid radius                    : ',xgrid(nx)
      if (cyl) then
         print 5000, ' Specular cylinder radius       : ',Rx
         print 5000, ' Reflected half-path length     : ',mx*dx
      end if
      print 5000, ' integration step               : ',dx
      print 5100, ' Accuracy                       : ',errtol
      print 5200, ' Nr. of trajs                   : ',tmax
      print 5200, ' Nr. of freqs                   : ',Ncmax
      print 5200, ' Ittyp (0=NN, 1=BR, 2=BB, AA=3) : ',ittyp
      print 5200, ' Max iterations                 : ',itmax
      print 5000, ' AA p_max                       : ',aa_pmax

      open(1,file='op_xyz',status='unknown')
      open(2,file='op_harm',status='unknown')
      open(3,file='curr',status='unknown')

      do ix=0,nx,1
         cmm(ix)=0.5*(dxx(ix)-dyy(ix)+w*(dxy(ix)+dyx(ix)))
         cmo(ix)=sqrth*(dxz(ix)+w*dyz(ix))
         cmp(ix)=0.5*(dxx(ix)+dyy(ix)-w*(dxy(ix)-dyx(ix)))

         com(ix)=sqrth*(dzx(ix)+w*dzy(ix))
         coo(ix)=dzz(ix)
         cop(ix)=sqrth*(dzx(ix)-w*dzy(ix))

         cpm(ix)=0.5*(dxx(ix)+dyy(ix)+w*(dxy(ix)-dyx(ix)))
         cpo(ix)=sqrth*(dxz(ix)-w*dyz(ix))
         cpp(ix)=0.5*(dxx(ix)-dyy(ix)-w*(dxy(ix)+dyx(ix)))
      end do

      do ix=0,nx,1
         x=xgrid(ix)

         write(1,1000) x,dxx(ix),dxy(ix),dxz(ix),dyx(ix),dyy(ix),dyz(ix),dzx(ix),dzy(ix),dzz(ix)
         write(2,1000) x,cpp(ix),cpo(ix),cpm(ix),cop(ix),coo(ix),com(ix),cmp(ix),cmo(ix),cmm(ix)

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
 5051 format(2(a,f8.4),a)
 5100 format(a,e10.2)
 5200 format(a,i5)
!
!---------------------------------------------------------------------
!
   end subroutine init_calc
!
!---------------------------------------------------------------------
!
   subroutine dc_op(x,xx,xy,xz,yx,yy,yz,zx,zy,zz)
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
   end subroutine dc_op
!
!---------------------------------------------------------------------
!
   subroutine ac_op(x,xx,xy,xz,yx,yy,yz,zx,zy,zz)
!
!---------------------------------------------------------------------
!
      real, intent(in) :: x
      complex, intent(out) :: xx, xy, xz, yx, yy, yz, zx, zy, zz

      real :: r, ph, dc4, xgl, ar, c, s, rz, az

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
   end subroutine ac_op
!
!---------------------------------------------------------------------
!
   subroutine nc_op(x,xx,xy,xz,yx,yy,yz,zx,zy,zz)
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
   end subroutine nc_op
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
