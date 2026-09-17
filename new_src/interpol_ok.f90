Module Interpolations
   use Global_variables
contains
!---------------------------------------------------------------------
!
!     This is aset of interoplation routines and functions that
!     delivers the self energies (sem) along a specified trajectory (it,ip)
!     that setends from a given spatial point (iy)
!
!---------------------------------------------------------------------
!
   subroutine intord_v(it,ip,iy,sem)  ! free vortex
!
!---------------------------------------------------------------------
!
      integer, intent(in) :: it, ip, iy
      complex, dimension(4,-mx:mx), intent(out) :: sem
      complex, dimension(10) :: intses

      real :: x, y, x0, y0, r, ar, ph, sx, sy, sz, py, px, d, arr, ar2, jproj
      complex :: dummy, cph
      integer :: ix

      sz=kz(ip)
      d=sqrt(1.0-sz*sz)
      sx=kx(it)*d
      sy=ky(it)*d

      ar=dx*nx     ! if you change her you must change in intpol
      y0=0.0

      x0=dx*float(iy)

      do ix=-mx,mx,1 
         d=dx*float(ix)
         x=x0+sx*d
         y=y0+sy*d
         r=sqrt(x*x+y*y)
         px=0.0
         py=0.0
         ph=0.0

         if(r.gt.1.0e-6) then
            px=-y/r
            py= x/r
            ph=atan2(y,x)
         endif

         jproj=aa0*(px*sx+py*sy)

         if(r.lt.ar) then
            dummy=intpol(r,ph,intses,1)
            sem(1,ix)=intses(1)*sx+intses(2)*sy+intses(3)*sz
            sem(2,ix)=intses(4)*sx+intses(5)*sy+intses(6)*sz
            sem(3,ix)=intses(7)*sx+intses(8)*sy+intses(9)*sz
            sem(4,ix)=intses(10)*jproj
         else
            dummy=intpol(ar,ph,intses,0)
            cph=exp(w*vort*ph)*deltat
            arr=ar/r
            ar2=arr*arr
            sem(1,ix)=(cph+intses(1)*ar2)*sx + intses(2)*ar2*sy + intses(3)*arr*sz
            sem(2,ix)=intses(4)*ar2*sx + (cph+intses(5)*ar2)*sy + intses(6)*arr*sz
            sem(3,ix)=intses(7)*arr*sx + intses(8)*arr*sy + (cph+intses(9)*ar2)*sz

            sem(4,ix)=intses(10)*arr*jproj
         endif
      enddo
!
!---------------------------------------------------------------------
!
   end subroutine intord_v
!
!---------------------------------------------------------------------
!
   complex function intpol(x,ph,ses,icomp)
!
!---------------------------------------------------------------------
!
      integer, intent(in) :: icomp
      real, intent(in) :: x, ph
      complex, dimension(10), intent(out) ::ses

      integer :: ix, iy, idx, i, ii
      real, dimension(-1:1) :: pw
      real :: p
      complex :: cpp, cpo, cpm, cop, coo, com, cmp, cmo, cmm, svy
      complex :: ixx, ixy, ixz, iyx, iyy, iyz, izx, izy, izz
      complex :: c00, cm1, cp1, cp2, cm2
      complex :: q00, qm1, qp1, qp2, qm2

      intpol=0.0
      c00=exp(w*vort*ph)
      cm1=exp(-w*ph)
      cm2=cm1*cm1
      cp1=1.0/cm1
      cp2=cp1*cp1

      if(icomp.eq.0) goto 10 ! Extrapolating for a free vortex

      iy=int(x/dx)
      ix=iy
      idx=1
      p=(x-dx*ix)/dx
!
! -- Take care of the endpoints
!
      if(iy.eq.nx) then
         ix=nx-1
         p=1.0
      endif
!     
! -- Second order iterpolation on an equidistant grid
!     
      pw(-1)=0.5*p*(p-1.0) 
      pw( 0)=1.0-p*p
      pw( 1)=0.5*p*(p+1.0)

      if(abs(p).gt.1) print *, ' intpol: p ',p,ix,x

      cmm=0.0
      cmo=0.0
      cmp=0.0

      com=0.0
      coo=0.0
      cop=0.0

      cpm=0.0
      cpo=0.0
      cpp=0.0

      svy=0.0
      if(iy.gt.0) then
         do ii=-1,1
            ixx=pw(ii)*dxx(ix+ii)
            ixy=pw(ii)*dxy(ix+ii)
            ixz=pw(ii)*dxz(ix+ii)

            iyx=pw(ii)*dyx(ix+ii)
            iyy=pw(ii)*dyy(ix+ii)
            iyz=pw(ii)*dyz(ix+ii)

            izx=pw(ii)*dzx(ix+ii)
            izy=pw(ii)*dzy(ix+ii)
            izz=pw(ii)*dzz(ix+ii)

            cmm=cmm+0.5*(ixx-iyy+w*(ixy+iyx))*cp2*c00
            cmo=cmo+sqrth*(ixz+w*iyz)*cp1*c00
            cmp=cmp+0.5*(ixx+iyy-w*(ixy-iyx))*c00

            com=com+sqrth*(izx+w*izy)*cp1*c00
            coo=coo+izz*c00
            cop=cop+sqrth*(izx-w*izy)*cm1*c00

            cpm=cpm+0.5*(ixx+iyy+w*(ixy-iyx))*c00
            cpo=cpo+sqrth*(ixz-w*iyz)*cm1*c00
            cpp=cpp+0.5*(ixx-iyy-w*(ixy+iyx))*cm2*c00

            svy=svy+pw(ii)*vy(ix+ii)
         enddo
      else
         q00=exp(w*vort*(ph+pi))
         qm1=exp(-w*(ph+pi))
         qm2=qm1*qm1
         qp1=1.0/qm1
         qp2=qp1*qp1

         do ii=-1,1
            i=abs(ii)
            ixx=pw(ii)*dxx(ix+i)
            ixy=pw(ii)*dxy(ix+i)
            ixz=pw(ii)*dxz(ix+i)

            iyx=pw(ii)*dyx(ix+i)
            iyy=pw(ii)*dyy(ix+i)
            iyz=pw(ii)*dyz(ix+i)

            izx=pw(ii)*dzx(ix+i)
            izy=pw(ii)*dzy(ix+i)
            izz=pw(ii)*dzz(ix+i)

            if(ii.lt.0) then
               cmm=cmm+0.5*(ixx-iyy+w*(ixy+iyx))*qp2*q00
               cmo=cmo+sqrth*(ixz+w*iyz)*qp1*q00
               cmp=cmp+0.5*(ixx+iyy-w*(ixy-iyx))*q00

               com=com+sqrth*(izx+w*izy)*qp1*q00
               coo=coo+izz*q00
               cop=cop+sqrth*(izx-w*izy)*qm1*q00

               cpm=cpm+0.5*(ixx+iyy+w*(ixy-iyx))*q00
               cpo=cpo+sqrth*(ixz-w*iyz)*qm1*q00
               cpp=cpp+0.5*(ixx-iyy-w*(ixy+iyx))*qm2*q00

               svy=svy-pw(ii)*vy(ix+i)
            else
               cmm=cmm+0.5*(ixx-iyy+w*(ixy+iyx))*cp2*c00
               cmo=cmo+sqrth*(ixz+w*iyz)*cp1*c00
               cmp=cmp+0.5*(ixx+iyy-w*(ixy-iyx))*c00

               com=com+sqrth*(izx+w*izy)*cp1*c00
               coo=coo+izz*c00
               cop=cop+sqrth*(izx-w*izy)*cm1*c00

               cpm=cpm+0.5*(ixx+iyy+w*(ixy-iyx))*c00
               cpo=cpo+sqrth*(ixz-w*iyz)*cm1*c00
               cpp=cpp+0.5*(ixx-iyy-w*(ixy+iyx))*cm2*c00
               svy=svy+pw(ii)*vy(ix+i)
            endif
         enddo
      endif

      ses(1)=  0.5*(cmp+cpm+cmm+cpp)
      ses(2)=w*0.5*(cmm-cpp+cmp-cpm)
      ses(3)=  (cmo+cpo)*sqrth

      ses(4)=w*0.5*(cpp-cmm-cmp+cpm)
      ses(5)=  0.5*(cmp+cpm-cmm-cpp)
      ses(6)=w*(cpo-cmo)*sqrth

      ses(7)=  (com+cop)*sqrth
      ses(8)=w*(cop-com)*sqrth
      ses(9)=   coo

      ses(10)=svy

 10   if(icomp.eq.0) then
         ixx=dxx(nx)-deltat
         ixy=dxy(nx)
         ixz=dxz(nx)

         iyx=dyx(nx)
         iyy=dyy(nx)-deltat
         iyz=dyz(nx)

         izx=dzx(nx)
         izy=dzy(nx)
         izz=dzz(nx)-deltat

         cmm=0.5*(ixx-iyy+w*(ixy+iyx))*cp2*c00
         cmo=sqrth*(ixz+w*iyz)*cp1*c00
         cmp=0.5*(ixx+iyy-w*(ixy-iyx))*c00

         com=sqrth*(izx+w*izy)*cp1*c00
         coo=izz*c00
         cop=sqrth*(izx-w*izy)*cm1*c00

         cpm=0.5*(ixx+iyy+w*(ixy-iyx))*c00
         cpo=sqrth*(ixz-w*iyz)*cm1*c00
         cpp=0.5*(ixx-iyy-w*(ixy+iyx))*cm2*c00

         ses(1)=  0.5*(cmp+cpm+cmm+cpp)
         ses(2)=w*0.5*(cmm-cpp+cmp-cpm)
         ses(3)=  (cmo+cpo)*sqrth

         ses(4)=w*0.5*(cpp-cmm-cmp+cpm)
         ses(5)=  0.5*(cmp+cpm-cmm-cpp)
         ses(6)=w*(cpo-cmo)*sqrth

         ses(7)=  (com+cop)*sqrth
         ses(8)=w*(cop-com)*sqrth
         ses(9)=   coo

         ses(10)=vy(nx)
      endif
!
!---------------------------------------------------------------------
!
   end function intpol
!
!---------------------------------------------------------------------
!
end MODULE Interpolations
!
!--------------------------------------------------------------------
