c---------------------------------------------------------------------
c
      complex function intpol_new(x,ph,icomp)
c
c--------------------------------------------------------------------
c
      include 'qcv.dat'

      integer ix,iy,ixp,ixm,icomp
      real x,p,ph,pm,p0,pp
      complex cpp,cpo,cpm,cop,coo,com,cmp,cmo,cmm
      complex ixx,ixy,ixz,iyx,iyy,iyz,izx,izy,izz
      complex cph,c00,cm1,cp1,cp2,cm2,cone,cdiv
      parameter (cone=cmplx(1.0,0.0))

      intpol_new=0.0
      c00=exp(w*vort*ph)
      cm1=exp(-w*ph)
      cm2=cm1*cm1
      cp1=cdiv(cone,cm1)
      cp2=cp1*cp1

      iy=int(x/dx)
      ixp=1
      ix=iy
      ixm=-1
      p=(x-dx*ix)/dx
      if(x.lt.dx) p=(x-dx)/dx
c
c -- Take care of the endpoints
c
      if(iy.eq.nx) then
         ix=nx-1
         p=1.0
      elseif(iy.eq.0) then
         ix=1
      endif
c
c -- Second order iterpolation on an equidistant grid
c
      pm=0.5*p*(p-1.0) 
      p0=1.0-p*p
      pp=0.5*p*(p+1.0)

      if(abs(p).gt.1) write(*,*) ' intpol: p ',p,ix,x

      if(icomp.gt.0) then
         ixx=pm*dxx(ix+ixm)+p0*dxx(ix)+pp*dxx(ix+ixp)
         ixy=pm*dxy(ix+ixm)+p0*dxy(ix)+pp*dxy(ix+ixp)
         ixz=pm*dxz(ix+ixm)+p0*dxz(ix)+pp*dxz(ix+ixp)

         iyx=pm*dyx(ix+ixm)+p0*dyx(ix)+pp*dyx(ix+ixp)
         iyy=pm*dyy(ix+ixm)+p0*dyy(ix)+pp*dyy(ix+ixp)
         iyz=pm*dyz(ix+ixm)+p0*dyz(ix)+pp*dyz(ix+ixp)

         izx=pm*dzx(ix+ixm)+p0*dzx(ix)+pp*dzx(ix+ixp)
         izy=pm*dzy(ix+ixm)+p0*dzy(ix)+pp*dzy(ix+ixp)
         izz=pm*dzz(ix+ixm)+p0*dzz(ix)+pp*dzz(ix+ixp)

      else  ! For the eventuality that extrapolation is needed

         ixx=dxx(nx)-deltat
         ixy=dxy(nx)
         ixz=dxz(nx)

         iyx=dyx(nx)
         iyy=dyy(nx)-deltat
         iyz=dyz(nx)

         izx=dzx(nx)
         izy=dzy(nx)
         izz=dzz(nx)-deltat
      endif

      cmm=0.5*(ixx-iyy+w*(ixy+iyx))*cp2*c00
      cmo=sqrth*(ixz+w*iyz)*cp1*c00
      cmp=0.5*(ixx+iyy-w*(ixy-iyx))*c00

      com=sqrth*(izx+w*izy)*cp1*c00
      coo=izz*c00
      cop=sqrth*(izx-w*izy)*cm1*c00

      cpm=0.5*(ixx+iyy+w*(ixy-iyx))*c00
      cpo=sqrth*(ixz-w*iyz)*cm1*c00
      cpp=0.5*(ixx-iyy-w*(ixy+iyx))*cm2*c00

      if(icomp.eq. 1) intpol_new=  0.5*(cmp+cpm+cmm+cpp)
      if(icomp.eq. 2) intpol_new=w*0.5*(cmm-cpp+cmp-cpm)
      if(icomp.eq. 3) intpol_new=  (cmo+cpo)*sqrth

      if(icomp.eq. 4) intpol_new=w*0.5*(cpp-cmm-cmp+cpm)
      if(icomp.eq. 5) intpol_new=  0.5*(cmp+cpm-cmm-cpp)
      if(icomp.eq. 6) intpol_new=w*(cpo-cmo)*sqrth

      if(icomp.eq. 7) intpol_new=  (com+cop)*sqrth
      if(icomp.eq. 8) intpol_new=w*(cop-com)*sqrth
      if(icomp.eq. 9) intpol_new=coo

      if(icomp.eq.10) 
     +   intpol_new=pm*vy(ix+ixm)+p0*vy(ix)+pp*vy(ix+ixp)

      if(icomp.eq. -1) intpol_new=0.5*(cmp+cpm+cmm+cpp)
      if(icomp.eq. -5) intpol_new=0.5*(cmp+cpm-cmm-cpp)
      if(icomp.eq. -9) intpol_new=coo

 1000 format(a,3(x,i5),2(x,f7.3))
      return
      end
c
c---------------------------------------------------------------------
