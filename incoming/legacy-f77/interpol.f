c---------------------------------------------------------------------
c
c     This is aset of interoplation routines and functions that
c     delivers the self energies (sem) along a specified trajectory (it,ip)
c     that setends from a given spatial point (iy)
c
c---------------------------------------------------------------------
c
      subroutine intord_v(it,ip,iy,sem)  ! free vortex
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'

      complex sem(4,-mx:mx),intses(10),dummy
      real x,y,x0,y0,r,ar,ph,sx,sy,sz,py,px,d,dd,arr,ar2,jproj
      complex intpol,mf,cph
      integer ix,iy,ii,it,ip,dix,diy

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
            sem(1,ix)=(cph+intses(1)*ar2)*sx
     +                    +intses(2)*ar2*sy
     +                    +intses(3)*arr*sz

            sem(2,ix)=     intses(4)*ar2*sx
     +               +(cph+intses(5)*ar2)*sy
     +                    +intses(6)*arr*sz

            sem(3,ix)=     intses(7)*arr*sx
     +                    +intses(8)*arr*sy
     +               +(cph+intses(9)*ar2)*sz

            sem(4,ix)=intses(10)*arr*jproj
         endif
      enddo
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine intord_c(it,ip,iy,sem)  ! vortex in cylinder
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'

      complex sem(4,-mx:mx),dummy,intses(10)
      real x,y,z,xo,yo,zo,r,ph,sx,sy,sz,py,px
      real ra,ra2,step,sd
      real bx(0:mx),by(0:mx),bz(0:mx)
      real fx(0:mx),fy(0:mx),fz(0:mx)
      real bpx(0:mx),bpy(0:mx),bpz(0:mx)
      real fpx(0:mx),fpy(0:mx),fpz(0:mx)
      complex intpol,mf,cph
      integer ix,iy,ii,it,ip,dix,diy,nsn

      sz=kz(ip)
      sd=sqrt(1.0-sz*sz)
      sx=kx(it)*sd
      sy=ky(it)*sd

      xo=dx*float(iy)
      zo=0.0
      yo=0.0
      ra=Rx

c     ra=Rx+dx*0.00001)  ! Seems to be a glitch in trajectories to treat
c                        ! a point on the radius of the cylinder
      if(iy.eq.nx) xo=xo-0.0000001*dx

      ra2=ra*ra
      step=dx
      nsn=mx

      call trajectories(xo,yo,zo,sx,sy,sz,ra,ra2,step,nsn,
     +               bx,by,bz,fx,fy,fz,bpx,bpy,bpz,fpx,fpy,fpz)
c
c steping along (px,py,pz)
c
      do ix=0,mx,1 
         x=fx(ix)
         y=fy(ix)
         z=fz(ix)
         r=sqrt(x*x+y*y)
         sx=fpx(ix)
         sy=fpy(ix)
         sz=fpz(ix)
c        if(r.ge.Rx) write(*,1000) '1:',iy,ix,r,sx,sy,sz
         px=0.0
         py=0.0
         ph=0.0
         if(r.gt.1e-6) then
            px=-y/r
            py= x/r
            ph=atan2(y,x)
         endif

         dummy=intpol(r,ph,intses,1)
         sem(1,ix)=intses(1)*sx+intses(2)*sy+intses(3)*sz
         sem(2,ix)=intses(4)*sx+intses(5)*sy+intses(6)*sz
         sem(3,ix)=intses(7)*sx+intses(8)*sy+intses(9)*sz
         sem(4,ix)=intses(10)*aa0*(px*sx+py*sy)
      enddo
c
c steping along (-px,-py,-pz)
c
      do ix=0,-mx,-1 
         ii=abs(ix)
         x=bx(ii)
         y=by(ii)
         z=bz(ii)
         r=sqrt(x*x+y*y)
         sx=bpx(ii)
         sy=bpy(ii)
         sz=bpz(ii)
c        if(r.ge.Rx) write(*,1000) '2:',iy,ix,r,sx,sy,sz
         px=0.0
         py=0.0
         ph=0.0
         if(r.gt.1e-6) then
            px=-y/r
            py= x/r
            ph=atan2(y,x)
         endif

         dummy=intpol(r,ph,intses,1)
         sem(1,ix)=intses(1)*sx+intses(2)*sy+intses(3)*sz
         sem(2,ix)=intses(4)*sx+intses(5)*sy+intses(6)*sz
         sem(3,ix)=intses(7)*sx+intses(8)*sy+intses(9)*sz
         sem(4,ix)=intses(10)*aa0*(px*sx+py*sy)
      enddo
 1000 format(a,2(x,i5),10(x,f9.3))
      return
      end
c
c---------------------------------------------------------------------
c
      complex function intpol(x,ph,ses,icomp)
c
c--------------------------------------------------------------------
c
      include 'qcv.dat'

      integer ix,iy,idx,i,ii,icomp
      real x,p,ph,pw(-1:1)
      complex cpp,cpo,cpm,cop,coo,com,cmp,cmo,cmm,svy
      complex ixx,ixy,ixz,iyx,iyy,iyz,izx,izy,izz
      complex c00,cm1,cp1,cp2,cm2,cone,cdiv,ses(10)
      complex q00,qm1,qp1,qp2,qm2
      parameter (cone=cmplx(1.0,0.0))

      intpol=0.0
      c00=exp(w*vort*ph)
      cm1=exp(-w*ph)
      cm2=cm1*cm1
      cp1=cdiv(cone,cm1)
      cp2=cp1*cp1

      if(icomp.eq.0) goto 10 ! Extrapolating for a free vortex

      iy=int(x/dx)
      ix=iy
      idx=1
      p=(x-dx*ix)/dx
c     if(x.lt.dx) p=(x-dx)/dx
c
c -- Take care of the endpoints
c
      if(iy.eq.nx) then
         ix=nx-1
         p=1.0
      endif
c     
c -- Second order iterpolation on an equidistant grid
c     
      pw(-1)=0.5*p*(p-1.0) 
      pw( 0)=1.0-p*p
      pw( 1)=0.5*p*(p+1.0)

      if(abs(p).gt.1) write(*,*) ' intpol: p ',p,ix,x


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
         qp1=cdiv(cone,qm1)
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

      return
      end
c
c---------------------------------------------------------------------
c
      subroutine newses
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'
      integer ix,iy,icomp,ex
      real errav,errmax,x
c
      integer icmp
      common /itercount/icmp
c
      ex=0
      errav=0.0
      errmax=0.0
      icomp=0
      do ix=0,nx,1
         x=abs(dxx(ix)-dnxx(ix))
         dxx(ix)=dnxx(ix)
         if(x.gt.errmax) then
            icomp=1
            errmax=x
            ex=ix
         endif
         errav=errav+x

         x=abs(dxy(ix)-dnxy(ix))
         dxy(ix)=dnxy(ix)
         if(x.gt.errmax) then
            icomp=2
            errmax=x
            ex=ix
         endif
         errav=errav+x

         x=abs(dxz(ix)-dnxz(ix))
         dxz(ix)=dnxz(ix)
         if(x.gt.errmax) then
            icomp=3
            errmax=x
            ex=ix
         endif
         errav=errav+x

         x=abs(dyx(ix)-dnyx(ix))
         dyx(ix)=dnyx(ix)
         if(x.gt.errmax) then
            icomp=4
            errmax=x
            ex=ix
         endif
         errav=errav+x

         x=abs(dyy(ix)-dnyy(ix))
         dyy(ix)=dnyy(ix)
         if(x.gt.errmax) then
            icomp=5
            errmax=x
            ex=ix
         endif
         errav=errav+x

         x=abs(dyz(ix)-dnyz(ix))
         dyz(ix)=dnyz(ix)
         if(x.gt.errmax) then
            icomp=6
            errmax=x
            ex=ix
         endif
         errav=errav+x

         x=abs(dzx(ix)-dnzx(ix))
         dzx(ix)=dnzx(ix)
         if(x.gt.errmax) then
            icomp=7
            errmax=x
            ex=ix
         endif
         errav=errav+x

         x=abs(dzy(ix)-dnzy(ix))
         dzy(ix)=dnzy(ix)
         if(x.gt.errmax) then
            icomp=8
            errmax=x
            ex=ix
         endif
         errav=errav+x

         x=abs(dzz(ix)-dnzz(ix))
         dzz(ix)=dnzz(ix)
         if(x.gt.errmax) then
            icomp=9
            errmax=x
            ex=ix
         endif
         errav=errav+x

         if(aa0.ne.0) then
            x=abs(vx(ix)-nvx(ix))
            if(x.gt.errmax) then
               icomp=10
               errmax=x
               ex=ix
            endif
            errav=errav+x

            x=abs(vy(ix)-nvy(ix))
            if(x.gt.errmax) then
               icomp=11
               errmax=x
               ex=ix
            endif
            errav=errav+x

            x=abs(vz(ix)-nvz(ix))
            if(x.gt.errmax) then
               icomp=12
               errmax=x
               ex=ix
            endif
            errav=errav+x
         endif

         vx(ix)=nvx(ix)
         vy(ix)=nvy(ix)
         vz(ix)=nvz(ix)
      enddo
      errav=errav/(nx+1)
      erreps=errmax
      write(30,1000) 
     + ' it nr',icmp,' errav = ',errav,' errmax = ',errmax,
     +   icomp,ex
      if(outerr) write(*,1000) 
     + ' it nr',icmp,' errav = ',errav,' errmax = ',errmax,
     +   icomp,ex
 1000 format(a,i4,a,e14.6,a,e14.6,2(x,i4))
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine dop(x,xx,xy,xz,yx,yy,yz,zx,zy,zz)
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'

      real x,r,ph,dc4,xgl,ar,c,s
      complex xx,xy,xz,yx,yy,yz,zx,zy,zz
                  
      xgl=1./sqrt(1-t**2)
      r=sqrt(x**2)
      ar=r/xgl/3.0
      if(r.eq.0) dc4=1.0
      if(r.ne.0) dc4=tanh(ar)/ar

      ph=0.0 !atan2(y,x)
      xx=exp(w*vort*ph)*tanh(ar)
      yy=xx
      zz=xx
      yx=0.0
      xy=0.0 
      zx=cmplx( dc4,0.0)
      zy=cmplx(0.0,0.07*dc4)
      xz=cmplx(-dc4,0.0)
      yz=cmplx(0.0,0.05*dc4)

      xx=deltat*xx
      xy=deltat*xy
      xz=deltat*xz
      yx=deltat*yx
      yy=deltat*yy
      yz=deltat*yz
      zx=deltat*zx
      zy=deltat*zy
      zz=deltat*zz
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine aop(x,xx,xy,xz,yx,yy,yz,zx,zy,zz)
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'

      real x,r,ph,dc4,xgl,ar,s,c,rz,az
      complex xx,xy,xz,yx,yy,yz,zx,zy,zz
                  
      xgl=1./sqrt(1.0-t**2)
      r=sqrt(x**2)
      ar=r/xgl/3.0
      if(r.eq.0) dc4=1.0
      if(r.ne.0) dc4=tanh(ar)/ar

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

 1000 format(30(x,e10.3))
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine nop(x,xx,xy,xz,yx,yy,yz,zx,zy,zz)
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'

      real x,r,xgl,ar,s,c,rz,az
      complex xx,xy,xz,yx,yy,yz,zx,zy,zz,cc
                  
      xgl=1./sqrt(1.0-t**2)
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

 1000 format(30(x,e10.3))
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine bulk(x,xx,xy,xz,yx,yy,yz,zx,zy,zz)
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'

      real x,r,ph,xgl,ar,s,c,rz,az
      complex xx,xy,xz,yx,yy,yz,zx,zy,zz,cc
                  
      xgl=1./sqrt(1.0-t**2)
      r=sqrt(x**2)
      ar=r/(3.*xgl)

      xx=deltat
      xy=0.0
      xz=0.0

      yx=0.0
      yy=deltat
      yz=0.0

      zx=0.0
      zy=0.0
      zz=deltat

 1000 format(30(x,e10.3))
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine bulkA(x,xx,xy,xz,yx,yy,yz,zx,zy,zz)
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'

      real x,r,ph,xgl,ar,s,c,rz,az
      complex xx,xy,xz,yx,yy,yz,zx,zy,zz,cc
                  
      xgl=1./sqrt(1.0-t**2)
      r=sqrt(x**2)
      ar=r/(3.*xgl)

      xx=0.0
      xy=0.0
      xz=0.0

      yx=0.0
      yy=0.0
      yz=0.0

      zx=sqrt2*deltat
      zy=w*sqrt2*deltat
      zz=0.0 

 1000 format(30(x,e10.3))
      return
      end
c
c--------------------------------------------------------------------
