c---------------------------------------------------------------------
c
      subroutine dop(x,y,xx,xy,xz,yx,yy,yz,zx,zy,zz)
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'

      real x,y,r,ph,dc4,xgl,ar,c,s
      complex xx,xy,xz,yx,yy,yz,zx,zy,zz
                  
      xgl=(1.0+aa0)/sqrt(1-t**2)
      r=sqrt(x**2+y**2)
c     ar=r/xgl/6.0
      ar=r/xgl/3.0  ! original version
      if(r.eq.0) dc4=1.0
      if(r.ne.0) dc4=tanh(ar)/ar

      ph=0.0
      if(y.eq.0.0 .and. x.eq.0.0) then
         xx=0.0
      else
         ph=atan2(y,x)
         xx=exp(w*ph)*tanh(ar)
      endif   
      yy=xx
      zz=xx
      yx=0.0
      xy=0.0 
      zx=cmplx( dc4,0.0)*cos(ph)**2
      xz=cmplx(-dc4,0.0)*cos(ph)**2
      zy=cmplx(0.0,0.0)
      yz=cmplx(0.0,0.0)

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
      subroutine aop(x,y,xx,xy,xz,yx,yy,yz,zx,zy,zz)
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'

      real x,y,r,ph,dc4,xgl,ar,s,c,rz,az
      complex xx,xy,xz,yx,yy,yz,zx,zy,zz
                  
      xgl=1./sqrt(1.0-t**2)
      r=sqrt(x**2+y**2)
      ar=r/xgl/3.0
      if(r.eq.0) dc4=1.0
      if(r.ne.0) dc4=tanh(ar)/ar

      if(y.eq.0.0 .and. x.eq.0.0) then
         ph=0.0
         xx=0.0
      else
         ph=atan2(y,x)
         xx=exp(w*ph)*tanh(ar)
      endif   

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
      subroutine nop(x,y,xx,xy,xz,yx,yy,yz,zx,zy,zz)
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'

      real x,y,r,ph,xgl,ar,s,c,rz,az
      complex xx,xy,xz,yx,yy,yz,zx,zy,zz,cc
                  
      xgl=1./sqrt(1.0-t**2)
      r=sqrt(x**2+y**2)
      ar=r/(3.*xgl)

      if(y.eq.0.0 .and. x.eq.0.0) then
         ph=0.0
         cc=0.0
      else
         ph=atan2(y,x)
         cc=exp(w*ph)*tanh(ar)
      endif   

      xx=deltat*cc
      yy=deltat*cc
      zz=deltat*cc

 1000 format(30(x,e10.3))
      return
      end
c
c--------------------------------------------------------------------
