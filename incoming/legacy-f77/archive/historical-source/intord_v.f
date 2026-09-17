c---------------------------------------------------------------------
c
      subroutine intord_v(it,ip,iy,sem)
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'

      complex sem(4,-mx:mx)
      real x,y,x0,y0,r,ar,a,ph,sx,sy,sz,py,px,d,arr,ar2
      complex intpol,mf,cph
      integer ix,iy,ii,it,ip,dix,diy

      a=aa0/(1.+aa0/3.)
     
      sz=kz(ip)
      d=sqrt(1.0-sz*sz)
      sx=kx(it)*d
      sy=ky(it)*d

      ar=dx*(nx-2.0)
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
         if(r.gt.1e-6) then
            px=-y/r
            py= x/r
            ph=atan2(y,x)
         endif

         if(r.le.ar) then
            sem(1,ix)=intpol(r,ph,1)*sx
     +               +intpol(r,ph,2)*sy
     +               +intpol(r,ph,3)*sz

            sem(2,ix)=intpol(r,ph,4)*sx
     +               +intpol(r,ph,5)*sy
     +               +intpol(r,ph,6)*sz

            sem(3,ix)=intpol(r,ph,7)*sx
     +               +intpol(r,ph,8)*sy
     +               +intpol(r,ph,9)*sz

            sem(4,ix)=intpol(r,ph,10)*a*(px*sx+py*sx)
         else
            cph=exp(w*vort*ph)*deltat
            arr=ar/r
            ar2=arr*arr
            sem(1,ix)=(cph+intpol(ar,ph,-1)*ar2)*sx
     +                    +intpol(ar,ph, 2)*ar2*sy
     +                    +intpol(ar,ph, 3)*arr*sz

            sem(2,ix)=    intpol(ar,ph, 4)*ar2*sx
     +              +(cph+intpol(ar,ph,-5)*ar2)*sy
     +                   +intpol(ar,ph, 6)*arr*sz

            sem(3,ix)=    intpol(ar,ph, 7)*arr*sx
     +                   +intpol(ar,ph, 8)*arr*sy
     +              +(cph+intpol(ar,ph,-9)*ar2)*sz

            mf=intpol(ar,ph,10)*arr*a
            sem(4,ix)=intpol(ar,ph,10)*arr*a*(px*sx+py*sy)
         endif
      enddo
      return
      end
c
c---------------------------------------------------------------------
