c---------------------------------------------------------------------
c
      subroutine greenR(it,ip,iy,tem,gamma10,gamma11,gamma12,gamma13,
     +                               gamma20,gamma21,gamma22,gamma23)
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'
      complex p0,p1,p2,p3,rh,ch,ga,gb,g,ngg,zintegrand
      complex normp(2,2),normm(2,2),ggt(2,2),gtg(2,2)
      complex tem(0:15)
      complex gamma10,gamma11,gamma12,gamma13
      complex gamma20,gamma21,gamma22,gamma23
      parameter (rh=cmplx(0.25,0.0),ch=cmplx(0.0,-0.25))
      real xd,px,py,pz,v,tpx,tpy,tpz,vpx,vpy,vpz,tlnt
      integer i,iy,it,ip
c
c-- Local varables
c
      complex sp(2,2),gp(2,2)
      complex Gg(2,2),Gt(2,2),prop(4,4)
c
c -- Make the matrix gammas
c
      p0=gamma10
      p1=gamma11
      p2=gamma12
      p3=gamma13

      Gg(1,1)=-p1+w*p2
      Gg(2,1)=-p0+  p3
      Gg(1,2)= p0+  p3
      Gg(2,2)= p1+w*p2
      ga=p0

      p0=gamma20
      p1=gamma21
      p2=gamma22
      p3=gamma23

      Gt(1,1)=-p1-w*p2
      Gt(2,1)=-p0+  p3
      Gt(1,2)= p0+  p3
      Gt(2,2)= p1-w*p2
      gb=p0

      g=-2.0*w*ga/(1.0+ga*gb)

      call AtimesB(ggt,Gg,Gt)
      ngg=(1.0-ggt(1,1))*(1.0-ggt(2,2))-ggt(2,1)*ggt(1,2)
      normp(1,1)=(1.0-ggt(2,2))/ngg
      normp(1,2)= ggt(1,2)/ngg
      normp(2,1)= ggt(2,1)/ngg
      normp(2,2)=(1.0-ggt(1,1))/ngg

      call AtimesB(gtg,Gt,Gg)
      ngg=(1.0-gtg(1,1))*(1.0-gtg(2,2))-gtg(2,1)*gtg(1,2)
      normm(1,1)=(1.0-gtg(2,2))/ngg
      normm(1,2)= gtg(1,2)/ngg
      normm(2,1)= gtg(2,1)/ngg
      normm(2,2)=(1.0-gtg(1,1))/ngg
c
c -- upper diagonal g
c
      sp(1,1)=-w*(1.0+ggt(1,1))
      sp(1,2)=-w*ggt(1,2)
      sp(2,1)=-w*ggt(2,1)
      sp(2,2)=-w*(1.0+ggt(2,2))
      call AtimesB(gp,normp,sp)
      prop(1,1)=gp(1,1)
      prop(1,2)=gp(1,2)
      prop(2,1)=gp(2,1)
      prop(2,2)=gp(2,2)
c
c --  f
c
      sp(1,1)=-2.*w*gg(1,1)
      sp(1,2)=-2.*w*gg(1,2)
      sp(2,1)=-2.*w*gg(2,1)
      sp(2,2)=-2.*w*gg(2,2)
      call AtimesB(gp,normp,sp)
      prop(1,3)=gp(1,1)
      prop(1,4)=gp(1,2)
      prop(2,3)=gp(2,1)
      prop(2,4)=gp(2,2)
c
c --  f-tilde
c
      sp(1,1)= 2.*w*gt(1,1)
      sp(1,2)= 2.*w*gt(1,2)
      sp(2,1)= 2.*w*gt(2,1)
      sp(2,2)= 2.*w*gt(2,2)
      call AtimesB(gp,normm,sp)
      prop(3,1)=gp(1,1)
      prop(3,2)=gp(1,2)
      prop(4,1)=gp(2,1)
      prop(4,2)=gp(2,2)
c
c -- lower diagonal g
c
      sp(1,1)= w*(1.0+gtg(1,1))
      sp(1,2)= w*gtg(1,2)
      sp(2,1)= w*gtg(2,1)
      sp(2,2)= w*(1.0+gtg(2,2))
      call AtimesB(gp,normm,sp)
      prop(3,3)=gp(1,1)
      prop(3,4)=gp(1,2)
      prop(4,3)=gp(2,1)
      prop(4,4)=gp(2,2)
c
c--  Spin-resolved DOS
c
      tem(0)=  0.5*(prop(1,1)+prop(2,2))
      tem(1)=  0.5*(prop(1,2)+prop(2,1))
      tem(2)=w*0.5*(prop(1,2)-prop(2,1))
      tem(3)=  0.5*(prop(1,1)-prop(2,2))

 1000 format(30(x,e10.3))
 1100 format(a,2(x,i3),30(x,e10.3))
 666  return
      end
c
c----------------------------------------------------------------------
