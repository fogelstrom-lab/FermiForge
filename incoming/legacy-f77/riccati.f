c--------------------------------------------------------------------
c
      subroutine makeprops(it,ip,iy,ien,en,sem,tem)
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'
      integer ndx,nn,it,ip,i,iy,ien,cstep
      parameter (ndx=4,cstep=20)
      real px,py,pz,d,dd,sd
      complex lam,cdiv,gg0
      complex en,sem(4,-mx:mx),tem(12)
      complex zen1,zen2,cen,zex1(3),zex2(3)
      complex del11(0:3),del12(0:3),del21(0:3),del22(0:3)
      complex gold(0:3),gnew(0:3)
      complex gamma10,gamma11,gamma12,gamma13
      complex gamma20,gamma21,gamma22,gamma23
c
      cen=w*en
      pz=kz(ip)
      sd=sqrt(1.0-pz*pz)
      px=kx(it)*sd
      py=ky(it)*sd
c
c-----------------  gamma1 ------------------------------------
c
      call bulk_gammas(gold,sem,cen,px,py,pz,-mx,iy,1)
      d=10*dx
      nn=int(cstep*ndx)
      call setenergy1(-mx,iy,px,py,pz,cen,del11,del21,zen1,zex1,sem)
      call ricc(gold,gnew,del11,del11,del21,del21,
     +                          zen1,zen1,zex1,zex1,d,nn)
      gold(0)=gnew(0)
      gold(1)=gnew(1)
      gold(2)=gnew(2)
      gold(3)=gnew(3)

      d= dx
      nn=int(ndx)*(1+int(abs(en)))
      do i=-mx+1,0,1
         call setenergy1(i,iy,px,py,pz,cen,del12,del22,zen2,zex2,sem)
         call ricc(gold,gnew,del11,del12,del21,del22,
     +                             zen1,zen2,zex1,zex2,d,nn)

         del11=del12
         del21=del22
         zen1=zen2
         zex1=zex2
         gold(0)=gnew(0)
         gold(1)=gnew(1)
         gold(2)=gnew(2)
         gold(3)=gnew(3)
      enddo

      gamma10=gnew(0)
      gamma11=gnew(1)
      gamma12=gnew(2)
      gamma13=gnew(3)
c
c-----------------  gamma2 ------------------------------------
c
      call bulk_gammas(gold,sem,cen,px,py,pz,mx,iy,2)
      d=10*dx
      nn=int(cstep*ndx)
      call setenergy2(mx,iy,px,py,pz,cen,del11,del21,zen1,zex1,sem)
      call ricc(gold,gnew,del11,del11,del21,del21,
     +                    zen1,zen1,zex1,zex1,d,nn)
      gold(0)=gnew(0)
      gold(1)=gnew(1)
      gold(2)=gnew(2)
      gold(3)=gnew(3)

      d= dx
      nn=int(ndx)*(1+int(abs(en)))
      do i=mx-1,0,-1    
         call setenergy2(i,iy,px,py,pz,cen,del12,del22,zen2,zex2,sem)
         call ricc(gold,gnew,del11,del12,del21,del22,
     +                             zen1,zen2,zex1,zex2,d,nn)

         del11=del12
         del21=del22
         zen1=zen2
         zex1=zex2
         gold(0)=gnew(0)
         gold(1)=gnew(1)
         gold(2)=gnew(2)
         gold(3)=gnew(3)
      enddo

      gamma20=gnew(0)
      gamma21=gnew(1)
      gamma22=gnew(2)
      gamma23=gnew(3)

      call green(it,ip,iy,ien,tem,gamma10,gamma11,gamma12,gamma13,
     +                            gamma20,gamma21,gamma22,gamma23)

 1000 format(a,4(x,i3),30(x,e11.4))
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine ricc(cold,cnew,ca1,ca2,cb1,cb2,
     +                          cep1,cep2,cex1,cex2,d,nsteps)
c
c-----------------------------------------------------------------------
c
c  solving i*y'=-D_2*y^2-2*E*y-D_1 with linearly varying
c  D_1 (ca1...ca2), D_2 (cb1...cb2), and E (cep1...cep2).
c  y evolves from cold to cnew in the intervall d.
c  uses nsteps steps.
c
c-----------------------------------------------------------------------
c
        implicit none
c
        complex cold(0:3),cnew(0:3)
	complex ca1(0:3),ca2(0:3),cb1(0:3),cb2(0:3)
	complex cep1,cep2
        complex cex1(3),cex2(3)

        real d
        integer nsteps
c
c       **** Functions ****
c
        intrinsic idint,dmax1,dcmplx,dfloat,dimag,dble,dabs
c
c       **** Temporary variables ****
c
	real  da1(2,0:3),da2(2,0:3),db1(2,0:3),db2(2,0:3)
        complex c1a1(0:3),c1a2(0:3),c1b1(0:3),c1b2(0:3)
        complex cca1(0:3),cca2(0:3),cca3(0:3)
	complex ccb1(0:3),ccb2(0:3),ccb3(0:3)
        complex cdela(0:3),cdelb(0:3)

	real  de1(2),de2(2)
	complex c1e1,c1e2,cce1,cce2,cce3,cdele

        real  dx1(2,3),dx2(2,3)
        complex c1x1(3),c1x2(3),ccx1(3),ccx2(3),ccx3(3),cdelx(3)

        complex cy1(0:3),ck1(0:3),cy2(0:3),ck2(0:3),ccy1(0:3)
	complex cy3(0:3),ck3(0:3),cy4(0:3),ck4(0:3),ccy3(0:3)

	complex cp0,cp1,cp2,cp3,csum1,csum2,cadd,cexpr,czero
        complex cvp,csp
        real six,h0,dn2,dn,dh,dh2
        parameter (six=0.16666666666666666667,czero=cmplx(0.0,1.0))
        integer i,ii
c
c       **** Commons and equivalences ****
c
	EQUIVALENCE (c1a1(0),da1(1,0))
	EQUIVALENCE (c1a2(0),da2(1,0))
	EQUIVALENCE (c1b1(0),db1(1,0))
	EQUIVALENCE (c1b2(0),db2(1,0))
	EQUIVALENCE (c1e1,de1(1))
	EQUIVALENCE (c1e2,de2(1))
        EQUIVALENCE (c1x1(1),dx1(1,1))
        EQUIVALENCE (c1x2(1),dx2(1,1))
c
c-----------------------------------------------------------------------
c
c       **** statements ****
c
	h0=d
	dn2=1.0/dfloat(nsteps+nsteps)
	dn=dn2+dn2
	dh=h0*dn
	dh2=dh+dh
c
c -- Diagonal S.E. multiplied by i and 2
c
	c1e1=cep1
	c1e2=cep2
	cce1=cmplx(-de1(2),de1(1))*dh2
	cce2=cmplx(-de2(2),de2(1))*dh2
	cdele=(cce2-cce1)*dn2
	cce3=cce1
c
c --  Exchange field S.E multiplied by i and 2
c
        do i=1,3
           c1x1(i)=cex1(i)
           c1x2(i)=cex2(i)
           ccx1(i)=cmplx(-dx1(2,i),dx1(1,i))*dh2
           ccx2(i)=cmplx(-dx2(2,i),dx2(1,i))*dh2
           cdelx(i)=(ccx2(i)-ccx1(i))*dn2
           ccx3(i)=ccx1(i)
        enddo
c
c -- Order parameter fields multiplied by i
c
	do i=0,3
           c1a1(i)=ca1(i)
	   c1a2(i)=ca2(i)
	   cca1(i)=cmplx(-da1(2,i),da1(1,i))*dh
	   cca2(i)=cmplx(-da2(2,i),da2(1,i))*dh
	   cca3(i)=cca1(i)
	   cdela(i)=(cca2(i)-cca1(i))*dn2

	   c1b1(i)=cb1(i)
	   c1b2(i)=cb2(i)
	   ccb1(i)=cmplx(-db1(2,i),db1(1,i))*dh
	   ccb2(i)=cmplx(-db2(2,i),db2(1,i))*dh
	   ccb3(i)=ccb1(i)
	   cdelb(i)=(ccb2(i)-ccb1(i))*dn2

	   ccy3(i)=cold(i)
	enddo

	do i=1,nsteps,1

           cce1=cce3
           cce2=cce1+cdele
           cce3=cce2+cdele

	   do ii=0,3
              cca1(ii)=cca3(ii)
              cca2(ii)=cca1(ii)+cdela(ii)
              cca3(ii)=cca2(ii)+cdela(ii)

              ccb1(ii)=ccb3(ii)
              ccb2(ii)=ccb1(ii)+cdelb(ii)
              ccb3(ii)=ccb2(ii)+cdelb(ii)

              ccy1(ii)=ccy3(ii)
              cy1(ii)=ccy1(ii)
	   enddo

	   do ii=1,3
              ccx1(ii)=ccx3(ii)
              ccx2(ii)=ccx1(ii)+cdelx(ii)
              ccx3(ii)=ccx2(ii)+cdelx(ii)
	   enddo

	   csp=ccb1(0)*cy1(0)
	   cvp=ccb1(1)*cy1(1)+ccb1(2)*cy1(2)+ccb1(3)*cy1(3)
	   csum1=cce1+2.0*cvp
	   csum2=csum1+2.0*csp
	   csp=cy1(0)*cy1(0)
	   cvp=cy1(1)*cy1(1)+cy1(2)*cy1(2)+cy1(3)*cy1(3)

           cp0=ccx1(1)*cy1(1)+ccx1(2)*cy1(2)+ccx1(3)*cy1(3)
           cp1=ccx1(1)*cy1(0)
           cp2=ccx1(2)*cy1(0)
           cp3=ccx1(3)*cy1(0)

	   ck1(0)=(csp+cvp)*ccb1(0)+cy1(0)*csum1+cca1(0)-cp0
           ck1(1)=(csp-cvp)*ccb1(1)+cy1(1)*csum2+cca1(1)-cp1
           ck1(2)=(csp-cvp)*ccb1(2)+cy1(2)*csum2+cca1(2)-cp2
           ck1(3)=(csp-cvp)*ccb1(3)+cy1(3)*csum2+cca1(3)-cp3

           cy2(0)=ccy1(0)+ck1(0)*0.5
           cy2(1)=ccy1(1)+ck1(1)*0.5
           cy2(2)=ccy1(2)+ck1(2)*0.5
           cy2(3)=ccy1(3)+ck1(3)*0.5

	   csp=ccb2(0)*cy2(0)
	   cvp=ccb2(1)*cy2(1)+ccb2(2)*cy2(2)+ccb2(3)*cy2(3)
	   csum1=cce2+2.0*cvp
	   csum2=csum1+2.0*csp
	   csp=cy2(0)*cy2(0)
	   cvp=cy2(1)*cy2(1)+cy2(2)*cy2(2)+cy2(3)*cy2(3)

           cp0=ccx2(1)*cy2(1)+ccx2(2)*cy2(2)+ccx2(3)*cy2(3)
           cp1=ccx2(1)*cy2(0)
           cp2=ccx2(2)*cy2(0)
           cp3=ccx2(3)*cy2(0)

	   ck2(0)=(csp+cvp)*ccb2(0)+cy2(0)*csum1+cca2(0)-cp0
           ck2(1)=(csp-cvp)*ccb2(1)+cy2(1)*csum2+cca2(1)-cp1
           ck2(2)=(csp-cvp)*ccb2(2)+cy2(2)*csum2+cca2(2)-cp2
           ck2(3)=(csp-cvp)*ccb2(3)+cy2(3)*csum2+cca2(3)-cp3

           cy3(0)=ccy1(0)+ck2(0)*0.5
           cy3(1)=ccy1(1)+ck2(1)*0.5
           cy3(2)=ccy1(2)+ck2(2)*0.5
           cy3(3)=ccy1(3)+ck2(3)*0.5

	   csp=ccb2(0)*cy3(0)
	   cvp=ccb2(1)*cy3(1)+ccb2(2)*cy3(2)+ccb2(3)*cy3(3)
	   csum1=cce2+2.0*cvp
	   csum2=csum1+2.0*csp
	   csp=cy3(0)*cy3(0)
	   cvp=cy3(1)*cy3(1)+cy3(2)*cy3(2)+cy3(3)*cy3(3)

           cp0=ccx2(1)*cy3(1)+ccx2(2)*cy3(2)+ccx2(3)*cy3(3)
           cp1=ccx2(1)*cy3(0)
           cp2=ccx2(2)*cy3(0)
           cp3=ccx2(3)*cy3(0)

	   ck3(0)=(csp+cvp)*ccb2(0)+cy3(0)*csum1+cca2(0)-cp0
           ck3(1)=(csp-cvp)*ccb2(1)+cy3(1)*csum2+cca2(1)-cp1
           ck3(2)=(csp-cvp)*ccb2(2)+cy3(2)*csum2+cca2(2)-cp2
           ck3(3)=(csp-cvp)*ccb2(3)+cy3(3)*csum2+cca2(3)-cp3

           cy4(0)=ccy1(0)+ck3(0)
           cy4(1)=ccy1(1)+ck3(1)
           cy4(2)=ccy1(2)+ck3(2)
           cy4(3)=ccy1(3)+ck3(3)

	   csp=ccb3(0)*cy4(0)
	   cvp=ccb3(1)*cy4(1)+ccb3(2)*cy4(2)+ccb3(3)*cy4(3)
	   csum1=cce3+2.0*cvp
	   csum2=csum1+2.0*csp
	   csp=cy4(0)*cy4(0)
	   cvp=cy4(1)*cy4(1)+cy4(2)*cy4(2)+cy4(3)*cy4(3)

           cp0=ccx3(1)*cy4(1)+ccx3(2)*cy4(2)+ccx3(3)*cy4(3)
           cp1=ccx3(1)*cy4(0)
           cp2=ccx3(2)*cy4(0)
           cp3=ccx3(3)*cy4(0)

	   ck4(0)=(csp+cvp)*ccb3(0)+cy4(0)*csum1+cca3(0)-cp0
           ck4(1)=(csp-cvp)*ccb3(1)+cy4(1)*csum2+cca3(1)-cp1
           ck4(2)=(csp-cvp)*ccb3(2)+cy4(2)*csum2+cca3(2)-cp2
           ck4(3)=(csp-cvp)*ccb3(3)+cy4(3)*csum2+cca3(3)-cp3

           ccy3(0)=ccy1(0)
     +            +(ck1(0)+ck2(0)+ck2(0)+ck3(0)+ck3(0)+ck4(0))*six
           ccy3(1)=ccy1(1)
     +            +(ck1(1)+ck2(1)+ck2(1)+ck3(1)+ck3(1)+ck4(1))*six
           ccy3(2)=ccy1(2)
     +            +(ck1(2)+ck2(2)+ck2(2)+ck3(2)+ck3(2)+ck4(2))*six
           ccy3(3)=ccy1(3)
     +            +(ck1(3)+ck2(3)+ck2(3)+ck3(3)+ck3(3)+ck4(3))*six
	enddo

	do ii=0,3
	   cnew(ii)=ccy3(ii)
	enddo
 1000   format(100(x,e14.6))
	return
	end 
c
c---------------------------------------------------------------------
c
      subroutine setenergy1(i,iy,px,py,pz,cen,del1,del2,zen,zex,sem)
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'
      integer i,iy
      real px,py,pz,hr,ha
      complex cen,zen,del1(0:3),del2(0:3),zex(3)
      complex opx,opy,opz,sem(4,-mx:mx)
 
      opx=sem(1,i)
      opy=sem(2,i)
      opz=sem(3,i)
      zen=cen-sem(4,i)

      del1(0)= 0.0 !hr+w*ha
      del2(0)= 0.0 !hr-w*ha

      del1(1)= opx
      del2(1)= conjg(opx)

      del1(2)= opy
      del2(2)= conjg(opy)

      del1(3)= opz    
      del2(3)= conjg(opz)

      zex(1)= 0.0 !ex1(i)
      zex(2)= 0.0 !ex2(i)
      zex(3)= 0.0 !ex3(i)
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine setenergy2(i,iy,px,py,pz,cen,del1,del2,zen,zex,sem)
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'
      integer i,iy
      real px,py,pz,hr,ha
      complex cen,zen,del1(0:3),del2(0:3),zex(3)
      complex opx,opy,opz,sem(4,-mx:mx)
 
      opx=sem(1,i)
      opy=sem(2,i)
      opz=sem(3,i)
      zen=cen-sem(4,i)

c     hr=r0s(i)
c     ha=a0s(i)
      del1(0)= 0.0 !hr-w*ha
      del2(0)= 0.0 !hr+w*ha

      del2(1)= opx
      del1(1)= conjg(opx)

      del2(2)= opy
      del1(2)= conjg(opy)

      del2(3)= opz    
      del1(3)= conjg(opz)

      zex(1)=0.0 !-ex1(i)
      zex(2)=0.0 !-ex2(i)
      zex(3)=0.0 !-ex3(i)
      return
      end
c
c---------------------------------------------------------------------
c
      complex function cdiv(d,n)
c
c---------------------------------------------------------------------
c
      implicit none
      complex  d,n,cn,cd
      real c
      cn=conjg(n)
      cd=d*cn
      c=real(n*cn)
      c=1.0/c
      cdiv=cd*c
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine green(it,ip,iy,ien,tem,gamma10,gamma11,gamma12,gamma13,
     +                                  gamma20,gamma21,gamma22,gamma23)
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'
      complex p0,p1,p2,p3,rh,ch,ga,gb,g,ngg,zintegrand
      complex normp(2,2),normm(2,2),ggt(2,2),gtg(2,2)
      complex tem(12)
      complex gamma10,gamma11,gamma12,gamma13
      complex gamma20,gamma21,gamma22,gamma23
      parameter (rh=cmplx(0.25,0.0),ch=cmplx(0.0,-0.25))
      real xd,px,py,pz,v,tpx,tpy,tpz,vpx,vpy,vpz,tlnt
      integer i,iy,it,ip,ien
c
c-- Local varables
c
      complex sp(2,2),gp(2,2)
      complex Gg(2,2),Gt(2,2),prop(4,4)
c
c-- Get going and set up the direction factors 
c
      pz=kz(ip)
      xd=sqrt(1.0-pz*pz)
      px=kx(it)*xd
      py=ky(it)*xd
      v=awei(it)*pwei(ip)
      tpx=0.5*t*v*px  ! 0.5 for the trace over spin
      tpy=0.5*t*v*py
      tpz=0.5*t*v*pz
      tlnt=3.0*v
      vpx=tlnt*px
      vpy=tlnt*py
      vpz=tlnt*pz
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

      zintegrand=(prop(2,4)-prop(1,3)
     +        -conjg(prop(4,2)-prop(3,1)))*rh

      tem(1)=tem(1)+vpx*zintegrand*Rp(ien)
      tem(2)=tem(2)+vpy*zintegrand*Rp(ien)
      tem(3)=tem(3)+vpz*zintegrand*Rp(ien)

      zintegrand=(prop(2,4)+prop(1,3)
     +        -conjg(prop(4,2)+prop(3,1)))*ch
c        write(*,1100)'OP_yk',it,ip,zintegrand,prop(2,4)+prop(1,3),
c    +                            -conjg(prop(4,2)+prop(3,1)),px,py,pz
c        OP_y_ki

      tem(4)=tem(4)+vpx*zintegrand*Rp(ien)
      tem(5)=tem(5)+vpy*zintegrand*Rp(ien)
      tem(6)=tem(6)+vpz*zintegrand*Rp(ien)

      zintegrand=(prop(1,4)+prop(2,3)
     +        -conjg(prop(4,1)+prop(3,2)))*rh
c        write(*,1100)'OP_zk',it,ip,zintegrand,prop(1,4)+prop(2,3),
c    +                            +conjg(prop(4,1)+prop(3,2)),px,py,pz
c        OP_z_ki

      tem(7)=tem(7)+vpx*zintegrand*Rp(ien)
      tem(8)=tem(8)+vpy*zintegrand*Rp(ien)
      tem(9)=tem(9)+vpz*zintegrand*Rp(ien)
c
c-- The current and the diagonal self energy
c
      xd=real(prop(1,1)+prop(2,2))
      tem(10)=tem(10)+tpx*xd*Rp(ien)
      tem(11)=tem(11)+tpy*xd*Rp(ien)
      tem(12)=tem(12)+tpz*xd*Rp(ien)

c     write(*,1100) 'OPx',it,ip,(tem(i),i=1,3)
c     write(*,1100) 'OPy',it,ip,(tem(i),i=4,6)
c     write(*,1100) 'OPz',it,ip,(tem(i),i=7,9)
 1000 format(30(x,e10.3))
 1100 format(a,2(x,i3),30(x,e10.3))
 666  return
      end
c
c----------------------------------------------------------------------
c
      subroutine bulk_gammas(gammaA,sem,cen,px,py,pz,i,iy,ityp)
c
c----------------------------------------------------------------------
c
      include 'qcv.dat'
      integer i,iy,ityp,k
      real px,py,pz,exf(3),exfl
      complex d2,lam,cdiv,sem(4,-mx:mx)
      complex zen,cen,zex(3)
      complex del1(0:3),del2(0:3)
      complex gammaA(0:3),gg,gg0,ggp,ggm
      
c
      if(ityp.eq.1) then
         call setenergy1(i,iy,px,py,pz,cen,del1,del2,zen,zex,sem)
         lam=del1(1)*del2(1)+del1(2)*del2(2)+del1(3)*del2(3)
         lam=w*sqrt(lam-zen*zen)
         gammaA(1)=cdiv(-del1(1),(zen+lam))
         gammaA(2)=cdiv(-del1(2),(zen+lam))
         gammaA(3)=cdiv(-del1(3),(zen+lam))
         gg=zen-gammaA(1)*del2(1)-gammaA(2)*del2(2)-gammaA(3)*del2(3)
         gg=(gammaA(1)*zex(1)+gammaA(2)*zex(2)+gammaA(3)*zex(3))/gg
         gammaA(0)=gg
      else
         call setenergy2(i,iy,px,py,pz,cen,del1,del2,zen,zex,sem)
         lam=del1(1)*del2(1)+del1(2)*del2(2)+del1(3)*del2(3)
         lam=w*sqrt(lam-zen*zen)
         gammaA(1)=cdiv(-del1(1),(zen+lam))
         gammaA(2)=cdiv(-del1(2),(zen+lam))
         gammaA(3)=cdiv(-del1(3),(zen+lam))
         gg=zen+gammaA(1)*del2(1)+gammaA(2)*del2(2)+gammaA(3)*del2(3)
         gg= (gammaA(1)*zex(1)+gammaA(2)*zex(2)+gammaA(3)*zex(3))/gg
         gammaA(0)=gg
      endif
c     write(*,1000) ityp,i,ien,(gammaA(k),k=0,3,1)
 1000 format(a,2(x,i4),30(x,e12.4))

      return
      end
c
c----------------------------------------------------------------------
c
      subroutine AtimesB(P,A,B)
c
c----------------------------------------------------------------------
c
      implicit none
      complex P(2,2),A(2,2),B(2,2)

      P(1,1)=A(1,1)*B(1,1)+A(1,2)*B(2,1)
      P(1,2)=A(1,1)*B(1,2)+A(1,2)*B(2,2)
      P(2,1)=A(2,1)*B(1,1)+A(2,2)*B(2,1)
      P(2,2)=A(2,1)*B(1,2)+A(2,2)*B(2,2)

      return
      end
c
c----------------------------------------------------------------------
