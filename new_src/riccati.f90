MODULE riccati
   use Global_variables
contains
!--------------------------------------------------------------------
!
   subroutine makeprops(it,ip,ien,en,sem,tem)
!
!---------------------------------------------------------------------
!
      integer, intent(in) :: it, ip, ien
      integer, parameter :: ndx=1, cstep=20
      complex, dimension(4,-mx:mx), intent(in) :: sem
      ! This array is an accumulator over momentum directions and Ozaki
      ! poles.  The historical intent(out) declaration was inconsistent with
      ! the additions performed below and made the program formally invalid.
      complex, dimension(12), intent(inout) :: tem
      complex, intent(in) :: en

      integer :: nn, i
      real :: px, py, pz, d, sd
      complex :: zen1, zen2, cen
      complex, dimension(3) :: zex1, zex2
      complex, dimension(0:3) :: del11, del12, del21, del22, gold, gnew
      complex :: gamma10,gamma11,gamma12,gamma13
      complex :: gamma20,gamma21,gamma22,gamma23
 
      cen = w * en
      pz=kz(ip)
      sd=sqrt(1.0-pz*pz)
      px=kx(it)*sd
      py=ky(it)*sd

!      print *, ien, it, ip, px,py,pz,en,cen
!
!-----------------  gamma1 ------------------------------------
!
      call bulk_gammas(gold,sem,cen,-mx,1)
!     print 1000, ' 1 ',-mx-1,gold,sem(1,-mx)
      d=10*dx
      nn=int(cstep*ndx)
      call setenergy1(-mx,cen,del11,del21,zen1,zex1,sem)
      call ricc(gold,gnew,del11,del11,del21,del21,zen1,zen1,zex1,zex1,d,nn)

      gold = gnew
!     print 1000, ' 1 ',-mx,gnew,sem(1,-mx)

      d = dx
      nn = int(ndx)*(1+int(abs(en)))
      do i=-mx+1,0,1
         call setenergy1(i,cen,del12,del22,zen2,zex2,sem)
         call ricc(gold,gnew,del11,del12,del21,del22,zen1,zen2,zex1,zex2,d,nn)

         del11=del12
         del21=del22
         zen1=zen2
         zex1=zex2

         gold = gnew
      end do
!     print 1000, ' 1 ',i,gnew,sem(1,i)

!     if( gnew(3) /= gnew(3)) stop
      gamma10=gnew(0)
      gamma11=gnew(1)
      gamma12=gnew(2)
      gamma13=gnew(3)
!
!-----------------  gamma2 ------------------------------------
!
      call bulk_gammas(gold,sem,cen,mx,2)
!     print 1000, ' 2 ',mx+1,gold,sem(1,mx)
      d=10*dx
      nn=int(cstep*ndx)
      call setenergy2(mx,cen,del11,del21,zen1,zex1,sem)
      call ricc(gold,gnew,del11,del11,del21,del21,zen1,zen1,zex1,zex1,d,nn)

      gold = gnew
!     print 1000, ' 2 ',mx,gnew,sem(1,mx)

      d= dx
      nn = int(ndx)*(1+int(abs(en)))

      do i=mx-1,0,-1    
         call setenergy2(i,cen,del12,del22,zen2,zex2,sem)
         call ricc(gold,gnew,del11,del12,del21,del22,zen1,zen2,zex1,zex2,d,nn)

         del11=del12
         del21=del22
         zen1=zen2
         zex1=zex2

         gold = gnew
      end do
!     print 1000, ' 2 ',0,gnew,sem(1,i)

      gamma20=gnew(0)
      gamma21=gnew(1)
      gamma22=gnew(2)
      gamma23=gnew(3)

      call green(it,ip,ien,tem,gamma10,gamma11,gamma12,gamma13,gamma20,gamma21,gamma22,gamma23)
!     print 1050, sem(1,0),sem(2,0),sem(3,0),sem(4,0)
!     print 1100, tem
!     print *, ' did ', ien,ip,it
!1000 format(a,1x,i4,10(1x,e12.4))
!1050 format(8(1x,f6.3))
!1100 format(24(1x,f6.3))
!
!---------------------------------------------------------------------
!
   end subroutine makeprops
!
!---------------------------------------------------------------------
!
   subroutine ricc(cold,cnew,ca1,ca2,cb1,cb2,cep1,cep2,cex1,cex2,d,nsteps)
!
!-----------------------------------------------------------------------
!
!  solving i*y'=-D_2*y^2-2*E*y-D_1 with linearly varying
!  D_1 (ca1...ca2), D_2 (cb1...cb2), and E (cep1...cep2).
!  y evolves from cold to cnew in the intervall d.
!  uses nsteps steps.
!
!-----------------------------------------------------------------------
!
      implicit none
 
      integer, intent(in) :: nsteps
      real, intent(in) ::  d
      complex, dimension(0:3), intent(in) :: cold, ca1, ca2, cb1, cb2 
      complex, dimension(3), intent(in) :: cex1, cex2
      complex, intent(in) :: cep1, cep2

      complex, dimension(0:3), intent(out) :: cnew
!
!       **** Temporary variables ****
!
      complex :: cce1, cce2, cce3, cdele
      complex, dimension(3) ::  ccx1, ccx2 ,ccx3 ,cdelx
      complex, dimension(0:3) ::  cca1, cca2, cca3, cdela
      complex, dimension(0:3) ::  ccb1, ccb2, ccb3, cdelb

      complex, dimension(0:3) :: cy1,ck1,cy2,ck2,ccy1
      complex, dimension(0:3) :: cy3,ck3,cy4,ck4,ccy3

      complex :: cp0, cp1, cp2, cp3, csum1, csum2
      complex :: cvp, csp
      real :: h0, dn2, dn, dh, dh2
      real, parameter :: six=1.0/6.0 
      complex, parameter :: cone =cmplx(0.0,1.0)
        
      integer :: i
!
!-----------------------------------------------------------------------
!
!       **** statements ****
!
      h0=d
      dn2=1.0/(nsteps+nsteps)
      dn=dn2+dn2
      dh=h0*dn
      dh2=dh+dh
!
! -- Diagonal S.E. multiplied by i and 2
!
      cce1 = cone * cep1 * dh2
      cce2 = cone * cep2 * dh2
      cce3 = cce1
      cdele = (cce2-cce1) * dn2
!
! --  Exchange field S.E multiplied by i and 2
!
      ccx1 = cone * cex1 * dh2
      ccx2 = cone * cex2 * dh2
      ccx3 = ccx1
      cdelx = (ccx2-ccx1) * dn2
!
! -- Order parameter fields multiplied by i
!
      cca1 = cone * ca1 * dh
      cca2 = cone * ca2 * dh
      cca3 = cca1
      cdela = (cca2-cca1) * dn2

      ccb1 = cone * cb1 * dh
      ccb2 = cone * cb2 * dh
      ccb3 = ccb1
      cdelb = (ccb2-ccb1) * dn2

      ccy3 = cold

      do i=1,nsteps,1

         cce1=cce3
         cce2=cce1+cdele
         cce3=cce2+cdele

         ccx1=ccx3
         ccx2=ccx1+cdelx
         ccx3=ccx2+cdelx

         cca1=cca3
         cca2=cca1+cdela
         cca3=cca2+cdela

         ccb1=ccb3
         ccb2=ccb1+cdelb
         ccb3=ccb2+cdelb

         ccy1=ccy3
         cy1=ccy1

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

         cy2 = ccy1 + ck1 * 0.5

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

         cy3 = ccy1 + ck2 * 0.5

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

         cy4 =  ccy1 + ck3

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

         ccy3 = ccy1+(ck1+ck2+ck2+ck3+ck3+ck4)*six
      end do

      cnew=ccy3
!
!---------------------------------------------------------------------
!
   end subroutine ricc
!
!---------------------------------------------------------------------
!
   subroutine setenergy1(i,cen,del1,del2,zen,zex,sem)
!
!---------------------------------------------------------------------
!
      integer, intent(in) :: i
      complex, dimension(4,-mx:mx), intent(in) :: sem
      complex, dimension(0:3), intent(out) :: del1, del2 
      complex, dimension(3), intent(out) :: zex 
      complex, intent(in) :: cen
      complex, intent(out) :: zen

      del1(0) = 0.0 
      del1(1) = sem(1,i)
      del1(2) = sem(2,i)
      del1(3) = sem(3,i)

      del2 = conjg(del1)

      zen=cen-sem(4,i)
      zex(1)=0.0 !-ex1(i)
      zex(2)=0.0 !-ex2(i)
      zex(3)=0.0 !-ex3(i)
!
!---------------------------------------------------------------------
!
   end subroutine setenergy1
!
!---------------------------------------------------------------------
!
   subroutine setenergy2(i,cen,del1,del2,zen,zex,sem)
!
!---------------------------------------------------------------------
!
      integer, intent(in) :: i
      complex, dimension(4,-mx:mx), intent(in) :: sem
      complex, dimension(0:3), intent(out) :: del1, del2 
      complex, dimension(3), intent(out) :: zex 
      complex, intent(in) :: cen
      complex, intent(out) :: zen

      del2(0) = 0.0 
      del2(1) = sem(1,i)
      del2(2) = sem(2,i)
      del2(3) = sem(3,i)

      del1 = conjg(del2)

      zen=cen-sem(4,i)
      zex(1)=0.0 !-ex1(i)
      zex(2)=0.0 !-ex2(i)
      zex(3)=0.0 !-ex3(i)
!
!---------------------------------------------------------------------
!
   end subroutine setenergy2
!
!---------------------------------------------------------------------
!
   subroutine green(it,ip,ien,tem,gamma10,gamma11,gamma12,gamma13,gamma20,gamma21,gamma22,gamma23)
!
!---------------------------------------------------------------------
!
      integer, intent(in) :: it, ip, ien
      ! Preserve the value accumulated by previous quadrature points.
      complex, dimension(12), intent(inout) :: tem
      complex, intent(in) :: gamma10, gamma11, gamma12, gamma13
      complex, intent(in) :: gamma20, gamma21, gamma22, gamma23
!
! -- some parameters
!
      real, dimension(2,2), parameter :: tau0 = reshape([1.0,0.0,0.0,1.0], [2,2])
      complex, parameter :: rh=cmplx(0.25,0.0), ch=cmplx(0.0,-0.25)
!
! -- Local varables
!
      integer :: i,j
      real :: xd, px, py, pz, v, tpx, tpy, tpz, vpx, vpy, vpz, tlnt
      complex :: ngg, zintegrand
      complex, dimension(2,2) :: Gg, Gt, normp ,normm ,ggt , gtg, ghelp, ahelp
      complex, dimension(4,4) :: prop
!
! -- Get going and set up the direction factors 
!
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
!     print *, ' green:', vpx,vpy,vpz
!
! -- Make the matrix gammas
!
      Gg(1,1)=-gamma11 + w * gamma12
      Gg(2,1)=-gamma10 +     gamma13
      Gg(1,2)= gamma10 +     gamma13
      Gg(2,2)= gamma11 + w * gamma12

      Gt(1,1)=-gamma21 - w * gamma22
      Gt(2,1)=-gamma20 +     gamma23
      Gt(1,2)= gamma20 +     gamma23
      Gt(2,2)= gamma21 - w * gamma22
!
! -- Make the normalisation matrices N_+, N_-
!
      ggt = matmul(Gg,Gt)
      ngg=(1.0-ggt(1,1))*(1.0-ggt(2,2))-ggt(2,1)*ggt(1,2)
      normp(1,1)=(1.0-ggt(2,2))/ngg
      normp(1,2)=     ggt(1,2) /ngg
      normp(2,1)=     ggt(2,1) /ngg
      normp(2,2)=(1.0-ggt(1,1))/ngg

      gtg = matmul(Gt,Gg)
      ngg=(1.0-gtg(1,1))*(1.0-gtg(2,2))-gtg(2,1)*gtg(1,2)
      normm(1,1)=(1.0-gtg(2,2))/ngg
      normm(1,2)=     gtg(1,2) /ngg
      normm(2,1)=     gtg(2,1) /ngg
      normm(2,2)=(1.0-gtg(1,1))/ngg
!
! -- upper diagonal g
!
      ahelp = 0.0
      ghelp = 0.0
      ahelp = tau0 + ggt
      ghelp = -w * matmul(normp,ahelp)
      do j=1,2
         do i=1,2
            prop(i,j) =  ghelp(i,j)
         end do
      end do
!
! --  f
!
      ghelp = 0.0
      ghelp = -2.0 * w * matmul(normp,Gg)
      do j=1,2
         do i=1,2
            prop(i,j+2) =  ghelp(i,j)
         end do
      end do
!
! --  f-tilde
!
      ghelp = 0.0
      ghelp = 2.0 * w * matmul(normm,Gt)
      do j=1,2
         do i=1,2
            prop(i+2,j) =  ghelp(i,j)
         end do
      end do
!
! -- lower diagonal g
!
      ahelp = 0.0
      ghelp = 0.0
      ahelp = tau0 + gtg
      ghelp =  w * matmul(normm,ahelp)
      do j=1,2
         do i=1,2
            prop(i+2,j+2) =  ghelp(i,j)
         end do
      end do

      zintegrand=(prop(2,4)-prop(1,3)-conjg(prop(4,2)-prop(3,1)))*rh

      tem(1)=tem(1)+vpx*zintegrand*Rp(ien)
      tem(2)=tem(2)+vpy*zintegrand*Rp(ien)
      tem(3)=tem(3)+vpz*zintegrand*Rp(ien)

      zintegrand=(prop(2,4)+prop(1,3)-conjg(prop(4,2)+prop(3,1)))*ch
      tem(4)=tem(4)+vpx*zintegrand*Rp(ien)
      tem(5)=tem(5)+vpy*zintegrand*Rp(ien)
      tem(6)=tem(6)+vpz*zintegrand*Rp(ien)

      zintegrand=(prop(1,4)+prop(2,3)-conjg(prop(4,1)+prop(3,2)))*rh

      tem(7)=tem(7)+vpx*zintegrand*Rp(ien)
      tem(8)=tem(8)+vpy*zintegrand*Rp(ien)
      tem(9)=tem(9)+vpz*zintegrand*Rp(ien)
!
!-- The current and the diagonal self energy
!
      xd=real(prop(1,1)+prop(2,2))
      tem(10)=tem(10)+tpx*xd*Rp(ien)
      tem(11)=tem(11)+tpy*xd*Rp(ien)
      tem(12)=tem(12)+tpz*xd*Rp(ien)
!
!----------------------------------------------------------------------
!
   end subroutine green
!
!----------------------------------------------------------------------
!
   subroutine bulk_gammas(gammaA,sem,cen,i,ityp)
!
!----------------------------------------------------------------------
!
      integer, intent(in) :: i, ityp
      complex, dimension(4,-mx:mx), intent(in) :: sem
      complex, intent(in) :: cen
      complex, dimension(0:3), intent(out) :: gammaA

      complex :: d2, lam, gg, zen
      complex, dimension(3) :: zex
      complex, dimension(0:3) :: del1, del2
      
      if(ityp.eq.1) then
         call setenergy1(i,cen,del1,del2,zen,zex,sem)
         lam=del1(1)*del2(1)+del1(2)*del2(2)+del1(3)*del2(3)
         lam=w*sqrt(lam-zen*zen)
         d2= 1.0/(zen+lam)
         gammaA(1)=-del1(1)*d2
         gammaA(2)=-del1(2)*d2
         gammaA(3)=-del1(3)*d2
         gg=zen-gammaA(1)*del2(1)-gammaA(2)*del2(2)-gammaA(3)*del2(3)
         gg=(gammaA(1)*zex(1)+gammaA(2)*zex(2)+gammaA(3)*zex(3))/gg
         gammaA(0)=gg
      else
         call setenergy2(i,cen,del1,del2,zen,zex,sem)
         lam=del1(1)*del2(1)+del1(2)*del2(2)+del1(3)*del2(3)
         lam=w*sqrt(lam-zen*zen)
         d2= 1.0/(zen+lam)
         gammaA(1)=-del1(1)*d2 
         gammaA(2)=-del1(2)*d2 
         gammaA(3)=-del1(3)*d2 
         gg=zen+gammaA(1)*del2(1)+gammaA(2)*del2(2)+gammaA(3)*del2(3)
         gg= (gammaA(1)*zex(1)+gammaA(2)*zex(2)+gammaA(3)*zex(3))/gg
         gammaA(0)=gg
      endif
!     write(*,1000) ityp,i,ien,(gammaA(k),k=0,3,1)
!1000 format(a,2(x,i4),30(x,e12.4))
!
!----------------------------------------------------------------------
!
   end subroutine bulk_gammas
!
!----------------------------------------------------------------------
!
end MODULE riccati
!
!----------------------------------------------------------------------
