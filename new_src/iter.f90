MODULE Iterate
   use Global_variables
   use MPI_variables
   use mpicalls 
   use newses
   use packses
   use anderson_iteration_adapter, only : iterateAA
   implicit none
contains
!
!---------------------------------------------------------------------
!
   subroutine iterateSE(tolA,tolL,p,md,AA,BR,BB,NN)
!
!---------------------------------------------------------------------
!
      integer, intent (inout) :: md
      real, intent (inout) :: tolA,tolL, p 
      logical,intent (in)  :: AA,BR,BB,NN
      integer :: ix
      real :: d, x

      md = 0
      if(AA) call iterateAA(tolA,tolL,p,md)        ! Anderson mixing scheme
      if(BR) call iterateBR(tolA,tolL)             ! Broyden update scheme
      if(BB) call iterateBB(tolA,tolL,p)           ! Barzilai-Borwein update scheme
      if(NN) call iterateNN(tolA,tolL)             ! simple update x_new = g(x_old)  

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

      call SES_bcast  ! Broadcast the new set of self energies

      if(myid == 0) then

         open(1,file='op_xyz',status='unknown')
         open(2,file='op_harm',status='unknown')
         open(3,file='curr',status='unknown')

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
      end if

 1000 format((1x,f8.3),30(1x,e15.9))
!
!---------------------------------------------------------------------
!
   end subroutine iterateSE
!
!---------------------------------------------------------------------
!
   include 'iter_NN.f90' ! This the simplest update scheme x_new = G(x_old)
!                        ! returns the update error in variable tol
!-----------------------------------------------------------
!
   include 'Iter_BR.f90' ! This the Broyden update scheme
!                        ! returns the update error in variable tol
!-----------------------------------------------------------
!
   include 'iter_BB.f90' ! This the Barzilai-Borwein update scheme
!                        ! returns the update error in variable tol
!                        ! and the accelaration factor in p
!---------------------------------------------------------------------
!
end module Iterate
!
!-----------------------------------------------------------------------
