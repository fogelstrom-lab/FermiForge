MODULE Iterate
   use Global_variables
   use MPI_variables
   use mpicalls 
   use newses
   use packses
contains
!
!---------------------------------------------------------------------
!
   subroutine iterateSE(tol,p,md,AA,BR,BB,NN)
!
!---------------------------------------------------------------------
!
      integer, intent (inout) :: md
      real, intent (inout) :: tol, p 
      logical,intent (in)  :: AA,BR,BB,NN

      complex :: cpp,cpo,cpm,cop,coo,com,cmp,cmo,cmm

      md = 0
      norm = 0.0
      if(AA) call iterateAA(tol,p,md)        ! Anderson mixing scheme
      if(BR) call iterateBR(tol)             ! Broyden update scheme
      if(BB) call iterateBB(tol,p)           ! Barzilai-Borwein update scheme
      if(NN) call iterateNN(tol)             ! simple update x_new = g(x_old)  

      call SES_bcast  ! Broadcast the new set of self energies

      if(myid == 0) then

         open(1,file='op_xyz',status='unknown')
         open(2,file='op_harm',status='unknown')
         open(3,file='curr',status='unknown')

         do ix=0,nx,1
            x=xgrid(ix)

            cmm=0.5*(dxx(ix)-dyy(ix)+w*(dxy(ix)+dyx(ix)))
            cmo=sqrth*(dxz(ix)+w*dyz(ix))
            cmp=0.5*(dxx(ix)+dyy(ix)-w*(dxy(ix)-dyx(ix)))

            com=sqrth*(dzx(ix)+w*dzy(ix))
            coo=dzz(ix)
            cop=sqrth*(dzx(ix)-w*dzy(ix))

            cpm=0.5*(dxx(ix)+dyy(ix)+w*(dxy(ix)-dyx(ix)))
            cpo=sqrth*(dxz(ix)-w*dyz(ix))
            cpp=0.5*(dxx(ix)-dyy(ix)-w*(dxy(ix)+dyx(ix)))

            write(1,1000) x,dxx(ix),dxy(ix),dxz(ix),dyx(ix),dyy(ix),dyz(ix),dzx(ix),dzy(ix),dzz(ix)
            write(2,1000) x,cpp,cpo,cpm,cop,coo,com,cmp,cmo,cmm
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

 1000 format((1x,f8.3),30(1x,e14.6))
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
   include 'iter_BR.f90' ! This the Broyden update scheme
!                        ! returns the update error in variable tol
!-----------------------------------------------------------
!
   include 'iter_BB.f90' ! This the Barzilai-Borwein update scheme
!                        ! returns the update error in variable tol
!                        ! and the accelaration factor in p
!---------------------------------------------------------------------
!
   include 'iter_AA.f90' ! This the Anderson-mixing update scheme
!                        ! returns the update error in variable tol,
!                        ! the accelaration factor in p,
!                        ! and the current numbers of vectors considered for mixing
!---------------------------------------------------------------------
!
end module Iterate
!
!-----------------------------------------------------------------------
