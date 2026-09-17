MODULE NewSES
   use global_variables
   use MPI_variables
   use MPI_variables
   use Interpolations
   use riccati
contains
!
!----------------------------------------------------------------------
!
   subroutine getnewop
!
!---------------------------------------------------------------------
!
      integer :: iy,ii,it,ip,ien,ierr
      real :: tlnt
      complex, parameter :: czero = cmplx(0.0,0.0)
      complex :: en
      complex, dimension(0:nx) :: dlxx, dlxy, dlxz
      complex, dimension(0:nx) :: dlyx, dlyy, dlyz
      complex, dimension(0:nx) :: dlzx, dlzy, dlzz, lvx, lvy, lvz
      complex, dimension(4,-mx:mx) :: sem
      complex, dimension(12) :: tem
 
      tlnt=esum !t/log(t)

      dlxx = czero
      dlxy = czero
      dlxz = czero
      dlyx = czero
      dlyy = czero
      dlyz = czero
      dlzx = czero
      dlzy = czero
      dlzz = czero
      lvx = czero
      lvy = czero
      lvz = czero

!     iy = 0
!     ip = 6
!     it = 12
      do iy = myid, nx , nproc

         tem = czero

         do it = 1, tmax, 1
            do ip = 1, 11, 1
               call intord_v(it,ip,iy,sem)
               do ien= 1, Ncmax, 1
                  en = zp(ien)
                  call makeprops(it,ip,ien,en,sem,tem)
               end do
            end do
         end do

         do ii=1,9
            tem(ii)=tem(ii)*tlnt
         end do
         dlxx(iy)=tem(1)
         dlxy(iy)=tem(2)
         dlxz(iy)=tem(3)

         dlyx(iy)=tem(4)
         dlyy(iy)=tem(5)
         dlyz(iy)=tem(6)

         dlzx(iy)=tem(7)
         dlzy(iy)=tem(8)
         dlzz(iy)=tem(9)

         lvx(iy)=tem(10)
         lvy(iy)=tem(11)
         lvz(iy)=tem(12)
      end do
! 1000 format(3(1x,i4),18(1x,e10.3))
!
! -- send out the new op suggestion
!         
      dxx = czero
      dxy = czero
      dxz = czero
      dyx = czero
      dyy = czero
      dyz = czero
      dzx = czero
      dzy = czero
      dzz = czero
      vx = czero
      vy = czero
      vz = czero

      call MPI_REDUCE(dlxx,dxx,pkgsz,MPICMPLX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(dlxy,dxy,pkgsz,MPICMPLX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(dlxz,dxz,pkgsz,MPICMPLX,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      call MPI_REDUCE(dlyx,dyx,pkgsz,MPICMPLX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(dlyy,dyy,pkgsz,MPICMPLX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(dlyz,dyz,pkgsz,MPICMPLX,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      call MPI_REDUCE(dlzx,dzx,pkgsz,MPICMPLX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(dlzy,dzy,pkgsz,MPICMPLX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(dlzz,dzz,pkgsz,MPICMPLX,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      call MPI_REDUCE(lvx,vx,pkgsz,MPICMPLX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(lvy,vy,pkgsz,MPICMPLX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(lvz,vz,pkgsz,MPICMPLX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!
!---------------------------------------------------------------------
!
   end subroutine getnewop
!
!---------------------------------------------------------------------
!
end MODULE newSES
!
!---------------------------------------------------------------------
