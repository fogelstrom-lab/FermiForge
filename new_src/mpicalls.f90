MODULE mpicalls
   use global_variables
   use mpi_variables
   implicit none
contains
!---------------------------------------------------------------------
!
   subroutine input_bcast   
!
!---------------------------------------------------------------------
!
      integer :: ierr
      
      call MPI_BCAST(pwave,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dwave,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(  cyl,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(   tmax, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(  ittyp, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(  itmax, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( istart, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(  Ncmax, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NcCutof, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(   icyl, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(  vort,   1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(     t,   1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(   aa0,   1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xgrid(0),size(xgrid),MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(    Rx,   1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(    tx,   1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(    dx,   1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(deltat,   1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(  esum,   1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(errtol,   1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(aa_pmax,  1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( kx(1),tmax,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( ky(1),tmax,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(awei(1),tmax,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( kz(1),   11,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(pwei(1),  11,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( zp(1),Ncmax,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( Rp(1),Ncmax,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
!
! -- Send out the initial self energies
!     
      call SES_bcast
!
!---------------------------------------------------------------------
!
   end subroutine input_bcast
!
!---------------------------------------------------------------------
!
   subroutine SES_bcast
!
!---------------------------------------------------------------------
!
      integer :: ierr

!     print *, 'sending out SES', myid
      call MPI_BCAST(dxx(0),size(dxx),MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dxy(0),size(dxy),MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dxz(0),size(dxz),MPICMPLX,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(dyx(0),size(dyx),MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dyy(0),size(dyy),MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dyz(0),size(dyz),MPICMPLX,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(dzx(0),size(dzx),MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dzy(0),size(dzy),MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dzz(0),size(dzz),MPICMPLX,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(cpp(0),size(cpp),MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cpo(0),size(cpo),MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cpm(0),size(cpm),MPICMPLX,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(cop(0),size(cop),MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(coo(0),size(coo),MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(com(0),size(com),MPICMPLX,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(cmp(0),size(cmp),MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cmo(0),size(cmo),MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cmm(0),size(cmm),MPICMPLX,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(vx(0),size(vx),MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vy(0),size(vy),MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vz(0),size(vz),MPICMPLX,0,MPI_COMM_WORLD,ierr)
!
!---------------------------------------------------------------------
!
   end subroutine SES_bcast
!
!---------------------------------------------------------------------
!
end MODULE mpicalls
!
!---------------------------------------------------------------------
