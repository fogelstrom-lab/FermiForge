MODULE Global_variables
   implicit none 

   integer, parameter :: nx = 99, mx = 800, mfrec = 400, vecl = nx+1, pkgsz = nx+1
   real, parameter :: pi = acos(-1.0), sqrt2 = sqrt(2.0), sqrth = 1.0/sqrt2
   complex, parameter :: w = cmplx(0.0,1.0)

   integer :: icyl, tmax, istart, itmax, ittyp
   real :: t, vort, aa0, errtol, aa_pmax

   logical :: pwave,dwave,cyl,outerr
   integer :: maxfrec, iiter
   real :: deltat, vc1, vc2, vc3, esum, Rx, dx, tx
   real, dimension(0:nx) :: xgrid
   real, dimension(0:400) :: kx, ky, awei
   real, dimension(0:40) :: kz, pwei
 
   integer, parameter :: nop=21  ! Number of selfenergies (real values)
   complex, dimension(0:nx) :: dxx, dxy, dxz, dyx, dyy, dyz, dzx, dzy, dzz, vx, vy, vz
   complex, dimension(0:nx) :: cpp, cpo, cpm, cop, coo, com, cmp, cmo, cmm

   integer ::  Ncmax, NcCutof
   real, dimension(100) :: zp, Rp
!
!------------------------------------------------------------  
!
end MODULE Global_variables
!
!------------------------------------------------------------  
!
MODULE MPI_variables
   use mpi
   implicit none
!
!------------------------------------------------------------
!
   integer :: myid,nproc

   integer, parameter  :: MPIFLOAT=MPI_DOUBLE_PRECISION
   integer, parameter  :: MPICMPLX=MPI_DOUBLE_COMPLEX
!
!-------------------------------------------------------------
!
end MODULE MPI_variables
!
!-------------------------------------------------------------
