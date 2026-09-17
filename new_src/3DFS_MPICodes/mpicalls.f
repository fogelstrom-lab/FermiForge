c---------------------------------------------------------------------
c
      subroutine input_bcast   
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'
      integer ierr

      call MPI_BCAST(pwave,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dwave,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(  cyl,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST( ienmax, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(maxfrec, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(   tmax, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(  itmax, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( istart, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(  Ncmax, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NcCutof, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(   irep, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(   icyl, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(  inmax,10,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(  vort,   1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(     t,   1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(   aa0,   1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(    Rx,   1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(    dx,   1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(deltat,   1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(  esum,   1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(   eps,   1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(  pmax, 10,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(    kx,400,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(    ky,400,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(  awei,400,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(    kz, 40,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(  pwei, 40,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(    zp,100,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(    Rp,100,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
c
c -- Send out the initial self energies
c     
      call dfunc_bcast_1

      return
      end
c
c---------------------------------------------------------------------
c
      subroutine runciterat_bcast(ndeltat)
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'
      integer ierr
      integer ndeltat 
      call MPI_BCAST(ndeltat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

c     write(*,*) ' was in runciterat_bcast',myid,ndeltat
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine dfunc_bcast_1
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'
      integer ierr

      call MPI_BCAST(dxx,pkgsz,MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dxy,pkgsz,MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dxz,pkgsz,MPICMPLX,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(dyx,pkgsz,MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dyy,pkgsz,MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dyz,pkgsz,MPICMPLX,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(dzx,pkgsz,MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dzy,pkgsz,MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dzz,pkgsz,MPICMPLX,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(vx,pkgsz,MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vy,pkgsz,MPICMPLX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vz,pkgsz,MPICMPLX,0,MPI_COMM_WORLD,ierr)

c     write(*,*) ' was in dfunc_bcast_1',myid
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine citerat_bcast_1(n,mmax)
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'
      integer ierr
      integer n,mmax
      call MPI_BCAST(   n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c     write(*,*) ' was in citerat_bcast_1',myid,mmax,n
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine citerat_bcast_2(k)
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'
      integer ierr
      integer k 
      call MPI_BCAST(k,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c     write(*,*) ' was in citerat_bcast_2',myid,k
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine citerat_bcast_3(eps1)
c
c---------------------------------------------------------------------
c
      include 'qcv.dat'
      integer ierr
      real eps1
      erreps=eps1
      call MPI_BCAST(  eps1,1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(erreps,1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
c     write(*,*) ' was in citerat_bcast_3',myid,eps1
      return
      end
c
c---------------------------------------------------------------------
