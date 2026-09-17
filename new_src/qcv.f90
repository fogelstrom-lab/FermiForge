!======================================================================
!
!                                        ! taken from  : 21.10.24
   program qcv                           ! last change : 05.01.25
                                         !
!
!======================================================================
!
      use Global_variables
      use MPI_variables
      use initialisation
      use iterate      
      implicit none
!
!---  Give the input and rig the calculation
!
      logical :: AA, BR, BB, NN
      integer :: ierr, md
      real :: p, epsA,epsL
      integer :: igo, istop, irate
      real :: start, finish, ctime, wtime
!
!--   MPI initialization
!
      epsA = 1.e10
      epsL = 1.e10
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierr )

!     print *, myid,nproc
!
! -- Initialise the calculation 
!
      igo = 0
      istop = 0
      irate = 1
      if(myid.eq.0) call init_calc
      call input_bcast

! -- Get the OPs done

      AA=.false.
      BR=.false.
      BB=.false.
      NN=.true.

      iiter=0
      call MPI_BCAST(iiter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      do while ( max(epsA,epsL) > errtol .and. iiter < 5) 
         if(myid == 0) then
            call cpu_time(start)
            call system_clock(igo, irate)
         end if

         iiter = iiter + 1

         call MPI_BCAST(iiter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call iterateSE(epsA,epsL,p,md,AA,BR,BB,NN) 
         call MPI_BCAST(epsA,1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(epsL,1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)

         if(myid == 0) then
            call cpu_time(finish)
            call system_clock(istop, irate)
            ctime=finish-start
            wtime=real(istop-igo)/real(irate)
 
            print 1212,'solved (NN) : ',iiter, ' time (c,w) [s]',ctime,wtime, &
                                  ' error(avg,max)= (',epsA,epsL,')'
         endif
      enddo

      NN=.false.
      if(ittyp == 0) NN=.true.
      if(ittyp == 1) BR=.true.
      if(ittyp == 2) BB=.true.
      if(ittyp == 3) AA=.true.

      iiter=0
      call MPI_BCAST(iiter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myid ==0 .and. NN) open(9,file='error_log_NN.dat',status='unknown')
      if(myid ==0 .and. BR) open(9,file='error_log_BR.dat',status='unknown')
      if(myid ==0 .and. BB) open(9,file='error_log_BB.dat',status='unknown')
      if(myid ==0 .and. AA) open(9,file='error_log_AA.dat',status='unknown')

      do while( max(epsA,epsL) > errtol .and.  iiter < itmax)

         if(myid == 0) then
            call cpu_time(start)
            call system_clock(igo, irate)
         end if

         iiter = iiter + 1
         call MPI_BCAST(iiter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call iterateSE(epsA,epsL,p,md,AA,BR,BB,NN) 

         call MPI_BCAST(epsA,1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(epsL,1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)

         if(myid == 0) then
            call cpu_time(finish)
            call system_clock(istop, irate)
            ctime=finish-start
            wtime=real(istop-igo)/real(irate)
 
            if(AA)  print 1210,'solved (AA) : ',iiter, ' time (c,w) [s]',ctime,wtime, &
                               ' error(avg,max)= (',epsA,epsL,') p =',p, ' md = ',md
            if(BR)  print 1212,'solved (BR) : ',iiter, ' time (c,w) [s]',ctime,wtime, &
                               ' error(avg,max)= (',epsA,epsL,')'
            if(BB)  print 1211,'solved (BB) : ',iiter, ' time (c,w) [s]',ctime,wtime, &
                               ' error(avg,max)= (',epsA,epsL,') alpha = ', p
            if(NN)  print 1212,'solved (NN) : ',iiter, ' time (c,w) [s]',ctime,wtime, &
                               ' error(avg,max)= (',epsA,epsL,')'
            write(9,*) iiter,epsA,epsL
         end if

      end do
      if(myid ==0) close(9)

 1210 format(a, (1x, i4), a, 2(1x, f7.1), a, 2(1x, e11.4), a, 1x, f5.2, a, 1x, i2)
 1211 format(a, (1x, i4), a, 2(1x, f7.1), a, 2(1x, e11.4), a, 1x, f6.2)
 1212 format(a, (1x, i4), a, 2(1x, f7.1), a, 2(1x, e11.4), a)
!
!---  By this point the calculation should be done
!
      call MPI_FINALIZE(ierr)
!
!----------------------------------------------------------------------
!
   end program qcv
!
!----------------------------------------------------------------------
