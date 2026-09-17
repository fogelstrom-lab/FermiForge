!-----------------------------------------------------------------------------!
i Coded by Tomas Löfwander, first verion, Evanston, July 2002                 !
!                     new implementation, Evanston, September 2003            !
!                     packaged into a module, Karlsruhe, January 2005         !
!                     updated, December 2005                                  !
!-----------------------------------------------------------------------------!
! Note: compile with flag -r8 (or equivalent) to get variables                !
!       with double precision.                                                !
!-----------------------------------------------------------------------------!
module scsolver

  implicit none
  
  ! modes:
  integer, parameter :: START=1
  integer, parameter :: STOPP=2
  integer, parameter :: WORK=3

  ! methods:
  integer, parameter :: RELAX=4 
  integer, parameter :: MOBESC=5

  ! itertypes:
  character(len=9), parameter :: IT_START =' IT_Start'
  character(len=9), parameter :: IT_STOP  ='  IT_Stop'
  character(len=9), parameter :: IT_BUILD =' IT_Build'
  character(len=9), parameter :: IT_NORMAL='IT_Normal'
  character(len=9), parameter :: IT_BAD   ='   IT_Bad'
  character(len=9), parameter :: IT_GOOD  ='  IT_Good'
  character(len=9), parameter :: IT_ERROR =' IT_Error'

contains

!-----------------------------------------------------------------------------!
subroutine scs(mode,method,m_max,itertyp,vopt,pm,pm_max,fnorm,epsG,epsL,nf,ng,fnewC)
!-----------------------------------------------------------------------------!

  implicit none

  integer, intent(in) :: mode                 ! mode flag
  integer, intent(in) :: method               ! method flag
  integer, intent(in) :: m_max                ! max number deviation vectors kept in memory
  character(len=9), intent(out) :: itertyp    ! result of throw-away procedure
  logical, intent(in) :: vopt                 ! verbose flag
  real, intent(inout) :: pm                   ! back-coupling strength
  real, intent(in)    :: pm_max               ! max allowed back-coupling strength
  real, intent(out) :: fnorm,epsG,epsL        ! accuracy measures (see above)
  integer, dimension(1), intent(out) :: nf,ng ! positions of biggest errors
  complex, dimension(:), pointer :: fnewC     ! double use:
  ! fnewC in mode=START: first guess of self-energies given as input
  ! fnewC in mode=WORK:  new deviation vector at input
  !                      new guess of self-energies at output
  ! fnewC in mode=STOP:  not used

!-----------------------------------------------------------------------------!
! General form of the main program that use the scs routine:                  !
! (the subroutine run_scsolver below)                                         !
!                                                                             !
! 1. guess SE (selfenergies)                                                  !
! 2. init scs by calling with mode=START and fnewC=SE                         !
! 3. do while( {fnorm,epsG,epsL} > required accuracy )                        !
!      compute fnewC (LHS-RHS of SE-equations)                                !
!      get new guess of SE by calling scs with mode=WORK, include fnewC       !
!      the new guess of SE is now in fnewC                                    !
!    end do                                                                   !
! 4. finalize scs by calling with mode=STOP                                   !
! 5. fnewC now contains the self-consistent SE                                ! 
!                                                                             !
!-----------------------------------------------------------------------------!
! Error estimates:                                                            !
!                         || dg ||                    max |dg |               !
!         global: epsG = ---------- ,  local: epsL = ----------,              !
!                         ||  g ||                    max | g |               !
!                                                                             !
!                   fnorm=||f||                                               !
!                                                                             !
! f_i(x) = g_i(x) - G[g_i(x)] (self-consistency equation at current iteration)!
! dg(x) = g_i(x)-g_{i-1}(x) (change of self-energy at current iteration)      !
!-----------------------------------------------------------------------------!
! scs use the LAPACK routine dsysv (http://www.netlib.org)                    !
! We include all necessary files for this to work in scsLA.f                  !
!-----------------------------------------------------------------------------!
  interface
     subroutine dsysv(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, DWORK, LWORK, INFO)
       CHARACTER UPLO
       INTEGER INFO, LDA, LDB, LWORK, N, NRHS
       INTEGER IPIV( * )
       real*8 A( LDA, * ), B( LDB, * ), DWORK( * )
     end subroutine dsysv
  end interface
!-----------------------------------------------------------------------------!

  logical :: swapped
  
  integer :: k,l,info,lwork
  integer, save :: m,n,itn
!  integer, dimension(1) :: nf,ng
  integer, dimension(m_max-1) :: IPIV

  real*8, dimension(:,:), allocatable :: A
  real*8, dimension(:), allocatable :: b,dwork

  real :: rndnr,pmt,rm,qm,fn0
  real, save :: fnorm_first

  real, dimension(2*size(fnewC)) :: achs,dfk,dfl,gs0,fs0
  real, dimension(2*size(fnewC)) :: fnew,fmin_old,gmin

  real, dimension(:), allocatable, save :: gnew,fmin,fnorms

  real, dimension(:,:), allocatable, save :: fs,gs

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
  if(mode.eq.START)then

     if(vopt) write(*,*) 'scs message: START - allocating memory'

     n=2*size(fnewC)
     
     if(allocated(gnew  )) deallocate(gnew  ); allocate(gnew(n))
     if(allocated(fmin  )) deallocate(fmin  ); allocate(fmin(n))
     if(allocated(fs    )) deallocate(fs    ); allocate(fs(n,m_max))
     if(allocated(gs    )) deallocate(gs    ); allocate(gs(n,m_max))
     if(allocated(fnorms)) deallocate(fnorms); allocate(fnorms(m_max))
     if(vopt) write(*,*) 'scs message: START - allocating memory --done'

     m=0 ! number of deviation vectors in memory
     itn=0 ! internal iteration number
     gnew(1:n/2)=real(fnewC(1:n/2)); gnew(n/2+1:n)=aimag(fnewC(1:n/2))

     itertyp=IT_START

     fmin = 0.0
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
  elseif(mode.eq.STOPP)then

     if(vopt) write(*,*) 'scs message: STOPP - deallocating memory'

     deallocate(gnew)
     deallocate(fmin)
     deallocate(fs)
     deallocate(gs)
     deallocate(fnorms)

     if(vopt) write(*,*) 'scs message: STOPP - deallocating memory --done'

     itertyp=IT_STOP

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
  elseif(mode.eq.WORK)then

     if(vopt) write(*,*) 'scs message: WORK'

     ! insert the new deviation vector:

     itn = itn+1

     fnew(1:n/2)=real(fnewC(1:n/2)); fnew(n/2+1:n)=aimag(fnewC(1:n/2))

     if(itn.eq.1)then
        fnorm_first=sqrt(dot_product(fnew,fnew))
        fnorm=1.0
     else
        fnorm=sqrt(dot_product(fnew,fnew))/fnorm_first
     end if

!-----------------------------------------------------------------------------!
     ! throw away procedure:
     
     if(m.eq.m_max)then

!        if(vopt) write(*,*) 'unsorted:','itn=',itn,'fnorms=',fnorms(1:m)

        ! sort the norms in descending order by bubble sort
        do
           swapped = .false.
           do k = m, 2, -1
              if( fnorms(k).gt.fnorms(k-1) )then ! swap

                 fn0 = fnorms(k)
                 gs0 = gs(1:n,k)
                 fs0 = fs(1:n,k)

                 fnorms(k) = fnorms(k-1)
                 gs(1:n,k) = gs(1:n,k-1)
                 fs(1:n,k) = fs(1:n,k-1)

                 fnorms(k-1) = fn0
                 gs(1:n,k-1) = gs0
                 fs(1:n,k-1) = fs0

                 swapped=.true.
              end if
           end do
           if(.not.swapped) exit
        end do

!!$        do k=2,m
!!$           fn0=fnorms(k)
!!$           gs0=gs(1:n,k)
!!$           fs0=fs(1:n,k)
!!$           l=k-1
!!$!write(*,*) 'k,l = ',k,l
!!$!           do while(l.gt.0.and.fnorms(l).lt.fn0)
!!$           do
!!$              fnorms(l+1)=fnorms(l)
!!$              gs(1:n,l+1)=gs(1:n,l)
!!$              fs(1:n,l+1)=fs(1:n,l)
!!$              if(l.eq.1.or.fnorms(l).lt.fn0) exit
!!$              l=l-1
!!$           end do
!!$           fnorms(l+1)=fn0
!!$           gs(1:n,l+1)=gs0
!!$           fs(1:n,l+1)=fs0
!!$        end do

!        if(vopt) write(*,*) '  sorted:','itn=',itn,'fnorms=',fnorms(1:m)

        if(fnorm.gt.fnorms(1).or.fnorm.lt.1.0e-2*fnorms(m))then

           ! a very bad step .or. a really good one
           ! throw away half of the vectors
           if(fnorm.gt.fnorms(1))then
              itertyp=IT_BAD
              if(vopt) write(*,*) 'scs message: a very bad step, '&
                   //'we throw away half of the vectors'
           else
              itertyp=IT_GOOD
              if(vopt) write(*,*) 'scs message: a very good step, '&
                   //'we throw away half of the vectors'
           end if
           do k=1,m/2
              gs(1:n,k)=gs(1:n,k+m/2)
              fs(1:n,k)=fs(1:n,k+m/2)
              fnorms(k)=fnorms(k+m/2)
           end do
           m=m/2+1

        else

           ! normal step, throw away the largest vector only
           itertyp=IT_NORMAL
           if(vopt) write(*,*) 'scs message: a normal step, we throw'&
                //'away the largest vector only'
           do k=1,m-1
              gs(1:n,k)=gs(1:n,k+1)
              fs(1:n,k)=fs(1:N,k+1)
              fnorms(k)=fnorms(k+1)
           end do

        end if

!!$        if(vopt)then
!!$           do k=1,m-1
!!$              write(*,*) 'k=',k,'norm=',fnorms(k),sqrt(dot_product(fs(1:n,k),fs(1:n,k)))/fnorm_first
!!$           end do
!!$           write(*,*) 'k=',m,'norm=',fnorm
!!$        end if

     else

        ! no throwing away, only saving in memory

        itertyp=IT_BUILD

        m=m+1

     end if

     gs(1:n,m)=gnew
     fs(1:n,m)=fnew
     fnorms(m)=fnorm

!-----------------------------------------------------------------------------!
     if(method.eq.RELAX)then

        if(vopt) write(*,*) 'scs message: RELAX'

        gnew=gs(1:n,m)+pm*fs(1:n,m)

        if(vopt) write(*,*) 'scs message: RELAX --done'

!-----------------------------------------------------------------------------!
     elseif(method.eq.MOBESC.and.m.gt.1)then

        if(vopt) write(*,*) 'scs message: MOBESC'

        if(m.eq.2)then

           dfk=fs(1:n,1)-fs(1:n,2)
           rm=-dot_product(dfk,fs(1:n,2))/dot_product(dfk,dfk)
           fmin=fs(1:n,2)+rm*(fs(1:n,1)-fs(1:n,2))
           gmin=gs(1:n,2)+rm*(gs(1:n,1)-gs(1:n,2))
           dfk=gmin-gs(1:n,2)
           pm=sqrt(dot_product(dfk,dfk)/dot_product(fs(1:n,2),fs(1:n,2)))
!!$           if(vopt) write(*,*) 'scs message: ','m=',m,'w1=',rm,'pm=',pm

        else

           info=0
           lwork=64*(m_max-1)
           allocate(dwork(lwork))
           allocate(A(m-1,m-1))
           allocate(b(m-1))

           do k=1,m-1  ! calculate the A matrix and b vector assuming that
              ! f_m is the origin on the hypersurface
              dfk=fs(1:n,k)-fs(1:n,m)
              b(k)=-dot_product(dfk,fs(1:n,m))
              do l=1,m-1
                 dfl=fs(1:n,l)-fs(1:n,m)
                 A(k,l)=dot_product(dfk,dfl)
              end do
           end do

           ! invert the symmetric matrix A:
           call dsysv('U',m-1,1,A,m-1,IPIV,b,m-1,dwork,lwork,info)

           if(info.eq.0)then  ! success

              if(vopt) write(*,*) 'scs message: m=',m,'weights=',b
              ! compute f_min and g_min:
              fmin_old=fmin
              fmin=fs(1:n,m)
              gmin=gs(1:n,m)
              do k=1,m-1
                 fmin=fmin+b(k)*(fs(1:n,k)-fs(1:n,m))
                 gmin=gmin+b(k)*(gs(1:n,k)-gs(1:n,m))
              end do

              ! update the back-coupling pm:
              dfk=gmin-gs(1:n,m)
              pmt=sqrt(dot_product(dfk,dfk)/dot_product(fs(1:n,m),fs(1:n,m)))
              rm=sqrt(dot_product(fmin_old,fmin_old)/&
                   &dot_product(fs(1:n,m),fs(1:n,m)))
              qm=min(1.0,rm)              
              if(qm>0.0) pm=pmt*qm

              if(pm.gt.pm_max) pm=pm_max
              
              if(vopt) write(*,*) 'scs message: ','m=',m,'pmt=',pmt,&
                   'rm=',rm,'qm=',qm,'pm=',pm

           elseif(info.gt.0)then  ! singular A matrix

              if(vopt) write(*,*) 'scs warning: stochastic step at m=',m,&
                   'info=',info
              gmin=gs(1:n,m)
              fmin=(0.0,0.0)
              do k=1,m
                 call random_number(rndnr)
                 fmin=fmin+rndnr*fs(1:n,m)
              end do
              fmin=fmin/(real(m,8))
              ! currently we use old pm, but pm could be tuned here

           else

              write(*,*) 'scs error: illegal argument info=',info
              itertyp=IT_ERROR

           end if

           deallocate(dwork)
           deallocate(A)
           deallocate(b)

        end if

        gnew=gmin+pm*fmin
!gnew=gmin-pm*fmin
        if(vopt) write(*,*) 'scs message: MOBESC --done'

!-----------------------------------------------------------------------------!
     elseif(method.eq.MOBESC.and.m.lt.2)then

        write(*,*) 'scs error: can not use MOBESC method, not enough vectors'

        itertyp=IT_ERROR

     else ! just in case

        write(*,*) 'scs error: unknown method - use RELAX or MOBESC'

        itertyp=IT_ERROR

!-----------------------------------------------------------------------------!
     end if ! method

     ! compute the error estimates:
     achs=gnew-gs(1:n,m)
     epsG=sqrt(dot_product(achs,achs))/sqrt(dot_product(gnew,gnew))
     epsL=maxval(abs(achs))/maxval(abs(gnew))
     nf=maxloc(abs(achs))
     ng=maxloc(abs(gnew))

     ! send back the new guess of self-energies
     fnewC=cmplx(gnew(1:n/2),gnew(n/2+1:n),8)

     if(vopt) write(*,*) 'scs message: WORK --done'

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
  else ! just in case...

     write(*,*) 'scs error: unknown mode - use START, WORK, or STOPP'

     itertyp=IT_ERROR

  end if ! mode

!-----------------------------------------------------------------------------!
end subroutine scs

!-----------------------------------------------------------------------------!
subroutine SCSflush(nout,thefname,FGtyp,xp,yp)
!-----------------------------------------------------------------------------!

  ! flush data with 17 decimals
  implicit none

  integer, intent(in) :: nout
  character(len=255), intent(in) :: thefname
  character(len=1), intent(in) :: FGtyp
  real, dimension(:), pointer :: xp
  complex, dimension(:), pointer :: yp

  integer k,n,KK,NN

  KK=size(xp)
  NN=size(yp)

  if(NN<KK)then
     write(*,*) 'SCSflush error: vectors with wrong sizes'
  else
     if(FGtyp.eq.'F')then
        open(nout,file=trim(thefname)//'_scsF.dat')
     elseif(FGtyp.eq.'G')then
        open(nout,file=trim(thefname)//'_scsG.dat')
     else
        open(nout,file=trim(thefname)//'.dat')
     end if
     do n=1,NN/KK
        do k=1,KK
           write(nout,"(e26.17e3,e26.17e3,e26.17e3)") &
                &(n-1)*(xp(KK)-xp(1))+xp(k),&
                &yp((n-1)*KK+k)
        end do
     end do
     close(nout)
  end if

!-----------------------------------------------------------------------------!
end subroutine SCSflush

!-----------------------------------------------------------------------------!
subroutine LEflush(nout,thefname,i,itertyp,pm,fnorm,eps1,eps2,nf,ng)
!-----------------------------------------------------------------------------!

  ! flush the last results of the scs-loop
  implicit none

  integer, intent(in) :: nout
  character(len=255), intent(in) :: thefname
  integer, intent(in) :: i
  character(len=9), intent(in) :: itertyp
  real, intent(in) :: pm,fnorm,eps1,eps2
  integer, dimension(1), intent(in) :: nf,ng

  if(nout.ne.6) open(nout,file=trim(thefname)//'_scsLE.dat')
  write(nout,"(A6,A10,A15,A15,A15,A15,A9,A9)") &
       &'i',&
       &'itertyp',&
       &'pm',&
       &'fnorm',&
       &'eps1',&
       &'eps2',&
       &'nf',&
       &'ng'
  write(nout,"(I6,A10,e15.6e3,e15.6e3,e15.6e3,e15.6e3,I9,I9)") &
       &i,&
       &itertyp,&
       &pm,&
       &fnorm,&
       &eps1,&
       &eps2,&
       &nf,&
       &ng
  if(nout.ne.6) close(nout)

!-----------------------------------------------------------------------------!
end subroutine LEflush

!-----------------------------------------------------------------------------!
subroutine Eflush(task,nout,thefname,i,itertyp,pm,fnorm,eps1,eps2,nf,ng)
!-----------------------------------------------------------------------------!

  ! flush the last results of the scs-loop
  implicit none

  character(len=1), intent(in) :: task
  integer, intent(in) :: nout
  character(len=255), intent(in) :: thefname
  integer, intent(in) :: i
  character(len=9), intent(in) :: itertyp
  real, intent(in) :: pm,fnorm,eps1,eps2
  integer, dimension(1), intent(in) :: nf,ng

  if(task.eq.'o')then
     if(nout.ne.6) open(nout,file=trim(thefname)//'_scsE.dat')
     write(nout,"(A6,A10,A15,A15,A15,A15,A9,A9)") &
          &'i',&
          &'itertyp',&
          &'pm',&
          &'fnorm',&
          &'eps1',&
          &'eps2',&
          &'nf',&
          &'ng'
     call flush(nout)
  elseif(task.eq.'c')then
     if(nout.ne.6) close(nout)
  elseif(task.eq.'w')then
     write(nout,"(I6,A10,e15.6e3,e15.6e3,e15.6e3,e15.6e3,I9,I9)") &
          &i,&
          &itertyp,&
          &pm,&
          &fnorm,&
          &eps1,&
          &eps2,&
          &nf,&
          &ng
     call flush(nout)
  else
     write(*,*) 'Eflush error: no task'
  end if

  call flush(nout)
!-----------------------------------------------------------------------------!
end subroutine Eflush

!-----------------------------------------------------------------------------!
subroutine E2flush(task,nout,thefname,i,fnorm,eps1,eps2)
!-----------------------------------------------------------------------------!
  character(len=1), intent(in) :: task
  integer, intent(in) :: nout
  character(len=255), intent(in) :: thefname
  integer, intent(in) :: i
  real, intent(in) :: fnorm,eps1,eps2

  if( task.eq.'o' )then
     if(nout.ne.6) open(nout,file=trim(thefname)//'_scsE2.dat')
  elseif( task.eq.'c' )then
     if(nout.ne.6) close(nout)
  elseif( task.eq.'w' )then
     write(nout,"(i5,3e15.6e3)") i,fnorm,eps1,eps2
  else
     write(*,*) 'E2flush error: no task'
  end if
  
  call flush(nout)
!-----------------------------------------------------------------------------!
end subroutine E2flush

!-----------------------------------------------------------------------------!
subroutine finalflush(nout,thefname,i,itertyp,pm,fnorm,eps1,eps2,nf,ng)
!-----------------------------------------------------------------------------!

  ! flush the final results of the scs-loop
  implicit none

  integer, intent(in) :: nout
  character(len=255), intent(in) :: thefname
  integer, intent(in) :: i
  character(len=9), intent(in) :: itertyp
  real, intent(in) :: pm,fnorm,eps1,eps2
  integer, dimension(1), intent(in) :: nf,ng

  write(nout,"(A20,A15,A15,A15,A9,A9)") &
       &'number of iterations',&
       &'fnorm',&
       &'eps1',&
       &'eps2',&
       &'nf',&
       &'ng'
  write(nout,"(I20,e15.6e3,e15.6e3,e15.6e3,I9,I9)") &
       &i,&
       &fnorm,&
       &eps1,&
       &eps2,&
       &nf,&
       &ng

!-----------------------------------------------------------------------------!
end subroutine finalflush

!-----------------------------------------------------------------------------!
subroutine run_scsolver(LHSmRHS,xp,SESp,&
     &itermax,m_max,nrr,pm,pm_max,epsL,vEopt,vopt,sopt,noutE,filnamn,success)
!subroutine run_scsolver(LHSmRHS,logfil,filnamn,mess,xp,SESp,&
!     &itermax,m_max,nrr,pm,fnorm,epsG,epsL,vopt,sopt,noutE,noutE2,noutG,noutF,success)
!-----------------------------------------------------------------------------!

  implicit none

  integer :: logfil ! logfile number where messages are flushed
  integer :: itermax ! max number of iterations
  integer :: m_max ! number of deviation vectors kept in memory
  integer :: nrr ! number of initial simple relaxation iterations to be done
  integer :: noutE,noutE2,noutG,noutF ! file numbers for flushing errors and vectors
  real :: pm ! initial back-coupling strength
  real :: pm_max ! max allowed back-coupling strength
  real :: fnorm,epsG,epsL ! required norm and errors
  logical :: vEopt,vopt,sopt ! verbose and save flags
  character(len=255) :: filnamn,mess ! basic filename, message used in logfile
  real, dimension(:), pointer :: xp ! "grid"
  complex, dimension(:), pointer :: SESp ! initial guess of self-energy
                                         ! --> self-consistent self-energy
  logical, intent(out) :: success ! flag
  
  interface
     subroutine LHSmRHS(SESp,fresh)
       implicit none
       complex, dimension(:), pointer, optional :: SESp ! selfenergies @ input, self-consistency eqs @ output
!       real,    optional :: eps
       logical, optional :: fresh
     end subroutine LHSmRHS
  end interface

  integer, dimension(1) :: nf,ng
  integer :: i,method
  real :: fnormW,epsGW,epsLW,pmW,eps_gamma
  character(len=9) :: itertyp

  mess   = 'N/A'
  fnorm  = 1.0
  epsG   = epsL
  logfil = noutE+1
  noutE2 = noutE+2
  noutG  = noutE+3
  noutF  = noutE+4
  
!  open(logfil,file=trim(filnamn)//'.log')

!  write(logfil,*) 'Iterate the self-consistency equations: '//trim(mess)
!  call flush(logfil)

  method=RELAX
  call scs(START,method,m_max,itertyp,vopt,pm,pm_max,fnorm,epsG,epsL,nf,ng,SESp)

  i=1

  fnormW=2*fnorm
  epsGW=0.0
  epsLW=0.0
  pmW=pm

  if(sopt) call SCSflush(noutG,filnamn,'G',xp,SESp)

  if(vEopt)then
     write(*,'(a21$)') 'eps_max'
     call Eflush('o',6,filnamn,i,itertyp,pmW,fnormW,epsGW,epsLW,nf,ng)
  end if
  call Eflush ('o',noutE, filnamn,i,itertyp,pmW,fnormW,epsGW,epsLW,nf,ng)
  call E2flush('o',noutE2,filnamn,i,fnormW,epsGW,epsLW)

!  do while((fnormW>fnorm.or.epsGW>epsG.or.epsLW>epsL).and.i<itermax)
!  do while((epsGW>epsG.or.epsLW>epsL).and.i<itermax)
  do
     call LHSmRHS(SESp,fresh=(i.eq.1))

!     if(sopt.and.mod(i,20).eq.0) call SCSflush(noutF,filnamn,'F',xp,SESp)
     if(sopt) call SCSflush(noutF,filnamn,'F',xp,SESp)

     if(i>nrr)then
        method=MOBESC
     else
        method=RELAX
     end if

     call scs(WORK,method,m_max,itertyp,vopt,pmW,pm_max,fnormW,epsGW,epsLW,nf,ng,SESp)

!     if(sopt.and.mod(i,20).eq.0) call SCSflush(noutG,filnamn,'G',xp,SESp)
     if(sopt) call SCSflush(noutG,filnamn,'G',xp,SESp)

     if(vEopt)then
        write(*,'(e21.12e3$)') eps_gamma
        call Eflush('w',6,filnamn,i,itertyp,pmW,fnormW,epsGW,epsLW,nf,ng)
     end if
     call Eflush ('w',noutE, filnamn,i,itertyp,pmW,fnormW,epsGW,epsLW,nf,ng)
     call E2flush('w',noutE2,filnamn,i,fnormW,epsGW,epsLW)

     if( (epsGW<epsG.and.epsLW<epsL).or.i.ge.itermax ) exit
     
     i=i+1
  end do

!  if(fnormW<fnorm.and.epsGW<epsG.and.epsLW<epsL)then
  if(epsGW<epsG.and.epsLW<epsL)then
     if(vopt) write(*,*) 'Reached self-consistency for: '//trim(mess)
!     write(logfil,*) 'Reached self-consistency for: '//trim(mess)
     success = .true.
  else
     if(vopt) write(*,*) 'Reached max number of iterations - NO self-consistency for: '//trim(mess)
!     write(logfil,*) 'Reached max number of iterations - NO self-consistency for: '//trim(mess)
     success = .false.
  end if

  if(vopt) call finalflush(6,filnamn,i,itertyp,pmW,fnormW,epsGW,epsLW,nf,ng)
!  call finalflush(logfil,filnamn,i,itertyp,pmW,fnormW,epsGW,epsLW,nf,ng)

  call SCSflush(noutG,filnamn,'G',xp,SESp)

  call Eflush ('c',noutE, filnamn,i,itertyp,pmW,fnormW,epsGW,epsLW,nf,ng)
  call E2flush('c',noutE2,filnamn,i,fnormW,epsGW,epsLW)

  call scs(STOPP,0,m_max,itertyp,vopt,pmW,pm_max,fnormW,epsGW,epsLW,nf,ng,SESp)

  ! return self-consisency parameters:
  itermax = i
  epsL = epsLW
!-----------------------------------------------------------------------------!
end subroutine run_scsolver

!-----------------------------------------------------------------------------!
subroutine run_mixed_scsolver(LHSmRHS,logfil,filnamn,mess,xp,SESp,&
     &itermax,m_max,nrr,pm,pm_max,fnorm,epsG,epsL,vopt,sopt,noutE,noutE2,noutG,noutF,&
     &success,iter_cycle_length,iter_cycle)
!-----------------------------------------------------------------------------!

  implicit none

  integer :: logfil ! logfile number where messages are flushed
  integer :: itermax ! max number of iterations
  integer :: m_max ! number of deviation vectors kept in memory
  integer :: nrr ! number of initial simple relaxation iterations to be done
  integer :: noutE,noutE2,noutG,noutF ! file numbers for flushing errors and vectors
  real :: pm ! initial back-coupling strength
  real :: pm_max ! max allowed back-coupling strength
  real :: fnorm,epsG,epsL ! required norm and errors
  logical :: vopt,sopt ! verbose and save flags
  character(len=255) :: filnamn,mess ! basic filename, message used in logfile
  real, dimension(:), pointer :: xp ! "grid"
  complex, dimension(:), pointer :: SESp ! initial guess of self-energy
                                         ! --> self-consistent self-energy
  logical, intent(out) :: success ! flag

  integer, intent(in) :: iter_cycle_length
  integer, dimension(iter_cycle_length), intent(in) :: iter_cycle
  
  interface
     subroutine LHSmRHS(SESp)
       implicit none
       complex, dimension(:), pointer :: SESp ! selfenergies @ input, self-consistency eqs @ output
     end subroutine LHSmRHS
  end interface

  integer, dimension(1) :: nf,ng
  integer :: i,method,i_cycle
  real :: fnormW,epsGW,epsLW,pmW
  character(len=9) :: itertyp

!  open(logfil,file=trim(filnamn)//'.log')

!  write(logfil,*) 'Iterate the self-consistency equations: '//trim(mess)
!  call flush(logfil)

  call scs(START,0,m_max,itertyp,vopt,pm,pm_max,fnorm,epsG,epsL,nf,ng,SESp)

  i=1
  i_cycle=1

  fnormW=2*fnorm
  epsGW=0.0
  epsLW=0.0
  pmW=pm

!  if(sopt) call SCSflush(noutG,filnamn,'G',xp,SESp)

  if(vopt) call Eflush('o',6,filnamn,i,itertyp,pmW,fnormW,epsGW,epsLW,nf,ng)
  call Eflush ('o',noutE, filnamn,i,itertyp,pmW,fnormW,epsGW,epsLW,nf,ng)
  call E2flush('o',noutE2,filnamn,i,fnormW,epsGW,epsLW)

!  do while((fnormW>fnorm.or.epsGW>epsG.or.epsLW>epsL).and.i<itermax)
!  do while((epsGW>epsG.or.epsLW>epsL).and.i<itermax)
  do

     call LHSmRHS(SESp)

     if(sopt.and.mod(i,5).eq.0) call SCSflush(noutF,filnamn,'F',xp,SESp)

     if(i>nrr)then
        method=MOBESC
     else
        method=RELAX
     end if

     if( iter_cycle(i_cycle).eq.RELAX ) pmW=pm

     call scs(WORK,iter_cycle(i_cycle),m_max,itertyp,vopt,pmW,pm_max,fnormW,epsGW,epsLW,nf,ng,SESp)

     if(sopt.and.mod(i,5).eq.0) call SCSflush(noutG,filnamn,'G',xp,SESp)

     if(vopt) call Eflush('w',6,filnamn,i,itertyp,pmW,fnormW,epsGW,epsLW,nf,ng)
     call Eflush ('w',noutE, filnamn,i,itertyp,pmW,fnormW,epsGW,epsLW,nf,ng)
     call E2flush('w',noutE2,filnamn,i,fnormW,epsGW,epsLW)

     if( (epsGW<epsG.and.epsLW<epsL).or.i>itermax ) exit

     i_cycle=i_cycle+1
     if( i_cycle.gt.iter_cycle_length ) i_cycle=1

     i=i+1
  end do

!  if(fnormW<fnorm.and.epsGW<epsG.and.epsLW<epsL)then
  if(epsGW<epsG.and.epsLW<epsL)then
     if(vopt) write(*,*) 'Reached self-consistency for: '//trim(mess)
!     write(logfil,*) 'Reached self-consistency for: '//trim(mess)
     success = .true.
  else
     if(vopt) write(*,*) 'Reached max number of iterations - NO self-consistency for: '//trim(mess)
!     write(logfil,*) 'Reached max number of iterations - NO self-consistency for: '//trim(mess)
     success = .false.
  end if

  if(vopt) call finalflush(6,filnamn,i,itertyp,pmW,fnormW,epsGW,epsLW,nf,ng)
!  call finalflush(logfil,filnamn,i,itertyp,pmW,fnormW,epsGW,epsLW,nf,ng)

  call SCSflush(noutG,filnamn,'G',xp,SESp)

  call Eflush ('c',noutE, filnamn,i,itertyp,pmW,fnormW,epsGW,epsLW,nf,ng)
  call E2flush('c',noutE2,filnamn,i,fnormW,epsGW,epsLW)

  call scs(STOPP,0,m_max,itertyp,vopt,pmW,pm_max,fnormW,epsGW,epsLW,nf,ng,SESp)

!-----------------------------------------------------------------------------!
end subroutine run_mixed_scsolver

!-----------------------------------------------------------------------------!
end module scsolver
