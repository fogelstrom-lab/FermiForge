!---------------------------------------------------------------------
!
   subroutine iterateAA(tolA,tolL,p,md)
!
!---------------------------------------------------------------------
!
   implicit none

   integer, intent (out) :: md
   real, intent (out) :: tolA, tolL, p
       
   integer, parameter :: mmax = 6
   integer :: N, i, j, ju, m, lda, info
   integer, dimension(nop) :: lb, ub

   logical :: swap, progress
   real :: q, pq, eps, pmax, prog, usub, norm

   real, dimension(:,:), allocatable :: v, f, a
   real, dimension(:), allocatable :: work, dv, df, b, epsq, h, d
   real, dimension(:), allocatable :: fmin, fmin_old
   save
!
!-----------------------------------------------------------------------
!
   if(myid == 0) then
      if(iiter == 1) then
         lb(1) = 1
         ub(1) = vecl
         do i = 2,  nop, 1
            lb(i)= lb(i-1)+ub(1)
            ub(i)= ub(i-1)+ub(1)
         end do
         N= ub(nop)

         allocate (v(N,0:mmax+1),f(N,0:mmax+1))
         allocate (dv(N),df(N),h(N),d(N),fmin_old(N),fmin(N))
         allocate (epsq(0:mmax+1))

         epsq = 0.0
         v = 0.0
         f = 0.0

         pmax = 1.0
         prog = 0.01
         pq=1.0
         p=0.01
         m=0
         usub=0.0

      end if

      call pacvector(lb,ub,h,1)
      v(:,m) = h
      norm = dot_product(h,h)   ! norm of the current self-energy vector
   end if

   call getnewop

   if(myid == 0) then
 
      call pacvector(lb,ub,h,2)
      f(:,m) = h
      eps = abs(dot_product(h,h))

      epsq(m) = eps 
      tolA = sqrt(eps)/norm
      tolL = float(N)*maxval(abs(h))/norm 
      ju = 0
!
!      case m <= 1
!
      if(m <= 1) then
         p = 1.0
         v(:,m) = v(:,m) + p * f(:,m)   ! Simple update
         usub = eps
         fmin = f(:,m)
         progress = .true.
!        print *, 'error vector (1) :', iiter,m,ju,sqrt(epsq(0:m))
      else
!
!     case m > 1
!
! (bubble) sort the vectors in decending order in the error 
!
         do 
            swap=.false.
            do i= m-1, 0, -1
               if (epsq(i) < epsq(i+1)) then
                  q = epsq(i+1)
                  dv = v(:,i+1)
                  df = f(:,i+1)
                  epsq(i+1) = epsq(i)
                  v(:,i+1) =  v(:,i)
                  f(:,i+1) =  f(:,i)
                  epsq(i) = q
                  v(:,i) = dv
                  f(:,i) = df

                  swap=.true.
               end if
            end do
            if(.not.swap) exit
         end do
!
! -- make the neccecary shifts
! 
         if(m == mmax+1) ju = 1
         if(eps > epsq(m))  ju=(m+1)/2             ! very bad step
         if(epsq(m) < 0.01*epsq(m-1)) ju=(m+1)/2   ! very good step

         if( ju /= 0) then

            do i = ju, m, 1 
               v(:,i-ju) = v(:,i)
               f(:,i-ju) = f(:,i)
               epsq(i-ju) = epsq(i)
            end do 
            m=m-ju
            progress = .true.
            if(m < 2) then        ! restart
               v(:,0) = v(:,m) + f(:,m)
               f(:,0) = 0.0
               epsq(0) = epsq(m)
               v(:,1:m) = 0.0
               f(:,1:m) = 0.0
               epsq(1:m) = 0.0
               m = 0
               progress = .false.
               goto 100
            endif
         end if
!        print *, 'error vector (2) :', iiter,m,ju,sqrt(epsq(0:m))

!     -----------------------------------------------------------
!
!      Build the matrix (f - f )(f - f ) and the vector -f (f - f ).
!                         i   m   k   m                   m  k   m 
!     -----------------------------------------------------------
 
         lda = m
         allocate(a(lda,lda),b(lda),work(lda))

         do i = 1, lda, 1
            d = f(:,i-1)-f(:,m)   ! Type-II Anderson
            do j = 1, lda, 1
               h = f(:,j-1)-f(:,m)
               a(i,j) = dot_product(d,h)
            end do
            a(i,i) = 2.0 * a(i,i)
            b(i)=-2.0 * dot_product(d,f(:,m))
         end do

         a = a / eps
         b = b / eps

!      ----------------------------------------------------------
!                                                    m-1
!      finding the minimum of |f|**2 for f(g) = f  - sum g *(f -f ).
!                                                m   i=0  i   i  m
!      ----------------------------------------------------------
 
         call dgesv(lda, 1, a, lda, work, b, lda, info)

         if(info /= 0) print *, 'Problem with the matrix a',info

         dv = 0.0
         df = 0.0
         do i = 1, lda, 1 
            dv=dv+(v(:,i-1)-v(:,m))*b(i)
            df=df+(f(:,i-1)-f(:,m))*b(i)
         end do


         deallocate(a,b,work)
 
!     -----------------------------------------------
!     finding the p_m  ! Replaced the origina p1 suboutine
!     -----------------------------------------------
 
         
         fmin_old = fmin
         usub = dot_product(fmin_old,fmin_old) 
         fmin = f(:,m) + df

         pq = min(1.0 , usub/eps)
         q=p
         p=sqrt(dot_product(dv,dv)/eps)
         p=pq*(p+q)*0.5
         p=min(pmax,p)
         p=max(p,0.01)
 
!     -----------------------------------------------
!     do the mixing
!     -----------------------------------------------
!
         v(:,m) = v(:,m) + dv + p * fmin

! --  exclude useless points

         if(usub/eps < prog .and. m > 2) then
            print *, ' throwing vectors ', usub,eps
            ju = (m+1)/2
            do i = ju, m, 1
               v(:,i-ju) = v(:,i)
               f(:,i-ju) = f(:,i)
               epsq(i-ju) = epsq(i)
            end do
            m=m-ju
            progress = .true.
         end if

      end if
 
! -- return the new self energies and add one to the counter m if needed

 100  md = m
      h= v(:,m)
      call pacvector(lb,ub,h,3)
      if(progress) m = m + 1
   endif

! --  send the error and end of the updating routine

!  call MPI_BCAST(tol,1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(norm,1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)

!-----------------------------------------------------------------------
!
   end subroutine iterateAA
!
!-----------------------------------------------------------------------
