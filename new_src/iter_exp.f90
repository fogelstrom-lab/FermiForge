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
   subroutine iterateSE(tol,norm,p,md,AA,BR,BB,NN)
!
!---------------------------------------------------------------------
!
      integer, intent (inout) :: md
      real, intent (inout) :: tol, norm, p 
      logical,intent (in)  :: AA,BR,BB,NN

      complex :: cpp,cpo,cpm,cop,coo,com,cmp,cmo,cmm

      md = 0
      norm = 0.0
      if(AA) call iterateAA(tol,norm,p,md)   ! Anderson mixing scheme
      if(BR) call iterateBR(tol)             ! Broyden update scheme
      if(BB) call iterateBB(tol,p)           ! Barzilai-Borwein update scheme
      if(NN) call iterateNN(tol)             ! simple update x_new = g(x_old)  

      call SES_bcast

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
   subroutine iterateNN(eps)   ! Most simple update x_new = g(x_old)
!
!---------------------------------------------------------------------
!
   real, intent (out) :: eps
   integer :: i
   integer, dimension(nop) :: lb, ub
   real, dimension(:), allocatable :: x, f, h
   save

   if(myid == 0) then
      if(iiter == 1) then
         lb(1) = 1
         ub(1) = vecl
         do i = 2,  nop, 1
            lb(i)= lb(i-1)+ub(1)
            ub(i)= ub(i-1)+ub(1)
         end do
         allocate (x(ub(nop)),f(ub(nop)),h(ub(nop)))
      endif 

      call pacvector(lb,ub,h,1)
      x = h
   end if

   call getnewop

   if(myid == 0) then
      call pacvector(lb,ub,h,2)
      f = h

      eps = sqrt(abs(dot_product(f,f)))

      x = x + f
      call pacvector(lb,ub,x,3)
   end if

! -- send the error and end of the updating routine

   call MPI_BCAST(eps,1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
!
!---------------------------------------------------------------------
!
   end subroutine iterateNN
!
!-----------------------------------------------------------
!
   subroutine iterateBR(eps)
!
!-----------------------------------------------------------
!
      real, intent(out) :: eps

      integer, dimension(nop) :: lb, ub
      real, dimension(:,:), allocatable :: jac,f
      real, dimension(:), allocatable :: v, dv, h, g, a
      real :: alpha
      integer :: i, j

      save 

      if(myid == 0) then

         if(iiter== 1) then

            lb(1) = 1
            ub(1) = vecl
            do i = 2,  nop, 1
               lb(i)= lb(i-1)+ub(1)
               ub(i)= ub(i-1)+ub(1)
            end do

            allocate (jac(ub(nop),ub(nop)),f(ub(nop),2))
            allocate (v(ub(nop)),dv(ub(nop)),h(ub(nop)),g(ub(nop)),a(ub(nop)))

            f = 0.0
            v = 0.0
            dv = 0.0
            h = 0.0
            a = 0.0
            g = 0.0

     ! Initial Jacobian approximation (identity matrix)

            jac = 0.0
            do i = 1, ub(nop)
               jac(i,i) = -1.0
            end do

         end if

         call pacvector(lb,ub,h,1)
         v = h
      endif  

      call getnewop
      
      if(myid == 0) then

         f(:,2) = f(:,1)
         call pacvector(lb,ub,h,2) 
         f(:,1) = h

         if(iiter > 1) then

     ! Update the Jacobian approximation using Broyden's formula 
     ! using simplest algorithm from Fang and Saad, Num. Lin. Algebra Appl. 2009: 16, 197-221

            h = f(:,1)-f(:,2)

            do i = 1, ub(nop), 1
               a(i) = dot_product(dv,jac(:,i))   ! Type-I Broyden 
!              a(i) = h(i)                       ! Type-II Broyden 
            end do

            alpha =  dot_product(a,h)
            alpha = 1.0 / alpha

            do i = 1, ub(nop), 1
               g(i) = dot_product(jac(i,:),h) 
            end do

            g = alpha * (dv - g)

            do i = lb(1), ub(nop)
               do j = lb(1), ub(nop)
                  jac(i,j) = jac(i,j) +  g(i) * a(j)
               end do 
            end do 
         end if

         do i = 1, ub(nop), 1 
            dv(i) = - dot_product(jac(i,:),f(:,1))
         end do

         eps = sqrt(dot_product(f(:,1),f(:,1)))
         v = v + dv
         call pacvector(lb,ub,v,3)

      end if

! -- send the error and end of the updating routine

      call MPI_BCAST(eps,1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
!
!-----------------------------------------------------------
!
    end subroutine iterateBR
!
!---------------------------------------------------------------------
!
   subroutine iterateBB(tol,alpha)
!
!---------------------------------------------------------------------
!
   real, intent(out) :: tol, alpha

   integer, parameter :: mvec=2
   integer :: i, ii
   integer, dimension(nop) :: lb, ub

   real ::  delta, eps1, eps2, eps3, eps4
   real, dimension(:,:), allocatable :: v, f
   real, dimension(:), allocatable :: d, epsq, h
   save 

   if(myid == 0) then
      if(iiter == 1) then
         lb(1) = 1
         ub(1) = vecl
         do i = 2,  nop, 1
            lb(i)= lb(i-1)+ub(1)
            ub(i)= ub(i-1)+ub(1)
         end do

         allocate (v(ub(nop),mvec),f(ub(nop),mvec),h(ub(nop)))
         allocate (epsq(mvec))
         epsq = 0.0
         v = 0.0
         f = 0.0

         allocate (d(ub(nop)))
         d = 0.0

         ii=1
         alpha = 1.0
      endif

      if(iiter > 1) then
! -- Save the latest relevant iterates

         v(:,2) = v(:,1)   ! The vector v
         f(:,2) = f(:,1)   ! The residue f(v)=g(v)-v
         epsq(2) = epsq(1)
      endif

      call pacvector(lb,ub,h,1)
      v(:,1) = h
   end if

   call getnewop

   if(myid == 0) then

      call pacvector(lb,ub,h,2)
      f(:,1) = h

      if(iiter == 1) then
         eps1 = abs(dot_product(f(:,1),f(:,1)))
         epsq(1) = eps1

! -- Make the first updates by simple update
 
         v(:,1) = v(:,1) + f(:,1)  ! simple update first time around
         tol = sqrt(eps1)
      else

         d = f(:,2)-f(:,1)  ! the difference in residues
         eps1=abs(dot_product(f(:,1),f(:,1)))
         epsq(1) = eps1

! -- Simple implementation of Barzilai–Borwein iteration algoritm

         eps2=abs(dot_product(f(:,2),f(:,2)))
         eps3=abs(dot_product(d,d))
         eps4=abs(dot_product(d,f(:,2)))

         if(mod(ii,2) /= 0) delta=eps2/eps4
         if(mod(ii,2) == 0) delta=eps4/eps3
         alpha=alpha*delta
         if(alpha > 100.0) alpha=100.0
         if(alpha < 0.001) alpha=0.001

         v(:,1) = v(:,1) + alpha * f(:,1)

      end if
!
! -- The next selfenergies
!
      call pacvector(lb,ub,v(:,1),3)

      ii=ii+1
      tol = sqrt(eps1)
   end if

! -- send the error and end of the updating routine

   call MPI_BCAST(tol,1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
!
!---------------------------------------------------------------------
!
   end subroutine iterateBB
!
!---------------------------------------------------------------------
!
   subroutine iterateAA(tol,norm,p,md)
!
!---------------------------------------------------------------------
!
   integer, intent (out) :: md
   real, intent (out) :: tol, p, norm
       
   integer, parameter :: mmax = 6    ! 6 seems optimal for now
   integer :: i, j, ju, m, lda, info
   integer, dimension(nop) :: lb, ub

   logical :: swap, flip
   real :: q, pq, eps, pmax, prog, eps_best

   complex, dimension(:,:), allocatable :: ca
   complex, dimension(:), allocatable :: cb

   real, dimension(:,:), allocatable :: v, f, a
   real, dimension(:), allocatable :: work, dv, df, b, epsq, h, d, v_best
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

         allocate (v(ub(nop),0:mmax+1),f(ub(nop),0:mmax+1))
         allocate (dv(ub(nop)),df(ub(nop)),h(ub(nop)),d(ub(nop)),v_best(ub(nop)))
         allocate (epsq(0:mmax+1))

         epsq = 0.0
         x = 0.0
         f = 0.0


         pmax = 1.0   ! 1.0 seems optimal for mmax=6
         prog = 0.01
         pq=1.0
         p=0.1
         m=0
         usub=0.0
         eps_best =  1.0e3
      end if

      call pacvector(lb,ub,h,1)
      v(:,m) = h
 
   end if

   call getnewop

   if(myid == 0) then
 
      call pacvector(lb,ub,h,2)
      f(:,m) = h
      eps = abs(dot_product(h,h))

      epsq(m) = eps 
      tol = sqrt(eps)
      ju = 0
!
! -- save the best vector v as a go-to vector
!
      if(eps < eps_best) then
         v_best = v(:,m)
         eps_best = eps
      end if
!
!      case m <= 1
!
      if(m <= 1) then
         v(:,m) = v(:,m) + p * f(:,m)   ! Simple update
         norm=1.0
         print *, 'error vector (1) :', iiter,m,ju,sqrt(epsq(0:m))
         usub = epsq(m)
         flip=.false.
      else
!
!     case m > 1
!
! sort the vectors in decending order in the error
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
         if(m == mmax+1)              ju = 1      ! just shift the solution index by one
         if(eps >  epsq(m))           ju=(m+1)/2  ! very bad step, throw away the worst half of the vectors
         if(epsq(m) < 0.01*epsq(m-1)) ju=(m+1)/2  ! very good step, throw away the worst half of the vectors
         if(ju == 0) flip = .false.

         if( ju /= 0) then

            if(epsq(m) < 0.01*epsq(m-1))    flip = .false.
            if(ju == 1 .and. eps <= epsq(m)) flip = .false.

            if(flip) then   ! after two bad updates in a row use the go-to vector and restart
               v = 0.0
               f = 0.0
               epsq = 0.0
               m = 0
               v(:,0) = v_best
               epsq(0) = eps_best
               p = 0.01
               print *, 'Bad updates, need to do a restart...',sqrt(eps),sqrt(eps_best)
               goto 100
            else
               do i = ju, m, 1 
                  v(:,i-ju) = v(:,i)
                  f(:,i-ju) = f(:,i)
                  epsq(i-ju) = epsq(i)
               end do 
               m=m-ju
               if(eps > epsq(m)) flip = .true.
            end if

         end if
         print *, 'error vector (2) :', iiter,m,ju,sqrt(epsq(0:m))

!     -----------------------------------------------------------
!
!      Build the matrix (f - f )(f - f ) and the vector -f (f - f ).
!                         i   m   k   m                   m  k   m 
!     -----------------------------------------------------------
 
         lda = m
         allocate(a(lda,lda),b(lda),work(lda))
         allocate(ca(lda,lda),cb(lda))

         do i = 1, lda, 1
!           d = v(:,i-1)-v(:,m)   ! Type-I Anderson (does not work/behave well)
            d = f(:,i-1)-f(:,m)   ! Type-II Anderson
            do j = 1, lda, 1
               h = f(:,j-1)-f(:,m)
               a(i,j) = dot_product(d,h)
            end do
            a(i,i) = 2.0 * a(i,i)
            b(i)= 2.0 * dot_product(d,f(:,m))
!           a(i,i) = a(i,i)
!           b(i)= dot_product(d,f(:,m))
         end do

         a = a / eps
         b = b / eps
         do i = 1, lda, 1
            cb(i) =  cmplx(b(i),0.0)
            do j=i, lda, 1 
               ca(i,j) =  cmplx(a(i,j),0.0)
            end do
         end do
 
!      ----------------------------------------------------------
!                                                    m-1
!      finding the minimum of |f|**2 for f(g) = f  - sum g *(f -f ).
!                                                m   i=0  i   i  m
!      ----------------------------------------------------------
 
         call dgesv(lda, 1, a, lda, work, b, lda, info)
!        call zposv('U', lda, 1, ca, lda, cb, lda, info)

         if(info /= 0) print *, 'Problem with the matrix a',info
!        b = real(cb)

         norm = sqrt(abs(dot_product(b,b)))

         dv = 0.0
         df = 0.0
         do i = 1, lda, 1 
            dv=dv+(v(:,i-1)-v(:,m))*b(i)
            df=df+(f(:,i-1)-f(:,m))*b(i)
         end do
         usub = dot_product(df,df)

!         pq=min(1.0,usub/eps)

         deallocate(a,b,ca,cb,work)
 
!     -----------------------------------------------
!     finding the p_m  ! Replaced the origina p1 suboutine
!     -----------------------------------------------
 
          q=p
          p=sqrt(dot_product(dv,dv)/eps)
          p=pq*(p+q)*0.5
!         p=pq*p
          p=min(pmax,p)
          p=max(p,0.01)
 
!     -----------------------------------------------
!     do the mixing
!     -----------------------------------------------
!
         v(:,m) = v(:,m) - dv + p * (f(:,m) - df)
  
      endif
!
! -- return the new self energies and add one to the counter m if needed

 100  md = m
      h= v(:,m)
      call pacvector(lb,ub,h,3)
      eps = tol 
      m=m+1
   endif

! -- send the error and end of the updating routine

   call MPI_BCAST(tol,1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(norm,1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)

!-----------------------------------------------------------------------
!
   end subroutine iterateAA
!
!-----------------------------------------------------------------------
!
end module Iterate
!
!-----------------------------------------------------------------------
