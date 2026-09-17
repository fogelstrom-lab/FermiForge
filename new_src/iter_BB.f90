!---------------------------------------------------------------------
!
   subroutine iterateBB(tolA,tolL,alpha)
!
!---------------------------------------------------------------------
!
   implicit none

   real, intent(out) :: tolA, tolL, alpha

   integer, parameter :: mvec=2
   integer :: N, i, ii
   integer, dimension(nop) :: lb, ub

   real ::  delta, eps1, eps2, eps3, eps4,norm
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
         N = ub(nop)

         allocate (v(N,mvec),f(N,mvec),h(N))
         allocate (epsq(mvec))
         epsq = 0.0
         v = 0.0
         f = 0.0

         allocate (d(N))
         d = 0.0

         ii=1
         alpha = 1.0
      endif

      if(iiter > 1) then
! -- Save the latest relevant iterates

         v(1:N,2) = v(1:N,1)   ! The vector v
         f(1:N,2) = f(1:N,1)   ! The residue f(v)=g(v)-v
         epsq(2) = epsq(1)
      endif

      call pacvector(lb,ub,h,1)
      norm = sqrt(abs(dot_product(h,h)))
      v(1:N,1) = h
   end if

   call getnewop

   if(myid == 0) then

      call pacvector(lb,ub,h,2)
      f(1:N,1) = h

      if(iiter == 1) then
         eps1 = abs(dot_product(f(1:N,1),f(1:N,1)))
         epsq(1) = eps1

! -- Make the first updates by simple update
 
         v(1:N,1) = v(1:N,1) + f(1:N,1)  ! simple update first time around
         
         tolA = sqrt(eps1)/norm
         tolL = N*maxval(abs(f(1:N,1)))/norm

      else

         d = f(1:N,2)-f(1:N,1)  ! the difference in residues
         eps1=abs(dot_product(f(1:N,1),f(1:N,1)))
         epsq(1) = eps1

! -- Simple implementation of Barzilai–Borwein iteration algoritm

         eps2=abs(dot_product(f(1:N,2),f(1:N,2)))
         eps3=abs(dot_product(d,d))
         eps4=abs(dot_product(d,f(1:N,2)))

         if(mod(ii,2) /= 0) delta=eps2/eps4
         if(mod(ii,2) == 0) delta=eps4/eps3
         alpha=alpha*delta
         if(alpha > 100.0) alpha=100.0
         if(alpha < 0.001) alpha=0.001

         v(1:N,1) = v(1:N,1) + alpha * f(1:N,1)

      end if
!
! -- The next selfenergies
!
      call pacvector(lb,ub,v(1:N,1),3)

      ii=ii+1
      
      tolA = sqrt(eps1)/norm
      tolL = N*maxval(abs(f(1:N,1)))/norm
   end if
!
!---------------------------------------------------------------------
!
   end subroutine iterateBB
!
!---------------------------------------------------------------------
