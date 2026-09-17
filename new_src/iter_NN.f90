!---------------------------------------------------------------------
!
   subroutine iterateNN(tolA,tolL)   ! Most simple update x_new = g(x_old)
!
!---------------------------------------------------------------------
!
   implicit none

   real, intent (out) :: tolA, tolL
   integer :: N, i
   real norm
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
         N = ub(nop)
         allocate (x(N),f(N),h(N))
      endif 

      call pacvector(lb,ub,h,1)
      x = h
      norm = sqrt(abs(dot_product(x,x)))
   end if

   call getnewop

   if(myid == 0) then
      call pacvector(lb,ub,h,2)
      f = h

      x = x + f
      call pacvector(lb,ub,x,3)

      tolA = sqrt(abs(dot_product(f,f)))/norm
      tolL = float(N)*maxval(abs(f))/norm
   end if
!
!---------------------------------------------------------------------
!
   end subroutine iterateNN
!
!-----------------------------------------------------------
