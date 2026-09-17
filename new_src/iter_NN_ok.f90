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
