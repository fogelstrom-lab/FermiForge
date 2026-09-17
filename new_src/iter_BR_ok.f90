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
