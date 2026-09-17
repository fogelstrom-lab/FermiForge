MODULE Iterate
   use Global_variables
   use MPI_variables
   use mpicalls 
   use newses
contains
!
!---------------------------------------------------------------------
!
   subroutine iterateSE(tol,norm,p,md,BR,BB,NN)
!
!---------------------------------------------------------------------
!
      integer, intent (inout) :: md
      real, intent (inout) :: tol, norm, p 
      logical,intent (in)  :: BR,BB,NN

      complex :: cpp,cpo,cpm,cop,coo,com,cmp,cmo,cmm

      md = 0
      norm = 0.0
      if(BR) call iterateBR(tol)             ! Broyden update scheme
      if(BB) call iterateBB(tol,p)           ! Barzilai-Borwein update scheme
      if(NN) call iterateNN(tol)             ! simple update x_new = g(x_old)  

      call SES_bcast

      if(myid == 0) then

         open(1,file='op_xyz',status='unknown')
         open(2,file='op_harm',status='unknown')
         open(3,file='curr',status='unknown')

         do ix=0,nx,1
            x=float(ix)*dx

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
   integer, parameter :: nop=21
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

      integer, parameter :: nop=21
      integer, dimension(nop) :: lb, ub
      real, dimension(:,:), allocatable :: jac,f
      real, dimension(:), allocatable :: x, dx, h, g, a
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
            allocate (x(ub(nop)),dx(ub(nop)),h(ub(nop)),g(ub(nop)),a(ub(nop)))

            f = 0.0
            x = 0.0
            dx = 0.0
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
         x = h
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
               a(i) = dot_product(dx,jac(:,i)) 
            end do

            alpha =  dot_product(a,h)
            alpha = 1.0 / alpha

            do i = 1, ub(nop), 1
               g(i) = dot_product(jac(i,:),h) 
            end do

            g = alpha * (dx - g)

            do i = lb(1), ub(nop)
               do j = lb(1), ub(nop)
                  jac(i,j) = jac(i,j) +  g(i) * a(j)
               end do 
            end do 
         end if

         do i = 1, ub(nop), 1 
            dx(i) = - dot_product(jac(i,:),f(:,1))
         end do

         eps = sqrt(dot_product(dx,dx))
         x = x + dx
         call pacvector(lb,ub,x,3)

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

   integer, parameter :: mvec=2, nop=21
   integer :: i, ii
   integer, dimension(nop) :: lb, ub

   real ::  delta, eps1, eps2, eps3, eps4
   real, dimension(:,:), allocatable :: x, f
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

         allocate (x(ub(nop),mvec),f(ub(nop),mvec),h(ub(nop)))
         allocate (epsq(mvec))
         epsq = 0.0
         x = 0.0
         f = 0.0

         allocate (d(ub(nop)))
         d = 0.0

         ii=1
         alpha = 1.0
      endif

      if(iiter > 1) then
! -- Save the latest relevant iterates

         x(:,2) = x(:,1)   ! The vector x
         f(:,2) = f(:,1)   ! The residue f(x)=g(x)-x
         epsq(2) = epsq(1)
      endif

      call pacvector(lb,ub,h,1)
      x(:,1) = h
   end if

   call getnewop

   if(myid == 0) then

      call pacvector(lb,ub,h,2)
      f(:,1) = h

      if(iiter == 1) then
         eps1 = abs(dot_product(f(:,1),f(:,1)))
         epsq(1) = eps1

! -- Make the first updates by simple update
 
         x(:,1) = x(:,1) + f(:,1)  ! simple update first time around
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

         x(:,1) = x(:,1) + alpha * f(:,1)

      end if
!
! -- The next selfenergies
!
      call pacvector(lb,ub,x(:,1),3)

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
! This subroutine does not work (20250106)
!
!---------------------------------------------------------------------
!
   integer, intent (out) :: md
   real, intent (out) :: tol, p, norm
       
   integer, parameter :: nop=21, mmax = 6
   integer :: i, j, m,  ju, lda, info
   integer, dimension(nop) :: lb, ub
   real :: q, pq, usub, duu, eps, v, cccc, pmax, prog, u1

   real, dimension(:,:), allocatable :: x, f, ga, a
   real, dimension(:), allocatable :: work, dx, df, b, g, epsq, h, d
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

         allocate (x(ub(nop),0:mmax),f(ub(nop),0:mmax),h(ub(nop)),d(ub(nop)))
         allocate (epsq(0:mmax))
         epsq = 0.0
         x = 0.0
         f = 0.0

         allocate(a(0:mmax,0:mmax),b(0:mmax))
         allocate (dx(ub(nop)),df(ub(nop)))

         pmax = 3.0
         prog = 0.01
         pq=1.0
         p=0.5
         m=0
         usub=0.0
      end if
!
!---------------------------------------------------------------------
!
!     storage x ==> x  and calculation of f  and of f *f , i=1,...,m
!                    m                     m         i  m
!
!-----------------------------------------------------------------------
!
      call pacvector(lb,ub,h,1)
      x(:,m) = h
   end if

   call getnewop

   if(myid == 0) then
      call pacvector(lb,ub,h,2)
      f(:,m) = h

      eps= abs(dot_product(f(:,m),f(:,m)))
      epsq(m)= eps
      tol = sqrt(eps)

      if(tol.le.errtol) then
!
!----------------------------------------------------------------------
!
!       convergence: restoration of x and return
!
!-----------------------------------------------------------------------
!
         x(:,m) = x(:,m) + f(:,m)
         call pacvector(lb,ub,x(:,m),3)
         eps = tol
         goto 100
      end if
!
!-----------------------------------------------------------------------
!
!     exclusion of useless points (note that m may be changed!)
!
!-----------------------------------------------------------------------
!
      ju=0
      if(m == mmax) ju=1
      if(mmax == 0) ju=0
      if(m > 1) then
         duu=epsq(0)
         if(usub/eps < prog .or. duu/eps < 1.0) ju=m/2
      endif
      if(ju.ne.0) then
         do j=ju,m
            x(:,j-ju)=x(:,j)
            f(:,j-ju)=f(:,j)
            epsq(j-ju)=epsq(j)
         enddo
         m=m-ju
         do j= m+1, mmax, 1 
            x(:,j) = 0.0
            f(:,j) = 0.0
            epsq(j) =0.0
         end do
      end if
      pq=min(1.0,sqrt(usub/eps))
!
!     -----------------------------------------------------------
!      case m=0
!     -----------------------------------------------------------
!
      if(m == 0) then

         x(:,m) = x(:,m) + p * f(:,m)   ! Simple update
         call pacvector(lb,ub,x(:,m),3)

         usub=epsq(m)
         md = m
         m = 1
         eps = tol
         norm=1.0
         goto 100
!
!     -----------------------------------------------------------
!      case m=1
!     -----------------------------------------------------------
!
      elseif(m == 1) then

         if (epsq(0) < epsq(1)) then
            v = epsq(1)
            dx = x(:,1)
            df = f(:,1)
            epsq(1) = epsq(0)
            x(:,1) =  x(:,0)
            f(:,1) =  f(:,0)
            epsq(0) = v
            x(:,0) = dx
            f(:,0) = df
         end if

         a(1,0)=dot_product(f(:,1),f(:,0))
         cccc=epsq(1)+epsq(0)-a(1,0)-a(1,0)
         u1= (max(cccc,0.0))/eps

         b(0)=(eps-a(1,0))/eps
         if(sqrt(u1) < 1e-15) stop 
         b(0)=b(0)/u1

         norm=abs(b(0))
      else
!
!     -----------------------------------------------------------
!      case m>1
!     -----------------------------------------------------------
!
! sort the vectors in decending order in the error
!                       
         do j= m-1, 0, -1
            if (epsq(j) < epsq(j+1)) then
               v = epsq(j+1)
               dx = x(:,j+1)
               df = f(:,j+1)
               epsq(j+1) = epsq(j)
               x(:,j+1) =  x(:,j)
               f(:,j+1) =  f(:,j)
               epsq(j) = v
               x(:,j) = dx
               f(:,j) = df
            end if
         end do
!
!     -----------------------------------------------------------
!
!      Build the matrix (f - f )(f - f ) and the vector -f (f - f ).
!                         i   m   k   m                   m  k   m 
!     -----------------------------------------------------------
!
         do i = 0, m-1, 1
            d = f(:,i)-f(:,m)
            do j = 0, m-1, 1
               h = f(:,j)-f(:,m)
               a(i,j) = dot_product(d,h)/eps
            end do
            a(i,i) = max(2.0*a(i,i),0.0)
            b(i) = -2.0 * dot_product(f(:,m),d)/eps
         end do
!
!      ----------------------------------------------------------
!                                                    m-1
!      finding the minimum of |f|**2 for f(g) = f  + sum g *(f -f ).
!                                                m   i=0  i   i  m
!      ----------------------------------------------------------
!      this is given by the eq.system (solve for g)
!        m-1             *  *            *  *
!        sum g *(f -f )(f -f )  =  - f (f -f )
!        i=0  i   i  m   k  m         m  k  m
!      ----------------------------------------------------------
!
         lda = m
         allocate(ga(lda,lda),g(lda),work(lda))
         do j= 0, m-1, 1 
            g(j+1) = b(j)
            do i= 0, m-1, 1 
               ga(i+1,j+1) = a(i,j)
            end do
         end do

         call dgesv(lda, 1, ga, lda, work, g, lda, info)

         norm=sqrt(abs(dot_product(g,g)))

         do i= 0, m-1
            b(i) = g(i+1)
         end do
         deallocate(ga,g,work)
      endif
!
!      ----------------------------------------------------------
!                      m-1
!     setting x = x  + sum g *(x -x ) + p *f   ,
!                  m   i=0  i   i  m     m  min
!
!                 m-1
!     where p =  |sum g *(x -x )| / |f |
!            m    i=0  i   i  m       m
!
!                      m-1
!     and   f   = f  + sum g *(f -f ) .
!            min   m   i=0  i   i  m 
!      ----------------------------------------------------------
!
      dx = 0.0
      do i = 0, m-1, 1 
         dx=dx+(x(:,i)-x(:,m))*b(i)
      enddo
!
!     -----------------------------------------------
!     finding the p_m  ! Replaced the origina p1 suboutine
!     -----------------------------------------------
!
      q=p
      p=sqrt(dot_product(dx,dx)/eps)
      p=pq*(p+q)*0.5
!     p=pq*p
      p=min(pmax,p)
      p=max(p,0.01)
!
!     -----------------------------------------------
!     calculate further the optimal czx and cx (here=f)
!     -----------------------------------------------
!
      df = 0.0
      do i = 0, m-1, 1 
         df=df+(f(:,i)-f(:,m))*b(i)
      end do
      df = f(:,m) + df
      usub=dot_product(df,df)
!
!     -----------------------------------------------
!     do the mixing
!     -----------------------------------------------
!
      x(:,m) = x(:,m) + dx + p * df
      md = m
!
!      print *, ' p=',p,' q=',q,' pq=',pq,' usub=',usub
!
!     -----------------------------------------------------------
!     exclusion of useless points (note that m may be changed!)
!     -----------------------------------------------------------
!
      if(usub/eps.lt.prog) then
         print *, usub/eps, prog
         ju=(m-1)/2
         do j= ju, m
            x(:,j-ju)=x(:,j)
            f(:,j-ju)=f(:,j)
            epsq(j-ju)=epsq(j)
         enddo
         m=m-ju
         do j= m+1, mmax, 1 
            x(:,j) = 0.0
            f(:,j) = 0.0
            epsq(j) =0.0
         end do
      end if

! -- return the new self energies and add one to the counter m

      call pacvector(lb,ub,x(:,m),3)
      m = m + 1
   endif

! -- send the error and end of the updating routine

 100 call MPI_BCAST(tol,1,MPIFLOAT,0,MPI_COMM_WORLD,ierr)
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
