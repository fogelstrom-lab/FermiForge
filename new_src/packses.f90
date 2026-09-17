MODULE packSES
   use global_variables
   use MPI_variables
contains
!
!---------------------------------------------------------------------
!
   subroutine pacvector(lb,ub,x,chr)
!
!---------------------------------------------------------------------
!
      integer, intent(in) :: chr
      integer :: i
      integer, dimension(nop), intent(in) :: lb, ub
      real, dimension(ub(nop)), intent(inout) :: x

      save 

      if(chr == 1) then
         i=0
         x(lb(1+i):ub(1+i)) =  real(dxx(0:nx))
         x(lb(2+i):ub(2+i)) = aimag(dxx(0:nx))
         x(lb(3+i):ub(3+i)) =  real(dxy(0:nx))
         x(lb(4+i):ub(4+i)) = aimag(dxy(0:nx))
         x(lb(5+i):ub(5+i)) =  real(dxz(0:nx))
         x(lb(6+i):ub(6+i)) = aimag(dxz(0:nx))
         i=6
         x(lb(1+i):ub(1+i)) =  real(dyx(0:nx))
         x(lb(2+i):ub(2+i)) = aimag(dyx(0:nx))
         x(lb(3+i):ub(3+i)) =  real(dyy(0:nx))
         x(lb(4+i):ub(4+i)) = aimag(dyy(0:nx))
         x(lb(5+i):ub(5+i)) =  real(dyz(0:nx))
         x(lb(6+i):ub(6+i)) = aimag(dyz(0:nx))
         i=12
         x(lb(1+i):ub(1+i)) =  real(dzx(0:nx))
         x(lb(2+i):ub(2+i)) = aimag(dzx(0:nx))
         x(lb(3+i):ub(3+i)) =  real(dzy(0:nx))
         x(lb(4+i):ub(4+i)) = aimag(dzy(0:nx))
         x(lb(5+i):ub(5+i)) =  real(dzz(0:nx))
         x(lb(6+i):ub(6+i)) = aimag(dzz(0:nx))
         i=18
         x(lb(1+i):ub(1+i)) =  real(vx(0:nx))
         x(lb(2+i):ub(2+i)) =  real(vy(0:nx))
         x(lb(3+i):ub(3+i)) =  real(vz(0:nx))
      end if

      if(chr == 2) then
         i=0
         x(lb(1+i):ub(1+i)) =  real(dxx(0:nx)) - x(lb(1+i):ub(1+i))
         x(lb(2+i):ub(2+i)) = aimag(dxx(0:nx)) - x(lb(2+i):ub(2+i))
         x(lb(3+i):ub(3+i)) =  real(dxy(0:nx)) - x(lb(3+i):ub(3+i))
         x(lb(4+i):ub(4+i)) = aimag(dxy(0:nx)) - x(lb(4+i):ub(4+i))
         x(lb(5+i):ub(5+i)) =  real(dxz(0:nx)) - x(lb(5+i):ub(5+i)) 
         x(lb(6+i):ub(6+i)) = aimag(dxz(0:nx)) - x(lb(6+i):ub(6+i))
         i=6
         x(lb(1+i):ub(1+i)) =  real(dyx(0:nx)) - x(lb(1+i):ub(1+i))
         x(lb(2+i):ub(2+i)) = aimag(dyx(0:nx)) - x(lb(2+i):ub(2+i))
         x(lb(3+i):ub(3+i)) =  real(dyy(0:nx)) - x(lb(3+i):ub(3+i))
         x(lb(4+i):ub(4+i)) = aimag(dyy(0:nx)) - x(lb(4+i):ub(4+i))
         x(lb(5+i):ub(5+i)) =  real(dyz(0:nx)) - x(lb(5+i):ub(5+i))
         x(lb(6+i):ub(6+i)) = aimag(dyz(0:nx)) - x(lb(6+i):ub(6+i))
         i=12
         x(lb(1+i):ub(1+i)) =  real(dzx(0:nx)) - x(lb(1+i):ub(1+i))
         x(lb(2+i):ub(2+i)) = aimag(dzx(0:nx)) - x(lb(2+i):ub(2+i))
         x(lb(3+i):ub(3+i)) =  real(dzy(0:nx)) - x(lb(3+i):ub(3+i))
         x(lb(4+i):ub(4+i)) = aimag(dzy(0:nx)) - x(lb(4+i):ub(4+i))
         x(lb(5+i):ub(5+i)) =  real(dzz(0:nx)) - x(lb(5+i):ub(5+i))
         x(lb(6+i):ub(6+i)) = aimag(dzz(0:nx)) - x(lb(6+i):ub(6+i))
         i=18
         x(lb(1+i):ub(1+i)) =  real(vx(0:nx)) - x(lb(1+i):ub(1+i))
         x(lb(2+i):ub(2+i)) =  real(vy(0:nx)) - x(lb(2+i):ub(2+i))
         x(lb(3+i):ub(3+i)) =  real(vz(0:nx)) - x(lb(3+i):ub(3+i))
      end if

      if(chr == 3) then
         i=0
         dxx(0:nx) = cmplx( x(lb(1+i):ub(1+i)) , x(lb(2+i):ub(2+i)))
         dxy(0:nx) = cmplx( x(lb(3+i):ub(3+i)) , x(lb(4+i):ub(4+i)))
         dxz(0:nx) = cmplx( x(lb(5+i):ub(5+i)) , x(lb(6+i):ub(6+i)))
         i=6
         dyx(0:nx) = cmplx( x(lb(1+i):ub(1+i)) , x(lb(2+i):ub(2+i)))
         dyy(0:nx) = cmplx( x(lb(3+i):ub(3+i)) , x(lb(4+i):ub(4+i)))
         dyz(0:nx) = cmplx( x(lb(5+i):ub(5+i)) , x(lb(6+i):ub(6+i)))
         i=12
         dzx(0:nx) = cmplx( x(lb(1+i):ub(1+i)) , x(lb(2+i):ub(2+i)))
         dzy(0:nx) = cmplx( x(lb(3+i):ub(3+i)) , x(lb(4+i):ub(4+i)))
         dzz(0:nx) = cmplx( x(lb(5+i):ub(5+i)) , x(lb(6+i):ub(6+i)))
         i=18
         vx(0:nx) = x(lb(1+i):ub(1+i))
         vy(0:nx) = x(lb(2+i):ub(2+i))
         vz(0:nx) = x(lb(3+i):ub(3+i))
      end if
!
!---------------------------------------------------------------------
!
   end subroutine pacvector
!
!---------------------------------------------------------------------
!
end MODULE packSES
!
!---------------------------------------------------------------------
