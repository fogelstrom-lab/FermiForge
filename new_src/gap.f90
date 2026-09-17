!---------------------------------------------------------------------
!
   program gap
!
!---------------------------------------------------------------------
!
!     Version edited 20241231
!
!---------------------------------------------------------------------
!
      use bulkop

      implicit none
      real :: t, gapM

      gapM = 0.1
      print *, ' Give temperature (temp), representation (irep)'
      read *, t
      read *, irep 

      print *, ' This is the result '
      do while (t < 1.0) 
         call bulkgap(t,gapM)
         print 1000,   t, gapM
         t = t +0.01
      end do

  1000 format(2(1x,f8.6))
!
end program gap
!
!---------------------------------------------------------------------

