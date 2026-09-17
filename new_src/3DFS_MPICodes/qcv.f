c======================================================================
c
c                                        ! taken from  : 21.10.24
      program qcv                        ! last change : 23.10.24
c                                        !
c
c======================================================================
c
      include 'qcv.dat'
c
c---  Give the input and rig the calculation
c
      integer i,ii,ien,inittyp,ierr,rc
      integer icmp
      common /itercount/icmp
c
c--   MPI initialization
c
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierr )
      write(*,*) myid,nproc
c
c-- Initialise the calculation 
c
      pi=acos(-1.)
      sqrt2=sqrt(2.0) 
      sqrth=1.0/sqrt2

      if(myid.eq.0) call init_calc
      call input_bcast
c
c--- Get the OPs done
c
      icmp=0
      outerr=.true.
      call MPI_BCAST(outerr,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      do i=1,5,1
         call dfunc_bcast_1
         call getnewop
         call citerat_bcast_3(erreps)
         if(erreps.lt.eps) goto 100
      enddo
      outerr=.false.
      call MPI_BCAST(outerr,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      do ii=1,5
         call runciterat(eps,pmax(ii),inmax(ii),itmax)
         if(erreps.lt.eps) goto 100
      enddo
 100  continue
c
c---  By this point the calculation should be done
c
 666  call MPI_FINALIZE(rc)
      end       
c
c
c----------------------------------------------------------------------
c----                                                            ------
c----   The following routines are called from the               ------
c----   Eschrig packet and they take care of producing           ------
c----   a new order parameter.                                   ------
c----                                                            ------
c----------------------------------------------------------------------
c
      subroutine getnewop
c
c----------------------------------------------------------------------
c
      include 'qcv.dat'
      integer ix,iy,ii,it,ip,ien,ierr
      real x,d,c,a,s2,ym,sx,sy,tlnt,elnt
      complex en,sem(4,-mx:mx),tem(12)
      complex cpp,cpo,cpm,cop,coo,com,cmp,cmo,cmm
      complex dlxx(0:nx),dlxy(0:nx),dlxz(0:nx),
     +        dlyx(0:nx),dlyy(0:nx),dlyz(0:nx),
     +        dlzx(0:nx),dlzy(0:nx),dlzz(0:nx),
     +        lvx(0:nx),lvy(0:nx),lvz(0:nx)
c
      integer icmp
      common /itercount/icmp
c
      icmp=icmp+1
      s2=sqrt(2.)
      tlnt=esum !t/log(t)
      elnt=0.0  !esum/log(t)

      do ix=0,nx,1
         dnxx(ix)=0.0
         dnxy(ix)=0.0
         dnxz(ix)=0.0
         dnyx(ix)=0.0
         dnyy(ix)=0.0
         dnyz(ix)=0.0
         dnzx(ix)=0.0
         dnzy(ix)=0.0
         dnzz(ix)=0.0
         nvx(ix)=0.0
         nvy(ix)=0.0
         nvz(ix)=0.0
         dlxx(ix)=0.0
         dlxy(ix)=0.0
         dlxz(ix)=0.0
         dlyx(ix)=0.0
         dlyy(ix)=0.0
         dlyz(ix)=0.0
         dlzx(ix)=0.0
         dlzy(ix)=0.0
         dlzz(ix)=0.0
         lvx(ix)=0.0
         lvy(ix)=0.0
         lvz(ix)=0.0
      enddo

      do iy=myid,nx,nproc

         do ii=1,12,1
            tem(ii)=0.0
         enddo 

         do it=1,tmax,1
            do ip=1,11,1
               if(cyl) then
                  call intord_c(it,ip,iy,sem)
               else
                  call intord_v(it,ip,iy,sem)
               endif
               do ien=1,Ncmax,1
                  en=zp(ien)
                  call makeprops(it,ip,iy,ien,en,sem,tem)
               enddo
            enddo
         enddo

         do ii=1,9
            tem(ii)=tem(ii)*tlnt
         enddo
         dlxx(iy)=tem(1)
         dlxy(iy)=tem(2)
         dlxz(iy)=tem(3)

         dlyx(iy)=tem(4)
         dlyy(iy)=tem(5)
         dlyz(iy)=tem(6)

         dlzx(iy)=tem(7)
         dlzy(iy)=tem(8)
         dlzz(iy)=tem(9)

         lvx(iy)=tem(10)
         lvy(iy)=tem(11)
         lvz(iy)=tem(12)
      enddo
c
c -- send out your result and loggit
c         
      call MPI_REDUCE(dlxx,dnxx,pkgsz,MPICMPLX,
     +                          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(dlxy,dnxy,pkgsz,MPICMPLX,
     +                          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(dlxz,dnxz,pkgsz,MPICMPLX,
     +                          MPI_SUM,0,MPI_COMM_WORLD,ierr)

      call MPI_REDUCE(dlyx,dnyx,pkgsz,MPICMPLX,
     +                          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(dlyy,dnyy,pkgsz,MPICMPLX,
     +                          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(dlyz,dnyz,pkgsz,MPICMPLX,
     +                          MPI_SUM,0,MPI_COMM_WORLD,ierr)

      call MPI_REDUCE(dlzx,dnzx,pkgsz,MPICMPLX,
     +                          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(dlzy,dnzy,pkgsz,MPICMPLX,
     +                          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(dlzz,dnzz,pkgsz,MPICMPLX,
     +                          MPI_SUM,0,MPI_COMM_WORLD,ierr)

      call MPI_REDUCE(lvx,nvx,pkgsz,MPICMPLX,
     +                          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(lvy,nvy,pkgsz,MPICMPLX,
     +                          MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(lvz,nvz,pkgsz,MPICMPLX,
     +                          MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if(myid.eq.0) then
         call newses
 
         open(1,file='op_xyz',status='unknown') 
         open(2,file='op_harm',status='unknown') 
         open(3,file='curr',status='unknown') 
         open(10,file='xaxis',status='unknown')
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

            write(1,1000) x,dxx(ix),dxy(ix),dxz(ix),
     +                      dyx(ix),dyy(ix),dyz(ix),
     +                      dzx(ix),dzy(ix),dzz(ix)
            write(2,1000) x,cpp,cpo,cpm,cop,coo,com,cmp,cmo,cmm
            d=abs(dxx(ix))**2+abs(dxy(ix))**2+abs(dxz(ix))**2
     +       +abs(dyx(ix))**2+abs(dyy(ix))**2+abs(dyz(ix))**2
     +       +abs(dzx(ix))**2+abs(dzy(ix))**2+abs(dzz(ix))**2
            d=sqrt(d/3.0)
            c=sqrt(real(vx(ix))**2+real(vy(ix))**2)
            write(3,1000) x,d,c,real(vx(ix)),real(vy(ix)),real(vz(ix))
            write(10,2000) x,real(vy(ix)),d
         enddo
         close(1) 
         close(2) 
         close(3) 
         close(10)
      endif

 1000 format((1x,f8.3),30(1x,e14.6))
 2000 format(1x,f8.3,30(1x,e14.6))
 1100 format()
      return
      end
c
c---------------------------------------------------------------------
