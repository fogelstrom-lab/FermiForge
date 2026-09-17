       program egrid
       integer nmax
       parameter (nmax=59) 
       real pi,x,y,Rx,dx,xgrid(0:nmax)
       pi=acos(-1.)
       Rx=8.0
       dx=0.42*pi/float(nmax)
       do i=0,nmax,1
          xgrid(i)=Rx*tan(dx*i)
          write(1,*) xgrid(i),i
       enddo

       do i=0,100,1
          x=0.5*i
          y =atan(x/Rx)/dx
          j=int(y)
          if(j.lt.nmax) write(*,*) i,j,xgrid(j),x,xgrid(j+1)
       enddo

       end

