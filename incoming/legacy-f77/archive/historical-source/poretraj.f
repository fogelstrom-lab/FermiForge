C     Up to the subroutine trajectories, the beginning
C     is a test program for the subroutine.  I have used
C     it indebugging.

      common/data/r,rsq
      common/bb/step,nsn
      common/reff/fxo,fx(0:1000),fyo,fy(0:1000),fzo,fz(0:1000)
      common/rerf/rxo,rx(0:1000),ryo,ry(0:1000),rzo,rz(0:1000)
      common/fpp/fpxo,fpx(0:1000),fpyo,fpy(0:1000),fpzo,fpz(0:1000)
      common/rpp/rpxo,rpx(0:1000),rpyo,rpy(0:1000),rpzo,rpz(0:1000)
CC      common/cfcc/nfc,fcx(0:1000),fcy(0:1000),fcz(0:1000)
CC      common/crcc/nrc,rcx(0:1000),rcy(0:1000),rcz(0:1000)

      write(*,*) ' give the radius '
      read(*,*) r
      write(*,*) ' give the step '
      read(*,*) step
      write(*,*) ' give the number of steps '
      read(*,*)  nsn
      rsq=r*r

      write(6,*) 'initial position'
      read(5,*) xo,yo,zo
      write(6,*) 'direction cosines'
      read(5,*) px,py,pz

      call trajectories(xo,yo,zo,px,py,pz)
 
      open(10,status='unknown',file='traj')
      write(10,200)xo,yo,zo
      write(10,200)px,py,pz
      write(10,*)
      do ii=1,nsn
       write(10,100) rx(ii),ry(ii),rz(ii),fx(ii),fy(ii),fz(ii)
      enddo
      write(10,*)
      do ii=1,nsn
       write(10,100) rpx(ii),rpy(ii),rpz(ii),
     +               -fpx(ii),-fpy(ii),-fpz(ii)
      enddo
100   format(3(x,f6.3),2x,3(x,f6.3))
200   format(3(x,f6.3))
      tp=atan2(ry(nsn),rx(nsn))
      tn=atan2(fy(nsn),fx(nsn))
      fp=atan2((-rpx(nsn)*sin(tp)+rpy(nsn)*cos(tp)),
     +         (rpx(nsn)*cos(tp)+rpy(nsn)*sin(tp)))
      fn=atan2((fpx(nsn)*sin(tn)-fpy(nsn)*cos(tn)),
     +         (-fpx(nsn)*cos(tn)-fpy(nsn)*sin(tn)))
      write(*,300) fp,tp,fn,tn 
      write(10,300) fp,tp,fn,tn 
300   format(2(x,f6.3),2x,2(x,f6.3))
      close(10)
      end

      subroutine trajectories(xo,yo,zo,px,py,pz)

C     computes the integration points on trajectories specified by
C     a coordinate point xo,yo,zo, the center, and the direction 
C     cosines px,py,pz. The desired number of integration points along 
C     the trajectory, nsn, and the step must be supplied in the common 
C     block common/bb/step,nsn  and the radius of the hole r, the
C     square of the radius rsq 
C     in the block common/data/r,rsq
C     The subroutine exits with the coordinates of the integration
C     points as well as the directions of the trajectories
C     along the half trajectories at the same points.  The
C     integration points will be in the arrays 
C     common/reff/fxo, fx(0:1000),fyo,fy(0:1000),fzo,fz(0:1000) forward
C     and common/rerf/rxo,rx(0:1000),ryo,ry(0:1000),rzo,rz(0:1000) backward
C     The direction of the trajectory at these points will be in
C     common/fpp/fpxo,fpx(0:1000),fpyo,fpy(0:1000),fpzo,fpz(0:1000) forward
C     and common/rpp/rpxo,rpx(0:1000),rpyo,rpy(0:1000),rpzo,rpz(0:1000)
C     backward.  The points
C     are indexed from the center outward so that the center point
C     figures in both fx(0),fy(0),fz(0) and in rx(0),ry,(0),rz(o).
C     which are the positions fxo,fyo,fzo,rxo,ryo,rzo on the common
C     list.  The same applies to the directions arrays. 
C     The rest of the common blocks are discussion between
C     this subroutine and the subroutine to it, zigzag.
C     The coordinate system has its origin somewhere in the hole,
C     and the z axis is along the axis of the cylinder.
C     Instructions made into comments with two C's refer
C     to lines writing down the collision points with
C     the walls.  The result will come in the common blocks
C     common/cfcc/nfc,fcx(0:1000),fcy(0:1000),fcz(0:1000) and
C     common/crcc/nrc,rcx(0:1000),fcy(0:1000),rcz(0:1000) where
C     the quantities nfc and nrc give the number of collisions
C     along each half trajectory.

      real xo,yo,zo,px,py,pz
      common/data/r,rsq
      common/inp/x1,y1,z1,x2,y2,z2,slong,q21
      common/bb/step,nsn
      common/resf/sxo,sx(0:1000),syo,sy(0:1000),szo,sz(0:1000)
      common/rerf/rxo,rx(0:1000),ryo,ry(0:1000),rzo,rz(0:1000)
      common/reff/fxo,fx(0:1000),fyo,fy(0:1000),fzo,fz(0:1000)
      common/spp/spxo,spx(0:1000),spyo,spy(0:1000),spzo,spz(0:1000)
      common/fpp/fpxo,fpx(0:1000),fpyo,fpy(0:1000),fpzo,fpz(0:1000)
      common/rpp/rpxo,rpx(0:1000),rpyo,rpy(0:1000),rpzo,rpz(0:1000)
CC      common/cfcc/nfc,fcx(0:1000),fcy(0:1000),fcz(0:1000)
CC      common/crcc/nrc,rcx(0:1000),rcy(0:1000),rcz(0:1000)
CC      common/ccp/nisec,scx(0:1000),scy(0:1000),scz(0:1000)

CC      nfc=0
CC      nrc=0

      eps=1.e-10
      ragn=1.e+10
      ssx=px*step
      ssy=py*step
      ssz=pz*step
C     a multiplication of the step with the direction cosines
C     done once and for all.

      sls=nsn*step
C     length of the half trajectory

      nspo=0
      fx(nspo)=xo
      fy(nspo)=yo
      fz(nspo)=zo
      rx(nspo)=xo
      ry(nspo)=yo
      rz(nspo)=zo
      fpx(nspo)=px
      fpy(nspo)=py
      fpz(nspo)=pz
      rpx(nspo)=px
      rpy(nspo)=py
      rpz(nspo)=pz
C     coordinates of the center in all the output files


C     prepare the stuff needed for the 
C     trajectories in the hole
      sth=sqrt(1.-pz*pz)
      if(abs(sth).lt.eps) then
        coph=1.
        siph=0.
        sth=eps
      else
        coph=px/sth
        siph=py/sth
      end if
C     sin(theta), cos(phi), sin(phi) from the direction
C     cosines.  Juglings take place in order to
C     avoid singularities. 

      qout=-xo*coph-yo*siph
      qin=(xo*siph-yo*coph)**2
      qroot=sqrt(rsq-qin)
      q2=qout+qroot
      q1=qout-qroot
      q21=2.*qroot
C     q21 is the distance between x2,y2 and x1,y1 with sign,
C     q2 is the trajectory coordinate ( with sign with respect to
C     the k-vector)  of the index 2 point and q1 of the index 1 point

      x1=xo+q1*coph
      x2=xo+q2*coph

      y1=yo+q1*siph
      y2=yo+q2*siph

      z1=zo+q1*pz/sth
      z2=zo+q2*pz/sth
C     x1,y1,z1 and x2,y2,z2 are the coordinates of the
C     points where the trajectory first meets the walls
C     on both sides of the center

      if(z2.gt.z1) then
         aa=q1
         q1=q2
         q2=aa
         aa=x1
         x1=x2
         x2=aa
         aa=y1
         y1=y2
         y2=aa
         aa=z1
         z1=z2
         z2=aa
      else
         continue
      end if
C     the points x1,y1,z2 and x2,y2,z2 are arranged so that
C     z2 is smaller than z1.  This is an unnecessary step
C     but it stems from the parent program and is not too
C     time consuming.

      slong=abs(q2/sth)
      q21=abs(q21/sth)
C     the definitions of the poins x1,y1,z1 and x2,y2,z2 for the
C     subroutine zigzag.  slong 
C     is the distance from the center to the above defined 
C     point 2, and q21 is the distance between points 2 and 1.   
C     z2 is taken as the point that lies deeper (deeper means
C     further in the negative direction) in the hole.
C     We will have to exchange the order of point 1 and point 
C     2 for computing  the upward zigzag.


 220  if(z1.eq.z2) then
         if(px.eq.0.) go to 132
         aaa=(x2-xo)/px
         go to 133
 132     aaa=(y2-yo)/py
 133     if(aaa.gt.0.) then
            pz=-eps
         else
            pz=eps
         end if
      else
         continue
      end if
C     this jugling is for the case that the trajectory is
C     strictly perpendicular to the cylinder axis within
C     the cylinder.  Then it is investigated whether the
C     direction cosine vector points to z2 or not.  If yes,
C     we go to the case two of the first choice below etc.
C     by assigning a negative value to pz.

      asit=slong
C     distance from zo to the wall

      if(pz.gt.0) then
C     we know that z2 lies below z1.  If the trajectory points
C     upward, the tail will point down in the direction of x2
   
 57      if(step.lt.asit) then
            nspo=nspo+1
            if(nspo.gt.nsn) go to 33
            rx(nspo)=rx(nspo-1)-ssx
            ry(nspo)=ry(nspo-1)-ssy
            rz(nspo)=rz(nspo-1)-ssz
            rpx(nspo)=px
            rpy(nspo)=py
            rpz(nspo)=pz
C           the integration points are written down

            asit=asit-step
            go to 57
C           even here there can be several steps before hitting
C           the wall

         else
            nspl=nspo+1
            call zigzag
            do 222 ii=nspl,nsn
            rx(ii)=sx(ii)
            ry(ii)=sy(ii)
            rz(ii)=sz(ii)
            rpx(ii)=-spx(ii)
            rpy(ii)=-spy(ii)
C           the integration points are in the array s when
C           they come out of the subroutine zigzag and they
C           must be written in the proper block according
C           to whether they belong to the forward part or
C           the tail.

            rpz(ii)=pz
 222        continue


CC            nrc=nisec
CC            do 887 ii=1,nisec
CC            rcx(ii)=scx(ii)
CC            rcy(ii)=scy(ii)
CC            rcz(ii)=scz(ii)
CC 887        continue

         end if
      else
C        the downward trajectory, the f-half ray is done
C        downward

 58      if(step.lt.asit) then
            nspo=nspo+1
            if(nspo.gt.nsn) go to 33
            fx(nspo)=fx(nspo-1)+ssx 
            fy(nspo)=fy(nspo-1)+ssy
            fz(nspo)=fz(nspo-1)+ssz
            fpx(nspo)=px
            fpy(nspo)=py
            fpz(nspo)=pz
            asit=asit-step
            go to 58
C           the steps before hitting the wall
         else
            nspl=nspo+1
            call zigzag
            do 62 ii=nspl,nsn
            fx(ii)=sx(ii)
            fy(ii)=sy(ii)
            fz(ii)=sz(ii)
            fpx(ii)=spx(ii)
            fpy(ii)=spy(ii)
            fpz(ii)=pz
 62         continue
 
CC            nfc=nisec
CC            do 885 ii=1,nisec
CC            fcx(ii)=scx(ii)
CC            fcy(ii)=scy(ii)
CC            fcz(ii)=scz(ii)
CC 885        continue

         end if 
      end if

C     remain the f-half of the upward trajectory and the r-half
C     of the downward trajectory.  The order of the points 1 and
C     2, the intersections of the trajectory and the cylinder,
C     have to be exchanged:

 33   aa=q2
      q2=q1
      q1=aa
      aa=x1
      x1=x2
      x2=aa
      aa=y1
      y1=y2
      y2=aa
      aa=z1
      z1=z2
      z2=aa
      slong=abs(q2/sth)
 
      asit=slong
      nspo=0
C     the distance to the wall and the s-index counter have
C     to be reset

      if(pz.gt.0.) then
 59      if(step.lt.asit) then
            nspo=nspo+1
            if(nspo.gt.nsn) return
            
            if(nspo.eq.1) write(6,*) fx(nspo-1)

            fx(nspo)=fx(nspo-1)+ssx 
            fy(nspo)=fy(nspo-1)+ssy
            fz(nspo)=fz(nspo-1)+ssz
            fpx(nspo)=px
            fpy(nspo)=py
            fpz(nspo)=pz
            asit=asit-step
            go to 59
C           the steps before hitting the wall
         else
            call zigzag
            nspl=nspo+1
            do 65 ii=nspl,nsn
            fx(ii)=sx(ii)
            fy(ii)=sy(ii)
            fz(ii)=sz(ii)
            fpx(ii)=spx(ii)
            fpy(ii)=spy(ii)
            fpz(ii)=pz
 65         continue

CC            nfc=nisec
CC            do 884 ii=1,nisec
CC            fcx(ii)=scx(ii)
CC            fcy(ii)=scy(ii)
CC            fcz(ii)=scz(ii)
CC 884        continue

         end if
      else
 56      if(step.lt.asit) then
            nspo=nspo+1
            if(nspo.gt.nsn) return
            rx(nspo)=rx(nspo-1)-ssx
            ry(nspo)=ry(nspo-1)-ssy
            rz(nspo)=rz(nspo-1)-ssz
            rpx(nspo)=px
            rpy(nspo)=py
            rpz(nspo)=pz
            asit=asit-step
            go to 56
C           even here there can be several steps before hitting
C           the wall
         else
            call zigzag
            nspl=nspo+1
            do 66 ii=nspl,nsn
            rx(ii)=sx(ii)
            ry(ii)=sy(ii)
            rz(ii)=sz(ii)
            rpx(ii)=-spx(ii)
            rpy(ii)=-spy(ii)
            rpz(ii)=pz
 66         continue

CC            nrc=nisec
CC            do 883 ii=1,nisec
CC            rcx(ii)=scx(ii)
CC            rcy(ii)=scy(ii)
CC            rcz(ii)=scz(ii)
CC 883        continue
         end if
      end if
C     the forward end of the upward trajectory and the
C     tail end of the downward trajectory have been completed 

      return
      end

      subroutine zigzag

C     this subroutine computes the zigzag path of a trajectory 
C     It needs the two intersection points,
C     indexed 1 and 2, of the cylinder and the    
C     trajectory.  It also needs the distance q21 along the trajectory
C     between these points  and the length of the trajectory up to the 
C     point x2,y2,z2.  The subroutine proceeds to look for the zigzag
C     path in the direction of x2,y2,z2 (as opposed to the direction of
C     x1,y1,z1. It exits with the coordinates of the integration points
C     along the trajectory at intervals specified
C     by the variable step plus the direction cosines at these 
C     points in the direction from point 1 toward point 2.
C     The former are in the array sx(0:1000), sy(0:1000), sz(0:1000) 
C     and the latter in the array spx(0:1000), spy(0:1000),spz(0:1000). 
C     The subroutine obviously also needs the radius of the hole =r
C     and the square of the same =rsq in particular.
C     The number of integration points required on the
C     half trajectory , nsn, and the step are clearly also needed.
C     These are given via the common blocks.
C     All variables beginning with the letter s refer to the integration
C     path and those beginning with q to the zigzag path of the
C     reflections from the cylinder (they are the same, of course,
C     but have different milestones, as it were).  The
C     intersection points themselves are an exception to the q-s
C     rule.
C     If the collision points with the walls are needed, they
C     are given written down in instructions that are made into
C     comments by two C's.


      common/data/ r,rsq
      common/inp/x1,y1,z1,x2,y2,z2,slong,q21
      common/bb/step,nsn
      common/resf/sxo,sx(0:1000),syo,sy(0:1000),szo,sz(0:1000)
      common/spp/spxo,spx(0:1000),spyo,spy(0:1000),spzo,spz(0:1000)
CC      common/ccp/nisec,scx(0:1000),scy(0:1000),scz(0:1000)

CC      nisec=0
      ragn=1.e+10
      ssa=step/q21
C     for the advancement along the trajectory, the integration
C     step is divided with the normalization of the direction cosines
C     determined by the successive intersection points.

      sinab=(x2*y1-y2*x1)/rsq
      cosab=(y1*y2+x1*x2)/rsq
C     sine and cosine of the angle between two intersections

      qdz=z2-z1
C     difference of the z-coordinate between consequtive intersections

      nspo=int(slong/step)
      spo=nspo*step
      qpo=slong
C     index of the integration point closest to x2,y2,z2, its coordinate
C     on trajectory, and the q-position of x2,y2,z2, i.e.
C     the position along the zigzag trajectory of intersection points

      xn=x2
      yn=y2
      zn=z2
      xcu=x1
      ycu=y1
      zcu=z1
C     the new and current positions of the intersection points
      
      saa=(qpo-spo)/q21
      sx(nspo)=(x1-x2)*saa+x2
      sy(nspo)=(y1-y2)*saa+y2
      sz(nspo)=-qdz*saa+z2
      spx(nspo)=(x2-x1)/q21
      spy(nspo)=(y2-y1)/q21
C     the integration point right before the first reflection
C     recorded in the output file.  The direction cosines come
C     with (2-1) although the position of the point is
C     computed backward from point 2.

      sfolpo=spo+step
C     the following integration point after the present spo
C     Precautions taken in the main subroutine guarantee 
C     that the present zn does not lie outside the passage.

 12   if(sfolpo.lt.qpo) then
C     will a new q-point have to be computed for the next
C     integration point? If not, then:

         nspo=nspo+1
         if(nspo.gt.nsn) return
C        if enough integration points have been calculated along
C        the half trajectory, return with the result in the array
C        sx,sy,sz from index int(slong/step) forward

         sx(nspo)=(xn-xcu)*ssa+sx(nspo-1)
         sy(nspo)=(yn-ycu)*ssa+sy(nspo-1)
         sz(nspo)=qdz*ssa+sz(nspo-1)
         spx(nspo)=(xn-xcu)/q21
         spy(nspo)=(yn-ycu)/q21
         sfolpo=sfolpo+step
C        new point on the integration trajectory.  The z-coordinate
C        is calculated simpler since the z-difference is always the same


         go to 12

      else

 10      xcu=xn
         ycu=yn
         zcu=zn

CC         nisec=nisec+1
CC         scx(nisec)=xn
CC         scy(nisec)=yn
CC         scz(nisec)=zn
C        the collision point with the wall is written down
C        if the double comment sign is missing in the above
C        instructions.

         xn=(xcu*cosab+ycu*sinab)
         yn=(ycu*cosab-xcu*sinab)
         zn=zcu+qdz
C        the new intersection point computed by the sine and cosine
C        addition formulas

         qpo=qpo+q21
C        the new q-position after incorporating one more intersection
         
         if(sfolpo.gt.qpo) go to 10
C        continue computing new q-points till you reach sfolpo
C        with qpo (as long as the trajectory collides with the
C        walls of the pore). The problem one grapples with here 
C        is that it is possible that the q-step will be small.
      end if

      nspo=nspo+1
      if(nspo.gt.nsn) return
C     again the half trajectory is complete

      saa=(-qpo+q21+sfolpo)/q21
 13   sx(nspo)=(xn-xcu)*saa+xcu
      sy(nspo)=(yn-ycu)*saa+ycu
      sz(nspo)=saa*qdz+zcu
      spx(nspo)=(xn-xcu)/q21
      spy(nspo)=(yn-ycu)/q21
      sfolpo=sfolpo+step
      go to 12

      end
c
c--------------------------------------------------------------------------
