c-----------------------------------------------------------------------
c
c----                     SUBROUTINE DFUNC1                        -----
c
c-----------------------------------------------------------------------
c
      subroutine dfunc1(clrarray,n)
c
c-----------------------------------------------------------------------
c-----                                                            ------
c-----     calls "iterat2d" to give the dfunc1 for CITERAT        ------
c-----     most of the iterat2d input enters via                  ------
c-----     the common/iiter/                                     ------
c-----         common/diter/                                     ------
c-----         common/self/                                    ------
c-----                                                            ------
c-----------------------------------------------------------------------
c-----                                                            ------
c-----     clrarray  = impurity self energies as long array       ------
c-----                 (input and output)                         ------
c-----     n         = dimension of long array (input)            ------
c-----                                                            ------
c-----------------------------------------------------------------------
c
      implicit real(a-h,o-z)
      include 'iter.dat'       ! includes parameter statements
c
c-----------------------------------------------------------------------
c     
      if(myid.eq.0) call makeself(clrarray)
c
c----------------------------------------------------------------------
c
      call dfunc_bcast_1          ! send out the current self-energy fields
      call getnewop               ! compute new the self-energy fields
c
c----------------------------------------------------------------------
c                                                                     |
c     this is the main iteration subroutine which iterates the        |
c     self consistency problem                                        |
c                                                                     |
c----------------------------------------------------------------------
c
      call makearray(clrarray,nn) 
c
c-----------------------------------------------------------------------
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine citerat(vector,czx,n,print,eps,prog,pmax,mmax,kmax)
c
c-----------------------------------------------------------------------
c
      implicit real (a-h,o-y)
      implicit complex (z)
      implicit integer (i-n)
c
c-----------------------------------------------------------------------
c
c     iterative solution of the n-dimensional non-linear equation system
c
c                    czx = vector(czx,n)
c
c     by searching in m-dimensional subspaces (m.le.mmax).
c     (for instance: mmax=5)
c
c    -if (print) the iteration is protocolled.
c    -convergence for n*max_i(|czx_i - vector_i(czx,n)|/||czx||) .le. eps.
c     see subroutine deps1 for condition of convergence!
c    -if the ratio of two increments ||czx - vector(czx,n)|| is
c     larger than 1/prog (or smaller than prog) then the
c     worst points are omitted.
c     (set for instance: prog=0.1d0 ... 0.005d0)
c    -iteration steps dx = p*f from the interpolated minimum of |f|**2
c     in the actual subspace are used (p self-controlled),
c     first p beeing 0.01, later for all p: p.le.pmax.
c     (for instance: pmax=0.1d0 ... 10d0)
c     if no convergence try smaller pmax; if good convergence
c     larger pmax may even accelerate convergence (there is a
c     'optimal' pmax).
c    -the maximum dimension of subspace is mmax.
c    -the maximum number of iterations is kmax.
c
c-----------------------------------------------------------------------
c
c    -uses LINPACK-subroutines which are given at the end of
c     the program in this file
c
c-----------------------------------------------------------------------
c
c  from A.Moebius and H.Eschrig, Dresden
c  (revised 1995/96 by M.Eschrig, Bayreuth)
c  reference: H. Eschrig, Optimized LCAO Method, Chap. 7.4,
c             Springer-Verlag, Berlin, 1989;
c
c-----------------------------------------------------------------------
c
c
      include 'iter.dat'      ! includes parameter statements
c
      logical print,update
      real    czx(nn0)
c
      external      vector
c
      equivalence 
     *     ( cx(1,0)  ,  c(1) ),
     *     ( cu(1,0)  ,  c(nn0*ldm+1) ),
     *     ( ca(0,0)  ,  c(2*nn0*ldm+1) ),
     *     ( cg(0)    ,  c((2*nn0+ldm)*ldm+1) )
      real
     *       cx( 1:nn0 , 0:mmax0 ),
     *       cu( 1:nn0 , 0:mmax0 ),
     *       ca( 0:mmax0 , 0:mmax0 ),
     *       cg( 0:mmax0 )
c
      integer kpvt(0:mmax0)
      real cn(0:mmax0)
      complex zhelp(0:mmax0),zca(0:mmax0,0:mmax0)
      equivalence ( cn(0),ca(0,mmax0) )
c
      common/terat/ c( ndim )
c
c-----------------------------------------------------------------------
c
      call citerat_bcast_1(n,mmax)
c
      if (n.gt.nn0) then
	write(*,*) 'dimension nn0 in citerat4 should be >',n-1
	return
      endif
      if (mmax.gt.n) then
	write(*,*) 'dimension of subspace (mmax) cannot be larger than'
	write(*,*) 'dimension of the whole space! I assume mmax=',n
	mmax=n
      endif
      if (mmax.gt.mmax0) then
	write(*,*) 'dimension mmax0 in citerat4 should be >',mmax-1
	return
      endif
c
c-----------------------------------------------------------------------
c
      pq=1.0
      p=0.01
      k=0
      m=0
      a=0.0
      usub=0.0
      eps1=1.e29
c
c---------------------------------------------------------------------
c
c     storage x ==> x  and calculation of f  and of f *f , i=1,...,m
c                    m                     m         i  m
c
c-----------------------------------------------------------------------
c
    1 continue

      if(myid.eq.0) k=k+1
      call citerat_bcast_2(k)
      if(k.gt.kmax) goto 666

      do i=1,n
         cx(i,m)=czx(i)
      enddo

      call vector(czx,n)
c
c-----------------------------------------------------------------------
c
c     This is the call for a new set of self-energies
c
c-----------------------------------------------------------------------
c
      update=.true.
      if(myid.eq.0) update=.true.
      if(update) then   
c
c-----------------------------------------------------------------------
c
         do i=1,n
           cu(i,m)=czx(i)-cx(i,m)
         enddo

         u=0.0
         do i=1,n
	    u=u+cu(i,m)**2
         enddo
         cn(m)=u

         call deps1(n,cu(1,m),czx(1),eps1,a)
         if(a.lt.eps) a=1.0
         if(eps1.le.eps) then
c
c----------------------------------------------------------------------
c
c       convergence: restoration of x and return
c
c-----------------------------------------------------------------------
c
	    do i=1,n
	       czx(i)=cx(i,m)
	    enddo
            if(myid.eq.0) write(*,101) k,m,sqrt(u/a),eps1,p
            goto 10
         endif
c
c-----------------------------------------------------------------------
c
c     exclusion of useless points (note that m may be changed!)
c
c-----------------------------------------------------------------------
c
         ju=0
         if(m.eq.mmax) ju=1
         if(mmax.eq.0) ju=0
         if(m.gt.1) then
            duu=cn(0)
            if(usub/u.lt.prog.or.duu/u.lt.1.0) ju=m/2
         endif
         if(ju.ne.0) call excluse(ju,n,m)
         pq=min(1.0,sqrt(usub/u))
         if(myid.eq.0) write(*,101) k,m+1,sqrt(u/a),eps1,p
c
c     -----------------------------------------------------------
c      case m=0
c     -----------------------------------------------------------
c
         if(m.eq.0) then
            do i=1,n
	       czx(i)=cx(i,0)+p*cu(i,0)
            enddo
            usub=u
            m=1
            goto 10
         endif
c
c     -----------------------------------------------------------
c      case m=1
c     -----------------------------------------------------------
c
         if(m.eq.1) then
            ca(1,0)=0.0
            do i=1,n
	       ca(1,0)=ca(1,0)+cu(i,1)*cu(i,0)
            enddo
            call rearrange(n,1,u)
            cccc=cn(1)+cn(0)-ca(1,0)-ca(1,0)
            u1= (max(cccc,0.0))/u
            ca(0,1)=(u-ca(1,0))/u
            uu=sqrt( u1 )
            if(uu.lt.1e-15) then
	       call stoch(czx,n,m,p)
   	       usub=1e29
	       m=m+1
	       goto 10
            else
               cg(0)=ca(0,1)/u1
            endif
         else
c
c     -----------------------------------------------------------
c      case m>1
c     -----------------------------------------------------------
c
            m1=m-1
            do j=0,m1
	       ca(m,j)=0.0
	       do i=1,n
	          ca(m,j)=ca(m,j)+cu(i,m)*cu(i,j) 
  	       enddo
            enddo
            call rearrange(n,m,u)
c
c     -----------------------------------------------------------
c
c      Build the matrix (f - f )(f - f ) and the vector -f (f - f ).
c                         i   m   k   m                   m  k   m 
c     -----------------------------------------------------------
c
            do i=0,m1
               cccc=(cn(m)+cn(i)-ca(m,i)-ca(m,i))
               u1=max(cccc,0.0)/u
               ca(i,i)=u1
               zca(i,i)=cmplx(ca(i,i),0.0)
               ca(i,m)=(u-ca(m,i))/u 
               zca(i,m)=cmplx(ca(i,m),0.0)
               i1=i+1
               do j=i1,m1
                  ca(i,j)=(ca(j,i)-ca(m,i)+u-ca(m,j))/u
                  zca(i,j)=cmplx(ca(i,j),0.0)
               enddo
            enddo
c  
c      ----------------------------------------------------------
c                                                    m-1
c      finding the minimum of |f|**2 for f(g) = f  + sum g *(f -f ).
c                                                m   i=0  i   i  m
c      ----------------------------------------------------------
c      this is given by the eq.system (solve for g)
c        m-1             *  *            *  *
c        sum g *(f -f )(f -f )  =  - f (f -f )
c        i=0  i   i  m   k  m         m  k  m
c      ----------------------------------------------------------
c
            call zhico(zca(0,0),ldm,m,kpvt(0),dcond,zhelp(0)) ! condition
c
c     ------------------------------------------------
c      stochastic step if condition of matrix too bad
c     ------------------------------------------------
c
            if(dcond+1.0.eq.1.0) then
  	       call stoch(czx,n,m,p)
    	       usub=1e29
	       m=m+1
	       goto 10
            else
	       m1=m-1
               call zhisl(zca(0,0),ldm,m,kpvt(0),zca(0,m)) ! solves system
	       do i=0,m1
	          cg(i)= real(zca(i,m))   ! solution of eqn.system
 	       enddo
            endif
         endif
c
c      ----------------------------------------------------------
c                      m-1
c     setting x = x  + sum g *(x -x ) + p *f   ,
c                  m   i=0  i   i  m     m  min
c
c                 m-1
c     where p =  |sum g *(x -x )| / |f |
c            m    i=0  i   i  m       m
c
c                      m-1
c     and   f   = f  + sum g *(f -f ) .
c            min   m   i=0  i   i  m 
c      ----------------------------------------------------------
c
         m1=m-1
         do j=1,n
            czx(j)=0.0
            do i=0,m1
	       czx(j)=czx(j)+(cx(j,i)-cx(j,m))*cg(i)
	    enddo
         enddo
c
c     -----------------------------------------------
c     finding the p_m
c     -----------------------------------------------
c
         call p1(n,czx(1),u,pq,pmax,p)
c
c     -----------------------------------------------
c     calculate further the optimal czx and cx (here=f)
c     -----------------------------------------------
c
         do j=1,n
            czx(j)=czx(j)+cx(j,m)
         enddo
         m1=m-1
         mp1=m+1
         do j=1,n
            cx(j,mp1)=cu(j,m)
            do i=0,m1
	       cx(j,mp1)=cx(j,mp1)+(cu(j,i)-cu(j,m))*cg(i)
            enddo
         enddo
         usub=0.0
         do j=1,n
	    usub=usub+cx(j,mp1)**2.
         enddo
c
c     -----------------------------------------------
c     do the mixing
c     -----------------------------------------------
c
         do j=1,n
            czx(j)=czx(j)+p*cx(j,mp1)
         enddo
c
c     -----------------------------------------------------------
c     exclusion of useless points (note that m may be changed!)
c     -----------------------------------------------------------
c
         if(usub/u.lt.prog) then
            ju=m/2
            call excluse(ju,n,m)
         endif

         m=m+1
c
c-- End of the updating routine
c
      endif
c
c-- Send the curreent error to all processors
c
 10   call citerat_bcast_3(eps1)
      erreps=eps1
      if(eps1.gt.eps) goto 1
 666  iterp=p
      continue
c
      return
c-----------------------------------------------------------------------
c
 101  format(1x,i3,1x,i3,1x,e14.5,1x,e14.5,1x,f8.4) 
 102  format(1x,i3,1x,i3,1x,i5,1x,e14.5,1x,e14.5,1x,f8.4) 
c
c-----------------------------------------------------------------------
c
      end
c
c
c-----------------------------------------------------------------------
c
      subroutine excluse(ju,n,m)
c
c-----------------------------------------------------------------------
c
      implicit real (a-h,o-z)
      integer ju,m,n
c
c     exclusion of useless points
c
      include 'iter.dat'
c
      equivalence 
     *     ( cx(1,0)  ,  c(1) ),
     *     ( cu(1,0)  ,  c(nn0*ldm+1) ),
     *     ( ca(0,0)  ,  c(2*nn0*ldm+1) )
      real 
     *       cx( 1:nn0 , 0:mmax0 ),
     *       cu( 1:nn0 , 0:mmax0 ),
     *       ca( 0:mmax0 , 0:mmax0 )
      common/terat/ c( ndim )
c
      do j=ju,m
	do i=1,n
	  cx(i,j-ju)=cx(i,j)
	  cu(i,j-ju)=cu(i,j)
	enddo
      enddo
      do i=ju,m
        do j=ju,m
          ca(i-ju,j-ju)=ca(i,j)
        enddo
	ca(i-ju,mmax0)=ca(i,mmax0)
      enddo
      m=m-ju

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine rearrange(n,m,u)
c
c-----------------------------------------------------------------------
c
      implicit real (a-h,o-z)
      integer n,m
c
c     rearranging of x , f , f *f  with descending order of f **2;
c                     i   i   i  j                           i
      include 'iter.dat'
c
      equivalence 
     *     ( cx(1,0)  ,  c(1) ),
     *     ( cu(1,0)  ,  c(nn0*ldm+1) ),
     *     ( ca(0,0)  ,  c(2*nn0*ldm+1) )
      real 
     *       cx( 1:nn0 , 0:mmax0 ),
     *       cu( 1:nn0 , 0:mmax0 ),
     *       ca( 0:mmax0 , 0:mmax0 )
      common/terat/ c( ndim )
c
      real   cn(0:mmax0)
      equivalence ( cn(0),ca(0,mmax0) )
c
      m1=m-1
      do 31 j=m1,0,-1
      j1=j+1
      uu= cn(j)
      ux= cn(j1)
      if(ux.le.uu) return
      do i=1,n
	cvar    =cx(i,j1)
	cx(i,j1)=cx(i,j)
	cx(i,j) =cvar
	cvar    =cu(i,j1)
	cu(i,j1)=cu(i,j)
	cu(i,j) =cvar
      enddo
      do i=j,m
	cvar    =ca(i,j1)
	ca(i,j1)=ca(i,j)
	ca(i,j) =cvar
      enddo
      do i=0,j1
         cvar=ca(j1,i)
         ca(j1,i)=ca(j,i)
         ca(j,i)=cvar
      enddo
      cn(j)=ux
      cn(j1)=uu
      ca(j1,j)=ca(j,j1)
   31 continue
      u= cn(m)
c
      return
      end

c
c-----------------------------------------------------------------------
c
      subroutine deps1(n,cum,czx,eps1,a)
c
c-----------------------------------------------------------------------
c
      implicit real (a-h,o-z)
      integer n
      real cum(n),czx(n)
c
c     calculating norm of x-vector ||x||:
c
      a=0.0
      do i=1,n
	a=a+czx(i)**2
      enddo
c
c     calculating maximal relative error:
c                                    /
c     eps1=n * max ( |G_i(x)-x_i| ) /  ||x||
c               i                  /
      eps1=0.0
      do j=1,n
         eps=cum(j)**2
         eps1=max(eps,eps1)
      enddo
      eps1=sqrt(eps1)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine p1(n,czx,u,pq,pmax,p)
c
c-----------------------------------------------------------------------
c
      implicit real (a-h,o-z)
      integer n
      real czx(n)
c
c     updating mixing strength p
c
      q=p
      dvar=0.0
      do i=1,n
        dvar=dvar+czx(i)**2
      enddo
      p=sqrt(dvar/u)
      p=pq*(p+q)*0.5
      p=min(pmax,p)
      p=max(p,0.01)
      return
      end

c
c-----------------------------------------------------------------------
c
      subroutine stoch(czx,n,m,p)
c
c-----------------------------------------------------------------------
c
      implicit real (a-h,o-z)
c
      parameter  (pi2=6.283185307179586)
      include 'iter.dat'
c
      real     czx(n)
c
      equivalence 
     *     ( cx(1,0)  ,  c(1) ),
     *     ( cu(1,0)  ,  c(nn0*ldm+1) )
      real 
     *       cx( 1:nn0 , 0:mmax0 ),
     *       cu( 1:nn0 , 0:mmax0 )
      common/terat/ c( ndim )
      save dzuf
      data dzuf/pi2/
c
c     creation of a stochastic step, if no minimum of u**2 is found
c
      write(99,*) 'stochastic step.'
      dzuf=(dzuf-nint(dzuf))*10+pi2*1e-10
      argmax1=pi2*1000*dzuf/float(n)
      do i=1,n
	di=float(i)
	p=sin(argmax1*di)
        czx(i)=cx(i,m)+p*cu(i,m)
      enddo
      p=1.0
c
      return
      end
c
c-----------------------------------------------------------------------
c     END OF CITERAT SUBROUTINES
c-----------------------------------------------------------------------
