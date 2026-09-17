c----------------------------------------------------------------------
c     BLOCK FOR USED LINPACK-ROUTINES ZHICO, ZHIFA AND ZHISL
c----------------------------------------------------------------------
c     first some level-1 blas routines:
c----------------------------------------------------------------------
      real function dcabs1(z)
      complex z,zz
      real t(2)
      equivalence (zz,t(1))
      zz = z
      dcabs1 = abs(t(1)) + abs(t(2))
      return
      end

      complex function zdotc(n,zx,incx,zy,incy)
c
c     forms the dot product of a vector.
c     jack dongarra, 3/11/78.
c
      complex zx(1),zy(1),ztemp
      ztemp = (0.0,0.0)
      zdotc = (0.0,0.0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = ztemp + conjg(zx(ix))*zy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      zdotc = ztemp
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        ztemp = ztemp + conjg(zx(i))*zy(i)
   30 continue
      zdotc = ztemp
      return
      end

      subroutine  zdscal(n,da,zx,incx)
c
c     scales a vector by a constant.
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      complex zx(1)
      real da
      integer i,incx,ix,n
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do 10 i = 1,n
        zx(ix) = cmplx(da,0.0)*zx(ix)
        ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        zx(i) = cmplx(da,0.0)*zx(i)
   30 continue
      return
      end

      real function dzasum(n,zx,incx)
c
c     takes the sum of the absolute values.
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      complex zx(1)
      real stemp,dcabs1
      integer i,incx,ix,n
c
      dzasum = 0.0
      stemp = 0.0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do 10 i = 1,n
        stemp = stemp + dcabs1(zx(ix))
        ix = ix + incx
   10 continue
      dzasum = stemp
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        stemp = stemp + dcabs1(zx(i))
   30 continue
      dzasum = stemp
      return
      end

      subroutine zaxpy(n,za,zx,incx,zy,incy)
c
c     constant times a vector plus a vector.
c     jack dongarra, 3/11/78.
c
      complex zx(1),zy(1),za
      real dcabs1
      if(n.le.0)return
      if (dcabs1(za) .eq. 0.0) return
      if (incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        zy(iy) = zy(iy) + za*zx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        zy(i) = zy(i) + za*zx(i)
   30 continue
      return
      end
      subroutine  zswap (n,zx,incx,zy,incy)
c
c     interchanges two vectors.
c     jack dongarra, 3/11/78.
c
      complex zx(1),zy(1),ztemp
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = zx(ix)
        zx(ix) = zy(iy)
        zy(iy) = ztemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
   20 do 30 i = 1,n
        ztemp = zx(i)
        zx(i) = zy(i)
        zy(i) = ztemp
   30 continue
      return
      end
      integer function izamax(n,zx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, 1/15/85.
c     modified 3/93 to return if incx .le. 0.
c
      complex zx(*)
      real smax,dx
      integer i,incx,ix,n
      real dcabs1
c
      izamax = 0
      if( n.lt.1 .or. incx.le.0 )return
      izamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = dcabs1(zx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dcabs1(zx(ix)).le.smax) go to 5
         izamax = i
         smax = dcabs1(zx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = dcabs1(zx(1))
      do 30 i = 2,n
	 dx=dcabs1(zx(i))
         if(dx.le.smax) go to 30
         izamax = i
         smax = dx
   30 continue
      return
      end
c
c-------------------------------------------------------------
c
      subroutine zhico(a,lda,n,kpvt,rcond,z)
c
c-------------------------------------------------------------
c
      integer lda,n,kpvt(1)
      complex a(lda,1),z(1)
      real rcond
c
c     zhico factors a complex hermitian matrix by elimination with
c     symmetric pivoting and estimates the condition of the matrix.
c
c     if  rcond  is not needed, zhifa is slightly faster.
c     to solve  a*x = b , follow zhico by zhisl.
c     to compute  inverse(a)*c , follow zhico by zhisl.
c     to compute  inverse(a) , follow zhico by zhidi.
c     to compute  determinant(a) , follow zhico by zhidi.
c     to compute  inertia(a), follow zhico by zhidi.
c
c     on entry
c
c        a       complex(lda, n)
c                the hermitian matrix to be factored.
c                only the diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     output
c
c        a       a block diagonal matrix and the multipliers which
c                were used to obtain it.
c                the factorization can be written  a = u*d*ctrans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , ctrans(u) is the
c                conjugate transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   real
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       complex(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack zhifa
c     blas zaxpy,zdotc,zdscal,dzasum
c     fortran abs,dmax1,cmplx,conjg,iabs
c
c     internal variables
c
      complex ak,akm1,bk,bkm1,zdotc,denom,ek,t
      real anorm,s,dzasum,ynorm
      integer i,info,j,jm1,k,kp,kps,ks
c
      complex zdum,zdum2,csign1
      real cabs1
      real real,aimag
      complex zdumr,zdumi
      real(zdumr) = zdumr
      aimag(zdumi) = (0.0,-1.0)*zdumi
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
      csign1(zdum,zdum2) = cabs1(zdum)*(zdum2/cabs1(zdum2))
c
c     find norm of a using only upper half
c
      do 30 j = 1, n
         z(j) = cmplx(dzasum(j,a(1,j),1),0.0)
         jm1 = j - 1
         if (jm1 .lt. 1) go to 20
         do 10 i = 1, jm1
            z(i) = cmplx(real(z(i))+cabs1(a(i,j)),0.0)
   10    continue
   20    continue
   30 continue
      anorm = 0.0
      do 40 j = 1, n
         anorm = max(anorm,real(z(j)))
   40 continue
c
c     factor
c
c-------------------------------------------------------------
c
      call zhifa(a,lda,n,kpvt,info)
c
c-------------------------------------------------------------
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  a*y = e .
c     the components of  e  are chosen to cause maximum local
c     growth in the elements of w  where  u*d*w = e .
c     the vectors are frequently rescaled to avoid overflow.
c
c     solve u*d*w = e
c
      ek = (1.0,0.0)
      do 50 j = 1, n
         z(j) = (0.0,0.0)
   50 continue
      k = n
   60 if (k .eq. 0) go to 120
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         kp = iabs(kpvt(k))
         kps = k + 1 - ks
         if (kp .eq. kps) go to 70
            t = z(kps)
            z(kps) = z(kp)
            z(kp) = t
   70    continue
         if (cabs1(z(k)) .ne. 0.0) ek = csign1(ek,z(k))
         z(k) = z(k) + ek
         call zaxpy(k-ks,z(k),a(1,k),1,z(1),1)
         if (ks .eq. 1) go to 80
            if (cabs1(z(k-1)) .ne. 0.0) ek = csign1(ek,z(k-1))
            z(k-1) = z(k-1) + ek
            call zaxpy(k-ks,z(k-1),a(1,k-1),1,z(1),1)
   80    continue
         if (ks .eq. 2) go to 100
            if (cabs1(z(k)) .le. cabs1(a(k,k))) go to 90
               s = cabs1(a(k,k))/cabs1(z(k))
               call zdscal(n,s,z,1)
               ek = cmplx(s,0.0)*ek
   90       continue
            if (cabs1(a(k,k)) .ne. 0.0) z(k) = z(k)/a(k,k)
            if (cabs1(a(k,k)) .eq. 0.0) z(k) = (1.0,0.0)
         go to 110
  100    continue
            ak = a(k,k)/conjg(a(k-1,k))
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = z(k)/conjg(a(k-1,k))
            bkm1 = z(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0
            z(k) = (akm1*bk - bkm1)/denom
            z(k-1) = (ak*bkm1 - bk)/denom
  110    continue
         k = k - ks
      go to 60
  120 continue
      s = 1.0/dzasum(n,z,1)
      call zdscal(n,s,z,1)
c
c     solve ctrans(u)*y = w
c
      k = 1
  130 if (k .gt. n) go to 160
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         if (k .eq. 1) go to 150
            z(k) = z(k) + zdotc(k-1,a(1,k),1,z(1),1)
            if (ks .eq. 2)
     *         z(k+1) = z(k+1) + zdotc(k-1,a(1,k+1),1,z(1),1)
            kp = iabs(kpvt(k))
            if (kp .eq. k) go to 140
               t = z(k)
               z(k) = z(kp)
               z(kp) = t
  140       continue
  150    continue
         k = k + ks
      go to 130
  160 continue
      s = 1.0/dzasum(n,z,1)
      call zdscal(n,s,z,1)
c
      ynorm = 1.0
c
c     solve u*d*v = y
c
      k = n
  170 if (k .eq. 0) go to 230
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         if (k .eq. ks) go to 190
            kp = iabs(kpvt(k))
            kps = k + 1 - ks
            if (kp .eq. kps) go to 180
               t = z(kps)
               z(kps) = z(kp)
               z(kp) = t
  180       continue
            call zaxpy(k-ks,z(k),a(1,k),1,z(1),1)
            if (ks .eq. 2) call zaxpy(k-ks,z(k-1),a(1,k-1),1,z(1),1)
  190    continue
         if (ks .eq. 2) go to 210
            if (cabs1(z(k)) .le. cabs1(a(k,k))) go to 200
               s = cabs1(a(k,k))/cabs1(z(k))
               call zdscal(n,s,z,1)
               ynorm = s*ynorm
  200       continue
            if (cabs1(a(k,k)) .ne. 0.0) z(k) = z(k)/a(k,k)
            if (cabs1(a(k,k)) .eq. 0.0) z(k) = (1.0,0.0)
         go to 220
  210    continue
            ak = a(k,k)/conjg(a(k-1,k))
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = z(k)/conjg(a(k-1,k))
            bkm1 = z(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0
            z(k) = (akm1*bk - bkm1)/denom
            z(k-1) = (ak*bkm1 - bk)/denom
  220    continue
         k = k - ks
      go to 170
  230 continue
      s = 1.0/dzasum(n,z,1)
      call zdscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve ctrans(u)*z = v
c
      k = 1
  240 if (k .gt. n) go to 270
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         if (k .eq. 1) go to 260
            z(k) = z(k) + zdotc(k-1,a(1,k),1,z(1),1)
            if (ks .eq. 2)
     *         z(k+1) = z(k+1) + zdotc(k-1,a(1,k+1),1,z(1),1)
            kp = iabs(kpvt(k))
            if (kp .eq. k) go to 250
               t = z(k)
               z(k) = z(kp)
               z(kp) = t
  250       continue
  260    continue
         k = k + ks
      go to 240
  270 continue
c     make znorm = 1.0
      s = 1.0/dzasum(n,z,1)
      call zdscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0) rcond = ynorm/anorm
      if (anorm .eq. 0.0) rcond = 0.0
      return
      end
c
c-------------------------------------------------------------
c
      subroutine zhisl(a,lda,n,kpvt,b)
c
c-------------------------------------------------------------
c
      integer lda,n,kpvt(1)
      complex a(lda,1),b(1)
c
c     zhisl solves the complex hermitian system
c     a * x = b
c     using the factors computed by zhifa.
c
c     on entry
c
c        a       complex(lda,n)
c                the output from zhifa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        kpvt    integer(n)
c                the pivot vector from zhifa.
c
c        b       complex(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero may occur if  zhico  has set rcond .eq. 0.0
c        or  zhifa  has set info .ne. 0  .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call zhifa(a,lda,n,kpvt,info)
c           if (info .ne. 0) go to ...
c           do 10 j = 1, p
c              call zhisl(a,lda,n,kpvt,c(1,j))
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas zaxpy,zdotc
c     fortran conjg,iabs
c
c     internal variables.
c
      complex ak,akm1,bk,bkm1,zdotc,denom,temp
      integer k,kp
c
c     loop backward applying the transformations and
c     d inverse to b.
c
      k = n
   10 if (k .eq. 0) go to 80
         if (kpvt(k) .lt. 0) go to 40
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 30
               kp = kpvt(k)
               if (kp .eq. k) go to 20
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
   20          continue
c
c              apply the transformation.
c
               call zaxpy(k-1,b(k),a(1,k),1,b(1),1)
   30       continue
c
c           apply d inverse.
c
            b(k) = b(k)/a(k,k)
            k = k - 1
         go to 70
   40    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 2) go to 60
               kp = iabs(kpvt(k))
               if (kp .eq. k - 1) go to 50
c
c                 interchange.
c
                  temp = b(k-1)
                  b(k-1) = b(kp)
                  b(kp) = temp
   50          continue
c
c              apply the transformation.
c
               call zaxpy(k-2,b(k),a(1,k),1,b(1),1)
               call zaxpy(k-2,b(k-1),a(1,k-1),1,b(1),1)
   60       continue
c
c           apply d inverse.
c
            ak = a(k,k)/conjg(a(k-1,k))
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = b(k)/conjg(a(k-1,k))
            bkm1 = b(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0
            b(k) = (akm1*bk - bkm1)/denom
            b(k-1) = (ak*bkm1 - bk)/denom
            k = k - 2
   70    continue
      go to 10
   80 continue
c
c     loop forward applying the transformations.
c
      k = 1
   90 if (k .gt. n) go to 160
         if (kpvt(k) .lt. 0) go to 120
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 110
c
c              apply the transformation.
c
               b(k) = b(k) + zdotc(k-1,a(1,k),1,b(1),1)
               kp = kpvt(k)
               if (kp .eq. k) go to 100
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  100          continue
  110       continue
            k = k + 1
         go to 150
  120    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 1) go to 140
c
c              apply the transformation.
c
               b(k) = b(k) + zdotc(k-1,a(1,k),1,b(1),1)
               b(k+1) = b(k+1) + zdotc(k-1,a(1,k+1),1,b(1),1)
               kp = iabs(kpvt(k))
               if (kp .eq. k) go to 130
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  130          continue
  140       continue
            k = k + 2
  150    continue
      go to 90
  160 continue
      return
      end
c
c-------------------------------------------------------------
c
      subroutine zhifa(a,lda,n,kpvt,info)
c
c-------------------------------------------------------------
c
      integer lda,n,kpvt(1),info
      complex a(lda,1)
c
c     zhifa factors a complex hermitian matrix by elimination
c     with symmetric pivoting.
c
c     to solve  a*x = b , follow zhifa by zhisl.
c     to compute  inverse(a)*c , follow zhifa by zhisl.
c     to compute  determinant(a) , follow zhifa by zhidi.
c     to compute  inertia(a) , follow zhifa by zhidi.
c     to compute  inverse(a) , follow zhifa by zhidi.
c
c     on entry
c
c        a       complex(lda,n)
c                the hermitian matrix to be factored.
c                only the diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       a block diagonal matrix and the multipliers which
c                were used to obtain it.
c                the factorization can be written  a = u*d*ctrans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , ctrans(u) is the
c                conjugate transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if the k-th pivot block is singular. this is
c                     not an error condition for this subroutine,
c                     but it does indicate that zhisl or zhidi may
c                     divide by zero if called.
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas zaxpy,zswap,izamax
c     fortran abs,dmax1,cmplx,conjg,dsqrt
c
c     internal variables
c
      complex ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
      real absakk,alpha,colmax,rowmax
      integer imax,imaxp1,j,jj,jmax,k,km1,km2,kstep,izamax
      logical swap
c
      complex zdum
      real cabs1
      real real,aimag
      complex zdumr,zdumi
      real(zdumr) = zdumr
      aimag(zdumi) = (0.0,-1.0)*zdumi
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
c     initialize
c
c     alpha is used in choosing pivot block size.
      alpha = (1.0 + sqrt(17.0))/8.0
c
      info = 0
c
c     main loop on k, which goes from n to 1.
c
      k = n
   10 continue
c
c        leave the loop if k=0 or k=1.
c
c     ...exit
         if (k .eq. 0) go to 200
         if (k .gt. 1) go to 20
            kpvt(1) = 1
            if (cabs1(a(1,1)) .eq. 0.0) info = 1
c     ......exit
            go to 200
   20    continue
c
c        this section of code determines the kind of
c        elimination to be performed.  when it is completed,
c        kstep will be set to the size of the pivot block, and
c        swap will be set to .true. if an interchange is
c        required.
c
         km1 = k - 1
         absakk = cabs1(a(k,k))
c
c        determine the largest off-diagonal element in
c        column k.
c
         imax = izamax(k-1,a(1,k),1)
         colmax = cabs1(a(imax,k))
         if (absakk .lt. alpha*colmax) go to 30
            kstep = 1
            swap = .false.
         go to 90
   30    continue
c
c           determine the largest off-diagonal element in
c           row imax.
c
            rowmax = 0.0
            imaxp1 = imax + 1
            do 40 j = imaxp1, k
               rowmax = max(rowmax,cabs1(a(imax,j)))
   40       continue
            if (imax .eq. 1) go to 50
               jmax = izamax(imax-1,a(1,imax),1)
               rowmax = max(rowmax,cabs1(a(jmax,imax)))
   50       continue
            if (cabs1(a(imax,imax)) .lt. alpha*rowmax) go to 60
               kstep = 1
               swap = .true.
            go to 80
   60       continue
            if (absakk .lt. alpha*colmax*(colmax/rowmax)) go to 70
               kstep = 1
               swap = .false.
            go to 80
   70       continue
               kstep = 2
               swap = imax .ne. km1
   80       continue
   90    continue
         if (max(absakk,colmax) .ne. 0.0) go to 100
c
c           column k is zero.  set info and iterate the loop.
c
            kpvt(k) = k
            info = k
         go to 190
  100    continue
         if (kstep .eq. 2) go to 140
c
c           1 x 1 pivot block.
c
            if (.not.swap) go to 120
c
c              perform an interchange.
c
               call zswap(imax,a(1,imax),1,a(1,k),1)
               do 110 jj = imax, k
                  j = k + imax - jj
                  t = conjg(a(j,k))
                  a(j,k) = conjg(a(imax,j))
                  a(imax,j) = t
  110          continue
  120       continue
c
c           perform the elimination.
c
            do 130 jj = 1, km1
               j = k - jj
               mulk = -a(j,k)/a(k,k)
               t = conjg(mulk)
               call zaxpy(j,t,a(1,k),1,a(1,j),1)
               a(j,j) = cmplx(real(a(j,j)),0.0)
               a(j,k) = mulk
  130       continue
c
c           set the pivot array.
c
            kpvt(k) = k
            if (swap) kpvt(k) = imax
         go to 190
  140    continue
c
c           2 x 2 pivot block.
c
            if (.not.swap) go to 160
c
c              perform an interchange.
c
               call zswap(imax,a(1,imax),1,a(1,k-1),1)
               do 150 jj = imax, km1
                  j = km1 + imax - jj
                  t = conjg(a(j,k-1))
                  a(j,k-1) = conjg(a(imax,j))
                  a(imax,j) = t
  150          continue
               t = a(k-1,k)
               a(k-1,k) = a(imax,k)
               a(imax,k) = t
  160       continue
c
c           perform the elimination.
c
            km2 = k - 2
            if (km2 .eq. 0) go to 180
               ak = a(k,k)/a(k-1,k)
               akm1 = a(k-1,k-1)/conjg(a(k-1,k))
               denom = 1.0 - ak*akm1
               do 170 jj = 1, km2
                  j = km1 - jj
                  bk = a(j,k)/a(k-1,k)
                  bkm1 = a(j,k-1)/conjg(a(k-1,k))
                  mulk = (akm1*bk - bkm1)/denom
                  mulkm1 = (ak*bkm1 - bk)/denom
                  t = conjg(mulk)
                  call zaxpy(j,t,a(1,k),1,a(1,j),1)
                  t = conjg(mulkm1)
                  call zaxpy(j,t,a(1,k-1),1,a(1,j),1)
                  a(j,k) = mulk
                  a(j,k-1) = mulkm1
                  a(j,j) = cmplx(real(a(j,j)),0.0)
  170          continue
  180       continue
c
c           set the pivot array.
c
            kpvt(k) = 1 - k
            if (swap) kpvt(k) = -imax
            kpvt(k-1) = kpvt(k)
  190    continue
         k = k - kstep
      go to 10
  200 continue
      return
      end
c
c-----------------------------------------------------------------------


