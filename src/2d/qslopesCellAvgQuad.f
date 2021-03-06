c
c ---------------------------------------------------------------------
c
       subroutine qslopesCellAvgQuad(qp,qx,qy,qxx,qxy,qyy,mitot,mjtot,
     &                  irr,lstgrd,lwidth,hx,hy,xlow,ylow,mptr,
     &                  nvar,istage) 

      use amr_module
      implicit double precision(a-h,o-z)

      dimension qp(nvar,mitot,mjtot),qx(nvar,mitot,mjtot),
     &          qy(nvar,mitot,mjtot),irr(mitot,mjtot)
      dimension qxx(nvar,mitot,mjtot),qyy(nvar,mitot,mjtot),
     &          qxy(nvar,mitot,mjtot)
      logical nolimiter
      common /order2/ ssw, quad, nolimiter
      include "cuserdt.i"


      dimension a(25,5),at(5,25),c(5,5),b(25,nvar),d(25,nvar)
      dimension rhs(25,nvar), w(5)
      dimension nlist(25,2)
      logical   prflag, quad, enufNbor
      logical all_nbors_exist
      logical IS_EXTERIOR_AND_BAD, IS_EXTERIOR
      data      prflag/.true./

      IS_EXTERIOR(x,y) = (x .lt. 0.d0 .or. x .gt. xprob .or. 
     &                    y .lt. 0.d0 .or. y .gt.  yprob)

      IS_EXTERIOR_AND_BAD(x,y,i,j) = (IS_EXTERIOR(x,y) .and.
     &                               (i .lt. 1+2*(istage-1) .or. 
     &                                j .lt. 1+2*(istage -1).or. 
     &                                i . gt. mitot-2*(istage-1) .or. 
     &                                j . gt. mjtot-2*(istage-1)))


c   ##########
c   #  compute slopes for cut cells using least squares approach
c   #  code for igradChoice  = 3 - fit quadratic using cell averages 
c   ##########

c
c   ssw = 0 for ctu (no slopes), 1 for muscl (second order).
c   qp contains primitive variables, regular cell slopes already set 
c   all slopes initialized to zero in regular slopes subr.
c
c
c     first save everything into unlim (even though limited already
c     on flow cells. will add unlimited cut gradients

      nTerms = 5 ! for quadratic through cell avg
      quad = .false.
      nco = 1

      do 110 ix0 = 1, mitot
      do 110 iy0 = 1, mjtot
         k = irr(ix0,iy0)
         if (k .eq. -1) then  ! set to 0 for easier debugging
            qx(:,ix0,iy0) = 0.d0
            qy(:,ix0,iy0) = 0.d0
            qxx(:,ix0,iy0) = 0.d0
            qxy(:,ix0,iy0) = 0.d0
            qyy(:,ix0,iy0) = 0.d0
            go to 110
         endif

c     
         if (ar(k)/ar(lstgrd) .lt. gradThreshold) then
            qx(:,ix0,iy0) = 0.d0
            qy(:,ix0,iy0) = 0.d0
            qxx(:,ix0,iy0) = 0.d0
            qxy(:,ix0,iy0) = 0.d0
            qyy(:,ix0,iy0) = 0.d0
            go to 110    ! leave 0 gradient:  more stable for teeny cells w/o slopes
         endif
c     
c      # need coords for exterior untrusted cells
         if (k .ne. lstgrd) then 
            x0 = xcirr(k)
            y0 = ycirr(k)
         else
            x0 = xlow + (ix0-.5d0)*hx
            y0 = ylow + (iy0-.5d0)*hy
         endif
         if (IS_EXTERIOR_AND_BAD(x0,y0,ix0,iy0)) then
            ! reset to zero whatever was computed in reg_slopes
            qx(:,ix0,iy0) = 0.d0
            qy(:,ix0,iy0) = 0.d0
            qxx(:,ix0,iy0) = 0.d0
            qxy(:,ix0,iy0) = 0.d0
            qyy(:,ix0,iy0) = 0.d0
            go to 110    
         endif

         if (all_nbors_exist(ix0,iy0,nco,irr,lstgrd,mitot,mjtot)) cycle
         
c      # this cell needs derivatives
         nlist(1,1) = ix0
         nlist(1,2) = iy0
         nst        = 1  
         nend       = 1  
         call addneighbors(irr,nlist,nst,nend,newend,mitot,mjtot,
     .        lstgrd,quad,xlow,ylow,hx,hy)

         if (newend .le. 2) then !could happen since starting at r/c=1
            qx(:,ix0,iy0) = 0.d0
            qy(:,ix0,iy0) = 0.d0
            qxx(:,ix0,iy0) = 0.d0
            qxy(:,ix0,iy0) = 0.d0
            qyy(:,ix0,iy0) = 0.d0
            go to 110    ! leave 0 gradient
         endif

         irow = 0
         do 22 n = 2, newend
            irow = irow + 1
            ixn = nlist(n,1)
            iyn = nlist(n,2)
            kn =  irr(ixn,iyn)
            if (kn .ne. lstgrd) then
               xn = xcirr(kn)
               yn = ycirr(kn)
            else
               xn = xlow + (ixn-.5d0)*hx
               yn = ylow + (iyn-.5d0)*hy
            endif

c     # code preserving integral average
              shiftxx = poly( 8,1,k)
              shiftxy = poly( 9,1,k)
              shiftyy = poly(10,1,k)
              a(irow,1:5) = 0.d0

              ! handle cut cells first
              if (kn .ne. lstgrd) then
                do 15 index = 1, 24
                  xp = points(index,1,kn)
                  yp = points(index,2,kn)
c     
                  a(irow,1) = wt(index,kn)*(xp-x0) + a(irow,1)
                  a(irow,2) = wt(index,kn)*(yp-y0) + a(irow,2)
                  a(irow,3) = .5*wt(index,kn)*(xp-x0)*(xp-x0)
     &                 + a(irow,3)
                  a(irow,4) = wt(index,kn)*(xp-x0)*(yp-y0)
     &                 + a(irow,4)
                  a(irow,5) = .5*wt(index,kn)*(yp-y0)*(yp-y0)
     &                 + a(irow,5)
 15             continue
              else
                do 17 index = 1, 5
                    wt(index,kn) = 1.d0/6.d0
 17             continue

                  wt(3,kn) = 1.d0/3.d0
                  points(1,1,lstgrd) = xn - hx/2.d0
                  points(1,2,lstgrd) = yn
                  points(2,1,lstgrd) = xn
                  points(2,2,lstgrd) = yn + hy/2.d0
                  points(3,1,lstgrd) = xn
                  points(3,2,lstgrd) = yn
                  points(4,1,lstgrd) = xn + hx/2.d0
                  points(4,2,lstgrd) = yn
                  points(5,1,lstgrd) = xn
                  points(5,2,lstgrd) = yn - hy/2.d0
                  do 19 index = 1, 5
                     xp = points(index,1,kn)
                     yp = points(index,2,kn)
                     a(irow,1) = wt(index,kn)*(xp-x0) + a(irow,1)
                     a(irow,2) = wt(index,kn)*(yp-y0) + a(irow,2)
                     a(irow,3) = .5*wt(index,kn)*(xp-x0)*(xp-x0)
     &                    + a(irow,3)
                     a(irow,4) = wt(index,kn)*(xp-x0)*(yp-y0)
     &                    + a(irow,4)
                     a(irow,5) = .5*wt(index,kn)*(yp-y0)*(yp-y0)
     &                    + a(irow,5)
 19               continue
            endif

c            #  shift to fit quadratic terms with mean 0 over cell k
             a(irow,3) = a(irow,3) - shiftxx
             a(irow,4) = a(irow,4) - shiftxy
             a(irow,5) = a(irow,5) - shiftyy
c
            do m = 1, nvar
               b(irow,m) = qp(m,ixn,iyn) - qp(m,ix0,iy0)
            end do
 22      continue
               
               
c
c        ## form normal equations and solve via cholesky
c
          do 30 it = 1, irow
          do 30 jt = 1, nTerms
             at(jt,it) = a(it,jt)
 30       continue

          do 50 it = 1, nTerms
             do 50 jt = 1, nTerms
                c(it,jt) = 0.d0
                do m = 1, nvar
                   d(it,m)  = 0.d0
                end do
                do 45 kt = 1, irow
                   c(it,jt) = c(it,jt) + at(it,kt)*a(kt,jt)
                   do m = 1, nvar
                      d(it,m) = d(it,m) + at(it,kt)*b(kt,m)
                   end do
 45             continue
 50       continue


c         # now solve C*w = d for least squares slopes. use cholesky
c         # put factors back in a
             a(1,1) = dsqrt(c(1,1))
             a(1,2) = c(1,2)/a(1,1)
             a(1,3) = c(1,3)/a(1,1)
             a(1,4) = c(1,4)/a(1,1)
             a(1,5) = c(1,5)/a(1,1)

             a(2,2) = dsqrt(c(2,2)-a(1,2)**2)
             a(2,3) = (c(2,3)-a(1,2)*a(1,3))/a(2,2)
             a(2,4) = (c(2,4)-a(1,2)*a(1,4))/a(2,2)
             a(2,5) = (c(2,5)-a(1,2)*a(1,5))/a(2,2)

             a(3,3) = dsqrt(c(3,3)-a(1,3)**2 - a(2,3)**2)
             a(3,4) = (c(3,4)-a(1,3)*a(1,4)-a(2,3)*a(2,4))/a(3,3)
             a(3,5) = (c(3,5)-a(1,3)*a(1,5)-a(2,3)*a(2,5))/a(3,3)

             a(4,4) = dsqrt(c(4,4)-a(1,4)**2-a(2,4)**2-a(3,4)**2)
             a(4,5) = (c(4,5)-a(1,4)*a(1,5)-a(2,4)*a(2,5)-a(3,4)*a(3,5))
     &            /a(4,4)
             a(5,5) = dsqrt(c(5,5)-a(1,5)**2-
     &                      a(2,5)**2-a(3,5)**2-a(4,5)**2)
c     
             do 60 m = 1, nvar
c     
c     # at*a = c. solve at*b = d, aw = b.  reuse b.
                b(1,m) = d(1,m) / a(1,1)
                b(2,m) = (d(2,m) - a(1,2)*b(1,m)) / a(2,2)
                b(3,m) = (d(3,m) - a(1,3)*b(1,m)-a(2,3)*b(2,m))/a(3,3)
                b(4,m) = (d(4,m) - a(1,4)*b(1,m)-a(2,4)*b(2,m)
     &                    -a(3,4)*b(3,m))/a(4,4)
                b(5,m) = (d(5,m) - a(1,5)*b(1,m)-a(2,5)*b(2,m)
     &                    -a(3,5)*b(3,m)-a(4,5)*b(4,m))/a(5,5)


                w(5) = b(5,m) / a(5,5)
                w(4) = (b(4,m)-a(4,5)*w(5))/a(4,4)
                w(3) = (b(3,m)-a(3,5)*w(5)-a(3,4)*w(4))/a(3,3)
                w(2) = (b(2,m)-a(2,5)*w(5)-a(2,4)*w(4)
     &               -a(2,3)*w(3))/a(2,2)
                w(1) = (b(1,m)-a(1,5)*w(5)-a(1,4)*w(4)
     &               -a(1,3)*w(3)-a(1,2)*w(2))/a(1,1)

                qyy(m,ix0,iy0)  =  w(5)  
                qxy(m,ix0,iy0)  =  w(4)
                qxx(m,ix0,iy0)  =  w(3) 
                qy(m,ix0,iy0)   =  w(2)
                qx(m,ix0,iy0)   =  w(1)
c
 60       continue
c  
 110  continue
c
 120  if (prflag) then
         write(21,*)' qx '
         call prDeriv(qx,irr,mitot,mjtot,nvar,lstgrd)
         write(21,*)' qy '
         call prDeriv(qy,irr,mitot,mjtot,nvar,lstgrd)
         write(21,*)' qxx'
         call prDeriv(qxx,irr,mitot,mjtot,nvar,lstgrd)
         write(21,*)' qxy'
         call prDeriv(qxy,irr,mitot,mjtot,nvar,lstgrd)
         write(21,*)' qyy'
         call prDeriv(qyy,irr,mitot,mjtot,nvar,lstgrd)

      endif
c
 99   return
      end
c
c ----------------------------------------------------------
c
      logical function all_nbors_exist(i,j,nco,irr,lstgrd,mitot,mjtot)

      dimension irr(mitot,mjtot)
     
      all_nbors_exist = .false.  ! initialize

      ! check if index out of bounds
      if (j-nco .lt. 1 .or. j+nco .gt. mjtot .or. 
     &    i-nco .lt. 1 .or. i+nco .gt. mitot) return

      ! check for regular stencil
      do joff = j-nco, j+nco
      do ioff = i-nco, i+nco
        k = irr(ioff,joff)
        if (k .eq. -1) return
        if (k .ne. lstgrd) return
      end do
      end do

      all_nbors_exist = .true. ! if made it to here
      return
      end
