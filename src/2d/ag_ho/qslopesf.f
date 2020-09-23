c
c ---------------------------------------------------------------------
c
       subroutine qslopesf(qp,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,
     &                    hx,hy,xlow,ylow,qxx,qxy,qyy,mptr,nvar, dlimit)
      implicit double precision(a-h,o-z)
      dimension qp(mitot,mjtot,nvar),qx(mitot,mjtot,nvar),
     &          qy(mitot,mjtot,nvar),irr(mitot,mjtot)
      dimension qxx(mitot,mjtot,nvar),qxy(mitot,mjtot,nvar)
      dimension qyy(mitot,mjtot,nvar)
      common /order2/ ssw, quad, nolimiter
      include "cirr.i"
      common   /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      dimension rhsmax(nvar), rhsmin(nvar)
      dimension dumax(nvar),dumin(nvar),phimin(nvar)
      dimension a(25,5),at(5,25),c(5,5),b(25,nvar),d(25,nvar)
      dimension rhs(25,nvar), w(5,nvar)
      dimension nlist(25,2)
      dimension dlimit(mitot,mjtot)
      logical   prflag, quad, enufNbor, outside, almost_outside
      data      prflag/.false./
      logical   nolimiter, ghost_ccg/.false./
      logical ALL_NBORS_EXIST
c
c     outside(x,y) = ((x.lt.0) .or. (y.lt.0.) .or. 
c    .                (x.gt.xprob) .or. (y.gt.yprob))

      ALL_NBORS_EXIST(i,j) = (i .gt. 1 .and. i .lt.mitot .and.
     .                        j .gt. 1 .and. j .lt. mjtot)


c
c   ##########
c   #  compute slopes for cut cells using least squares approach
c   ##########
c   
c   ssw = 0 for ctu (no slopes), 1 for muscl (second order).
c   qp contains primitive variables, regular cell slopes already set 
c   all slopes initialized to zero in regular slopes subr.
c
c
      if (quad) then 
         nterms = 5 
      else 
         nterms = 2 
      endif


c      do 110 ix0 = lwidth-2, mitot-lwidth+3
c      do 110 iy0 = lwidth-2, mjtot-lwidth+3
      do 110 ix0 = 1, mitot
      do 110 iy0 = 1, mjtot
         k = irr(ix0,iy0)
         if (k .eq. -1) go to 110
c     ::: next line is so only cut cells done quadratically,
c     ::: otherwise do their neighbors too
c     if (k .eq. lstgrd) go to 110
         if (k .eq. lstgrd .and. ALL_NBORS_EXIST(ix0,iy0)) then
           if (irr(ix0+1,iy0) .eq. lstgrd .and. 
     .         irr(ix0,iy0+1) .eq. lstgrd .and. 
     .         irr(ix0,iy0-1) .eq. lstgrd .and. 
     .         irr(ix0-1,iy0) .eq. lstgrd) go to 110
         endif
c     
         if (ar(k)/ar(lstgrd) .lt. gradThreshold) then
            go to 110           ! leave 0 gradient:  more stable for teeny cells w/o slopes
         endif
c     
c     ::: use one-sided 2nd order slopes 
c     if (k .eq. lstgrd .and. .not. quad) go to 110
c     
c      # this cell needs a slope
         if (k .ne. lstgrd) then
            x0 = xcirr(k)
            y0 = ycirr(k)
            do 37 kside = 1, 6
               if (poly(kside+2,1,k).eq.-11)then
                  sidex2 = poly(kside,1,k)
                  sidey2 = poly(kside,2,k)
                  sidex1 = poly(kside+1,1,k)
                  sidey1 = poly(kside+1,2,k)
                  go to 39
               endif
 37         continue
 39         rlenb = dsqrt((sidey1-sidey2)**2 + (sidex1-sidex2)**2)
            alf  = (sidey1-sidey2)/rlenb
            beta = (sidex2-sidex1)/rlenb
            bxpt = .5d0*(sidex1+sidex2)
            bypt = .5d0*(sidey1+sidey2)
         else
            x0 = xlow + (ix0-.5d0)*hx
            y0 = ylow + (iy0-.5d0)*hy
         endif
         
c     if (outside(x0,y0)) go to 110   ! to match cart3d, no gradients in ghost cells
         nlist(1,1) = ix0
         nlist(1,2) = iy0
         nst        = 1  
         nend       = 1  
         call addneighbors(irr,nlist,nst,nend,newend,mitot,mjtot,
     .        lstgrd,quad,xlow,ylow,hx,hy)
         if (.not. quad .and. newend .gt. 3) go to 16
c     in the linear case only 2 nbors. could be ill conditioned. add diagonal nbor
c     first have to find diag nbor. this code doesnt test for thin bodies,etc
         ixn = nlist(2,1)
         iyn = nlist(2,2)
         if (nlist(2,1) .eq. ix0) ixn = nlist(3,1)
         if (nlist(2,2) .eq. iy0) iyn = nlist(3,2)
         newend = newend + 1
         nlist(newend,1) = ixn 
         nlist(newend,2) = iyn
         go to 16

c     cell itself in 1st row, 2-6 is 5 more rows. see if sufficient.
         if (quad .and. newend .gt. 6) then ! want more than 5 nbors (pos. 1 is cell itself)
            enufNbor = .TRUE.
            go to 16
         endif

c        ::: not enough neighbors for quadratic fit
c        ::: add neighbors of neighbors    
         nst        = 2  
         nend       = newend  
         call addneighbors(irr,nlist,nst,nend,newend,mitot,mjtot,
     .                     lstgrd,quad,xlow,ylow,hx,hy)
         if (newend .ge. 6)then
             enufNbor = .TRUE.
             go to 16
         endif
c        write(6,*)" problem finding neighbors for quadratic fit"
c        write(6,*) ix0, iy0, " only have ", newend-1
c        write(6,*) "not sufficient "
c        stop
         enufNbor = .FALSE.
c
 16      irow = 0
         if (quad .and. enufNbor) then
            shiftxx = poly( 8,1,k)
            shiftxy = poly( 9,1,k)
            shiftyy = poly(10,1,k)
         endif
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
c     # code treating q as pt.wise values
            if ((.not. quad) .or. (.not. enufNbor)) then
               a(irow,1) = (xn - x0)
               a(irow,2) = (yn - y0)
               go to 21
            endif
c     a(irow,3) = .5*a(irow,1)*a(irow,1)
c     a(irow,4) =    a(irow,1)*a(irow,2)
c     a(irow,5) = .5*a(irow,2)*a(irow,2)

c     # code to do quadratic reconstruction preserving integral avgs.
            a(irow,1) = 0.d0
            a(irow,2) = 0.d0
            a(irow,3) = 0.d0
            a(irow,4) = 0.d0
            a(irow,5) = 0.d0

c            # handle cut cells first
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
 15            continue
            else
               do 17 index = 1, 5
 17               wt(index,kn) = 1.d0/6.d0
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

 21            do m = 1, nvar
                  b(irow,m) = qp(ixn,iyn,m) - qp(ix0,iy0,m)
               end do
 22         continue
c
c  ## in linear case try adding ghost cut cell into gradient (if is a cut cell)
!--          if (ghost_ccg .and.  k .ne. lstgrd) then  ! same name as in flowCart
!--            irow = irow + 1
!--            ddotn = (bxpt-x0)*alf + (bypt-y0)*beta  ! dist from centroid to bndry in normal dir
!--            distx = 2.*ddotn*alf
!--            disty = 2.*ddotn*beta
!--            a(irow,1) = distx
!--            a(irow,2) = disty
!--            b(irow,1) = 0.
!--            b(irow,4) = 0.
!--            velnorm = qp(ix0,iy0,2)*alf + qp(ix0,iy0,3)*beta
!--            ubnd = qp(ix0,iy0,2)-2.d0*velnorm*alf   ! subtract twice normal velocity at bndry to get 
!--            vbnd = qp(ix0,iy0,3)-2.d0*velnorm*beta  ! to ghost cell. leave tangential vel alone
!--            b(irow,2) = ubnd - qp(ix0,iy0,2)
!--            b(irow,3) = vbnd - qp(ix0,iy0,3)
!--c fake it as a test of accuracy of implementaitno
!--            ptx = x0 + distx
!--            pty = y0 + disty
!--            call p19old(ptx,pty,rho,u,v,pr)
!--            b(irow,1) = rho-qp(ix0,iy0,1)
!--            b(irow,2) = u-qp(ix0,iy0,2)
!--            b(irow,3) = v-qp(ix0,iy0,3)
!--            b(irow,4) = pr-qp(ix0,iy0,4)
!--          endif

          if (((.not.quad .or. .not.enufNbor)) .and. (irow .lt. 2)) then
c          write(*,*) "no slope for cut cell ",ix0,iy0
           go to 110    ! not enough neighbors for gradient recon
          endif
c 
          do 30 it = 1, irow
          do 30 jt = 1, nterms
             at(jt,it) = a(it,jt)
 30       continue
c     
          do 40  m = 1, nvar
             rhsmax(m) = b(1,m)
             rhsmin(m) = b(1,m)
             do 40 it = 1, irow
                rhs(it,m) = b(it,m)
                rhsmax(m) = dmax1(rhsmax(m),b(it,m))
                rhsmin(m) = dmin1(rhsmin(m),b(it,m))
 40          continue
c     
             do 50 it = 1, nterms
                do 50 jt = 1, nterms
                   c(it,jt) = 0.d0
                   do m = 1, nvar
                      d(it,m)  = 0.d0
                   end do
                   do 45 kt = 1, irow
                      c(it,jt) = c(it,jt) + at(it,kt)*a(kt,jt)
                      do m = 1, nvar
                      d(it,m) = d(it,m) + at(it,kt)*b(kt,m)
                      end do
 45                continue
 50             continue

          if (quad .and. enufNbor) then
c
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



                w(5,m) = b(5,m) / a(5,5)
                w(4,m) = (b(4,m)-a(4,5)*w(5,m))/a(4,4)
                w(3,m) = (b(3,m)-a(3,5)*w(5,m)-a(3,4)*w(4,m))/a(3,3)
                w(2,m) = (b(2,m)-a(2,5)*w(5,m)-a(2,4)*w(4,m)
     &               -a(2,3)*w(3,m))/a(2,2)
                w(1,m) = (b(1,m)-a(1,5)*w(5,m)-a(1,4)*w(4,m)
     &               -a(1,3)*w(3,m)-a(1,2)*w(2,m))/a(1,1)

                qyy(ix0,iy0,m) =  w(5,m)
                qxy(ix0,iy0,m) =  w(4,m)
                qxx(ix0,iy0,m) =  w(3,m)
                qy(ix0,iy0,m)  =  w(2,m)
                qx(ix0,iy0,m)  =  w(1,m)
c
 60       continue

        else
c
c         here do linear fit

c     # solve C*w = d for least squares slopes. use cholesky
c     # put factors back in a
           a(1,1) = dsqrt(c(1,1))
           a(1,2) = c(1,2)/a(1,1)
           a(2,2) = dsqrt(c(2,2)-a(1,2)**2)
c     
           do 61 m = 1, nvar
c     
c     # at*a = c. solve at*b = d, aw = b.  reuse b.
              b(1,m) = d(1,m) / a(1,1)
              b(2,m) = (d(2,m) - a(1,2)*b(1,m)) / a(2,2)
              w2 =   b(2,m)/a(2,2)
              qy(ix0,iy0,m) =  w2
              qx(ix0,iy0,m) = (b(1,m)-a(1,2)*w2)/a(1,1)
c     
              qyy(ix0,iy0,m) = 0.
              qxy(ix0,iy0,m) = 0.
              qxx(ix0,iy0,m) = 0.
 61        continue

        endif

c     turn off limiting 
        if (nolimiter)  go to 110

c
c
c  check how far from zero is normal velocity at the boundary
c
c      ub  = qp(ix0,iy0,2) + ((bxpt-x0)*qx(ix0,iy0,2) +
c    &                       (bypt-y0)*qy(ix0,iy0,2))
c      vb  = qp(ix0,iy0,3) + ((bxpt-x0)*qx(ix0,iy0,3) +
c    &                       (bypt-y0)*qy(ix0,iy0,3))

c      resid = alf*ub + beta*vb
c      write(21,*) " pt.",ix0,iy0," resid = ",resid
c
c   doesn't seem to help - skip this (go on to limiting)
       go to 63
c
c   want normal velocity at the wall to be zero. solve minimal norm
c   underdetermined least squares problem
!--       c1 =  alf*(bxpt - x0)
!--       c2 =  alf*(bypt - y0)
!--       c3 = beta*(bxpt - x0)
!--       c4 = beta*(bypt - y0)
!--       rn =  sqrt((bxpt-x0)**2 + (bypt-y0)**2)
!--       if (c1 .ge. 0.) then
!--          bet = rn
!--       else
!--
!--          bet = -rn
!--       endif
!--       v1 = c1 + bet
!--       v2 = c2 
!--       v3 = c3 
!--       v4 = c4 
!--c     solve
!--       z1 = resid / rn
!--c   
!--c    recover change in slopes, delta, by delta = householder * z
!--       vtv = v1*v1+v2*v2+v3*v3+v4*v4
!--       v1z1 = v1 * z1
!--       coeff = 2. * v1z1 / vtv
!--       delta1 = z1 - coeff * v1
!--       delta2 =     -coeff * v2
!--       delta3 =     -coeff * v3
!--       delta4 =     -coeff * v4
!--       delnorm2 = (delta1**2+delta2**2+delta3**2+delta4)
!--       gnorm2 = (qx(ix0,iy0,2)**2+qy(ix0,iy0,2)**2+
!--     &      qx(ix0,iy0,3)**2+qy(ix0,iy0,3)**2)
!--       if (gnorm .ne. 0.) then
!--          reldel = delnorm / gnorm
!--       else
!--          reldel = delnorm
!--       endif
!--
!--       qx(ix0,iy0,2) = qx(ix0,iy0,2) - delta1
!--       qy(ix0,iy0,2) = qy(ix0,iy0,2) - delta2
!--       qx(ix0,iy0,3) = qx(ix0,iy0,3) - delta3
!--       qy(ix0,iy0,3) = qy(ix0,iy0,3) - delta4
!--
!--c     check boundary residual again
!--       ub  = qp(ix0,iy0,2) + ((bxpt-x0)*qx(ix0,iy0,2) +
!--     &      (bypt-y0)*qy(ix0,iy0,2))
!--       vb  = qp(ix0,iy0,3) + ((bxpt-x0)*qx(ix0,iy0,3) +
!--     &      (bypt-y0)*qy(ix0,iy0,3))
!--       resid = alf*ub + beta*vb
!--c     write(21,*) " pt.",ix0,iy0," delnrm2 ",delnorm2, " gnrm2", gnorm2
!--
!--c
!--c      go to 110 

c      ::: now some kind of limiting of more accurate slope computed above
c      ::: fix up matrix first - transpose remained
c     
 63    do 65 i = 1, irow
         a(i,1) = at(1,i)
         a(i,2) = at(2,i)
 65    continue
c
c      ::::: check if solution feasible (doesn't need limiting)
c            see if exceed neighboring values at corners of cell
c
       if (k .eq. lstgrd) then  ! put in poly for this regular cell
          call makep(poly(1,1,lstgrd),ix0,iy0,xlow,ylow,hx,hy)
       endif

       if (.false.) then !limit at midpoints of cell edges
       do 95 m = 1, nvar
          frac = 1.d0
          do 66 iside = 1, 10
             if (poly(iside+1,1,k) .eq. -11) go to 67
             x = poly(iside,1,k)
             y = poly(iside,2,k)
             xdif = (x-x0)
             ydif = (y-y0)
             if (quad) then
                dd = xdif*qx(ix0,iy0,m)+ydif*qy(ix0,iy0,m)
     .               + (.5d0*xdif*xdif-shiftxx)*qxx(ix0,iy0,m)
     .               + (.5d0*ydif*ydif-shiftyy)*qyy(ix0,iy0,m)
     .               +      (xdif*ydif-shiftxy)*qxy(ix0,iy0,m)
             else
                dd = xdif*qx(ix0,iy0,m)+ydif*qy(ix0,iy0,m)
             endif
             if (dd .gt. rhsmax(m)) frac=dmin1(frac,rhsmax(m)/dd)
             if (dd .lt. rhsmin(m)) frac=dmin1(frac,rhsmin(m)/dd)
             if ((dd .lt. 0.d0 .and. rhsmin(m) .gt. 0.d0) .or.
     .            (dd .gt. 0.d0 .and. rhsmax(m) .lt. 0.d0)) frac = 0.d0
 66       continue
 67       continue
c        qx(ix0,iy0,m) = frac*qx(ix0,iy0,m)
c        qy(ix0,iy0,m) = frac*qy(ix0,iy0,m)
c
c      ::: check for peaks, i.e. that doesn't predict neighboring
c      ::: states in the wrong direction
c      do 68 i = 1, irow
c	       if (quad) then
c	        dd = a(i,1)*qx(ix0,iy0,m) + a(i,2)*qy(ix0,iy0,m)
c     .                             + a(i,3)*qxx(ix0,iy0,m)
c     .                             + a(i,4)*qxy(ix0,iy0,m)
c     .                             + a(i,5)*qyy(ix0,iy0,m)
c	       else
c	        dd = a(i,1)*qx(ix0,iy0,m) + a(i,2)*qy(ix0,iy0,m)
c	       endif
c	       if (dd*rhs(i,m) .lt. 0.d0) frac = 0.d0
c 68   continue
       qx(ix0,iy0,m) = frac*qx(ix0,iy0,m)
       qy(ix0,iy0,m) = frac*qy(ix0,iy0,m)
       if (quad) then
          qxx(ix0,iy0,m) = frac*qxx(ix0,iy0,m)
          qxy(ix0,iy0,m) = frac*qxy(ix0,iy0,m)
          qyy(ix0,iy0,m) = frac*qyy(ix0,iy0,m)
       endif
c     
 95   continue
      endif
c   
c do it again to limit against neighboring cell centers. (as in minmod vs previous mc limiter)
       ! original version
       if (.false.) then
       do 96 m = 1, nvar
          frac = 1.d0
          do 69 it = 1, irow
             ixn = nlist(1+it,1)
             iyn = nlist(1+it,2)

             xdif = a(it,1)  ! (x-x0)
             ydif = a(it,2)  ! (y-y0)

             dd = xdif*qx(ix0,iy0,m)+ydif*qy(ix0,iy0,m)
             qdif = qp(ixn,iyn,m) - qp(ix0,iy0,m)

             if (dd*qdif .le. 0.) then
                frac = 0.d0
             else
                frac = dmin1(frac, qdif/dd)
             endif
 69       continue
          qx(ix0,iy0,m) = frac*qx(ix0,iy0,m)
          qy(ix0,iy0,m) = frac*qy(ix0,iy0,m)
 96   continue
      endif

       ! BJ limiting
       dumax = 0.0d0
       dumin = 0.0d0
       phimin = 1.d0
       do it = 1, irow
          ixn = nlist(1+it,1)
          iyn = nlist(1+it,2)
          dumax = max(dumax(:),qp(ixn,iyn,:) - qp(ix0,iy0,:))
          dumin = min(dumin(:),qp(ixn,iyn,:) - qp(ix0,iy0,:))
       end do
       do 107 it = 1, irow
         ixn = nlist(1+it,1)
         iyn = nlist(1+it,2)
         xdif = a(it,1)  ! (x-x0)
         ydif = a(it,2)  ! (y-y0)

         do 97 m = 1, nvar
             dd = xdif*qx(ix0,iy0,m)+ydif*qy(ix0,iy0,m)

             if (dd .lt. 0.d0) then
                frac = min(1.d0, dumin(m)/dd)
             else if (dd .gt. 0.d0) then
                frac = min(1.d0, dumax(m)/dd)
             else
                frac = 1.d0
             endif
           phimin(m) = min(frac,phimin(m))
 97        continue
 107      continue
          qx(ix0,iy0,:) = qx(ix0,iy0,:)*phimin(:)
          qy(ix0,iy0,:) = qy(ix0,iy0,:)*phimin(:)

          if(dsqrt(qx(ix0,iy0,1)*qx(ix0,iy0,1)
     .             + qy(ix0,iy0,1)*qy(ix0,iy0,1)) > 1e-10) then
            dlimit(ix0,iy0) = 1.d0-phimin(1)
          else
            dlimit(ix0,iy0) = 0.d0
          endif

 110  continue
c
      if (prflag) then
         write(21,*)' qx '
         do 180 i = 2, mitot-1
         do 180 j = 2, mjtot-1
            if (irr(i,j) .ne. -1) then
               write(21,190)i,j,(qx(i,j,m),m=1,nvar)
 190           format('  i,j  ',2i4,4e14.6)
            endif
 180     continue
         write(21,*)' qy '
         do 181 i = 2, mitot-1
         do 181 j = 2, mjtot-1
            if (irr(i,j) .ne. -1) then
               write(21,190)i,j,(qy(i,j,m),m=1,nvar)
            endif
 181     continue
         if (quad) then
            write(21,*)' qxx '
            do 182 i = 2, mitot-1
            do 182 j = 2, mjtot-1
               if (irr(i,j) .ne. -1) then
                  write(21,190)i,j,(qxx(i,j,m),m=1,nvar)
               endif
 182        continue
            write(21,*)' qxy '
            do 183 i = 2, mitot-1
            do 183 j = 2, mjtot-1
               if (irr(i,j) .ne. -1) then
                  write(21,190)i,j,(qxy(i,j,m),m=1,nvar)
               endif
 183        continue
            write(21,*)' qyy '
            do 184 i = 2, mitot-1
            do 184 j = 2, mjtot-1
               if (irr(i,j) .ne. -1) then
                  write(21,190)i,j,(qyy(i,j,m),m=1,nvar)
               endif
 184        continue
         endif
      endif
c
 99    return
       end
