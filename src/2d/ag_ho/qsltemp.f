c
c ---------------------------------------------------------------------
c
       subroutine qsltemp(qp,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,
     &                    hx,hy,xlow,ylow,maxip1,maxjp1)
       implicit double precision(a-h,o-z)
       parameter (msize = 264)
c      parameter (msize = 564)
       dimension qp(maxip1,maxjp1,4),qx(mitot,mjtot,4),
     &           qy(mitot,mjtot,4),irr(mitot,mjtot)
       common /order2/ ssw, quad
      include "cirr.i"
       dimension rhsmax(4), rhsmin(4)
       dimension a(25,5),at(5,25),c(5,5),b(25,4),d(25,4),rhs(25,4)
       dimension w(5,4)
       dimension nlist(25,2)
       logical   prflag, quad
       data      prflag/.false./

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

       do 110 ix0 = lwidth+1, mitot-lwidth
       do 110 iy0 = lwidth+1, mjtot-lwidth
          k = irr(ix0,iy0)
 	  if (k .eq. -1) go to 110
c         ::: next line is so only cut cells done quadratically,
c         ::: otherwise do their neighbors too
c	  if (k .eq. lstgrd) go to 110
 	  if (k .eq. lstgrd .and. 
     .        irr(ix0+1,iy0) .eq. lstgrd .and. 
     .        irr(ix0,iy0+1) .eq. lstgrd .and. 
     .        irr(ix0,iy0-1) .eq. lstgrd .and. 
     .        irr(ix0-1,iy0) .eq. lstgrd) go to 110
c         ::: use one-sided 2nd order slopes 
 	  if (k .eq. lstgrd .and. .not. quad) go to 110
c
c         # this cell needs a slope
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
 39	    rlenb = dsqrt((sidey1-sidey2)**2 + (sidex1-sidex2)**2)
	    alf  = (sidey1-sidey2)/rlenb
	    beta = (sidex2-sidex1)/rlenb
	    bxpt = .5d0*(sidex1+sidex2)
	    bypt = .5d0*(sidey1+sidey2)
	  else
	    x0 = xlow + (ix0-.5d0)*hx
	    y0 = ylow + (iy0-.5d0)*hy
	  endif

          nlist(1,1) = ix0
          nlist(1,2) = iy0
          nst        = 1  
          nend       = 1  
          call addtemp(irr,nlist,nst,nend,newend,mitot,mjtot,
     .                      lstgrd,quad,lwidth)
	  if (.not. quad) then
	    if (newend .gt. 2) go to 16
c           not enough points for slope - set to 0
	    qx(ix0,iy0,1) = 0.d0
	    qx(ix0,iy0,2) = 0.d0
	    qx(ix0,iy0,3) = 0.d0
	    qx(ix0,iy0,4) = 0.d0
	    qy(ix0,iy0,1) = 0.d0
	    qy(ix0,iy0,2) = 0.d0
	    qy(ix0,iy0,3) = 0.d0
	    qy(ix0,iy0,4) = 0.d0
	    go to 110
	  endif
          if (newend .gt. 7) go to 16
c         ::: not enough neighbors for quadratic fit
c         ::: add neighbors of neighbors    
          nst        = 2  
          nend       = newend  
          call addtemp(irr,nlist,nst,nend,newend,mitot,mjtot,
     .                      lstgrd,quad,lwidth)
          if (newend .gt. 7) go to 16
c           not enough points for slope - set to 0
	    do 543 iv = 1, 4
	      qx(ix0,iy0,iv) = 0.d0
 543          qy(ix0,iy0,iv) = 0.d0
            go to 110
c
 16	  irow = 0
	  if (quad) then
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
c            # code treating q as pt.wise values
             if (.not. quad) then
	       a(irow,1) = (xn - x0)/hx
	       a(irow,2) = (yn - y0)/hy
	       go to 21
	     endif
c	     a(irow,3) = .5*a(irow,1)*a(irow,1)
c	     a(irow,4) =    a(irow,1)*a(irow,2)
c	     a(irow,5) = .5*a(irow,2)*a(irow,2)

c            # code to do quadratic reconstruction preserving integral avgs.
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
	       a(irow,1) = wt(index,kn)*(xp-x0)/hx + a(irow,1)
	       a(irow,2) = wt(index,kn)*(yp-y0)/hy + a(irow,2)
	       a(irow,3) = .5*wt(index,kn)*(xp-x0)*(xp-x0)/(hx*hx)
     &                      + a(irow,3)
	       a(irow,4) = wt(index,kn)*(xp-x0)*(yp-y0)/(hx*hy)
     &                      + a(irow,4)
	       a(irow,5) = .5*wt(index,kn)*(yp-y0)*(yp-y0)/(hy*hy)
     &                      + a(irow,5)
 15          continue
	     else
	     do 17 index = 1, 5
 17           wt(index,kn) = 1.d0/6.d0
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
	       a(irow,1) = wt(index,kn)*(xp-x0)/hx + a(irow,1)
	       a(irow,2) = wt(index,kn)*(yp-y0)/hy + a(irow,2)
	       a(irow,3) = .5*wt(index,kn)*(xp-x0)*(xp-x0)/(hx*hx)
     &                      + a(irow,3)
	       a(irow,4) = wt(index,kn)*(xp-x0)*(yp-y0)/(hx*hy)
     &                      + a(irow,4)
	       a(irow,5) = .5*wt(index,kn)*(yp-y0)*(yp-y0)/(hy*hy)
     &                      + a(irow,5)
 19          continue

	     endif

c            #  shift to fit quadratic terms with mean 0 over cell k
	     a(irow,3) = a(irow,3) - shiftxx
	     a(irow,4) = a(irow,4) - shiftxy
	     a(irow,5) = a(irow,5) - shiftyy

 21          b(irow,1) = qp(ixn-lwidth+1,iyn-lwidth+1,1) - 
     .                   qp(ix0-lwidth+1,iy0-lwidth+1,1)
	     b(irow,2) = qp(ixn-lwidth+1,iyn-lwidth+1,2) - 
     .                   qp(ix0-lwidth+1,iy0-lwidth+1,2)
	     b(irow,3) = qp(ixn-lwidth+1,iyn-lwidth+1,3) -
     .                   qp(ix0-lwidth+1,iy0-lwidth+1,3)
	     b(irow,4) = qp(ixn-lwidth+1,iyn-lwidth+1,4) - 
     .                   qp(ix0-lwidth+1,iy0-lwidth+1,4)
22        continue
c 
          do 30 it = 1, irow
	  do 30 jt = 1, nterms
	    at(jt,it) = a(it,jt)
 30       continue
c
	  do 40  m = 1, 4
	  rhsmax(m) = b(1,m)
	  rhsmin(m) = b(1,m)
          do 40 it = 1, irow
	    rhs(it,m) = b(it,m)
	    rhsmax(m) = dmax1(rhsmax(m),b(it,m))
	    rhsmin(m) = dmin1(rhsmin(m),b(it,m))
 40       continue
c
	  do 50 it = 1, nterms
	  do 50 jt = 1, nterms
	     c(it,jt) = 0.d0
	     d(it,1)  = 0.d0
	     d(it,2)  = 0.d0
	     d(it,3)  = 0.d0
	     d(it,4)  = 0.d0
	     do 45 kt = 1, irow
		c(it,jt) = c(it,jt) + at(it,kt)*a(kt,jt)
		d(it,1) = d(it,1) + at(it,kt)*b(kt,1)
		d(it,2) = d(it,2) + at(it,kt)*b(kt,2)
		d(it,3) = d(it,3) + at(it,kt)*b(kt,3)
		d(it,4) = d(it,4) + at(it,kt)*b(kt,4)
 45          continue
 50       continue

          if (quad) then
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
     &               /a(4,4)
          a(5,5) = dsqrt(c(5,5)-a(1,5)**2-a(2,5)**2-a(3,5)**2-a(4,5)**2)
c
          do 60 m = 1, 4
c
c            # at*a = c. solve at*b = d, aw = b.  reuse b.
             b(1,m) = d(1,m) / a(1,1)
	     b(2,m) = (d(2,m) - a(1,2)*b(1,m)) / a(2,2)
	     b(3,m) = (d(3,m) - a(1,3)*b(1,m)-a(2,3)*b(2,m))/a(3,3)
	     b(4,m) = (d(4,m) - a(1,4)*b(1,m)-a(2,4)*b(2,m)
     &                        -a(3,4)*b(3,m))/a(4,4)
	     b(5,m) = (d(5,m) - a(1,5)*b(1,m)-a(2,5)*b(2,m)
     &                        -a(3,5)*b(3,m)-a(4,5)*b(4,m))/a(5,5)



             w(5,m) = b(5,m) / a(5,5)
	     w(4,m) = (b(4,m)-a(4,5)*w(5,m))/a(4,4)
	     w(3,m) = (b(3,m)-a(3,5)*w(5,m)-a(3,4)*w(4,m))/a(3,3)
	     w(2,m) = (b(2,m)-a(2,5)*w(5,m)-a(2,4)*w(4,m)
     &                       -a(2,3)*w(3,m))/a(2,2)
	     w(1,m) = (b(1,m)-a(1,5)*w(5,m)-a(1,4)*w(4,m)
     &                       -a(1,3)*w(3,m)-a(1,2)*w(2,m))/a(1,1)

 	      qy(ix0,iy0,m) =  w(2,m)
              qx(ix0,iy0,m) =  w(1,m)
c
 60       continue

        else
c
c         here do linear fit

c         # solve C*w = d for least squares slopes. use cholesky
c         # put factors back in a
          a(1,1) = dsqrt(c(1,1))
          a(1,2) = c(1,2)/a(1,1)
          a(2,2) = dsqrt(c(2,2)-a(1,2)**2)
c
          do 61 m = 1, 4
c
c            # at*a = c. solve at*b = d, aw = b.  reuse b.
             b(1,m) = d(1,m) / a(1,1)
             b(2,m) = (d(2,m) - a(1,2)*b(1,m)) / a(2,2)
             w2 =   b(2,m)/a(2,2)
             qy(ix0,iy0,m) =  w2
             qx(ix0,iy0,m) = (b(1,m)-a(1,2)*w2)/a(1,1)
c
 61       continue

        endif

c  turn off limiting 
c         go to 110

c
c
c  check how far from zero is normal velocity at the boundary
c
c      ub  = qp(ix0,iy0,2) + ((bxpt-x0)*qx(ix0,iy0,2) +
c    &                       (bypt-y0)*qy(ix0,iy0,2))/hx
c      vb  = qp(ix0,iy0,3) + ((bxpt-x0)*qx(ix0,iy0,3) +
c    &                       (bypt-y0)*qy(ix0,iy0,3))/hx

c      resid = alf*ub + beta*vb
c      write(21,*) " pt.",ix0,iy0," resid = ",resid
c
c   doesn't seem to help - skip this (go on to limiting)
       go to 63
c
c   want normal velocity at the wall to be zero. solve minimal norm
c   underdetermined least squares problem
       c1 =  alf*(bxpt - x0)
       c2 =  alf*(bypt - y0)
       c3 = beta*(bxpt - x0)
       c4 = beta*(bypt - y0)
       rn =  sqrt((bxpt-x0)**2 + (bypt-y0)**2)
       if (c1 .ge. 0.) then
	  bet = rn
       else
	  bet = -rn
       endif
       v1 = c1 + bet
       v2 = c2 
       v3 = c3 
       v4 = c4 
c      solve
       z1 = resid / rn
c   
c    recover change in slopes, delta, by delta = householder * z
	vtv = v1*v1+v2*v2+v3*v3+v4*v4
	v1z1 = v1 * z1
	coeff = 2. * v1z1 / vtv
	delta1 = z1 - coeff * v1
	delta2 =     -coeff * v2
	delta3 =     -coeff * v3
	delta4 =     -coeff * v4
	delnorm2 = (delta1**2+delta2**2+delta3**2+delta4)
	gnorm2 = (qx(ix0,iy0,2)**2+qy(ix0,iy0,2)**2+
     &          qx(ix0,iy0,3)**2+qy(ix0,iy0,3)**2)
	if (gnorm .ne. 0.) then
	     reldel = delnorm / gnorm
	else
	     reldel = delnorm
	endif

	qx(ix0,iy0,2) = qx(ix0,iy0,2) - delta1
	qy(ix0,iy0,2) = qy(ix0,iy0,2) - delta2
	qx(ix0,iy0,3) = qx(ix0,iy0,3) - delta3
	qy(ix0,iy0,3) = qy(ix0,iy0,3) - delta4

c     check boundary residual again
       ub  = qp(ix0,iy0,2) + ((bxpt-x0)*qx(ix0,iy0,2) +
     &                       (bypt-y0)*qy(ix0,iy0,2))/hx
       vb  = qp(ix0,iy0,3) + ((bxpt-x0)*qx(ix0,iy0,3) +
     &                       (bypt-y0)*qy(ix0,iy0,3))/hx
       resid = alf*ub + beta*vb
c      write(21,*) " pt.",ix0,iy0," delnrm2 ",delnorm2, " gnrm2", gnorm2

c
c      go to 110 

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
       do 95 m = 1, 4
       frac = 1.d0
       do 66 iside = 1, 10
c	  if (poly(iside+1,1,k) .eq. -1) go to 67
	  if (poly(iside+1,1,k) .eq. -11) go to 67
	  x = poly(iside,1,k)
	  y = poly(iside,2,k)
          dd = ((x-x0)*qx(ix0,iy0,m)+(y-y0)*qy(ix0,iy0,m))/hx
          if (dd .gt. rhsmax(m)) frac=dmin1(frac,rhsmax(m)/dd)
	  if (dd .lt. rhsmin(m)) frac=dmin1(frac,rhsmin(m)/dd)
	  if ((dd .lt. 0.d0 .and. rhsmin(m) .gt. 0.d0) .or.
     .        (dd .gt. 0.d0 .and. rhsmax(m) .lt. 0.d0)) frac = 0.d0
 66    continue
 67    qx(ix0,iy0,m) = frac*qx(ix0,iy0,m)
       qy(ix0,iy0,m) = frac*qy(ix0,iy0,m)
c
c      ::: check for peaks, i.e. that doesn't predict neighboring
c      ::: states in the wrong direction
       do 68 i = 1, irow
	  dd = a(i,1)*qx(ix0,iy0,m) + a(i,2)*qy(ix0,iy0,m)
	  if (dd*rhs(i,m) .lt. 0.d0) frac = 0.d0
 68    continue
       qx(ix0,iy0,m) = frac*qx(ix0,iy0,m)
       qy(ix0,iy0,m) = frac*qy(ix0,iy0,m)
c
 95    continue
c
 110   continue
c
       if (prflag) then
          write(21,*)' qx '
          do 180 i = 2, mitot-1
          do 180 j = 2, mjtot-1
          if (irr(i,j) .ne. -1) then
	    write(21,190)i,j,(qx(i,j,m),m=1,4)
 190        format('  i,j  ',2i4,4e14.6)
          endif
 180      continue
          write(21,*)' qy '
          do 181 i = 2, mitot-1
          do 181 j = 2, mjtot-1
          if (irr(i,j) .ne. -1) then
	    write(21,190)i,j,(qy(i,j,m),m=1,4)
          endif
 181      continue
       endif
c
 99    return
       end
