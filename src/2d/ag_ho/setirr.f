c
c ------------------------------------------------------
c
      subroutine setirr(irr,mitot,mjtot,mptr,quad,gradThreshold)
c
      implicit double precision (a-h,o-z)
      include "./quadrature.i"
      dimension irr(mitot,mjtot) 


      include "cirr.i"
      include "./reconinfo.i"
      common  /cloops/ xloops(10),yloops(10),nloops
      logical            graf, quad, done
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      common /order2/ ssw, temp2, temp3
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThresholdp, ilts,ibc
      common /center/ sx,sy
      logical IS_REAL
      IS_REAL(i,j) = (i > 0 .and. i < mitot+1 .and.
     .                j > 0 .and. j < mjtot+1)
c
c lstgrd = starting index into the linked list of irregular info.
c    for grid mptr. In fact, the first entry is the "regular"
c    cell info, and stores the area of the regular cells.
c
c set the irregular cell array for the augmented grid as follows:
c -1 = cell is completely exterior to the domain
c  lstgrd = cell is a regular cell
c  k <> lstgrd but > 0 = cell is irregular, with k = index into the
c                     linked list for that cell's info.
c
c determine this by calling fbody = returns postiive if in the domain,
c  negative if outside. when a cell has a sign change, a boundary
c  segment has been found. must be careful about orientation -
c  exterior of the domain is on the right. when such a change
c  has been found,  call mkbdry to make the list for that
c  boundary piece. must also take care of completely closed
c  loops inside the grid.
c
       node(17,mptr) = lstget(dummy)
       lstgrd          = node(17,mptr)
       ar(lstgrd)    = rnode(9,mptr)*rnode(10,mptr)
       nxtirr(lstgrd)  = 0
       done            = .false.
c
c initialize irr array. mkbdry will later modify one of the cells
c at a sign change. use the loc. of the cell center for this first pass
c
       hx    = rnode(9,mptr)
       hy    = rnode(10, mptr)
       xlow  = rnode(1,mptr) - lwidth*hx
       ylow  = rnode(2,mptr) - lwidth*hy
       centx = xlow + .5*hx
       centy = ylow + .5*hy
       do 10 i = 1, mitot
       do 10 j = 1, mjtot 
          x = centx + (i-1)*hx
          y = centy + (j-1)*hy
	  d = fbody(x,y)
	  if (d .lt. 0.) then
	     irr(i,j) = -1
	  else
	     irr(i,j) = lstgrd
	  endif
 10    continue
c
c look for sign changes in irr as an indication of the starting
c or ending location of a boundary segment. only for starting
c segments, call mkbdry for the info. on that bndry segment
c as it affects this grid.  only look at the perimeter to
c find starting and ending segments. only other case is
c segment is completely contained in interior. must be more
c careful here and look at outside cell edge not center loc.
c this is an example why:
c            |   |   |   |        |___|___|___|___|___|
c            |__/|\__|___|        |___|___|__|=|__|___|
c                                 |___|___|__| |__|___|  
c                                            | |
c
c must catch boundary segments that only glancingly touch the
c fine grid, but don't cover cell center. once it's caught,
c mkbdry will fix it up in the interior, as in the ex. on the right.
c do the bottom, top, left, and right sides.
c
       kloc = lstgrd
       fpre = fbody(xlow,ylow)
       do 20 i = 1, mitot
	  fpost = fbody(xlow+i*hx,ylow)
          if ((fpost .lt. 0) .and. (fpre .gt. 0)) 
     1	       call mkbdry(kloc,i,1,irr,mitot,mjtot,xlow,ylow,
     2                     hx,hy,lstgrd,mptr,gradThreshold)
	  fpre  = fpost
 20    continue
      yhigh = ylow + mjtot*hy
      fpre  = fbody(xlow,yhigh)
       do 30 i = 1, mitot
	  fpost = fbody(xlow+i*hx,yhigh)
          if ((fpost .gt. 0) .and. (fpre .lt. 0)) 
     1	    call mkbdry(kloc,i,mjtot,irr,mitot,mjtot,xlow,ylow,
     2                  hx,hy,lstgrd,mptr,gradThreshold)
	  fpre = fpost
 30    continue
       fpre = fbody(xlow,ylow)
       do 40 j = 1, mjtot
	  fpost = fbody(xlow,ylow+j*hy)
          if ((fpost .gt. 0) .and. (fpre .lt. 0)) 
     1	       call mkbdry(kloc,1,j,irr,mitot,mjtot,xlow,ylow,
     2                     hx,hy,lstgrd,mptr,gradThreshold)
	  fpre = fpost
 40    continue
       xhigh = xlow + mitot*hx
       fpre = fbody(xhigh,ylow)
       do 50 j = 1, mjtot
	  fpost = fbody(xhigh,ylow+j*hy)
          if ((fpost .lt. 0) .and. (fpre .gt. 0)) 
     1	    call mkbdry(kloc,mitot,j,irr,mitot,mjtot,xlow,ylow,
     2                  hx,hy,lstgrd,mptr,gradThreshold)
	  fpre = fpost
 50    continue
c
c last job is to check for closed loops completely contained inside
c the grid, so they might not have intersected the boundary.
c also, double check that every boundary segment accounted for.
c
       do 60 iloops = 1, nloops
          x = xloops(iloops)
          y = yloops(iloops)
          if ((x .lt. xlow) .or. (x .gt. xhigh) .or. (y .lt. ylow) .or.
     1      (y .gt. yhigh)) go to 60
c         point inside the augmented grid. has it been acounted for yet?
          i = (x - xlow)/hx + 1.00d0
          j = (y - ylow)/hy + 1.00d0
	  if ((irr(i,j) .eq. lstgrd) .or. (irr(i,j) .eq. -1))
     1	  call mkbdry(kloc,i,j,irr,mitot,mjtot,xlow,ylow,
     2                hx,hy,lstgrd,mptr,gradThreshold)
 60    continue
c
c  for each cell with a boundary segment, calculate the weight
c  needed for quadratic reconstruction if nec.
c
c     if (quad) then
c     ###  compute weights all the time, even if not used
c
	  do 70 i = 1, mitot
	  do 70 j = 1, mjtot
	   k = irr(i,j)
	   if (k .eq. -1) go to 70
	   if (k .ne. lstgrd) then
              call weights(poly(1,1,k),xcirr(k),ycirr(k),totxx,totxy,
     &       	           totyy,wt(1,k),ar(k),points(1,1,k),hx,hy,k)
c             save tots in rest of poly
	      poly( 8,1,k) = totxx
	      poly( 9,1,k) = totxy
	      poly(10,1,k) = totyy
	   else if (.not. done) then
	       call makep(poly(1,1,lstgrd),i,j,xlow,ylow,hx,hy)
	       xc = (poly(1,1,lstgrd) + poly(3,1,lstgrd))/2.d0
	       yc = (poly(1,2,lstgrd) + poly(2,2,lstgrd))/2.d0
	       call weights(poly(1,1,lstgrd),xc,yc,totxx,totxy,totyy,wt,
     .                  ar(lstgrd),points(1,1,lstgrd),hx,hy,k)
	       poly(8,1,lstgrd)  = totxx
	       poly(9,1,lstgrd)  = totxy
	       poly(10,1,lstgrd) = totyy
	       done = .true.
	    endif
 70      continue
c     endif




           poly(8,1,lstgrd)  = 1.d0/12.d0
           poly(9,1,lstgrd)  = 0.d0
           poly(10,1,lstgrd) = 1.d0/12.d0
           dcubicshifts(1,lstgrd) = 0.d0
           dcubicshifts(2,lstgrd) = 0.d0
           dcubicshifts(3,lstgrd) = 0.d0
           dcubicshifts(4,lstgrd) = 0.d0

      ! recompute the moments
      do 71 i = 1, mitot
      do 71 j = 1, mjtot
       kirr = irr(i,j)
       if (kirr .eq. -1 .or. kirr .eq. lstgrd) go to 71

           shiftxx = 0.d0
           shiftxy = 0.d0
           shiftyy = 0.d0

           shiftxxx = 0.d0
           shiftxxy = 0.d0
           shiftxyy = 0.d0
           shiftyyy = 0.d0

           x0 = xcirr(kirr)
           y0 = ycirr(kirr)

          arr = ar(kirr)
          ivert = 1
          do 23 while (poly(ivert+1,1,kirr) .ne. -11.)
            ivert = ivert + 1
  23     continue




          itri = ivert - 3
          idx1 = 1

          do 21 it = 1, itri ! for each  triangle
            idx2 = it + 1
            idx3 = it + 2

            x1 = poly(idx1,1,kirr)
            y1 = poly(idx1,2,kirr)

            x2 = poly(idx2,1,kirr)
            y2 = poly(idx2,2,kirr)

            x3 = poly(idx3,1,kirr)
            y3 = poly(idx3,2,kirr)

            artri = area(x1, x2, x3, y1, y2, y3)

            do 22 itq = 1,ntriquad

                xval = x1 * rtri(itq) + x2 * stri(itq)
     .              +  x3 * (1.d0-rtri(itq)-stri(itq))
                yval = y1 * rtri(itq) + y2 * stri(itq)
     .              +  y3 * (1.d0-rtri(itq)-stri(itq))

      shiftxx = shiftxx + (artri/arr)*wtri(itq)*(xval-x0)**2 /(hx**2)
      shiftxy = shiftxy
     .          + (artri/arr)*wtri(itq)*(xval-x0)*(yval-y0)/(hx * hy)
      shiftyy = shiftyy + (artri/arr)*wtri(itq)*(yval-y0)**2 /(hy**2)


      shiftxxx = shiftxxx + (artri/arr)*wtri(itq)*(xval-x0)**3 /(hx**3)
      shiftxxy = shiftxxy +
     .          (artri/arr)*wtri(itq)*(yval-y0)*(xval-x0)**2 /(hy*hx**2)
      shiftxyy = shiftxyy +
     .          (artri/arr)*wtri(itq)*(xval-x0)*(yval-y0)**2 /(hx*hy**2)
      shiftyyy =shiftyyy + (artri/arr)*wtri(itq)*(yval-y0)**3 /(hy**3)



  22        continue ! for each quadrature point on each triangle
  21      continue ! for each triangle

           poly(8,1,kirr)  = shiftxx
           poly(9,1,kirr)  = shiftxy
           poly(10,1,kirr) = shiftyy
           dcubicshifts(1,kirr) = shiftxxx
           dcubicshifts(2,kirr) = shiftxxy
           dcubicshifts(3,kirr) = shiftxyy
           dcubicshifts(4,kirr) = shiftyyy

!            if(i .eq. 9 .and. j .eq. 24) then
!            print *,shiftxx, shiftxy, shiftyy, shiftxxx, shiftxxy,
!     . shiftxyy, shiftyyy
!            print *, x0,y0, arr
!
!          ivert = 1
!          do 2329 while (poly(ivert+1,1,kirr) .ne. -11.)
!            ivert = ivert + 1
! 2329     continue
!
!          itri = ivert - 3
!          idx1 = 1
!             write(*, 91,advance="no") '{'
! 91          format(A)
!            do ii = 1,itri+2
!      write(*, 90,advance="no") '{', poly(ii,1,kirr),',',
!     .                                      poly(ii,2,kirr),'}'
! 90          format(A,e25.16,A,e25.16,A)
!            end do
!             write(*, 92,advance="no") '};'
! 92          format(A)
!            print *, ""
!
!            endif

 71      continue









































      ! make merging hoods
       areaTOL = 0.50d0*hx*hy
       call makeHood(irr,mitot,mjtot,lwidth,lstgrd,xlow,ylow,hx,hy,
     .areaTOL)

      ! make reconstruction hoods
      if(ssw .eq. 0) then
      recontolx = 0.5d0 *hx
      recontoly = 0.5d0 *hy
      initval = 1
      ! standard reconstruction hoods
       call makeReconHood(irr, mitot, mjtot, lwidth, lstgrd, hx,hy,
     . recontolx, recontoly, xlow, ylow, initval)
       call makeFVReconHood(irr, mitot, mjtot, lwidth, lstgrd,
     . hx,hy,reconTOLx,reconTOLy, xlow, ylow, initval)

      elseif(ssw .eq. 1) then
      recontolx = 0.5d0 *hx
      recontoly = 0.5d0 *hy
      initval = 1 ! originally was 1
      ! standard reconstruction hoods
       call makeReconHood(irr, mitot, mjtot, lwidth, lstgrd, hx,hy,
     . recontolx, recontoly, xlow, ylow, initval)
       call makeFVReconHood(irr, mitot, mjtot, lwidth, lstgrd,
     . hx,hy,reconTOLx,reconTOLy, xlow, ylow, initval)

      elseif(ssw .eq. 2) then
!      recontolx = 1.5d0 *hx
!      recontoly = 1.5d0 *hy
!      initval = 2
!      ! standard reconstruction hoods
!       call makeReconHood(irr, mitot, mjtot, lwidth, lstgrd, hx,hy,
!     . recontolx, recontoly, xlow, ylow, initval)
!       call makeFVReconHood(irr, mitot, mjtot, lwidth, lstgrd,
!     . hx,hy,reconTOLx,reconTOLy, xlow, ylow, initval)
      ! quadratic slopes on merging neighborhoods, use 5x5
      recontolx = 1.5d0 *hx
      recontoly = 1.5d0 *hy
      initval = 2
      ! standard reconstruction hoods for merging neighborhoods
       call makeReconHood(irr, mitot, mjtot, lwidth, lstgrd, hx,hy,
     . recontolx, recontoly, xlow, ylow, initval)

       call makeFVReconHood_m1(irr, mitot, mjtot, lwidth, lstgrd,
     . hx,hy,
     . xlow, ylow)

      elseif(ssw .eq. 3 ) then
      recontolx = 2.5d0 *hx
      recontoly = 2.5d0 *hy
      initval = 3
      ! standard reconstruction hoods
       call makeReconHood(irr, mitot, mjtot, lwidth, lstgrd, hx,hy,
     . recontolx, recontoly, xlow, ylow, initval)
       call makeFVReconHood(irr, mitot, mjtot, lwidth, lstgrd,
     . hx,hy,reconTOLx,reconTOLy, xlow, ylow, initval)

      elseif(ssw .eq. -1) then
      ! quadratic slopes on merging neighborhoods, use 5x5
      recontolx = 1.5d0 *hx
      recontoly = 1.5d0 *hy
      initval = 2
      ! standard reconstruction hoods for merging neighborhoods
       call makeReconHood(irr, mitot, mjtot, lwidth, lstgrd, hx,hy,
     . recontolx, recontoly, xlow, ylow, initval)

       call makeFVReconHood_m1(irr, mitot, mjtot, lwidth, lstgrd,
     . hx,hy,
     . xlow, ylow)

      elseif(ssw .eq. -10) then ! pointwise quadratics
      ! quadratic slopes on merging neighborhoods, use 5x5
      recontolx = 1.5d0 *hx
      recontoly = 1.5d0 *hy
      initval = 2
      ! standard reconstruction hoods for merging neighborhoods
       call makeReconHood(irr, mitot, mjtot, lwidth, lstgrd, hx,hy,
     . recontolx, recontoly, xlow, ylow, initval)

       call makeFVReconHood_m1(irr, mitot, mjtot, lwidth, lstgrd,
     . hx,hy,
     . xlow, ylow)

      elseif(ssw .eq. -2) then
      ! cubically accurate second derivatives on merging neighborhoods, use 7x7
      recontolx = 2.5d0 *hx
      recontoly = 2.5d0 *hy
      initval = 3
      ! standard reconstruction hoods for merging neighborhoods
       call makeReconHood(irr, mitot, mjtot, lwidth, lstgrd, hx,hy,
     . recontolx, recontoly, xlow, ylow, initval)

       call makeFVReconHood_m2(irr, mitot, mjtot, lwidth, lstgrd,
     . hx,hy,
     . xlow, ylow)



!      recontolx = 2.5d0 *hx
!      recontoly = 2.5d0 *hy
!      initval = 3
!      ! standard reconstruction hoods
!       call makeReconHood(irr, mitot, mjtot, lwidth, lstgrd, hx,hy,
!     . recontolx, recontoly, xlow, ylow, initval)
!       call makeFVReconHood(irr, mitot, mjtot, lwidth, lstgrd,
!     . hx,hy,reconTOLx,reconTOLy, xlow, ylow, initval)
      elseif(ssw .eq. -3) then
      recontolx = 2.5d0 *hx
      recontoly = 2.5d0 *hy
      initval = 3
      ! standard reconstruction hoods for merging neighborhoods
       call makeReconHood(irr, mitot, mjtot, lwidth, lstgrd, hx,hy,
     . recontolx, recontoly, xlow, ylow, initval)

       call makeFVReconHood_m2(irr, mitot, mjtot, lwidth, lstgrd,
     . hx,hy,
     . xlow, ylow)
      endif








      ! noooooow compute the merged cell shifts

      ! recompute the moments
      do 719 j = lwidth+1, mjtot-lwidth
      do 719 i = lwidth+1, mitot-lwidth
       kirr = irr(i,j)
       if (kirr .eq. -1 .or. kirr .eq. lstgrd) go to 719



           shiftmxx = 0.d0
           shiftmxy = 0.d0
           shiftmyy = 0.d0

           shiftmxxx = 0.d0
           shiftmxxy = 0.d0
           shiftmxyy = 0.d0
           shiftmyyy = 0.d0

           x0 = xcentMerge(kirr)
           y0 = ycentMerge(kirr)

           arr = volMerge(kirr)

!        if(i .eq. 28 .and. j .eq. 7) then
!        print *, "here"
!        endif

      do 340 ic = 1, ncount(kirr)
            ioff = iidx(kirr,ic)
            joff = jidx(kirr,ic)

            koff = irr(ioff, joff)

      if(koff .eq. lstgrd) then ! full cell
            x1 = xlow + (dfloat(ioff)-1.d0)*hx
            y1 = ylow + (dfloat(joff)-1.d0)*hy

            x3 = xlow + (dfloat(ioff)-0.d0)*hx
            y3 = ylow + (dfloat(joff)-0.d0)*hy

            arcell = hx * hy / numhoods(ioff, joff)

      do 11 iq = 1, nquadquad
          xval = x1 * (1.d0 - rquad(iq))/2. + x3 * (1.d0 + rquad(iq))/2.
          yval = y1 * (1.d0 - squad(iq))/2. + y3 * (1.d0 + squad(iq))/2.

      shiftmxx = shiftmxx
     .       + (arcell/arr)*wquad(iq)*(xval-x0)**2 /(hx**2)
      shiftmxy = shiftmxy
     .       + (arcell/arr)*wquad(iq)*(xval-x0)*(yval-y0)/(hx * hy)
      shiftmyy = shiftmyy
     .       + (arcell/arr)*wquad(iq)*(yval-y0)**2 /(hy**2)


      shiftmxxx = shiftmxxx
     .       + (arcell/arr)*wquad(iq)*(xval-x0)**3 /(hx**3)
      shiftmxxy = shiftmxxy
     .       + (arcell/arr)*wquad(iq)*(yval-y0)*(xval-x0)**2/(hx*hx*hy)
      shiftmxyy = shiftmxyy
     .       + (arcell/arr)*wquad(iq)*(xval-x0)*(yval-y0)**2/(hx*hy*hy)
      shiftmyyy = shiftmyyy
     .       + (arcell/arr)*wquad(iq)*(yval-y0)**3/(hy**3)
 11         continue ! for each quadrature point on full cell




      else                      ! cut cell




          ivert = 1
          do 239 while (poly(ivert+1,1,koff) .ne. -11.)
            ivert = ivert + 1
 239     continue


          itri = ivert - 3
          idx1 = 1

          do 219 it = 1, itri ! for each  triangle
            idx2 = it + 1
            idx3 = it + 2

            x1 = poly(idx1,1,koff)
            y1 = poly(idx1,2,koff)

            x2 = poly(idx2,1,koff)
            y2 = poly(idx2,2,koff)

            x3 = poly(idx3,1,koff)
            y3 = poly(idx3,2,koff)

            artri = area(x1, x2, x3, y1, y2, y3) / numhoods(ioff, joff)

            do 229 itq = 1,ntriquad

                xval = x1 * rtri(itq) + x2 * stri(itq)
     .              +  x3 * (1.d0-rtri(itq)-stri(itq))
                yval = y1 * rtri(itq) + y2 * stri(itq)
     .              +  y3 * (1.d0-rtri(itq)-stri(itq))

      shiftmxx = shiftmxx + (artri/arr)*wtri(itq)*(xval-x0)**2 /(hx**2)
      shiftmxy = shiftmxy
     .          + (artri/arr)*wtri(itq)*(xval-x0)*(yval-y0)/(hx * hy)
      shiftmyy = shiftmyy + (artri/arr)*wtri(itq)*(yval-y0)**2 /(hy**2)


      shiftmxxx = shiftmxxx
     .                   + (artri/arr)*wtri(itq)*(xval-x0)**3 /(hx**3)
      shiftmxxy = shiftmxxy
     .        + (artri/arr)*wtri(itq)*(yval-y0)*(xval-x0)**2/(hy*hx**2)
      shiftmxyy = shiftmxyy
     .        + (artri/arr)*wtri(itq)*(xval-x0)*(yval-y0)**2 /(hx*hy**2)
      shiftmyyy = shiftmyyy
     .        + (artri/arr)*wtri(itq)*(yval-y0)**3 /(hy**3)

 229         continue ! for each quadrature point on each triangle
 219      continue ! for each triangle

      endif

 340     continue


           dmergeshifts(1,kirr)  = shiftmxx
           dmergeshifts(2,kirr)  = shiftmxy
           dmergeshifts(3,kirr)  = shiftmyy

           dmergeshifts(4,kirr)  = shiftmxxx
           dmergeshifts(5,kirr)  = shiftmxxy
           dmergeshifts(6,kirr)  = shiftmxyy
           dmergeshifts(7,kirr)  = shiftmyyy

 719      continue




!      call merge_shifts_ho(irr,mitot,mjtot,lwidth,hx,hy,xlow,ylow,
!     . lstgrd,
!     . volmerge, xcentmerge, ycentmerge, numhoods,
!     . iidx, jidx, ncount,
!     . qmshift)

      if(ihob .eq. 1) then
      call ho_boundary(mitot, mjtot, lstgrd, ssw, hx, hy, irr,sx,sy)
      call merge_shifts_ho(irr,mitot,mjtot,lwidth,hx,hy,xlow,ylow,
     . lstgrd, ssw)
      endif





























!     computing info for FV reconstruction


      inuf = 0
      inuforder = 0

      if(ssw .eq. 1) inuforder = 1
      if(ssw .eq.-1) inuforder = 1
      if(ssw.eq.-10) inuforder = 1
      if(ssw .eq. 2) inuforder = 1
      if(ssw .eq.-2) inuforder = 2
      if(ssw .eq. 3) inuforder = 2
      if(ssw .eq.-3) inuforder = 2



       do 33 i = 1, mitot
       do 33 j = 1, mjtot
            k = irr(i,j)
            if (k .ne. lstgrd) goto 33
            inuf(i,j) = 1

            do 333 ioff = -inuforder,inuforder
            do 333 joff = -inuforder,inuforder

                if(ssw .eq. 1 .or. ssw .eq. -1 .or. ssw .eq. -10) then
                    if(abs(ioff) .eq. 1 .and. abs(joff) .eq. 1) then
                        cycle ! skip the diagonal elements for p = 1
                    endif
                endif

                if(.not. IS_REAL(i+ioff, j+joff) ) then
                    inuf(i,j) = 0
                elseif(irr(i+ioff, j+joff) .ne. lstgrd ) then
                    inuf(i,j) = 0
                endif
 333        continue
 33    continue

       return
       end

      subroutine makeFVReconHood(irr, mitot, mjtot, lwidth, lstgrd,
     . dx,dy,reconTOLx,reconTOLy, xlow, ylow, initval)
      implicit double precision (a-h, o-z)
      include "cirr.i"
      include "reconinfo.i"
      dimension irr(mitot,mjtot)
      logical IS_REAL
      IS_REAL(i,j) = (i > 0 .and. i < mitot+1 .and.
     .                j > 0 .and. j < mjtot+1)

      iir = initval
      jjr = initval

!      do 10 j = lwidth, mjtot-lwidth+1
!      do 10 i = lwidth, mitot-lwidth+1
      do 10 j = 1, mjtot
      do 10 i = 1, mitot
         k = irr(i,j)

         if( k .eq. -1) goto 10 ! solid so do nothing

         if( k .eq. lstgrd) then
             xi = xlow + (i-.5d0)*dx
             yi = ylow + (j-.5d0)*dy
         else
             xi = xcirr(k)
             yi = ycirr(k)
         endif

!         if( i .eq. 10 .and. j .eq. 25) then
!         print *, "here"
!         endif

         if(inuf(i,j) .eq. 1 .or. k .eq. -1) cycle
!         if(i .eq. 75 .and. j .eq. 45) then
!         print *,mioff(i,j), mjoff(i,j)
!         endif
         icont = 1

         do while(icont .eq. 1)
         icont = 0

         do 20 mi = 1, iir(i,j)
         do 20 mj = 1, jjr(i,j)

            diffx = -1.d0
            diffy = -1.d0

            do 30 ioff = -mi, mi
            do 30 joff = -mj, mj

                if (.not. IS_REAL(i+ioff,j+joff)) go to 30
                koff = irr(i+ioff, j+joff)

                if (koff .eq. -1) goto 30

!                if(abs(ioff) .ne. mi .or. abs(joff) .ne. mj) goto 30

                if( koff .eq. lstgrd) then
                    xoff = xlow + (i+ioff-.5d0)*dx
                    yoff = ylow + (j+joff-.5d0)*dy
                else
                    xoff = xcirr(koff)
                    yoff = ycirr(koff)
                endif

                diffx = max(diffx, dabs(xoff- xi) )
                diffy = max(diffy, dabs(yoff- yi) )
30          continue
20       continue

         if(diffx < reconTOLx - 1.d-10) then
             iir(i,j) = iir(i,j) + 1
             icont = 1
         endif
         if(diffy < reconTOLy - 1.d-10) then
             jjr(i,j) = jjr(i,j) + 1
             icont = 1
         endif

         end do

!         print *,iir(i,j), jjr(i,j)

!         if(i .eq. 75 .and. j .eq. 45) then
!         print *,iir(i,j), jjr(i,j)
!         endif

10    continue


      end subroutine









      subroutine makeFVReconHood_m3(irr, mitot, mjtot, lwidth, lstgrd,
     . dx,dy,
     . xlow, ylow)
      implicit double precision (a-h, o-z)
      include "cirr.i"
      include "reconinfo.i"
      dimension irr(mitot,mjtot)
      logical IS_REAL, IS_GHOST, IS_BORDER
      IS_REAL(i,j) = (i > 0 .and. i < mitot +1 .and.
     .                j > 0 .and. j < mjtot +1)
      IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                 j .le. lwidth .or. j .gt. mjtot-lwidth)
      IS_BORDER(i,j) = (i .eq. 1 .or. i .eq. mitot .or.
     .                  j .eq. 1 .or. j .eq. mjtot)

!      if(.true.) then ! quadratic slopes on cc (5x5 neigh)
                      !  and eneigh of cc (3x3 neigh)
      recontolx_cc = 1.5d0 *dx
      recontoly_cc = 1.5d0 *dy

      recontolx_eneigh = 0.5d0 *dx
      recontoly_eneigh = 0.5d0 *dy


      iir = 1 ! cc eneighs are 3x3
      jjr = 1 ! cc eneighs are 3x3

      do 19 j = 1, mjtot
      do 19 i = 1, mitot
        if(inuf(i,j) .eq. 1 .or. irr(i,j) .eq. -1) cycle

        if(irr(i,j) .ne. lstgrd .or. IS_BORDER(i,j) ) then ! I'm on a cut cell or on the border
            iir(i,j) = 2 ! cc and ghost cells are 5x5
            jjr(i,j) = 2 ! cc and ghost cells are 5x5
        endif

 19   continue
!      else
!      recontolx_cc = 1.5d0 *dx
!      recontoly_cc = 1.5d0 *dy
!
!      recontolx_eneigh = 1.5d0 *dx
!      recontoly_eneigh = 1.5d0 *dy
!
!      iir = 2 ! cc and eneighs are 5x5
!      jjr = 2 ! cc and eneighs are 5x5
!      endif


      do 10 j = 1, mjtot
      do 10 i = 1, mitot
         k = irr(i,j)

         if( k .eq. -1) goto 10 ! solid so do nothing

         if( k .eq. lstgrd) then
             xi = xlow + (i-.5d0)*dx
             yi = ylow + (j-.5d0)*dy
         else
             xi = xcirr(k)
             yi = ycirr(k)
         endif

!         if( i .eq. 10 .and. j .eq. 25) then
!         print *, "here"
!         endif

         if(inuf(i,j) .eq. 1 .or. k .eq. -1) cycle
!         if(i .eq. 75 .and. j .eq. 45) then
!         print *,mioff(i,j), mjoff(i,j)
!         endif
         icont = 1

         do while(icont .eq. 1)
         icont = 0

         do 20 mi = 1, iir(i,j)
         do 20 mj = 1, jjr(i,j)

            diffx = -1.d0
            diffy = -1.d0

            do 30 ioff = -mi, mi
            do 30 joff = -mj, mj

                if (.not. IS_REAL(i+ioff,j+joff)) go to 30
                koff = irr(i+ioff, j+joff)

                if (koff .eq. -1) goto 30

!                if(abs(ioff) .ne. mi .or. abs(joff) .ne. mj) goto 30

                if( koff .eq. lstgrd) then
                    xoff = xlow + (i+ioff-.5d0)*dx
                    yoff = ylow + (j+joff-.5d0)*dy
                else
                    xoff = xcirr(koff)
                    yoff = ycirr(koff)
                endif

                diffx = max(diffx, dabs(xoff- xi) )
                diffy = max(diffy, dabs(yoff- yi) )
30          continue
20       continue

         if(irr(i,j) .ne. lstgrd .or. IS_BORDER(i,j) ) then ! I'm a cut cell, or ghost cell so use larger neighborhood
             if(diffx < reconTOLx_cc - 1.d-10) then
                 iir(i,j) = iir(i,j) + 1
                 icont = 1
             endif
             if(diffy < reconTOLy_cc - 1.d-10) then
                 jjr(i,j) = jjr(i,j) + 1
                 icont = 1
             endif
         else ! I'm a cut cell edge neighbor
             if(diffx < reconTOLx_eneigh - 1.d-10) then
                 iir(i,j) = iir(i,j) + 1
                 icont = 1
             endif
             if(diffy < reconTOLy_eneigh - 1.d-10) then
                 jjr(i,j) = jjr(i,j) + 1
                 icont = 1
             endif
         endif

         end do



10    continue


      end subroutine





















      subroutine makeFVReconHood_m2(irr, mitot, mjtot, lwidth, lstgrd,
     . dx,dy,
     . xlow, ylow)
      implicit double precision (a-h, o-z)
      include "cirr.i"
      include "reconinfo.i"
      dimension irr(mitot,mjtot)
      logical IS_REAL, IS_GHOST, IS_BORDER
      IS_REAL(i,j) = (i > 0 .and. i < mitot+1 .and.
     .                j > 0 .and. j < mjtot+1)
      IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                 j .le. lwidth .or. j .gt. mjtot-lwidth)
      IS_BORDER(i,j) = (i .eq. 1 .or. i .eq. mitot .or.
     .                  j .eq. 1 .or. j .eq. mjtot)

!     ! second order accurate accurate second derivatives on cc (7x7 neigh)
      ! and eneigh of cc (5x5 neigh)
      recontolx_cc = 2.5d0 *dx
      recontoly_cc = 2.5d0 *dy

      recontolx_eneigh = 1.5d0 *dx
      recontoly_eneigh = 1.5d0 *dy


      iir = 2 ! cc eneighs are 5x5
      jjr = 2 ! cc eneighs are 5x5

      do 19 j = 1, mjtot
      do 19 i = 1, mitot
        if(inuf(i,j) .eq. 1 .or. irr(i,j) .eq. -1) cycle

        if(irr(i,j) .ne. lstgrd  .or. IS_BORDER(i,j) ) then ! I'm on a cut cell
            iir(i,j) = 3 ! cc and ghost cells are 7x7
            jjr(i,j) = 3 ! cc and ghost cells are 7x7
        endif

 19   continue
!      else
!      recontolx_cc = 1.5d0 *dx
!      recontoly_cc = 1.5d0 *dy
!
!      recontolx_eneigh = 1.5d0 *dx
!      recontoly_eneigh = 1.5d0 *dy
!
!      iir = 2 ! cc and eneighs are 5x5
!      jjr = 2 ! cc and eneighs are 5x5
!      endif


      do 10 j = 1, mjtot
      do 10 i = 1, mitot
         k = irr(i,j)

         if( k .eq. -1) goto 10 ! solid so do nothing

         if( k .eq. lstgrd) then
             xi = xlow + (i-.5d0)*dx
             yi = ylow + (j-.5d0)*dy
         else
             xi = xcirr(k)
             yi = ycirr(k)
         endif

!         if( i .eq. 10 .and. j .eq. 25) then
!         print *, "here"
!         endif

         if(inuf(i,j) .eq. 1 .or. k .eq. -1) cycle
!         if(i .eq. 75 .and. j .eq. 45) then
!         print *,mioff(i,j), mjoff(i,j)
!         endif
         icont = 1

         do while(icont .eq. 1)
         icont = 0

         do 20 mi = 1, iir(i,j)
         do 20 mj = 1, jjr(i,j)

            diffx = -1.d0
            diffy = -1.d0

            do 30 ioff = -mi, mi
            do 30 joff = -mj, mj

                if (.not. IS_REAL(i+ioff,j+joff)) go to 30
                koff = irr(i+ioff, j+joff)

                if (koff .eq. -1) goto 30

!                if(abs(ioff) .ne. mi .or. abs(joff) .ne. mj) goto 30

                if( koff .eq. lstgrd) then
                    xoff = xlow + (i+ioff-.5d0)*dx
                    yoff = ylow + (j+joff-.5d0)*dy
                else
                    xoff = xcirr(koff)
                    yoff = ycirr(koff)
                endif

                diffx = max(diffx, dabs(xoff- xi) )
                diffy = max(diffy, dabs(yoff- yi) )
30          continue
20       continue

         if(irr(i,j) .ne. lstgrd  .or. IS_BORDER(i,j)  ) then ! I'm a cut cell, or ghost cell so use larger neighborhood
             if(diffx < reconTOLx_cc - 1.d-10) then
                 iir(i,j) = iir(i,j) + 1
                 icont = 1
             endif
             if(diffy < reconTOLy_cc - 1.d-10) then
                 jjr(i,j) = jjr(i,j) + 1
                 icont = 1
             endif
         else ! I'm a cut cell edge neighbor
             if(diffx < reconTOLx_eneigh - 1.d-10) then
                 iir(i,j) = iir(i,j) + 1
                 icont = 1
             endif
             if(diffy < reconTOLy_eneigh - 1.d-10) then
                 jjr(i,j) = jjr(i,j) + 1
                 icont = 1
             endif
         endif

         end do



10    continue


      end subroutine















      subroutine makeFVReconHood_m1(irr, mitot, mjtot, lwidth, lstgrd,
     . dx,dy,
     . xlow, ylow)
      implicit double precision (a-h, o-z)
      include "cirr.i"
      include "reconinfo.i"
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold, ilts,ibc
      dimension irr(mitot,mjtot)
      logical IS_REAL, IS_GHOST, IS_BORDER, IS_REAL_BORDER
      IS_REAL(i,j) = (i > 0 .and. i < mitot+1 .and.
     .                j > 0 .and. j < mjtot+1)
      IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                 j .le. lwidth .or. j .gt. mjtot-lwidth)
      IS_BORDER(i,j) = (i .eq. 1 .or. i .eq. mitot .or.
     .                  j .eq. 1 .or. j .eq. mjtot)
      IS_REAL_BORDER(i,j) = (i .eq. lwidth+1 .or.
     .                       i .eq. mitot-lwidth .or.
     .                       j .eq. lwidth+1 .or.
     .                       j .eq. mjtot-lwidth)
!      if(.true.) then ! quadratic slopes on cc (5x5 neigh)
                      !  and eneigh of cc (3x3 neigh)
      recontolx_cc = 1.5d0 *dx
      recontoly_cc = 1.5d0 *dy

      recontolx_eneigh = 0.5d0 *dx
      recontoly_eneigh = 0.5d0 *dy


      iir = 1 ! cc eneighs are 3x3
      jjr = 1 ! cc eneighs are 3x3

      do 19 j = 1, mjtot
      do 19 i = 1, mitot
        if(inuf(i,j) .eq. 1 .or. irr(i,j) .eq. -1) cycle

        if(irr(i,j) .ne. lstgrd .or. IS_BORDER(i,j) ) then ! I'm on a cut cell
            iir(i,j) = 2 ! cc and ghost cells are 5x5
            jjr(i,j) = 2 ! cc and ghost cells are 5x5
        endif

        if((irr(i,j) .ne. lstgrd .or. IS_REAL_BORDER(i,j)) .and.
     .     ibc .eq. 1 ) then ! I'm on a cut cell
            iir(i,j) = 2 ! cc and ghost cells are 5x5
            jjr(i,j) = 2 ! cc and ghost cells are 5x5
        endif

 19   continue
!      else
!      recontolx_cc = 1.5d0 *dx
!      recontoly_cc = 1.5d0 *dy
!
!      recontolx_eneigh = 1.5d0 *dx
!      recontoly_eneigh = 1.5d0 *dy
!
!      iir = 2 ! cc and eneighs are 5x5
!      jjr = 2 ! cc and eneighs are 5x5
!      endif


      do 10 j = 1, mjtot
      do 10 i = 1, mitot
         k = irr(i,j)

         if( k .eq. -1) goto 10 ! solid so do nothing

         if( k .eq. lstgrd) then
             xi = xlow + (i-.5d0)*dx
             yi = ylow + (j-.5d0)*dy
         else
             xi = xcirr(k)
             yi = ycirr(k)
         endif

!         if( i .eq. 10 .and. j .eq. 25) then
!         print *, "here"
!         endif

         if(inuf(i,j) .eq. 1 .or. k .eq. -1) cycle
!         if(i .eq. 75 .and. j .eq. 45) then
!         print *,mioff(i,j), mjoff(i,j)
!         endif
         icont = 1

         do while(icont .eq. 1)
         icont = 0

         do 20 mi = 1, iir(i,j)
         do 20 mj = 1, jjr(i,j)

            diffx = -1.d0
            diffy = -1.d0

            do 30 ioff = -mi, mi
            do 30 joff = -mj, mj

                if (ibc .eq. 1 .and. IS_GHOST(i+ioff,j+joff)) go to 30


                if (.not. IS_REAL(i+ioff,j+joff)) go to 30
                koff = irr(i+ioff, j+joff)

                if (koff .eq. -1) goto 30

!                if(abs(ioff) .ne. mi .or. abs(joff) .ne. mj) goto 30

                if( koff .eq. lstgrd) then
                    xoff = xlow + (i+ioff-.5d0)*dx
                    yoff = ylow + (j+joff-.5d0)*dy
                else
                    xoff = xcirr(koff)
                    yoff = ycirr(koff)
                endif

                diffx = max(diffx, dabs(xoff- xi) )
                diffy = max(diffy, dabs(yoff- yi) )
30          continue
20       continue

         if(irr(i,j) .ne. lstgrd .or. IS_BORDER(i,j) ) then ! I'm a cut cell, or border ghost cell so use larger neighborhood
             if(diffx < reconTOLx_cc - 1.d-10) then
                 iir(i,j) = iir(i,j) + 1
                 icont = 1
             endif
             if(diffy < reconTOLy_cc - 1.d-10) then
                 jjr(i,j) = jjr(i,j) + 1
                 icont = 1
             endif
         else ! I'm a cut cell edge neighbor
             if(diffx < reconTOLx_eneigh - 1.d-10) then
                 iir(i,j) = iir(i,j) + 1
                 icont = 1
             endif
             if(diffy < reconTOLy_eneigh - 1.d-10) then
                 jjr(i,j) = jjr(i,j) + 1
                 icont = 1
             endif
         endif

         end do



10    continue


      end subroutine


c
c ----------------------------------------------------------------------------
c


      subroutine makeReconHood(irr, mitot, mjtot, lwidth, lstgrd, dx,dy,
     .reconTOLx,reconTOLy, xlow, ylow, initval)
      implicit double precision (a-h, o-z)
      include "cirr.i"
      dimension irr(mitot,mjtot)
      logical IS_REAL, IS_GHOST, IS_BORDER, IS_REAL_BORDER
      IS_REAL(i,j) = (i > 0 .and. i < mitot+1 .and.
     .                j > 0 .and. j < mjtot+1)
      IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                 j .le. lwidth .or. j .gt. mjtot-lwidth)
      IS_BORDER(i,j) = (i .eq. 1 .or. i .eq. mitot .or.
     .                  j .eq. 1 .or. j .eq. mjtot)
      IS_REAL_BORDER(i,j) = (i .eq. lwidth+1 .or.
     .                       i .eq. mitot-lwidth .or.
     .                       j .eq. lwidth+1 .or.
     .                       j .eq. mjtot-lwidth)


      mioff = initval
      mjoff = initval

      do 10 j = lwidth+1, mjtot-lwidth
      do 10 i = lwidth+1, mitot-lwidth
         k = irr(i,j)

         if( k .eq. -1 .or. k .eq. lstgrd) goto 10

         icont = 1
!         if(i .eq. 75 .and. j .eq. 45) then
!         print *,mioff(i,j), mjoff(i,j)
!         endif

         do while(icont .eq. 1)
         icont = 0
         do 20 mi = 1, mioff(i,j)
         do 20 mj = 1, mjoff(i,j)
            k    = irr(i,j)
            diffx = -1.d0
            diffy = -1.d0
            do 30 ioff = -mi, mi
            do 30 joff = -mj, mj
                if (IS_GHOST(i+ioff,j+joff)) go to 30
                koff = irr(i+ioff, j+joff)

                if (koff .eq. -1) goto 30

!                if(abs(ioff) .ne. mi .or. abs(joff) .ne. mj) goto 30

                if( koff .eq. lstgrd) then
                    xoff = xlow + (i+ioff-.5d0)*dx
                    yoff = ylow + (j+joff-.5d0)*dy
                else
                    xoff = xcentmerge(koff)
                    yoff = ycentmerge(koff)
                endif

                diffx = max(diffx, dabs(xoff- xcentmerge(k)) )
                diffy = max(diffy, dabs(yoff- ycentmerge(k)) )
!         if(i .eq. 75 .and. j .eq. 45) then
!         print *,xoff, yoff, xcentmerge(k), ycentmerge(k), recontolx
!         endif
30          continue
20          continue
            if(diffx < reconTOLx - 1.d-10) then
                mioff(i,j) = mioff(i,j) + 1
                icont = 1
            endif
            if(diffy < reconTOLy - 1.d-10) then
                mjoff(i,j) = mjoff(i,j) + 1
                icont = 1
            endif

!         if(i .eq. 75 .and. j .eq. 45) then
!         print *,mioff(i,j), mjoff(i,j)
!         endif
         end do


10    continue


      end subroutine



c
c ----------------------------------------------------------------------------
c
      subroutine makeHood(irr,mitot,mjtot,lwidth,lstgrd,xlow,ylow,dx,dy,
     .areaTOL)
      implicit double precision (a-h, o-z)
      include "cirr.i"
      dimension irr(mitot,mjtot)

      ! there are different strategies, how about only merging in y direction.
      ! on simple geometries there will not be any overlaps
      numHoods = 1
      xcentmerge = 0.d0
      ycentmerge = 0.d0
      volmerge = 0.d0
      xcentmerge_ho = 0.d0
      ycentmerge_ho = 0.d0
      volmerge_ho = 0.d0
      ncount = 1
      nsize = 1
      iidx = -1
      jidx = -1

      do 10 j = lwidth+1, mjtot-lwidth
      do 10 i = lwidth+1, mitot-lwidth
         k = irr(i,j)

         if( k .eq. -1 .or. k .eq. lstgrd) goto 10

         iidx(k, 1 ) = i
         jidx(k, 1 ) = j

         icurr = i
         jcurr = j

         vqmerge = ar(k)
         if (k .eq. -1 .or. k .eq. lstgrd) goto 10
         if (ar(k) > areaTOL) goto 10

        call determineDirection(areaTOL,irr,mitot,mjtot,lwidth,nvar,i,
     .j,idir)

         if(idir .eq. -1) then


        print *, "couldn't find an appropriate merging direction - need
     . to properly implement a centered merge"
        call exit(-1)
!         ir = 1
!         do while (vqmerge < areaTOL) ! merging strategy, merge up
!         do 33 ioff = -ir, ir
!         do 33 joff = -ir, ir
!            if(abs(ioff) .ne. ir .or. abs(joff) .ne. ir) goto 33
!            kcurr = irr(i+ioff, j+joff)
!            if( kcurr .eq. -1 ) goto 33
!
!            ncount(k) = ncount(k)+1
!            iidx(k, ncount(k) ) = i+ioff
!            jidx(k, ncount(k) ) = j+joff
!            numHoods(i+ioff,j+joff) = numHoods(i+ioff,j+joff) + 1
!            vqmerge = vqmerge + ar(kcurr)
! 33      continue
!         end do


         else

         do while (vqmerge < areaTOL) ! merging strategy, merge up

            if(idir .eq. 1) then
                icurr = icurr - 1
            elseif(idir .eq. 2) then
                icurr = icurr + 1
            elseif(idir .eq. 3) then
                jcurr = jcurr - 1
            elseif(idir .eq. 4) then
                jcurr = jcurr + 1
            endif

            kcurr = irr(icurr,jcurr)
            if(kcurr .eq. -1) then
                print *, "uh oh gone into a solid wall.",i,j
            return
            endif


            ncount(k) = ncount(k)+1
            iidx(k, ncount(k) ) = icurr
            jidx(k, ncount(k) ) = jcurr
            numHoods(icurr,jcurr) = numHoods(icurr,jcurr) +1
            vqmerge = vqmerge + ar(kcurr)
         end do

         endif

 10   continue

      do 30 j = lwidth+1, mjtot-lwidth
      do 30 i = lwidth+1, mitot-lwidth
         k = irr(i,j)

         if(k .eq. -1 .or. k .eq. lstgrd) goto 30

         do 20 ic = 1, ncount(k)
            icurr = iidx(k,ic)
            jcurr = jidx(k,ic)
            kcurr = irr(icurr,jcurr)
            call getCellCentroid_lo(lstgrd,icurr,jcurr,xc,yc,xlow,ylow,
     .                           dx,dy,kcurr)
            volMerge(k) = volMerge(k)
     .                    + ar(kcurr) / numHoods(icurr,jcurr)
            xcentmerge(k) = xcentmerge(k)
     .                    + xc*ar(kcurr) / numHoods(icurr,jcurr)
            ycentmerge(k) = ycentmerge(k)
     .                    + yc*ar(kcurr) / numHoods(icurr,jcurr)



 20      continue

         xcentmerge(k) = xcentmerge(k) / volMerge(k)
         ycentmerge(k) = ycentmerge(k) / volMerge(k)

 30   continue
      end


c     MOST NORMAL DIRECTION
c ----------------------------------------------------------------------------
c
      subroutine determineDirection(areaTol,irr,mitot,mjtot,lwidth,
     .nvar,i,j,idir)
      implicit double precision (a-h, o-z)
      include "cirr.i"
      dimension irr(mitot,mjtot)
      dimension area(4), inum(4)
      logical IS_GHOST
            IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                 j .le. lwidth .or. j .gt. mjtot-lwidth)

!        if( i .eq. 38 .and. j .eq. 93) then
!        print *, "here"
!        endif


        ! determine irregular face normal
        ! compute volume in the four directions
        ! ignore volumes that are < tolerance

         k = irr(i,j)
         do 20 kside=1,6
            if (poly(kside+2,1,k).eq.-11.) then
               x1 = poly(kside,1,k)
               y1 = poly(kside,2,k)
               x2 = poly(kside+1,1,k)
               y2 = poly(kside+1,2,k)
               go to 25
            endif
 20      continue
 25      continue
c
c        # boundary segment face:
c        ------------------------

c     # compute (vx,vy) = unit normal to boundary pointing in.
         hsx1 = x2
         hsx2 = x1
         hsy1 = y2
         hsy2 = y1
         rlen = dsqrt((hsy1-hsy2)**2 + (hsx1-hsx2)**2)
         vx = (hsy1-hsy2)/rlen
         vy = (hsx2-hsx1)/rlen




        do 30 id = 1,4
            inum(id) = 0
            icurr = i
            jcurr = j
            kcurr = irr(icurr,jcurr)
            area(id) = ar(kcurr)
            do while (area(id) < areaTol)
                if(id .eq. 1) then
                    icurr = icurr - 1
                elseif(id .eq. 2) then
                    icurr = icurr + 1
                elseif(id .eq. 3) then
                    jcurr = jcurr - 1
                elseif(id .eq. 4) then
                    jcurr = jcurr + 1
                endif

                if(icurr < 1 .or. jcurr < 1) then
                idir = -1
                exit
                endif

                if( icurr > mitot .or. jcurr > mjtot) then
                idir = -1
                exit
                endif


                kcurr = irr(icurr,jcurr)

                if(kcurr .eq. -1 .or. IS_GHOST(icurr, jcurr)) then
                    idir = -1
                    exit
                endif

                area(id) = area(id) + ar(kcurr)
                inum(id) = inum(id) + 1
            end do
  30    continue







        ! find the direction that has max v dot x  such that
        ! area(id) > areaTOL.

        closest = -1.d0
        idxclosest = -1
        do 50 id = 1,4
            if(area(id) > areaTOL) then
                if(id .eq. 1) then
                    dirx = -1.d0
                    diry =  0.d0
                elseif(id .eq. 2) then
                    dirx = 1.d0
                    diry =  0.d0
                elseif(id .eq. 3) then
                    dirx =  0.d0
                    diry = -1.d0
                else
                    dirx =  0.d0
                    diry =  1.d0
                endif

                temp = vx*dirx + vy*diry
                if(temp > closest) then
                    idxclosest = id
                    closest = temp
                endif

            endif
  50    continue

        if(idxclosest .eq. -1) print *, "PROBLEM with direction."
        idir = idxclosest


      return
      end

      function rs2xy_2(x1,x2,x3,x4, i)
      implicit double precision (a-h, o-z)
      include "./quadrature.i"
      r = rtri_ho(i)
      s = stri_ho(i)
      rs2xy_2 = (1.d0-r-s)*x1
     . +x2*r*(1.d0-2.d0*s)+x3*s*(1.d0-2.d0*r)+4.d0*x4*r*s
      return
      end function


      function rs2xy_3(x1,x2,x3,x6,x7, i)
      implicit double precision (a-h, o-z)
      include "./quadrature.i"
      r = rtri_ho(i)
      s = stri_ho(i)
      rs2xy_3 = ((1 - r - s) * x1) + (-(9 * r * r * s)
     .+ ((-9 * s * s + 9 * s + 2) * r) / 0.2d1) * x2
     .+ (-0.9d1 / 0.2d1 * (r * r) * s
     .+ ((-18 * s * s + 9 * s) * r) / 0.2d1 + s) * x3
     .+ (0.27d2 / 0.2d1 * (r * r) * s - 0.9d1 / 0.2d1 * r * s)
     .* x6 + ((27 * s * s - 9 * s) * r * x7) / 0.2d1

      return
      end function


      function dJ2(x1,x2,x3,x4,y1,y2,y3,y4, i)
      implicit double precision (a-h, o-z)
      include "./quadrature.i"
      r = rtri_ho(i)
      s = stri_ho(i)

      dJ2 = dabs( ((-x1 + x2) * (4.d0 * y4 - 2.d0 * y2 - 2.d0 * y3)
     .- (-y1 + y2) * (-2.d0 * x2 - 2.d0 * x3 + 4.d0 * x4)) * r
     .+ ((-2.d0 * x2 - 2.d0 * x3 + 4.d0 * x4) * (-y1 + y3)
     .- (4.d0 * y4 - 2.d0 * y2 - 2.d0 * y3) * (-x1 + x3)) * s
     .+ (-x1 + x2) * (-y1 + y3) - (-y1 + y2) * (-x1 + x3) )

!      temp = ((-x1 + x2) * (4.d0 * y4 - 2.d0 * y2 - 2.d0 * y3)
!     .- (-y1 + y2) * (-2.d0 * x2 - 2.d0 * x3 + 4.d0 * x4)) * r
!     .+ ((-2.d0 * x2 - 2.d0 * x3 + 4.d0 * x4) * (-y1 + y3)
!     .- (4.d0 * y4 - 2.d0 * y2 - 2.d0 * y3) * (-x1 + x3)) * s
!     .+ (-x1 + x2) * (-y1 + y3) - (-y1 + y2) * (-x1 + x3)
!      if(temp > 0.d0) print *, temp

      return
      end function

      function dJ3(x1,x2,x3,x6,x7,y1,y2,y3,y6,y7, i)
      implicit double precision (a-h, o-z)
      include "./quadrature.i"
      r = rtri_ho(i)
      s = stri_ho(i)

      dJ3 = dabs((-x1 + (- (18 * r * s)
     .- 0.9d1 / 0.2d1 *  s *  s + 0.9d1 / 0.2d1 *  s + 0.1d1) * x2
     .+ (- (9 * r * s) -  (9 * s * s) + 0.9d1 / 0.2d1 *  s) * x3
     .+ ( (27 * r * s) - 0.9d1 / 0.2d1 *  s) * x6
     .+  ((27 * s * s - 9 * s) * x7) / 0.2d1) * (-y1 + (- (9 * r * r)
     .+  ((-18 * s + 9) * r) / 0.2d1) * y2
     .+ (-0.9d1 / 0.2d1 *  r *  r
     .+  ((-36 * s + 9) * r) / 0.2d1 + 0.1d1) * y3
     .+ (0.27d2 / 0.2d1 *  r *  r - 0.9d1 / 0.2d1 *  r) * y6
     .+  ((54 * s - 9) * r * y7) / 0.2d1)
     .- (-y1 + (- (18 * r * s)
     .- 0.9d1 / 0.2d1 *  s *  s + 0.9d1 / 0.2d1 *  s + 0.1d1) * y2
     .+ (- (9 * r * s) -  (9 * s * s) + 0.9d1 / 0.2d1 *  s) * y3
     .+ ( (27 * r * s) - 0.9d1 / 0.2d1 *  s) * y6
     .+  ((27 * s * s - 9 * s) * y7) / 0.2d1) * (-x1 + (- (9 * r * r)
     .+  ((-18 * s + 9) * r) / 0.2d1) * x2
     .+ (-0.9d1 / 0.2d1 *  r *  r
     .+  ((-36 * s + 9) * r) / 0.2d1 + 0.1d1) * x3
     .+ (0.27d2 / 0.2d1 *  r *  r - 0.9d1 / 0.2d1 *  r) * x6
     .+  ((54 * s - 9) * r * x7) / 0.2d1) )




!      temp = (-x1 + (- (18 * r * s)
!     .- 0.9d1 / 0.2d1 *  s *  s + 0.9d1 / 0.2d1 *  s + 0.1d1) * x2
!     .+ (- (9 * r * s) -  (9 * s * s) + 0.9d1 / 0.2d1 *  s) * x3
!     .+ ( (27 * r * s) - 0.9d1 / 0.2d1 *  s) * x6
!     .+  ((27 * s * s - 9 * s) * x7) / 0.2d1) * (-y1 + (- (9 * r * r)
!     .+  ((-18 * s + 9) * r) / 0.2d1) * y2
!     .+ (-0.9d1 / 0.2d1 *  r *  r
!     .+  ((-36 * s + 9) * r) / 0.2d1 + 0.1d1) * y3
!     .+ (0.27d2 / 0.2d1 *  r *  r - 0.9d1 / 0.2d1 *  r) * y6
!     .+  ((54 * s - 9) * r * y7) / 0.2d1)
!     .- (-y1 + (- (18 * r * s)
!     .- 0.9d1 / 0.2d1 *  s *  s + 0.9d1 / 0.2d1 *  s + 0.1d1) * y2
!     .+ (- (9 * r * s) -  (9 * s * s) + 0.9d1 / 0.2d1 *  s) * y3
!     .+ ( (27 * r * s) - 0.9d1 / 0.2d1 *  s) * y6
!     .+  ((27 * s * s - 9 * s) * y7) / 0.2d1) * (-x1 + (- (9 * r * r)
!     .+  ((-18 * s + 9) * r) / 0.2d1) * x2
!     .+ (-0.9d1 / 0.2d1 *  r *  r
!     .+  ((-36 * s + 9) * r) / 0.2d1 + 0.1d1) * x3
!     .+ (0.27d2 / 0.2d1 *  r *  r - 0.9d1 / 0.2d1 *  r) * x6
!     .+  ((54 * s - 9) * r * x7) / 0.2d1)
!      return
      end function

      subroutine ho_boundary(mitot, mjtot, lstgrd, ssw, hx,hy,irr,sx,sy)
      implicit double precision (a-h, o-z)
      include "./quadrature.i"
      include "./cirr.i"

      dimension irr(mitot,mjtot)
      pi = 3.14159265358979d0


      ar_ho(lstgrd) = hx*hy
!     compute extra boundary points for ssv
      nbdry = 0
!      if(ihob .eq. 1) then
      do 17 i = 1, mitot
      do 17 j = 1, mjtot

       kirr = irr(i,j)
       if ((kirr .eq. -1) .or. (kirr .eq. lstgrd)) cycle

       ivert = 1
       do 213 while (poly(ivert+1,1,kirr) .ne. -11.)
         ivert = ivert + 1
  213  continue

      x1 = poly(ivert-1,1,kirr)
      y1 = poly(ivert-1,2,kirr)

      x2 = poly(ivert  ,1,kirr)
      y2 = poly(ivert  ,2,kirr)

      nbdry = 2
      if(ssw .eq. 2 .or. ssw .eq. -2) then
        nbdry = 3
      elseif(ssw .eq. 3 .or. ssw .eq. -3) then
        nbdry = 4
      endif
      bdry(1,1,kirr) = x1
      bdry(1,2,kirr) = y1

!      if(kirr .eq. 2 .or. kirr .eq. 50) then
!      print *, "here"
!      endif

      do 99 nb = 2,nbdry-1
      angle1 = atan(y1-sy,x1-sx)
      angle2 = atan(y2-sy,x2-sx)
      da = angle2-angle1

      if(da >  pi) da = da - 2.d0*pi
      if(da < -pi) da = da + 2.d0*pi

      curr_angle = angle1+da*(nb-1.d0)/(nbdry-1.d0)
      bdry(nb,1,kirr) = dsqrt((x1-sx)**2 + (y1-sy)**2) * cos(curr_angle)
     .+sx
      bdry(nb,2,kirr) = dsqrt((x1-sx)**2 + (y1-sy)**2) * sin(curr_angle)
     .+sy
  99  continue

      bdry(nbdry,1,kirr) = x2
      bdry(nbdry,2,kirr) = y2

!      if(ssw .eq. 2) then
!      bdry(2,1,kirr) = (bdry(1,1,kirr) + bdry(3,1,kirr))/2.d0
!      bdry(2,2,kirr) = (bdry(1,2,kirr) + bdry(3,2,kirr))/2.d0
!      endif

!      if(ssw .eq. -2) then
!      bdry(2,1,kirr) =  (2.d0/9.d0)*bdry(1,1,kirr)
!     .-(1/9.d0)*bdry(nbdry,1,kirr)
!     .+(8.d0/9.d0)*dsqrt(x1**2 + y1**2) * cos(angle1+da/2.d0)
!      bdry(3,1,kirr) = -(1.d0/9.d0)*bdry(1,1,kirr)
!     .+(2/9.d0)*bdry(nbdry,1,kirr)
!     .+(8.d0/9.d0)*dsqrt(x1**2 + y1**2) * cos(angle1+da/2.d0)
!
!      bdry(2,2,kirr) =  (2.d0/9.d0)*bdry(1,2,kirr)
!     .-(1/9.d0)*bdry(nbdry,2,kirr)
!     .+(8.d0/9.d0)*dsqrt(x1**2 + y1**2) * sin(angle1+da/2.d0)
!      bdry(3,2,kirr) = -(1.d0/9.d0)*bdry(1,2,kirr)
!     .+(2/9.d0)*bdry(nbdry,2,kirr)
!     .+(8.d0/9.d0)*dsqrt(x1**2 + y1**2) * sin(angle1+da/2.d0)
!      endif



      ! compute normals on the fly @ the quadrature points
      ! recompute areas


      ar_ho(kirr) = 0.d0
      xcirr_ho(kirr) = 0.d0
      ycirr_ho(kirr) = 0.d0
      do iv = 1, ivert-4
      x1 = poly(1,1,kirr)
      y1 = poly(1,2,kirr)

      x2 = poly(iv + 1,1,kirr)
      y2 = poly(iv + 1,2,kirr)

      x3 = poly(iv + 2,1,kirr)
      y3 = poly(iv + 2,2,kirr)

      xc = (x1+x2+x3)/3.d0
      yc = (y1+y2+y3)/3.d0
      da = area(x1,x2,x3,y1,y2,y3)

      ar_ho(kirr) = ar_ho(kirr) + da
      xcirr_ho(kirr) = xcirr_ho(kirr) + xc*da
      ycirr_ho(kirr) = ycirr_ho(kirr) + yc*da
      end do

      ! area of the second order/third order triangle

      if(ssw .eq. 2 .or. ssw .eq. -2) then
      x1 = poly(ivert-2,1,kirr)
      y1 = poly(ivert-2,2,kirr)

      x2 = bdry(1,1,kirr)
      y2 = bdry(1,2,kirr)

      x3 = bdry(3,1,kirr)
      y3 = bdry(3,2,kirr)

      x4 = bdry(2,1,kirr)
      y4 = bdry(2,2,kirr)


!      Aq = 0.d0
      do nq = 1,ntriquad_ho
      da = dJ2(x1,x2,x3,x4,y1,y2,y3,y4,nq)
      ar_ho(kirr) = ar_ho(kirr) + wtri_ho(nq)*da

      xc = rs2xy_2(x1,x2,x3,x4,nq)
      yc = rs2xy_2(y1,y2,y3,y4,nq)

      xcirr_ho(kirr) = xcirr_ho(kirr) + wtri_ho(nq)*xc*da
      ycirr_ho(kirr) = ycirr_ho(kirr) + wtri_ho(nq)*yc*da
      enddo


!      print *, ar_ho(kirr), ar(kirr)
!      print *, xcirr_ho(kirr), xcirr(kirr)
!      print *, ycirr_ho(kirr), ycirr(kirr)


      elseif(ssw .eq. 3 .or. ssw .eq. -3) then
      x1 = poly(ivert-2,1,kirr)
      y1 = poly(ivert-2,2,kirr)

      x2 = bdry(1,1,kirr)
      y2 = bdry(1,2,kirr)

      x3 = bdry(4,1,kirr)
      y3 = bdry(4,2,kirr)

      x4 = bdry(2,1,kirr)
      y4 = bdry(2,2,kirr)

      x5 = bdry(3,1,kirr)
      y5 = bdry(3,2,kirr)

!      ar_ho(kirr) = 0.d0
!      xcirr_ho(kirr) = 0.d0
!      ycirr_ho(kirr) = 0.d0

      do nq = 1,ntriquad_ho
      da = dJ3(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,nq)

      xc = rs2xy_3(x1,x2,x3,x4,x5,nq)
      yc = rs2xy_3(y1,y2,y3,y4,y5,nq)

      ar_ho(kirr)=ar_ho(kirr)+wtri_ho(nq)*da
      xcirr_ho(kirr) = xcirr_ho(kirr) + wtri_ho(nq)*xc*da
      ycirr_ho(kirr) = ycirr_ho(kirr) + wtri_ho(nq)*yc*da
      end do

      endif

      xcirr_ho(kirr) = xcirr_ho(kirr)/ar_ho(kirr)
      ycirr_ho(kirr) = ycirr_ho(kirr)/ar_ho(kirr)

  17  continue




!      endif





         poly_ho(8,1,lstgrd)  = 1.d0/12.d0
         poly_ho(9,1,lstgrd)  = 0.d0
         poly_ho(10,1,lstgrd) = 1.d0/12.d0
         dcubicshifts_ho(1,lstgrd) = 0.d0
         dcubicshifts_ho(2,lstgrd) = 0.d0
         dcubicshifts_ho(3,lstgrd) = 0.d0
         dcubicshifts_ho(4,lstgrd) = 0.d0

      ! recompute the moments
      do 171 i = 1, mitot
      do 171 j = 1, mjtot
       kirr = irr(i,j)
       if (kirr .eq. -1 .or. kirr .eq. lstgrd) go to 171

           shiftxx = 0.d0
           shiftxy = 0.d0
           shiftyy = 0.d0

           shiftxxx = 0.d0
           shiftxxy = 0.d0
           shiftxyy = 0.d0
           shiftyyy = 0.d0

           x0 = xcirr_ho(kirr)
           y0 = ycirr_ho(kirr)

          arr = ar_ho(kirr)
          ivert = 1
          do 123 while (poly(ivert+1,1,kirr) .ne. -11.)
            ivert = ivert + 1
  123     continue



          itri = ivert - 3
          idx1 = 1

          do 121 it = 1, itri-1 ! for each  triangle except ho one.
            idx2 = it + 1
            idx3 = it + 2

            x1 = poly(idx1,1,kirr)
            y1 = poly(idx1,2,kirr)

            x2 = poly(idx2,1,kirr)
            y2 = poly(idx2,2,kirr)

            x3 = poly(idx3,1,kirr)
            y3 = poly(idx3,2,kirr)

            artri = area(x1, x2, x3, y1, y2, y3)

            do 122 itq = 1,ntriquad

                xval = x1 * rtri(itq) + x2 * stri(itq)
     .              +  x3 * (1.d0-rtri(itq)-stri(itq))
                yval = y1 * rtri(itq) + y2 * stri(itq)
     .              +  y3 * (1.d0-rtri(itq)-stri(itq))

      shiftxx = shiftxx + (artri/arr)*wtri(itq)*(xval-x0)**2 /(hx**2)
      shiftxy = shiftxy
     .          + (artri/arr)*wtri(itq)*(xval-x0)*(yval-y0)/(hx * hy)
      shiftyy = shiftyy + (artri/arr)*wtri(itq)*(yval-y0)**2 /(hy**2)


      shiftxxx = shiftxxx + (artri/arr)*wtri(itq)*(xval-x0)**3 /(hx**3)
      shiftxxy = shiftxxy +
     .          (artri/arr)*wtri(itq)*(yval-y0)*(xval-x0)**2 /(hy*hx**2)
      shiftxyy = shiftxyy +
     .          (artri/arr)*wtri(itq)*(xval-x0)*(yval-y0)**2 /(hx*hy**2)
      shiftyyy =shiftyyy + (artri/arr)*wtri(itq)*(yval-y0)**3 /(hy**3)



  122        continue ! for each quadrature point on each triangle
  121      continue ! for each triangle



      if(ssw .eq. 2  .or. ssw .eq. -2) then
      x1 = poly(ivert-2,1,kirr)
      y1 = poly(ivert-2,2,kirr)

      x2 = bdry(1,1,kirr)
      y2 = bdry(1,2,kirr)

      x3 = bdry(3,1,kirr)
      y3 = bdry(3,2,kirr)

      x4 = bdry(2,1,kirr)
      y4 = bdry(2,2,kirr)

      do nq = 1,ntriquad_ho
      artri = dJ2(x1,x2,x3,x4,y1,y2,y3,y4,nq)

      xval = rs2xy_2(x1,x2,x3,x4,nq)
      yval = rs2xy_2(y1,y2,y3,y4,nq)

      shiftxx = shiftxx
     .          + (artri/arr)*wtri_ho(nq)*(xval-x0)**2 /(hx**2)
      shiftxy = shiftxy
     .          + (artri/arr)*wtri_ho(nq)*(xval-x0)*(yval-y0)/(hx * hy)
      shiftyy = shiftyy
     .          + (artri/arr)*wtri_ho(nq)*(yval-y0)**2 /(hy**2)


      shiftxxx = shiftxxx +
     .     (artri/arr)*wtri_ho(nq)*(xval-x0)**3 /(hx**3)
      shiftxxy = shiftxxy +
     .     (artri/arr)*wtri_ho(nq)*(yval-y0)*(xval-x0)**2 /(hy*hx**2)
      shiftxyy = shiftxyy +
     .     (artri/arr)*wtri_ho(nq)*(xval-x0)*(yval-y0)**2 /(hx*hy**2)
      shiftyyy =shiftyyy
     .   + (artri/arr)*wtri_ho(nq)*(yval-y0)**3 /(hy**3)

      enddo
      elseif(ssw .eq. 3 .or. ssw .eq. -3) then
      x1 = poly(ivert-2,1,kirr)
      y1 = poly(ivert-2,2,kirr)

      x2 = bdry(1,1,kirr)
      y2 = bdry(1,2,kirr)

      x3 = bdry(4,1,kirr)
      y3 = bdry(4,2,kirr)

      x4 = bdry(2,1,kirr)
      y4 = bdry(2,2,kirr)

      x5 = bdry(3,1,kirr)
      y5 = bdry(3,2,kirr)

      do nq = 1,ntriquad_ho
      artri = dJ3(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,nq)

      xval = rs2xy_3(x1,x2,x3,x4,x5,nq)
      yval = rs2xy_3(y1,y2,y3,y4,y5,nq)

      shiftxx = shiftxx
     .       + (artri/arr)*wtri_ho(nq)*(xval-x0)**2 /(hx**2)
      shiftxy = shiftxy
     .       + (artri/arr)*wtri_ho(nq)*(xval-x0)*(yval-y0)/(hx * hy)
      shiftyy = shiftyy
     .       + (artri/arr)*wtri_ho(nq)*(yval-y0)**2 /(hy**2)


      shiftxxx = shiftxxx
     .     + (artri/arr)*wtri_ho(nq)*(xval-x0)**3 /(hx**3)
      shiftxxy = shiftxxy
     .     + (artri/arr)*wtri_ho(nq)*(yval-y0)*(xval-x0)**2 /(hy*hx**2)
      shiftxyy = shiftxyy
     .     + (artri/arr)*wtri_ho(nq)*(xval-x0)*(yval-y0)**2 /(hx*hy**2)
      shiftyyy =shiftyyy
     .     + (artri/arr)*wtri_ho(nq)*(yval-y0)**3 /(hy**3)

      end do

      endif







           poly_ho(8,1,kirr)  = shiftxx
           poly_ho(9,1,kirr)  = shiftxy
           poly_ho(10,1,kirr) = shiftyy
           dcubicshifts_ho(1,kirr) = shiftxxx
           dcubicshifts_ho(2,kirr) = shiftxxy
           dcubicshifts_ho(3,kirr) = shiftxyy
           dcubicshifts_ho(4,kirr) = shiftyyy


 171      continue
      end subroutine


       subroutine getCellCentroid_lo(lstgrd,i,j,xc,yc,xlow,ylow,dx,dy,k)

      implicit double precision (a-h,o-z)
      include "cirr.i"

      if (k .eq. lstgrd) then
         xc = xlow + (i-0.5d0)*dx
         yc = ylow + (j-0.5d0)*dy
      else if (k .eq. -1) then
         xc = 0.d0
         yc = 0.d0
      else
           xc = xcirr(k)
           yc = ycirr(k)
      endif

      return
      end

