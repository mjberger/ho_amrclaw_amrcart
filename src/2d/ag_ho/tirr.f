c
c ---------------------------------------------------------------
c
       subroutine tirr(irr,nrow,ncol,xl,xr,yb,yt,level)
c
c set the portion of the irregular array at the given level
c bounded by the given coordinates
c
      implicit double precision (a-h,o-z)
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
      dimension ialloc(allocsize)
      equivalence (alloc, ialloc)
c
       logical   xint,yint
       dimension irr(nrow,ncol)
c
       iirr(i,j) = 2*(locirr-1)+1+i-1+mitot*(j-1)
c
       hx = hxposs(level)
       hy = hyposs(level)
c
       mptr = lstart(level)
c
c check if mptr and patch intersect. recall that irr array
c is enlarged by lwidth around grid, so use effective size.
c
 10    maxi   = node(5,mptr)
       maxj   = node(6,mptr)
       maxip1 = maxi + 1
       maxjp1 = maxj + 1
       locirr = node(14,mptr)
       mitot  = maxi-1+2*lwidth
       mjtot  = maxj-1+2*lwidth
       xst    = rnode(1,mptr)-lwidth*hx
       yst    = rnode(2,mptr)-lwidth*hy
       xend   = rnode(7,mptr)+lwidth*hx
       yend   = rnode(4,mptr)+lwidth*hy
       xint   = .false.
       yint   = .false.

      if (((xl .le. xst) .and. (xr .gt. xst)) .or.
     1    ((xl .ge. xst) .and. (xl .lt. xend)))   xint = .true.
      if (((yb .le. yst) .and. (yt .gt. yst)) .or.
     1    ((yb .ge. yst) .and. (yb .lt. yend)))   yint = .true.
      if (.not. (xint .and. yint)) go to 90
c
c  calculate starting and ending indices for copying values from
c  mptr to patch
c
      ist  = max(1, idint((xst-xl)/hx + 1.1))
      jst  = max(1, idint((yst-yb)/hy + 1.1))
      iend = min(nrow,idint( (xend-xl)/hx +  .1))
      jend = min(ncol,idint( (yend-yb)/hy +  .1))
c     play it safe
      if ((ist .gt. iend) .or. (jst .gt. jend)) go to 90
c
c  calculate starting index for mptr as source updater
c
      isrc = max(1, idint((xl-xst)/hx + 1.1))
      jsrc = max(1, idint((yb-yst)/hy + 1.1)) 
c
      jdonor = jsrc
      do 55 j = jst, jend
         idonor = isrc
         do 45 i = ist, iend
            irr(i,j) = ialloc(iirr(idonor,jdonor))
            idonor = idonor + 1
 45      continue
         jdonor = jdonor + 1
 55    continue
c
 90    mptr = node(10, mptr)
       if (mptr .ne. 0) go to 10
c
       return
       end
