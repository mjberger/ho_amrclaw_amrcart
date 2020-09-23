c
c ------------------------------------------------------
c
       subroutine qsintfil(xl,xr,yb,yt,level,nrowst,nrowend,
     1                    ncolst,ncolend,nvar,valbig,mitot,mjtot,
     2                   time,intime,ifill,locuse,qx,qy)
c
c  INTF: interpolates values for a patch at the specified level and
c  location, using values from grids at LEVEL and coarser, if nec.
c
c  take the intersection of a grid patch with corners at XL,XR,YB,YT
c  and all grids mptr at LEVEL.  If there is a non-null intersection
c  copy the solution vaues from mptr (at TIME) into VAL array.
c  assumes patch at same level so do straight copy, not skipping
c  every intrat or doing any interpolation here,
c  assume called in correct order of levels, so that when copying
c  is ok to overwrite. there are no dummy points around patch, since
c  this is not an official grid.
c
c ifill =  0 : calling program is fill for coarser level.
c              in this case, all values should be set by
c              interp. from grids at same level. if used array has
c              any remaining zeroes it is an error. this is because
c              grid patch being interpolated will be a new finer
c              level, and proper nesting should be in force.
c          1:  calling program is fill for finer level. ignore
c              any unfilled values and do not do consistency check.
c
c  used array marks when point filled. if points left over after
c  rectangle intersections at specified level, for option 3,use
c  getsrc/incent to finish it off.  switch between treating used
c  array as 2d or 1d for cray vectorization purposes.
c
      implicit double precision (a-h,o-z)
      parameter  (msize = 264)
      parameter  (maxgr = 192, maxlv=12)
      logical    graf
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
c
      logical     xint, yint
      dimension   valbig(mitot,mjtot,nvar), intinc(maxlv)
      dimension   qx(msize,msize,nvar), qy(msize,msize,nvar)
c
      iadd(i,j,ivar)   = loc    + i - 1 + maxip1*((ivar-1)*maxjp1+j-1)
      iadnew(i,j,ivar) = locnew + i - 1 + maxip1*((ivar-1)*maxjp1+j-1)
      iadold(i,j,ivar) = locold + i - 1 + maxip1*((ivar-1)*maxjp1+j-1)
      iaduse(i,j)      = locuse + i - 1 + nrow*(j-1)

      iaddsx(i,j,ivar)   = locqx  + i - 1 + misrc*((ivar-1)*mjsrc+j-1)
      iaddsy(i,j,ivar)   = locqy  + i - 1 + misrc*((ivar-1)*mjsrc+j-1)
c
      if (level.eq.0) return
 
      dt     = possk(level)
      hx     = hxposs(level)
      hy     = hyposs(level)
      dt10   = dt / 10.
      intinc(mxnest) = 1
      do 3 i = 1, mxnest-1
         ii = mxnest - i
 3       intinc(ii) = intinc(ii+1)*intrat(ii)
      nrow = nrowend-nrowst+1
      ncol = ncolend-ncolst+1
      ntot   = nrow *ncol
      do 5 i = 1, ntot
5        alloc(locuse+i-1) = 0.0
      mptr   = lstart(level)
      dt = rnode(11,mptr)
c
 10   if (mptr .eq. 0) go to 105
c
c  check if grid mptr and patch intersect
c
      maxi   = node(5,mptr)
      maxj   = node(6,mptr)
      maxip1 = maxi + 1
      maxjp1 = maxj + 1
      xst    = rnode(1,mptr)
      yst    = rnode(2,mptr)
      xend   = rnode(7,mptr)
      yend   = rnode(4,mptr)
      misrc  = maxi-1+2*lwidth
      mjsrc  = maxj-1+2*lwidth
c
c  there are no tolerances in the following arithmetic - the
c  borderline cases should all work out ok, and there is a check
c  later in case something snuck through the cracks (where xr = xst)
c
      xint = .false.
      yint = .false.
      if (((xl .le. xst) .and. (xr .gt. xst)) .or.
     1    ((xl .ge. xst) .and. (xl .lt. xend)))   xint = .true.
      if (((yb .le. yst) .and. (yt .gt. yst)) .or.
     1    ((yb .ge. yst) .and. (yb .lt. yend)))   yint = .true.
      if (.not. (xint .and. yint)) go to 90
c
c  calculate starting and ending indices for copying values from
c  mptr to patch
c
      ist  = max(nrowst, nrowst+idint((xst-xl)/hx + .1))
      jst  = max(ncolst, ncolst+idint((yst-yb)/hy + .1))
      iend = min(nrowend,nrowst-1+idint( (xend-xl)/hx +  .1))
      jend = min(ncolend,ncolst-1+idint( (yend-yb)/hy +  .1))
c     play it safe
      if ((ist .gt. iend) .or. (jst .gt. jend)) go to 90
c
c  calculate starting index for mptr as source updater
c
      isrc = max(1, idint((xl-xst)/hx + 1.1)) + 1
      jsrc = max(1, idint((yb-yst)/hy + 1.1)) + 1
c
c  figure out at which time solution values are wanted.
c  if time interpolation needed, use separate loop, for speed
c  alphai = 1 for new time; 0 for old time
c
      alphac = dfloat(intcnt(level)-intime)
     .          / dfloat(intinc(level))
      alphai  = 1.d0-alphac
c     alphai = (time - rnode(12, mptr) + dt) / dt
      if ((alphai .lt. -dt10) .or. (alphai .gt. 1.+dt10)) then
          write(6,900) time, mptr, level
 900      format(' time wanted ',e15.7,' not available from grid ',i4,
     1           'level',i4)
          write(6,901)xl,xr,yb,yt,mptr,level,time,rnode(12,mptr),
     .                alphai,dt10,ifill
          call outtre(mstart,.false.,nvar)
          stop
      endif
c
      if (dabs(alphai - 1.) .lt. dt10) then
          loc = node(7,mptr)
	  locqx = node(3,mptr)
	  locqy = locqx + 4*misrc*mjsrc
      else if (dabs(alphai) .lt. dt10) then
          loc = node(8,mptr)
	  locqx = node(3,mptr)
	  locqy = locqx + 4*misrc*mjsrc
          if (level .eq. mxnest) then
              write(6,901) xl,xr,yb,yt,mptr,level,time,rnode(12,mptr),
     .                     alphai,dt10,ifill
              stop
           endif
      else
          locold = node(8,mptr)
          locnew = node(7,mptr)
          loc    = node(7,mptr)
	  locqx = node(3,mptr)
	  locqy = locqx + 4*misrc*mjsrc
          if (level .eq. mxnest) then
             write(6,901) xl,xr,yb,yt,mptr,level,time,rnode(12,mptr),
     .                    alphai,dt10,ifill
             stop
          endif
c         go to 60
      endif
 901  format(' trying to interpolate from previous time values ',/,
     .       ' for a patch with corners xl,xr,yb,yt:'
     .       ,/,2x,4e15.7,/,
     .       ' from source grid ',i4,' at level ',i4,/,
     .       ' time wanted ',e15.7,' source time is ',e15.7,/,
     .       ' alphai, dt10, ifill ',2e15.7,i4)
c
c  copy the solution values.  no time interp. in this loop version.
c
      do 45 ivar = 1, nvar
      jdonor  = jsrc
          do 35 j = jst, jend
          idonor  = isrc
              do 20 i = ist, iend
              idonor  = (isrc + i - ist)
              valbig(i,j,ivar) = alloc(iadd(idonor, jdonor, ivar))
                  qx(i+lwidth-1,j+lwidth-1,ivar) = 
     &                      alloc(iaddsx(idonor+lwidth-1, 
     &                                   jdonor+lwidth-1, ivar))
                  qy(i+lwidth-1,j+lwidth-1,ivar) = 
     &                      alloc(iaddsy(idonor+lwidth-1, 
     &                                   jdonor+lwidth-1, ivar))
              alloc(iaduse(i-nrowst+1,j-ncolst+1)) = 1.0
 20           continue
          jdonor = jdonor + 1
 35       continue
 45   continue
      go to 90
c
c  time interpolation version of the loop
c
 60   do 85 ivar = 1, nvar
      jdonor  = jsrc
          do 75 j = jst, jend
          idonor  = isrc
              do 65 i = ist, iend
              idonor  = (isrc + i - ist)
              valbig(i,j,ivar) = alloc(iadnew(idonor,jdonor,ivar))*
     1              alphai + alphac*alloc(iadold(idonor,jdonor,ivar))
              alloc(iaduse(i-nrowst+1,j-ncolst+1)) = 1.0
 65           continue
          jdonor = jdonor + 1
 75       continue
 85   continue
c
 90   mptr = node(10, mptr)
      go to 10
c
 105  continue
 
c  set used array points which intersect boundary to be equal to 1;
c  they will be set in physbd.
 
      hxmarg = hx*.01
      hymarg = hy*.01
      if (yt.gt.yprob + hymarg) then
        nyt = (yt - yprob + hymarg)/hy
        jbeg = max0(ncolend - nyt + 1,ncolst)
        do 1000 j= jbeg, ncolend
        do 1000 i= nrowst,nrowend
           alloc(iaduse(i-nrowst+1,j-ncolst+1)) = 1.
1000    continue
      endif
 
      if (yb.lt.-hymarg) then
        nyb = (hymarg - yb)/hy
        nyb = min0(nyb,ncolend)
        do 1200 j= ncolst, nyb
        do 1200 i= nrowst,nrowend
            alloc(iaduse(i-nrowst+1,j-ncolst+1)) = 1.
1200    continue
      endif
 
      if (xl.lt.-hxmarg) then
        nxl = (hxmarg - xl)/hx
        nxl = min0(nxl,nrowend)
        do 1400 i= nrowst,nxl
        do 1400 j= ncolst,ncolend
           alloc(iaduse(i-nrowst+1,j-ncolst+1)) = 1.
1400    continue
      endif
 
      if (xr.gt.xprob + hxmarg) then
        nxr = (xr - xprob + hxmarg)/hx
        ibeg = max0(nrowst,nrowend - nxr + 1)
        do 1600 i= ibeg,nrowend
        do 1600 j= ncolst,ncolend
           alloc(iaduse(i-nrowst+1,j-ncolst+1)) = 1.
1600    continue
      endif
c
c  check if any unfilled cells = any 0.'s left in used array.
c we ignore corner zones if ifill = 0, since intfil is then being called
c from filval, and the corner zones are not used in the
c interpolation of new values.
c
      if (ifill.eq.0) then
         alloc(iaduse(1,1))     = 1.
         alloc(iaduse(nrow,1))  = 1.
         alloc(iaduse(1,ncol))  = 1.
         alloc(iaduse(nrow,ncol)) = 1.
      endif
c 
      usemin  = alloc(iaduse(1,1))
      do 110 i = 1, ntot
 110     usemin  = dmin1(usemin, alloc(locuse+i-1))
      if (usemin .eq. 1.) go to 99
c
      if (ifill .eq. 0) then
         write(6,902) ifill,level,xl,xr,yb,yt,
     1                nrowst,nrowend,ncolst,ncolend
 902     format(' intfil has points not filled with ifill ',i4,/,1x,
     1          ' trying to fill patch at level ',i5,/,1x,
     2          ' with corners: ',4e15.7,/,1x,
     2          ' nrowst,nrowend,ncolst,ncolend ',4i3)
         stop
      endif
c
99    continue
      return
      end
