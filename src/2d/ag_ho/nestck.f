c
c ---------------------------------------------------------
c
      logical function nestck(mnew,lbase,badpts,npts,
     1                        numptc,icl,nclust,hxbase,hybase)
c
      implicit double precision (a-h,o-z)
      logical   phys, inbase, ylong
      dimension      badpts(2,npts)
      integer   numptc(maxcl)
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      logical    nprint
      data       nprint/.false./
c
c ******************************************************************
c
c nestck - check that the potential grid mnew is completely
c          contained in the (coarser) finest grid which stays
c          fixed, at level lbase. projec algo. will guarantee
c          containment in all finer grids twixt them.
c          if grid not contained in some coarse grid,  then
c          bisect in long direction.
c          EVENTUALLY this has to work.
c
c input parameter:
c   mnew    - grid descriptor of potential grid
c   lbase   - level which stays fixed during regridding
c   badpts  - only the flagged pts. in this cluster (# icl)
c ******************************************************************
c
c  go around all sides checking for point containment.skip if a
c  physical boundary pt. it is possible that first and last point
c  on a side are physical, and not others. order of sides are
c     |-2-|
c     1   3
c     |-4-|
c
       nestck = .true.
c
c  with prebuffering, grid can not exceed coarsest grid (since
c  assuming convex domain), so don't bother to nestck
c
       if (lbase .eq. 1)  go to 99
       levnew = node(4,mnew)
       mi     = max(idint((rnode(7,mnew)-rnode(1,mnew))/hxbase+.001),1)
       mj     = max(idint((rnode(4,mnew)-rnode(2,mnew))/hybase+.001),1)
c  side 1
       if (phys(rnode(1,mnew),
     1         (rnode(2,mnew)+rnode(4,mnew))/2.)) go to 15
c
       x = rnode(1,mnew) - hxbase / 4.
       y = rnode(2,mnew) + hybase / 4.
       do 10 j = 1, mj
       if (inbase(x,y,lbase)) go to 10
       if (phys(x,y)) go to 10
       go to 50
 10    y = y + hybase
c  side 3
 15    if (phys(rnode(7,mnew),
     1         (rnode(6,mnew)+rnode(8,mnew))/2.)) go to 25
c
       x = rnode(7,mnew) + hxbase / 4.
       y = rnode(6,mnew) - hybase / 4.
       do 20 j = 1, mj
       if (inbase(x,y,lbase)) go to 20
       if (phys(x,y)) go to 20
       go to 50
 20    y = y - hybase
c  side 2
 25    if (phys((rnode(3,mnew)+rnode(5,mnew))/2.,
     1           rnode(4,mnew))) go to 35
c
       x = rnode(3,mnew) + hxbase / 4.
       y = rnode(4,mnew) + hybase / 4.
       do 30 i = 1, mi
       if (inbase(x,y,lbase)) go to 30
       if (phys(x,y)) go to 30
       go to 50
 30    x = x + hxbase
c  side 4
 35    if (phys((rnode(7,mnew)+rnode(1,mnew))/2.,
     1           rnode(8,mnew))) go to 99
c
       x = rnode(7,mnew) - hxbase / 4.
       y = rnode(8,mnew) - hybase / 4.
       do 40 i = 1, mi
       if (inbase(x,y,lbase)) go to 40
       if (phys(x,y)) go to 40
       go to 50
 40    x = x - hxbase
c
       go to 99
c
c  grid not properly nested. bisect in long direction, and return
c  two clusters instead of 1.
c
 50    if (npts .gt. 1) go to 55
           write(6,101) levnew
 101       format(' nestck: 1 pt. cluster at level ',i5,' still not',
     1       ' nested',/,'          pt. too close to boundary')
           call outtre(mstart, .false.,nvar)
           call outmsh(mnew, .false.,nvar)
           stop
 55    if (nclust .lt. maxcl) go to 60
           write(6,102) maxcl
 102       format(' too many clusters: > ',i5,' (from nestck) ')
           stop
 60   if (nprint) write(6,103) icl, npts
 103  format(' bisecting cluster ',i5,' with ',i5,' pts. in nestck')
      if (rnode(7,mnew)-rnode(1,mnew) .gt.
     1     rnode(4,mnew) - rnode(2,mnew)) go to 70
       ymid  = (rnode(4,mnew) + rnode(2,mnew) ) / 2.
       ylong = .true.
       go to 80
 70    xmid  = (rnode(7,mnew) + rnode(1,mnew) ) / 2.
       ylong = .false.
c
 80    ipt  = 1
       ntop = npts
 90    if (ylong) then
          if (badpts(2,ipt) .lt. ymid) go to 100
        else
          if (badpts(1,ipt) .lt. xmid) go to 100
       endif
c
c  swap with a point in top half not yet tested. keep smaller
c  half of rect. in bottom half
c
       temp           = badpts(1,ipt)
       badpts(1,ipt)  = badpts(1,ntop)
       badpts(1,ntop) = temp
       temp           = badpts(2,ipt)
       badpts(2,ipt)  = badpts(2,ntop)
       badpts(2,ntop) = temp
       ntop           = ntop - 1
       if (ipt .le. ntop) go to 90
       go to 110
 100   ipt = ipt +1
       if (ipt .le. ntop) go to 90
c
c  ntop points to top of 1st cluster (= no. of points in 1st cluster)
c
 110   numptc(icl)     = npts - ntop
       do 120 i        = icl, nclust
       nmove           = nclust + icl - i
 120   numptc(nmove+1) = numptc(nmove)
       numptc(icl)     = ntop
       nclust          = nclust + 1
       nestck          = .false.
c
 99    return
       end
