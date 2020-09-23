c
c -------------------------------------------------------------
c
      subroutine flgbuf(rbuff,iflow,numbuf,level)
c
      implicit double precision (a-h,o-z)
      logical            graf, phys, cntain
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
      dimension rbuff(3,iflow)
c
      iadd(i,j) = locerr + i - 1 + maxip1*(j-1)
c
c  every point in the rbuff array should either be outside the
c  physical domain, or have a neighboring grid to flag.
c  otherwise, flagged points have come too close to the boundary
c
      nbad = 0
      do 30 ifind = 1, iflow
        xcent = rbuff(1,ifind)
        ycent = rbuff(2,ifind)
        mptr  = idint(rbuff(3,ifind))
        mtry  = lstart(level)
c
 10     if (mptr .eq. mtry) go to 20
        if (.not. cntain(xcent,ycent,mtry)) go to 20
c
c pt. (xcent,ycent) is in grid mtry. flag mtry's error plane if
c not already flagged
c
        locerr = node(8,mtry)
        maxip1 = node(5,mtry) + 1
        maxjp1 = node(6,mtry) + 1
        hx = rnode(9,  mtry)
        hy = rnode(10, mtry)
        i = idint((xcent-rnode(1,mtry))/hx + 1.51)
        j = idint((ycent-rnode(2,mtry))/hy + 1.51)
        if (alloc(iadd(i,j)) .eq. 0) then
            alloc(iadd(i,j)) = 1.0
            numbuf = numbuf + 1
        endif
        go to 30
c
 20     mtry = node(10,mtry)
        if (mtry .ne. 0) go to 10
c
c  see if point is outside domain
c
       if (phys(xcent,ycent)) go to 30
c
c  output error msg. - no home found for this point
c
c      write(6,100) xcent,ycent,mptr
 100   format(' bad news: no grid to flag pt. ',2e14.7,' from grid ',i4)
       nbad = nbad + 1
c
 30    continue
c
       if (nbad .gt. 0) write(6,101) nbad,mptr
 101   format(' bad news: ',i4,' buffer pts. ignored for grid ',i4)
c
       return
       end
