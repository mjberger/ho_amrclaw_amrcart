c
c -----------------------------------------------------------
c
      subroutine regrid  (nvar,lbase,cut,cdist,work,steady,quad)
c
      implicit double precision (a-h,o-z)
      dimension          work(nvar)
      logical            steady, quad
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
c
c ******************************************************************
c  regrid = flag points on each grid with a level > = lbase.
c  cluster them, and fit new subgrids around the clusters.
c  the lbase grids stay fixed during regridding operation.
c  when a parent grid has its error estimated, add its kid grid
c  information to the error grid before clustering. (project)
c  order of grid examination - all grids at the same level, then
c  do the next coarser level.
c
c input parameters:
c     lbase  = highest level that stays fixed during regridding
c     cutoff = criteria for measuring goodness of rect. fit.
c local variables:
c     lcheck = the level being examined.
c     lfnew  = finest grid to be. will replace lfine.
c global
c    mstart  = start of very coarsest grids.
c *****************************************************************
c
      lcheck    = min0(lfine,mxnest-1)
      lfnew     = lbase
      do 10 i   = 1, mxnest
 10   newstl(i) = 0
      time      = rnode(12, lstart(lbase))
c
 20   if (lcheck .lt. lbase) go to 50
         call grdfit(lbase,lcheck,nvar,cut,cdist,lfnew,time,
     .               steady)
          if (newstl(lcheck+1) .eq. 0) go to 40
          lfnew = max0(lcheck + 1,lfnew)
 40       continue
          lcheck = lcheck - 1
c
      go to 20
 50   continue
c
c  end of level loop
c
c  remaining tasks left in regridding:
c
c  interpolate storage for the new grids.  the starting pointers
c  for each level are in newstl. also reclaim some space before new
c  allocations.
      call fill(lbase, lfnew, nvar, work, quad)
c
c  merge data structures (newstl and lstart )
c  finish storage allocation, reclaim space, etc. set up boundary
c  flux conservation arrays
c
      call join (lbase, nvar)
      if (lbase .eq. lfine) go to 70
      do 60 level = lbase, lfine-1
      call prepf(level+1,nvar)
      call prepc(level,nvar)
 60   continue
c
 70   continue
c
      return
      end
