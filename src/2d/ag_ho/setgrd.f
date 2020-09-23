c
c -----------------------------------------------------------------
c
      subroutine setgrd (nvar,cut,cdist,lfix,quad)
c
      implicit double precision (a-h,o-z)
      logical            graf,quad
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
      common   /stats/ evol, rvol, rvoll(maxlv), lentot, lenmax
      logical  vtime, steady
      data     vtime/.false./, steady/.false./
c     # may as well not bother to calculate time step for error est.
c
c  set up the entire tree/grid structure.  only at this time t = 0
c  can we take advantage of initialization routines.
c  remember that regridding/error estimation needs to have two
c  time steps of soln. values.
c
      if (lfix .eq. mxnest) go to 60
c
c set up fixed levels so can continue error est. to make new levels
c
      do 5 lcheck = 1, lfix-1
         call advanc(lcheck,nvar,dtlev,vtime,steady)
         evol = evol + rvol
         rvol = 0.d0
 5    continue
c
      levnew = lfix+1
      lbase  = lfix
      time   = 0.0d0
      mxinit = 6
c
 10   if (levnew .gt. mxnest) go to 30
c10   if (levnew .gt. mxinit) go to 30
          levold = levnew - 1
          if (lstart(levold) .eq. 0) go to 30
          lbase  = levold
          lfnew  = lbase
c
c  set up level to be flagged
c
         call advanc(levold,nvar,dtlev,vtime,steady)
         evol = evol + rvol
         rvol = 0.d0
 
         do 20 level=1,mxnest
 20         rvoll(level) = 0.d0
c
c  flag, cluster, and make new grids
c
         call grdfit(lbase,levold,nvar,cut,cdist,lfnew,time,
     .               steady)
         if (newstl(levnew) .ne. 0) lfnew = levnew
c
c  init new level. after each iteration. fix the data structure
c  also reinitalize coarser grids so fine grids can be advanced
c  and interpolate correctly for their bndry vals from coarser grids.
c
         call ginit(newstl(levnew),.true., nvar,quad)
         lstart(levnew) = newstl(levnew)
         lfine = lfnew
         call join(levold, nvar)
         call ginit(lstart(levold),.false., nvar,quad)
c
         levnew = levnew + 1
      go to 10
 30   continue
c
c  switch location of 7 and 8 vals, and reset time to 0.0
c
      if (mxnest .eq. 1) go to 99
c
      lev = 1
 40   if ((lev .eq. lfine) .or. (lev .gt. lfine))  go to 60
        mptr = lstart(lev)
 50        temp           = node(7,mptr)
           node(7,mptr)   = node(8,mptr)
           node(8,mptr)   = temp
           rnode(12,mptr) = 0.0d0
           mptr = node(10,mptr)
           if (mptr .ne. 0) go to 50
       lev = lev + 1
       go to 40
 60   continue
c set up boundary flux conservation arrays
      if (lfine .eq. 1) go to 99
      do 70 level = 1, lfine-1
      call prepf(level+1,nvar)
      call prepc(level,nvar)
 70   continue
c
 99   continue
c
      return
      end
