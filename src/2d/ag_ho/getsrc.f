c
c -----------------------------------------------------------
c
      subroutine getsrc( x,y, levrec, gptr, srcptr, srclev)
c
      implicit double precision (a-h,o-z)
      integer  muse, srcptr, srclev, levrec, gptr
      logical  cntain
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
c
c *****************************************************************
c getsrc = find the best (finest) source grid to use for
c          interpolating and updating the point (x,y).
c
c input parameters:
c     x,y   = coordinates of receiving point.
c   levrec  =  level of receiving pt = max. source level to check.
c              doesn't pay to look for the absolutely finest source
c              grid, if (x,y) at a coarser level.
c     gptr  = initial guess for the best source grid
c output parameters:
c     srcptr    = best source grid
c     srclevel  = level of the aforementioned source
c
c local variables:
c     muse = contains a possible source grid not nec. the best.
c     potentialptr, trykid = possible sources
c note: if source grid cannot be determined from gptr, then the
c       entire tree is searched (using the global mstart)
c *****************************************************************
c
      if (gptr .eq. 0) gptr = mstart
      if (.not. (cntain(x,y,gptr))) go to 15
         muse = gptr
         go to 40
 15   continue
c
c   find coarser grid containing pt. to start
c
       muse = mstart
 20    if (cntain(x,y,muse)) go to 40
           muse = node(10,muse)
           if (muse .ne. 0) go to 20
                write(6,100)  x,y
100             format(' pt.',2e12.5,' not contained in a coarse grid')
                call outlev(levrec,nvar)
                stop
 40    continue
c
c muse is a source grid.  check for a finer grid for better source.
c
      level = node(4, muse)  + 1
 50   if (level .gt. levrec) go to 90
          mptr = lstart(level)
 60       if (mptr .eq. 0) go to 80
              if (.not.(cntain(x,y,mptr))) go to 70
                  muse = mptr
                  go to 80
 70           mptr = node(10, mptr)
             go to 60
 80       level = level + 1
          go to 50
c
c
 90   srcptr     =   muse
      srclev     =   node(4, srcptr)
c
      return
      end
