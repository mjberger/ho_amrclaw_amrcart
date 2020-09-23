c
c  -----------------------------------------------------------
c
      subroutine fill (lbase, lfnew, nvar, work, quad)
c
      implicit double precision (a-h,o-z)
      dimension  work(nvar)
      logical            graf, quad
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
      logical     fprint
      data        fprint /.false./
c
c fill = interpolate initial values for the newly created grids.
c        the start of each level is located in newstl array.
c        since only levels greater than lbase were examined, start
c        looking there.
c
c   reclaim old storage (position 8) and list space 15 and 16
c   before allocating new storage. remember, finest level grids
c  (if level = mxnest so that error never estimated) don't have
c  2 copies of solution values at old and new times.
c
c
      call putsp(lbase,lbase,nvar)
      level = lbase + 1
 1    if (level .gt. lfine) go to 4
      call putsp(lbase,level,nvar)
          mptr = lstart(level)
 2        if (mptr .eq. 0) go to 3
              maxi   = node(5, mptr)
              maxj   = node(6, mptr)
              maxip1 = maxi + 1
              maxjp1 = maxj + 1
              nwords        = maxip1*maxjp1*nvar
              if (level .lt. mxnest) call reclam(node(8, mptr), nwords)
              node(8, mptr) = 0
              maxipb = maxi + lwidth - 1
              maxjpb = maxj + lwidth - 1
c             call reclam(node(2,mptr), 2*nvar*lwidth*(maxipb+maxjpb))
              mptr          = node(10, mptr)
          go to 2
 3        level   = level + 1
          go to 1
c
 4    lcheck = lbase + 1
 5    if (lcheck .gt. maxlv) go to 99
c
c  interpolate level lcheck
c
          mptr   = newstl(lcheck)
 10       if (mptr .eq. 0) go to 80
              maxi   = node(5,mptr)
              maxj   = node(6,mptr)
              maxip1 = maxi + 1
              maxjp1 = maxj + 1
              loc    = igetsp(maxip1 * maxjp1 * nvar)
              node(7, mptr)  = loc
              time   = rnode(12, mptr)
              intime = intcnt(lbase)
              mitot = maxi-1 + 2*lwidth
              mjtot = maxj-1 + 2*lwidth
              !  3 * to include ncount and numHoods
              node(14, mptr)  = igetsp(3*mitot*mjtot)
c
              call setirr(alloc(node(14,mptr)),mitot,mjtot,mptr,quad)
              lstgrd = node(17,mptr)
c
c      We now fill in the values for grid mptr using filvar. It uses
c      piecewise linear interpolation to obtain values from the
c      (lcheck - 1) grid, then overwrites those with whatever (lcheck)
c      grids are available. We take advantage of the fact that the
c      (lcheck - 1) grids have already been set, and that the distance
c      between any point in mptr and a (lcheck - 1) and (lcheck - 2)
c      interface is at least one (lcheck - 1) cell wide.
c
 
           mxcip1 = (maxi - 1)/intrat(lcheck-1) + 2
           mxcjp1 = (maxj - 1)/intrat(lcheck-1) + 2
           ivalc  = igetsp(mxcip1*mxcjp1*(nvar+1))
           islope = igetsp(2*mxcip1)
           xl = rnode(1,mptr)
           xr = rnode(7,mptr)
           yb = rnode(2,mptr)
           yt = rnode(4,mptr)
           hx = hxposs(lcheck)
           hy = hyposs(lcheck)
           hxcrse = hxposs(lcheck-1)
           hycrse = hyposs(lcheck-1)
           irc    = ivalc + mxcip1*mxcjp1*nvar
           lstnew = -2
           call tirrfil(alloc(irc),mxcip1,mxcjp1,xl-hxcrse,xr+hxcrse,
     1               yb-hycrse,yt+hycrse,lcheck-1,lstnew)
 
           call filval(alloc(loc),maxip1,maxjp1,hx,hy,lcheck,time,
     1                 intime,alloc(ivalc),mxcip1,mxcjp1,
     2                 alloc(node(14,mptr)),mitot,mjtot,
     3                 xl,xr,yb,yt,nvar,intrat(lcheck-1),
     4                 alloc(islope),alloc(islope+mxcip1),lwidth,lstgrd,
     5                 alloc(irc),mptr,lstnew)
 
           call reclam(ivalc,mxcip1*mxcjp1*(nvar+1))
           call reclam(islope,2*mxcip1)
 
           mptr = node(10, mptr)
           go to 10
c
c  done filling new grids at level. move them into lstart from newstl
c  (so can use as source grids for filling next level). can also
c  get rid of loc. 7 storage for old level.
c
 80   mptr = lstart(lcheck)
 85   if (mptr .eq. 0) go to 90
          maxip1 = node(5,mptr)+1
          maxjp1 = node(6,mptr)+1
          call reclam(node(7,mptr),maxip1*maxjp1*nvar)
          call lstput(node(17,mptr))
          mitot = maxip1-2 + 2*lwidth
          mjtot = maxjp1-2 + 2*lwidth
          !  3 * to include ncount and numHoods
          call reclam(node(14,mptr),3*mitot*mjtot)
          mold   = mptr
          mptr   = node(10,mptr)
          call putnod(mold)
          go to 85
 90   lstart(lcheck) = newstl(lcheck)
      lcheck = lcheck + 1
      go to 5
c
 99   lfine = lfnew
c
c  grid structure now complete again. safe to print, etc. assuming
c  things initialized to zero in nodget.
c
      return
      end
