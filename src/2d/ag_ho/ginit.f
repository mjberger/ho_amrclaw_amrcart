c
c -------------------------------------------------------------
c
      subroutine ginit(msave, first, nvar, quad, gradThreshold)
c
      implicit double precision (a-h,o-z)
      logical first, quad
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
c
c ginit = init's all grids at 'level'  by calling to finit
c         if first = true, (first call to init), then allocate a
c         storage area, else use one already allocated.
c
      mptr  = msave
 10   if (mptr .eq. 0) go to 99
          maxi    = node(5,mptr)
          maxj    = node(6, mptr)
          maxip1  = maxi + 1
          maxjp1  = maxj + 1
          mitot   = maxi-1 + 2*lwidth
          mjtot   = maxj-1 + 2*lwidth
          corn1   = rnode(1,mptr)
          corn2   = rnode(2,mptr)
          hx      = rnode(9,mptr)
          hy      = rnode(10,mptr)
          if(.not. (first)) go to 20
              loc             = igetsp(maxip1*maxjp1*nvar)
              node(7,mptr)    = loc
              rnode(12, mptr) = 0.0
c             get space for irr array, and for shifted values for quad.recon.
              !! get space for ncount and numHoods too
              node(14, mptr) = igetsp(3*mitot*mjtot)
              call setirr(alloc(node(14,mptr)),mitot,mjtot,mptr,quad,
     .                    gradThreshold)
              go to 30
 20       continue
c
c  if 2nd time through, put initial values in 8 so finer grids
c  can be advanced with interpolation of their boundary values.
c  new time soln should still be in location 7.
c
          loc     = node(8,mptr)
c
   30     continue
c         # added irr,mitot,mjtot,lstgrd to finit call
c         # is this lstgrd value right???
          lstgrd  = node(17,mptr)
          call finit(alloc(loc),nvar,maxip1,maxjp1,
     1               corn1,corn2,hx,hy,alloc(node(14,mptr)),
     2               mitot,mjtot,lstgrd)
c
          mptr  = node(10, mptr)
      go to 10
c
c
 99   return
      end
