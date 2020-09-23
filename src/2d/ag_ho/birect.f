c
c --------------------------------------------------
c
      subroutine birect(mptr1)
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
      data     maxrow,maxcol/1200,1200/
c     data     maxrow,maxcol/140,140/
c     data     maxrow,maxcol/58,58/
c
c check each grid, starting with mptr1 (either newstl or lstart)
c to see that it has no more than maxrow and maxcol points.
c needed so that scratch array space in method not exceeded.
c1x
         return !! only one large grid



         mptr  = mptr1
         level = node(4,mptr)
         hx    = hxposs(level)
         hy    = hyposs(level)
c
10    continue
      c1x = rnode(1,mptr)
      c4x = rnode(7,mptr)
      c1y = rnode(2,mptr)
      c2y = rnode(4,mptr)
      nx  = 1 + int((c4x - c1x + .0001*hx)/hx)
      ny  = 1 + int((c2y - c1y + .0001*hy)/hy)
c
c check number of rows first - if too many, bisect grid with vertical
c line down the middle. make sure new grid corners are anchored
c on coarse grid point.
c
      if (nx .gt. maxrow) then
 
        nxl    = (nx+1)/2
	if (level .gt. 1) nxl = (nxl/intrat(level-1))*intrat(level-1)+1
        nxr    = nx - nxl + 1
        cxmid  = (nxl-1)*hx + c1x
 
        mptrnx = nodget(dummy)
        node(10,mptrnx) = node(10,mptr)
        node(10,mptr)   = mptrnx
 
        rnode(5,mptr) = cxmid
        rnode(7,mptr) = cxmid
        rnode(1,mptrnx) = cxmid
        rnode(2,mptrnx) = c1y
        rnode(3,mptrnx) = cxmid
        rnode(4,mptrnx) = c2y
        rnode(5,mptrnx) = c4x
        rnode(6,mptrnx) = c2y
        rnode(7,mptrnx) = c4x
        rnode(8,mptrnx) = c1y
        rnode(12,mptrnx)= rnode(12,mptr)
        rnode(11,mptrnx)= rnode(11,mptr)
        node(4,mptrnx)  = node(4,mptr)
        go to 10
c
c check number of columns next - if too many, bisect grid with horizontal
c line down the middle
c
      else if (ny .gt. maxcol) then
 
        nyl    = (ny+1)/2
	if (level .gt. 1) nyl = (nyl/intrat(level-1))*intrat(level-1)+1
        nyr    = ny - nyl + 1
        cymid  = (nyl-1)*hy + c1y
 
        mptrnx = nodget(dummy)
        node(10,mptrnx) = node(10,mptr)
        node(10,mptr)   = mptrnx
 
        rnode(4,mptr) = cymid
        rnode(6,mptr) = cymid
        rnode(1,mptrnx) = c1x
        rnode(2,mptrnx) = cymid
        rnode(3,mptrnx) = c1x
        rnode(4,mptrnx) = c2y
        rnode(5,mptrnx) = c4x
        rnode(6,mptrnx) = c2y
        rnode(7,mptrnx) = c4x
        rnode(8,mptrnx) = cymid
        node(4,mptrnx)  = node(4,mptr)
        rnode(12,mptrnx)= rnode(12,mptr)
        rnode(11,mptrnx)= rnode(11,mptr)
        go to 10
c
c  grid ok - check the next
c
      else
        mptr = node(10,mptr)
        if (mptr.ne.0) go to 10
 
      endif
 
      return
      end
