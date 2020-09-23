 
c
c -----------------------------------------------------------
c
      logical function ovrlap(m1ptr,m2ptr)
      implicit double precision (a-h,o-z)
      logical            xint, yint
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,
     *  mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
c
c ******************************************************************
c ovrlap - determine if two grids intersect
c input paramters:
c    m1ptr, m2ptr - ptrs to grid descriptors
c returns:
c    true if grids described by m1ptr and m2ptr overlap
c note:
c    this routine will return true if one of the grids is wholly
c    contained in the other. also, grids that share an edge (top of
c    one is the bottom of the other) are considered overlapping.
c ******************************************************************
c
      xm     = rnode(9, m1ptr) / 100.
      ym     = rnode(10,m1ptr) / 100.
      xst    = rnode(1,m1ptr)
      xend   = rnode(7,m1ptr)
      yst    = rnode(2,m1ptr)
      yend   = rnode(4,m1ptr)
      xleft  = rnode(1,m2ptr)
      xright = rnode(7,m2ptr)
      ybot   = rnode(2,m2ptr)
      ytop   = rnode(4,m2ptr)
      xint   = .false.
      yint   = .false.
      ovrlap = .false.
c
      if  ( ((xleft .le. xst + xm) .and. (xright .ge. xst - xm)) .or.
     1      ((xleft .ge. xst - xm) .and. (xleft  .le. xend + xm)) )
     2    xint = .true.
      if  ( ((ybot .le. yst + ym) .and. (ytop .ge. yst - ym)) .or.
     1      ((ybot .ge. yst - ym) .and. (ybot  .le. yend + ym)) )
     2    yint = .true.
c
      if (xint .and. yint) ovrlap = .true.
c
      return
      end
