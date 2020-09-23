c
c ------------------------------------------------------------
c
      logical function cntain(x,y, mptr)
c
      implicit double precision (a-h,o-z)
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
c
c  cntain =  true if pt. (x,y) is cntained in the
c             rectangular grid mptr, else false.
c
      hxmarg = rnode(9,mptr) / 10.d0
      hymarg = rnode(10,mptr) / 10.d0
c
      cntain = .false.
      if ( x .lt. rnode( 7,mptr) + hxmarg .and.
     1     x .gt. rnode( 1,mptr) - hxmarg .and.
     2     y .gt. rnode( 2,mptr) - hymarg .and.
     3     y .lt. rnode( 4,mptr) + hymarg )
     4   cntain = .true.
c
      return
      end
