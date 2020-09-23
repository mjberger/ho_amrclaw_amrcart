c
c -------------------------------------------------------------
c
       logical function inbase(x,y,lbase)
c
c  inbase = true if pt. (x,y) in any grid at level lbase, else false.
c   do in-line cntain checking here, for speed and less subr. overhead.
c   no tolerances (hxmarg) since called with cell centered points.
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
      inbase = .true.
      mpar   = lstart(lbase)
c
 10   if ( x .lt. rnode( 7,mpar) .and.
     1     x .gt. rnode( 1,mpar) .and.
     2     y .gt. rnode( 2,mpar) .and.
     3     y .lt. rnode( 4,mpar) )
     4   go to 99
c
      mpar  = node(10,mpar)
      if (mpar .ne. 0) go to 10
      inbase = .false.
c
 99   return
      end
