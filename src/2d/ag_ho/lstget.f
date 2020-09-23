c
c -------------------------------------------------------------
c
      integer function lstget(dummy)
c
      implicit double precision (a-h,o-z)
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "cirr.i"
c
c first free spot in linked list of irregular info. pointed to
c by lhead.
c
      if (lhead .ne. 0) go to 10
	 write(6,100)
 100     format(' out of irregular list space ')
	 stop
c
 10   lstget = lhead
      lhead = iabs(nxtirr(lhead))
c
      return
      end
