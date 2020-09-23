c
c -------------------------------------------------------------
c
      subroutine lstput(lstgrd)
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
c lstgrd = start of linked list for the grid. is never 0, since
c the "regular cell" info. (at least area) uses at least
c one spot on the info. list.
c
c find the end of the grid's list
c
       lptr = lstgrd
 10    if (nxtirr(lptr) .eq. 0) go to 20
          lptr = iabs(nxtirr(lptr))
          go to 10
c
c now lptr points to end of the grid's list. attach to rest of list
c
 20    nxtirr(lptr) = lhead
       lhead = lstgrd
c
       return
       end
