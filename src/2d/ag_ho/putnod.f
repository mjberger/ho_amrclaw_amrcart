c
c -------------------------------------------------------------
c
      subroutine putnod (mptr)
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
c putnod = return mptr node to the linked list kept in node array.
c
      node(2, mptr) = ndfree
      ndfree        = mptr
c
      return
      end
