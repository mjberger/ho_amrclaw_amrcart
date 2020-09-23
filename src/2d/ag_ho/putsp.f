c
c ----------------------------------------------------------
c
      subroutine putsp(lbase,level,nvar)
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
c
c reclaim list space in nodes 15 and 16 for grids at level
c
c first compute max. space allocated in node 15.
c
      if (level .eq. lfine) go to 30
c
      mptr  = lstart(level)
 20      call reclam(node(15,mptr), 5*listsp(level))
         node(15,mptr) = 0
      mptr  = node(10,mptr)
      if (mptr .ne. 0) go to 20
c
 30    if (level .eq. lbase) go to 99
      mptr = lstart(level)
 40       maxi  = node(5, mptr)
          maxj  = node(6, mptr)
          ikeep = (maxi-1)/intrat(level-1)
          jkeep = (maxj-1)/intrat(level-1)
          lenbc = nvar*2*(ikeep+jkeep)
          call reclam(node(16,mptr), lenbc)
          mptr = node(10,mptr)
          if (mptr .ne. 0) go to 40
c
 99    return
       end
