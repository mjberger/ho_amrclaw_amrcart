c
c ----------------------------------------------------------
c
      subroutine prepf(level,nvar)
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
c prepare new fine grids to save fluxes after each integration step
c for future conservative fixing of coarse grids
c save all boundary fluxes of fine grid (even if phys. bndry.) -
c but only save space for every intrat. (remember - 4 fluxes).
c
      mkid = lstart(level)
 10   if (mkid .eq. 0) go to 99
          maxi  = node(5, mkid)
          maxj  = node(6, mkid)
          ikeep = (maxi-1)/intrat(level-1)
          jkeep = (maxj-1)/intrat(level-1)
          lenbc = nvar*2*(ikeep+jkeep)
          node(16,mkid) = igetsp(lenbc)
	  ist   = node(16,mkid)
	  do 20 i = 1, lenbc
 20          alloc(ist+i-1) = 0.d0
          mkid = node(10,mkid)
          go to 10
 99    return
       end
