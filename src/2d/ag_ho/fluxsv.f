c
c ----------------------------------------------------------
c
      subroutine fluxsv(mptr,xflux,yflux,listbc,ndimx,ndimy,
     1		        nvar,maxsp)
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
      include "cirr.i"
      dimension xflux(ndimx,ndimy,nvar), yflux(ndimx,ndimy,nvar)
      dimension listbc(5,maxsp)
c
c  listbc holds info for the coarse grid mptr to add its fluxes (in
c  xflux and yflux) into a fine grid.
c  xflux holds 'f' fluxes, yflux holds 'g' fluxes.
c  they are dimensioned at ndimx for other reasons - the fluxes
c  are actually stored starting at 1,1 = lower left cell.
c  cell variables and fluxes look like this (tricky indexing):
c                        g i-1,j
c            f i-1,j-1   u i,j     f i,j-1
c                        g i-1,j-1
c
c  there is at least 1 item in listbc or else would not be called from advanc.
c
c    The above statement is not true. We return if there are no entries
c     in listbc .
 
      if (listbc(1,1).eq.0) return
 
      ispot   = 1
      level   = node(4,mptr)
      refrat  = dfloat(intrat(level))
c     refrat  = 1  ! for steady calcs, no ref in time
c
 10      mkid   = listbc(4,ispot)
         intopl = listbc(5,ispot)
         kidlst = node(16,mkid)
         i      = listbc(1,ispot)
         j      = listbc(2,ispot)
         inlist = kidlst + nvar*(intopl-1) - 1
         if (listbc(3,ispot) .eq. 1) then
c	    ::::: Right side
	    do 100 ivar = 1, nvar
               alloc(inlist + ivar) = -xflux(i-1,j-1,ivar)*refrat
100	    continue
         endif
         if (listbc(3,ispot) .eq. 2) then
c	    ::::: Bottom side
	    do 200 ivar = 1, nvar
               alloc(inlist + ivar) = -yflux(i-1,j,ivar)*refrat
200	    continue
         endif
         if (listbc(3,ispot) .eq. 3) then
c	    ::::: Left side
	    do 300 ivar = 1, nvar
               alloc(inlist + ivar) = -xflux(i,j-1,ivar)*refrat
300	    continue
         endif
         if (listbc(3,ispot) .eq. 4) then
c	    ::::: Top side
	    do 400 ivar = 1, nvar
               alloc(inlist + ivar) = -yflux(i-1,j-1,ivar)*refrat
400	    continue
         endif
c
      ispot = ispot + 1
      if (ispot .gt. maxsp) go to 99
      if (listbc(1,ispot) .ne. 0) go to 10
c
 99   return
      end
