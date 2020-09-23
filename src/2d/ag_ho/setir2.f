c
c ------------------------------------------------------
c
       subroutine setir2(irreg2,mi2tot,mj2tot,mptr,lstgrd)
c
       implicit double precision (a-h,o-z)
       dimension  irreg2(mi2tot,mj2tot)
       integer srcptr, srclev, gptr
       logical outside

      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
      include "cirr.i"
      dimension ialloc(2*allocsize)
      equivalence(alloc,ialloc)

       outside(x,y) = (x .lt. tolb) .or. (y .lt. tolb) .or.
     *                (x .gt. xprob-tolb) .or. (y .gt. yprob-tolb)
      iirrc(i,j)      = 2*(locirc-1)+1 + i - 1 + misrc*(j-1)
c
c make the irregular array for the (interior) coarsened grid using the 
c irregular array for next coarser level grid.
c
       tolb  = 1.d-8
       hx   = rnode(9,mptr)
       hy   = rnode(10,mptr)
       hx2  = 2.d0*hx
       hy2  = 2.d0*hy
       xlow = rnode(1,mptr) - lwidth*hx2
       ylow = rnode(2,mptr) - lwidth*hy2
       gptr = 0
       level = node(4,mptr) - 1
       if (level .lt. 1) go to 99

       do 10 i = 1, mi2tot
       do 10 j = 1, mj2tot
          x = xlow + (i-.5)*hx2
          y = ylow + (j-.5)*hy2
	  if (outside(x,y)) then
	    irreg2(i,j) = lstgrd
	  else
             call getsrc(x,y,level,gptr, srcptr, srclev)
             if (srclev .lt. level) then
               irreg2(i,j) = -1
             else
               locirc = node(14, srcptr)
               misrc  = node(5,srcptr)-1+2*lwidth
               xsrc   = rnode(1,srcptr)-lwidth*hx2
               ysrc   = rnode(2,srcptr)-lwidth*hy2
	       isrc   = (x-xsrc)/hx2 + 1
	       jsrc   = (y-ysrc)/hy2 + 1
               ksrc = ialloc(iirrc(isrc,jsrc))
	       lstsrc = node(17,srcptr)
	       if (ksrc .ne. lstsrc) then
	         irreg2(i,j) = ksrc
	       else
	         irreg2(i,j) = lstgrd
	       endif
               gptr = srcptr
             endif
          endif
 10    continue
c
 99    return
       end
