c
c ----------------------------------------------------------------
c
       subroutine trimgr(mnew,hx,hy)
       implicit double precision (a-h,o-z)
       logical            graf
       parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
c
c nclust grids have been formed by bisecting and merging.
c remove some inefficiency here by reducing the amount of overlap
c in the simplest case. check that grid doesn't completely disappear.
c
       hx2 = hx / 2.d0
       hy2 = hy / 2.d0
       hxmarg = hx / 100.
       hymarg = hy / 100.
c
       iprev = 0
       itry  = mnew
 5     icl   = mnew
 10    if (icl .eq. itry) go to 20
c
c do grids intersect in this special way
c
       xst  = rnode(1,icl)
       yst  = rnode(2,icl)
       xend = rnode(7,icl)
       yend = rnode(4,icl)
c
c check right side first
c
       if ((rnode(7,itry) .ge. xst+hx2) .and.
     .    (rnode(7,itry) .le. xend+hxmarg) .and.
     .    (rnode(2,itry) .ge. yst-hymarg) .and.
     .    (rnode(4,itry) .le. yend+hymarg)) then
		  rnode(7,itry) = rnode(1,icl)
		  rnode(5,itry) = rnode(1,icl)
                  if (rnode(7,itry)-rnode(1,itry) .lt. hx2) go to 30
       endif
c
check left side next
c
       if ((rnode(1,itry) .ge. xst-hxmarg) .and.
     .    (rnode(1,itry) .le. xend-hx2) .and.
     .    (rnode(2,itry) .ge. yst-hymarg) .and.
     .    (rnode(4,itry) .le. yend+hymarg)) then
		  rnode(1,itry) = rnode(7,icl)
		  rnode(3,itry) = rnode(7,icl)
                  if (rnode(7,itry)-rnode(1,itry) .lt. hx2) go to 30
       endif
c
c check top
c
       if ((rnode(4,itry) .ge. yst+hy2) .and.
     .    (rnode(4,itry) .le. yend+hymarg) .and.
     .    (rnode(1,itry) .ge. xst-hxmarg) .and.
     .    (rnode(7,itry) .le. xend+hxmarg)) then
		  rnode(4,itry) = rnode(2,icl)
		  rnode(6,itry) = rnode(2,icl)
                  if (rnode(4,itry)-rnode(2,itry) .lt. hy2) go to 30
       endif
c
c check bottom
c
       if ((rnode(2,itry) .ge. yst-hymarg) .and.
     .    (rnode(2,itry) .le. yend-hy2) .and.
     .    (rnode(1,itry) .ge. xst-hxmarg) .and.
     .    (rnode(7,itry) .le. xend+hxmarg)) then
		  rnode(2,itry) = rnode(4,icl)
		  rnode(8,itry) = rnode(4,icl)
                  if (rnode(4,itry)-rnode(2,itry) .lt. hy2) go to 30
       endif
c
 20   icl = node(10,icl)
      if (icl .ne. 0) go to 10
c
      iprev = itry
      itry  = node(10, itry)
      if (itry .ne. 0) go to 5
      go to 99
c
c grid itry has been completely eliminated.
c
 30   continue
      inext = node(10,itry)
      if (iprev .eq. 0) then
 	 mnew = node(10,itry)
      else
         node(10,iprev) = node(10,itry)
      endif
      call putnod(itry)
      itry = inext
      if (itry .ne. 0) go to 5
c
 99   return
      end
