c
c -----------------------------------------------------------
c
      subroutine unpack(mptr, bot, top, left, right)
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
      logical   phys, bot, top, left, right
c
c  set the 4 logical variables describing the boundaries for grid mptr
c
       x = (rnode(1,mptr)+rnode(7,mptr))/2.
       y = rnode(2,mptr)
       if (phys(x,y)) then
          bot = .true.
        else
          bot = .false.
       endif
       x = rnode(1,mptr)
       y = (rnode(2,mptr)+rnode(4,mptr))/2.
       if (phys(x,y)) then
          left = .true.
        else
          left = .false.
       endif
       x = (rnode(1,mptr) + rnode(7,mptr))/2.
       y = rnode(4,mptr)
       if (phys(x,y)) then
          top = .true.
        else
          top = .false.
       endif
       x = rnode(7,mptr)
       y = (rnode(2,mptr)+rnode(4,mptr))/2.
       if (phys(x,y)) then
          right = .true.
        else
          right = .false.
       endif
       return
       end
