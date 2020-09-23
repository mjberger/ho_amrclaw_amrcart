c
c -----------------------------------------------------------
c
      subroutine sethk(mptr1)
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
c ******************************************************************
c sethk  - set the h's, k, and maxnumrow, col for the grids
c          starting at msave.
c set the prevlevel (backwards pointing) grid pointers
c input parameter:
c    mptr1   - starting ptr to grid to find hrow, hcol, k, maxi, maxj for
c ******************************************************************
c
      mptr   = mptr1
      mlevel = node(4,mptr)
      iprev  = 0
c
 10   if (mptr .eq. 0) go to 99
        rnode(11,mptr) = possk(mlevel)
        rnode( 9,mptr) = hxposs(mlevel)
        rnode(10,mptr) = hyposs(mlevel)
c
        node(5,mptr) = (rnode(7,mptr)-rnode(1,mptr)+.01*hxposs(mlevel))
     1  / hxposs(mlevel)  + 1
        node(6,mptr) = (rnode(4,mptr)-rnode(2,mptr)+.01*hyposs(mlevel))
     1  / hyposs(mlevel)  + 1
c
        node(12,mptr) = iprev
        iprev         = mptr
c
        mptr = node(10,mptr)
      go to 10
c
 99   return
      end
