c
c --------------------------------------------------------------
c
      subroutine stst1
c
      implicit double precision (a-h,o-z)
      common   /space/  lfree(150,2), lenf, idimf
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
      include "cirr.i"
c
c     intialize a few variables needed before calling user set up
c     routine domain.
c     the spatial and temporal stepsizes are set. the node array
c     is kept as a linked list of free nodes.  "ndfree" points to the
c    head of the list, i.e.-first free node.  use first row of each
c    col to hold this pointer, set by the macro "nextfree".
c    the free space list, managed in lfree, will have first and
c    last positions filled with an allocation of zero words,
c    to avoid boundary cases.
c
      ndfree = 1
      do 10 i   = 1, maxgr
      node(2,i) = i+1
 10   continue
c
c the last free node will have a null pointer
 
      node(2, maxgr) = 0
c
      lfine = 1
      do 20 i  = 1,allocsize
      alloc(i) = 0.0d0
 20   continue
c
c  initialize linked list of alloc storage as well.
c  first and last locations are dummy placeholders of zero words
c  of allocation each, to avoid boundary cases.
      idimf = 150
      do  40 i  = 1, idimf
        lfree(i,1) = 0
        lfree(i,2) = 0
 40   continue
      lfree(3,1) = allocsize + 2
      lfree(2,1) = 1
      lfree(2,2) = allocsize
      lenf       = 3
c
c the linked list of irregular areas and indices will be
c threaded through the ix/iy arrays.
c
      do 50 i = 1, irrsize-1
 50   nxtirr(i) = i+1
      nxtirr(irrsize) = 0
      lhead = 1
c
c after kcheck integrations of parent grid, move its refinements.
c finest level grid never needs to have its finer subgrids moved.
c
      do 60 i   = 1, maxlv
      lstart(i) = 0
 60   icheck(i) = 0
c
c finish initializing spatial and counting arrays
c
      level      = 2
      rr         = dfloat(intrat(1))
 70   if (level .gt. mxnest) go to 80
          hxposs(level) = hxposs(level-1) / rr
          hyposs(level) = hyposs(level-1) / rr
          rr            =  intrat(level)
          level         = level + 1
      go to 70
 80   continue

      return
      end
