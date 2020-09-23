c
c ---------------------------------------------------------
c
      subroutine projec (level, numpro, nvar)
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
      logical   pprint, ovrlap
      data      pprint/.false./, goodpt/0.0/, badpt/3.0/
c
      iadd(i,j,ivar) = loc + i - 1 + maxip1*((ivar-1)*maxjp1 + j - 1)
c
c  projec *****************************************************
c  for all newly created fine grids, project area onto a coarser
c  grid 2 levels down. Used to recreate grids 1 level down, and
c  insure proper level nesting.
c
c  on entry, all coarse grids have already had error estimated, so
c  add bad flags.   count number of 'added' flags only.
c
c input parameters:
c    level = project all fine subgrids onto grids at this level.
c output parameters:
c  numpro = number of additional flagged pts. at 'level'.
c           (initialized to 0 in flglvl)
c local variables:
c     mptr    = holds coarser grid projected into
c     mkid    = grid doing the projecting
c  *************************************************************
c
      hx      = hxposs(level)
      hy      = hyposs(level)
      hxmid   = hxposs(level+1)
      hymid   = hyposs(level+1)
      levpro  =  level + 2
      buffx = bzonex
      buffy = bzoney
c
      mkid = newstl(levpro)
 10   if (mkid .eq. 0) go to 90
       xst  = rnode(1,mkid) - buffx * hxmid
       yst  = rnode(2,mkid) - buffy * hymid
       xend = rnode(7,mkid) + buffx * hxmid
       yend = rnode(4,mkid) + buffy * hymid
c
c  project onto any overlapping coarse grid. faster than insuring
c  that one pt. not projected twice onto 2 difference coarse
c  (intersecting) grids.
c
      mptr    = lstart(level)
 20   if (.not. ovrlap(mkid,mptr)) go to 70
      maxi    = node(5, mptr)
      maxj    = node(6, mptr)
      maxip1  = maxi + 1
      maxjp1  = maxj + 1
      loc     = node(8, mptr)
      xpar    = rnode(1,mptr) - hx/2.d0
      ypar    = rnode(2,mptr) - hy/2.d0
c
      ist  = max(idint((xst -xpar)/hx +1.1), 2)
      jst  = max(idint((yst -ypar)/hy +1.1), 2)
      iend = min(idint((xend-xpar)/hx +1.1), maxi)
      jend = min(idint((yend-ypar)/hy +1.1), maxj)
c
       do 60 j = jst, jend
       do 60 i = ist, iend
           if (alloc(iadd(i,j,1)) .ne. goodpt) go to 60
               alloc(iadd(i,j,1)) = badpt
               numpro                 = numpro + 1
           if (pprint) write(6,101) i,j,mptr,mkid
101        format(' pt.',2i5,' of grid ',i5,' projected from grid',i5)
 60    continue
c
c  done with gridpt. loop for grid mkid.
c
 70     mptr = node(10,mptr)
        if (mptr .ne. 0) go to 20
c
 80     mkid = node(10, mkid)
        go to 10
c
 90   if (numpro .eq. 0) go to 99
      write(6,102) numpro,level
 102  format(i7,' more pts. projected to level ',i5)
      mptr = lstart(level)
c
 95   write(6,103) mptr
 103  format(/,'  flagged pts. for grid ',i4,':')

      loc     = node(8, mptr)
      maxip1  = node(5, mptr) + 1
      maxjp1  = node(6, mptr) + 1
c
      do 110 jj = 2, maxj
      j        = maxj + 2 - jj
      write(6,104) (nint(alloc(iadd(i,j,1))),i=2,maxi)
104   format(1h ,80i1)
 110  continue
      mptr = node(10,mptr)
      if (mptr .ne. 0) go to 95
c
 99   return
      end
