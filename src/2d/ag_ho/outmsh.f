c
c ---------------------------------------------------------
c
      subroutine outmsh(mptr,outgrd,nvar)
      implicit double precision (a-h,o-z)
      logical  outgrd
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
c
c ******************************************************************
c outmsh - output the grid descriptor and optionally the values on
c          the grid (for a single grid - see "outtre" for outputing
c          a subtree)
c input parameters:
c    mptr   - ptr to grid descriptor
c    outgrd - if true, output value on grid
c special case
c     if grid has level < 1, nothing is printed. (forest pointer
c has level = 0).  this simplifies logic of outtre; also, any grid
c with level <= 0 is not properly set (the forest pointer
c is used only to provide pointers into the tree of coarsest grids)
c ******************************************************************
100   format(1x,47h+----------------------------------------------,
     *30h-----------------------------+)
c
      if (node(4,mptr) .lt. 1) go to 99
c
      write(6,100)
      write(6,101) mptr
101   format(1x,10h! grid no:,i3,63x,1h!)
      write(6,102)node(4,mptr),node(5,mptr),node(6,mptr),rnode(12,mptr)
102   format(1x,1h!,1x,11h nestlevel=,i3,9h, maxrow=,i6,9h, maxcol=,i6
     *,12h, time mult=,f11.5,7x,1h!)
      mitot    = node(5,mptr) - 1 + 2*lwidth
      mjtot    = node(6,mptr) - 1 + 2*lwidth
      write(6,202) mitot,mjtot
202   format(1x,1h!,1x,14h              ,9h  mitot =,i6,9h, mjtot =,i6
     *,' incl. ghost cells ',11x,1h!)
      write(6,103) node(7,mptr),node(8,mptr),node(15,mptr),
     1             node(16,mptr)
 103  format(1x,'! storage locs =',2i9,'    bndry locs =',2i7,12x,1h!)
      write(6,110) node(10,mptr), node(14, mptr)
110   format(1x,13h! level ptr =,i4,11h, irreg ptr,i10,38x,1h!)
      write(6,104)
      write(6,111) rnode(3,mptr),rnode(4,mptr),rnode(5,mptr),
     1             rnode(6,mptr)
      write(6,111) rnode(1,mptr),rnode(2,mptr),rnode(7,mptr),
     1             rnode(8,mptr)
c104   format(1x,23h! corners of rectangle:,53x,1h!)
104   format(1x,51h! corners of rectangle: (not including ghost cells),
     .25x,1h!)
111   format(1x,2h! ,18x,2(1h(,f11.5,2h, ,f11.5,1h)),4x,1h!)
      write(6,105) rnode(9,mptr),rnode(10,mptr),rnode(11,mptr)
105   format(1x,7h! hrow=,f8.5,7h, hcol=,f8.5,8h, ktime=,f8.5,30x,1h!)
      write(6,100)
c
      if (.not. outgrd) go to 99
c output the grid
 20   maxip1   = node(5,mptr) + 1
      maxjp1   = node(6,mptr) + 1
      loc      = node(7,mptr)
      locirr   = node(14,mptr)
      mitot    = maxip1 - 2 + 2*lwidth
      mjtot    = maxjp1 - 2 + 2*lwidth
      lstgrd   = node(17,mptr)
         call outval(alloc(loc),maxip1,maxjp1,nvar,
     1               alloc(locirr),mitot,mjtot,mptr,outgrd,lstgrd)
 99   return
      end
