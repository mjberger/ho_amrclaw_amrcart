c
c -----------------------------------------------------------
c
      subroutine colate (badpts, len, lcheck, nvar, ignore, lbase)
c
      implicit double precision (a-h,o-z)
      dimension          badpts(2,len)
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
      integer  gptr, srcptr, srclev
      logical  top, bot, left, right
      data     goodpt/0.0/
c
      iadd(i,j,ivar) = loc + i - 1 + maxip1*((ivar-1)*maxjp1 + j - 1)
c
c *************************************************************
c
c colate = takes all error planes with flagged pts at level lcheck
c          and puts the (x,y) coords into badpts. array.
c          to insure proper nesting, do not count flagged points
c          on grid boundary, unless another grid neighbors it, or
c          underlying grid will also be regridded (enlarged) to
c          accomodate it.
c
c *************************************************************
c
      mptr = lstart(lcheck)
      time = rnode(12, mptr)
      index  = 1
      ignore = 0
 10   if (mptr .eq. 0) go to 99
          maxi   = node(5,mptr)
          maxj   = node(6,mptr)
          maxip1 = maxi + 1
          maxjp1 = maxj + 1
          hx     = rnode(9, mptr)
          hy     = rnode(10,mptr)
          xst    = rnode(1,mptr)
          yst    = rnode(2,mptr)
          loc    = node(8,mptr)
          gptr   = mstart
          call unpack(mptr,bot,top,left,right)
c
          do 20 i   = 2, maxi
          do 20 j   = 2, maxj
          if (alloc(iadd(i,j,1)) .eq. goodpt) go to 20
c
c  check if flagged point on the boundary would create properly
c  properly nested grid
c
          if (i .eq. 2) then
             if (left) then
               go to 16
             else
               x = rnode(1,mptr) - hx/2.d0
               y = rnode(2,mptr) + (j-1.5d0)*hy
               call getsrc(x,y,lcheck,gptr,srcptr,srclev)
               if (srclev .ge. lcheck) then
                  go to 16
               else
                  ignore = ignore + 1
                  go to 20
               endif
             endif
       endif
          if (i .eq. maxi) then
             if (right) then
               go to 16
             else
               x = rnode(7,mptr) + hx/2.d0
               y = rnode(2,mptr) + (j-1.5d0)*hy
               call getsrc(x,y,lcheck,gptr,srcptr,srclev)
               if (srclev .ge. lcheck) then
                  go to 16
               else
                  ignore = ignore + 1
                  go to 20
               endif
             endif
          endif
c
 16       if (j .eq. 2) then
             if (bot) then
               go to 19
             else
               x = rnode(1,mptr) + (i-1.5d0)*hx
               y = rnode(2,mptr) - hy/2.d0
               call getsrc(x,y,lcheck,gptr,srcptr,srclev)
               if (srclev .ge. lcheck) then
                  go to 19
               else
                  ignore = ignore + 1
                  go to 20
               endif
            endif
       endif
          if (j .eq. maxj) then
             if (top) then
               go to 19
             else
               x = rnode(1,mptr) + (i-1.5d0)*hx
               y = rnode(4,mptr) + hy/2.d0
               call getsrc(x,y,lcheck,gptr,srcptr,srclev)
               if (srclev .ge. lcheck) then
                  go to 19
               else
                  ignore = ignore + 1
                  go to 20
               endif
             endif
          endif
c
c flagged point OK - count it
c
 19           badpts(1,index) = xst + hx*(dfloat(i)-1.5d0)
              badpts(2,index) = yst + hy*(dfloat(j)-1.5d0)
              index = index + 1
 20       continue
c
          mptr = node(10, mptr)
      go to 10
c
 99   if (ignore .ne. 0) write(6,900) ignore, time, lcheck
 900  format( i5,' flagged points ignored at time ',e15.7,
     1        ' on level ',i4,' grids ')
      return
      end
