c
c --------------------------------------------------------------
c
      subroutine outtre(mlev,outgrd,nvar)
c
      implicit double precision (a-h,o-z)
      logical  outgrd
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
c
c ******************************************************************
c outtre - output subtree
c input parameters:
c    mlev   - ptr to subtree to output i.e., start at level(mlev)
c    outgrd - if true, output the values on the grid
c tree is output from 'level' to finest level.
c ******************************************************************
c
      write (6,1)
1     format(1x,14hthe subtree is)
c
      level = node(4, mlev)
10    if (level .gt. lfine) go to 99
          mptr    = lstart(level)
 20       if (mptr .eq. 0) go to 30
              call outmsh(mptr,outgrd,nvar)
              mptr = node(10, mptr)
          go to 20
 30       continue
          level = level + 1
      go to 10
c
 99   return
      end
