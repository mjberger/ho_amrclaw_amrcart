c
c --------------------------------------------------------------
c
      subroutine outlev(levnew,nvar)
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
      data     outgrd /.false./
c
c ******************************************************************
c outlev - output a level in newstl  (grid being formed).
c input parameters:
c    levnew   - ptr to level to output
c ******************************************************************
c
      write (6,1)
1     format(1x,'the new level is')
c
          mptr    = newstl(levnew)
 20       if (mptr .eq. 0) go to 99
              call outmsh(mptr,outgrd,nvar)
              mptr = node(10, mptr)
          go to 20
c
 99   return
      end
