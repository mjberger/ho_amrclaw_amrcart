c
c --------------------------------------------------------------
c
      subroutine join (lbase, nvar)
c
      implicit double precision (a-h,o-z)
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
c
c  join = allocate 2nd storage block for each grid.
c         subr. fill now does some data structure manipulation
c         (it joins newstl into lstart a level at a time. it also
c         reclaims node 7 storage from old grids as soon as it can,
c         and resets lfnew to lfine)
c
c allocate rest of storage for new grids.
c
      level = lbase + 1
 70   if (level .gt. lfine) go to 99
          mptr = lstart(level)
 80       if (mptr .eq. 0) go to 90
              maxi   = node(5, mptr)
              maxj   = node(6, mptr)
              maxip1 = maxi + 1
              maxjp1 = maxj + 1
              maxipb = maxi + lwidth - 1
              maxjpb = maxj + lwidth - 1
              nwords       = maxip1*maxjp1*nvar
	      if (level .lt. mxnest) node(8,mptr) = igetsp(nwords)
c             node(2,mptr) = igetsp(2*nvar*lwidth*(maxipb + maxjpb))
          mptr         =  node(10, mptr)
          go to 80
 90       continue
          level = level + 1
      go to 70
c
 99   return
      end
