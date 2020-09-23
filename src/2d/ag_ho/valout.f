c
c -----------------------------------------------------
c
      subroutine valout (mgrid, lst, lend, time, nvar)
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
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      logical         flag
c
c valout = graphics output of soln values for contour or surface plots.
c    if mgrid <> 0,  output only that grid, (should have
c                          lst = lend here)
c
      call basic (time, lst, lend, 3 )
c
      write(3,100) mgrid,  nvar
100   format(10h*VALS     ,2i10)
c
      flag = .false.
      if (lst .eq. lend .and. mgrid .ne. 0) flag = .true.
c
      level = lst
10    if (level .gt. lend) go to 99
          mptr = lstart(level)
          if (flag) mptr = mgrid
20        if (mptr .eq. 0) go to 50
              maxi   = node(5,mptr)
              maxj   = node(6,mptr)
	      xlow = rnode(1,mptr)
	      ylow = rnode(2,mptr)
              maxip1 = maxi + 1
              maxjp1 = maxj + 1
              loc = node(7,mptr)
	      locirr = node(14,mptr)
	      mitot = maxi-1+2*lwidth
	      mjtot = maxj-1+2*lwidth
	      locsm = igetsp(mitot*mjtot*nvar)
              call outvar(alloc(loc),maxip1,maxjp1,nvar,mptr,
     1                    alloc(locirr),mitot,mjtot,lwidth,
     2                    node(17,mptr),rnode(9,mptr),
     3                    rnode(10,mptr),xlow,ylow,alloc(locsm))
	      call reclam(locsm,mitot*mjtot*nvar)
              mptr = node(10,mptr)
              if (flag) go to 50
          go to 20
50        continue
          level = level + 1
      go to 10
c
 99   return
      end
