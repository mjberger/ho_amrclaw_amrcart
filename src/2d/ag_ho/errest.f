c
c -------------------------------------------------------------
c
      subroutine errest (nvar,numbad,lcheck,rbuff,nbuff,iflow,
     .                   steady)
c
      implicit double precision (a-h,o-z)
      logical    graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "cirr.i"
      include "calloc.i"
      common   /stats/  evol, rvol, rvoll(maxlv), lentot, lenmax
      dimension rbuff(3,nbuff), intinc(maxlv)
      logical   rflag, vtime, steady
      logical lastout/.false./
      data      rflag/.false./, vtime/.false./
c     #   no sense computing new time step if just for error estimation,
c     #   so set vtime false here.
 
c for all grids at level lcheck:
c  estimate the error by taking a large (2h,2k) step based on the
c  values in the old storage loc., then take one regular (and for
c  now wasted) step based on the new info.   compare using an
c  error relation for a pth order  accurate integration formula.
c  flag error plane as either bad (needs refinement), or good.
c
c  call the regular integrator on a grid twice as coarse.
c  initialize such a grid directly, instead of trick dimensioning.
c  this is to make other l1 type estimates easier to program.
c
       numbad = 0
       intinc(mxnest) = 1
       do 3 i = 1, mxnest-1
          level = mxnest - i
 3        intinc(level) = intinc(level+1)*intrat(level)
c
       mptr = lstart(lcheck)
 5     continue
          maxi   = node(5, mptr)
          maxj   = node(6, mptr)
          maxip1 = maxi + 1
          maxjp1 = maxj + 1
	  mitot  = maxi-1 + 2*lwidth
	  mjtot  = maxj-1 + 2*lwidth
          locnew = node(7,mptr)
          locold = node(8,mptr)
          mi2tot = (maxi-1)/2  + 2*lwidth
          mj2tot = (maxj-1)/2  + 2*lwidth
          time   = rnode(12,mptr)
	  intime = intcnt(lcheck)
          dt     = rnode(11,mptr)
          tpre   = time - dt
	  intpre = intime - intinc(lcheck)
          hx     = rnode(9,mptr)
          hy     = rnode(10,mptr)
c
          xl = rnode(1,mptr)
          yb = rnode(2,mptr)
          xr = rnode(7,mptr)
          yt = rnode(4,mptr)
c 
c     prepare double the stencil size worth of boundary values,
c            then coarsen them for the giant step integration.
c 
	  midub = maxi-1+4*lwidth
	  mjdub = maxj-1+4*lwidth
          locdub = igetsp(midub*mjdub*nvar)
	  do 10 i = 1, midub*mjdub*nvar
 10       alloc(locdub+i-1) = 1.d0
	  locbgc = igetsp(mi2tot*mj2tot*nvar)
	  node(13,mptr) = locbgc
          lw2 = 2*lwidth
	  call prem(alloc(locdub),alloc(locold),nvar,maxip1,maxjp1,
     1              midub,mjdub,lw2)
          call bdcrse(xl,xr,yb,yt,hx,hy,tpre,intpre,lcheck,nvar,lw2,
     1                alloc(locdub),midub,mjdub,alloc(locbgc),mi2tot,
     2                mj2tot,alloc(node(14,mptr)),mitot,mjtot)
	  call reclam(locdub,midub*mjdub*nvar)
c
c
c We now call bound at time t = time, in preparation
c for calculating the solution on the grid mptr for error estimation.
c 
          locbig = igetsp(mitot*mjtot*nvar)
	  node(2,mptr) = locbig
          call prem(alloc(locbig),alloc(locnew),nvar,maxip1,maxjp1,
     1              mitot,mjtot,lwidth)
c         call bound(xl,xr,yb,yt,hx,hy,time,intime,lcheck,nvar,lwidth,
c    1               alloc(locbig),mitot,mjtot,alloc(node(14,mptr)),
c    2               node(17,mptr))
c
       mptr = node(10,mptr)
       if (mptr .ne. 0) go to 5
c
       mptr = lstart(lcheck)
 25    continue
	  maxi   = node(5, mptr)
	  maxj   = node(6, mptr)
	  maxip1 = maxi + 1
	  maxjp1 = maxj + 1
	  mitot  = maxi-1 + 2*lwidth
	  mjtot  = maxj-1 + 2*lwidth
	  mi2    = (maxi-1) / 2 + 1
	  mj2    = (maxj-1) / 2 + 1
	  mi2p1  = mi2 + 1
	  mj2p1  = mj2 + 1
          mi2tot = mi2-1 + 2*lwidth
          mj2tot = mj2-1 + 2*lwidth
	  dt     = rnode(11,mptr)
	  dt2    = 2. * dt
	  hx     = rnode(9, mptr)
	  hy     = rnode(10,mptr)
	  hx2    = 2.*hx
	  hy2    = 2.*hy
c
          locsc1 = igetsp(mi2tot*mj2tot*nvar)
          locsc5 = igetsp(mi2tot*mj2tot*nvar)
          locqx  = igetsp(mi2tot*mj2tot*nvar)
          locqy  = igetsp(mi2tot*mj2tot*nvar)
	  locbgc = node(13,mptr)
	  locirc = igetsp(mi2tot*mj2tot)
	  lst2   = node(17,mptr)
	  savenx = nxtirr(lst2)
	  nxtirr(lst2) = 0
	  savear = ar(lst2)
	  ar(lst2) = hx2*hy2
	  if (lcheck .gt. 1 .and. intrat(lcheck-1) .eq. 2) then
	    call setir2(alloc(locirc),mi2tot,mj2tot,
     1             mptr,lst2)
	  else
	     call setir21(alloc(locirc),mi2tot,mj2tot,
     1                    alloc(node(14,mptr)),mitot,mjtot,
     2                    lwidth,lst2)
	  endif
c
          evol = evol + (mi2p1 - 2)*(mj2p1 - 2)
	  xlow = rnode(1,mptr) - lwidth*hx2
	  ylow = rnode(2,mptr) - lwidth*hy2
          time   = rnode(12,mptr)
          dt     = rnode(11,mptr)
          tpre   = time - dt
          call method(alloc(locbgc),alloc(locsc1),alloc(locsc5),
     1                alloc(locirc),mi2tot,mj2tot,lwidth,
     2                dt2,dtnew2,lst2,hx2,hy2,rflag,
     3                iorder,xlow,ylow,mptr,vtime,steady,
     4                alloc(locqx),alloc(locqy),level,difmax,
     5                lastout,nvar,tpre)
	  nxtirr(lst2) = savenx
	  ar(lst2)   = savear
c
          call reclam(locsc1,mi2tot*mj2tot*nvar)
          call reclam(locsc5,mi2tot*mj2tot*nvar)
          call reclam(locqx, mi2tot*mj2tot*nvar)
          call reclam(locqy, mi2tot*mj2tot*nvar)
c
c  the one giant step based on old values is done. now take
c  one regular step based on new values. 
c
      locsc1 = igetsp(mitot*mjtot*nvar)
      locsc5 = igetsp(mitot*mjtot*nvar)
      locqx  = igetsp(mitot*mjtot*nvar)
      locqy  = igetsp(mitot*mjtot*nvar)
c
      evol   = evol + (maxip1 - 2)*(maxjp1 - 2)
      xlow   = rnode(1,mptr) - lwidth*hx
      ylow   = rnode(2,mptr) - lwidth*hy
      locbig = node(2,mptr)
      loctmp = node(8,mptr)
c     call method(alloc(locbig),alloc(locsc1),alloc(locsc5),
c    1            alloc(node(14,mptr)),mitot,mjtot,lwidth,
c    2            dt,dtnew,node(17,mptr),hx,hy,.true.,
c    3            iorder,xlow,ylow,mptr,vtime,steady,
c    4            alloc(locqx),alloc(locqy),level,difmax,lastout)
c     call postm(alloc(locbig),alloc(loctmp),nvar,maxip1,maxjp1,
c    1               mitot,mjtot,lwidth,nvar)
c
      call reclam(locqx, mitot*mjtot*nvar)
      call reclam(locqy, mitot*mjtot*nvar)
      call reclam(locsc1,mitot*mjtot*nvar)
      call reclam(locsc5,mitot*mjtot*nvar)
c
c need 2 errfs since rbuff list won't parallelize easily.
c 
      call errf1(alloc(locbig),alloc(node(14,mptr)),
     1          nvar,maxi,maxj,maxip1,maxjp1,
     1          alloc(locbgc),mptr,
     2          alloc(locirc),mi2tot,mj2tot,node(17,mptr),
     3          mitot,mjtot,alloc(locbig))
      call reclam(locirc,mi2tot*mj2tot)
      call reclam(locbgc,mi2tot*mj2tot*nvar)
      do 42 i = 1, maxip1*maxjp1*nvar
 42   alloc(loctmp+i-1) = alloc(locbig+i-1)
      call reclam(locbig,mitot*mjtot*nvar)
c
      mptr = node(10, mptr)
      if (mptr .ne. 0) go to 25
c
c final loop needed for the part of errf1 that didn't parallelize
c
      mptr = lstart(lcheck)
 45   continue
         maxi   = node(5, mptr)
         maxj   = node(6, mptr)
         maxip1 = maxi + 1
	 maxjp1 = maxj + 1
	 mitot  = maxi + 2*lwidth - 1
	 mjtot  = maxj + 2*lwidth - 1
	 loctmp = node(8, mptr)
	 numflg = 0
         call errf2(alloc(loctmp),nvar,maxi,maxj,maxip1,maxjp1,
     1              numflg,mptr,rbuff,nbuff,iflow,
     2              alloc(node(14,mptr)),mitot,mjtot)
	 numbad = numbad + numflg
      mptr = node(10,mptr)
      if (mptr .ne. 0) go to 45
c
      return
      end
