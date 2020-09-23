ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c common macros for all routines
c define(good,0.0)
c define(goodpt,0.0)
c define(badpt,2.0)
c define(bad,2.0)
c define(yes,1)
c define(no,0)
c define(macheps,.0001)
c define(rinfinity, 10.e32)
c define(iinfinity, 1000000000)
c define(maxlv,12) - max. permitted levels of nesting
c define(maxcls,192) - max. number of clusters (corner dimensions)
c define(maxgr,192) - max. number of grids (node dimensions)
c define(rsize,12)
c define(nsize,17)
c define(maxstorage,allocsize)    allocsize set in calloc.i
 
c macros for the nodes of the grid tree
c for node (integer part of record)
c define(parent,1)
c define(bndryptr,2) - for perimeter pieces - phasing out, use as
c                      temp. storage for enlarged grids in adv. and ee.
c define(owner,3) - for use in multiprocessing version
c define(nestlevel,4)
c define(maxnumrow,5)
c define(maxnumcol,6)
c define(store1,7)
c define(store2,8)
c define(-------,9) - not used now
c define(levelptr,10)
c define(errptr,11)
c define(prevlevel,12)
c define(cbndryptr,13) - temporary storage for for coarsened perimeter 
c                       pieces and errest temp.  storage
c define(permstore,14) - pointer to irreg array for each grid
c define(cfluxptr,15)
c define(ffluxptr,16)
c define (lstptr,17) - pointer to start of irregular list for that grid.
 
c for rnode (real part of record)
c define(corn1x,1)
c define(corn1y,2)
c define(corn2x,3)
c define(corn2y,4)
c define(corn3x,5)
c define(corn3y,6)
c define(corn4x,7)
c define(corn4y,8)
c define(hrow,9)
c define(hcol,10)
c define(ktime,11)
c define(timemult,12)
 
c tree - related macros:
c define(nextfree,2)
c define(null, 0)
c define(noinsect, 1)
c define(nil,  0)
c
c
c --------------------------------------------------------------
c
      subroutine advanc (level,nvar,dtlev,vtime,steady,lastout)
c
      implicit double precision (a-h,o-z)
      logical    rflag,vtime,steady
      logical    graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder, mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      common   /stats/   evol, rvol, rvoll(maxlv), lentot, lenmax
      include "calloc.i"
      data    rflag/.true./
      logical first/.false./, lastout
c
c  fill all boundary pieces of grids at this level. done
c  in this order to (a) help parallelize and (b) don't need
c  2 copies of finest level storage, since only needed
c  to fill boundary pieces for grids not yet advanced. 
c  integrate all grids at 'level' by one step of its delta(t)
c
      mptr = lstart(level)
 3    continue
          maxi   = node(5,mptr)
          maxj   = node(6,mptr)
          maxip1 = maxi + 1
          maxjp1 = maxj + 1
          mitot  = maxi-1 + 2*lwidth
          mjtot  = maxj-1 + 2*lwidth
          locbig = igetsp(mitot*mjtot*nvar)
          node(2,mptr) = locbig
          locnew = node(7,mptr)
          xl   = rnode(1, mptr)
          yb   = rnode(2, mptr)
          xr   = rnode(7, mptr)
          yt   = rnode(4, mptr)
          hx   = rnode(9, mptr)
          hy   = rnode(10,mptr)
          time = rnode(12,mptr)
          intime = intcnt(level)
          call prem(alloc(locbig),alloc(locnew),nvar,maxip1,maxjp1,
     1              mitot,mjtot,lwidth)

c  #### read cart3d solution for restart
         if (first) then
            call readcart(alloc(locbig),mitot,mjtot)
         endif
          call bound(xl,xr,yb,yt,hx,hy,time,intime,level,nvar,lwidth,
     1               alloc(locbig),mitot,mjtot,alloc(node(14,mptr)),
     2               node(17,mptr))
        mptr = node(10, mptr)
        if (mptr .ne. 0) go to 3
c
      dtlev = 1.e32
      mptr  = lstart(level)
 5    continue
          locold = node(8, mptr)
          locnew = node(7, mptr)
          maxi   = node(5, mptr)
          maxj   = node(6, mptr)
          maxip1 = maxi + 1
          maxjp1 = maxj + 1
          delt   = possk(level)
          rnode(11,mptr) = delt
          hx     = rnode(9,mptr)
          hy     = rnode(10,mptr)
          time   = rnode(12,mptr)
c
          mitot  = maxi-1 + 2*lwidth
          mjtot  = maxj-1 + 2*lwidth
          locbig = node(2, mptr)
          locsc1 = igetsp(mitot*mjtot*nvar)
          locsc5 = igetsp(mitot*mjtot*nvar)
          locqx  = igetsp(mitot*mjtot*nvar)
          locqy  = igetsp(mitot*mjtot*nvar)
c
c  copy old soln. values into  next time steps soln. values
c  since integrator will overwrite it. only for grids not at
c  the finest level. finest level grids do not maintain copies
c  of old and new time solution values.
c
          if (level .lt. mxnest) then
             ntot   = maxip1 * maxjp1 * nvar
cdir$ ivdep
             do 10 i = 1, ntot
 10          alloc(locold + i - 1) = alloc(locnew + i - 1)
          endif
c
      xlow = rnode(1,mptr) - lwidth*hx
      ylow = rnode(2,mptr) - lwidth*hy
      rvol = rvol + (maxip1 - 2)*(maxjp1 - 2)
      rvoll(level) = rvoll(level) + (maxip1 - 2)*(maxjp1 - 2)
      locirr = node(14,mptr)
      locncount = locirr + mitot*mjtot
      locnumHoods = locncount + mitot*mjtot
c
      call method(alloc(locbig),alloc(locsc1),alloc(locsc5),
     1            alloc(locirr),mitot,mjtot,lwidth,
     2            delt,dtnew,node(17,mptr),hx,hy,rflag,
     3            iorder,xlow,ylow,mptr,vtime,steady,
     4            alloc(locqx),alloc(locqy),level,difmax,
     5            lastout,nvar,time)
       call postm(alloc(locbig),alloc(locnew),nvar,maxip1,maxjp1,
     1            mitot,mjtot,lwidth)
      if (node(15,mptr) .ne. 0)
     1   call fluxsv(mptr,alloc(locsc1),alloc(locsc5),
     2               alloc(node(15,mptr)),mitot,mjtot,
     3               nvar,listsp(level))
      if (node(16,mptr) .ne. 0) then
         lenbc = 2*((maxi-1)/intrat(level-1)+(maxj-1)/intrat(level-1))
         call fluxad(alloc(locsc1),alloc(locsc5),
     1               alloc(node(16,mptr)),mptr,mitot,mjtot,maxi,maxj,
     2               nvar,lenbc,intrat(level-1))
      endif
c
          call reclam(locsc1,mitot*mjtot*nvar)
          call reclam(locsc5,mitot*mjtot*nvar)
          call reclam(locqx, mitot*mjtot*nvar)
          call reclam(locqy, mitot*mjtot*nvar)
          call reclam(locbig,mitot*mjtot*nvar)
c
          if (vtime) then
             dtlev = dmin1(dtlev,dtnew)
          endif
            
c
          rnode(12,mptr)  = rnode(12,mptr)+delt
          mptr            = node(10, mptr)
          if (mptr .ne. 0) go to 5
c
      return
      end
c
c    
      subroutine readcart(val,mitot,mjtot)
      implicit double precision(a-h,o-z)
      dimension val(mitot,mjtot,4)
      dimension q(4)

10    continue
      read(8,100,end=99) i,j,(q(k),k=1,4)
 100  format(23x,i3,3x,i3,4x,4e16.7)
c     write(6,*)i,j,(q(k),k=1,4)
c     q is in primitive vars
      val(i,j,1) = q(1)
      val(i,j,2) = q(1)*q(2)
      val(i,j,3) = q(1)*q(3)
      val(i,j,4) = q(4)/.4 + .5d0*q(1)*(q(2)*q(2)+q(3)*q(3))


      go to 10 

99    return
      end
