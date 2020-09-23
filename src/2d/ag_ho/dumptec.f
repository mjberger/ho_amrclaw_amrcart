c
c -----------------------------------------------------
c
      subroutine dumptec (lst,lend,nvar,nplot,time)
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
      character*23  filename , filename2
c
c dumptec = make tecplot file for finest level soln values over entire domain
c (but for now assume only 1 level  -- output finest level)
c
c nplot == plotnum (0 for initial conditions, 1 (or more) for final conditions
c
      filename = 'disjointCutPlanesxx.dat'
      filename(18:18) = '0'
      filename(19:19) = char(ichar('0')+nplot)
      filename2 = 'Uxx.pos'
      filename2(2:2) = '0'
      filename2(3:3) = char(ichar('0')+nplot)
      nplot = nplot+1

      open(14,file=filename,status='unknown',form='formatted')
      open(231,file=filename2,status='unknown',form='formatted')
      write(*,*)" writing graphics file ",filename," at time ",time
      write(14,100) 
 100  format('TITLE = "Extracted cutting Planes through mesh"' )

10    level = lfine
          mptr = lstart(level)
20        if (mptr .eq. 0) go to 50
              maxi   = node(5,mptr)
              maxj   = node(6,mptr)
              xlow = rnode(1,mptr)
              ylow = rnode(2,mptr)
              maxip1 = maxi + 1
              maxjp1 = maxj + 1
              mitot = maxi-1+2*lwidth
              mjtot = maxj-1+2*lwidth
              locnew = node(7,mptr)
              locbig = igetsp(mitot*mjtot*nvar)
              xl   = rnode(1, mptr)
              yb   = rnode(2, mptr)
              xr   = rnode(7, mptr)
              yt   = rnode(4, mptr)
              hx   = rnode(9, mptr)
              hy   = rnode(10,mptr)
              intime = intcnt(level)

             call prem(alloc(locbig),alloc(locnew),nvar,maxip1,maxjp1,
     1              mitot,mjtot,lwidth)
             call bound(xl,xr,yb,yt,hx,hy,time,intime,level,nvar,lwidth,
     1               alloc(locbig),mitot,mjtot,alloc(node(14,mptr)),
     2               node(17,mptr))

              ! remember got 3 times for irr to include other arrays
              locirr = node(14,mptr)
              locncount = locirr + mitot*mjtot
              locnumHoods = locncount + mitot*mjtot
              locqx  = igetsp(mitot*mjtot*nvar)
              locqy  = igetsp(mitot*mjtot*nvar)

              call outtec(alloc(locbig),nvar,mptr,
     1                    alloc(locirr),mitot,mjtot,
     2                    alloc(locqx),alloc(locqy),lwidth,
     2                    node(17,mptr),hx,
     3                    hy,xlow,ylow,time)
              call outgmsh(alloc(locbig),nvar,mptr,
     1                    alloc(locirr),mitot,mjtot,
     2                    alloc(locqx),alloc(locqy),lwidth,
     2                    node(17,mptr),hx,
     3                    hy,xlow,ylow,time)
              call reclam(locqx,mitot*mjtot*nvar)
              call reclam(locqy,mitot*mjtot*nvar)
              call reclam(locbig,mitot*mjtot*nvar)
c
              mptr = node(10,mptr)
          go to 20
50        continue
c
      close(14)
c
 99   return
      end
