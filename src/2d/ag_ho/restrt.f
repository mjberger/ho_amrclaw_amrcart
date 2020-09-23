c
c ---------------------------------------------------------
c
      subroutine restrt(nsteps,time,nvar,lfix)
c
      implicit double precision (a-h,o-z)
      logical            graf,graf1,ee
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
      include "cirr.i"
      common   /space/  lfree(150,2),lenf,idimf
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      common  /stats/  evol, rvol, rvoll(maxlv), lentot, lenmax
      common  /cloops/  xloops(10),yloops(10),nloops
      dimension        intrt1(maxlv)
c
c alloc array could have been written with a smaller size at checkpoint
c 
      read(9) lenmax, lentot, isize

c      ichunk = 500000
c      ist = 1
c5    continue
c     iend = min(ist + ichunk-1, lenmax)
c     read(9)  (alloc(i),i=ist,iend)
c     ist = iend+1
c     if (iend .lt. lenmax) go to 5

c     read(9) (alloc(i),i=1,lenmax)

      read(9) alloc

      read(9) hxposs,hyposs,possk,icheck,intcnt
      read(9) lfree,lenf,idimf
      read(9) rnode,node,lstart,newstl,listsp,tl,
     1  bzonex,bzoney,mstart,ndfree,lfine,iorder,mxnold,intrt1,
     2        kcheck1,maxcl,graf1,lhead,nsteps,time
      read(9) cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,ismp,lfix
      read(9) evol, rvol, rvoll
      read(9) poly,ar,xcirr,ycirr,points
      read(9) wt,ix,iy,nxtirr
      read(9) xloops,yloops,nloops
      write(6,100) nsteps,time
 100  format(' restarting the calculation after ',i5,' steps',
     1        /,'  (time = ',e15.7,')')
c
c adjust free list of storage in case size has changed.
c
      idif = allocsize - isize
      if (idif .gt. 0) then
          lfree(lenf,1) = allocsize + 2
	  call reclam(isize+1,idif)
      else if (idif .lt. 0) then
	    write(6,900) isize,allocsize
 900        format(' size of alloc not allowed to shrink with restart ',
     .             /,' old size ',i10,' current size',i10)
	    stop
      endif
c
c adjust storage in case mxnest has changed - only allow it to increase,
c and only at non-multiples of error estimation on old mxnest.
c
       if (mxnest .eq. mxnold) go to 99

       if (mxnest .lt. mxnold) then
         if (lfine .lt. mxnest) then
             go to 99
         else
             write(6,901) mxnold, mxnest
901          format(' only allow mxnest to increase: ',/,
     &            '  old mxnest ',i4, ' new mxnest ',i4)
             stop
	 endif
       endif

c      see if simple enough situation to allow changing mxnest
	ee = .false.
	do 10 level = 1, mxnold
	   if (icheck(level) .ge. kcheck) then
	      ee = .true.
	      kmust = icheck(level)
	   endif
10      continue
	if (ee) then
	   write(6,902) mxnold, mxnest, kmust
902        format(/,' only allow changes in mxnest (from ',
     &                i4,' to ',i4,')',/,
     &            ' when not time to error estimate: ',/,
     &            ' please run a few more steps before changing ',/,
     &            ' so that # of steps not greater then kcheck',/,
     &            ' or make kcheck > ',i4 )
	     stop
	else
c          #  add second storage location to previous mxnest level
	   mptr = lstart(mxnold)
15	   if (mptr .eq. 0) go to 25 
	      maxip1 = node(5,mptr)+1
	      maxjp1 = node(6,mptr)+1
	      node(8,mptr) = igetsp(maxip1*maxjp1*nvar)
	      mptr = node(10,mptr)
	      go to 15
25         continue
	endif
c
c          # add new info. to spatial and counting arrays
 99	   level = lfine + 1
	   rr = dfloat(intrat(lfine))
35         if (level .gt. mxnest) go to 45
	     hxposs(level) = hxposs(level-1) / rr
	     hyposs(level) = hyposs(level-1) / rr
	     possk (level) = possk (level-1) / rr
	     rr            = intrat(level)
	     level         = level + 1
	     go to 35
45         continue
c
c
      return
      end
