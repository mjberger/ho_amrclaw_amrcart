c
c --------------------------------------------------------------
c
      subroutine errf1(rctold,irr,nvar,maxi,maxj,maxip1,maxjp1,
     1                 rcterr,mptr,irr2,mi2tot,mj2tot,kreg,
     2                 mitot,mjtot,rctsm)
c
c  compare error estimates in rctold, rcterr. if exceed tol, flag.
c  We put in a hook which allows us to ignore the error estimates
c  and not refine if we are outside some prescribed region.
c  rcterr and rctold are enlarged grids.
 
      implicit double precision (a-h,o-z)
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
      include "cirr.i"
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
 
      dimension rctold(mitot,mjtot,nvar),   irr(mitot,mjtot)
      dimension rcterr(mi2tot,mj2tot,nvar), irr2(mi2tot,mj2tot)
      dimension rctsm(maxip1,maxjp1,nvar)
      logical   eprint,  edebug
      data      eprint /.false./, edebug/.false./, nres/60/
c
      time  = rnode(12, mptr)
      xleft = rnode(1,mptr)
      levm  = node(4, mptr)
      hx    = hxposs(levm)
      ybot  = rnode(2,mptr)
      hy    = hyposs(levm)
      dt    = rnode(11,mptr)
 
      errmax = 0.0d0
      err2   = 0.0d0
c     order  = dt*dfloat(2**(iorder+1) - 2)
      order  = dfloat(2**(iorder+1) - 2)
c
      if (.not. (edebug)) go to 20
         write(6,107)
 107     format(//,' coarsened grid residuals')
         do 10 jj = lwidth+1, mj2tot-lwidth
            j = mj2tot - jj + 1
            write(6,101) (rcterr(i,j,1),i = lwidth+1, mi2tot-lwidth)
10       continue
101      format(' ',10e8.1)
c
c zero out the exterior locations so they don't affect err.est.
c
 20   continue
      ar(-1) = 0.d0
      jfine = lwidth+1
      do 35  j = lwidth+1, mj2tot-lwidth
c     yofj  = ybot + (dfloat(jfine - 1) - .5)*hy
      ifine = lwidth+1
c
      do 30  i  = lwidth+1, mi2tot-lwidth
          rflag = 0.0
          xofi  = xleft + (dfloat(ifine - 1) - .5)*hx
	  est = dabs(rcterr(i,j,1))
c	  if (irr2(i,j) .ne. kreg) est = 0.d0
          if (est .gt. errmax) errmax = est
	  err2 = err2 + est*est
c         write(6,102) i,j,est
 102      format(' i,j,est ',2i5,e12.5)
c         rcterr(i,j,2) = est
          if (est .ge. tol)   rflag  = 2.0
          rcterr(i,j,1) = rflag
          ifine = ifine + 2
 30   continue
      jfine = jfine + 2
 35   continue
c
c
c  transfer flagged points on cell centered coarse grid
c  to cell centered fine grid. count flagged points.
c
c  (re) initialize rctold to 0.0 before flagging
c
      do 40 i = 1, mitot
      do 40 j = 1, mjtot
 40      rctold(i,j,1) = 0.0
c
c  print out intermediate flagged rcterr (for debugging)
c
      if (eprint) then
	 err2 = dsqrt(err2/dfloat((mi2tot-2*lwidth)*(mj2tot-2*lwidth)))
         write(6,103) mptr, levm, errmax, err2
 103     format(' grid ',i4,' level ',i4,
     .          ' max. error = ',e15.7,' err2 = ',e15.7,/,
     1          ' flagged points on coarsened grid ')
         do 45 jj = lwidth+1, mj2tot-lwidth
         j = mj2tot - jj + 1
         write(6,106) (nint(rcterr(i,j,1)),i=lwidth+1,mi2tot-lwidth)
106      format(1h ,80i1)
45       continue
      endif
c
c rctsm is the same space as rctold, but dimensioned as the smaller
c grid, so that can be copied into difference storage in errest
c
      jfine   = 2
      do 70 j = lwidth+1, mj2tot-lwidth
      ifine   = 2
      do 60 i = lwidth+1, mi2tot-lwidth
      if (rcterr(i,j,1) .eq. 0.0d0) go to 55
      if (levm .eq. 1) then
         if (irr(ifine+lwidth-1,jfine+lwidth-1).ne.-1) 
     .       rctsm(ifine,jfine,1)=2.d0
         if (irr(ifine+lwidth,jfine+lwidth-1).ne.-1) 
     .       rctsm(ifine+1,jfine,1)=2.d0
         if (irr(ifine+lwidth-1,jfine+lwidth).ne.-1) 
     .       rctsm(ifine,jfine+1,1)=2.d0
         if (irr(ifine+lwidth,jfine+lwidth).ne.-1) 
     .       rctsm(ifine+1,jfine+1,1)=2.d0
 	 go to 55
      else
         rctsm(ifine,jfine,1)     = 2.0d0
         rctsm(ifine+1,jfine,1)   = 2.0d0
         rctsm(ifine,jfine+1,1)   = 2.0d0
         rctsm(ifine+1,jfine+1,1) = 2.0d0
      endif
 55   ifine   = ifine + 2
 60   continue
      jfine   = jfine + 2
 70   continue
c
      if (eprint) then
         write(6,105) mptr
 105     format(' grid ',i4,' before res flagging ',/)
         do 47 jj = 2, maxj
         j = maxjp1 - jj + 1
         write(6,106) (nint(rctsm(i,j,1)),i=2,maxi)
47       continue
      endif

c     refine if not enough points around body (put curvature here too)
      nresadd = 0
      write(6,*) ' checking level ',levm
      if (levm .lt. 9 .and. iprob .eq. 6) then ! for naca airfoil
      do 38 i = 2, maxi
      do 38 j = 2, maxj
	ifine = i + lwidth - 1
	jfine = j + lwidth - 1
        k = irr(ifine,jfine)
        if (k .ne.kreg .and. k .ne.-1) then
c          find length of arc of solid body
           do 39 kside = 1, 6
             if (poly(kside+2,1,k) .eq. -11) then
               sidex2 = poly(kside,1,k)
               sidey2 = poly(kside,2,k)
               sidex1 = poly(kside+1,1,k)
               sidey1 = poly(kside+1,2,k)
               go to 41
             endif
 39        continue
   41      rlen = dsqrt((sidey1-sidey2)**2 + (sidex1-sidex2)**2)
c try always refining around body
	   rctsm(i,j,1) = 2.0
           if (rlen .gt. .1/nres) then
	      rctsm(i,j,1) = 2.0
	      nresadd = nresadd + 1
	   endif
        endif
 38   continue
      endif

      if (iprob .eq. 2) then  ! flat plate, refine at geom
      nresadd = 0
      do 49 i = 2, maxi
      do 48 j = 2, maxj
      	ifine = i + lwidth - 1
	jfine = j + lwidth - 1
        k = irr(ifine,jfine)
        if (k .ne.kreg .and. k .ne.-1) then
	   numtoflag = 12    ! flag more on coarser grids
	   do jfl = j,min(j+numtoflag,maxj)  ! found first cut cell
              rctsm(i,jfl,1) = 2.0  ! flag a few above it
           nresadd = nresadd + 1
	   end do
	   go to 49   ! for flat plate, no more geometry in this row
        endif
   
 48   continue
 49   continue
      endif


      if (eprint) then
         write(6,104) mptr, nresadd
 104     format(' grid ',i4,' after res flagging of',i4,' pts. ',/)
         do 46 jj = 2, maxj
         j = maxjp1 - jj + 1
         write(6,106) (nint(rctsm(i,j,1)),i=2,maxi)
46       continue
      endif
c

      return
      end
