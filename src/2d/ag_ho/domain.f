c
c  ----------------------------------------------------------
c
      subroutine domain (nvar,lfix,vtime,quad)
c
      implicit double precision (a-h,o-z)
      logical    graf, vtime, ok, quad
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "cirr.i"
      include "calloc.i"
      common /userdt/cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .               ismp,gradThreshold
c
c  allocate initial coarse grid domain. set node info & initialize grid
c  initial space and time step set here too
c
      mstart = nodget(dummy)
c
c code assumes in many places that lower left corner at (0,0)
c this initial code sets the domain - assumed rectangular
c if it is too large, birect will chop it up into several rectangular
c pieces
c
      rnode(1,mstart)= 0.0
      rnode(2,mstart)= 0.0
      rnode(3,mstart)= 0.0
      rnode(4,mstart)= yprob
      rnode(5,mstart)= xprob
      rnode(6,mstart)= yprob
      rnode(7,mstart)= xprob
      rnode(8,mstart)= 0.0
      node(4,mstart) = 1
      node(10,mstart)= 0
      lstart(1) = mstart
      rnode(12,mstart) = 0.d0
      call sethk(mstart)
      if (mxnest .eq. 1) then
	  lfix = 1
	  go to 35
      endif
c
c     ymid = 12.319431966
c     ymid = 11.121709413750
      ymid = yprob/2.d0
      mnext = nodget(dummy)
      rnode(1,mnext) = 0.
      rnode(2,mnext) = 0.
      rnode(3,mnext) = 0.
      rnode(4,mnext) = ymid
      rnode(5,mnext) = xprob
      rnode(6,mnext) = ymid
      rnode(7,mnext) = xprob
      rnode(8,mnext) = 0.
      node(4,mnext)  = 2
      node(10,mnext) = 0
      lstart(2) = mnext
      call sethk(mnext)
      if (mxnest .eq. 2) then
        lfix = 2
	go to 35
      endif
c
c     ymid =  8.21295464400
      ymid =  yprob/4.d0
      mnext = nodget(dummy)
      rnode(1,mnext) = 0.
      rnode(2,mnext) = 0.
      rnode(3,mnext) = 0.
      rnode(4,mnext) = ymid
      rnode(5,mnext) = xprob
      rnode(6,mnext) = ymid
      rnode(7,mnext) = xprob
      rnode(8,mnext) = 0.
      node(4,mnext)  = 3
      node(10,mnext) = 0
      lstart(3) = mnext
      call sethk(mnext)
      lfix = 3
      go to 35
c
      ymid = 2.05323866100
      mnext = nodget(dummy)
      rnode(1,mnext) = 0.
      rnode(2,mnext) = 0.
      rnode(3,mnext) = 0.
      rnode(4,mnext) = ymid
      rnode(5,mnext) = xprob/2.
      rnode(6,mnext) = ymid
      rnode(7,mnext) = xprob/2.
      rnode(8,mnext) = 0
      node(4,mnext) = 4
      node(10,mnext)=0
      lstart(4) = mnext
      mold = mnext

      ymid = 1.7110322175
      ytop = 4.106477322000
      mnext = nodget(dummy)
      rnode(10,mold) = mnext
      rnode(1,mnext) = xprob/2.
      rnode(2,mnext) = ymid
      rnode(3,mnext) = xprob/2.
      rnode(4,mnext) = ytop
      rnode(5,mnext) = xprob
      rnode(6,mnext) = ytop
      rnode(7,mnext) = xprob
      rnode(8,mnext) = ymid
      node(4,mnext) = 4
      node(10,mnext)=0
      lstart(4) = mnext

      call sethk(lstart(4))
c
c
c  check that user specified stepsize really does divide user specified domain
c  and that number of points is odd in both directions.
c
 34    lfix = mxnest
 35    continue
       do 30 lcheck = 1, lfix
         mcheck = lstart(lcheck) 
         ok = .true.
         if (dabs(rnode(1,mcheck)+hxposs(lcheck)*
     1       dfloat(node(5,mcheck)-1)-rnode(7,mcheck)).gt..0001) 
     2       ok = .false.
         if (dabs(rnode(2,mcheck)+hyposs(lcheck)*
     1       dfloat(node(6,mcheck)-1)-rnode(4,mcheck)).gt..0001) 
     2       ok = .false.
         if (ok) go to 30
            write(6,100)
 100        format('  user specified stepsize doesn''t divide domain ')
            write(6,102) xprob,yprob
 102        format("domain edges are ",2e25.15)
            computedXedge = dabs(rnode(1,mcheck)+hxposs(lcheck)*
     1       dfloat(node(5,mcheck)-1)-rnode(7,mcheck))
            computedYEdge = dabs(rnode(2,mcheck)+hyposs(lcheck)*
     1       dfloat(node(6,mcheck)-1)-rnode(4,mcheck))
            write(6,103) computedXEdge,computedYedge
 103        format("computed edges are ",2e25.15)
            stop
 30   continue
c
c  for this code, must have odd number gridpoints
c
      do 40 lcheck = 1, lfix
         mcheck = lstart(lcheck)
         ok = .true.
         if (.not. ( ((node(5,mcheck)/2)*2 .ne. node(5,mcheck)) .and.
     1        ((node(6,mcheck)/2)*2 .ne. node(6,mcheck)) .and.
     2        ((intrat(lcheck)/2)*2 .eq. intrat(lcheck)) )) ok = .false.

      ! commented out starting here - AG
!         if (ok) go to 40
!           write(6,101)
! 101       format(' mesh widths must give odd number grid points '
!     1           /,' and even refinement ratio')
!           stop
!     commented out ending here - AG

 40   continue
c
      do 50 lcheck = 1, lfix
         mcheck = lstart(lcheck)
         call birect(mcheck)
         call sethk(mcheck)
         call ginit (mcheck, .true.,nvar,quad,gradThreshold)
 50   continue
c
      lfine = lfix
      call join(0,nvar)
c
c  set stable initial time step using coarse grid data
c
      if (vtime) then
         dtgrid = 1.e+20
         mptr = lstart(1)
         dx = hxposs(1)
         dy = hyposs(1)
 60         maxip1 = node(5,mptr) + 1
            maxjp1 = node(6,mptr) + 1
	          mitot  = maxip1 - 2 + 2*lwidth
	          mjtot  = maxjp1 - 2 + 2*lwidth
	          call estdt(alloc(node(7,mptr)),maxip1,maxjp1,nvar,
     1                       dx,dy,dt,alloc(node(14,mptr)),mitot,mjtot,
     2                       lwidth)
            dtgrid = dmin1(dt,dtgrid)
            mptr = node(10,mptr)
            if (mptr .ne. 0) go to 60
         possk(1) = dtgrid
      endif
c
c set rest of possk array for refined timesteps
c
      rr = intrat(1)
      do 70 level = 2, mxnest
         possk(level) = possk(level-1)/rr
c        possk(level) = possk(level-1)  ! for steady calcs pretend same timestep
         rr = intrat(level)
 70   continue
c
      return
      end
