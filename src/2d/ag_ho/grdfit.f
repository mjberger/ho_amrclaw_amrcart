c
c ---------------------------------------------------------
c
      subroutine grdfit (lbase,lcheck,nvar,cut,cdist,lfnew,time,
     .                   steady)
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
c
      parameter  (maxcls = 192)
      dimension  corner(12,maxcls)
      integer    numptc(maxcls),  prvptr
      logical    cprint, fit, nestck, steady, newalg
      data       cprint/.false. /, newalg/.true./
c
c  grdfit called by setgrd and regrid to actually fit the new grids
c         on each level. lcheck is the level being error estimated
c         so that lcheck+1 will be the level of the new grids.
c
      call flglvl (nvar, lcheck, nptmax, index, ignore,
     .             lbase,steady)
      npts = nptmax - ignore
      if (npts .eq. 0) go to 99
      hxbase = hxposs(lbase)
      hybase = hyposs(lbase)
      maxcl  = maxcls
c
      levnew    = lcheck + 1
      hx        = hxposs(lcheck)
      hy        = hyposs(lcheck)
      buffx     = .5 * hx
      buffy     = .5 * hy
c
c  initially all flagged points are put in 1 cluster
c
      nclust    = 1
      numptc(1) = npts
      lratio = intrat(lcheck)
c
c -- bisect the clusters till each cluster ok ---
c
      if (newalg) then
c         this assumes lower left corner is (0,0)
	  idim = xprob/hx + .1
	  jdim = yprob/hy + .1
	  lociscr = igetsp(idim)
	  locjscr = igetsp(jdim)
	  call smartbis(alloc(index),npts,cut,buffx,buffy,hx,hy,
     1                  numptc,nclust,lbase,corner,maxcl,
     2                  alloc(lociscr),alloc(locjscr),idim,jdim)
	  call reclam(lociscr,idim)
	  call reclam(locjscr,jdim)
       else
          call bisect(alloc(index),npts,cut,buffx,buffy,hx,hy,
     1                numptc,nclust,lbase,corner,maxcl)
      endif
      if (cprint) write(6,103) nclust
 103  format(' ',i4,' clusters after bisect ')
c
c  merge clusters back together to make bigger ones if possible
c
      if (nclust .eq. 1) go to 50
      if (.not. newalg) then
         isp  =  igetsp(2*npts)
         call merge(alloc(isp),alloc(index),npts,cut,buffx,buffy,hx,hy,
     1              numptc,nclust,corner,lratio,maxcl)
         call reclam(isp, 2*npts)
         if (cprint) write(6,104) nclust
 104     format(' ',i4,' clusters after merge')
      endif
c
c  for each cluster, fit the actual grid, set up some data structures
c
 50   ibase   = 0
      prvptr  = 0
      icl     = 1
c
 70   mnew      = nodget(dummy)
 75   call  moment(rnode(1,mnew),alloc(index+2*ibase),numptc(icl),
     *                 usage,buffx,buffy,hx,hy)
      if (cprint) write(6,100) icl,mnew,usage,numptc(icl)
100   format('cluster ',i5,' new rect.,',i5,
     1       ' usage ',e12.5,' with ',i5,' pts.')
      node(4,mnew)   = levnew
c
c  if new grid doesn't fit in base grid, nestck bisect it
c  and returns 2 clusters where there used to be 1.
c
c  could put in corner usage into nestck and remove moment calling
c  from this segment.  make moment a lower level routine, call from
c  nestck when bisect.
c
      fit = nestck(mnew,lbase,alloc(index+2*ibase),numptc(icl),
     1             numptc,icl,nclust,hxbase,hybase)
      if (.not. fit) go to 75
c
      if (prvptr .ne. 0) node(10,prvptr) = mnew
      if (prvptr .eq. 0) newstl(levnew)  = mnew
      rnode(12,mnew) = time
      prvptr = mnew
      ibase  = ibase + numptc(icl)
      icl = icl + 1
 80   if (icl .le. nclust) go to 70
c
c     call drawrg(time,lcheck,newstl(levnew),
c    1            nclust,numptc,npts,alloc(index))
      call birect(newstl(levnew))
      if (.not. newalg) call trimgr(newstl(levnew),hx,hy)
      call reclam(index, 2*nptmax)
c
c  finish setting fields for new grids and some data structures
c
      call sethk(newstl(levnew))
c
 99   return
      end
