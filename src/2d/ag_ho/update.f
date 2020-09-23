c
c -----------------------------------------------------------
c
      subroutine update (level, nvar, work)
c
      implicit double precision (a-h,o-z)
      dimension work (nvar)
      logical  uprint, ovrlap
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
      dimension ialloc(2*allocsize)
      equivalence(alloc,ialloc)
      include "cirr.i"
      data     uprint/.false./
      iadd(i,j,ivar)  = loc + i - 1 + maxip1*((ivar-1)*maxjp1+j-1)
      iaddf(i,j,ivar) = locf + i - 1 + mip1*((ivar-1)*mjp1+j-1)
      iirreg(i,j)     = 2*(locirr-1)+1 + i - 1 + mitot*(j-1)
      iirrc(i,j)      = 2*(locirc-1)+1 + i - 1 + miget*(j-1)
c
c ******************************************************************
c update - update all grids at level 'level'.
c          this routine assumes cell centered variables.
c          the update is done from 1 level finer meshes under it.
c input parameter:
c    level  - ptr to the only level to be updated. levels coarser than
c             this will be at a diffeent time.
c ******************************************************************
c
      lget = level
      if (uprint) write(6,100) lget
100   format(19h    updating level ,i5)
c
c  grid loop for each level
c
      mptr = lstart(lget)
 20   if (mptr .eq. 0) go to 80
         loc    = node(7,mptr)
	 locirc = node(14,mptr)
         maxi   = node(5,mptr)
         maxj   = node(6,mptr)
         maxip1 = maxi + 1
         maxjp1 = maxj + 1
	 miget  = maxi-1+2*lwidth
	 mjget  = maxj-1+2*lwidth
         xst    = rnode(1,mptr)
         yst    = rnode(2,mptr)
         hx     = rnode(9,mptr)
         hy     = rnode(10,mptr)
         dt     = rnode(11,mptr)
         dtfine = possk(lget+1)
c
         if (node(15,mptr) .eq. 0) go to 25
c         locuse = igetsp(maxip1*maxjp1)
c         call upbnd(alloc(node(15,mptr)),alloc(loc),nvar,
c     1              maxip1,maxjp1,listsp(lget),alloc(locuse),dtfine,
c     2              alloc(node(14,mptr)),miget,mjget)
c         call reclam(locuse,maxip1*maxjp1)
c
c  loop through all intersecting fine grids as source updaters.
c
 25      mkid = lstart(lget+1)
 30      if (.not. ovrlap(mkid,mptr)) go to 75
         mi     = node(5,mkid)
         mj     = node(6,mkid)
         mip1   = mi + 1
         mjp1   = mj + 1
         locf   = node(7,mkid)
	 locirr = node(14,mkid)
	 mitot  = mi-1+2*lwidth
	 mjtot  = mj-1+2*lwidth
         hxkid  = rnode(9,mkid)
         hykid  = rnode(10,mkid)
c
c  calculate starting and ending indices for coarse grid update
c
         ist  = max(idint((rnode(1,mkid)-xst)/hx + 2.1), 2)
         jst  = max(idint((rnode(2,mkid)-yst)/hy + 2.1), 2)
         iend = min(idint((rnode(7,mkid)-xst)/hx + 1.1), maxi)
         jend = min(idint((rnode(4,mkid)-yst)/hy + 1.1), maxj)
         if ((ist .gt. iend) .or. (jst. gt. jend)) go to 75
c
c  calculate starting index for fine grid source pts.
c
         ikid = max(idint((xst-rnode(1,mkid))/hxkid + 2.1), 2)
         jkid = max(idint((yst-rnode(2,mkid))/hykid + 2.1), 2)
         iff  = ikid
         jff   = jkid
 
         do 71 i = ist, iend
         do 70 j = jst, jend
           if (uprint) write(6,101) i,j,mptr,iff,jff,mkid
 101       format(' updating pt. ',2i4,' of grid ',i3,' using ',2i4,
     1            ' of grid ',i4)
c
c  update using intrat fine points in each direction
c
           do 40 ivar = 1, nvar
 40        work(ivar) = 0.d0
	   warea      = 0.d0
c
           do 50 jco  = 1, intrat(lget)
           do 50 ico  = 1, intrat(lget)
           index = ialloc(iirreg(iff+ico-2+lwidth,jff+jco-2+lwidth))
           if (index .eq. -1) go to 50
           farea = ar(index)
           do 45 ivar = 1, nvar
              work(ivar) = work(ivar)+
     1                     farea*alloc(iaddf(iff+ico-1,jff+jco-1,ivar))
 45           continue
	   warea = warea + farea
 50        continue
c
            if (uprint) write(6,102)(alloc(iadd(i,j,ivar)),ivar=1,nvar)
 102        format(' old vals: ',4e12.4)
c
c only update if cell was really in the domain, not exterior. could
c double check here that coarse grid was also exterior, if warea=0
c
	   if (warea .ne. 0.d0) then
               do 65 ivar = 1, nvar
 65            alloc(iadd(i,j,ivar)) = work(ivar)/warea
	   else if (ialloc(iirrc(i+lwidth-1,j+lwidth-1)) .ne. -1) then
	          write(6,900) i,j,iff,jff
 900              format(' update error: fine grid cells ',2i4, 
     1                                  'are exterior ',
     1                   /,'             coarse cell',2i4,' isn"t ')
		  stop
           endif
	       
            if (uprint) write(6,103)(alloc(iadd(i,j,ivar)),ivar=1,nvar)
 103        format(' new vals: ',4e12.4)
c
           jff = jff + intrat(lget)
 70        continue
           iff = iff + intrat(lget)
           jff = jkid
 71        continue
c
 75         mkid = node(10,mkid)
            if (mkid .ne. 0) go to 30
c
            mptr = node(10, mptr)
            go to 20
c
 80       continue
c
 99   return
      end
