c
c ---------------------------------------------------------------
c
        subroutine filpc(xlp,xrp,ybp,ytp,level,nvar,valdub,time,
     &                      intime,midub,mjdub,nrowst,nrowend,ncolst,
     &                      ncolend,irr,mitot,mjtot)
 
c This routine calculates the values for a patch at time time, 
c integer offset time intime, and level level,
c with lower left and upper right corners located at (xlp,ybp), (xrp,ytp),
c and inserting them in the enlarged array valdub. It first fills in
c those values in the interior which are obtainable from the level level
c grids, and then recusively obtaining the remaining values from
c coarser levels. This routine calls the following routines:
c
c intfil: finds all of the available values of the solution at a given
c         level contained in a patch, copies them and marks a used array;
c
c physbd: sets all the values in a patch corresponding to values exterior
c         to the physical domain;
c
c trimbd: given a used array, finds the limits of the smallest subarray
c         containing all the points which have not been marked in the
c         used array.
c
c The implementation of this routine is done in two steps. First, we
c go down in levels starting from the level at which patch is being
c called: at each level we find out which points can be filled from that level
c (using intfil) and calculate the smallest rectangle which will contain
c the points required to interpolate the remaining values. The second step is
c to perform the interpolations, beginning form the bottom up. The temporary
c storage required for this routine is allocated from alloc: the pointers
c to the used arrays marking which points have already been filled from a given
c level are contained in locuse; the pointers to the values themselves are
c contained in locvar.
 
      implicit double precision (a-h,o-z)
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
      dimension ialloc(2*allocsize)
      equivalence (alloc,ialloc)
 
        logical flag, allin, reset, outside
        dimension valdub(midub,mjdub,nvar), irr(mitot,mjtot)
        dimension xl(maxlv),xr(maxlv),yb(maxlv),yt(maxlv),
     &            locuse(maxlv),locvar(maxlv),locirr(maxlv),il(maxlv),
     &            ir(maxlv),jb(maxlv),jt(maxlv),nrow(maxlv),ncol(maxlv)
 
        ival(ii,jj,ll,nn) = locvar(ll) + (ii - 1) + nrow(ll)*(jj - 1)
     &                      + nrow(ll)*ncol(ll)*(nn - 1)
c       int2(xx) = int(xx) - int(.5 - sign(.5,xx))
        int2(xx) = idint(xx) - idint(.5d0 - dsign(.5d0,xx))
	iirr(i,j,l) = 2*(locirr(l)-1)+1 + i - 1 + nrow(l)*(j-1)
	outside(i,j) = ((i .lt. 1)     .or. (j .lt. 1) .or.
     &                  (i .gt. mitot) .or. (j .gt. mjtot))
c
c We begin by filling values for grids at level level. If all values can be
c filled in this way, we return; otherwise, we recursively descend through
c as many levels as required to interpolate the values in the patch. We
c know that the process must terminate eventually, since the original patch
c is either contained in a level 1 grid or exterior to the computational domain.
 
	nrowp = nrowend-nrowst+1
	ncolp = ncolend-ncolst+1
        locuse(level) = igetsp((nvar+2)*nrowp*ncolp)
 
	ifill = 1
        call intfil
     &  (xlp,xrp,ybp,ytp,level,nrowst,nrowend,ncolst,ncolend,nvar,
     &  valdub,midub,mjdub,time,intime,ifill,locuse(level))
 
c If level = 1 then return. intfil + physbd will fill any
c level 1 patch. physbd called from bound, after all pieces filled.
 
        if (level.eq.1) then
           call reclam(locuse(level),(nvar+2)*nrowp*ncolp)
           return
        endif
 
c Trimbd returns flag = true if all of the entries of alloc are 0.,
c flag = false, otherwise. If flag = true, then no other levels are
c are required to interpolate, and we return. Note that the used array is
c filled entirely in intfil, i.e. the marking done there also takes
c into account the points filled by the boundary conditions.
c physbd will be called from bound, after all 4 boundary pieces filled.
 
        call trimbd(alloc(locuse(level)),nrowp,ncolp,flag,
     &              il(level),ir(level),jb(level),jt(level))
        if (flag) then
           call reclam(locuse(level),(nvar+2)*nrowp*ncolp)
           return
        endif
 
c flag = false. we will have to interpolate some values from coarser
c levels. We begin by initializing the level level arrays, so that we can use
c purely recursive formulation for interpolating. 
 
        locvar(level) = locuse(level) + nrowp*ncolp
        nrow(level) = nrowp
        ncol(level) = ncolp
 
        levc = level - 1
        levf = level
        xl(level) = xlp
        xr(level) = xrp
        yb(level) = ybp
        yt(level) = ytp
 
c begin loop going down through levels.
 
10      continue
 
        hxc = hxposs(levc)
        hyc = hyposs(levc)
        hxf = hxposs(levf)
        hyf = hyposs(levf)
 
        xlf = xl(levf) + (il(levf) - 1)*hxf + .5*hxf
        xrf = xl(levf) + (ir(levf) - 1)*hxf + .5*hxf
        ybf = yb(levf) + (jb(levf) - 1)*hyf + .5*hyf
        ytf = yb(levf) + (jt(levf) - 1)*hyf + .5*hyf
c 
c WARNING: The following formula uses in an essential way that the absolute
c          origin of the coordinates is (0,0).
c
        xlc = int2(-.5 + ( xlf + .0001*hxf)/hxc)*hxc
        xrc = int2(-.5 + ( xrf - .0001*hxf)/hxc)*hxc+ 2.*hxc
        ybc = int2(-.5 + ( ybf + .0001*hyf)/hyc)*hyc
        ytc = int2(-.5 + ( ytf - .0001*hyf)/hyc)*hyc+ 2.*hyc
c       arg = (-.5 + ( xlf + .0001*hxf)/hxc)
c       xlc = (int(arg) - int(.5 - sign(.5,arg)))*hxc
c       arg = (-.5 + ( xrf - .0001*hxf)/hxc)
c       xrc = (int(arg) - int(.5 - sign(.5,arg)))*hxc+ 2.*hxc
c       arg = (-.5 + ( ybf + .0001*hyf)/hyc)
c       ybc = (int(arg) - int(.5 - sign(.5,arg)))*hyc
c       arg = (-.5 + ( ytf - .0001*hyf)/hyc)
c       ytc = (int(arg) - int(.5 - sign(.5,arg)))*hyc+2.*hyc
 
        xl(levc) = xlc
        xr(levc) = xrc
        yb(levc) = ybc
        yt(levc) = ytc
 
        nrow(levc) = idint((xrc - xlc)/hxc + .0001)
        ncol(levc) = idint((ytc - ybc)/hyc + .0001)
 
        ntot = nrow(levc)*ncol(levc)
        locuse(levc) = igetsp(ntot*(nvar+2))
        locvar(levc) = locuse(levc) + ntot
        locirr(levc) = locvar(levc) + nvar*ntot
 
	ifill = 1
        call intfil
     1  (xl(levc),xr(levc),yb(levc),yt(levc),levc,
     2  1,nrow(levc),1,ncol(levc),nvar,alloc(locvar(levc)),
     3  nrow(levc),ncol(levc),time,intime,ifill,locuse(levc))
 
	call tirr(alloc(locirr(levc)),nrow(levc),ncol(levc),
     &            xl(levc),xr(levc),yb(levc),yt(levc),levc)

        call trimbd(alloc(locuse(levc)),nrow(levc),ncol(levc),flag,
     &              il(levc),ir(levc),jb(levc),jt(levc))
 
        if (flag) then
            call physbd(xl(levc),xr(levc),yb(levc),yt(levc),levc,
     &                  nrow(levc),ncol(levc),nvar,alloc(locvar(levc)),
     &                  time,hxposs(levc),hyposs(levc))
        else
            levf = levc
            levc = levc - 1
            go to 10
        endif
 
c       go back up the levels, interpolating as you go.
 
20      continue
 
        do 100 iff = 1,nrow(levf)
        do 100 jf  = 1,ncol(levf)

                flag =
     &  (alloc(locuse(levf) + iff - 1 + (jf - 1)*nrow(levf)).eq.0.)
 
                if (flag) then
 
                xif = xl(levf) + (.5 + dfloat(iff - 1))*hxposs(levf)
                yjf = yb(levf) + (.5 + dfloat(jf - 1))*hyposs(levf)
 
                ic=idint((xif-xl(levc)+.5*hxposs(levc))/hxposs(levc))
		jc=idint((yjf-yb(levc)+.5*hyposs(levc))/hyposs(levc))
 
                xc = xl(levc) + (dfloat(ic) - .5)*hxposs(levc)
                yc = yb(levc) + (dfloat(jc) - .5)*hyposs(levc)
 
                eta1 = (xif - xc)/hxposs(levc)
                eta2 = (yjf - yc)/hyposs(levc)
c
c coarse irr arrays should always be fully set in tirr, even if
c lwidth = 1 and refine by 2
c
		if ((ialloc(iirr(ic,jc,levc)) .ne. -1) .and.
     .        	    (ialloc(iirr(ic+1,jc,levc)) .ne. -1) .and.
     .       	    (ialloc(iirr(ic,jc+1,levc)) .ne. -1) .and.
     .       	    (ialloc(iirr(ic+1,jc+1,levc)) .ne. -1)) then
		  allin = .true.
		else 
		    allin = .false.
		    if (eta1 .lt. .5d0) then
		      iclose = ic
		    else
		      iclose = ic+1
		    endif
		    if (eta2 .lt. .5d0) then
		      jclose = jc
		    else
		      jclose = jc+1
		    endif
c
c finest level irr array not set from lwidth+1 to 2*lwidth all around grid
c
		    if (ialloc(iirr(iclose,jclose,levc)).ne.-1) go to 30

                    if (levf .eq. level) then
                        ieff = iff + nrowst-1-lwidth
                        jeff = jf  + ncolst-1-lwidth
c                       # if no irr info. available for extra lwidth border,
c                       # just ignore for now, assume the best (=cell is ext.)
   		        if (outside(ieff,jeff)) go to 30
                         if (irr(ieff,jeff) .ne. -1) then
		 	     reset = .true.
			 else
			     reset = .false.
			 endif
		    else 
			 if (ialloc(iirr(iff,jf,levf)) .ne. -1) then
			   reset = .true.
			 else
			   reset = .false.
			 endif
		    endif

		    if (reset) then
c                    # reset iclose,jclose to extrapolate soln
                        if (ialloc(iirr(ic,jc,levc)) .ne. -1) then
			   iclose = ic
			   jclose = jc
		        else if (ialloc(iirr(ic+1,jc,levc)).ne.-1) then
		           iclose = ic+1
			   jclose = jc
		        else if (ialloc(iirr(ic,jc+1,levc)).ne.-1) then
			   iclose = ic
			   jclose = jc+1
		        else if (ialloc(iirr(ic+1,jc+1,levc)).ne.-1)then
			   iclose = ic+1
			   jclose = jc+1
		        else
			   write(6,900)
 900                       format(' can"t find value to extrapolate',
     .                            ' in filpc ')
			   stop
		       endif
		    endif
		
 30                continue
		endif
 
                do 101 ivar = 1,nvar
 
                   valc00 = alloc(ival(ic,jc,levc,ivar))
                   valc10 = alloc(ival(ic+1,jc,levc,ivar))
                   valc01 = alloc(ival(ic,jc+1,levc,ivar))
                   valc11 = alloc(ival(ic+1,jc+1,levc,ivar))
                      
		   if (allin) then
                     valint = (1. - eta2)*
     &               ((1. - eta1)*valc00 + eta1*valc10)
     &               + eta2*((1. - eta1)*valc01 + eta1*valc11)
		   else
		     valint = alloc(ival(iclose,jclose,levc,ivar))
		   endif
 
		   if (levf .lt. level) then
                      alloc(ival(iff,jf,levf,ivar)) = valint
		   else
		      valdub(iff+nrowst-1,jf+ncolst-1,ivar) = valint
		   endif
 
101             continue
 
                endif
 
100     continue
 
        call reclam(locuse(levc),nrow(levc)*ncol(levc)*(nvar+2))
 
 
 
        if (levf .lt. level) then
             call physbd(xl(levf),xr(levf),yb(levf),yt(levf),levf,
     1                  nrow(levf),ncol(levf),nvar, alloc(locvar(levf)),
     2                  time,hxposs(levf),hyposs(levf))
             levc = levf
             levf = levf + 1
             go to 20
        endif
 
        call reclam(locuse(level),nrow(level)*ncol(level)*(nvar + 2))

        return
        end
