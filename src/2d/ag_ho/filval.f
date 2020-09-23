c
c ------------------------------------------------------------------
c
      subroutine filval(val,maxip1,maxjp1,hx,hy,lev,time,intime,
     1                  valc,mxcip1,mxcjp1,
     2                  irr,mitot,mjtot,xleft,xright,ybot,ytop,
     3                  nvar,lratio,dudx,dudy,lwdth,lstgrd,
     4                  irc,mptr,lstnew)
 
      parameter (msize = 1700)
      implicit double precision (a-h,o-z)
      dimension val(maxip1,maxjp1,nvar),valc(mxcip1,mxcjp1,nvar),
     1          irr(mitot,mjtot),dudx(mxcip1),dudy(mxcip1),
     2          irc(mxcip1,mxcjp1)
      dimension qx(msize,msize,nvar)
      dimension qy(msize,msize,nvar)
      logical    graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder, mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "calloc.i"
      include "cirr.i"
c
c create and fill coarser (lev-1) patch with one extra coarse cell all
c around. will interpolate from this patch to grid mptr 
c without needing special boundary code. coarsened irr (=irc) array
c for the coarser patch already set in fill. it is same size as valc,
c not bigger to include lwidth border.
c
      hxcrse = hx*lratio
      hycrse = hy*lratio
      xl     = xleft - hxcrse
      xr     = xright + hxcrse
      yb     = ybot - hycrse
      yt     = ytop + hycrse
      maxci  = mxcip1 - 1
      maxcj  = mxcjp1 - 1
c
c first fill slopes in coarser grid for use in qsintfil
      mpar = lstart(lev-1)
 21   if (mpar .eq. 0) go to 22
	  mipar = node(5,mpar)-1+2*lwidth
	  mjpar = node(6,mpar)-1+2*lwidth
	  mparip1 = node(5,mpar)+1
	  mparjp1 = node(6,mpar)+1
	  node(3,mpar) = igetsp(8*mipar*mjpar)
	  locqx = node(3,mpar)
	  locqy = locqx+4*mipar*mjpar
	  call qsltemp(alloc(node(7,mpar)),alloc(locqx),alloc(locqy),
     .          mipar,mjpar,alloc(node(14,mpar)),node(17,mpar),
     .          lwidth,hxcrse,hycrse,
     .          rnode(1,mpar)-lwidth*hxcrse,
     .          rnode(2,mpar)-lwidth*hycrse,mparip1,mparjp1)
      mpar = node(10,mpar)
      go to 21
 22   continue
c
c  these 2 calls completely fill valc.  physbd (possibly) needed to fill
c  extra cell only.  note: intfil uses a used array. we get the
c  largest size needed (for the 2nd call), so don't have to get it twice.
c
      locuse = igetsp(maxip1*maxjp1)
      ifill  = 0
c qsintfil fills slopes qx,qy as well as valc. based on enlarged grid
      call qsintfil(xl,xr,yb,yt,lev-1,1,mxcip1,1,mxcjp1,nvar,valc,
     1            mxcip1,mxcjp1,time,intime,ifill,locuse,qx,qy)
      call physbd(xl,xr,yb,yt,lev-1,mxcip1,mxcjp1,nvar,valc,time
     1            ,hxcrse,hycrse)
 
      do 30 j = 2,maxcj
 
      do 30 ivar = 1,nvar
c
c slope calculations for entire row of coarse cells at a time
c
       do 10 i=2,maxci

        slp = valc(i+1,j,ivar) - valc(i,j,ivar)
        slm = valc(i,j,ivar)   - valc(i-1,j,ivar)
c
c really want to say "if (irc(i,j) .ne. regular)" but irc array was 
c set using several grids irr array, each with a diff. regular marker.
c
 	if (irc(i,j).eq.-1) then
  	    slopex = 0.d0
 	else if (irc(i,j).ne.lstnew) then
  	    slopex = qx(i+lwidth-1,j+lwidth-1,ivar)
 	else
 	    slopex = dmin1(dabs(slp),dabs(slm))*
     .               dsign(1.0d0,valc(i+1,j,ivar) - valc(i-1,j,ivar))
            if ( slm*slp .lt. 0.d0) then
               slopex = 0.d0
            endif
        endif
	dudx(i) = slopex
c
c same for y slopes
c
          slp = valc(i,j+1,ivar) - valc(i,j,ivar)
          slm = valc(i,j,ivar)   - valc(i,j-1,ivar)
c
c again, really want to say "if (irc(i,j) .ne. lstgrd)" but...
c
          if (irc(i,j).eq.-1) then
 	     slopey = 0.d0
          else  if (irc(i,j) .ne. lstnew) then
 	     slopey = qy(i+lwidth-1,j+lwidth-1,ivar)
	  else
 	     slopey = dmin1(dabs(slp),dabs(slm))*
     .                 dsign(1.0d0,valc(i,j+1,ivar) - valc(i,j-1,ivar))
             if ( slm*slp .lt. 0.d0) then
 	       slopey = 0.d0
             endif
	  endif
	  dudy(i) = slopey

  10   continue
c
c interp. from coarse cell to fine grid
c

      do 20 ico = 1,lratio

         do 20 jco = 1,lratio
         jfine = (j - 2)*lratio + 1 + jco

            do 20 i = 2,maxci
	      kc = irc(i,j)
	      if (kc .eq. lstnew) then
		xc =  xl + (i-.5d0)*hxcrse
		yc =  yb + (j-.5d0)*hycrse
	      else if (kc .ne. -1) then
		xc = xcirr(kc)
		yc = ycirr(kc)
	      endif
	      ifine   = (i - 2)*lratio + 1 + ico
	      kfine = irr(ifine+lwidth-1,jfine+lwidth-1)
	      if (kfine .eq. lstgrd .or. kfine .eq. -1) then
		xfine = xleft + (ifine-1.5d0)*hx
		yfine = ybot  + (jfine-1.5d0)*hy
	      else
		xfine = xcirr(kfine)
		yfine = ycirr(kfine)
	      endif
	      xoff = (xfine - xc)/hxcrse
	      yoff = (yfine - yc)/hycrse
	      val(ifine,jfine,ivar) = valc(i,j,ivar)
     1                            + xoff*dudx(i) + yoff*dudy(i)
 20         continue
c
 30   continue

c reclaim slope space before moving on

      mpar = lstart(lev-1)
 41   if (mpar .eq. 0) go to 42
	  mipar = node(5,mpar)-1+2*lwidth
	  mjpar = node(6,mpar)-1+2*lwidth
	  call reclam(node(3,mpar),8*mipar*mjpar)
      mpar = node(10,mpar)
      go to 41
 42   continue
c
c
c  overwrite interpolated values with fine grids values, if available.
c
      xl = xleft  - hx
      xr = xright + hx
      yb = ybot   - hy
      yt = ytop   + hy
      ifill = 1
      call intfil(xl,xr,yb,yt,lev,1,maxip1,1,maxjp1,
     1            nvar,val,maxip1,maxjp1,time,intime,ifill,locuse)
      call reclam(locuse, maxip1*maxjp1)
 
      return
      end
