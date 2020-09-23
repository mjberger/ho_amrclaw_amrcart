      subroutine transfer(qbest,mifine,mjfine,q,maxip1,maxjp1, 
     .                    mitot,mjtot,
     .                    irrfine,irr,lwidth,lstgrd,hx,hy)
  
      implicit double precision (a-h,o-z)
      include "cirr.i"

      dimension  pfine(10,2,irrsize),arfine(-1:irrsize),
     .		 xcfine(irrsize),ycfine(irrsize),ixfine(irrsize),
     .           iyfine(irrsize),nxtfine(irrsize)
      dimension  scratch(10,2)
      common /cirr2/  rlink(9,2,irrsize)

      dimension irrfine(mifine,mjfine), irr(mitot,mjtot)
      dimension qbest(mifine,mjfine,4), q(maxip1,maxjp1,4)
      logical border
      border(i,j) = (i.le.lwidth .or. j.le.lwidth .or.
     .               i .gt. mifine-lwidth .or. j .gt. mjfine-lwidth)


      read(43) irrfine
      read(43) pfine,arfine
      read(43) xcfine,ycfine,ixfine,iyfine,nxtfine
      read(43) qbest

c
c transfer info from finer to coarser grid
c if a coarse cell is regular, so are all fine cells within
c only fix coarse cut cells, but use both cut and full fine cells
c to do it
c
      lratio = 4
      rl2    = dfloat(lratio * lratio)
      lstfine = 1
      hxfine = hx / dfloat(lratio)
      hyfine = hy / dfloat(lratio)

      do 10 i = lwidth+1, mitot-lwidth
      do 10 j = lwidth+1, mjtot-lwidth
c     ## handle lwidth border somewhere ##

      if ((irr(i,j) .eq. lstgrd) .or. (irr(i,j) .eq. -1)) go to 10
      k = irr(i,j)

      ifine = (i-lwidth)*lratio + 1
      jfine = (j-lwidth)*lratio + 1

      asum = 0.d0
      xsum = 0.d0
      ysum = 0.d0
      qsum1 = 0.d0
      qsum2 = 0.d0
      qsum3 = 0.d0
      qsum4 = 0.d0
      shiftxx = 0.d0
      shiftxy = 0.d0
      shiftyy = 0.d0

      do 20 ico = 1, lratio
      do 20 jco = 1, lratio

	kfine = irrfine(ifine+ico-1,jfine+jco-1)
	if (kfine .eq. lstfine) then
	   area = hxfine*hyfine
	   xcf =  (dfloat(ifine+ico-1-lwidth)-.5)*hxfine
	   ycf =  (dfloat(jfine+jco-1-lwidth)-.5)*hyfine
	else
	   area = arfine(kfine)
	   if (kfine .ne. -1) then
	     xcf = xcfine(kfine)
	     ycf = ycfine(kfine)
	   endif
	endif

	asum = asum + area
	xsum = xsum + area*xcf
	ysum = ysum + area*ycf
	shiftxx = shiftxx +  pfine( 8,1,kfine)*area
	shiftxy = shiftxy +  pfine( 9,1,kfine)*area
	shiftyy = shiftyy +  pfine(10,1,kfine)*area

c       # qbest already multiplied by arfine
	qsum1 = qsum1 + qbest(ifine+ico-1,jfine+jco-1,1)
	qsum2 = qsum2 + qbest(ifine+ico-1,jfine+jco-1,2)
	qsum3 = qsum3 + qbest(ifine+ico-1,jfine+jco-1,3)
	qsum4 = qsum4 + qbest(ifine+ico-1,jfine+jco-1,4)

 20   continue

      ar(k)    = asum
      xcirr(k)    = xsum / asum
      ycirr(k)    = ysum / asum

c     q(i-lwidth+1,j-lwidth+1,1) = qsum1/asum
c     q(i-lwidth+1,j-lwidth+1,2) = qsum2/asum
c     q(i-lwidth+1,j-lwidth+1,3) = qsum3/asum
c     q(i-lwidth+1,j-lwidth+1,4) = qsum4/asum
c     poly( 8,1,k) = shiftxx/asum
c     poly( 9,1,k) = shiftxy/asum
c     poly(10,1,k) = shiftyy/asum

 10   continue

c
c hard part next - transfer the boundary segments
c
      kfine = lstfine
 21   kfine = iabs(nxtfine(kfine))
      if (kfine .eq. 0) go to 99

      ixuse = ixfine(kfine)
      iyuse = iyfine(kfine)
 25   if (border(ixuse,iyuse))  then
	go to 21
      else
	ixc = (ixuse-lwidth-1)/lratio + lwidth + 1
	iyc = (iyuse-lwidth-1)/lratio + lwidth + 1
      endif
      kco = irr(ixc,iyc)
      ilink = 1

c     copy the link. copy link bckwards, (as in hs) and
c     reverse everything later
 30   continue
c   find ivert
      do 31 ivert = 1, 10
      if (pfine(ivert+1,1,kfine) .eq. -11) go to 32
 31   continue

 32   continue
c     rlink(ilink,1,kco) = hsfine(1,1,kfine)
c     rlink(ilink,2,kco) = hsfine(1,2,kfine)
      rlink(ilink,1,kco) = pfine(ivert,1,kfine)
      rlink(ilink,2,kco) = pfine(ivert,2,kfine)
c     hxsave = hsfine(2,1,kfine)
c     hysave = hsfine(2,2,kfine)
      hxsave = pfine(ivert-1,1,kfine)
      hysave = pfine(ivert-1,2,kfine)
      ilink = ilink + 1
c     get next one
      kfine = iabs(nxtfine(kfine))
      if (kfine .eq. 0) go to 99
      ixuse = ixfine(kfine)
      iyuse = iyfine(kfine)
      if (border(ixuse,iyuse)) go to 35
      ixcn = (ixuse-lwidth-1)/lratio + lwidth + 1
      iycn = (iyuse-lwidth-1)/lratio + lwidth + 1
      if (ixcn.eq.ixc .and. iycn.eq.iyc) go to 30
      
c     finish off previous link and get ready for next one.
 35   rlink(ilink,1,kco) = hxsave
      rlink(ilink,2,kco) = hysave
      rlink(ilink+1,1,kco) = -11
      if (kfine .ne. 0) go to 25
c
c     reverse
 99   continue
      do 40 i = lwidth+1, mitot-lwidth
      do 40 j = lwidth+1, mjtot-lwidth

      kco = irr(i,j)
      if (kco .eq. lstgrd .or. kco .eq. -1) go to 40
        do 50 ilink = 1, 10
	  scratch(ilink,1) = rlink(ilink,1,kco)
	  scratch(ilink,2) = rlink(ilink,2,kco)
	  if (scratch(ilink,1) .eq. -11) then
	     isave = ilink - 1
	     go to 60
	  endif
 50    continue
 60    do 65 ilink = 1, isave
	rlink(isave+1-ilink,1,kco) = scratch(ilink,1)
	rlink(isave+1-ilink,2,kco) = scratch(ilink,2)
 65    continue
       rlink(isave+1,1,kco) = -11
c	      write(6,*) " links for cell kco = ",kco
c	      do 84 ilink = 1, isave+1
c	      write(6,929) rlink(ilink,1,kco),rlink(ilink,2,kco)
c 929          format(2e15.7)
c 84           continue
c	      write(6,*)

 40    continue

c fix border cells - specific for circular vortex problem.
c i=1,4 two segments, and j=1,4 2 segments. nothing at the other end.
c
       do 70 i = 1, lwidth
	  do 70 j = 1, mjtot
	     if (irr(i,j) .eq. lstgrd .or. irr(i,j).eq.-1)go to 70
c            copy reflection - interior cell  info
             isrc     = 9-i
	     k        = irr(i,j)
	     ksrc     = irr(isrc,j)
	     ar(k)    = ar(ksrc)
	     xcirr(k) = -xcirr(ksrc)
	     ycirr(k) =  ycirr(ksrc)
c            reverse the list then copy 
	      do 51 ilink = 1, 10
		scratch(ilink,1) = rlink(ilink,1,ksrc)
		scratch(ilink,2) = rlink(ilink,2,ksrc)
		if (scratch(ilink,1) .eq. -11) then
		   isave = ilink - 1
		   go to 61
		endif
 51             continue
 61           do 52 ilink = 1, isave
		rlink(isave+1-ilink,1,k) = -scratch(ilink,1)
		rlink(isave+1-ilink,2,k) =  scratch(ilink,2)
 52           continue
	      rlink(isave+1,1,k) = -11
 70    continue

       do 80 j = 1, lwidth
	  do 80 i = 1, mitot
	     if (irr(i,j) .eq. lstgrd .or. irr(i,j).eq.-1)go to 80
c            copy reflection - interior cell  info
             jsrc     = 9-j
	     k        = irr(i,j)
	     ksrc     = irr(i,jsrc)
	     ar(k)    = ar(ksrc)
	     xcirr(k) =  xcirr(ksrc)
	     ycirr(k) = -ycirr(ksrc)
c            reverse the list then copy 
	      do 81 ilink = 1, 10
		scratch(ilink,1) = rlink(ilink,1,ksrc)
		scratch(ilink,2) = rlink(ilink,2,ksrc)
		if (scratch(ilink,1) .eq. -11) then
		   isave = ilink - 1
		   go to 91
		endif
 81             continue
 91           do 82 ilink = 1, isave
		rlink(isave+1-ilink,1,k) =  scratch(ilink,1)
		rlink(isave+1-ilink,2,k) = -scratch(ilink,2)
 82           continue
	      rlink(isave+1,1,k) = -11
 80    continue

      continue
      return
      end
