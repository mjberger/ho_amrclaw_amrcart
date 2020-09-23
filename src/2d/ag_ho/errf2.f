c
c --------------------------------------------------------------
c
      subroutine errf2(rctold,nvar,maxi,maxj,maxip1,maxjp1,
     1                 numflg,mptr,rbuff,nbuff,iflow,irr,mitot,mjtot)
c
c Finish by:
c 1. flagging points that need to go on buffer list. count
c    number of flagged points.
c 2. unflag (if necesssary) flagged points exterior to the domain
c    that were accidentally flagged due to the buffer flagging.
c 3. repeat entire code from errf1, saving points to be flagged that lie
c    outside grid - might be on neighboring grid.
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
      include "cirr.i"
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      dimension rctold(maxip1,maxjp1,nvar),rbuff(3,nbuff)
      dimension irr(mitot,mjtot)
      logical   eprint,  edebug
      logical   ingrid
      data      eprint /.true./, edebug/.false./
      ingrid(iflag,jflag) = ((iflag .ge. 2) .and. (iflag .le. maxi) 
     .               .and. (jflag .ge. 2) .and. (jflag .le. maxj))
c
c flag buffer zone (prebuffering), with a different no. so doesn't
c get recursively out of hand. colate will account for this.
c take special care if flag buffer point outside the grid, but there
c is another grid which could be flagged.
c
c don't flag points inside body. this may cause a problem with finer
c grids being level-nested. projec routine should fix this, except on
c a restart run, until all levels are regridded.
c
       level  = node(4, mptr)
c no buffer zone when flagging for geometry
c      if (level .eq. 1) go to 86
c
       do 85 j = 2, maxj
       do 85 i = 2, maxi
         if (rctold(i,j,1) .ne. 2.0) go to 85
	 if (irr(i+lwidth-1,j+lwidth-1) .eq. -1) go to 85
         ist  = i-bzonex
         iend = i+bzonex
         jst  = j-bzoney
         jend = j+bzoney
         do 84 iflag = ist, iend
         do 84 jflag = jst, jend
         if (ingrid(iflag,jflag)) then
	     if ((irr(iflag+lwidth-1,jflag+lwidth-1) .eq. -1) .and.
     .            (level+1 .eq. mxnest)) then
                 rctold(iflag,jflag,1) = 0.d0
	     else
		 if (rctold(iflag,jflag,1) .eq. 0.d0) 
     .               rctold(iflag,jflag,1) = 1.d0
	     endif
	 go to 84
	 endif
	 if (irr(iflag+lwidth-1,jflag+lwidth-1) .eq. -1) go to 84
               iflow = iflow + 1
               if (iflow .gt. nbuff) then
                  write(6,900) nbuff
 900              format(' out of flagging buffer zone space in rbuff:',/,
     1                   ' dimensioned at nbuff = ',i7)
                  stop
               endif
            rbuff(1,iflow)=rnode(1,mptr)+rnode( 9,mptr)
     1                          *(dfloat(iflag)-1.5)
            rbuff(2,iflow)=rnode(2,mptr)+rnode(10,mptr)
     1                          *(dfloat(jflag)-1.5)
            rbuff(3,iflow)=dfloat(mptr)
 84      continue
 85      continue
c
 86   continue
      do 90 i = 2, maxi
      do 90 j = 2, maxj
      if (rctold(i,j,1) .ne. 0) numflg = numflg + 1
 90   continue
c
      if (eprint) then
         write(6,105) mptr,level,numflg
 105     format(' grid',i4,' level',i4,' has ',i4,
     1              ' pts flagged in total ')
	 if (numflg .gt. 0) then
            do 95 jj = 2, maxj
            j        = maxj + 2 - jj
            write(6,106) (nint(rctold(i,j,1)),i=2,maxi)
106         format(1h ,80i1)
 95         continue
	 endif
      endif
c
 99   return
      end
