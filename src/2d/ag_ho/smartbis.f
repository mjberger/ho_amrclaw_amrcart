c
c ---------------------------------------------------------
c
      subroutine smartbis(badpts,npts,cutoff,buffx,buffy,hx,hy,
     1                  newnum,nclust,lbase,corner,maxcl,
     2                  iscr,jscr,idim,jdim)
c
      implicit double precision (a-h,o-z)
      dimension     badpts(2,npts),corner(12,maxcl)
      dimension     iscr(idim), jscr(jdim)
      integer       nclust, newnum(maxcl)
      integer       horizontal, vertical
      parameter     (horizontal = 1)
      parameter     (vertical = 2)
      logical       gprint
      data          gprint/.true./
c
c smart bisect rectangles until cutoff reached for each.
c replaces old bisect/merge procedure
c
      if (gprint) write(6,100) nclust
 100  format(' starting smart bisection with ',i5,' clusters')
c
      icl         = 1
      ist         = 1
      iend        = newnum(icl)
c
 10   call moment(corner(1,icl),badpts(1,ist),newnum(icl),usenew,
     1             buffx,buffy,hx,hy)
      if (gprint) write(6,101) icl,newnum(icl),usenew
 101  format(' testing cluster ',i4,' with ',i5,' pts. use ',e12.4)
c
      if (usenew .lt. cutoff) go to 20
c
c  this cluster ok - on to next
c
      if (.not. gprint) go to 15
         write(6,102) icl,newnum(icl),usenew
 102     format(' accepting smart bisected cluster',i4,' with ',i5,
     1          ' pts. use = ',e10.3)
 15   icl   = icl + 1
      if (icl .gt. nclust) go to 200
      ist   = iend + 1
      iend  = ist + newnum(icl) - 1
      go to 10
c
c  smart bisect rectangle (and its cluster) in best location
c
 20   if (nclust .lt. maxcl) go to 25
          write(6,900) maxcl
 900      format('  too many clusters:  > ',i5)
          stop
 25   continue
c
c smart bisection computes signatures, finds best cut and splits there
c
      call signs(badpts,npts,iscr,jscr,idim,jdim,
     &           ist,iend,ilo,ihi,jlo,jhi,hx,hy)
      call findcut(icl,iscr,jscr,idim,jdim,index,iside,
     &             ilo,ihi,jlo,jhi)
      if (index .eq. 0) then
	 icl = icl + 1
	 if (icl .gt. nclust) go to 200
	 ist = iend + 1
	 iend = ist + newnum(icl) - 1
	 go to 10
      endif
c
      if (iside .eq. vertical) then
         fmid = (index-.5)*hy
	 idir = 2
      else
         fmid = (index-.5)*hx
	 idir = 1
      endif
c
      itop = ist - 1
      ibot = iend + 1
      i    = ist
 50   if (badpts(idir,i) .lt. fmid) go to 60
c
c  point in top half. let it stay, increment counter
c
        itop = itop + 1
        if (itop+1 .ge. ibot) go to 80
             i = i + 1
             go to 50
c
c  point in bottom half. switch with a bottom point that's not yet
c  checked, and increment bot. pointer
c
 60    ibot           = ibot - 1
       temp           = badpts(1,ibot)
       badpts(1,ibot) = badpts(1,i)
       badpts(1,i)    = temp
       temp           = badpts(2,ibot)
       badpts(2,ibot) = badpts(2,i)
       badpts(2,i)    = temp
       if (itop+1 .lt. ibot) go to 50
c
c done smartbisecting icl'th clusters. adjust counts, repeat bisect stage .
c
 80   newnum(icl) = itop - ist + 1
      icl         = icl + 1
c
c  bump down remaining clusters to make room for the new half of one.
c
      if (icl .gt. nclust) go to 120
      do 90 ico         = icl,nclust
      nmove             = nclust - ico + icl
 90   newnum(nmove + 1) = newnum(nmove)
 120  newnum(icl)       = iend - ibot + 1
      nclust            = nclust + 1
      iend              = itop
      icl               = icl - 1
      go to 10
c
c  done: there are nclust  clusters.
c
 200  continue
c
      return
      end
