c
c ---------------------------------------------------------
c
      subroutine bisect(badpts,npts,cutoff,buffx,buffy,hx,hy,
     1                  newnum,nclust,lbase,corner,maxcl)
c
      implicit double precision (a-h,o-z)
      dimension          badpts(2,npts),corner(12,maxcl)
      integer       nclust, newnum(maxcl)
      logical       gprint
      data          gprint/.false./
c
c bisect rectangles until cutoff reached for each.
c this will be followed by a merge step
c
      if (gprint) write(6,100) nclust
 100  format(' starting bisection with ',i5,' clusters')
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
 102     format(' accepting bisected cluster',i4,' with ',i5,
     1          ' pts. use = ',e10.3)
 15   icl   = icl + 1
      if (icl .gt. nclust) go to 200
      ist   = iend + 1
      iend  = ist + newnum(icl) - 1
      go to 10
c
c  bisect rectangle (and its cluster) in long direction
c
 20   if (nclust .lt. maxcl) go to 25
          write(6,900) maxcl
 900      format('  too many clusters:  > ',i5)
          stop
 25   rlen1  = corner(7,icl) - corner(1,icl)
      rlen2  = corner(4,icl) - corner(2,icl)
      if (rlen1 .gt. rlen2) go to 30
         evx   = 0.0
         evy   = 1.0
         rhalf = (corner(4,icl) + corner(2,icl) ) / 2.
         go to 40
 30   evx   = 1.0
      evy   = 0.0
      rhalf = (corner(7,icl) + corner(1,icl) ) / 2.
c
c  do in place bisect of the icl'th cluster
c
  40  itop = ist - 1
      ibot = iend + 1
      i    = ist
 50   prod = badpts(1,i) * evx + badpts(2,i) * evy
      if (prod .lt. rhalf) go to 60
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
c done bisecting icl'th clusters. adjust counts, repeat bisect stage .
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
c  done with bisect stage. there are nclust 'initial' clusters.
c  next step is to repeatedly merge neighboring clustering until
c  no more successful merges can be done.
c
 200  continue
c
      return
      end
