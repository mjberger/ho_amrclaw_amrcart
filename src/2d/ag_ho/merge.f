c
c ---------------------------------------------------------
c
      subroutine merge(work,badpts,npts,cutoff,buffx,buffy,
     1                 hx,hy,numptc,nclust,corner,lratio,maxcl)
c
      implicit double precision (a-h,o-z)
      dimension          badpts(2,npts), work(2,npts), corner(12,maxcl)
      integer       numptc(maxcl)
      logical       mprint, mdebug
      data          mprint/.true./,mdebug/.false./
c
c process  into rectangles
c
c data structures :
c   work(2,npts)    : used for work backup storage only
c   corner(12,maxcl): has the 12 fields for the rectangles of
c                     the upto maxcl clusters.
c   there are nclust 'initial' clusters. (nclust>1 or wouldn't call)
c   repeatedly merge clusters until no more successful merges left.
c
c  for now trying merging with all other clusters. try for one
c  succession of merges, rather than bottom up, two at a time.
c
c  will make one pass through trying all possible merges.
c  will also bump down points in linear list style.
c  all pointers into the badpts array point to one less than
c  then starting location, to avoid subtracting 1 each index op.
c
c  try merging cluster  ist  with all the rest.
c
      if (mprint) write(6,100) nclust
 100  format(i5,' clusters before merge')
      ist   = 1
      nist  = 0
 10   itry   =  1
      itrypt = 0
      if (itry .eq. ist) then
         itry   = itry + 1
         itrypt = itrypt + numptc(ist)
      endif
      if (itry .gt. nclust) go to 110
 20   c1x    = dmin1(corner(1,ist),corner(1,itry))
      c1y    = dmin1(corner(2,ist),corner(2,itry))
      c2y    = dmax1(corner(4,ist),corner(4,itry))
      c4x    = dmax1(corner(7,ist),corner(7,itry))
      iside1 = (hx + c4x - c1x - 2.0*buffx + .00001) / hx
      iside2 = (hy + c2y - c1y - 2.0*buffy + .00001) / hy
      gpall  = iside1 * iside2
      usage  = dfloat(numptc(ist) + numptc(itry)) / gpall
c
      side1  = lratio*(corner(7,ist)-corner(1,ist))/hx  + 1
      side2  = lratio*(corner(4,ist)-corner(2,ist))/hy + 1
      side3  = lratio*(corner(7,itry)-corner(1,itry))/hx + 1
      side4  = lratio*(corner(4,itry)-corner(2,itry))/hy + 1
      side5  = lratio*(c4x-c1x)/hx + 1
      side6  = lratio*(c2y-c1y)/hy + 1
      cost1 = side1*side2 + side1+side2
      cost2 = side3*side4 + side3+side4
      costc = side5*side6 + side5+side6
c
      if (mdebug) write(6,101) ist,itry,numptc(ist),numptc(itry),usage,
     1                         cost1,cost2,costc
 101  format(' trying cluster ',i3,' with cluster ',i3,' :',
     1      i4,' pts + ',i4,' pts',/,
     2       '  use= ',e11.3,' cost1,2,c ',3f10.3)
      if (mdebug) write(6,102) (corner(i,ist),i=1,8),
     1                         (corner(i,itry),i=1,8)
 102  format(' corners: ',4('(',2f7.4,')'), /, 10x, 4('(',2f7.4,')'))
      if (( usage .gt. cutoff) .or. (costc .lt. cost1+cost2))
     1   go to 30
c
c  merge rejected. try next
c
       itrypt = itrypt + numptc(itry)
       itry = itry + 1
       if (itry .eq. ist)  then
          itry   = itry + 1
          itrypt = itrypt + numptc(ist)
        endif
       if (itry .le. nclust) go to 20
       go to 110
c
c  merge accepted. merge ist and itry clusters.
c
c  save itry cluster in work. shift pts. down. move work to be
c  near ist pts. need to save pts. for nest-checking later.
c  nmove = number of points to be moved
c
c if points adjacent (most likely case), don't bother
c
 30   if ((itry .eq. ist+1) .or. (itry .eq. ist-1)) go to 72
      do 50 i   = 1, numptc(itry)
      work(1,i) = badpts(1,itrypt+i)
 50   work(2,i) = badpts(2,itrypt+i)
c
      if (itry .gt. ist) then
         nmove   = itrypt - numptc(ist) - nist
         nblank  = itrypt + numptc(itry) +1
         do 60 i = 1, nmove
         badpts(1,nblank-i) = badpts(1,itrypt+1-i)
 60      badpts(2,nblank-i) = badpts(2,itrypt+1-i)
      else
         nmove   = nist - itrypt - numptc(itry)
         nbump   = itrypt + numptc(itry)
         do 62 i = 1, nmove
         badpts(1,itrypt+i) = badpts(1,nbump+i)
 62      badpts(2,itrypt+i) = badpts(2,nbump+i)
      endif
c
      if (itry .gt. ist) then
         nafter = nist + numptc(ist)
      else
         nafter = nist - numptc(itry)
      endif
      do 70 i            = 1, numptc(itry)
      badpts(1,nafter+i) = work(1,i)
 70   badpts(2,nafter+i) = work(2,i)
c
 72   numptc(ist) = numptc(ist)+numptc(itry)
      if (itry .gt. ist) itrypt = itrypt + numptc(itry)
      if (itry .lt. ist) then
         ist  = ist - 1
         nist = nist - numptc(itry)
       endif
      nclust      = nclust - 1
      do 80  i    = itry, nclust
      do 75  ipos    = 1, 12
 75   corner(ipos,i) = corner(ipos,i+1)
 80   numptc(i)   = numptc(i+1)
c
c  fix rectangle to reflect merger
c
      corner(1,ist)  = c1x
      corner(2,ist)  = c1y
      corner(3,ist)  = c1x
      corner(4,ist)  = c2y
      corner(5,ist)  = c4x
      corner(6,ist)  = c2y
      corner(7,ist)  = c4x
      corner(8,ist)  = c1y
c
      if (mdebug) write(6,103) numptc(ist)
 103  format(' merge accepted, new # points ',i5)
c
      if (itry .eq. ist) then
         itry   = itry + 1
         itrypt = itrypt + numptc(ist)
      endif
      if (itry .le. nclust) go to 20
c
c  no more clusters can merge with cluster ist. try next cluster.
c
 110  nist = nist + numptc(ist)
      ist  = ist + 1
      if (ist .le. nclust) go to 10
c
c  no more clusters left to try merging. accept the last cluster by
c  itself.
c
      if (mprint) write(6,104) nclust
 104  format(i4,' clusters after merge ')
      return
      end
