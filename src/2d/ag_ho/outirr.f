c
c -------------------------------------------------------------
c
      subroutine outirr(irr,mitot,mjtot,lhptr)
c
      implicit double precision (a-h,o-z)
      dimension irr(mitot,mjtot)
      include "cirr.i"
      common  /cloops/  xloops(10),yloops(10),nloops
c
      write(3,101) ((irr(i,j),i=1,mitot),j=1,mjtot)
 101  format(20i4)
      write(3,104) nloops,(xloops(i),yloops(i),i=1,nloops)
 104  format(i2,5(2e13.5))
      write(3,102) lhptr
 102  format(i6,'  start of irregular list info ')
      lptr = lhptr
 10   write(3,103) lptr,ix(lptr),iy(lptr),nxtirr(lptr),
     1             xcirr(lptr),ycirr(lptr),
     1             ((poly(i,j,lptr),i=1,10),j=1,2)
 103  format(4i6,2e15.7,4(/,5e15.8))
      if (nxtirr(lptr) .eq. 0) go to 99
	  lptr = iabs(nxtirr(lptr))
	  go to 10
c
 99   return
      end
