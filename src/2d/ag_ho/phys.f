c
c ----------------------------------------------------------
c
      logical function phys(x,y)
      implicit double precision (a-h,o-z)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
c
c  true = pt. (x,y) is on the physical boundary of the region or
c         exterior to the domain.
c  false (is an interior pt. of the domain)
c
c  set for the square on (0,xprob) X (0,yprob)
c
      phys = .false.
c
      if ( xprob-x .lt..00001 .or. yprob-y .lt. .00001 .or.
     1     x .lt. .00001 .or. y .lt.  .00001)
     2 phys = .true.
c
c
      return
      end
