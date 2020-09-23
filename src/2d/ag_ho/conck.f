c
c -----------------------------------------------------------
c
      subroutine conck(level, nvar, time)
c
      implicit double precision (a-h,o-z)

      logical    rflag,vtime,steady
      logical    graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder, mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "cirr.i"
      include "calloc.i"

      dimension ialloc(2*allocsize)
      equivalence(alloc,ialloc)

      iadd(i,j,ivar)  = loc + i - 1 + maxip1*((ivar-1)*maxjp1+j-1)
      iirreg(i,j)     = 2*(locirr-1)+1 + i - 1 + mitot*(j-1)
c
c ******************************************************************
c conck - conservation check  for specified level
c         mostly a debugging tool
c         this assumes grids dont overlap
c ******************************************************************
c
c
c  grid loop for given level
c
      hx      = hxposs(level)
      hy      = hyposs(level)
      dt      = possk(level)
      totmass = 0.d0

      mptr = lstart(level)
 20   if (mptr .eq. 0) go to 85
         loc    = node(7,mptr)
         maxi   = node(5,mptr) 
         maxj   = node(6,mptr)
         maxip1 = maxi + 1
         maxjp1 = maxj + 1
         mitot  = maxi-1+2*lwidth
         mjtot  = maxj-1+2*lwidth
         locirr = node(14,mptr)
         lstgrd = node(17,mptr)
         ar(lstgrd) = hx*hy
         ar(-1)     = 0.d0
c
         do 50 j  = 2, maxip1-1
         do 50 i  = 2, maxjp1-1
            iadj = i + lwidth - 1
            jadj = j + lwidth - 1
            k = ialloc(iirreg(iadj,jadj))
            totmass = totmass + ar(k)*alloc(iadd(i,j,1)) 
 50      continue
c
       mptr = node(10,mptr)
       go to 20
c
 85    write(6,777) time, totmass
       write(*,777) time, totmass
 777   format('At time t = ',e15.7,',   total mass = ',e30.20)
c
 99   return
      end
