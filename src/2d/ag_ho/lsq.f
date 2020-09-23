c
c -------------------------------------------------------------------
c
       subroutine lsq(ix0,iy0,qMerged,qxmerge,qymerge,q,irr,mitot,mjtot,
     .                nvar,xlow,ylow,hx,hy,lstgrd,x0,y0,lwidth,istage)

       implicit double precision (a-h, o-z)
       include "cirr.i"
       
       dimension q(mitot,mjtot,nvar), irr(mitot,mjtot)
       dimension qxmerge(nvar), qymerge(nvar),qMerged(nvar)
       dimension rhsmax(nvar),rhsmin(nvar)

       integer asize
       parameter (nsize=25)  ! also in addMoreNeighbors
       dimension a(nsize,2),at(2,nsize),c(2,2)
       dimension b(nsize,nvar), d(nsize,nvar),rhs(nsize,nvar)
       dimension nlist(nsize,2)
       dimension qxexact(4), qyexact(4)

       logical quad,  nolimiter, debug/.false./, print/.true./

      common   /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      common /order2/ ssw, quad, nolimiter
       
c
c======================================================
c
c for now set to 0, see if code maintains pw const solution
c
!--       do ivar = 1, nvar
!--          qxmerge(ivar) = 0.d0
!--          qymerge(ivar) = 0.d0
!--       end do
!--       go to 99
c
c for linear test case put in exact gradient for rho
!--          qxmerge(1) = -.1d0
!--          qymerge(1) = 1.d0
!--          go to 99
c
c======================================================
c
c  do least square fit to get gradients using merged values and its neighbors
c  i,j is cut cell whose neighborhood made the merged cell
c  qmerge is merged soln in that cell
c  q is regular grid values in primitive variables
c
c  not sure who to include in stencil right now.
c
c  cant use qslopes routine since merged cell not a cartesian cell, and will choose
c  stencil differently
c
       nterms = 2
       nlist(1,1) = ix0
       nlist(1,2) = iy0
       nst  = 1
       nend = 1
c
c try using addneighbors as if merged cell were cut cell, may have to enlarge tho or
c not use neighbors or add weights but ...
c
       nhood = 1   ! size of neighboorhod to search, to start
       call addMoreNeighbors(irr,nlist,nst,nend,newend,mitot,mjtot,
     .                 lstgrd,quad,xlow,ylow,hx,hy,lwidth,nhood,istage)

       if (debug) then
         write(*,*)" for cell ",ix0,iy0," using ",newend," nbors "
       endif

       if (newend .ge. 5) go to 16  ! more stable with more cells
c
c      didnt get enough cells first time around, try again.
c

       
       nhood = 2
       nend  = 1  !reset and search again
       if (debug) then
          write(*,*)" second search bigger nhood for ",ix0,iy0
       endif
       call addMoreNeighbors(irr,nlist,nst,nend,newend,mitot,mjtot,
     .                lstgrd,quad,xlow,ylow,hx,hy,lwidth,nhood,istage)

       if (debug) then
          write(*,*)"   now have ",newend-1,
     .              " not counting cut cell itself"
       endif
c
       if (newend .ge. 5) go to 16  ! changed since trying to not use cut cell
        write(*,*)" NOT ENOUGH NEIGHBORS IN LSQ"
        write(*,*)" for cell ",ix0,iy0
        write(*,*)" only have ",newend-1
        do ivar = 1, nvar
         qxmerge(ivar)=0.d0
         qymerge(ivar)=0.d0
        end do
        go to 99
c       stop   !for debugging just stop and check instead of setting to 0
c
 16    irow = 0
c testing using all cut cells but with weights
c      :: only use cut cell itself (in position 1, seeding the list) if its area is large enough
c       if (ar(irr(ix0,iy0)) .gt. .5*ar(lstgrd)) then   !easy way to turn off using 1st cut cell
c          nst = 1
c       else
c          nst = 2
c       endif
       if (print) then
        write(*,919) ix0,iy0,newend
 919    format(" merged Cell gradient for ",2i5," used ",i5," nbors")
       endif

       nst = 1
       do 22 n = nst, newend 

          irow = irow + 1
          if (irow .gt. nsize) then
            write(*,*)" need larger arrays in lsq.f"
            stop
          endif
          ixn = nlist(n,1)
          iyn = nlist(n,2)
          kn =  irr(ixn,iyn)
          if (kn .ne. lstgrd) then
             xn = xcirr(kn)
             yn = ycirr(kn)
          else  !no solid cells on list
             xn = xlow + (ixn-.5d0)*hx
             yn = ylow + (iyn-.5d0)*hy
          endif
c
c   use volume weighting (sqrt(ar)) in this case so small cut cells dont count much
c
c          wtVolFrac = (ar(kn)/ar(lstgrd))**4
          wtVolFrac = 1.

          a(irow,1) = (xn - x0) * wtVolFrac
          a(irow,2) = (yn - y0) * wtVolFrac

          do m = 1, nvar
            b(irow,m) = (q(ixn,iyn,m) - qMerged(m))*wtVolFrac
          end do
 22       continue

c 
          do 30 it = 1, irow
          do 30 jt = 1, nterms
             at(jt,it) = a(it,jt)
 30       continue

          do 40  m = 1, nvar
             rhsmax(m) = b(1,m)
             rhsmin(m) = b(1,m)
             do 40 it = 1, irow
                rhs(it,m) = b(it,m)
                rhsmax(m) = dmax1(rhsmax(m),b(it,m))
                rhsmin(m) = dmin1(rhsmin(m),b(it,m))
 40          continue
c     
             do 50 it = 1, nterms
                do 50 jt = 1, nterms
                   c(it,jt) = 0.d0
                   do m = 1, nvar
                      d(it,m)  = 0.d0
                   end do
                   do 45 kt = 1, irow
                      c(it,jt) = c(it,jt) + at(it,kt)*a(kt,jt)
                      do m = 1, nvar
                      d(it,m) = d(it,m) + at(it,kt)*b(kt,m)
                      end do
 45                continue
 50             continue

c         here do linear fit

c     # solve C*w = d for least squares slopes. use cholesky
c     # put factors back in a
           a(1,1) = dsqrt(c(1,1))
           a(1,2) = c(1,2)/a(1,1)
           a(2,2) = dsqrt(c(2,2)-a(1,2)**2)
c     
           do 61 m = 1, nvar
c     
c     # at*a = c. solve at*b = d, aw = b.  reuse b.
              b(1,m) = d(1,m) / a(1,1)
              b(2,m) = (d(2,m) - a(1,2)*b(1,m)) / a(2,2)
              w2 =   b(2,m)/a(2,2)
              qymerge(m) =  w2
              qxmerge(m) = (b(1,m)-a(1,2)*w2)/a(1,1)
61         continue

c      write(*,900) (qxmerge(m),qymerge(m),m=1,nvar)
 900  format("      ",2e15.7)

c     turn off limiting 
        if (nolimiter)  go to 90
c
c LIMIT HERE

 90    if (print) then
         do ivar=2,nvar
           qxexact(ivar)=0.d0
           qyexact(ivar)=0.d0
         end do
         qxexact(1) = -.0d0
         qyexact(1) = 0.d0
         write(*,901) ix0,iy0,((qxmerge(ivar)-qxexact(ivar)),
     .          (qymerge(ivar)-qyexact(ivar)),ivar=1,nvar)
 901     format(" cell ",2i5, 4(2e30.20,/,16x))
       endif

c
 99    return
       end
