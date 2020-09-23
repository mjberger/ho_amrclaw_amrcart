c
c ------------------------------------------------------------------
c
      subroutine irrlen(f,g,firreg,irr,mitot,mjtot,dt,dx,dy,lstgrd,
     &                  mptr,meqn,ffluxlen,gfluxlen,msize)

      implicit double precision (a-h,o-z)
      include "cirr.i"

      dimension f(mitot,mjtot,meqn),g(mitot,mjtot,meqn),irr(mitot,mjtot)
      dimension firreg(-1:irrsize,meqn)
      dimension ffluxlen(msize,msize)
      dimension gfluxlen(msize,msize)
c
c     # irregular side modifications to fluxes for length only
c     # this used to be part of irreg1/2.
c
      logical done(irrsize),edgein(4)
c
c     # initialize fluxes at boundary to zero and set done array to false.
c     # done(k) is set to true when k'th boundary cell is handled and is
c     # used to insure that interfaces between boundary cells
c     # are only handled once.
c
      do 20 ik = 1, irrsize
 20      done(ik) = .false.
      do 30 m = 1, meqn
      do 30 ik = -1, irrsize
 30       firreg(ik,m) = 0.d0
c
c     # modify fluxes on sides of irregular cells to reflect the fact that
c     # the length of each side may be less than the full dx or dy.
c     # (Recall that the fluxes f and g are the actual fluxes multiplied by
c     # the length of the side).
c
       k = lstgrd
       do 65 ik=1,irrsize
          done(k) = .true.
          k = iabs(nxtirr(k))
          if (k.eq.0) go to 66
c     # determine which edge of cell the first side is on (1=left,
c     # 2=top, 3=right, 4=bottom):
          if (poly(1,1,k).eq.poly(2,1,k)) then
c         # vertical edge (1 or 3)
             if (poly(1,2,k).gt.poly(2,2,k)) then
                iedge = 3
             else
                iedge = 1
             endif
          else
c     # horizontal edge (2 or 4)
             if (poly(1,1,k).gt.poly(2,1,k)) then
                iedge = 4
             else
                iedge = 2
             endif
          endif
c     
c     # keep track of which cell edges are present in this polygon:
          do 35 ied=1,4
 35          edgein(ied) = .false.
c     
c     # march around the sides of this cell:
          do 60 kside=1,6
             if (poly(kside+2,1,k).eq.-11) go to 62
             edgein(iedge) = .true.
c     # determine coordinates of adjacent cell:
             ixadj = ix(k)
             iyadj = iy(k)
             go to (41,42,43,44) iedge
 41          ixadj = ixadj-1
             go to 45
 42          iyadj = iyadj+1
             go to 45
 43          ixadj = ixadj+1
             go to 45
 44          iyadj = iyadj-1
             go to 45
 45          continue
c     
c     # if the adjacent cell has already been dealt with (or if it
c     # is regular) then we need not modify the flux for this side:
             if (ixadj.eq.0 .or. ixadj.gt.mitot .or.
     &           iyadj.eq.0 .or. iyadj.gt.mjtot) go to 58
             kirr = irr(ixadj,iyadj) ! protect against zero length edge
             if (kirr .ne. -1) then ! zero out flux 
                    if (done(kirr)) go to 58
             endif
c     
c     # compute the length of this side and modify the appropriate
c     # flux:
             hside = dsqrt((poly(kside,1,k)-poly(kside+1,1,k))**2 +
     &               (poly(kside,2,k)-poly(kside+1,2,k))**2)
             do 55 m=1,meqn
                go to (51,52,53,54) iedge
 51             f(ix(k),iy(k),m) = hside/dy * f(ix(k),iy(k),m)
                ffluxlen(ix(k),iy(k)) =  hside
                go to 55
 52             g(ixadj,iyadj,m) = hside/dx * g(ixadj,iyadj,m)
                gfluxlen(ixadj,iyadj) = hside
                go to 55
 53             f(ixadj,iyadj,m) = hside/dy * f(ixadj,iyadj,m)
                ffluxlen(ixadj,iyadj) = hside
                go to 55
 54             g(ix(k),iy(k),m) = hside/dx * g(ix(k),iy(k),m)
                gfluxlen(ix(k),iy(k)) = hside                
                go to 55
 55          continue
c     
 58          continue
             iedge = iedge+1
             if (iedge.gt.4) iedge = iedge-4
 60       continue
 62    continue
c     
c     # set fluxes to zero for cell edges that are not in this poly:
c     
          do 63 iedge=1,4
              if (.not. edgein(iedge)) then
                do 64 m=1,meqn
                   go to (601,602,603,604) iedge
 601               f(ix(k),iy(k),m) = 0.d0
                   ffluxlen(ix(k),iy(k)) = 0.d0
                   go to 64
 602               if (iy(k).eq.mjtot) go to 64
                   g(ix(k),iy(k)+1,m) = 0.d0
                   gfluxlen(ix(k),iy(k)+1) = 0.d0
                   go to 64
 603               if (ix(k).eq.mitot) go to 64
                   f(ix(k)+1,iy(k),m) = 0.d0
                   ffluxlen(ix(k)+1,iy(k)) = 0.d0
                   go to 64
 604               g(ix(k),iy(k),m) = 0.d0
                   gfluxlen(ix(k),iy(k)) = 0.d0
                   go to 64
 64             continue
              endif
 63       continue
 65     continue
c     
c     
 66    continue
c     
       return
       end
