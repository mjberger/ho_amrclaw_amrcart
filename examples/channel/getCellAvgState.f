c
c
c     =====================================================
       subroutine getCellAvgState(meqn,xlow,ylow,
     &                   dx,dy,qij,lstgrd,kirr,i,j)
c     =====================================================
c
c     # Set initial conditions for q.
c     # rotated channel problem
c     # i,j in 1, mitot coordinates
c
       use amr_module
       implicit double precision (a-h,o-z)
       include "quadrature.i"

       dimension qij(meqn),state(meqn)
       dimension fakeState(4)
       integer kirr
       data fakestate/0.5d0,0.d0,0.d0,2.5d0/

c
       if (kirr .eq. -1) then
             qij = fakeState
             return
       endif

       if (kirr .eq. lstgrd) then
           ! put in poly to evaluate the subtriangles
          call makep(poly(1,1,lstgrd),i,j,xlow,ylow,dx,dy)
       endif
       arr = ar(kirr)

        !LOOP over tris in poly to compute cell avg
       ivert = 3  ! first tri uses verts 1,2,3
       ntriquad = 3 ! num quadrature pts in each triangle 
       qij = 0.d0
       do while (poly(ivert,1,kirr) .ne. -11)
         call setVerts(poly(1,1,kirr),artri,x1,y1,x2,y2,x3,y3,ivert)
         do itq = 1, ntriquad
            xval = x1 * rtri(itq) + x2 * stri(itq)
     &           +  x3 * (1.d0-rtri(itq)-stri(itq))
            yval = y1 * rtri(itq) + y2 * stri(itq)
     &           +  y3 * (1.d0-rtri(itq)-stri(itq))
            call channelInit(xval,yval,state)
            qij(:) = qij(:) + (artri/arr)*wtri(itq) * state(:)
         end do
         ivert = ivert + 1  ! next tri is verts 1 3 4
       end do

       return
       end
