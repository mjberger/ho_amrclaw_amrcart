c
c ---------------------------------------------------------------------
c
       subroutine SRD_cellMerge(q,nvar,irr,mitot,mjtot,lstgrd,
     .                      dx,dy,lwidth,xlow,ylow,istage,
     .                      numHoods,mptr)

       use amr_module
       implicit double precision (a-h, o-z)
       include "cuserdt.i"
       include "quadrature.i"

       dimension q(nvar,mitot,mjtot),  irr(mitot,mjtot)
       dimension gradmx(nvar,irrsize), gradmy(nvar,irrsize)
       dimension gradmxx(nvar,irrsize), gradmyy(nvar,irrsize)
       dimension gradmxy(nvar,irrsize)

       dimension valnew(nvar,mitot,mjtot)
       dimension qMerge(nvar,mitot,mjtot), numHoods(mitot,mjtot)
       dimension mioff(mitot,mjtot),mjoff(mitot,mjtot)
       dimension fakeStateCons(nvar), qm(nvar), rhs(5,nvar)
       dimension a(5,5),b(5),db(nvar)
       dimension nborList(35,2)
       character c2

       logical OUT_OF_RANGE
       logical quad, nolimiter
       common /order2/ ssw, quad, nolimiter
       integer omp_get_max_threads
       integer maxthreads/1/


c  xlow,ylow   refers to the corner of grid (includes ghost cells)
c  xlower,ylower refers to the domain

c NEW WAY - do one stage at a time then copy in ghost cells at intermediate
c stages for neighboring grids.  
c never use cells from exterior to the domain for SRD


      OUT_OF_RANGE(i,j) = (i .lt. 1 .or. i .gt. mitot .or.
     &                     j .lt. 1 .or. j .gt. mjtot)

c :::::::::::::;:::
c
c   try cell merging and reconstruction to stabilize updates in small cut cells
c   calculate all updates using provisional values, then actually do the update
c   (as in, jacobi rather than gauss seidel)
c   
c   the setup for this is called from setirr, and is in routine makeMergeHood
c ::::::::::::::::

c

       call countCellType(irr,mitot,mjtot,lwidth,numSolid,numCut,
     &                    numFull,lstgrd)
       nx = mitot - 2*lwidth 
       ny = mjtot - 2*lwidth 
       if (numSolid .eq. nx*ny) return

c      some initializations
       ar(lstgrd) = dx*dy   ! area of regular grid cell 
       areaMin = areaFrac*dx*dy
       qMerge   = 0.d0
         

       ! put in 'normal' values to prevent errors e.g. in converting to prim vars
       ! these fake vals are in conserved variables
       fakeStateCons(1) = 1.4d0
       fakeStateCons(2) = 0.d0
       fakeStateCons(3) = 0.d0
       fakeStateCons(4) = 2.5d0 !(p/gm1, corresponding to fakestate in method)

c   merging nhoods make in makeMergeHood, called from setirr

      if (igradChoice .eq. 3 .or. igradChoice .eq. 4) then
        call merge_shifts(irr,mitot,mjtot,lwidth,dx,dy,xlow,ylow,lstgrd,
     .                  numHoods)
      endif


c       form qMerge vals 
        do 11 j = 1, mjtot
        do 10 i = 1, mitot
            k = irr(i,j)
            if (k .eq. -1) go to 10 ! no solid cells
            ! with multistep RK, or if not then adjust do loop indices
            if (OUT_OF_RANGE(i,j)) then 
              !qMerge(:,i,j) = rinfinity ! to make sure we dont use it
              qMerge(:,i,j) = q(:,i,j) ! to make sure we dont use it
              go to 10 
            endif
            if (k.eq.lstgrd) then ! full cell is its own nhood
              qMerge(:,i,j) = q(:,i,j)
              go to 10
            endif
            if (ar(k) .gt. areaMin) then
              qMerge(:,i,j) = q(:,i,j)  ! cut cell is large enough
              go to 10
            endif

            ! my ncount doesn't include cell itself so initialize qMerge to it, not 0
            qMerge(:,i,j) = ar(k)*q(:,i,j)/numHoods(i,j)
            do ic = 1, ncount(k)
               icurr = iidx(ic,k)
               jcurr = jidx(ic,k)
               koff = irr(icurr,jcurr)
               qMerge(:,i,j) = qMerge(:,i,j) + ar(koff)*
     &                         q(:,icurr,jcurr)/numHoods(icurr,jcurr)
            end do
            qMerge(:,i,j) = qMerge(:,i,j) / volMerge(k)

 10     continue
 11     continue

        ! gradient of merging neighborhoods, initialized to 0. set using neighboring tiles
        ! local variables so ok to zero out
        gradmx = 0.d0 
        gradmy = 0.d0
        gradmxx = 0.d0
        gradmxy = 0.d0
        gradmyy = 0.d0

        ! compute stable neighborhood (mioff,mjoff) for gradients on merging tiles
        call makeMergeGradHood(irr,lwidth,mitot,mjtot,lstgrd,dx,dy,
     &                      xlow,ylow,mptr,mioff,mjoff)
        if (igradChoice .eq. 3) then
           call qmslopes(irr,mitot,mjtot,lwidth,dx,dy,xlow,ylow,lstgrd,
     &                numHoods,mioff,mjoff,qMerge,nvar,
     &                gradmx,gradmy,gradmxx,gradmxy,gradmyy)
        else if (igradChoice .eq. 1 .or. igradChoice .eq. 2) then
          call mnslopes(irr,mitot,mjtot,lwidth,dx,dy,xlow,ylow,lstgrd,
     &                numHoods,mioff,mjoff,
     &                qMerge,nvar,gradmx,gradmy)
        endif
c
c      apply limiter if requested. Go over all neighbors, do BJ
        if (nolimiter) go to 20

       if (limitTile .eq. 1) then
c         CELL based version
c         call limitTileBJ2(qmerge,gradmx,gradmy,
c    &                         xlow,ylow,dx,dy,
c    &                         irr,nvar,mitot,mjtot,lstgrd,
c    &                         i,j,k,mioff,mjoff,lwidth) 
c         limit all of grid all at once
          call limitTileGradientBJ(qmerge,gradmx,gradmy,
     &                         xlow,ylow,dx,dy,
     &                         mioff,mjoff, 
     &                         irr,nvar,mitot,mjtot,lstgrd,lwidth)

        else
           call limitTileGradientLP(qmerge,gradmx,gradmy,
     &                           xlow,ylow,dx,dy,irr,lwidth,
     &                           nvar,mitot,mjtot,lstgrd,mptr,
     &                           lpChoice,mioff,mjoff,
     &                           nborList,nborCount)
       endif
        
 20     continue

c
      ! redo neighborhood calc as in makeNHoods but putting merged vals INTO cells
      ! instead of getting FROm cells. have to do it this way because a cell doesn't 
      ! know which cells contribute TO it, only the giving cell knows.


      valnew = 0.d0  !  all cells initialized to 0

c     need to look at some ghost cells because they may distribute
c     to the last real cell
c     dont use first and last cell since no good update in method
      do 51 j = 1, mjtot
      do 50 i = 1, mitot
          k = irr(i,j)
          if (k .eq. -1) then
             ! set valnew to 'robust' fake state
             valnew(:,i,j) = fakeStateCons
             go to 50  ! does not contribute
          endif
          if (k .eq. lstgrd .or. ar(k) .gt. areaMin) then
             valnew(:,i,j) = valnew(:,i,j)+qMerge(:,i,j)/numHoods(i,j)
             go to 50
          endif
          call getCellCentroid(lstgrd,i,j,xc,yc,
     &                         xlow,ylow,dx,dy,k)
          if (OUT_OF_RANGE(i,j)) then
             valnew(:,i,j) = qMerge(:,i,j)  ! copy what came in  
             go to 50 
          endif

          ! have a cut cell to fix
          do ic = 0, ncount(k)   ! check if should be +1,  AG/and I differ
             if (ic .eq. 0) then
               ioff = i
               joff = j
               koff = k
             else
               ioff = iidx(ic,k)
               joff = jidx(ic,k)
               koff = irr(ioff,joff)
             endif

             qm = 0.d0  
             if (koff .eq. lstgrd) then
               arr = dx*dy
               call makep(poly(1,1,lstgrd),ioff,joff,xlow,ylow,dx,dy)
               ntris = 2
             else
               arr = ar(koff)
               ivert = 1
               do while (poly(ivert+1,1,koff).ne.-11)
                  ivert = ivert + 1
               end do
               ntris = ivert - 3
             endif

             indx1 = 1
             x1 = poly(indx1,1,koff)
             y1 = poly(indx1,2,koff)
             do it = 1, ntris
               indx2 = it + 1
               indx3 = it + 2

               x2 = poly(indx2,1,koff)
               y2 = poly(indx2,2,koff)

               x3 = poly(indx3,1,koff)
               y3 = poly(indx3,2,koff)

               artri = triangle_area(x1,x2,x3,y1,y2,y3)

               do itq = 1, ntriquad
                  xval = x1*rtri(itq) + x2*stri(itq) +
     &                   x3*(1.d0-rtri(itq)-stri(itq))
                  yval = y1*rtri(itq) + y2*stri(itq) +
     &                   y3*(1.d0-rtri(itq)-stri(itq))

                  diffx = xval - xcentMerge(k)
                  diffy = yval - ycentMerge(k)

                  qm = qm + (artri/arr)*wtri(itq)*(qMerge(:,i,j) +
     &             diffx * gradmx(:,k) +
     &             diffy * gradmy(:,k) +
     &           0.5d0*gradmxx(:,k)*(diffx**2 - (dx**2)*qmshifts(1,k)) +
     &              gradmxy(:,k)*(diffx*diffy - dx*dy*qmshifts(2,k)) +
     &           0.5d0*gradmyy(:,k)*(diffy**2 - (dy**2)*qmshifts(3,k)))

               end do ! loop over num terms in quadrature rule 
             end do ! loop over ntris per polyhedra

             valnew(:,ioff,joff) = valnew(:,ioff,joff) + 
     &                             qm(:)/numHoods(ioff,joff)

          end do

 50    continue
 51    continue
c
       q = valnew
c
c      check positivity
       call checkPhys(q,irr,mitot,mjtot,mptr,istage,lstgrd,
     .                   'from SRD_merge',lwidth+1,mitot-lwidth,
     .                   lwidth+1,mjtot-lwidth)

       return
       end
