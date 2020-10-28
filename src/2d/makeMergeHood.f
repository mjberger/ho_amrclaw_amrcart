c
c   -------------------------------------------------------------
c
      subroutine makeMergeHood(irr,numHoods,mitot,mjtot,lwidth,lstgrd,
     &                         xlow,ylow,hx,hy)

      use amr_module

      implicit real*8 (a-h, o-z)
      dimension irr(mitot,mjtot), numHoods(mitot,mjtot)
      include "cuserdt.i"

c
c setup merging neighborhoods for cut cells - how far in each direction
c
       areaMin = areaFrac*hx*hy
       numHoods = 1  !initialization  - every cell is its own nhood, at least, maybe more

       do j = 1,mjtot
       do i = 1,mitot
          k = irr(i,j)
          if (k .eq. -1 .or. k .eq. lstgrd) cycle

          iidx(1,k) = i
          jidx(1,k) = j

          icurr = i
          jcurr = j

          vqmerge = ar(k)  ! initial vol for merged hood is this cell
          if (ar(k) .ge. areaMin) cycle

          call determineDirection(irr,mitot,mjtot,lwidth,nvar,i,j,
     &         idir,areaMin)
          if (idir .eq. -1) then
             write(*,*)"Could not find merging direction. "
             write(*,*)"Should implement centered merge "
             stop
          endif
      
         ncount(k) = 0 !count how many cells BESIDES you are in your own nhood
         do while (vqmerge < areaMin) 
            if (idir .eq. 1) then
               icurr = icurr - 1
            else if (idir .eq. 2) then
               icurr = icurr + 1
            else if (idir .eq. 3) then
               jcurr = jcurr - 1
            else if (idir .eq. 4) then
             jcurr = jcurr + 1
            endif

            kcurr = irr(icurr,jcurr)
            if (kcurr .eq. -1) then
               write(*,*)"oops - hit solid wall"
               return
            endif

            ncount(k) = ncount(k) + 1
            iidx(ncount(k),k) = icurr
            jidx(ncount(k),k) = jcurr
            numHoods(icurr,jcurr) = numHoods(icurr,jcurr) + 1
            vqmerge = vqmerge + ar(kcurr)
         end do

       end do
       end do

c
c  now that have neighborhood, compute vol and centroids
c
       do j = 1, mjtot
       do i = 1, mitot
          k = irr(i,j)
          if (k .eq. -1 .or. k .eq. lstgrd) cycle

          call getCellCentroid(lstgrd,i,j,xc,yc,xlow,ylow,
     &                     hx,hy,k)

          volMerge(k) = ar(k)/numHoods(i,j)
          xcentMerge(k) = xc*ar(k)/numHoods(i,j)
          ycentMerge(k) = yc*ar(k)/numHoods(i,j)
          do ic = 1, ncount(k)
             icurr = iidx(ic,k)
             jcurr = jidx(ic,k)
             kcurr = irr(icurr,jcurr)
             call getCellCentroid(lstgrd,icurr,jcurr,xc,yc,xlow,ylow,
     &                        hx,hy,kcurr)
             volMerge(k) = volMerge(k) + ar(kcurr)/numHoods(icurr,jcurr)
             xcentMerge(k) = xcentMerge(k) + 
     &                       xc*ar(kcurr)/numHoods(icurr,jcurr)
             ycentMerge(k) = ycentMerge(k) + 
     &                       yc*ar(kcurr)/numHoods(icurr,jcurr)
          end do

          xcentMerge(k) = xcentMerge(k) / volMerge(k)
          ycentMerge(k) = ycentMerge(k) / volMerge(k)

       end do
       end do

       ncount(lstgrd) = 0   ! for regular cells, remember doesnt include cell itself
       volMerge(1) = hx * hy

       return
       end

