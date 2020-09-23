      subroutine limiter_srd(qp, qx, qy, mitot, mjtot, irr, nvar, hx,
     . hy, lwidth, lstgrd, xlow, ylow)
       implicit double precision(a-h,o-z)


       common /order2/ ssw, primitive, nolimiter
       dimension nlist(25,2)
       logical primitive
       include "./reconinfo.i"
       include "cirr.i"
       dimension qp(mitot,mjtot,4),qx(irrsize,4),
     &           qy(irrsize,4),irr(mitot,mjtot)



       dimension dumax(4),dumin(4),phimin(4)
       dimension graddot(4),alpha(4),recon(4)
       logical IS_GHOST, verbose/.true./
       IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                  j .le. lwidth .or. j .gt. mjtot-lwidth)




!        qx = 0.d0
!        qy = 0.d0
!        return

       if(primitive .eqv. .false.) then
        print *, "PROBLEM: THE CURRENT LIMITER IMPLEMENTATION ONLY WORKS
     .            WHEN RECONSTRUCTING IN PRIMITIVE VARIABLES."
       endif



      do 30 j = lwidth+1, mjtot-lwidth
      do 30 i = lwidth+1, mitot-lwidth

!          if( i .eq. 12 .and. j .eq. 14) then
!          print *, "here"
!          endif


          k = irr(i,j)
          if (k .eq. -1 .or. k .eq. lstgrd) goto 30 ! solid/full cells have no need for a gradient

          ! find max and min needed for BJ limiting
          dumax = 0.d0
          dumin = 0.d0
          do 31 joff = -mjoff(i,j), mjoff(i,j)
          do 31 ioff = -mioff(i,j), mioff(i,j)
            if (ioff .eq. 0 .and. joff .eq. 0) go to 31
            if (IS_GHOST(i+ioff,j+joff)) go to 31

            koff = irr(i+ioff,j+joff)
            if (koff .eq. -1) go to 31
            dumax = max(dumax,qp(i+ioff,j+joff,:)-qp(i,j,:))
            dumin = min(dumin,qp(i+ioff,j+joff,:)-qp(i,j,:))
  31      continue

          phimin = 1.d0
          do 32 joff = -mjoff(i,j), mjoff(i,j)
          do 32 ioff = -mioff(i,j), mioff(i,j)
             if (ioff .eq. 0 .and. joff .eq. 0) go to 32 ! no eqn to solve
             if (IS_GHOST(i+ioff,j+joff)) go to 32
             koff = irr(i+ioff,j+joff)
             if (koff .eq. -1) go to 32

             if(koff .eq. lstgrd) then
              call getCellCentroid(lstgrd,i+ioff,j+joff,xc,yc,xlow,
     .                       ylow,hx,hy,koff)
             else
                 xc = xcentMerge(koff)
                 yc = ycentMerge(koff)
             endif

              diffx = xc-xcentMerge(k)
              diffy = yc-ycentMerge(k)
              graddot  = qx(k,:)*diffx + qy(k,:)*diffy
              recon = qp(i,j,:) + graddot
                do m = 1,4
                   if (graddot(m) > 0.d0) then
                      alpha(m) = min(1.d0, dumax(m)/graddot(m))
                   else if (graddot(m) < 0.d0) then
                      alpha(m) = min(1.d0, dumin(m)/graddot(m))
                   else
                      alpha(m) = 1.d0
                   endif
                end do
                ! one last check for positivity
!                if (recon(1) .le. 0.d0) alpha = 0.d0
!                velsq = recon(2)**2+recon(3)**2
!                press = .4d0*(recon(4)-0.5d0*velsq/recon(1))
!                if (press .le. 0.d0) alpha = 0.d0
                phimin = min(phimin, alpha)
  32      continue
          qx(k,:) = qx(k,:)*phimin(:)
          qy(k,:) = qy(k,:)*phimin(:)

  30  continue

110   continue









      ! final check for positivity
      do 310 j = lwidth+1, mjtot-lwidth
      do 310 i = lwidth+1, mitot-lwidth
          k = irr(i,j)
          if (k .eq. -1 .or. k .eq. lstgrd) goto 310 ! solid/full cells have no need for a gradient

          alpha = 1.d0
          do 27 ic = 1, ncount(k)
             icurr = iidx(k,ic)
             jcurr = jidx(k,ic)
             koff = irr(icurr,jcurr)
             call getCellCentroid(lstgrd,icurr,jcurr,xc,yc,xlow,
     .                            ylow,hx,hy,koff)
             diffx = xc-xcentMerge(k)
             diffy = yc-ycentMerge(k)
             graddot  = qx(k,:)*diffx + qy(k,:)*diffy
             recon = qp(i,j,:) + graddot
             ! one last check for positivity
             if (recon(1) .le. 0.d0) alpha = 0.d0
             velsq = recon(2)**2+recon(3)**2
             press = .4d0*(recon(4)-0.5d0*velsq/recon(1))
             if (press .le. 0.d0) alpha = 0.d0
  27      continue
          qx(k,:) = qx(k,:)*alpha(:)
          qy(k,:) = qy(k,:)*alpha(:)
  310 continue









      ! final check for positivity
      do 111 j = lwidth+1, mjtot-lwidth
      do 111 i = lwidth+1, mitot-lwidth
          k = irr(i,j)
          if (k .eq. -1 .or. k .eq. lstgrd) goto 111 ! solid/full cells have no need for a gradient

          alpha = 1.d0
          do 17 ic = 1, ncount(k)
             icurr = iidx(k,ic)
             jcurr = jidx(k,ic)
             koff = irr(icurr,jcurr)
             call getCellCentroid(lstgrd,icurr,jcurr,xc,yc,xlow,
     .                            ylow,hx,hy,koff)
             diffx = xc-xcentMerge(k)
             diffy = yc-ycentMerge(k)
             recon = qp(i,j,:) + qx(k,:)*diffx + qy(k,:)*diffy
             ! one last check for positivity
             if (recon(1) .le. 0.d0) alpha = 0.d0
             velsq = recon(2)**2+recon(3)**2
             press = .4d0*(recon(4)-0.5d0*velsq/recon(1))
             if (press .le. 0.d0) alpha = 0.d0
  17      continue

          if(alpha(1) .eq. 0) then
          print *, "FINAL POSITIVITY CHECK FAILED."
          endif

  111 continue







      return
      end



