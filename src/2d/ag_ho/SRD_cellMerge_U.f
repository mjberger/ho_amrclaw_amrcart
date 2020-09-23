c
c ---------------------------------------------------------------------
c
       subroutine SRD_cellMerge_U(q,dlimit,nvar,irr,mitot,mjtot,qx,qy,
     .                       lstgrd,
     .                       dx,dy,lwidth,xlow,ylow,istage,
     .                       ncount,numHoods)

       implicit double precision (a-h, o-z)
       include "cirr.i"

       dimension q(mitot,mjtot,nvar)
       dimension irr(mitot,mjtot)
       dimension qx(mitot,mjtot,nvar), qy(mitot,mjtot,nvar)
       dimension gradmx(mitot,mjtot,nvar), gradmy(mitot,mjtot,nvar)
       dimension dlimit(mitot,mjtot,nvar)
       dimension perimeter(mitot,mjtot)


       dimension valnew(mitot,mjtot,nvar), volMerge(mitot,mjtot)
       dimension qMerge(mitot,mjtot,nvar), numHoods(mitot,mjtot)
       dimension nCount(mitot,mjtot)
       dimension xcentMerge(mitot,mjtot),ycentMerge(mitot,mjtot)
       dimension fakeState(nvar), qm(nvar), rhs(2,nvar)
       dimension dumax(nvar),dumin(nvar),phimin(nvar)
       dimension graddot(nvar),alpha(nvar),recon(nvar)
       dimension a(2,2),b(2)
       real*8 minmod

       logical IS_GHOST, IS_FAR_GHOST, verbose/.true./
       logical quad, nolimiter
       logical BAD_STATE
       common /order2/ ssw, quad, nolimiter

c  this next statement considers a ghost cell to be anything beyond the
c  loop indices, since they change according to the stage
       IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                  j .le. lwidth .or. j .gt. mjtot-lwidth)
C        BAD_STATE = ((qm(1) .le. 0.d0) .or. 
C      .  (.4d0*(qm(4)- 0.5d0*(qm(2)**2+qm(3)**2)/qm(1)) .le.0.d0))

c :::::::::::::;:::
c
c   try cell merging and reconstruction to stabilize updates in small cut cells
c   calculate all updates using provisional values, then actually do the update
c   (as in, jacobi rather than gauss seidel)
c ::::::::::::::::

c
       ar(-1) = 0.d0        ! zero out area of solid cells for this loop
       ar(lstgrd) = dx*dy   ! area of regular grid cell 
       qMerge   = 0.d0


       ! put in 'normal' values to prevent errors e.g. in converting to prim vars
       fakeState(1) = 1.d0
       fakeState(2) = 0.d0
       fakeState(3) = 0.d0
       fakeState(4) = 2.5d0

c     first make neighborhoods - need count for each cells, and width (nhood above)
c     nCount is size of neighborhood, numHoods is number of merged nhoods each cells is in
      call makeNHood(volMerge,xcentMerge,ycentMerge,ncount,irr,numHoods,
     .               mitot,mjtot,lwidth,lstgrd,xlow,ylow,dx,dy)   

      call getPerimeter(lstgrd,mitot,mjtot, lwidth, ncount, irr, dx, dy,
     .                  xlow,ylow, perimeter)







       if (verbose) then
          totmass =  bigconck(q,irr,mitot,mjtot,lwidth,nvar)
          write(*,911) totmass
 911      format(/,"         SRDmass before redistribution is ",e30.20)
       endif

c       form qMerge vals 
c       do 10 j = lwidth+1, mjtot-lwidth
c       do 10 i = lwidth+1, mitot-lwidth
        do 10 j = 1, mjtot
        do 10 i = 1, mitot
            if (irr(i,j) .eq. -1) go to 10 ! no solid cells
            if (irr(i,j) .eq. lstgrd .or. IS_GHOST(i,j)) then 
              qMerge(i,j,:) = q(i,j,:)
              go to 10 
            endif
c
c           sum in blocks of size 2*ncount on a size using only valid cells 
c
             do 27 joff = -ncount(i,j), ncount(i,j)
             do 27 ioff = -ncount(i,j), ncount(i,j)
                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 27  ! solid cells dont contribute to neighborhood
                if (IS_GHOST(i+ioff,j+joff)) go to 27  ! nor ghost cells
c               count this cell 
                qMerge(i,j,:) = qMerge(i,j,:) + ar(koff)*
     .                        q(i+ioff,j+joff,:)/numHoods(i+ioff,j+joff)
 27          continue
             qMerge(i,j,:) = qMerge(i,j,:) / volMerge(i,j)

 10     continue

        ! gradient of merge neighborhoods, initialized to 0. set using neighboring merged tiels
        gradmx = 0.d0
        gradmy = 0.d0

        do 20 j = lwidth+1, mjtot-lwidth
        do 20 i = lwidth+1, mitot-lwidth



            k = irr(i,j)
            if (k .eq. -1) go to 20 ! solid cells have no gradient
            !! numHood=1 is its own neighborhood, no gradient needed
c            if (numHoods(i,j) .eq. 1) go to 20 ! regular cell doesnt participate in anything
            !!if (k .eq. lstgrd .and. ncount(i,j) .eq. 0) then 
            if (ncount(i,j) .eq. 0) then 
            ! these should be all regular cells that are left
               nco = 1 ! need to use some kind of nhood to make gradient
            else 
               nco = ncount(i,j)
            endif

            ! can reconstruct high order poly here, but for now do linear reconstruction
            ! solving least squares problem with all neighboring merged vals. Should preprocess the matrix and save
            rhs = 0.d0 ! initialize for accumulation
            a = 0.d0
            x0 = xcentMerge(i,j)
            y0 = ycentMerge(i,j)
            do 22 joff = -nco, nco
            do 22 ioff = -nco, nco
                if (ioff .eq. 0 .and. joff .eq. 0) go to 22 ! no eqn to solve
                if (IS_GHOST(i+ioff,j+joff)) go to 22
                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 22
               deltax = xcentMerge(i+ioff,j+joff) - x0
               deltay = ycentMerge(i+ioff,j+joff) - y0
               a(1,1) = a(1,1) + deltax*deltax
               a(1,2) = a(1,2) + deltax*deltay
               a(2,2) = a(2,2) + deltay*deltay
               rhs(1,:) = rhs(1,:) + deltax * 
     .                    (qMerge(i+ioff,j+joff,:) - qMerge(i,j,:))
               rhs(2,:) = rhs(2,:) + deltay * 
     .                    (qMerge(i+ioff,j+joff,:) - qMerge(i,j,:))
 22          continue

             ! solve a^t*a  * grad = a^t rhs.  First form cholesky factors.
             ! will have to add robustness checks
             c11 = sqrt(a(1,1))
             c12 = a(1,2) / c11
             c22 = sqrt(a(2,2) - c12**2)

             ! now back solve (C^t C = rhs of A^tdu ) to get x and y gradient for all variables
             do m = 1, nvar
               b(1) = rhs(1,m) / c11
               b(2) = (rhs(2,m) - c12*b(1))/c22
               gradmy(i,j,m) = b(2) / c22
               gradmx(i,j,m) = (b(1) - c12*gradmy(i,j,m))/c11
             end do




 20     continue

c
c      apply limiter if requested. Go over all neighbors, do BJ
        if (.false.) go to 35

        do 30 j = lwidth+1, mjtot-lwidth
        do 30 i = lwidth+1, mitot-lwidth
            k = irr(i,j)
            if (k .eq. -1) go to 30 ! solid cells have no gradient
c            if (numHoods(i,j) .eq. 1) go to 30  ! CHECK THAT NOTHING TO DO AND VAL NOT CHANGED
            nco = ncount(i,j)

            ! find max and min needed for BJ limiting
            dumax = 0.d0 
            dumin = 0.d0
            do 31 joff = -nco, nco
            do 31 ioff = -nco, nco
              if (ioff .eq. 0 .and. joff .eq. 0) go to 31
              if (IS_GHOST(i+ioff,j+joff)) go to 31
              koff = irr(i+ioff,j+joff)
              if (koff .eq. -1) go to 31
              dumax = max(dumax,qmerge(i+ioff,j+joff,:)-qmerge(i,j,:))
              dumin = min(dumin,qmerge(i+ioff,j+joff,:)-qmerge(i,j,:))
 31         continue

            phimin = 1.d0
            do 32 joff = -nco, nco
            do 32 ioff = -nco, nco
                if (ioff .eq. 0 .and. joff .eq. 0) go to 32 ! no eqn to solve
                if (IS_GHOST(i+ioff,j+joff)) go to 32
                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 32
                ! use merged val at ioff,joff to limit
                !diffx = xcentMerge(i+ioff,j+joff)-xcentMerge(i,j)
                !diffy = ycentMerge(i+ioff,j+joff)-ycentMerge(i,j)
                ! these next lines limit at cell where will evaluate
                ! previous lines limit at cells used in making gradient
                call getCellCentroid(lstgrd,i+ioff,j+joff,xc,yc,xlow,
     .                          ylow,dx,dy,koff)
                diffx = xc-xcentMerge(i,j)
                diffy = yc-ycentMerge(i,j)
                graddot  = gradmx(i,j,:)*diffx + gradmy(i,j,:)*diffy
                recon = qMerge(i,j,:) + graddot
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
                  if (recon(1) .le. 0.d0) alpha = 0.d0
                  velsq = recon(2)**2+recon(3)**2
                  press = .4d0*(recon(4)-0.5d0*velsq/recon(1))
                  if (press .le. 0.d0) alpha = 0.d0
                  phimin = min(phimin, alpha)
 32         continue
            gradmx(i,j,:) = gradmx(i,j,:)*phimin(:)
            gradmy(i,j,:) = gradmy(i,j,:)*phimin(:)
            dlimit(i,j,:) = phimin(:)
 30     continue

 
c      gradmx = 0.d0  
c      gradmy = 0.d0 
c
      ! redo neighborhood calc as in makeNHoods but putting merged vals INTO cells
      ! instead of getting FROm cells. have to do it this way because a cell doesn't 
      ! know which cells contribute TO it, only the giving cell knows.

 35   continue


      valnew = 0.d0  !  all cells initialized to 0

      do 50 j = 1, mjtot
      do 50 i = 1, mitot

c       if(i  .eq. 5 .and. j .eq. 19) then
c       print *,"here"
c       endif

          k = irr(i,j)
          if (k .eq. -1) then
             ! set valnew to 'robust' fake state
             valnew(i,j,:) = fakeState
             go to 50  ! does not contribute
          endif
          if (IS_GHOST(i,j)) then
             valnew(i,j,:) = qMerge(i,j,:)  ! copy from what came in
             go to 50 
          endif
          !!if (k .eq. lstgrd .and. numHoods(i,j) .eq. 1) then ! doesnt participate in anything
          !!if (numHoods(i,j) .eq. 1) then ! doesnt participate in anything
          !!   valnew(i,j,:) = qMerge(i,j,:) ! takes its own state
          !!   go to 50  
          !!endif


          do 40 joff = -ncount(i,j), ncount(i,j)
          do 40 ioff = -ncount(i,j), ncount(i,j)
             if (IS_GHOST(i+ioff,j+joff)) go to 40  
             koff = irr(i+ioff,j+joff)
             if (koff .eq. -1) go to 40  

             call getCellCentroid(lstgrd,i+ioff,j+joff,xc,yc,xlow,
     .                            ylow,dx,dy,koff)
             qm(:) = qMerge(i,j,:) + (xc-xcentMerge(i,j))*gradmx(i,j,:) 
     .                             + (yc-ycentMerge(i,j))*gradmy(i,j,:)






             !!if (BAD_STATE) then
c              pr = .4d0*(qm(4)- 0.5d0*(qm(2)**2+qm(3)**2)/qm(1))
c              if ((qm(1) .le. 0.d0) .or. (pr .le. 0.d0)) then
c                 write(*,*)" should not happen"
c              endif
             valnew(i+ioff,j+joff,:) = valnew(i+ioff,j+joff,:) + 
     .              qm(:)/numHoods(i+ioff,j+joff)
 40       continue

 50    continue

c
       q = valnew
c
       if (verbose) then
          totmass2 =  bigconck(q,irr,mitot,mjtot,lwidth,nvar)
          dif = totmass2 - totmass

          write(*,912) totmass2,dif
 912      format("         SRDmass after  redistribution is ",e30.20,
     .           "  dif is ",e15.7)
       endif

       return
       end
c
c ----------------------------------------------------------------------------
c
      function bigconck(q,irr,mitot,mjtot,lwidth,nvar)
c
      implicit double precision (a-h,o-z)
      include "cirr.i"

      dimension q(mitot,mjtot,nvar), irr(mitot,mjtot)

      totmass = 0.d0
      ar(-1)  = 0.d0

      do 10 j = lwidth+1, mjtot-lwidth
      do 10 i = lwidth+1, mitot-lwidth

         k = irr(i,j)
         totmass = totmass + q(i,j,1)*ar(k)

 10   continue

      bigconck = totmass

      return
      end

c
c ----------------------------------------------------------------------------
c
      subroutine makeNHood(volMerge,xcentMerge,ycentMerge,ncount,irr,
     .                     numHoods,mitot,mjtot,lwidth,lstgrd,xlow,ylow,
     ,                     dx,dy)

      implicit double precision (a-h, o-z)
      include "cirr.i"

      dimension numHoods(mitot,mjtot), volMerge(mitot,mjtot)
      dimension xcentMerge(mitot,mjtot), ycentMerge(mitot,mjtot)
      dimension ncount(mitot,mjtot), irr(mitot,mjtot)

      logical IS_GHOST, firstTimeThru
      IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                 j .le. lwidth .or. j .gt. mjtot-lwidth)

      ! merge until vqmerge at least this big (analogous to 1d where left and right nhoods each dx
      areaMin = 0.5d0*ar(lstgrd)
      numHoods = 0  ! initialize, loop below will add each cell to its own nhood
      ncount = 0
      ar(-1) = 0.d0  ! reset here to remind us

      do 10 j = lwidth+1, mjtot-lwidth
      do 10 i = lwidth+1, mitot-lwidth

         k = irr(i,j)  
         if (k .eq. -1) go to 10
         if (k .eq. lstgrd) then
            numHoods(i,j) =  numHoods(i,j) + 1
            go to 10 ! a full  flow cell is its own merged neighborhood
         endif
         vqmerge = 0.d0
         firstTimeThru = .true.
         nco = 0   ! initial size of neighborhood, from -1 to 1 square centered on cell

            !do while (dcell_h < dlenMin)
            do while (vqmerge < areaMin)
               do 15 joff = -nco, nco
               do 15 ioff = -nco, nco
                   if (IS_GHOST(i+ioff,j+joff)) go to 15  
                   koff = irr(i+ioff,j+joff)
                   if (koff .eq. -1) go to 15  ! solid cells dont help
                   vqmerge = vqmerge + ar(koff)
                   if (firstTimeThru) then ! count everybody
                      numHoods(i+ioff,j+joff)=numHoods(i+ioff,j+joff)+1
                   else ! only add to new cells-on newly enlarged nhood border
                     if (abs(ioff).eq. nco .or. abs(joff).eq.nco)
     .                numHoods(i+ioff,j+joff)=numHoods(i+ioff,j+joff)+1
                   endif
 15             continue

                !if (dcell_h .ge. dlenMin-1.e-10) then
                if (vqmerge > areaMin) then
                   ncount(i,j) = nco
                   go to 10
                else   ! redo with larger neighborhood
                   nco = nco + 1
                   vqmerge = 0.d0
                   firstTimeThru = .false.
                endif
            end do
             
 10   continue      


c     needed number of neighbhoods to compute volMerge = which is not
c     the real volume of the merging neighborhood
c

!   initialize array with  most common case, overwritten below as needed
      volMerge = ar(lstgrd) 
      xcentMerge = 0.d0      
      ycentMerge = 0.d0

      do 20 j = 1, mjtot
      do 20 i = 1, mitot
         k = irr(i,j)  
         if (k .eq. lstgrd) then
             call getCellCentroid(lstgrd,i,j,xcentMerge(i,j),
     .                            ycentMerge(i,j),
     .                            xlow,ylow,dx,dy,k)
             go to 20 ! a full  flow cell is its own merged neighborhood
         endif
         if (IS_GHOST(i,j)) then
            volMerge(i,j) = 0.d0
            go to 20
         endif
         if (k .eq. -1) then
            volMerge(i,j) = 0.d0
            go to 20
         endif

         vmerge = 0.d0  ! diff variable than above, weighted by numhoods
         xcent = 0.d0
         ycent = 0.d0
         nco = ncount(i,j)
         do 25 joff = -nco, nco
         do 25 ioff = -nco, nco
            if (IS_GHOST(i+ioff,j+joff)) go to 25  
            koff = irr(i+ioff,j+joff)
            if (koff .eq. -1) go to 25 ! solid cells dont help
            call getCellCentroid(lstgrd,i+ioff,j+joff,xc,yc,xlow,ylow,
     .                           dx,dy,koff)
            vmerge = vmerge +  ar(koff)/numHoods(i+ioff,j+joff)
            xcent = xcent + xc*ar(koff)/numHoods(i+ioff,j+joff)
            ycent = ycent + yc*ar(koff)/numHoods(i+ioff,j+joff)

 25      continue
         volMerge(i,j) = vmerge
         xcentMerge(i,j) = xcent/vmerge
         ycentMerge(i,j) = ycent/vmerge
c         print *,ncount(i,j)
 20   continue



      return
      end
c
c -------------------------------------------------------------------------
c
      subroutine getCellCentroid(lstgrd,i,j,xc,yc,xlow,ylow,dx,dy,k)

      implicit double precision (a-h,o-z)
      include "cirr.i"

      if (k .eq. lstgrd) then
         xc = xlow + (i-0.5d0)*dx
         yc = ylow + (j-0.5d0)*dy
      else if (k .eq. -1) then
         xc = 0.d0
         yc = 0.d0
      else
         xc = xcirr(k)
         yc = ycirr(k)
      endif

      return
      end
c
c -------------------------------------------------------------------
c
      double precision function minmod(a,b)
      double precision a, b, ans


      if (a*b <= 0.d0) then
        ans = 0.d0
      else if (a .gt. 0.d0 .and. b .gt. 0.d0) then
        ans = min(a,b)
      else
        ans = max(a,b)
      endif

      minmod = ans
      return
      end


c
c -------------------------------------------------------------------------
c
      subroutine getPerimeter(lstgrd, mitot, mjtot, lwidth, ncount, irr,
     .                        dx, dy, xlow, ylow, perimeter)

      implicit double precision (a-h,o-z)
      include "cirr.i"
      dimension ncount(mitot,mjtot), irr(mitot,mjtot)
      logical IS_GHOST, IS_SAME
      IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth) .or.
     .                (j .le. lwidth .or. j .gt. mjtot-lwidth)
      IS_SAME(x1,x2,x3,x4,y1,y2,y3,y4) =
     .                          ( (abs( x1 - x3 ) .le. 1.e-10)   .and.
     .                            (abs( y1 - y3 ) .le. 1.e-10)   .and.
     .                            (abs( x2 - x4 ) .le. 1.e-10)   .and.
     .                            (abs( y2 - y4 ) .le. 1.e-10) ) .or.
     .                          ( (abs( x1 - x4 ) .le. 1.e-10)   .and.
     .                            (abs( y1 - y4 ) .le. 1.e-10)   .and.
     .                            (abs( x2 - x3 ) .le. 1.e-10)   .and.
     .                            (abs( y2 - y3 ) .le. 1.e-10) )
      dimension faces(4,100)
      dimension idup(100)
      dimension perimeter(mitot,mjtot)

      perimeter = 0.d0

      do 20 j = 1, mjtot
      do 20 i = 1, mitot
        nco = ncount(i,j)
        k = irr(i, j)

        if(IS_GHOST(i,j) .or. k .eq. -1) goto 20

        ifaces = 0
        faces = -1

        do 25 joff = -nco, nco
        do 25 ioff = -nco, nco
          k = irr(i + ioff, j + joff)

          ! do not include ghost and completely skip solid cells
          if(IS_GHOST(i+ioff,j+joff) .or. k .eq. -1)   goto 25

          ! if regular cell
          if(k .eq. lstgrd) then
            x1 = xlow + dfloat(i + ioff - 1)*dx
            y1 = ylow + dfloat(j + joff - 1)*dy
            x2 = xlow + dfloat(i + ioff - 1)*dx
            y2 = ylow + dfloat(j + joff)*dy
            x3 = xlow + dfloat(i + ioff)*dx
            y3 = ylow + dfloat(j + joff)*dy
            x4 = xlow + dfloat(i + ioff)*dx
            y4 = ylow + dfloat(j + joff - 1)*dy

            ifaces = ifaces + 1
            faces(1,ifaces) = x1
            faces(2,ifaces) = y1
            faces(3,ifaces) = x2
            faces(4,ifaces) = y2

            ifaces = ifaces + 1
            faces(1,ifaces) = x2
            faces(2,ifaces) = y2
            faces(3,ifaces) = x3
            faces(4,ifaces) = y3

            ifaces = ifaces + 1
            faces(1,ifaces) = x3
            faces(2,ifaces) = y3
            faces(3,ifaces) = x4
            faces(4,ifaces) = y4

            ifaces = ifaces + 1
            faces(1,ifaces) = x4
            faces(2,ifaces) = y4
            faces(3,ifaces) = x1
            faces(4,ifaces) = y1
            goto 25
           endif

          ! if irregular cell
          kside = 1
          do 28 while (poly(kside+1,1,k) .ne. -11.)
            ifaces = ifaces + 1
            faces(1,ifaces) = poly(kside,1,k)  !x1
            faces(2,ifaces) = poly(kside,2,k)  !y1
            faces(3,ifaces) = poly(kside+1,1,k)!x2
            faces(4,ifaces) = poly(kside+1,2,k)!y2
            kside = kside + 1
 28       continue

 25      continue





        idup = 0
              ! faces found
         do 26 icount = 1,ifaces
         do 26 jcount = 1,ifaces
             if (icount .eq. jcount) goto 26

             x1 = faces(1,icount)
             y1 = faces(2,icount)
             x2 = faces(3,icount)
             y2 = faces(4,icount)

             x3 = faces(1,jcount)
             y3 = faces(2,jcount)
             x4 = faces(3,jcount)
             y4 = faces(4,jcount)
             if(IS_SAME(x1,x2,x3,x4,y1,y2,y3,y4)) then
              idup(icount) =  1
              idup(jcount) =  1
             endif
 26       continue


        do 27 icount = 1,ifaces
            if(idup(icount) .eq. 0) then
                 x1 = faces(1,icount)
                 y1 = faces(2,icount)
                 x2 = faces(3,icount)
                 y2 = faces(4,icount)
                 dlength = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) )
                 perimeter(i,j) = perimeter(i,j) + dlength
            endif
 27     continue

 20   continue
      ! loop over grid


      end

c
c -------------------------------------------------------------------------
c
      function cellPerimeter(lstgrd, mitot, mjtot, lwidth, nco,
     . irr,dx, dy, xlow, ylow, i,j, poly2)

      implicit double precision (a-h,o-z)
      include "cirr.i"
      dimension ncount(mitot,mjtot), irr(mitot,mjtot)
      logical IS_GHOST, IS_SAME, ARE_EQUAL
      IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth) .or.
     .                (j .le. lwidth .or. j .gt. mjtot-lwidth)
      IS_SAME(x1,x2,x3,x4,y1,y2,y3,y4) =
     .                          ( (abs( x1 - x3 ) .le. 1.e-10)   .and.
     .                            (abs( y1 - y3 ) .le. 1.e-10)   .and.
     .                            (abs( x2 - x4 ) .le. 1.e-10)   .and.
     .                            (abs( y2 - y4 ) .le. 1.e-10) ) .or.
     .                          ( (abs( x1 - x4 ) .le. 1.e-10)   .and.
     .                            (abs( y1 - y4 ) .le. 1.e-10)   .and.
     .                            (abs( x2 - x3 ) .le. 1.e-10)   .and.
     .                            (abs( y2 - y3 ) .le. 1.e-10) )

      ARE_EQUAL(x1,x2,y1,y2) =  ( (abs( x1 - x2 ) .le. 1.e-10)   .and.
     .                            (abs( y1 - y2 ) .le. 1.e-10) )

      dimension faces(4,100)
      dimension idup(100)
      dimension poly2(100,2,2)
      dimension poly3(100,2,2)
      poly2 = -11
      poly3 = -11

      k = irr(i, j)

      if(IS_GHOST(i,j) .or. k .eq. -1) then
        cellPerimeter = -1
        return
      endif

      ifaces = 0
      faces = -1

      do 25 joff = -nco, nco
      do 25 ioff = -nco, nco
        k = irr(i + ioff, j + joff)

        ! do not include ghost and completely skip solid cells
        if(IS_GHOST(i+ioff,j+joff) .or. k .eq. -1)   goto 25

        ! if regular cell
        if(k .eq. lstgrd) then
          x1 = xlow + dfloat(i + ioff - 1)*dx
          y1 = ylow + dfloat(j + joff - 1)*dy
          x2 = xlow + dfloat(i + ioff - 1)*dx
          y2 = ylow + dfloat(j + joff)*dy
          x3 = xlow + dfloat(i + ioff)*dx
          y3 = ylow + dfloat(j + joff)*dy
          x4 = xlow + dfloat(i + ioff)*dx
          y4 = ylow + dfloat(j + joff - 1)*dy

          ifaces = ifaces + 1
          faces(1,ifaces) = x1
          faces(2,ifaces) = y1
          faces(3,ifaces) = x2
          faces(4,ifaces) = y2

          ifaces = ifaces + 1
          faces(1,ifaces) = x2
          faces(2,ifaces) = y2
          faces(3,ifaces) = x3
          faces(4,ifaces) = y3

          ifaces = ifaces + 1
          faces(1,ifaces) = x3
          faces(2,ifaces) = y3
          faces(3,ifaces) = x4
          faces(4,ifaces) = y4

          ifaces = ifaces + 1
          faces(1,ifaces) = x4
          faces(2,ifaces) = y4
          faces(3,ifaces) = x1
          faces(4,ifaces) = y1
          goto 25
         endif

        ! if irregular cell
        kside = 1
        do 28 while (poly(kside+1,1,k) .ne. -11.)
          ifaces = ifaces + 1
          faces(1,ifaces) = poly(kside,1,k)  !x1
          faces(2,ifaces) = poly(kside,2,k)  !y1
          faces(3,ifaces) = poly(kside+1,1,k)!x2
          faces(4,ifaces) = poly(kside+1,2,k)!y2
          kside = kside + 1
 28     continue

 25    continue





      idup = 0
            ! faces found
       do 26 icount = 1,ifaces
       do 26 jcount = 1,ifaces
           if (icount .eq. jcount) goto 26

           x1 = faces(1,icount)
           y1 = faces(2,icount)
           x2 = faces(3,icount)
           y2 = faces(4,icount)

           x3 = faces(1,jcount)
           y3 = faces(2,jcount)
           x4 = faces(3,jcount)
           y4 = faces(4,jcount)
           if(IS_SAME(x1,x2,x3,x4,y1,y2,y3,y4)) then
            idup(icount) =  1
            idup(jcount) =  1
           endif
 26     continue

      cellPerimeter = 0.d0
      icount2 = 1
      do 27 icount = 1,ifaces
          if(idup(icount) .eq. 0) then
               x1 = faces(1,icount)
               y1 = faces(2,icount)
               x2 = faces(3,icount)
               y2 = faces(4,icount)
               dlength = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) )
               cellPerimeter = cellPerimeter + dlength
               poly2(icount2,1,1) = x1
               poly2(icount2,2,1) = y1
               poly2(icount2,1,2) = x2
               poly2(icount2,2,2) = y2
               icount2 = icount2+1
          endif
 27   continue
      icount2 = icount2-1



      ! sort the faces
      do 29 icount = 1,icount2

c      if(icount2 .eq. 10) then
c      write(*,*) "HERE !!!!!!!!!!! SORTING INFO"
c      do 229 ip = 1,icount2
c        write(*,*) poly2(ip,1,1),poly2(ip,2,1),
c     .             poly2(ip,1,2),poly2(ip,2,2)
c 229   continue
c      endif

          x1 = poly2(icount,1,2)
          y1 = poly2(icount,2,2)
      do 30 jcount = icount+1, icount2
          ! find (x1,y1) in the first coordinate
          x2 = poly2(jcount,1,1)
          y2 = poly2(jcount,2,1)

          ! find (x1,y1) in the second coordinate
          x3 = poly2(jcount,1,2)
          y3 = poly2(jcount,2,2)
          if(ARE_EQUAL(x1,x2,y1,y2) .and. jcount .ne. icount+1) then ! swap icount + 1 with jcount
            temp1_x = poly2(jcount,1,1)
            temp1_y = poly2(jcount,2,1)
            temp2_x = poly2(jcount,1,2)
            temp2_y = poly2(jcount,2,2)

            poly2(jcount,1,1) = poly2(icount+1,1,1)
            poly2(jcount,2,1) = poly2(icount+1,2,1)
            poly2(jcount,1,2) = poly2(icount+1,1,2)
            poly2(jcount,2,2) = poly2(icount+1,2,2)

            poly2(icount+1,1,1) =  temp1_x
            poly2(icount+1,2,1) =  temp1_y
            poly2(icount+1,1,2) =  temp2_x
            poly2(icount+1,2,2) =  temp2_y
!            write(*,*) "swap",icount+1, jcount
            goto 29
          elseif(ARE_EQUAL(x1,x3,y1,y3) ) then ! flip jcount, then swap icount + 1 with jcount
          ! flip
            temp2_x =  poly2(jcount,1,2)
            temp2_y =  poly2(jcount,2,2)
            poly2(jcount,1,2)  = poly2(jcount,1,1)
            poly2(jcount,2,2)  = poly2(jcount,2,1)
            poly2(jcount,1,1)  = temp2_x
            poly2(jcount,2,1)  = temp2_y
!            write(*,*) "flip",jcount
          ! swap only if icount+1 .neq. jcount
            if(jcount .ne. icount+1) then
                temp1_x = poly2(jcount,1,1)
                temp1_y = poly2(jcount,2,1)
                temp2_x = poly2(jcount,1,2)
                temp2_y = poly2(jcount,2,2)

                poly2(jcount,1,1) = poly2(icount+1,1,1)
                poly2(jcount,2,1) = poly2(icount+1,2,1)
                poly2(jcount,1,2) = poly2(icount+1,1,2)
                poly2(jcount,2,2) = poly2(icount+1,2,2)

                poly2(icount+1,1,1) =  temp1_x
                poly2(icount+1,2,1) =  temp1_y
                poly2(icount+1,1,2) =  temp2_x
                poly2(icount+1,2,2) =  temp2_y
!                write(*,*) "swap",icount+1, jcount
            endif
            goto 29
          endif

!          write(*,*) "not found"
 30   continue
 29   continue

c      if(icount2 .eq. 10) then
c      write(*,*) "HERE !!!!!!!!!!! SORTING INFO"
c      do 229 ip = 1,icount2
c        write(*,*) poly2(ip,1,1),poly2(ip,2,1),
c     .             poly2(ip,1,2),poly2(ip,2,2)
c 229   continue
c      endif
      ! check if two edges are collinear
 999  continue
      do 31 icount = 1,icount2
c      write(*,*) "HERE !!!!!!!!!!! SORTING INFO"
c      do 229 ip = 1,icount2
c        write(*,*) poly2(ip,1,1),poly2(ip,2,1),
c     .             poly2(ip,1,2),poly2(ip,2,2)
c 229   continue


        iv1 = icount
        iv2 = icount+1
        if(iv2 == icount2+1) iv2 = 1

        da = poly2(iv1,1,1)
        db = poly2(iv1,2,1)
        dc = poly2(iv1,1,2)
        dd = poly2(iv1,2,2)
        de = poly2(iv2,1,2)
        df = poly2(iv2,2,2)

      area = 0.5d0 * abs(-db*dc + da * dd + db * de - dd * de - da * df
     .                      + dc*df)
        if(area < 1.e-10) then ! merge edges 1 and 2 and remove edge 2
c        write(*,*) "collinear",iv1,iv2
            if(iv1 .ne. icount2) then
                poly2(iv1,1,2) = poly2(iv2,1,2)
                poly2(iv1,2,2) = poly2(iv2,2,2)
                ! shift everything below upward
                do 33 ii = iv2,icount2-1
                    poly2(ii,1,1) = poly2(ii+1,1,1)
                    poly2(ii,2,1) = poly2(ii+1,2,1)
                    poly2(ii,1,2) = poly2(ii+1,1,2)
                    poly2(ii,2,2) = poly2(ii+1,2,2)
  33            continue
            else
                poly2(iv2,1,1) = poly2(iv1,1,1)
                poly2(iv2,2,1) = poly2(iv1,2,1)
            endif
            poly2(icount2,1,1) = -11
            poly2(icount2,2,1) = -11
            poly2(icount2,1,2) = -11
            poly2(icount2,2,2) = -11
            icount2 = icount2 - 1
            goto 999
        endif
  31  continue

c      poly3 = -11
c      irun = 1
c      do 32 icount = 1,icount2
c        if(idup(icount) .eq. 0) then
c            poly3(irun,1,1) = poly2(icount,1,1)
c            poly3(irun,2,1) = poly2(icount,2,1)
c            poly3(irun,1,2) = poly2(icount,1,2)
c            poly3(irun,2,2) = poly2(icount,2,2)
c            irun = irun + 1
c        endif
c  32  continue
c      poly2 = poly3
      end
