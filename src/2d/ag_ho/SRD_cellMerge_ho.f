c
c ---------------------------------------------------------------------
c
       subroutine SRD_cellMerge_ho(q,dlimit,nvar,irr,mitot,mjtot,qx,qy,
     .                       lstgrd,
     .                       dx,dy,lwidth,xlow,ylow,istage)
!
       implicit double precision (a-h, o-z)
       include "cirr.i"
       include "./quadrature.i"
       common /order2/ ssw, quad, nolimiter
       logical    quad, nolimiter
       dimension q(mitot,mjtot,nvar)
       dimension irr(mitot,mjtot)
       dimension gradmx(irrsize,nvar), gradmy(irrsize,nvar)
       dimension gradmxx(irrsize,nvar), gradmxy(irrsize,nvar)
       dimension gradmyy(irrsize,nvar)
       dimension gradmxxx(irrsize,nvar), gradmxxy(irrsize,nvar)
       dimension gradmxyy(irrsize,nvar), gradmyyy(irrsize,nvar)
       dimension qMerge(mitot,mjtot, nvar), qm(nvar)

       dimension a(9,9), rhs(9,nvar), b(9), f(9), w(9,nvar)
       dimension G(9,9)
       dimension valnew(mitot,mjtot,nvar)
       common   /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                  ismp,gradThreshold

       logical IS_GHOST, verbose/.true./
       IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                  j .le. lwidth .or. j .gt. mjtot-lwidth)


      TOL = 1.d-10

       if (verbose) then
          totmass =  bigconck_ho(q,irr,mitot,mjtot,lwidth,nvar)
          write(*,911) totmass
 911      format(/,"         SRDmass before redistribution is ",e30.20)
       endif






      qMerge(:,:,:) = 0.d0
      qMerge(1:lwidth,:,:) =  q(1:lwidth,:,:)
      qMerge(mitot-lwidth+1:mitot,:,:) =  q(mitot-lwidth+1:mitot,:,:)
      qMerge(:,1:lwidth,:) =  q(:,1:lwidth,:)
      qMerge(:,mjtot-lwidth+1:mjtot,:) =  q(:,mjtot-lwidth+1:mjtot,:)


c       form qMerge vals, these are the same regardless of the order of accuracy!
      do 10 j = lwidth+1, mjtot-lwidth
      do 10 i = lwidth+1, mitot-lwidth
            k = irr(i,j)

            if (k .eq. -1 .or. k .eq. lstgrd) then
              qMerge(i,j,:) = q(i,j,:)
              go to 10
            endif


             do 27 ic = 1, ncount(k)
                icurr = iidx(k,ic)
                jcurr = jidx(k,ic)
                koff = irr(icurr,jcurr)
                qMerge(i,j,:) = qMerge(i,j,:) + ar_ho(koff)*
     .                        q(icurr,jcurr,:)/numHoods(icurr,jcurr)

 27          continue
             qMerge(i,j,:) = qMerge(i,j,:) / volMerge_ho(k)
 10     continue





         ! gradient of merge neighborhoods, initialized to 0. set using neighboring merged tiels
         gradmx = 0.d0
         gradmy = 0.d0

         gradmxx = 0.d0
         gradmxy = 0.d0
         gradmyy = 0.d0

         gradmxxx = 0.d0
         gradmxxy = 0.d0
         gradmxyy = 0.d0
         gradmyyy = 0.d0

        if (ssw .eq. 2) then
        do 120 j = lwidth+1, mjtot-lwidth
        do 120 i = lwidth+1, mitot-lwidth
             k = irr(i,j)
             if (k .eq. -1 .or. k .eq. lstgrd) go to 120
             ! solid/regular cell doesnt participate in anything
             ! can reconstruct high order poly here, but for now do linear reconstruction
             ! solving least squares problem with all neighboring merged vals. Should preprocess the matrix and save

             rhs = 0.d0 ! initialize for accumulation
             a = 0.d0
             x0 = xcentMerge_ho(k)
             y0 = ycentMerge_ho(k)


            do 334 ioff = -mioff(i,j), mioff(i,j)
            do 334 joff = -mjoff(i,j), mjoff(i,j)
                 if (IS_GHOST(i+ioff,j+joff)) go to 334
                 koff = irr(i+ioff,j+joff)
                 if (koff .eq. -1) goto 334

                 if(koff .eq. lstgrd) then
                     xcoff = xlow + (i+ioff-0.5d0)*dx
                     ycoff = ylow + (j+joff-0.5d0)*dy
                     deltax = xcoff - x0
                     deltay = ycoff - y0
                 else
                     deltax = xcentMerge_ho(koff) - x0
                     deltay = ycentMerge_ho(koff) - y0
                 endif

                 f(1) = deltax/dx
                 f(2) = deltay/dy

                 f(3) = -qmshifts_ho(1,k)
                 f(4) = -qmshifts_ho(2,k)
                 f(5) = -qmshifts_ho(3,k)

                 if(koff .eq. lstgrd) then ! set up data on whole cell
                    ncount(koff)  = 1
                    iidx(koff, 1) = i+ioff
                    jidx(koff, 1) = j+joff
                  volmerge_ho(koff) = dx * dy / numhoods(i+ioff, j+joff)
                 endif


                 cvm = volmerge_ho(koff) ! off vol merge
      do 297 ic = 1, ncount(koff) ! compute weighted inner product of monomials on this neighborhood
                    icurr = iidx(koff,ic)
                    jcurr = jidx(koff,ic)
                    kcurr = irr(icurr,jcurr)
                    nhc = numhoods(icurr, jcurr) ! num hoods current
                    if(kcurr .eq. lstgrd) then
                      itri = 2 + 1
                      poly(1,1,kcurr) = xlow + (dfloat(icurr)-1.d0)*dx
                      poly(1,2,kcurr) = ylow + (dfloat(jcurr)-1.d0)*dy

                      poly(2,1,kcurr) = xlow + (dfloat(icurr)-0.d0)*dx
                      poly(2,2,kcurr) = ylow + (dfloat(jcurr)-1.d0)*dy

                      poly(3,1,kcurr) = xlow + (dfloat(icurr)-0.d0)*dx
                      poly(3,2,kcurr) = ylow + (dfloat(jcurr)-0.d0)*dy

                      poly(4,1,kcurr) = xlow + (dfloat(icurr)-1.d0)*dx
                      poly(4,2,kcurr) = ylow + (dfloat(jcurr)-0.d0)*dy
                    else
                     ivert = 1
                      do 290 while (poly(ivert+1,1,kcurr) .ne. -11.)
                        ivert = ivert + 1
  290                   continue
                       itri = ivert - 3
                    endif

      ! computing the inner product on each triangle of neighborhood member
                  idx1 = 1
                  do 291 it = 1, itri-1 ! for each  triangle except the last one
                    idx2 = it + 1
                    idx3 = it + 2

                    x1 = poly(idx1,1,kcurr)
                    y1 = poly(idx1,2,kcurr)

                    x2 = poly(idx2,1,kcurr)
                    y2 = poly(idx2,2,kcurr)

                    x3 = poly(idx3,1,kcurr)
                    y3 = poly(idx3,2,kcurr)

                    artri = area(x1, x2, x3, y1, y2, y3)

                    do 292 itq = 1,ntriquad
                        xval = x1 * rtri(itq) + x2 * stri(itq)
     .                     +  x3 * (1.d0-rtri(itq)-stri(itq))
                        yval = y1 * rtri(itq) + y2 * stri(itq)
     .                     +  y3 * (1.d0-rtri(itq)-stri(itq))

            f(3) = f(3) + (artri/nhc/cvm) * wtri(itq) *( (xval-x0)**2 )
     .                                                / (dx**2)

            f(4) = f(4) + (artri/nhc/cvm) * wtri(itq) *( (xval-x0)/dx )
     .                                                *( (yval-y0)/dy )

            f(5) = f(5) + (artri/nhc/cvm) * wtri(itq) *( (yval-y0)**2 )
     .                                                / (dy**2)
 292        continue ! for each quadrature point on each triangle
 291        continue ! for each triangle

            if(kcurr .ne. lstgrd) then
      nhc = numhoods(icurr, jcurr)
      cvm = volmerge_ho(koff)
      x1 = poly(ivert-2,1,kcurr)
      y1 = poly(ivert-2,2,kcurr)
      x2 = bdry(1,1,kcurr)
      y2 = bdry(1,2,kcurr)
      x3 = bdry(3,1,kcurr)
      y3 = bdry(3,2,kcurr)
      x4 = bdry(2,1,kcurr)
      y4 = bdry(2,2,kcurr)

      do nq = 1,ntriquad_ho
      artri = dJ2(x1,x2,x3,x4,y1,y2,y3,y4,nq)
      xval = rs2xy_2(x1,x2,x3,x4,nq)
      yval = rs2xy_2(y1,y2,y3,y4,nq)
      f(3) = f(3) + (artri/nhc/cvm) * wtri_ho(nq) *( (xval-x0)**2 )
     .                                                / (dx**2)

      f(4) = f(4) + (artri/nhc/cvm) * wtri_ho(nq) *( (xval-x0)/dx )
     .                                                *( (yval-y0)/dy )

      f(5) = f(5) + (artri/nhc/cvm) * wtri_ho(nq) *( (yval-y0)**2 )
     .                                                / (dy**2)
      enddo
            endif

 297        continue ! for each subelement of the current neighborhood
        ! integration is complete, let's accumulate

            do ii = 1,5
            do jj = 1,5
                a(ii,jj) = a(ii,jj) + f(ii)*f(jj)
            enddo

                rhs(ii,:) = rhs(ii,:)
     .               + f(ii) * (qMerge(i+ioff,j+joff,:) - qMerge(i,j,:))
            enddo

 2334         continue


 334           continue



      ! SOLVE THE NORMAL EQUATIONS.
             call cholesky(9, 5, a, G)
             do mm = 1, nvar
                call solve(9,5,G,rhs(:,mm))

                gradmyy(k,mm) =  2.d0*rhs(5,mm)/(dy**2)
                gradmxy(k,mm) =  rhs(4,mm)/(dx*dy)
                gradmxx(k,mm) =  2.d0*rhs(3,mm)/(dx**2)
                gradmy(k,mm)  =  rhs(2,mm)/dy
                gradmx(k,mm)  =  rhs(1,mm)/dx
            end do




 120    continue ! iterate over merging tiles on entire grid

        elseif(ssw .eq. 3) then


        do 620 j = lwidth+1, mjtot-lwidth
        do 620 i = lwidth+1, mitot-lwidth
             k = irr(i,j)
             if (k .eq. -1 .or. k .eq. lstgrd) go to 620
             ! solid/regular cell doesnt participate in anything
             ! can reconstruct high order poly here, but for now do linear reconstruction
             ! solving least squares problem with all neighboring merged vals. Should preprocess the matrix and save

             rhs = 0.d0 ! initialize for accumulation
             a = 0.d0
             x0 = xcentMerge_ho(k)
             y0 = ycentMerge_ho(k)


            do 634 ioff = -mioff(i,j), mioff(i,j)
            do 634 joff = -mjoff(i,j), mjoff(i,j)
                 if (IS_GHOST(i+ioff,j+joff)) go to 634
                 koff = irr(i+ioff,j+joff)
                 if (koff .eq. -1) goto 634

                 if(koff .eq. lstgrd) then
                     xcoff = xlow + (i+ioff-0.5d0)*dx
                     ycoff = ylow + (j+joff-0.5d0)*dy
                     deltax = xcoff - x0
                     deltay = ycoff - y0
                 else
                     deltax = xcentMerge_ho(koff) - x0
                     deltay = ycentMerge_ho(koff) - y0
                 endif

                 f(1) = deltax/dx
                 f(2) = deltay/dy
                 f(3) = -qmshifts_ho(1,k)
                 f(4) = -qmshifts_ho(2,k)
                 f(5) = -qmshifts_ho(3,k)
                 f(6) = -qmshifts_ho(4,k)
                 f(7) = -qmshifts_ho(5,k)
                 f(8) = -qmshifts_ho(6,k)
                 f(9) = -qmshifts_ho(7,k)

                 if(koff .eq. lstgrd) then ! set up data on whole cell
                    ncount(koff)  = 1
                    iidx(koff, 1) = i+ioff
                    jidx(koff, 1) = j+joff
                  volmerge_ho(koff) = dx * dy / numhoods(i+ioff, j+joff)
                 endif


                 cvm = volmerge_ho(koff) ! off vol merge
      do 697 ic = 1, ncount(koff) ! compute weighted inner product of monomials on this neighborhood
                    icurr = iidx(koff,ic)
                    jcurr = jidx(koff,ic)
                    kcurr = irr(icurr,jcurr)
                    nhc = numhoods(icurr, jcurr) ! num hoods current
                    if(kcurr .eq. lstgrd) then
                      itri = 2 + 1
                      poly(1,1,kcurr) = xlow + (dfloat(icurr)-1.d0)*dx
                      poly(1,2,kcurr) = ylow + (dfloat(jcurr)-1.d0)*dy

                      poly(2,1,kcurr) = xlow + (dfloat(icurr)-0.d0)*dx
                      poly(2,2,kcurr) = ylow + (dfloat(jcurr)-1.d0)*dy

                      poly(3,1,kcurr) = xlow + (dfloat(icurr)-0.d0)*dx
                      poly(3,2,kcurr) = ylow + (dfloat(jcurr)-0.d0)*dy

                      poly(4,1,kcurr) = xlow + (dfloat(icurr)-1.d0)*dx
                      poly(4,2,kcurr) = ylow + (dfloat(jcurr)-0.d0)*dy
                    else
                     ivert = 1
                      do 690 while (poly(ivert+1,1,kcurr) .ne. -11.)
                        ivert = ivert + 1
  690                   continue
                       itri = ivert - 3
                    endif

      ! computing the inner product on each triangle of neighborhood member
                  idx1 = 1
                  do 691 it = 1, itri - 1! for each  triangle except the last one
                    idx2 = it + 1
                    idx3 = it + 2

                    x1 = poly(idx1,1,kcurr)
                    y1 = poly(idx1,2,kcurr)

                    x2 = poly(idx2,1,kcurr)
                    y2 = poly(idx2,2,kcurr)

                    x3 = poly(idx3,1,kcurr)
                    y3 = poly(idx3,2,kcurr)

                    artri = area(x1, x2, x3, y1, y2, y3)

                    do 692 itq = 1,ntriquad
                        xval = x1 * rtri(itq) + x2 * stri(itq)
     .                     +  x3 * (1.d0-rtri(itq)-stri(itq))
                        yval = y1 * rtri(itq) + y2 * stri(itq)
     .                     +  y3 * (1.d0-rtri(itq)-stri(itq))

            f(3) = f(3) + (artri/nhc/cvm) * wtri(itq) *( (xval-x0)**2 )
     .                                                / (dx**2)

            f(4) = f(4) + (artri/nhc/cvm) * wtri(itq) *( (xval-x0)/dx )
     .                                                *( (yval-y0)/dy )

            f(5) = f(5) + (artri/nhc/cvm) * wtri(itq) *( (yval-y0)**2 )
     .                                                / (dy**2)





            f(6) = f(6) + (artri/nhc/cvm) * wtri(itq) *( (xval-x0)**3 )
     .                                                / (dx**3)

            f(7) = f(7) + (artri/nhc/cvm) * wtri(itq)
     .                   *( (xval-x0)**2 / dx**2)*( (yval-y0)/dy )

            f(8) = f(8) + (artri/nhc/cvm) * wtri(itq)
     .                   *( (xval-x0) / dx)*( (yval-y0)**2/dy**2 )

            f(9) = f(9) + (artri/nhc/cvm) * wtri(itq) *( (yval-y0)**3 )
     .                                                / (dy**3)




 692        continue ! for each quadrature point on each triangle
 691        continue ! for each triangle

             if(kcurr .ne. lstgrd) then
      nhc = numhoods(icurr, jcurr)
      cvm = volmerge_ho(koff)
      x1 = poly(ivert-2,1,kcurr)
      y1 = poly(ivert-2,2,kcurr)

      x2 = bdry(1,1,kcurr)
      y2 = bdry(1,2,kcurr)

      x3 = bdry(4,1,kcurr)
      y3 = bdry(4,2,kcurr)

      x4 = bdry(2,1,kcurr)
      y4 = bdry(2,2,kcurr)

      x5 = bdry(3,1,kcurr)
      y5 = bdry(3,2,kcurr)

      do nq = 1,ntriquad_ho
      artri = dJ3(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,nq)

      xval = rs2xy_3(x1,x2,x3,x4,x5,nq)
      yval = rs2xy_3(y1,y2,y3,y4,y5,nq)

      f(3) = f(3) +
     .       + (artri/nhc/cvm)*wtri_ho(nq)*(xval-x0)**2 /(dx**2)
      f(4) = f(4) +
     .       + (artri/nhc/cvm)*wtri_ho(nq)*(xval-x0)*(yval-y0)/(dx * dy)
      f(5) = f(5) +
     .       + (artri/nhc/cvm)*wtri_ho(nq)*(yval-y0)**2 /(dy**2)


      f(6) = f(6) +
     .  + (artri/nhc/cvm)*wtri_ho(nq)*(xval-x0)**3 /(dx**3)
      f(7) = f(7) +
     .  + (artri/nhc/cvm)*wtri_ho(nq)*(yval-y0)*(xval-x0)**2 /(dy*dx**2)
      f(8) = f(8) +
     .  + (artri/nhc/cvm)*wtri_ho(nq)*(xval-x0)*(yval-y0)**2 /(dx*dy**2)
      f(9) = f(9) +
     .  + (artri/nhc/cvm)*wtri_ho(nq)*(yval-y0)**3 /(dy**3)

      end do
            endif



 697        continue ! for each subelement of the current neighborhood
        ! integration is complete, let's accumulate

            do ii = 1,9
            do jj = 1,9
                a(ii,jj) = a(ii,jj) + f(ii)*f(jj)
            enddo

                rhs(ii,:) = rhs(ii,:)
     .               + f(ii) * (qMerge(i+ioff,j+joff,:) - qMerge(i,j,:))
            enddo

 6334         continue


 634           continue



      ! SOLVE THE NORMAL EQUATIONS.
             call cholesky(9, 9, a, G)
             do mm = 1, nvar
                call solve(9,9,G,rhs(:,mm))
                gradmyyy(k,mm) = 6.d0*rhs(9,mm)/(dy**3)
                gradmxyy(k,mm) = 2.d0*rhs(8,mm)/(dx*dy**2)
                gradmxxy(k,mm) = 2.d0*rhs(7,mm)/(dy*dx**2)
                gradmxxx(k,mm) = 6.d0*rhs(6,mm)/(dx**3)
                gradmyy(k,mm) =  2.d0*rhs(5,mm)/(dy**2)
                gradmxy(k,mm) =  rhs(4,mm)/(dx*dy)
                gradmxx(k,mm) =  2.d0*rhs(3,mm)/(dx**2)
                gradmy(k,mm)  =  rhs(2,mm)/dy
                gradmx(k,mm)  =  rhs(1,mm)/dx
            end do


 620    continue ! iterate over merging tiles on entire grid












        endif


















      if(nolimiter .eqv. .false.) then
        call limiter_srd(qMerge,gradmx,gradmy,mitot,mjtot,irr, nvar,
     .                   dx, dy,lwidth, lstgrd, xlow, ylow)
      endif















      valnew = 0.d0  !  all cells initialized to 0

      if(ssw .eq. 2) then


      do 510 j = lwidth+1, mjtot-lwidth
      do 510 i = lwidth+1, mitot-lwidth
          k = irr(i,j)
          if (k .eq. -1) then
             valnew(i,j,:) = q(i,j,:)
             go to 510  ! does not contribute
          endif
          if(k .eq. lstgrd) then
           valnew(i,j,:) = valnew(i,j,:) + qMerge(i,j,:) / numHoods(i,j)
           cycle
          endif

      ! if here, I am on a cut cell

      do 410 ic = 1, ncount(k)
                 ioff = iidx(k,ic)
                 joff = jidx(k,ic)
                 koff = irr(ioff,joff)
                 call getCellCentroid_ho(lstgrd,ioff,joff,xc,yc,xlow,
     .                            ylow,dx,dy,koff)
      qm = 0.d0
      if(koff .eq. lstgrd) then
          arr = dx*dy
          itri = 2 + 1
          poly(1,1,koff) = xlow + (dfloat(ioff)-1.d0)*dx
          poly(1,2,koff) = ylow + (dfloat(joff)-1.d0)*dy

          poly(2,1,koff) = xlow + (dfloat(ioff)-0.d0)*dx
          poly(2,2,koff) = ylow + (dfloat(joff)-1.d0)*dy

          poly(3,1,koff) = xlow + (dfloat(ioff)-0.d0)*dx
          poly(3,2,koff) = ylow + (dfloat(joff)-0.d0)*dy

          poly(4,1,koff) = xlow + (dfloat(ioff)-1.d0)*dx
          poly(4,2,koff) = ylow + (dfloat(joff)-0.d0)*dy
      else
          arr = ar_ho(koff)
          ivert = 1
          do 240 while (poly(ivert+1,1,koff) .ne. -11.)
            ivert = ivert + 1
  240      continue
           itri = ivert - 3
      endif

      idx1 = 1
      do 21 it = 1, itri -1 ! for each  triangle except the last one
            idx2 = it + 1
            idx3 = it + 2

            x1 = poly(idx1,1,koff)
            y1 = poly(idx1,2,koff)

            x2 = poly(idx2,1,koff)
            y2 = poly(idx2,2,koff)

            x3 = poly(idx3,1,koff)
            y3 = poly(idx3,2,koff)

            artri = area(x1, x2, x3, y1, y2, y3)

            do 225 itq = 1,ntriquad
                xval = x1 * rtri(itq) + x2 * stri(itq)
     .              +  x3 * (1.d0-rtri(itq)-stri(itq))
                yval = y1 * rtri(itq) + y2 * stri(itq)
     .              +  y3 * (1.d0-rtri(itq)-stri(itq))

                diffx = xval-xcentMerge_ho(k)
                diffy = yval-ycentMerge_ho(k)

            qm(:) = qm(:) + (artri/arr) * wtri(itq) *
     .  (qMerge(i,j,:)+ gradmx(k,:) * diffx + gradmy(k,:) * diffy
     . + 0.5d0*gradmxx(k,:)*(diffx*diffx - (dx**2)*qmshifts_ho(1,k))
     . +       gradmxy(k,:)*(diffx*diffy - dx*dy*qmshifts_ho(2,k))
     . + 0.5d0*gradmyy(k,:)*(diffy*diffy - (dy**2)*qmshifts_ho(3,k))  )
  225   continue ! for each quadrature point on each triangle
  21   continue ! for each triangle
            if(koff .ne. lstgrd) then
      x1 = poly(ivert-2,1,koff)
      y1 = poly(ivert-2,2,koff)
      x2 = bdry(1,1,koff)
      y2 = bdry(1,2,koff)
      x3 = bdry(3,1,koff)
      y3 = bdry(3,2,koff)
      x4 = bdry(2,1,koff)
      y4 = bdry(2,2,koff)

      do nq = 1,ntriquad_ho
      artri = dJ2(x1,x2,x3,x4,y1,y2,y3,y4,nq)
      xval = rs2xy_2(x1,x2,x3,x4,nq)
      yval = rs2xy_2(y1,y2,y3,y4,nq)

                diffx = xval-xcentMerge_ho(k)
                diffy = yval-ycentMerge_ho(k)

            qm(:) = qm(:) + (artri/arr) * wtri_ho(nq) *
     .  (qMerge(i,j,:)+ gradmx(k,:) * diffx + gradmy(k,:) * diffy
     . + 0.5d0*gradmxx(k,:)*(diffx*diffx - (dx**2)*qmshifts_ho(1,k))
     . +       gradmxy(k,:)*(diffx*diffy -   dx*dy*qmshifts_ho(2,k))
     . + 0.5d0*gradmyy(k,:)*(diffy*diffy - (dy**2)*qmshifts_ho(3,k))  )
      enddo
            endif


        valnew(ioff,joff,:) = valnew(ioff,joff,:) +
     .              qm(:)/numHoods(ioff,joff)

 410       continue ! cycle through all the subcells



 510    continue ! onto the next cell in the grid







      elseif(ssw .eq. 3) then


      do 1510 j = lwidth+1, mjtot-lwidth
      do 1510 i = lwidth+1, mitot-lwidth
          k = irr(i,j)
          if (k .eq. -1) then
             valnew(i,j,:) = q(i,j,:)
             go to 1510  ! does not contribute
          endif
          if(k .eq. lstgrd) then
           valnew(i,j,:) = valnew(i,j,:) + qMerge(i,j,:) / numHoods(i,j)
           cycle
          endif

      ! if here, I am on a cut cell

      do 1410 ic = 1, ncount(k)
                 ioff = iidx(k,ic)
                 joff = jidx(k,ic)
                 koff = irr(ioff,joff)
                 call getCellCentroid_ho(lstgrd,ioff,joff,xc,yc,xlow,
     .                            ylow,dx,dy,koff)
      qm = 0.d0
      if(koff .eq. lstgrd) then
          arr = dx*dy
          itri = 2 + 1
          poly(1,1,koff) = xlow + (dfloat(ioff)-1.d0)*dx
          poly(1,2,koff) = ylow + (dfloat(joff)-1.d0)*dy

          poly(2,1,koff) = xlow + (dfloat(ioff)-0.d0)*dx
          poly(2,2,koff) = ylow + (dfloat(joff)-1.d0)*dy

          poly(3,1,koff) = xlow + (dfloat(ioff)-0.d0)*dx
          poly(3,2,koff) = ylow + (dfloat(joff)-0.d0)*dy

          poly(4,1,koff) = xlow + (dfloat(ioff)-1.d0)*dx
          poly(4,2,koff) = ylow + (dfloat(joff)-0.d0)*dy
      else
          arr = ar_ho(koff)
          ivert = 1
          do 1240 while (poly(ivert+1,1,koff) .ne. -11.)
            ivert = ivert + 1
1240      continue
           itri = ivert - 3
      endif

      idx1 = 1
      do 121 it = 1, itri - 1 ! for each  triangle except the last one
            idx2 = it + 1
            idx3 = it + 2

            x1 = poly(idx1,1,koff)
            y1 = poly(idx1,2,koff)

            x2 = poly(idx2,1,koff)
            y2 = poly(idx2,2,koff)

            x3 = poly(idx3,1,koff)
            y3 = poly(idx3,2,koff)

            artri = area(x1, x2, x3, y1, y2, y3)

            do 1225 itq = 1,ntriquad
                xval = x1 * rtri(itq) + x2 * stri(itq)
     .              +  x3 * (1.d0-rtri(itq)-stri(itq))
                yval = y1 * rtri(itq) + y2 * stri(itq)
     .              +  y3 * (1.d0-rtri(itq)-stri(itq))

                diffx = xval-xcentMerge_ho(k)
                diffy = yval-ycentMerge_ho(k)

            qm(:) = qm(:) + (artri/arr) * wtri(itq) *
     .  (qMerge(i,j,:)+ gradmx(k,:) * diffx + gradmy(k,:) * diffy
     . + 0.5d0*gradmxx(k,:)*(diffx*diffx - (dx**2)*qmshifts_ho(1,k))
     . +       gradmxy(k,:)*(diffx*diffy - dx*dy*qmshifts_ho(2,k))
     . + 0.5d0*gradmyy(k,:)*(diffy*diffy - (dy**2)*qmshifts_ho(3,k))
     .+(1.d0/6.d0)*gradmxxx(k,:)*(diffx**3 - (dx**3)*qmshifts_ho(4,k))
     .+ 0.5d0*gradmxxy(k,:)*(diffy*diffx**2-dx**2*dy*qmshifts_ho(5,k))
     .+ 0.5d0*gradmxyy(k,:)*(diffy**2*diffx-dx*dy**2*qmshifts_ho(6,k))
     .+(1.d0/6.d0)*gradmyyy(k,:)*(diffy**3 - (dy**3)*qmshifts_ho(7,k))
     .)
 1225   continue ! for each quadrature point on each triangle
 121   continue ! for each triangle
             if(koff .ne. lstgrd) then
      x1 = poly(ivert-2,1,koff)
      y1 = poly(ivert-2,2,koff)

      x2 = bdry(1,1,koff)
      y2 = bdry(1,2,koff)

      x3 = bdry(4,1,koff)
      y3 = bdry(4,2,koff)

      x4 = bdry(2,1,koff)
      y4 = bdry(2,2,koff)

      x5 = bdry(3,1,koff)
      y5 = bdry(3,2,koff)

      do nq = 1,ntriquad_ho
      artri = dJ3(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,nq)

      xval = rs2xy_3(x1,x2,x3,x4,x5,nq)
      yval = rs2xy_3(y1,y2,y3,y4,y5,nq)


                diffx = xval-xcentMerge_ho(k)
                diffy = yval-ycentMerge_ho(k)

            qm(:) = qm(:) + (artri/arr) * wtri_ho(nq) *
     .  (qMerge(i,j,:)+ gradmx(k,:) * diffx + gradmy(k,:) * diffy
     . + 0.5d0*gradmxx(k,:)*(diffx*diffx - (dx**2)*qmshifts_ho(1,k))
     . +       gradmxy(k,:)*(diffx*diffy - dx*dy*qmshifts_ho(2,k))
     . + 0.5d0*gradmyy(k,:)*(diffy*diffy - (dy**2)*qmshifts_ho(3,k))
     .+(1.d0/6.d0)*gradmxxx(k,:)*(diffx**3 - (dx**3)*qmshifts_ho(4,k))
     .+ 0.5d0*gradmxxy(k,:)*(diffy*diffx**2-dx**2*dy*qmshifts_ho(5,k))
     .+ 0.5d0*gradmxyy(k,:)*(diffy**2*diffx-dx*dy**2*qmshifts_ho(6,k))
     .+(1.d0/6.d0)*gradmyyy(k,:)*(diffy**3 - (dy**3)*qmshifts_ho(7,k))
     .)
      enddo
            endif



        valnew(ioff,joff,:) = valnew(ioff,joff,:) +
     .              qm(:)/numHoods(ioff,joff)



 1410       continue ! cycle through all the subcells



 1510    continue ! onto the next cell in the grid






      endif

      valnew(1:lwidth,:,:) =  q(1:lwidth,:,:)
      valnew(mitot-lwidth+1:mitot,:,:) =  q(mitot-lwidth+1:mitot,:,:)
      valnew(:,1:lwidth,:) =  q(:,1:lwidth,:)
      valnew(:,mjtot-lwidth+1:mjtot,:) =  q(:,mjtot-lwidth+1:mjtot,:)

      q = valnew



       if (verbose) then
          totmassafter =  bigconck_ho(q,irr,mitot,mjtot,lwidth,nvar)
          write(*,921) totmassafter
 921      format(/,"         SRDmass after redistribution is  ",e30.20)
          print *, "diff is ",totmassafter-totmass
       endif







      end
!
!
!
!
!
!c
!c -------------------------------------------------------------------------
!c
      subroutine getCellCentroid_ho(lstgrd,i,j,xc,yc,xlow,ylow,dx,dy,k)

      implicit double precision (a-h,o-z)
      include "cirr.i"

      if (k .eq. lstgrd) then
         xc = xlow + (i-0.5d0)*dx
         yc = ylow + (j-0.5d0)*dy
      else if (k .eq. -1) then
         xc = 0.d0
         yc = 0.d0
      else
         xc = xcirr_ho(k)
         yc = ycirr_ho(k)
      endif

      return
      end




      function bigconck_ho(q,irr,mitot,mjtot,lwidth,nvar)
c
      implicit double precision (a-h,o-z)
      include "cirr.i"

      dimension q(mitot,mjtot,nvar), irr(mitot,mjtot)

      totmass = 0.d0
      ar_ho(-1)  = 0.d0


      do 10 iv = 1, nvar
      do 10 j = lwidth+1, mjtot-lwidth
      do 10 i = lwidth+1, mitot-lwidth

         k = irr(i,j)
         totmass = totmass + q(i,j,iv)*ar_ho(k)

 10   continue

      bigconck_ho = totmass

      return
      end

