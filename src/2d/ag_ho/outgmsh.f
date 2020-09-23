c
c -------------------------------------------------------------
c
      subroutine outgmsh(qin,nvar,mptr,irr,
     1                  mitot,mjtot,qx,qy,lwidth,lstgrd,
     2                  dx,dy,xlow,ylow,time)
c
      implicit double precision (a-h,o-z)
      include "cirr.i"
      dimension irr(mitot,mjtot), q(mitot,mjtot,nvar)
      dimension qin(mitot,mjtot,nvar)
      dimension qx(mitot,mjtot,nvar),qy(mitot,mjtot,nvar)
      dimension qxx(mitot,mjtot,nvar),qyy(mitot,mjtot,nvar)
      dimension qxy(mitot,mjtot,nvar)
      dimension qxxx(mitot,mjtot,nvar),qxyy(mitot,mjtot,nvar)
      dimension qxxy(mitot,mjtot,nvar),qyyy(mitot,mjtot,nvar)
      dimension valprim(4)
      dimension exactsoln(1)
      dimension exactval(mitot,mjtot,nvar)
      dimension dlimit(mitot,mjtot)
      dimension recon(4,nvar), corners(4,2)
      dimension numvert(irrsize)
      dimension dnumsol(nvar)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold,ilts,ibc,imodal
      common /order2/ ssw, quad, nolimiter
      logical  eout, quad, nolimiter
      logical pwconst

      data      spval/-1.e10/, eout/.true./
c      character*23  filename
      common /dump/ uvel,vvel
      dimension uvel(849,848),vvel(848,849)


      q = qin

!      open(231, file = 'out.pos')
      ibuff = 4

      xlowb = xlow - lwidth*dx
      ylowb = ylow - lwidth*dy

      qx = 0.d0
      qy = 0.d0
      qxx = 0.d0
      qxy = 0.d0
      qyy = 0.d0
      qxxx = 0.d0
      qxxy = 0.d0
      qxyy = 0.d0
      qyyy = 0.d0
!
!      call vctoprm(q,q,mitot,mjtot)
c
c     ### call for exterior bcs at each stage so can use slopes
c    ## NOTE THAT BNDRY CELLS FROM OTHER GRIDS NOT SET
            xhigh = xlowb + mitot*dx
            yhigh = ylowb + mjtot*dy
            call pphysbdlin(xlowb,xhigh,ylowb,yhigh,level,mitot,mjtot,
     &                   nvar,q,time,dx,dy,qx,qy,irr,lstgrd)





      if(quad .and. (nolimiter .eqv. .true.) ) then ! yes in primitive variables, and without limiter (this is an accurate projection between the two
      call slopes(q,qx,qy,qxx,qxy,qyy,
     &               qxxx, qxxy, qxyy, qyyy,
     &               mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &               xlowb,ylowb,nvar)

      call vctoprm(mitot,mjtot,nvar,
     . dx, dy, lstgrd, xlowb, ylowb, irr,
     . q,q,qx,qy,qxx, qxy, qyy, qxxx, qxxy, qxyy, qyyy)

      call slopes(q,qx,qy,qxx,qxy,qyy,
     &               qxxx, qxxy, qxyy, qyyy,
     &               mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &               xlowb,ylowb,nvar)




      elseif(quad .and. (nolimiter .eqv. .false.) ) then ! yes in primitive variables, and with limiter
                                                      ! only works for p = 1
      call slopes(q,qx,qy,qxx,qxy,qyy,
     &               qxxx, qxxy, qxyy, qyyy,
     &               mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &               xlowb,ylowb,nvar)

      call vctoprm(mitot,mjtot,nvar,
     . dx, dy, lstgrd, xlowb, ylowb, irr,
     . q,q,qx,qy,qxx, qxy, qyy, qxxx, qxxy, qxyy, qyyy)


      call slopes(q,qx,qy,qxx,qxy,qyy,
     &               qxxx, qxxy, qxyy, qyyy,
     &               mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &               xlowb,ylowb,nvar)

      ! limiter is only applied to the primitive variables.
!      if(imodal .eq. 0) then
      call limiter(q,  qx, qy,  mitot, mjtot, irr, nvar,dx, dy,
     .               lwidth, lstgrd, xlowb, ylowb)
!      else
!      call limiter_modal(q,  qx, qy,  qxx,qxy,qyy,
!     .             qxxx, qxxy, qxyy, qyyy,
!     .               mitot, mjtot, irr, nvar,dx, dy,
!     .               lwidth, lstgrd, xlowb, ylowb)
!      endif

      elseif (ssw .ne. 0.d0) then   ! conserved variables
         call slopes(q,qx,qy,qxx,qxy,qyy,
     &               qxxx, qxxy, qxyy, qyyy,
     &               mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &               xlowb,ylowb,nvar)
      endif

      ivar = 1
      write(231,*) 'View "U0 - T = ',time,'" {'
      do 15 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 15 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 15
         if (kirr .eq. lstgrd ) then
            xcen = xlowb + (dfloat(i)-.5d0)* dx
            ycen = ylowb + (dfloat(j)-.5d0)* dy

             xcorner = xlowb + (dfloat(i)-1.)* dx
             ycorner = ylowb + (dfloat(j)-1.)* dy
             do itimes = 1,4
                if (itimes .eq. 1) then ! get all 4 corners of mesh, pw constant sol
                   xc = xcorner
                   yc = ycorner
                else if (itimes .eq. 2) then
                   xc = xcorner
                   yc = ycorner + dy
                else if (itimes .eq. 3) then
                   xc = xcorner + dx
                   yc = ycorner + dy
                else if (itimes .eq. 4) then
                   xc = xcorner + dx
                   yc = ycorner
                endif

            corners(itimes,1) = xc
            corners(itimes,2) = yc


       deltax = xc-xcen
       deltay = yc-ycen

       do mmm = 1,nvar
       call evalU(dnumsol,deltax, deltay, dx, dy,
     .            q, qx, qy,
     .            qxx, qxy, qyy,
     .            qxxx, qxxy, qxyy, qyyy,
     .            i,j, kirr, mitot, mjtot, nvar, mmm)
       end do
       recon(:,itimes) = dnumsol(:)

!            recon(:,itimes) = q(i,j,:)
!     .   +(xc-xcen)*qx(i,j,:)
!     .   +(yc-ycen)*qy(i,j,:)
!     .   + 0.5d0*qxx(i,j,:)*( (xc-xcen)**2-(dx**2)/12.d0)
!     .   + (xc-xcen)*(yc-ycen)*qxy(i,j,:)
!     .   + 0.5d0*qyy(i,j,:)*( (yc-ycen)**2-(dy**2)/12.d0)
              end do

       write(231,*)  'SQ ( ', corners(1,1),",",
     .                        corners(1,2),",0,",
     .                        corners(2,1),",",
     .                        corners(2,2),",0,",
     .                        corners(3,1),",",
     .                        corners(3,2),",0,",
     .                        corners(4,1),",",
     .                        corners(4,2),",0)",
     .                        "{",
     .                        recon(ivar,1),",",
     .                        recon(ivar,2),",",
     .                        recon(ivar,3),",",
     .                        recon(ivar,4),
     .                        "};"


         else
            ivert = 1
            do 290 while (poly(ivert+1,1,kirr) .ne. -11.)
               ivert = ivert + 1
  290       continue
             if(ihob .eq. 0 ) then
                itri = ivert - 3
             else
                itri = ivert - 4 ! don't do the last one
             endif
                  idx1 = 1
                  do 291 it = 1, itri ! for each  triangle
                    idx2 = it + 1
                    idx3 = it + 2

                    x1 = poly(idx1,1,kirr)
                    y1 = poly(idx1,2,kirr)

                    x2 = poly(idx2,1,kirr)
                    y2 = poly(idx2,2,kirr)

                    x3 = poly(idx3,1,kirr)
                    y3 = poly(idx3,2,kirr)


                xcen = xcirr(kirr)
                ycen = ycirr(kirr)

            xc = x1
            yc = y1
            recon(:,1) = q(i,j,:)
     .   +(xc-xcen)*qx(i,j,:)
     .   +(yc-ycen)*qy(i,j,:)
     .   + 0.5d0*qxx(i,j,:)*( (xc-xcen)**2-(dx**2)/12.d0)
     .   + (xc-xcen)*(yc-ycen)*qxy(i,j,:)
     .   + 0.5d0*qyy(i,j,:)*( (yc-ycen)**2-(dy**2)/12.d0)

            xc = x2
            yc = y2
            recon(:,2) = q(i,j,:)
     .   +(xc-xcen)*qx(i,j,:)
     .   +(yc-ycen)*qy(i,j,:)
     .   + 0.5d0*qxx(i,j,:)*( (xc-xcen)**2-(dx**2)/12.d0)
     .   + (xc-xcen)*(yc-ycen)*qxy(i,j,:)
     .   + 0.5d0*qyy(i,j,:)*( (yc-ycen)**2-(dy**2)/12.d0)

            xc = x3
            yc = y3
            recon(:,3) = q(i,j,:)
     .   +(xc-xcen)*qx(i,j,:)
     .   +(yc-ycen)*qy(i,j,:)
     .   + 0.5d0*qxx(i,j,:)*( (xc-xcen)**2-(dx**2)/12.d0)
     .   + (xc-xcen)*(yc-ycen)*qxy(i,j,:)
     .   + 0.5d0*qyy(i,j,:)*( (yc-ycen)**2-(dy**2)/12.d0)

       write(231,*)  'ST ( ', x1,",",
     .                        y1,",0,",
     .                        x2,",",
     .                        y2,",0,",
     .                        x3,",",
     .                        y3,",0)",
     .                        "{",
     .                        recon(ivar,1),",",
     .                        recon(ivar,2),",",
     .                        recon(ivar,3),
     .                        "};"
291   continue






                  idx1 = ivert-2
                  x1 = poly(idx1,1,kirr)
                  y1 = poly(idx1,2,kirr)
                  do 11 it = 1, nbdry-1 ! for each  triangle
                    idx2 = it + 0
                    idx3 = it + 1

                    x2 = bdry(idx2,1,kirr)
                    y2 = bdry(idx2,2,kirr)

                    x3 = bdry(idx3,1,kirr)
                    y3 = bdry(idx3,2,kirr)


                xcen = xcirr(kirr)
                ycen = ycirr(kirr)

            xc = x1
            yc = y1
            recon(:,1) = q(i,j,:)
     .   +(xc-xcen)*qx(i,j,:)
     .   +(yc-ycen)*qy(i,j,:)
     .   + 0.5d0*qxx(i,j,:)*( (xc-xcen)**2-(dx**2)/12.d0)
     .   + (xc-xcen)*(yc-ycen)*qxy(i,j,:)
     .   + 0.5d0*qyy(i,j,:)*( (yc-ycen)**2-(dy**2)/12.d0)

            xc = x2
            yc = y2
            recon(:,2) = q(i,j,:)
     .   +(xc-xcen)*qx(i,j,:)
     .   +(yc-ycen)*qy(i,j,:)
     .   + 0.5d0*qxx(i,j,:)*( (xc-xcen)**2-(dx**2)/12.d0)
     .   + (xc-xcen)*(yc-ycen)*qxy(i,j,:)
     .   + 0.5d0*qyy(i,j,:)*( (yc-ycen)**2-(dy**2)/12.d0)

            xc = x3
            yc = y3
            recon(:,3) = q(i,j,:)
     .   +(xc-xcen)*qx(i,j,:)
     .   +(yc-ycen)*qy(i,j,:)
     .   + 0.5d0*qxx(i,j,:)*( (xc-xcen)**2-(dx**2)/12.d0)
     .   + (xc-xcen)*(yc-ycen)*qxy(i,j,:)
     .   + 0.5d0*qyy(i,j,:)*( (yc-ycen)**2-(dy**2)/12.d0)

       write(231,*)  'ST ( ', x1,",",
     .                        y1,",0,",
     .                        x2,",",
     .                        y2,",0,",
     .                        x3,",",
     .                        y3,",0)",
     .                        "{",
     .                        recon(ivar,1),",",
     .                        recon(ivar,2),",",
     .                        recon(ivar,3),
     .                        "};"
11    continue


         endif
 15   continue
      write(231,*) "};"













!      write(231,*) 'View "error" {'
!      do 915 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
!      do 915 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
!         kirr = irr(i,j)
!         if (kirr .eq. -1) go to 915
!         if (kirr .eq. lstgrd ) then
!            xcen = xlowb + (dfloat(i)-.5d0)* dx
!            ycen = ylowb + (dfloat(j)-.5d0)* dy
!
!             xcorner = xlowb + (dfloat(i)-1.)* dx
!             ycorner = ylowb + (dfloat(j)-1.)* dy
!             do itimes = 1,4
!                if (itimes .eq. 1) then ! get all 4 corners of mesh, pw constant sol
!                   xc = xcorner
!                   yc = ycorner
!                else if (itimes .eq. 2) then
!                   xc = xcorner
!                   yc = ycorner + dy
!                else if (itimes .eq. 3) then
!                   xc = xcorner + dx
!                   yc = ycorner + dy
!                else if (itimes .eq. 4) then
!                   xc = xcorner + dx
!                   yc = ycorner
!                endif
!
!            corners(itimes,1) = xc
!            corners(itimes,2) = yc
!
!              end do
!
!
!                call computed_error(i,j,dx,dy,xlow,ylow,
!     .          time, irr,
!     .          q, qx, qy, qxx, qxy, qyy,
!     .          qxxx, qxxy, qxyy, qyyy,
!     .          mitot, mjtot, nvar,
!     .          ivar, iprob, lstgrd,lwidth,
!     .          douterr, doutmaxerr, doutemag, doutmaxemag,
!     .          douterr2, doutmaxerr2, doutemag2, doutmaxemag2)
!
!       write(231,*)  'SQ ( ', corners(1,1),",",
!     .                        corners(1,2),",0,",
!     .                        corners(2,1),",",
!     .                        corners(2,2),",0,",
!     .                        corners(3,1),",",
!     .                        corners(3,2),",0,",
!     .                        corners(4,1),",",
!     .                        corners(4,2),",0)",
!     .                        "{",
!     .                        douterr,",",
!     .                        douterr,",",
!     .                        douterr,",",
!     .                        douterr,
!     .                        "};"
!
!
!         else
!            ivert = 1
!            do 9290 while (poly(ivert+1,1,kirr) .ne. -11.)
!               ivert = ivert + 1
! 9290       continue
!            itri = ivert - 3
!
!                  idx1 = 1
!                  do 991 it = 1, itri ! for each  triangle
!                    idx2 = it + 1
!                    idx3 = it + 2
!
!                    x1 = poly(idx1,1,kirr)
!                    y1 = poly(idx1,2,kirr)
!
!                    x2 = poly(idx2,1,kirr)
!                    y2 = poly(idx2,2,kirr)
!
!                    x3 = poly(idx3,1,kirr)
!                    y3 = poly(idx3,2,kirr)
!
!
!                xcen = xcirr(kirr)
!                ycen = ycirr(kirr)
!
!                call computed_error(i,j,dx,dy,xlow,ylow,
!     .          time, irr,
!     .          q, qx, qy, qxx, qxy, qyy,
!     .          qxxx, qxxy, qxyy, qyyy,
!     .          mitot, mjtot, nvar,
!     .          ivar, iprob, lstgrd,lwidth,
!     .          douterr, doutmaxerr, doutemag, doutmaxemag,
!     .          douterr2, doutmaxerr2, doutemag2, doutmaxemag2)
!
!       write(231,*)  'ST ( ', x1,",",
!     .                        y1,",0,",
!     .                        x2,",",
!     .                        y2,",0,",
!     .                        x3,",",
!     .                        y3,",0)",
!     .                        "{",
!     .                        douterr,",",
!     .                        douterr,",",
!     .                        douterr,
!     .                        "};"
!991   continue
!
!
!
!         endif
!915   continue
!      write(231,*) "};"
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!      write(231,*) 'View "numhoods" {'
!      do 215 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
!      do 215 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
!         kirr = irr(i,j)
!         if (kirr .eq. -1) go to 215
!         if (kirr .eq. lstgrd ) then
!            xcen = xlowb + (dfloat(i)-.5d0)* dx
!            ycen = ylowb + (dfloat(j)-.5d0)* dy
!
!             xcorner = xlowb + (dfloat(i)-1.)* dx
!             ycorner = ylowb + (dfloat(j)-1.)* dy
!             do itimes = 1,4
!                if (itimes .eq. 1) then ! get all 4 corners of mesh, pw constant sol
!                   xc = xcorner
!                   yc = ycorner
!                else if (itimes .eq. 2) then
!                   xc = xcorner
!                   yc = ycorner + dy
!                else if (itimes .eq. 3) then
!                   xc = xcorner + dx
!                   yc = ycorner + dy
!                else if (itimes .eq. 4) then
!                   xc = xcorner + dx
!                   yc = ycorner
!                endif
!
!            corners(itimes,1) = xc
!            corners(itimes,2) = yc
!
!            recon(:,itimes) = q(i,j,:)
!     .   +(xc-xcen)*qx(i,j,:)
!     .   +(yc-ycen)*qy(i,j,:)
!     .   + 0.5d0*qxx(i,j,:)*( (xc-xcen)**2-(dx**2)/12.d0)
!     .   + (xc-xcen)*(yc-ycen)*qxy(i,j,:)
!     .   + 0.5d0*qyy(i,j,:)*( (yc-ycen)**2-(dy**2)/12.d0)
!              end do
!
!       write(231,*)  'SQ ( ', corners(1,1),",",
!     .                        corners(1,2),",0,",
!     .                        corners(2,1),",",
!     .                        corners(2,2),",0,",
!     .                        corners(3,1),",",
!     .                        corners(3,2),",0,",
!     .                        corners(4,1),",",
!     .                        corners(4,2),",0)",
!     .                        "{",
!     .                        numhoods(i,j),",",
!     .                        numhoods(i,j),",",
!     .                        numhoods(i,j),",",
!     .                        numhoods(i,j),
!     .                        "};"
!
!
!         else
!            ivert = 1
!            do 2290 while (poly(ivert+1,1,kirr) .ne. -11.)
!               ivert = ivert + 1
! 2290       continue
!            itri = ivert - 3
!
!                  idx1 = 1
!                  do 2291 it = 1, itri ! for each  triangle
!                    idx2 = it + 1
!                    idx3 = it + 2
!
!                    x1 = poly(idx1,1,kirr)
!                    y1 = poly(idx1,2,kirr)
!
!                    x2 = poly(idx2,1,kirr)
!                    y2 = poly(idx2,2,kirr)
!
!                    x3 = poly(idx3,1,kirr)
!                    y3 = poly(idx3,2,kirr)
!
!
!                xcen = xcirr(kirr)
!                ycen = ycirr(kirr)
!
!            xc = x1
!            yc = y1
!            recon(:,1) = q(i,j,:)
!     .   +(xc-xcen)*qx(i,j,:)
!     .   +(yc-ycen)*qy(i,j,:)
!     .   + 0.5d0*qxx(i,j,:)*( (xc-xcen)**2-(dx**2)/12.d0)
!     .   + (xc-xcen)*(yc-ycen)*qxy(i,j,:)
!     .   + 0.5d0*qyy(i,j,:)*( (yc-ycen)**2-(dy**2)/12.d0)
!
!            xc = x2
!            yc = y2
!            recon(:,2) = q(i,j,:)
!     .   +(xc-xcen)*qx(i,j,:)
!     .   +(yc-ycen)*qy(i,j,:)
!     .   + 0.5d0*qxx(i,j,:)*( (xc-xcen)**2-(dx**2)/12.d0)
!     .   + (xc-xcen)*(yc-ycen)*qxy(i,j,:)
!     .   + 0.5d0*qyy(i,j,:)*( (yc-ycen)**2-(dy**2)/12.d0)
!
!            xc = x3
!            yc = y3
!            recon(:,3) = q(i,j,:)
!     .   +(xc-xcen)*qx(i,j,:)
!     .   +(yc-ycen)*qy(i,j,:)
!     .   + 0.5d0*qxx(i,j,:)*( (xc-xcen)**2-(dx**2)/12.d0)
!     .   + (xc-xcen)*(yc-ycen)*qxy(i,j,:)
!     .   + 0.5d0*qyy(i,j,:)*( (yc-ycen)**2-(dy**2)/12.d0)
!
!       write(231,*)  'ST ( ', x1,",",
!     .                        y1,",0,",
!     .                        x2,",",
!     .                        y2,",0,",
!     .                        x3,",",
!     .                        y3,",0)",
!     .                        "{",
!     .                        numhoods(i,j),",",
!     .                        numhoods(i,j),",",
!     .                        numhoods(i,j),
!     .                        "};"
!2291   continue
!
!
!
!
!         endif
! 215   continue
!      write(231,*) "};"
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!      write(231,*) 'View "ncount" {'
!      do 2159 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
!      do 2159 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
!         kirr = irr(i,j)
!         if (kirr .eq. -1) go to 2159
!         if (kirr .eq. lstgrd ) then
!            xcen = xlowb + (dfloat(i)-.5d0)* dx
!            ycen = ylowb + (dfloat(j)-.5d0)* dy
!
!             xcorner = xlowb + (dfloat(i)-1.)* dx
!             ycorner = ylowb + (dfloat(j)-1.)* dy
!             do itimes = 1,4
!                if (itimes .eq. 1) then ! get all 4 corners of mesh, pw constant sol
!                   xc = xcorner
!                   yc = ycorner
!                else if (itimes .eq. 2) then
!                   xc = xcorner
!                   yc = ycorner + dy
!                else if (itimes .eq. 3) then
!                   xc = xcorner + dx
!                   yc = ycorner + dy
!                else if (itimes .eq. 4) then
!                   xc = xcorner + dx
!                   yc = ycorner
!                endif
!
!            corners(itimes,1) = xc
!            corners(itimes,2) = yc
!
!            recon(:,itimes) = q(i,j,:)
!     .   +(xc-xcen)*qx(i,j,:)
!     .   +(yc-ycen)*qy(i,j,:)
!     .   + 0.5d0*qxx(i,j,:)*( (xc-xcen)**2-(dx**2)/12.d0)
!     .   + (xc-xcen)*(yc-ycen)*qxy(i,j,:)
!     .   + 0.5d0*qyy(i,j,:)*( (yc-ycen)**2-(dy**2)/12.d0)
!              end do
!
!       write(231,*)  'SQ ( ', corners(1,1),",",
!     .                        corners(1,2),",0,",
!     .                        corners(2,1),",",
!     .                        corners(2,2),",0,",
!     .                        corners(3,1),",",
!     .                        corners(3,2),",0,",
!     .                        corners(4,1),",",
!     .                        corners(4,2),",0)",
!     .                        "{",
!     .                        0,",",
!     .                        0,",",
!     .                        0,",",
!     .                        0,
!     .                        "};"
!
!
!         else
!            ivert = 1
!            do 22909 while (poly(ivert+1,1,kirr) .ne. -11.)
!               ivert = ivert + 1
!22909       continue
!            itri = ivert - 3
!
!                  idx1 = 1
!                  do 2299 it = 1, itri ! for each  triangle
!                    idx2 = it + 1
!                    idx3 = it + 2
!
!                    x1 = poly(idx1,1,kirr)
!                    y1 = poly(idx1,2,kirr)
!
!                    x2 = poly(idx2,1,kirr)
!                    y2 = poly(idx2,2,kirr)
!
!                    x3 = poly(idx3,1,kirr)
!                    y3 = poly(idx3,2,kirr)
!
!
!                xcen = xcirr(kirr)
!                ycen = ycirr(kirr)
!
!            xc = x1
!            yc = y1
!            recon(:,1) = q(i,j,:)
!     .   +(xc-xcen)*qx(i,j,:)
!     .   +(yc-ycen)*qy(i,j,:)
!     .   + 0.5d0*qxx(i,j,:)*( (xc-xcen)**2-(dx**2)/12.d0)
!     .   + (xc-xcen)*(yc-ycen)*qxy(i,j,:)
!     .   + 0.5d0*qyy(i,j,:)*( (yc-ycen)**2-(dy**2)/12.d0)
!
!            xc = x2
!            yc = y2
!            recon(:,2) = q(i,j,:)
!     .   +(xc-xcen)*qx(i,j,:)
!     .   +(yc-ycen)*qy(i,j,:)
!     .   + 0.5d0*qxx(i,j,:)*( (xc-xcen)**2-(dx**2)/12.d0)
!     .   + (xc-xcen)*(yc-ycen)*qxy(i,j,:)
!     .   + 0.5d0*qyy(i,j,:)*( (yc-ycen)**2-(dy**2)/12.d0)
!
!            xc = x3
!            yc = y3
!            recon(:,3) = q(i,j,:)
!     .   +(xc-xcen)*qx(i,j,:)
!     .   +(yc-ycen)*qy(i,j,:)
!     .   + 0.5d0*qxx(i,j,:)*( (xc-xcen)**2-(dx**2)/12.d0)
!     .   + (xc-xcen)*(yc-ycen)*qxy(i,j,:)
!     .   + 0.5d0*qyy(i,j,:)*( (yc-ycen)**2-(dy**2)/12.d0)
!
!       write(231,*)  'ST ( ', x1,",",
!     .                        y1,",0,",
!     .                        x2,",",
!     .                        y2,",0,",
!     .                        x3,",",
!     .                        y3,",0)",
!     .                        "{",
!     .                        ncount(kirr),",",
!     .                        ncount(kirr),",",
!     .                        ncount(kirr),
!     .                        "};"
!2299   continue
!
!
!
!
!         endif
! 2159   continue
!      write(231,*) "};"















      close(231)


!      if(iprob .eq. 11) then
!      call p11plot(q, qx, qy,lwidth, mitot, mjtot, nvar, lstgrd, irr)
!      endif


      return
      end
