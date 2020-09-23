c
c -------------------------------------------------------------
c
      subroutine outtec(qin,nvar,mptr,irr,
     1                  mitot,mjtot,qx,qy,lwidth,lstgrd,
     2                  dx,dy,xlow,ylow,time)
c
      implicit double precision (a-h,o-z)
      include "cirr.i"
      include "./reconinfo.i"
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
      dimension recon(nvar)
      dimension numvert(irrsize)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .         ismp,gradThreshold, ilts,ibc,imodal,
     .id2,res(ijsize,ijsize,4)
      common /order2/ ssw, quad, nolimiter
      logical  eout, quad, nolimiter
      logical pwconst

      dimension dir(2)

      data      spval/-1.e10/, eout/.true./
c      character*23  filename 
      common /dump/ uvel,vvel
      dimension uvel(849,848),vvel(848,849)
c
      ibuff = 4
      if(ibc .eq. 1) ibuff = 0

      xlowb = xlow - dfloat(lwidth)*dx
      ylowb = ylow - dfloat(lwidth)*dy

      qx = 0.d0
      qy = 0.d0
      qxx = 0.d0
      qxy = 0.d0
      qyy = 0.d0
      qxxx = 0.d0
      qxxy = 0.d0
      qxyy = 0.d0
      qyyy = 0.d0

      q = qin





!
!      call vctoprm(q,q,mitot,mjtot)
c
c     ### call for exterior bcs at each stage so can use slopes
c    ## NOTE THAT BNDRY CELLS FROM OTHER GRIDS NOT SET
            xhigh = xlowb + mitot*dx
            yhigh = ylowb + mjtot*dy
            call pphysbdlin(xlowb,xhigh,ylowb,yhigh,level,mitot,mjtot,
     &                   nvar,q,time,dx,dy,qx,qy,irr,lstgrd)


!      print *,q(4,23,1),q(5,23,1)

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


!      print *,q(4,23,1),q(5,23,1)

      elseif(quad .and. (nolimiter .eqv. .false.) ) then ! yes in primitive variables, and with limiter
                                                      ! only works for p = 1

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
!     .             mitot, mjtot, irr, nvar,dx, dy,
!     .             lwidth, lstgrd, xlowb, ylowb)
!      endif

      elseif (ssw .ne. 0.d0) then   ! conserved variables
         call slopes(q,qx,qy,qxx,qxy,qyy,
     &               qxxx, qxxy, qxyy, qyyy,
     &               mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &               xlowb,ylowb,nvar)
      endif




      do 115 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 115 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
      kirr = irr(i,j)
      if (kirr .eq. -1 .or. kirr .eq. lstgrd) go to 115
          ivert = 1
          do while (poly(ivert+1,1,kirr) .ne. -11.)
            ivert = ivert + 1
          end do
          numvert(kirr) = ivert-1
 115    continue

!      if(iprob .eq. 28) then
!
!        if(ssw .eq. 2) then
!            numvert(kirr) = numvert(kirr) + 1
!        elseif(ssw .eq. 3) then
!            numvert(kirr) = numvert(kirr) + 2
!        endif
!
!      endif


c  count needed for unstructured tec format (so dont have to look up new format)
      numelem = 0
      numnodes = 0
      numfaces = 0
      do 10 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 10 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) then
            cycle
         elseif (kirr .eq. lstgrd) then
            numelem = numelem+1
            numnodes = numnodes + 4
            numfaces = numfaces + 4
         else
            numelem = numelem +1
            numnodes = numnodes + numvert(kirr)
            numfaces = numfaces + numvert(kirr)

          if(ihob .eq. 1) then

            if(ssw .eq. 2 .or. ssw .eq .-2) then
                numnodes = numnodes + 1
                numfaces = numfaces + 1
            elseif(ssw .eq. 3  .or. ssw .eq. -3)  then
                numnodes = numnodes + 2
                numfaces = numfaces + 2
            endif

          endif

         endif
 10   continue
c
      write(14,99) numnodes, numelem,numfaces
 99   format('VARIABLES = "X","Y","Rho","error","Xcent","Ycent",
     ."i", "j", "q", "qx", "qy", "qxx","qxy", "qyy", "exact",
     . "numhoods", "inuf", "mioff", "mjoff" , "iir", "jjr", "nc",
     . "kirr",  "Xcent_ho","Ycent_ho", "shift_xx","shift_xy","shift_yy",
     . "res"'
     .,/,'ZONE T="CUTCELLS" ',/,
     .'ZONETYPE=FEPOLYGON',/,'Nodes    =',i8,/,'Elements =',i8,/,
     .'Faces    =',i8,/,'NumConnectedBoundaryFaces=0',/,
     .'TotalNumBoundaryConnections=0')

      write(14,88)
 88   format('# X Data')
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
                write(14,141, advance="no") xc
 141            format(8e25.15)
              end do
         else
            do iv = 1, numvert(kirr)
                write(14,981, advance="no") poly(iv,1,kirr)
 981            format(8e25.15)
            end do
            do iv = 2,nbdry-1
                write(14,981, advance="no") bdry(iv,1,kirr)
            enddo

         endif
      write(14,*)
 15   continue

      write(14,*)
      write(14,77)
 77   format('# Y Data')
      do 16 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 16 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 16
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
                write(14,142, advance="no") yc
 142            format(8e25.15)
              end do
        else
            do iv = 1, numvert(kirr)
                write(14,91, advance="no") poly(iv,2,kirr)
 91            format(8e25.15)
            end do
            do iv = 2,nbdry-1
                write(14,981, advance="no") bdry(iv,2,kirr)
            enddo

         endif

      write(14,*)
 16   continue




      ivar = 1
      write(14,*)
      write(14,66)
 66   format('# Rho Data')
      do 17 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 17 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 17
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

       deltax = xc-xcen
       deltay = yc-ycen
       call evalU(valprim,deltax, deltay, dx, dy,
     .            q, qx, qy,
     .            qxx, qxy, qyy,
     .            qxxx, qxxy, qxyy, qyyy,
     .            i,j, lstgrd, mitot, mjtot, nvar, ivar)
!            valprim(ivar) = q(i,j,ivar)
!     .   +(xc-xcen)*qx(i,j,ivar)
!     .   +(yc-ycen)*qy(i,j,ivar)
!     .   + 0.5d0*qxx(i,j,ivar)*( (xc-xcen)**2-(dx**2)/12.d0)
!     .   + (xc-xcen)*(yc-ycen)*qxy(i,j,ivar)
!     .   + 0.5d0*qyy(i,j,ivar)*( (yc-ycen)**2-(dy**2)/12.d0)

                write(14,143, advance="no") valprim(ivar)
 143            format(8e25.15)
              end do
      else
            do iv = 1, numvert(kirr)
                if(ihob .eq. 0) then
                xcen = xcirr(kirr)
                ycen = ycirr(kirr)
                else
                xcen = xcirr_ho(kirr)
                ycen = ycirr_ho(kirr)
                endif

                xc = poly(iv,1,kirr)
                yc = poly(iv,2,kirr)

               deltax = xc-xcen
               deltay = yc-ycen
       call evalU(valprim,deltax, deltay, dx, dy,
     .            q, qx, qy,
     .            qxx, qxy, qyy,
     .            qxxx, qxxy, qxyy, qyyy,
     .            i,j, kirr, mitot, mjtot, nvar, ivar)
                write(14,144, advance="no") valprim(ivar)
 144            format(8e25.15)
            end do
            do iv = 2,nbdry-1
                if(ihob .eq. 0) then
                xcen = xcirr(kirr)
                ycen = ycirr(kirr)
                else
                xcen = xcirr_ho(kirr)
                ycen = ycirr_ho(kirr)
                endif

                xc = bdry(iv,1,kirr)
                yc = bdry(iv,2,kirr)

                deltax = xc-xcen
                deltay = yc-ycen
                call evalU(valprim,deltax, deltay, dx, dy,
     .                     q, qx, qy,
     .                     qxx, qxy, qyy,
     .                     qxxx, qxxy, qxyy, qyyy,
     .                     i,j, kirr, mitot, mjtot, nvar, ivar)
                write(14,144, advance="no") valprim(ivar)
            enddo



      endif
              write(14,*)


 17   continue




      dmag = 0.d0
      dl1error = 0.d0
      dmaxerror = -1.d0
      dmaxmag = -1.d0
      vol = 0.d0

      dmag2 = 0.d0
      dl1error2 = 0.d0
      dmaxerror2 = -1.d0
      dmaxmag2 = -1.d0

      fv_vol_err = 0.d0


      write(14,*)
      write(14,661)
 661   format('# error Data')
      do 171 i = lwidth+1-ibuff, mitot-lwidth+ibuff
      do 171 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) cycle

                call computed_error(i,j,dx,dy,xlow,ylow,
     .          time, irr,
     .          q, qx, qy, qxx, qxy, qyy,
     .          qxxx, qxxy, qxyy, qyyy,
     .          mitot, mjtot, nvar,
     .          ivar, iprob, lstgrd,lwidth,
     .          douterr, doutmaxerr, doutmag, doutmaxmag,
     .          doutvol, dfv_vol_err)

            if( i .ge. lwidth+1 .and. i .le. mitot-lwidth .and.
     .          j .ge. lwidth+1 .and. j .le. mjtot-lwidth) then

            fv_vol_err = fv_vol_err + dfv_vol_err

            dl1error = dl1error+douterr
            dmag = dmag+doutmag
            vol = vol + doutvol
            if(dmaxerror < doutmaxerr) dmaxerror = doutmaxerr
            if(dmaxmag < doutmaxmag) dmaxmag = doutmaxmag
        endif

         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
!                write(14,431, advance="no") doutmaxerr
                write(14,431, advance="no") douterr/(dx*dy)
 431            format(8e25.15)
              end do

        else
            do iv = 1, numvert(kirr)
!               write(14,441, advance="no") doutmaxerr
                write(14,431, advance="no") douterr/ar(kirr)
 441            format(8e25.15)
            end do
            do iv = 2,nbdry-1
                write(14,431, advance="no") douterr/ar(kirr)
            enddo



      endif
              write(14,*)


 171   continue


      write(14,*)
      write(14,929)
 929   format('# Xcent Data')
      do 420 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 420 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 420
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,921, advance="no") xlowb + (dfloat(i)-.5d0)* dx
 921            format(8e25.15)
              end do
      else
            do iv = 1, numvert(kirr)
                write(14,922, advance="no") xcirr(kirr)
 922            format(8e25.15)
            end do
            do iv = 2,nbdry-1
                write(14,922, advance="no") xcirr(kirr)
            enddo

      endif
              write(14,*)


 420   continue





      write(14,*)
      write(14,939)
 939   format('# Ycent Data')
      do 430 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 430 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 430
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,931, advance="no") ylowb + (dfloat(j)-.5d0)* dy
 931            format(8e25.15)
              end do
      else
            do iv = 1, numvert(kirr)
                write(14,932, advance="no") ycirr(kirr)
 932            format(8e25.15)
            end do
            do iv = 2,nbdry-1
                write(14,932, advance="no") ycirr(kirr)
            enddo


      endif
              write(14,*)


 430   continue

















      write(14,*)
      write(14,199)
 199   format('# i Data')
      do 11 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 11 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 11
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,1431, advance="no") i
 1431            format(i12)
              end do
      else
            do iv = 1, numvert(kirr)
                write(14,1441, advance="no") i
 1441            format(i12)
            end do
            do iv = 2,nbdry-1
                write(14,1441, advance="no") i
            enddo
      endif
              write(14,*)


 11   continue



      write(14,*)
      write(14,299)
 299   format('# j Data')
      do 111 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 111 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 111
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,113, advance="no") j
 113            format(i12)
              end do
      else
            do iv = 1, numvert(kirr)
                write(14,124, advance="no") j
 124            format(i12)
            end do
            do iv = 2,nbdry-1
                write(14,124, advance="no") j
            enddo

      endif
              write(14,*)


 111   continue






      write(14,*)
      write(14,919)
 919   format('# q Data')
      do 410 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 410 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
!      if( i .eq. 4 .and. j .eq. 24) then
!      print *, "herE"
!      endif
!      if( i .eq. 5 .and. j .eq. 24) then
!      print *, "herE"
!      endif



         kirr = irr(i,j)
         if (kirr .eq. -1) go to 410
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,991, advance="no") q(i,j,ivar)
 911            format(8e25.15)
              end do
      else
            do iv = 1, numvert(kirr)
                write(14,992, advance="no") q(i,j,ivar)
 912            format(8e25.15)
            end do
            do iv = 2,nbdry-1
                write(14,992, advance="no") q(i,j,ivar)
            enddo

      endif
              write(14,*)


 410   continue








      write(14,*)
      write(14,999)
 999   format('# qx Data')
      do 110 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 110 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 110
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,991, advance="no") qx(i,j,ivar)
 991            format(8e25.15)
              end do
      else
            do iv = 1, numvert(kirr)
                write(14,992, advance="no") qx(i,j,ivar)
 992            format(8e25.15)
            end do
            do iv = 2,nbdry-1
                write(14,992, advance="no") qx(i,j,ivar)
            enddo

      endif
              write(14,*)


 110   continue




      write(14,*)
      write(14,12)
 12    format('# qy Data')
      do 109 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 109 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 109
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,881, advance="no") qy(i,j,ivar)
 881            format(8e25.15)
              end do
      else
            do iv = 1, numvert(kirr)
                write(14,882, advance="no") qy(i,j,ivar)
 882            format(8e25.15)
            end do
            do iv = 2,nbdry-1
                write(14,882, advance="no") qy(i,j,ivar)
            enddo

      endif
              write(14,*)


 109   continue


      write(14,*)
      write(14,13)
 13    format('# qxx Data')
      do 131 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 131 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 131
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,883, advance="no") qxx(i,j,ivar)
 883            format(8e25.15)
              end do
      else
            do iv = 1, numvert(kirr)
                write(14,884, advance="no") qxx(i,j,ivar)
 884            format(8e25.15)
            end do
            do iv = 2,nbdry-1
                write(14,884, advance="no") qxx(i,j,ivar)
            enddo
      endif
              write(14,*)


 131   continue


      write(14,*)
      write(14,21)
 21    format('# qxy Data')
      do 118 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 118 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 118
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,116, advance="no") qxy(i,j,ivar)
 116            format(8e25.15)
              end do
      else
            do iv = 1, numvert(kirr)
                write(14,117, advance="no") qxy(i,j,ivar)
 117            format(8e25.15)
            end do
            do iv = 2,nbdry-1
                write(14,117, advance="no") qxy(i,j,ivar)
            enddo
      endif
              write(14,*)


 118   continue



      write(14,*)
      write(14,119)
 119    format('# qyy Data')
      do 122 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 122 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 122
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,120, advance="no") qyy(i,j,ivar)
 120            format(8e25.15)
              end do
      else
            do iv = 1, numvert(kirr)
                write(14,121, advance="no") qyy(i,j,ivar)
 121            format(8e25.15)
            end do
            do iv = 2,nbdry-1
                write(14,121, advance="no") qyy(i,j,ivar)
            enddo
      endif
              write(14,*)


 122   continue



      ivar = 1
      write(14,*)
      write(14,626)
 626   format('# Exact Data')
      do 127 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 127 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 127
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


                call f(valprim, xc, yc, time, iprob)
                write(14,1243, advance="no") valprim(ivar)
 1243            format(8e25.15)
              end do
      else
            do iv = 1, numvert(kirr)
                xc = poly(iv,1,kirr)
                yc = poly(iv,2,kirr)

                call f(valprim, xc, yc, time, iprob)
                write(14,1244, advance="no") valprim(ivar)
 1244            format(8e25.15)
            end do
            do iv = 2,nbdry-1
                xc = bdry(iv,1,kirr)
                yc = bdry(iv,2,kirr)

                call f(valprim, xc, yc, time, iprob)
                write(14,1244, advance="no") valprim(ivar)
            enddo
      endif
              write(14,*)


 127   continue












      write(14,*)
      write(14,715)
 715    format('# numhoods Data')
      do 209 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 209 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 209
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,114, advance="no") numhoods(i,j)

              end do
      else
            do iv = 1, numvert(kirr)
                write(14,114, advance="no") numhoods(i,j)
            end do
            do iv = 2,nbdry-1
                write(14,114, advance="no") numhoods(i,j)
            enddo
      endif
              write(14,*)


 209   continue









      write(14,*)
      write(14,915)
 915    format('# inuf Data')
      do 2994 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 2994 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 2994
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,114, advance="no") inuf(i,j)

              end do
      else
            do iv = 1, numvert(kirr)
                write(14,114, advance="no") inuf(i,j)
            end do
            do iv = 2,nbdry-1
                write(14,114, advance="no") inuf(i,j)
            enddo
      endif
              write(14,*)


 2994   continue









      write(14,*)
      write(14,995)
 995    format('# mioff Data')
      do 24 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 24 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 24
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,114, advance="no") mioff(i,j)

              end do
      else
            do iv = 1, numvert(kirr)
                write(14,114, advance="no") mioff(i,j)
            end do
            do iv = 2,nbdry-1
                write(14,114, advance="no") mioff(i,j)
            enddo
      endif
              write(14,*)


 24   continue

      write(14,*)
      write(14,935)
 935    format('# mjoff Data')
      do 23 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 23 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 23
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,114, advance="no") mjoff(i,j)

              end do
      else
            do iv = 1, numvert(kirr)
                write(14,114, advance="no") mjoff(i,j)
            end do
            do iv = 2,nbdry-1
                write(14,114, advance="no") mjoff(i,j)
            enddo
      endif
              write(14,*)


 23   continue
















      write(14,*)
      write(14,9435)
 9435    format('# iir Data')
      do 246 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 246 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 246
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,114, advance="no") iir(i,j)

              end do
      else
            do iv = 1, numvert(kirr)
                write(14,114, advance="no") iir(i,j)
            end do
            do iv = 2,nbdry-1
                write(14,114, advance="no") iir(i,j)
            enddo
      endif
              write(14,*)


 246   continue

      write(14,*)
      write(14,1414)
 1414    format('# jjr Data')
      do 63 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 63 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 63
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,114, advance="no") jjr(i,j)

              end do
      else
            do iv = 1, numvert(kirr)
                write(14,114, advance="no") jjr(i,j)
            end do
            do iv = 2,nbdry-1
                write(14,114, advance="no") jjr(i,j)
            enddo
      endif
              write(14,*)


 63   continue

















      write(14,*)
      write(14,4)
 4    format('# nc Data')
      do 613 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 613 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 613
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,114, advance="no") 1

              end do
      else
            do iv = 1, numvert(kirr)
                write(14,114, advance="no") ncount(kirr)
            end do
            do iv = 2,nbdry-1
                write(14,114, advance="no") ncount(kirr)
            enddo
      endif
              write(14,*)


 613   continue



      write(14,*)
      write(14,8)
 8    format('# kirr Data')
      do 6813 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 6813 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 6813
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,114, advance="no") lstgrd

              end do
      else
            do iv = 1, numvert(kirr)
                write(14,114, advance="no") kirr
            end do
            do iv = 2,nbdry-1
                write(14,114, advance="no") kirr
            enddo
      endif
              write(14,*)


 6813   continue














      write(14,*)
      write(14,949)
 949   format('# Xcent_ho Data')
      do 440 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 440 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 440
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,941, advance="no") xlowb + (dfloat(i)-.5d0)* dx
 941            format(8e25.15)
              end do
      else
            do iv = 1, numvert(kirr)
                write(14,942, advance="no") xcirr_ho(kirr)
 942            format(8e25.15)
            end do
            do iv = 2,nbdry-1
                write(14,942, advance="no") xcirr_ho(kirr)
            enddo

      endif
              write(14,*)


 440   continue





      write(14,*)
      write(14,959)
 959   format('# Ycent_ho Data')
      do 450 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 450 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 450
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,951, advance="no") ylowb + (dfloat(j)-.5d0)* dy
 951            format(8e25.15)
              end do
      else
            do iv = 1, numvert(kirr)
                write(14,952, advance="no") ycirr_ho(kirr)
 952            format(8e25.15)
            end do
            do iv = 2,nbdry-1
                write(14,952, advance="no") ycirr_ho(kirr)
            enddo


      endif
              write(14,*)


 450   continue








      write(14,*)
      write(14,962)
 962   format('# shift_xx Data')
      do 453 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 453 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 453
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,951, advance="no") ylowb + (dfloat(j)-.5d0)* dy
              end do
      else
            do iv = 1, numvert(kirr)

            if(ihob .eq. 0) then
                write(14,952, advance="no") poly(8 ,1,kirr)
            else
                write(14,952, advance="no") poly_ho(8 ,1,kirr)
            endif

            end do
            do iv = 2,nbdry-1

            if(ihob .eq. 0) then
                write(14,952, advance="no") poly(8 ,1,kirr)
            else
                write(14,952, advance="no") poly_ho(8 ,1,kirr)
            endif

            enddo


      endif
              write(14,*)


 453   continue


      write(14,*)
      write(14,961)
 961   format('# shift_xy Data')
      do 451 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 451 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 451
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,951, advance="no") ylowb + (dfloat(j)-.5d0)* dy
              end do
      else
            do iv = 1, numvert(kirr)
            if(ihob .eq. 0) then
                write(14,952, advance="no") poly(9 ,1,kirr)
            else
                write(14,952, advance="no") poly_ho(9 ,1,kirr)
            endif
            end do
            do iv = 2,nbdry-1
            if(ihob .eq. 0) then
                write(14,952, advance="no") poly(9 ,1,kirr)
            else
                write(14,952, advance="no") poly_ho(9 ,1,kirr)
            endif
            enddo


      endif
              write(14,*)


 451   continue



      write(14,*)
      write(14,960)
 960   format('# shift_yy Data')
      do 452 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 452 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 452
         if (kirr .eq. lstgrd ) then
             do itimes = 1,4
                write(14,951, advance="no") ylowb + (dfloat(j)-.5d0)* dy
              end do
      else
            do iv = 1, numvert(kirr)
            if(ihob .eq. 0) then
                write(14,952, advance="no") poly(10 ,1,kirr)
            else
                write(14,952, advance="no") poly_ho(10 ,1,kirr)
            endif
            end do
            do iv = 2,nbdry-1
            if(ihob .eq. 0) then
                write(14,952, advance="no") poly(10 ,1,kirr)
            else
                write(14,952, advance="no") poly_ho(10 ,1,kirr)
            endif
            enddo


      endif
              write(14,*)


 452   continue












      dmaxresfull = -1.d0
      write(14,*)
      write(14,969)
 969   format('# res Data')
      do 454 i = lwidth+1-ibuff, mitot-lwidth+ibuff
      do 454 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 454



         if (kirr .eq. lstgrd ) then
             dmaxresfull = max(dmaxresfull, dabs(res(i,j,1)))
             do itimes = 1,4
                write(14,951, advance="no") dabs(res(i,j,1))
              end do
         else
            do iv = 1, numvert(kirr)
            if(ihob .eq. 0) then
                write(14,952, advance="no") dabs(res(i,j,1))
            else
                write(14,952, advance="no") dabs(res(i,j,1))
            endif
            end do
            do iv = 2,nbdry-1
            if(ihob .eq. 0) then
                write(14,952, advance="no") dabs(res(i,j,1))
            else
                write(14,952, advance="no") dabs(res(i,j,1))
            endif
            enddo


      endif
              write(14,*)


 454   continue
       print *, "Maximum residual is: ", dmaxresfull
!      ivar = 1
!      write(14,*)
!      write(14,67)
! 67   format('# u Data')
!      do 51 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
!      do 51 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
!         kirr = irr(i,j)
!         if (kirr .eq. -1) go to 51
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
!                call getdir(dir, xc, yc)
!
!
!                write(14,143, advance="no") dir(1)
!              end do
!      else
!            do iv = 1, numvert(kirr)
!                if(ihob .eq. 0) then
!                xcen = xcirr(kirr)
!                ycen = ycirr(kirr)
!                else
!                xcen = xcirr_ho(kirr)
!                ycen = ycirr_ho(kirr)
!                endif
!
!                xc = poly(iv,1,kirr)
!                yc = poly(iv,2,kirr)
!
!                call getdir(dir, xc, yc)
!                write(14,144, advance="no") dir(1)
!            end do
!            do iv = 2,nbdry-1
!                if(ihob .eq. 0) then
!                xcen = xcirr(kirr)
!                ycen = ycirr(kirr)
!                else
!                xcen = xcirr_ho(kirr)
!                ycen = ycirr_ho(kirr)
!                endif
!
!                xc = bdry(iv,1,kirr)
!                yc = bdry(iv,2,kirr)
!                call getdir(dir, xc, yc)
!                write(14,144, advance="no") dir(1)
!            enddo
!
!
!
!      endif
!              write(14,*)
!
!
! 51   continue
!
!      ivar = 1
!      write(14,*)
!      write(14,68)
! 68   format('# v Data')
!      do 41 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
!      do 41 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
!         kirr = irr(i,j)
!         if (kirr .eq. -1) go to 41
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
!                call getdir(dir, xc, yc)
!
!
!                write(14,143, advance="no") dir(2)
!              end do
!      else
!            do iv = 1, numvert(kirr)
!                if(ihob .eq. 0) then
!                xcen = xcirr(kirr)
!                ycen = ycirr(kirr)
!                else
!                xcen = xcirr_ho(kirr)
!                ycen = ycirr_ho(kirr)
!                endif
!
!                xc = poly(iv,1,kirr)
!                yc = poly(iv,2,kirr)
!
!                call getdir(dir, xc, yc)
!                write(14,144, advance="no") dir(2)
!            end do
!            do iv = 2,nbdry-1
!                if(ihob .eq. 0) then
!                xcen = xcirr(kirr)
!                ycen = ycirr(kirr)
!                else
!                xcen = xcirr_ho(kirr)
!                ycen = ycirr_ho(kirr)
!                endif
!
!                xc = bdry(iv,1,kirr)
!                yc = bdry(iv,2,kirr)
!                call getdir(dir, xc, yc)
!                write(14,144, advance="no") dir(2)
!            enddo
!
!
!
!      endif
!              write(14,*)
!
!
! 41   continue























      idx = 1

      write(14,55)
 55   format('# Face Nodes List')
      do 18 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 18 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 18
         if (kirr .eq. lstgrd ) then


            write(14, *) idx + 0, idx + 1
            write(14, *) idx + 1, idx + 2
            write(14, *) idx + 2, idx + 3
            write(14, *) idx + 3, idx + 0
            idx = idx+4
         else
            idxorig = idx
            do iv = 1, numvert(kirr)-1
                write(14,*) idx+0, idx+1
                idx = idx + 1
            end do
            do iv = 2,nbdry-1
                write(14,*) idx+0, idx+1
                idx = idx + 1
            enddo
            write(14,*) idx, idxorig
            idx = idx + 1
         endif

 18   continue

      icount = 1
      write(14,*)
      write(14,44)
 44   format('# Face Left Elements')
      do 19 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 19 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 19
         if (kirr .eq. lstgrd ) then


            write(14, *) icount, icount, icount, icount
            icount = icount+1
         else
            do iv = 1, numvert(kirr)
                write(14,114, advance="no") icount
 114            format(i12)
            end do
            do iv = 2,nbdry-1
                write(14,114, advance="no") icount
            enddo

            write(14, *)
            icount = icount + 1
         endif

 19   continue

      write(14,*)
      write(14,33)
 33   format('# Face Right Elements')
      do 20 i = lwidth+1 - ibuff, mitot-lwidth+ibuff
      do 20 j = lwidth+1-ibuff, mjtot-lwidth+ibuff
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 20
         if (kirr .eq. lstgrd ) then


            write(14, *) 0,0,0,0
         else
            do iv = 1, numvert(kirr)
                write(14,104, advance="no") 0
 104            format(i12)
            end do
            do iv = 2,nbdry-1
                write(14,104, advance="no") 0
            enddo
            write(14, *)
         endif

 20   continue






      dboundl1err = 0.d0
      dboundl1mag = 0.d0
      dboundmaxerr = -1.d0
      dboundmaxmag = -1.d0
      perimeter = 0.d0
      fv_boundary_err = 0.d0


!      dbounderr2 = 0.d0
!      dboundmaxerr2 = -1.d0
!      dboundmag2 = 0.d0
!      dboundmaxmag2 = -1.d0

      do 1711 i = lwidth+1, mitot-lwidth
      do 1711 j = lwidth+1, mjtot-lwidth
         kirr = irr(i,j)
         if (kirr .eq. -1) cycle
         if (kirr .eq. lstgrd ) cycle

                if(ihob .eq. 0) then
                call computed_boundary_lo(i,j,dx,dy,xlow,ylow,
     .          time, irr,
     .          q, qx, qy, qxx, qxy, qyy,
     .          qxxx, qxxy, qxyy, qyyy,
     .          mitot, mjtot, nvar,
     .          ivar, iprob, lstgrd,lwidth,
     .          douterr, doutmaxerr, doutmag, doutmaxmag,
     .          dlen, dfv_bound_err)

                fv_boundary_err = fv_boundary_err + dfv_bound_err
                else
                call computed_boundary_ho(i,j,dx,dy,xlow,ylow,
     .          time, irr,
     .          q, qx, qy, qxx, qxy, qyy,
     .          qxxx, qxxy, qxyy, qyyy,
     .          mitot, mjtot, nvar,
     .          ivar, iprob, lstgrd,lwidth,
     .          douterr, doutmaxerr, doutmag, doutmaxmag,
     .          dlen)
                endif
                dboundl1err = dboundl1err+douterr

                if(dboundmaxerr < doutmaxerr) dboundmaxerr = doutmaxerr
                perimeter = perimeter + dlen


 1711   continue


      print *,""
      print *,""
      print *,"VOLUME ERROR"
      print *,"------------"
      print *, "L1 volume error = ", dl1error
      print *, "maximum volume error = ", dmaxerror
      print *, "area = ", vol

      print *,""
      print *,""
      print *, "BOUNDARY ERROR"
      print *, "--------------"
      print *, "L1 boundary error = ", dboundl1err
      print *, "maximum boundary error = ", dboundmaxerr
      print *, "perimeter = ", perimeter
      print *, ""

      if(ihob .eq. 0) then
      print *,""
      print *,""
      print *,"FV ERRORs"
      print *,"------------"
      print *, "L1 volume error = ", fv_vol_err
      print *, "L1 boundary error = ", fv_boundary_err
      endif

!      print *,""
!      print *,""
!      print *,"RELATIVE VOLUME ERROR"
!      print *, "--------------------"
!      print *, "relative L1 volume error = ", dl1error/dmag

!                dboundl1mag = dboundl1mag+doutmag
!                if(dboundmaxmag < doutmaxmag) dboundmaxmag = doutmaxmag
!
!                dbounderr2 = dbounderr2+douterr2
!                dboundmag2 = dboundmag2+doutmag2
!                if(dboundmaxerr2 < doutmaxerr2) then
!                dboundmaxerr2 = doutmaxerr2
!                endif
!                if(dboundmaxmag2 < doutmaxmag2) then
!                dboundmaxmag2 = doutmaxmag2
!                endif






!      print *, "relative maximum error = ", dmaxerror/dmaxmag
!      print *, "dmag = ", dmag
!      print *, "dmaxmag = ", dmaxmag
!      if(dboundmag > 0.d0) then
!      print *, "L1 relative boundary error = ", dbounderr/dboundmag
!      print *, "maximum relative boundary error = ",
!     . dboundmaxerr/dboundmaxmag
!      print *, "dboundmag = ", dboundmag
!      print *, "dboundmaxmag = ", dboundmaxmag
!      endif



!      print *,""
!      print *,""
!      print *,""
!      print *,"Aftomis's version of relative norm"
!      print *, "L1 error weighted = ", dl1error2
!      print *, "maximum error weighted = ", dmaxerror2
!      print *, "L1 boundary error weighted = ", dbounderr2
!      print *, "maximum boundary error weighted = ", dboundmaxerr2
!      print *,""

!      print *, "relative L1 error weighted = ", dl1error2/dmag2
!      print *, "relative maximum error weighted = ", dmaxerror2/dmaxmag2
!      print *, "dmag weighted (area) = ", dmag2
!      print *, "dmaxmag weighted = ", dmaxmag2
!      print *,""
!
!      print *, "L1 relative boundary error weighted = ",
!     . dbounderr2/dboundmag2
!      print *, "maximum relative boundary error weighted = ",
!     . dboundmaxerr2/dboundmaxmag2
!      if(dboundmag2 > 0.d0) then
!      print *, "dboundmag2 weighted (perimeter) = ", dboundmag2
!      print *, "dboundmaxmag2 weighted = ", dboundmaxmag2
!      endif

      close(14)






      if(iprob .eq. 11) then
      print *, "Plotting Density on Cylinder!"
      do 215 i = lwidth+1, mitot-lwidth
      do 215 j = lwidth+1, mjtot-lwidth

         kirr = irr(i,j)
         if (kirr .eq. -1) go to 215
         if (kirr .eq. lstgrd) go to 215

         do 210 kside=1,6
            if (poly(kside+2,1,kirr).eq.-11.) then
               x1 = poly(kside,1,kirr)
               y1 = poly(kside,2,kirr)
               x2 = poly(kside+1,1,kirr)
               y2 = poly(kside+1,2,kirr)
               go to 25
            endif
 210      continue
 25      continue
         xcenter = 0.5d0
         ycenter = 0.5d0

         xm = 0.5d0*(x1+x2)
         ym = 0.5d0*(y1+y2)


         if(ihob .eq. 0) then
         xcen = xcirr(kirr)
         ycen = ycirr(kirr)
         else
         xcen = xcirr_ho(kirr)
         ycen = ycirr_ho(kirr)
         endif
!         xcen = xcirr(kirr)
!         ycen = ycirr(kirr)

         diffx = xm - xcen
         diffy = ym - ycen

         recon(:) = q(i,j,:) +(xm-xcen)*qx(i,j,:)+
     .                        (ym-ycen)*qy(i,j,:)

         val = atan2((ym-ycenter),(xm-xcenter))
         print *,val,recon

 215   continue
      endif



      if(iprob .eq. 21) then
      print *, "Plotting Density on wedge boundary!"
      do 225 i = lwidth+1, mitot-lwidth
      do 225 j = lwidth+1, mjtot-lwidth

         kirr = irr(i,j)
         if (kirr .eq. -1) go to 225
         if (kirr .eq. lstgrd) go to 225

         do 220 kside=1,6
            if (poly(kside+2,1,kirr).eq.-11.) then
               x1 = poly(kside,1,kirr)
               y1 = poly(kside,2,kirr)
               x2 = poly(kside+1,1,kirr)
               y2 = poly(kside+1,2,kirr)
               go to 211
            endif
 220      continue
 211      continue
         xcenter = 0.5d0
         ycenter = 0.5d0

         xm = 0.5d0*(x1+x2)
         ym = 0.5d0*(y1+y2)


         if(ihob .eq. 0) then
         xcen = xcirr(kirr)
         ycen = ycirr(kirr)
         else
         xcen = xcirr_ho(kirr)
         ycen = ycirr_ho(kirr)
         endif

         diffx = xm - xcen
         diffy = ym - ycen

         recon(:) = q(i,j,:) +(xm-xcen)*qx(i,j,:)+
     .                        (ym-ycen)*qy(i,j,:)

         print *,xm,recon

 225   continue
      endif








      return
      end


      subroutine computed_error(i,j,dx,dy,xlow,ylow,time, irr,
     .          q, qx, qy, qxx, qxy, qyy,
     .          qxxx, qxxy, qxyy, qyyy,
     .          mitot, mjtot, nvar,
     .          ivar, iprob, lstgrd, lwidth,
     .          d_error, d_maxerror, d_mag, d_maxmag,
     .          dvol, dfv_vol_err)
      implicit double precision (a-h,o-z)
      include "cirr.i"
      include "quadrature.i"
      dimension state(4)
      dimension err(4), dnumsol(4), abssol(4)
      dimension q(mitot,mjtot, nvar)
      dimension qx(mitot,mjtot, nvar),qy(mitot,mjtot, nvar)
      dimension qxx(mitot,mjtot, nvar),qxy(mitot,mjtot, nvar)
      dimension qyy(mitot,mjtot, nvar)
      dimension qxxx(mitot,mjtot, nvar),qxxy(mitot,mjtot, nvar)
      dimension qxyy(mitot,mjtot, nvar),qyyy(mitot,mjtot, nvar)
      dimension irr(mitot,mjtot)
      dimension dl1error(4), dmaxerror(4), d_exact_mag(4)
      dimension d_exact_max_mag(4)


      common /order2/ ssw, quad, nolimiter
      logical quad

      dvol = 0.d0
      dl1error = 0.d0
      dmaxerror = -1.d0
      d_exact_mag = 0.d0
      d_exact_max_mag = -1.d0

      kirr = irr(i,j)


      if (kirr .eq. lstgrd ) then
            xcen = xlow + (dfloat(i-lwidth)-.5d0)* dx
            ycen = ylow + (dfloat(j-lwidth)-.5d0)* dy
      else

                if(ihob .eq. 0) then
                xcen = xcirr(kirr)
                ycen = ycirr(kirr)
                else
                xcen = xcirr_ho(kirr)
                ycen = ycirr_ho(kirr)
                endif
!            xcen = xcirr(kirr)
!            ycen = ycirr(kirr)
      endif


          if (kirr .eq. lstgrd) then ! only full cells for now
             goto 1
          elseif (kirr .ne. -1 .and. kirr .ne. lstgrd
     .           .and. ihob .eq. 0) then ! cut cells lo
             goto 2
          elseif (kirr .ne. -1 .and. kirr .ne. lstgrd
     .           .and. ihob .eq. 1) then ! cut cells ho
             goto 3
          else
            goto 821
          endif

1         continue ! whole cell projection
          val = 0.d0



          xcen1 = xlow + (dfloat(i-lwidth)-1.d0)*dx
          ycen1 = ylow + (dfloat(j-lwidth)-1.d0)*dy

          xcen3 = xlow + (dfloat(i-lwidth)-0.d0)*dx
          ycen3 = ylow + (dfloat(j-lwidth)-0.d0)*dy






        xc = xlow + (dfloat(i-lwidth)-0.5d0)*dx
        yc = ylow + (dfloat(j-lwidth)-0.5d0)*dy
        call f(state, xc, yc, time, iprob)
        if(quad) then
             gamma1 = 0.4d0
             rhot = state(1)
             ut   = state(2)/state(1)
             vt   = state(3)/state(1)
             pt   = gamma1*( state(4)
     .    - 0.5d0 * (state(2)**2 + state(3)**2)/state(1) )
          err(1) = dabs( rhot - q(i,j,1) )
          err(2) = dabs( ut   - q(i,j,2) )
          err(3) = dabs( vt   - q(i,j,3) )
          err(4) = dabs( pt   - q(i,j,4) )
        else
          err(:) = dabs(state(:)- q(i,j,:))
        endif
        dfv_vol_err = dx*dy*err(1)
        err = 0.d0







        do 11 iq = 1, nquadquad

      xc = xcen1 * (1.d0 - rquad(iq))/2. + xcen3 * (1.d0 + rquad(iq))/2.
      yc = ycen1 * (1.d0 - squad(iq))/2. + ycen3 * (1.d0 + squad(iq))/2.


          call f(state, xc, yc, time, iprob)

       deltax = xc-xcen
       deltay = yc-ycen

       do mmm = 1,nvar
       call evalU(dnumsol,deltax, deltay, dx, dy,
     .            q, qx, qy,
     .            qxx, qxy, qyy,
     .            qxxx, qxxy, qxyy, qyyy,
     .            i,j, kirr, mitot, mjtot, nvar, mmm)
       end do


          if(quad) then
             gamma1 = 0.4d0
             rhot = state(1)
             ut   = state(2)/state(1)
             vt   = state(3)/state(1)
             pt   = gamma1*( state(4)
     .    - 0.5d0 * (state(2)**2 + state(3)**2)/state(1) )
          err(1) = abs(rhot - dnumsol(1))
          err(2) = abs(ut- dnumsol(2))
          err(3) = abs(vt- dnumsol(3))
          err(4) = abs(pt- dnumsol(4))

          abssol(1) = abs(rhot )
          abssol(2) = abs(ut)
          abssol(3) = abs(vt)
          abssol(4) = abs(pt)
          else
          err(1) = abs(state(1)- dnumsol(1))
          err(2) = abs(state(2)- dnumsol(2))
          err(3) = abs(state(3)- dnumsol(3))
          err(4) = abs(state(4)- dnumsol(4))

          abssol(1) = abs(state(1))
          abssol(2) = abs(state(2))
          abssol(3) = abs(state(3))
          abssol(4) = abs(state(4))
          endif


      dvol = dvol + wquad(iq)  * dx * dy
      do jj=1,4
       dl1error(jj) = dl1error(jj) + wquad(iq) * err(jj)  * dx * dy
       d_exact_mag(jj) = d_exact_mag(jj) + wquad(iq) * abssol(jj)*dx*dy

      if (dmaxerror(jj) < err(jj)) then
        dmaxerror(jj) = err(jj)
      endif

      if (d_exact_max_mag(jj) < abssol(jj) ) then
            d_exact_max_mag(jj) =  abssol(jj)
      endif
      end do




 11   continue

      goto 821 ! go to compute error now







 2    continue ! cut cell projection now lo
          val = 0.d0
          arr = ar(kirr)
          ivert = 1
          do 20 while (poly(ivert+1,1,kirr) .ne. -11.)
            ivert = ivert + 1
  20     continue


        xc = xcirr(kirr)
        yc = ycirr(kirr)
        call f(state, xc, yc, time, iprob)
        if(quad) then
             gamma1 = 0.4d0
             rhot = state(1)
             ut   = state(2)/state(1)
             vt   = state(3)/state(1)
             pt   = gamma1*( state(4)
     .    - 0.5d0 * (state(2)**2 + state(3)**2)/state(1) )
          err(1) = abs( rhot - q(i,j,1) )
          err(2) = abs( ut   - q(i,j,2) )
          err(3) = abs( vt   - q(i,j,3) )
          err(4) = abs( pt   - q(i,j,4) )
          else
          err(:) = abs(state(:)- q(i,j,:))
          endif
        dfv_vol_err = arr*err(1)
        err = 0.d0







          itri = ivert - 3
          idx1 = 1

          tempar = 0.d0
          do 21 it = 1, itri ! for each  triangle
            idx2 = it + 1
            idx3 = it + 2

            x1 = poly(idx1,1,kirr)
            y1 = poly(idx1,2,kirr)

            x2 = poly(idx2,1,kirr)
            y2 = poly(idx2,2,kirr)

            x3 = poly(idx3,1,kirr)
            y3 = poly(idx3,2,kirr)

            artri = area(x1, x2, x3, y1, y2, y3)




            do 22 itq = 1,3

                xc = x1 * rtri(itq) + x2 * stri(itq)
     .              +  x3 * (1.d0-rtri(itq)-stri(itq))
                yc = y1 * rtri(itq) + y2 * stri(itq)
     .              +  y3 * (1.d0-rtri(itq)-stri(itq))


                  call f(state, xc, yc, time, iprob)


       deltax = xc-xcen
       deltay = yc-ycen

       do mmm = 1,nvar
       call evalU(dnumsol,deltax, deltay, dx, dy,
     .            q, qx, qy,
     .            qxx, qxy, qyy,
     .            qxxx, qxxy, qxyy, qyyy,
     .            i,j, kirr, mitot, mjtot, nvar, mmm)
       end do


          if(quad) then
             gamma1 = 0.4d0
             rhot = state(1)
             ut   = state(2)/state(1)
             vt   = state(3)/state(1)
             pt   = gamma1*( state(4)
     .    - 0.5d0 * (state(2)**2 + state(3)**2)/state(1) )
          err(1) = abs(rhot - dnumsol(1))
          err(2) = abs(ut- dnumsol(2))
          err(3) = abs(vt- dnumsol(3))
          err(4) = abs(pt- dnumsol(4))

          abssol(1) = abs(rhot )
          abssol(2) = abs(ut)
          abssol(3) = abs(vt)
          abssol(4) = abs(pt)
          else
          err(1) = abs(state(1)- dnumsol(1))
          err(2) = abs(state(2)- dnumsol(2))
          err(3) = abs(state(3)- dnumsol(3))
          err(4) = abs(state(4)- dnumsol(4))

          abssol(1) = abs(state(1))
          abssol(2) = abs(state(2))
          abssol(3) = abs(state(3))
          abssol(4) = abs(state(4))
          endif

      dvol = dvol + artri*wtri(itq)
      dl1error = dl1error + artri*wtri(itq) * err
      d_exact_mag = d_exact_mag+ artri*wtri(itq) * abssol




      do jj=1,4

      if (dmaxerror(jj) < err(jj)) then
        dmaxerror(jj) = err(jj)
      endif

      if (d_exact_max_mag(jj) < abssol(jj) ) then
            d_exact_max_mag(jj) =  abssol(jj)
      endif

      end do




  22        continue ! for each quadrature point on each triangle

  21      continue ! for each triangle

      goto 821












 3    continue ! cut cell projection now ho
          val = 0.d0
          arr = ar(kirr)
          ivert = 1
          do 15 while (poly(ivert+1,1,kirr) .ne. -11.)
            ivert = ivert + 1
  15     continue




          itri = ivert - 3
          idx1 = 1

          tempar = 0.d0
          do 13 it = 1, itri-1 ! for each  triangle
            idx2 = it + 1
            idx3 = it + 2

            x1 = poly(idx1,1,kirr)
            y1 = poly(idx1,2,kirr)

            x2 = poly(idx2,1,kirr)
            y2 = poly(idx2,2,kirr)

            x3 = poly(idx3,1,kirr)
            y3 = poly(idx3,2,kirr)

            artri = area(x1, x2, x3, y1, y2, y3)




            do 14 itq = 1,3

                xc = x1 * rtri(itq) + x2 * stri(itq)
     .              +  x3 * (1.d0-rtri(itq)-stri(itq))
                yc = y1 * rtri(itq) + y2 * stri(itq)
     .              +  y3 * (1.d0-rtri(itq)-stri(itq))


                  call f(state, xc, yc, time, iprob)


       deltax = xc-xcirr_ho(kirr)
       deltay = yc-ycirr_ho(kirr)

       do mmm = 1,nvar
       call evalU(dnumsol,deltax, deltay, dx, dy,
     .            q, qx, qy,
     .            qxx, qxy, qyy,
     .            qxxx, qxxy, qxyy, qyyy,
     .            i,j, kirr, mitot, mjtot, nvar, mmm)
       end do


          if(quad) then
             gamma1 = 0.4d0
             rhot = state(1)
             ut   = state(2)/state(1)
             vt   = state(3)/state(1)
             pt   = gamma1*( state(4)
     .    - 0.5d0 * (state(2)**2 + state(3)**2)/state(1) )
          err(1) = abs(rhot - dnumsol(1))
          err(2) = abs(ut- dnumsol(2))
          err(3) = abs(vt- dnumsol(3))
          err(4) = abs(pt- dnumsol(4))

          abssol(1) = abs(rhot )
          abssol(2) = abs(ut)
          abssol(3) = abs(vt)
          abssol(4) = abs(pt)
          else
          err(1) = abs(state(1)- dnumsol(1))
          err(2) = abs(state(2)- dnumsol(2))
          err(3) = abs(state(3)- dnumsol(3))
          err(4) = abs(state(4)- dnumsol(4))

          abssol(1) = abs(state(1))
          abssol(2) = abs(state(2))
          abssol(3) = abs(state(3))
          abssol(4) = abs(state(4))
          endif

      dvol = dvol + artri*wtri(itq)
      dl1error = dl1error + artri*wtri(itq) * err
      d_exact_mag = d_exact_mag+ artri*wtri(itq) * abssol




      do jj=1,4

      if (dmaxerror(jj) < err(jj)) then
        dmaxerror(jj) = err(jj)
      endif

      if (d_exact_max_mag(jj) < abssol(jj) ) then
            d_exact_max_mag(jj) =  abssol(jj)
      endif

      end do




  14        continue ! for each quadrature point on each triangle

  13      continue ! for each triangle


      if(ssw .eq. 2 .or. ssw .eq. -2) then
      x1 = poly(ivert-2,1,kirr)
      y1 = poly(ivert-2,2,kirr)
      x2 = bdry(1,1,kirr)
      y2 = bdry(1,2,kirr)
      x3 = bdry(3,1,kirr)
      y3 = bdry(3,2,kirr)
      x4 = bdry(2,1,kirr)
      y4 = bdry(2,2,kirr)

      do nq = 1,ntriquad_ho
      artri = dJ2(x1,x2,x3,x4,y1,y2,y3,y4,nq)
      xc = rs2xy_2(x1,x2,x3,x4,nq)
      yc = rs2xy_2(y1,y2,y3,y4,nq)
      call f(state, xc, yc, time, iprob)


       deltax = xc-xcirr_ho(kirr)
       deltay = yc-ycirr_ho(kirr)

       do mmm = 1,nvar
       call evalU(dnumsol,deltax, deltay, dx, dy,
     .            q, qx, qy,
     .            qxx, qxy, qyy,
     .            qxxx, qxxy, qxyy, qyyy,
     .            i,j, kirr, mitot, mjtot, nvar, mmm)
       end do


          if(quad) then
             gamma1 = 0.4d0
             rhot = state(1)
             ut   = state(2)/state(1)
             vt   = state(3)/state(1)
             pt   = gamma1*( state(4)
     .    - 0.5d0 * (state(2)**2 + state(3)**2)/state(1) )
          err(1) = abs(rhot - dnumsol(1))
          err(2) = abs(ut- dnumsol(2))
          err(3) = abs(vt- dnumsol(3))
          err(4) = abs(pt- dnumsol(4))

          abssol(1) = abs(rhot )
          abssol(2) = abs(ut)
          abssol(3) = abs(vt)
          abssol(4) = abs(pt)
          else
          err(1) = abs(state(1)- dnumsol(1))
          err(2) = abs(state(2)- dnumsol(2))
          err(3) = abs(state(3)- dnumsol(3))
          err(4) = abs(state(4)- dnumsol(4))

          abssol(1) = abs(state(1))
          abssol(2) = abs(state(2))
          abssol(3) = abs(state(3))
          abssol(4) = abs(state(4))
          endif

      dvol = dvol + artri*wtri_ho(nq)
      dl1error = dl1error + artri*wtri_ho(nq) * err
      d_exact_mag = d_exact_mag+ artri*wtri_ho(nq) * abssol




      do jj=1,4

      if (dmaxerror(jj) < err(jj)) then
        dmaxerror(jj) = err(jj)
      endif

      if (d_exact_max_mag(jj) < abssol(jj) ) then
            d_exact_max_mag(jj) =  abssol(jj)
      endif

      end do

      enddo



      elseif(ssw .eq. 3) then
      x1 = poly(ivert-2,1,kirr)
      y1 = poly(ivert-2,2,kirr)
      x2 = bdry(1,1,kirr)
      y2 = bdry(1,2,kirr)
      x3 = bdry(4,1,kirr)
      y3 = bdry(4,2,kirr)
      x4 = bdry(2,1,kirr)
      y4 = bdry(2,2,kirr)
      x5 = bdry(3,1,kirr)
      y5 = bdry(3,2,kirr)

!        du = 0.d0
      do nq = 1,ntriquad_ho
      artri = dJ3(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,nq)
      xc = rs2xy_3(x1,x2,x3,x4,x5,nq)
      yc = rs2xy_3(y1,y2,y3,y4,y5,nq)

      call f(state, xc, yc, time, iprob)


       deltax = xc-xcirr_ho(kirr)
       deltay = yc-ycirr_ho(kirr)

       do mmm = 1,nvar
       call evalU(dnumsol,deltax, deltay, dx, dy,
     .            q, qx, qy,
     .            qxx, qxy, qyy,
     .            qxxx, qxxy, qxyy, qyyy,
     .            i,j, kirr, mitot, mjtot, nvar, mmm)
       end do


          if(quad) then
             gamma1 = 0.4d0
             rhot = state(1)
             ut   = state(2)/state(1)
             vt   = state(3)/state(1)
             pt   = gamma1*( state(4)
     .    - 0.5d0 * (state(2)**2 + state(3)**2)/state(1) )
          err(1) = abs(rhot - dnumsol(1))
          err(2) = abs(ut- dnumsol(2))
          err(3) = abs(vt- dnumsol(3))
          err(4) = abs(pt- dnumsol(4))

          abssol(1) = abs(rhot )
          abssol(2) = abs(ut)
          abssol(3) = abs(vt)
          abssol(4) = abs(pt)
          else
          err(1) = abs(state(1)- dnumsol(1))
          err(2) = abs(state(2)- dnumsol(2))
          err(3) = abs(state(3)- dnumsol(3))
          err(4) = abs(state(4)- dnumsol(4))

          abssol(1) = abs(state(1))
          abssol(2) = abs(state(2))
          abssol(3) = abs(state(3))
          abssol(4) = abs(state(4))
          endif

      dvol = dvol + artri*wtri_ho(nq)
      dl1error = dl1error + artri*wtri_ho(nq) * err
      d_exact_mag = d_exact_mag+ artri*wtri_ho(nq) * abssol




      do jj=1,4

      if (dmaxerror(jj) < err(jj)) then
        dmaxerror(jj) = err(jj)
      endif

      if (d_exact_max_mag(jj) < abssol(jj) ) then
            d_exact_max_mag(jj) =  abssol(jj)
      endif

      enddo

      enddo
      endif




 821  continue
      d_error = dl1error(ivar)
      d_maxerror = dmaxerror(ivar)
      d_mag = d_exact_mag(ivar)
      d_maxmag = d_exact_max_mag(ivar)

      return
      end









































      subroutine computed_boundary_lo(i,j,dx,dy,xlow,ylow,time, irr,
     .          q, qx, qy, qxx, qxy, qyy,
     .          qxxx, qxxy, qxyy, qyyy,
     .          mitot, mjtot, nvar,
     .          ivar, iprob2, lstgrd, lwidth,
     .          dl1error_out, dmaxerror_out, dl1mag_out, dmaxmag_out,
     .          dlen,
     .          dfv_err)
      implicit double precision (a-h,o-z)
      include "./cirr.i"
      include "./quadrature.i"
      common /order2/ ssw, quad, nolimiter
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .       ismp,gradThreshold, ilts,ibc,id,id2,res(ijsize,ijsize,4)
      logical quad
      dimension Uout(4), state(4), err(4),abssol(4)
      dimension dl1error(4), dmaxerror(4), dl1mag(4), dmaxmag(4)
      dimension dnumsol(4)
      dimension irr(mitot,mjtot)
      dimension q(mitot,mjtot, nvar)
      dimension qx(mitot,mjtot, nvar),qy(mitot,mjtot, nvar)
      dimension qxx(mitot,mjtot, nvar),qxy(mitot,mjtot, nvar)
      dimension qyy(mitot,mjtot, nvar)
      dimension qxxx(mitot,mjtot, nvar),qxxy(mitot,mjtot, nvar)
      dimension qxyy(mitot,mjtot, nvar),qyyy(mitot,mjtot, nvar)



         dl1error = 0.d0
         dmaxerror = -1.d0
         dl1mag = 0.d0
         dmaxmag = -1.d0
         dlen = 0.d0

         k = irr(i,j)
         do 20 kside=1,6
            if (poly(kside+2,1,k).eq.-11.) then
               x1 = poly(kside,1,k)
               y1 = poly(kside,2,k)
               x2 = poly(kside+1,1,k)
               y2 = poly(kside+1,2,k)
               go to 25
            endif
 20      continue
 25      continue
c
c        # boundary segment face:
c        ------------------------
c
         ix0 = ix(k)
         iy0 = iy(k)
c     # compute (alf,beta) = unit normal to boundary pointing out.
         hsx1 = x2
         hsx2 = x1
         hsy1 = y2
         hsy2 = y1
         rlen = dsqrt((hsy1-hsy2)**2 + (hsx1-hsx2)**2)
         alf = -(hsy1-hsy2)/rlen
         beta = -(hsx2-hsx1)/rlen
c
c     # compute boundary data  - q already has primitive variables
c



        x0 = xcirr(k)
        y0 = ycirr(k)

        call f(state, x0, y0, time, iprob2)
        if(quad) then
           gamma1 = 0.4d0
           rhot = state(1)
           ut   = state(2)/state(1)
           vt   = state(3)/state(1)
           pt   = gamma1*( state(4)
     .  - 0.5d0 * (state(2)**2.d0 + state(3)**2.d0)/state(1) )

          err(1) = abs(rhot - q(ix0, iy0, 1))
          err(2) = abs(ut   - q(ix0, iy0, 2))
          err(3) = abs(vt   - q(ix0, iy0, 3))
          err(4) = abs(pt   - q(ix0, iy0, 4))

          else
          err(1) = abs(state(1)- q(ix0, iy0, 1))
          err(2) = abs(state(2)- q(ix0, iy0, 2))
          err(3) = abs(state(3)- q(ix0, iy0, 3))
          err(4) = abs(state(4)- q(ix0, iy0, 4))
          endif

          dfv_err = err(1) * rlen

          err = 0.d0

      do 33 nn = 1,nlinequad
            bxpt = 0.5d0*(1+rline(nn))*hsx1+0.5d0*(1-rline(nn))*hsx2
            bypt = 0.5d0*(1+rline(nn))*hsy1+0.5d0*(1-rline(nn))*hsy2


         xdif = (bxpt-x0)
         ydif = (bypt-y0)

        do m = 1,4
        call evalU(dnumsol,xdif, ydif, dx, dy,
     .                 q, qx, qy,
     .                 qxx, qxy, qyy,
     .                 qxxx, qxxy, qxyy, qyyy,
     .                 ix0, iy0, k, mitot, mjtot, 4, m)
        end do

          call f(state, bxpt, bypt, time, iprob2)

          if(quad) then
             gamma1 = 0.4d0
             rhot = state(1)
             ut   = state(2)/state(1)
             vt   = state(3)/state(1)
             pt   = gamma1*( state(4)
     .    - 0.5d0 * (state(2)**2.d0 + state(3)**2.d0)/state(1) )

          err(1) = abs(rhot - dnumsol(1))
          err(2) = abs(ut   - dnumsol(2))
          err(3) = abs(vt   - dnumsol(3))
          err(4) = abs(pt   - dnumsol(4))

          abssol(1) = abs(rhot )
          abssol(2) = abs(ut   )
          abssol(3) = abs(vt   )
          abssol(4) = abs(pt   )
          else
          err(1) = abs(state(1)- dnumsol(1))
          err(2) = abs(state(2)- dnumsol(2))
          err(3) = abs(state(3)- dnumsol(3))
          err(4) = abs(state(4)- dnumsol(4))

          abssol(1) = abs(state(1) )
          abssol(2) = abs(state(2) )
          abssol(3) = abs(state(3) )
          abssol(4) = abs(state(4) )
          endif

          dl1error = dl1error + wline(nn)* err * rlen
          dl1mag = dl1mag + wline(nn)* abssol * rlen
          dlen = dlen + wline(nn) * rlen
          do jj = 1,4
            if (dmaxerror(jj) < err(jj)) then
                dmaxerror(jj) = err(jj)
            endif
            if (dmaxmag(jj) < abssol(jj)) then
                dmaxmag(jj) = abssol(jj)
            endif
          end do



 33   continue
      dl1error_out = dl1error(ivar)
      dmaxerror_out = dmaxerror(ivar)
      dl1mag_out = dl1mag(ivar)
      dmaxmag_out = dmaxmag(ivar)


      return
      end


      subroutine computed_boundary_ho(i,j,dx,dy,xlow,ylow,time, irr,
     .          q, qx, qy, qxx, qxy, qyy,
     .          qxxx, qxxy, qxyy, qyyy,
     .          mitot, mjtot, nvar,
     .          ivar, iprob2, lstgrd, lwidth,
     .          dl1error_out, dmaxerror_out, dl1mag_out, dmaxmag_out,
     .          dlen)
      implicit double precision (a-h,o-z)
      common /order2/ ssw, quad, nolimiter
      include "./cirr.i"
      include "./quadrature.i"
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .        ismp,gradThreshold, ilts,ibc,id,id2,res(ijsize,ijsize,4)
      logical quad
      dimension Uout(4), state(4), err(4),abssol(4)
      dimension dl1error(4), dmaxerror(4), dl1mag(4), dmaxmag(4)
      dimension dnumsol(4)
      dimension irr(mitot,mjtot)
      dimension q(mitot,mjtot, nvar)
      dimension qx(mitot,mjtot, nvar),qy(mitot,mjtot, nvar)
      dimension qxx(mitot,mjtot, nvar),qxy(mitot,mjtot, nvar)
      dimension qyy(mitot,mjtot, nvar)
      dimension qxxx(mitot,mjtot, nvar),qxxy(mitot,mjtot, nvar)
      dimension qxyy(mitot,mjtot, nvar),qyyy(mitot,mjtot, nvar)



         dl1error = 0.d0
         dmaxerror = -1.d0
         dl1mag = 0.d0
         dmaxmag = -1.d0
         dlen = 0.d0

         k = irr(i,j)
         do 20 kside=1,6
            if (poly(kside+2,1,k).eq.-11.) then
               x1 = poly(kside,1,k)
               y1 = poly(kside,2,k)
               x2 = poly(kside+1,1,k)
               y2 = poly(kside+1,2,k)
               go to 25
            endif
 20      continue
 25      continue
c
c        # boundary segment face:
c        ------------------------
c
         ix0 = ix(k)
         iy0 = iy(k)


       ivert = 1
       do 213 while (poly(ivert+1,1,k) .ne. -11.)
         ivert = ivert + 1
  213  continue

      x0 = xcirr_ho(k)
      y0 = ycirr_ho(k)

      if(ssw .eq. 2  .or. ssw .eq. -2) then
      x1 = poly(ivert-2,1,k)
      y1 = poly(ivert-2,2,k)

      x2 = bdry(1,1,k)
      y2 = bdry(1,2,k)

      x3 = bdry(3,1,k)
      y3 = bdry(3,2,k)

      x4 = bdry(2,1,k)
      y4 = bdry(2,2,k)
      elseif(ssw .eq. 3) then
      x1 = poly(ivert-2,1,k)
      y1 = poly(ivert-2,2,k)

      x2 = bdry(1,1,k)
      y2 = bdry(1,2,k)

      x3 = bdry(4,1,k)
      y3 = bdry(4,2,k)

      x6 = bdry(2,1,k)
      y6 = bdry(2,2,k)

      x7 = bdry(3,1,k)
      y7 = bdry(3,2,k)
      else
        print *,"problem in irregfluxgauss"
      endif

         ix0 = ix(k)
         iy0 = iy(k)

         x0   = xcirr_ho(k)
         y0   = ycirr_ho(k)

      do 33 nn = 1,nlinequad_ho
         r = (rline_ho(nn)+1.d0)/2.d0
         if(ssw .eq. 2  .or. ssw .eq. -2) then

         bxpt = x2*r*(-1.d0+2*r)
     .         +x3*(1.d0-r)*(1.d0-2.d0*r)
     .         +4.d0*x4*r*(1.d0-r)
         bypt = y2*r*(-1.d0+2*r)
     .         +y3*(1.d0-r)*(1.d0-2.d0*r)
     .         +4.d0*y4*r*(1.d0-r)
         Tx = x2*(-1.d0+2.d0*r)+2.d0*x2*r-x3*(1.d0-2.d0*r)
     .-2.d0*x3*(1.d0-r)+4.d0*x4*(1.d0-r)-4.d0*x4*r
         Ty = y2*(-1.d0+2.d0*r)+2.d0*y2*r-y3*(1.d0-2.d0*r)
     .-2.d0*y3*(1.d0-r)+4.d0*y4*(1.d0-r)-4.d0*y4*r
         dl = dsqrt(Tx**2 + Ty**2)
         dNx =  Ty/dl
         dNy = -Tx/dl
         elseif(ssw .eq. 3) then
      bxpt = (-dble(9 * r ** 2 * (1 - r)) + dble((-9 * (1 - r) ** 2 + 11
     . - 9 * r) * r) / 0.2D1) * x2 + (-0.9D1 / 0.2D1 * dble(r ** 2) * db
     .le(1 - r) + dble((-18 * (1 - r) ** 2 + 9 - 9 * r) * r) / 0.2D1 + 0
     ..1D1 - dble(r)) * x3 + (0.27D2 / 0.2D1 * dble(r ** 2) * dble(1 - r
     .) - 0.9D1 / 0.2D1 * dble(r) * dble(1 - r)) * x6 + dble((27 * (1 -
     .r) ** 2 - 9 + 9 * r) * r * x7) / 0.2D1
      bypt = (-dble(9 * r ** 2 * (1 - r)) + dble((-9 * (1 - r) ** 2 + 11
     . - 9 * r) * r) / 0.2D1) * y2 + (-0.9D1 / 0.2D1 * dble(r ** 2) * db
     .le(1 - r) + dble((-18 * (1 - r) ** 2 + 9 - 9 * r) * r) / 0.2D1 + 0
     ..1D1 - dble(r)) * y3 + (0.27D2 / 0.2D1 * dble(r ** 2) * dble(1 - r
     .) - 0.9D1 / 0.2D1 * dble(r) * dble(1 - r)) * y6 + dble((27 * (1 -
     .r) ** 2 - 9 + 9 * r) * r * y7) / 0.2D1

      Tx = (-dble(18 * r * (1 - r)) + dble(9 * r ** 2) + dble((9 - 18 *
     .r) * r) / 0.2D1 - 0.9D1 / 0.2D1 * dble((1 - r) ** 2) + 0.11D2 / 0.
     .2D1 - 0.9D1 / 0.2D1 * dble(r)) * x2 + (-dble(9 * r * (1 - r)) + 0.
     .9D1 / 0.2D1 * dble(r ** 2) + dble((27 - 36 * r) * r) / 0.2D1 - dbl
     .e(9 * (1 - r) ** 2) + 0.7D1 / 0.2D1 - 0.9D1 / 0.2D1 * dble(r)) * x
     .3 + (dble(27 * r * (1 - r)) - 0.27D2 / 0.2D1 * dble(r ** 2) - 0.9D
     .1 / 0.2D1 + dble(9 * r)) * x6 + dble((-45 + 54 * r) * r * x7) / 0.
     .2D1 + dble((27 * (1 - r) ** 2 - 9 + 9 * r) * x7) / 0.2D1

      Ty = (-dble(18 * r * (1 - r)) + dble(9 * r ** 2) + dble((9 - 18 *
     .r) * r) / 0.2D1 - 0.9D1 / 0.2D1 * dble((1 - r) ** 2) + 0.11D2 / 0.
     .2D1 - 0.9D1 / 0.2D1 * dble(r)) * y2 + (-dble(9 * r * (1 - r)) + 0.
     .9D1 / 0.2D1 * dble(r ** 2) + dble((27 - 36 * r) * r) / 0.2D1 - dbl
     .e(9 * (1 - r) ** 2) + 0.7D1 / 0.2D1 - 0.9D1 / 0.2D1 * dble(r)) * y
     .3 + (dble(27 * r * (1 - r)) - 0.27D2 / 0.2D1 * dble(r ** 2) - 0.9D
     .1 / 0.2D1 + dble(9 * r)) * y6 + dble((-45 + 54 * r) * r * y7) / 0.
     .2D1 + dble((27 * (1 - r) ** 2 - 9 + 9 * r) * y7) / 0.2D1

         dl = dsqrt(Tx**2 + Ty**2)
         dNx =  Ty/dl
         dNy = -Tx/dl

         else
         print *, "problem in irregflux_gauss"
         endif

         xdif = (bxpt-x0)
         ydif = (bypt-y0)

        do m = 1,4
        call evalU(dnumsol,xdif, ydif, dx, dy,
     .                 q, qx, qy,
     .                 qxx, qxy, qyy,
     .                 qxxx, qxxy, qxyy, qyyy,
     .                 ix0, iy0, k, mitot, mjtot, 4, m)
        end do

          call f(state, bxpt, bypt, time, iprob2)

          if(quad) then
             gamma1 = 0.4d0
             rhot = state(1)
             ut   = state(2)/state(1)
             vt   = state(3)/state(1)
             pt   = gamma1*( state(4)
     .    - 0.5d0 * (state(2)**2.d0 + state(3)**2.d0)/state(1) )

          err(1) = abs(rhot - dnumsol(1))
          err(2) = abs(ut   - dnumsol(2))
          err(3) = abs(vt   - dnumsol(3))
          err(4) = abs(pt   - dnumsol(4))

          abssol(1) = abs(rhot )
          abssol(2) = abs(ut   )
          abssol(3) = abs(vt   )
          abssol(4) = abs(pt   )
          else
          err(1) = abs(state(1)- dnumsol(1))
          err(2) = abs(state(2)- dnumsol(2))
          err(3) = abs(state(3)- dnumsol(3))
          err(4) = abs(state(4)- dnumsol(4))

          abssol(1) = abs(state(1) )
          abssol(2) = abs(state(2) )
          abssol(3) = abs(state(3) )
          abssol(4) = abs(state(4) )
          endif

          dl1error = dl1error + wline_ho(nn)* err * dl
          dl1mag = dl1mag + wline_ho(nn)* abssol * dl
          dlen = dlen + wline_ho(nn) * dl
          do jj = 1,4
            if (dmaxerror(jj) < err(jj)) then
                dmaxerror(jj) = err(jj)
            endif
            if (dmaxmag(jj) < abssol(jj)) then
                dmaxmag(jj) = abssol(jj)
            endif
          end do


  33  continue




      dl1error_out = dl1error(ivar)
      dmaxerror_out = dmaxerror(ivar)
      dl1mag_out = dl1mag(ivar)
      dmaxmag_out = dmaxmag(ivar)


      return
      end

