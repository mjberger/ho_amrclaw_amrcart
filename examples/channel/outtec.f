c
c -------------------------------------------------------------
c
      subroutine outtec(q,nvar,mptr,irr,mitot,mjtot,
     1                  lstgrd,dx,dy,xlow,ylow,time,
     2                  numHoods,ibunit)
c
      use amr_module
      implicit double precision (a-h,o-z)

      dimension q(nvar,mitot,mjtot) 
      dimension state(nvar)
      ! use temporary array for qp to avoid converting back
      dimension qp(nvar,mitot,mjtot)  
      integer numHoods(mitot,mjtot) 
      integer irr(mitot,mjtot) 
      dimension qx(nvar,mitot,mjtot),qy(nvar,mitot,mjtot)
      dimension qxx(nvar,mitot,mjtot),qyy(nvar,mitot,mjtot)
      dimension qxy(nvar,mitot,mjtot)
      dimension valprim(4), errprim(4)
      dimension exactsoln(1)
      include "cuserdt.i"
      common /order2/ ssw, quad, nolimiter
      logical  quad, nolimiter
      logical ghostOut
      logical isAllSolid,checkIfAllSolid
      character ch
c
c     first count to see if any non solid cells on grid
c     tecplot doesnt like zeroes
      isAllSolid = checkIfAllSolid(irr,mitot,mjtot,nghost)
      call countCellType(irr,mitot,mjtot,nghost,numSolid,numCut,
     &                         numFull,lstgrd)
      if (isAllSolid) then
        if (numSolid .ne. (mitot-2*nghost)*(mjtot-2*nghost)) then
          write(*,*)"count off in outtec"
          stop
        else
          return
        endif
      endif


 8    xlowb = xlow - nghost*dx
      ylowb = ylow - nghost*dy

      ! output primitive variables, not conserved  LATER
      !call vctoprm(q,qp,mitot,mjtot,nvar)
      qp = q
c
       qx = 0.d0
       qy = 0.d0
       qxx = 0.d0
       qxy = 0.d0
       qyy = 0.d0

c  if want to output ghost cells too change this flag
       !ghostOut = .false.
       ghostOut = .true.
       if (ghostOut) then
          ist = 1
          iend = mitot
          jst = 1
          jend = mjtot
       else
          ist = nghost+1
          iend = mitot-nghost
          jst = nghost+1
          jend = mjtot-nghost
       endif


c  set pwconst true for piecewise constant plots, set to false for slopes in tec output
c     pwconst =  .true.
c     pwconst =  .false.
      if (pwconst) go to 9

      if (ssw .ne. 0.d0) then
        istage = 0 ! signifies called from outtec NOT USED ATM
        call qslopes(qp,qx,qy,qxx,qxy,qyy,mitot,mjtot,irr,lstgrd,
     &               nghost,dx,dy,xlowb,ylowb,mptr,nvar)
      endif

 9    continue

c  count needed for unstructured tec format (so dont have to look up new format)
      nCellsinPlane = 0  
       do 11 i = ist,iend
       do 10 j = jst,jend
         if (irr(i,j) .ne. -1) then
            nCellsinPlane = nCellsinPlane+1
         endif
 10   continue
 11   continue
c
      write(14,103) 4*nCellsinPlane,nCellsinPlane
      write(13,1033) 4*nCellsinPlane,nCellsinPlane
 103  format('VARIABLES = x,y,Rho,U,V,Pressure,Xcent,Ycent,',
     .                   'ncount,numHoods,i,j,k,volFrac,mptr',/,
     .          'Zone T="Cut",N =',i10,' E= ',i10,' F=FEPOINT')
1033  format('VARIABLES = x,y,ErrRho,ErrU,ErrV,ErrPressure,Xcent,',
     .                   'Ycent,ncount,numHoods,i,j,k,volFrac,mptr',/,
     .          'Zone T="Cut",N =',i10,' E= ',i10,' F=FEPOINT')


c     initialize for error computation
      volDenErrorL1   = 0.d0
      volExactDenL1   = 0.d0
      exactVol        = 0.d0


c  only output real rows and cols, no ghost cells 
c
      do 16 j = jst,jend
      do 15 i = ist, iend
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 15

c        # test again for kirr -1 in case want to view in tecplot
         if (kirr .eq. lstgrd .or. kirr .eq. -1) then
            xcen = xlowb + (dfloat(i)-.5d0)* dx
            ycen = ylowb + (dfloat(j)-.5d0)* dy
            ch = ' '
         else 
            xcen = xcirr(kirr)
            ycen = ycirr(kirr)
            ch = '+'
c           reconstruct to midpt of solid bndry and output to special file for cylinder case
            call dumpBndry(qp,qx,qy,qxx,qxy,qyy,irr,mitot,mjtot,i,j,
     &                     nvar,ibunit,mptr) 
            bLength = getBlength(irr,kirr,mitot,mjtot)
         endif

         volFrac = ar(kirr)/ar(lstgrd)

62       continue

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
c
c  reconstruct to corners so can contour through disjoint dataset
c
         do ivar = 1, nvar
            valprim(ivar) = qp(ivar,i,j)+(xc-xcen)*qx(ivar,i,j) +
     .                                   (yc-ycen)*qy(ivar,i,j)
         end do

         write(14,102) xc,yc,(valprim(ivar),ivar=1,nvar),
     &                 xcen,ycen,ncount(kirr),numHoods(i,j),i,j,
     &                 kirr,volFrac,mptr
         ! this computes and  outputs cell centered error for 2nd order scheme
         !call channelInit(xcen,ycen,state) 
         ! now does cell average
         call getCellAvgState(nvar,xlowb,ylowb,dx,dy,state,
     &                        lstgrd,kirr,i,j,time)
         !errprim(:) = qp(:,i,j) - state
         errprim(:) = q(:,i,j) - state ! altho for now in cons vars
         rhot = state(1)
         write(13,102) xc,yc,(errprim(ivar),ivar=1,nvar),
     &                  xcen,ycen,ncount(kirr),numHoods(i,j),i,j,
     &                  kirr,volFrac,mptr
 102    format(8e25.15,5i8,1e10.2,i5)
        end do

c       add this point to error computation, get exact solution at bndry too
        volDenErrorL1 = volDenErrorL1 + ar(kirr)*abs(errprim(1))
        volExactDenL1 = volExactDenL1 + ar(kirr)*rhot
        exactVol =  exactVol + ar(kirr)

 15    continue
 16    continue
c
c write mesh
c
       ico = 0
       do i = 1, nCellsinPlane
          write(14,104) ico+1,ico+2,ico+3,ico+4
          write(13,104) ico+1,ico+2,ico+3,ico+4
 104      format(4i10)
          ico = ico + 4
       end do

c output errors
      write(outunit,600) volDenErrorL1,volExactDenL1,
     .                   volDenErrorL1/volExactDenL1,exactVol
 600  format("L1 density volume error ", e15.7,/,
     .       "L1 density exact soln   ", e15.7,/,
     .       "L1 Relative density error",e15.7,/,
     .       "Computed volume          ",e15.7,//)

      return
      end
c
c -----------------------------------------------------------
c
      logical function checkIfAllSolid(irr,mitot,mjtot,nghost)

      dimension irr(mitot,mjtot)

      checkIfAllSolid = .true.

      do 6 j = nghost+1, mjtot-nghost
      do 5 i = nghost+1, mitot-nghost
         if (irr(i,j) .ne. -1) then
           checkIfAllSolid = .false.
           go to 99
         endif
 5    continue
 6    continue

 99   return 
      end
c
c -----------------------------------------------------------
c
      subroutine countCellType(irr,mitot,mjtot,nghost,numSolid,numCut,
     &                         numFull,lstgrd)
      
      dimension irr(mitot,mjtot)

      numSolid = 0
      numCut = 0
      numFull = 0

      do 11 j = nghost+1,mjtot-nghost
      do 10 i = nghost+1,mitot-nghost
         if (irr(i,j) .eq. -1) then
            numSolid = numSolid + 1
         else if (irr(i,j) .ne. lstgrd) then
            numCut = numCut + 1
         else
            numFull = numFull + 1
         endif
 10   continue
 11   continue

      return
      end
c
c -----------------------------------------------------------------
c
      double precision function getBlength(irr,k,mitot,mjtot)

      use amr_module
      implicit double precision (a-h,o-z)
      dimension irr(mitot,mjtot)

      do 20 kside=1,6
         if (poly(kside+2,1,k).eq.-11.) then
            x1 = poly(kside,1,k)
            y1 = poly(kside,2,k)
            x2 = poly(kside+1,1,k)
            y2 = poly(kside+1,2,k)
            go to 25
         endif
 20   continue
 25   continue
      rlen = dsqrt((y2-y1)**2+(x2-x1)**2)

      getBlength = rlen
      return
      end
