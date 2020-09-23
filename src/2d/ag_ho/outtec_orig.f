c
c -------------------------------------------------------------
c
      subroutine outtec(q,nvar,mptr,irr,
     1                  mitot,mjtot,qx,qy,lwidth,lstgrd,
     2                  dx,dy,xlow,ylow,time,ncount,numHoods)
c
      implicit double precision (a-h,o-z)
      include "cirr.i"
      dimension irr(mitot,mjtot), q(mitot,mjtot,nvar) 
      dimension ncount(mitot,mjtot), numHoods(mitot,mjtot) 
      dimension qx(mitot,mjtot,nvar),qy(mitot,mjtot,nvar)
      dimension qxx(mitot,mjtot,nvar),qyy(mitot,mjtot,nvar)
      dimension qxy(mitot,mjtot,nvar)
      dimension valprim(4)
      dimension exactsoln(1)
      dimension exactval(mitot,mjtot,nvar)
      dimension dlimit(mitot,mjtot)
      dimension recon(nvar)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      common /order2/ ssw, quad, nolimiter
      logical  eout, quad, nolimiter
      logical pwconst

      data      spval/-1.e10/, eout/.true./
c      character*23  filename 
      common /dump/ uvel,vvel
      dimension uvel(849,848),vvel(848,849)
c
      dlimit(:,:) = 1.d0
      xlowb = xlow - lwidth*dx
      ylowb = ylow - lwidth*dy

      if (iprob .ne. 20)  call vctoprm(q,q,mitot,mjtot)
c
c     ### call for exterior bcs at each stage so can use slopes
c    ## NOTE THAT BNDRY CELLS FROM OTHER GRIDS NOT SET
            xhigh = xlowb + mitot*dx
            yhigh = ylowb + mjtot*dy
            call pphysbdlin(xlowb,xhigh,ylowb,yhigh,level,mitot,mjtot,
     &                   nvar,q,time,dx,dy,qx,qy,irr,lstgrd)
 
       do 5 ivar = 1, nvar
       do 5 i = 1, mitot
       do 5 j = 1, mjtot
         qx(i,j,ivar) = 0.d0
         qy(i,j,ivar) = 0.d0
         qxx(i,j,ivar) = 0.d0
         qxy(i,j,ivar) = 0.d0
         qyy(i,j,ivar) = 0.d0
 5     continue

c  set pwconst true for piecewise constant plots, set to false for slopes in tec output
c       pwconst =  .true.
      pwconst =  .false.
      if (pwconst) go to 8

      if (ssw .ne. 0.d0) then
         call slopes(q,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &               xlowb,ylowb,nvar)
         call qslopes(q,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &                 xlowb,ylowb,qxx,qxy,qyy,mptr,nvar)
      endif

 8    continue



c  count needed for unstructured tec format (so dont have to look up new format)
      nCellsinPlane = 0  
      do 10 i = lwidth+1, mitot-lwidth
      do 10 j = lwidth+1, mjtot-lwidth
         if (irr(i,j) .ne. -1) then
            nCellsinPlane = nCellsinPlane+1
         endif
 10   continue
c
      if (iprob .eq. 20 ) then
         write(14,101) 4*nCellsinPlane,nCellsinPlane
 101     format('VARIABLES = x,y,Rho,U,V,Pressure,Xcent,Ycent,Err',/,
     .          'Zone T="Cut",N =',i8,' E= ',i8,' F=FEPOINT')
      else
         write(14,103) 4*nCellsinPlane,nCellsinPlane
 103     format('VARIABLES = x,y,Rho,U,V,Pressure,Xcent,Ycent,',
     .                      'ncount,numHoods,i,j,err',/,
     .          'Zone T="Cut",N =',i10,' E= ',i10,' F=FEPOINT')
      endif

c  only output real rows and cols, no ghost cells 
c
      dmaxerror = 0.d0
      dl1error = 0.d0
      do 15 i = lwidth+1, mitot-lwidth
      do 15 j = lwidth+1, mjtot-lwidth
c       do 15 i = lwidth-1,mitot-lwidth+2
c       do 15 j = lwidth-1,mjtot-lwidth+2
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 15


         if(iprob .eq. 11 .and. kirr .ne. lstgrd) then

             do 20 kside=1,6
                if (poly(kside+2,1,kirr).eq.-11.) then
                   x1 = poly(kside,1,kirr)
                   y1 = poly(kside,2,kirr)
                   x2 = poly(kside+1,1,kirr)
                   y2 = poly(kside+1,2,kirr)
                   go to 25
                endif
 20          continue
 25          continue
              xcenter = 2.d0
              ycenter = 2.d0

              xm = 0.5d0*(x1+x2)
              ym = 0.5d0*(y1+y2)

              xcen = xcirr(kirr)
              ycen = ycirr(kirr)

              diffx = xm - xcen
              diffy = ym - ycen

              recon(:) = q(i,j,:) +(xm-xcen)*qx(i,j,:)+
     .                             (ym-ycen)*qy(i,j,:)

              val = atan2((ym-ycenter),(xm-xcenter))
              print *,val,recon
         endif







c        # test again for kirr -1 in case want to view in tecplot
         if (kirr .eq. lstgrd .or. kirr .eq. -1) then
            xcen = xlowb + (dfloat(i)-.5d0)* dx
            ycen = ylowb + (dfloat(j)-.5d0)* dy
         else 
            xcen = xcirr(kirr)
            ycen = ycirr(kirr)
         endif

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
c  reconstruct to corners to can contour through disjoint dataset
c
         do ivar = 1, nvar
            valprim(ivar) = q(i,j,ivar)+(xc-xcen)*qx(i,j,ivar)+
     .                                  (yc-ycen)*qy(i,j,ivar)
         end do


         if (iprob .eq. 20)then
            call p20fn(xcen,ycen,exactsoln,time)
             errprim =  q(i,j,1) - exactsoln(1)
         endif
         if(iprob .eq. 24 .or. iprob .eq. 25) then
            call raref(xlowb,xhigh,ylowb,yhigh,level,mitot,mjtot,
     &             nvar,exactval,time,dx,dy,qx,qy,irr,lstgrd, i, j)
            currerror = dabs(exactval(i,j,1) - q(i,j,1))
            errprim = currerror


            k = irr(i,j)
            if(k .eq. lstgrd) then
                dl1error = dl1error + currerror*dx*dy
            else
                dl1error = dl1error + currerror*ar(k)
            endif

            if(currerror > dmaxerror) dmaxerror = currerror
         endif
         if (quad) then   ! add terms. 

            shiftxx = poly(8,1,k)
            shiftxy = poly(9,1,k)
            shiftyy = poly(10,1,k)

            do ivar = 1, nvar
              valprim(ivar) = valprim(ivar)   
     .              + (.5d0*(xc-xcen)*(xc-xcen)-shiftxx)*qxx(i,j,ivar)
     .              + (.5d0*(yc-ycen)*(yc-ycen)-shiftyy)*qyy(i,j,ivar)
     .              + (     (xc-xcen)*(yc-ycen)-shiftxy)*qxy(i,j,ivar)
            end do
         endif

          if (iprob .eq. 20) then  ! output error and soln, nvar=1
             write(14,102) xc,yc,(valprim(ivar),ivar=1,nvar),
     &                  uvel(i,j),vvel(i,j),1.d0,xcen,ycen,errprim
          else if (iprob .eq. 16) then  ! only output soln
c for debugging output error instead
             rhoex = ycen - .1d0*xcen + .5d0
c            rhoex = 1.d0
c            uncomment next line to output error
c            stuffed into w field
c            valprim(3) = valprim(1) - rhoex
c            uncomment next line to output density
c            valprim(1) = valprim(1) 
c            valprim(2) = valprim(2)
c            valprim(3) = valprim(3)
c            valprim(4) = valprim(4)
             write(14,1022) xc,yc,(valprim(ivar),ivar=1,nvar),
     &                  xcen,ycen,ncount(i,j),numHoods(i,j),i,j
          else if (iprob .eq. 19) then  
            if (kirr .eq. lstgrd) then
               call makep(poly(1,1,kirr),i,j,xlowb,ylowb,dx,dy)
            endif
            call p19tru(xcen,ycen,rhot,ut,vt,pt,poly(1,1,kirr),kirr)
            valprim(1) = valprim(1) - rhot
            valprim(2) = valprim(2) - ut
            valprim(3) = valprim(3) - vt
            valprim(4) = valprim(4) - pt
             write(14,1022) xc,yc,(valprim(ivar),ivar=1,nvar),
     &                  xcen,ycen
          else
            write(14,1022) xc,yc,(valprim(ivar),ivar=1,nvar),
     &                  xcen,ycen,ncount(i,j),numHoods(i,j),i,j,
     &                  errprim
          endif
 102      format(9e18.9)
 1022      format(8e25.15,4i7,1e25.15)
        end do

 15    continue

c
c write mesh
c
       ico = 0
       do i = 1, nCellsinPlane
          write(14,104) ico+1,ico+2,ico+3,ico+4
 104      format(4i10)
          ico = ico + 4
       end do


      print *,"l1 error at time ", time, " is : ", dl1error
      print *,"max error at time ", time, " is : ", dmaxerror
      return
      end
