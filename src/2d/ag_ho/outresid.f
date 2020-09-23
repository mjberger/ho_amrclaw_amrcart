c
c -------------------------------------------------------------
c
      subroutine outresid(res,mitot,mjtot,nvar,msize,irr,
     1                     lwidth,lstgrd,
     2                     dx,dy,xlowb,ylowb,nplot)
c
      implicit double precision (a-h,o-z)
      include "cirr.i"
      dimension irr(mitot,mjtot), res(mitot,mjtot,nvar)
      dimension valprim(4)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      common /order2/ ssw, quad, nolimiter
      logical  eout, quad, nolimiter
      logical pwconst
      data      spval/-1.e10/, eout/.true./
      character*23  filename 
c
c nplot == plotnum (2 for residual)
c
c  only output real rows and cols, no ghost cells 
c
c  will need to change for more than 1 grid for smallest h
      filename = 'disjointCutPlanesxx.dat'
      filename(18:18) = '0'
      filename(19:19) = char(ichar('0')+nplot)

      open(14,file=filename)
      write(14,100) 
 100  format('TITLE = "Extracted cutting Planes through mesh"' )
c

 8    continue

c  count needed for unstructured tec format (so dont have to look up new format)
      nCellsinPlane = 0  
c      do 10 i = lwidth+1, mitot-lwidth
c      do 10 j = lwidth+1, mjtot-lwidth
       do 10 i = lwidth-1, mitot-lwidth+2
       do 10 j = lwidth-1, mjtot-lwidth+2
         if (irr(i,j) .ne. -1) then
           nCellsinPlane = nCellsinPlane+1
         endif
 10   continue
c
      write(14,101) 4*nCellsinPlane,nCellsinPlane
 101  format('VARIABLES = x,y,Rho,U,V,Pressure,Xcent,Ycent',/,
     .       'Zone T="Cut",N =',i8,' E= ',i8,' F=FEPOINT')

c      do 15 i = lwidth+1, mitot-lwidth
c      do 15 j = lwidth+1, mjtot-lwidth
       do 15 i = lwidth-1, mitot-lwidth+2
       do 15 j = lwidth-1, mjtot-lwidth+2

         kirr = irr(i,j)
         if (kirr .eq. -1) go to 15

c        # test again for kirr -1 in case want to view in tecplot
        if (kirr .eq. lstgrd .or. kirr .eq. -1) then
            xcen = xlowb + (dfloat(i)-.5d0)* dx
            ycen = ylowb + (dfloat(j)-.5d0)* dy
        else 
           xcen = xcirr(kirr)
           ycen = ycirr(kirr)
        endif

         xcorner = xlowb + (dfloat(i)-1.d0)* dx    
         ycorner = ylowb + (dfloat(j)-1.d0)* dy
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
              valprim(ivar)=res(i,j,ivar)
           end do
             
           write(14,102) xc,yc,(valprim(ivar),ivar=1,4),xcen,ycen
 102       format(8e18.9)
        end do

 15   continue
c
c write mesh
c
      ico = 0
      do i = 1, nCellsinPlane
         write(14,104) ico+1,ico+2,ico+3,ico+4
 104     format(4i8)
         ico = ico + 4
      end do

      return
      end
