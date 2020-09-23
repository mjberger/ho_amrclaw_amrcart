c     
c     -------------------------------------------------------------
c     
      subroutine errout20(rect,maxip1,maxjp1,nvar,mptr,irr,
     1     mitot,mjtot,lwidth,lstgrd,
     2     dx,dy,xlow,ylow,time)
c     
      implicit double precision (a-h,o-z)
      include "cirr.i"
      dimension rect(maxip1,maxjp1,nvar),irr(mitot,mjtot)
      dimension err(20), soln(20), berr(20), bsoln(20)
      dimension binner(20), bouter(20), sinner(20), souter(20)
      dimension errmax(20), solnmax(20)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .     ismp,gradThresholdx

      data xshift/1.5d0/, yshift/1.5d0/

c     
c     compute error in final solution in rect. 
c     note different dimensions for rect and locirr - this version does not
c     use enlarged grid for all amr stuff
c     
      do ivar = 1, nvar
         err(ivar)   = 0.d0
         soln(ivar)  = 0.d0
         berr(ivar)  = 0.d0
         bsoln(ivar) = 0.d0
         errmax(ivar) = 0.d0
         binner(ivar) = 0.d0
         bouter(ivar) = 0.d0
         sinner(ivar) = 0.d0
         souter(ivar) = 0.d0
      end do
      imax = 1
      jmax = 1

      maxi = maxip1 - 1
      maxj = maxjp1 - 1
      xlowb = xlow - lwidth*dx
      ylowb = ylow - lwidth*dy
      write(6,100) mptr
 100  format(//,'L1 error/soln/rel.err for grid ',i5)
      do 15 i = 2, maxi
      do 15 j = 2, maxj
         iadj = i - 1 + lwidth
         jadj = j - 1 + lwidth
         kirr = irr(iadj,jadj)
         if (kirr .eq. lstgrd) then
            call makep(poly(1,1,kirr),iadj,jadj,xlowb,ylowb,dx,dy)
            xcen = xlow + (dfloat(i)-1.5d0)* dx
            ycen = ylow + (dfloat(j)-1.5d0)* dy
         else if (kirr .ne. -1) then
            xcen = xcirr(kirr)
            ycen = ycirr(kirr)
         endif
         if (kirr .ne. -1) then
c           call p20tru(xcen,ycen,rho,poly(1,1,kirr),kirr,nvar,time)
            call p20fn(xcen,ycen,rho,time)
         endif
         if (irr(iadj,jadj) .ne. -1) then
            rhocomp = rect(i,j,1)
            err(1)  = err(1) + abs(rhocomp - rho)*ar(kirr)
            soln(1) = soln(1) + abs(rho)  *ar(kirr)
            if (kirr .ne. lstgrd) then ! add to bndry err calc
c           rl1d = sqrt(ar(kirr))
            do kside=1,10
              if (poly(kside+2,1,kirr) .eq. -11) then
                 x1 = poly(kside,1,kirr)
                 y1 = poly(kside,2,kirr)
                 x2 = poly(kside+1,1,kirr)
                 y2 = poly(kside+1,2,kirr)
                 go to 25
              endif
           end do
 25        rl1d = sqrt((x1-x2)**2 + (y1-y2)**2)  
           berr(1)  = berr(1) + abs(rhocomp - rho)*rl1d
           bsoln(1) = bsoln(1) + abs(rho)*rl1d
c          also compute error separately for inner and outer radius cut cells
           rad = sqrt((x1-xshift)**2 + (y1-yshift)**2)  ! shouldnt be close. no need for care
           if (abs(rad-.75) .lt. abs(rad-1.25)) then
             binner(1) = binner(1) + abs(rhocomp - rho)*rl1d
             sinner(1) = sinner(1) + abs(rho)*rl1d
           else
             bouter(1) = bouter(1) + abs(rhocomp - rho)*rl1d
             souter(1) = souter(1) + abs(rho)*rl1d
           endif
         endif
c     
          if (errmax(1) .lt.  abs(rhocomp-rho)) then
             errmax(1) = abs(rhocomp-rho)
             imax = iadj
             jmax = jadj
          endif
          solnmax(1) = max(solnmax(1),abs(rho))
       endif
 15   continue
c     
       write(6,101) err(1),soln(1),err(1)/soln(1)
 101   format("rho: err ",e15.7," soln ",e15.7," rel err ",e15.7)
c     
       write(6,102)  mptr
 102   format(//,'L1 bndry. error/soln/rel.err for grid ',i5,
     .        ' wghted by bndry seg')
       write(6,101) berr(1),bsoln(1),berr(1)/bsoln(1)
       write(6,104) binner(1), bouter(1)
 104   format("  inner boundary/ outer boundary error",/,2e15.7)
       write(6,105) sinner(1),souter(1)
 105   format("  inner boundary / outer boundary soln",/,2e15.7)
       write(6,114) binner(1)/sinner(1), bouter(1)/souter(1)
 114   format("  inner boundary/ outer boundary rel error",/,2e15.7)

       write(6,103) errmax(1),imax,jmax,solnmax(1)
 103   format(/,"max error ",e15.7," at ",2i5," soln max ",e15.7)
       kmax = irr(imax,jmax)
       if (kmax .eq. lstgrd) then
           xmax = xlowb + (imax-.5d0)*dx
           ymax = ylowb + (jmax-.5d0)*dy
           write(6,*)" this is a full cell, located at ",xmax,ymax
       else
           write(6,*)" this is a cut cell located at",
     &                      xcirr(kmax),ycirr(kmax)
      endif
c     
      return
      end
