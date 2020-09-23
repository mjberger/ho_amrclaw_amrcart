c     
c     -------------------------------------------------------------
c     
      subroutine errout19(rect,maxip1,maxjp1,nvar,mptr,irr,
     1     mitot,mjtot,lwidth,lstgrd,
     2     dx,dy,xlow,ylow)
c     
      implicit double precision (a-h,o-z)
      include "cirr.i"
      dimension rect(maxip1,maxjp1,nvar),irr(mitot,mjtot)
      dimension err(20), soln(20), berr(20), bsoln(20)
      dimension errmax(20), solnmax(20)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .     ismp,gradThresholdx

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
      end do

      maxi = maxip1 - 1
      maxj = maxjp1 - 1
      xlowb = xlow - lwidth*dx
      ylowb = ylow - lwidth*dy
      write(6,100) mptr
 100  format(//,'L1 error/soln/rel.err for grid ',i5)
      do 15 i = 2, maxi
         do 15 j = 2, maxj
            if (iprob .eq. 19 .or. iprob .eq. 15) then
               iadj = i - 1 + lwidth
               jadj = j - 1 + lwidth
               kirr = irr(iadj,jadj)
               if (kirr .eq. lstgrd) then
                  call makep(poly(1,1,kirr),iadj,jadj,xlowb,ylowb,dx,dy)
                  xcen = xlow + (dfloat(i)-1.5)* dx
                  ycen = ylow + (dfloat(j)-1.5)* dy
               else if (kirr .ne. -1) then
                  xcen = xcirr(kirr)
                  ycen = ycirr(kirr)
               endif
               if (kirr .ne. -1) then
                  if (iprob .eq. 19) then
                    call p19tru(xcen,ycen,rho,u,v,p,poly(1,1,kirr),kirr)
                  elseif (iprob .eq. 15) then
                     call p15tru(xcen,ycen,rho,u,v,p,poly(1,1,kirr))
                  endif
               endif
            endif
            if (irr(i+lwidth-1,j+lwidth-1) .ne. -1) then
              pcomp =.4*(rect(i,j,4)-.5*(rect(i,j,2)**2+rect(i,j,3)**2)/
     .              rect(i,j,1))
              rhocomp = rect(i,j,1)
              err(1)  = err(1) + abs(rhocomp - rho)*ar(kirr)
              err(2)  = err(2) + abs(rect(i,j,2)/rhocomp - u)*ar(kirr)
              err(3)  = err(3) + abs(rect(i,j,3)/rhocomp - v)*ar(kirr)
              err(4)  = err(4) + abs(pcomp - p)*ar(kirr)
              soln(1) = soln(1) + abs(rho)  *ar(kirr)
              soln(2) = soln(2) + abs(u)  *ar(kirr)
              soln(3) = soln(3) + abs(v)*ar(kirr)
              soln(4) = soln(4) + abs(p)  *ar(kirr)
              if (irr(i+lwidth-1,j+lwidth-1) .ne. lstgrd) then ! add to bndry err calc
c     rl1d = sqrt(ar(kirr))
                  do kside=1,10
                     if (poly(kside+2,1,kirr) .eq. -11) then
                        x1 = poly(kside,1,kirr)
                        y1 = poly(kside,2,kirr)
                        x2 = poly(kside+1,1,kirr)
                        y2 = poly(kside+1,2,kirr)
                        go to 25
                     endif
                  end do
 25               rl1d = sqrt((x1-x2)**2 + (y1-y2)**2)  
                  berr(1)  = berr(1) + abs(rhocomp - rho)*rl1d
                  berr(2)  = berr(2) + abs(rect(i,j,2)/rhocomp - u)*rl1d
                  berr(3)  = berr(3) + abs(rect(i,j,3)/rhocomp - v)*rl1d
                  berr(4)  = berr(4) + abs(pcomp - p)*rl1d
                  bsoln(1) = bsoln(1) + abs(rho)*rl1d
                  bsoln(2) = bsoln(2) + abs(u)*rl1d
                  bsoln(3) = bsoln(3) + abs(v)*rl1d
                  bsoln(4) = bsoln(4) + abs(p)*rl1d
               endif
c     
               errmax(1) = max(errmax(1),abs(rhocomp-rho))
               errmax(2) = max(errmax(2),abs(rect(i,j,2)/rhocomp - u))
               errmax(3) = max(errmax(3),abs(rect(i,j,3)/rhocomp - v))
               errmax(4) = max(errmax(4),abs(pcomp - p))
               solnmax(1) = max(solnmax(1),abs(rho))
               solnmax(2) = max(solnmax(2),abs(u))
               solnmax(3) = max(solnmax(3),abs(v))
               solnmax(4) = max(solnmax(4),abs(p))

            endif
 15      continue
c     
         write(6,101) (err(m),soln(m),err(m)/soln(m),m=1,nvar)
 101     format("rho   ",3e15.7,/,"u vel ",3e15.7,/,"v vel ",3e15.7,/,
     .        "press ",3e15.7)
c     
         write(6,102)  mptr
 102     format(//,'L1 bndry. error/soln/rel.err for grid ',i5,
     .        ' wghted by bndry seg')
         write(6,101) (berr(m),bsoln(m),err(m)/bsoln(m),m=1,nvar)

         write(6,103) (errmax(m),solnmax(m),m=1,nvar)
 103     format(/,"max         error         soln",/,
     .        " rho   ",2e15.7,/,"  u    ",2e15.7,/,
     .        "  v    ",2e15.7,/," press ",2e15.7,//)
c     
         return
         end
