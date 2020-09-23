c
c -------------------------------------------------------------
c
c
       subroutine vprmtoc(q,mitot,mjtot,midim,mjdim,nvar)
       implicit double precision(a-h,o-z)
c      dimension q(mitot,mjtot,4)
       dimension q(midim,mjdim,nvar)
       common   /userdt/  gamma,dgamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
c
c NB: new dimensioning so can use on scrath arrays and actual soln arrays
c
       gamma1 = .4d0
       do 20 i = 1, mitot
       do 20 j = 1, mjtot
         rho  = q(i,j,1)
         u    = q(i,j,2)
         v    = q(i,j,3)
         pr   = q(i,j,4)
         ei   = pr/(rho*gamma1)
         en   = ei + .5d0*(u**2+v**2)
         q(i,j,2) = u *rho
         q(i,j,3) = v *rho
         q(i,j,4) = en*rho
 20       continue
c     
       return
       end
