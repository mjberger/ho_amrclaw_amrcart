c
c ----------------------------------------------------------
c
      subroutine prem(valbig,val,nvar,maxip1,maxjp1,
     1          mitot,mjtot,lwidth)
c
      implicit double precision (a-h,o-z)
      dimension valbig(mitot,mjtot,nvar), val(maxip1,maxjp1,nvar)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      data pi/3.1415926535d0/
c
c copy solution and boundary patches into one extended solution array
c
c for now, initialize to freestream so even solid cells have real vals
c then copy in from real solution array
c
       do 10 ivar = 1, nvar
       do 10 j = 2, maxjp1-1
       do 10 i = 2, maxip1-1
          valbig(i+lwidth-1,j+lwidth-1,ivar) = val(i,j,ivar)
 10    continue
c
       return
       end
