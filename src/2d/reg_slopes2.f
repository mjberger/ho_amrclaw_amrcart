c
c------------------------------------------------------------
c
       subroutine reg_slopes2(qp,qx,qy,qxx,qxy,qyy,mitot,mjtot,irr,
     &                   lstgrd,lwidth,hx,hy,xlow,ylow,mptr,nvar,istage)

       use amr_module
       implicit double precision(a-h,o-z)

       dimension qp(nvar,mitot,mjtot),qx(nvar,mitot,mjtot),
     &           qy(nvar,mitot,mjtot),irr(mitot,mjtot)
       dimension qxx(nvar,mitot,mjtot), qyy(nvar,mitot,mjtot)
       dimension qxy(nvar,mitot,mjtot)
       dimension qxxx(nvar,mitot,mjtot),qyyy(nvar,mitot,mjtot)
       dimension qxxy(nvar,mitot,mjtot),qxyy(nvar,mitot,mjtot)
       logical  quad, nolimiter
       include "cuserdt.i"
       common /order2/ ssw, quad, nolimiter


c
c  Set all slopes, even for irregular or solid.
c  Will be fixed in slope routines for irregular/solid cells called afterwards

c      # ssw = slope switch (1. for slopes, 0 for donor cell 0 slopes)
c      # now set in amrcart
c
c   ssw = 0 for ctu (no slopes), 1 for muscl (second order).
c   compute slopes using muscl limiter
c   q contains either conserved or primitive variables
c
c  Here do more accurate quadratic reconstruction. 
c  Construct cubic but only use quadratic part. Limiting done elsewhere
c  initialized in method (in case no slopes at all)

c
      do j = 3, mjtot-2
      do i = 3, mitot-2

      qyyy(:,i,j) = (2*qp(:,i,j-1) - qp(:,i,j-2) - 2*qp(:,i,j+1)
     .+ qp(:,i,j+2) + 2*qp(:,i-1,j-1) - qp(:,i-1,j-2) - 2*qp(:,i-1,j+1)
     .+ qp(:,i-1,j+2) + 2*qp(:,i-2,j-1) - qp(:,i-2,j-2)
     . - 2*qp(:,i-2,j+1) + qp(:,i-2,j+2) +  2*qp(:,i+1,j-1)
     . - qp(:,i+1,j-2) - 2*qp(:,i+1,j+1) + qp(:,i+1,j+2)
     .+  2*qp(:,i+2,j-1) - qp(:,i+2,j-2) - 2*qp(:,i+2,j+1)
     . + qp(:,i+2,j+2))/(10.*hy**3)

      qxyy(:,i,j) = (2*qp(:,i-1,j) + qp(:,i-1,j-1) - 2*qp(:,i-1,j-2)
     . + qp(:,i-1,j+1) - 2*qp(:,i-1,j+2) + 4*qp(:,i-2,j)
     .+ 2*qp(:,i-2,j-1) - 4*qp(:,i-2,j-2) + 2*qp(:,i-2,j+1)
     .- 4*qp(:,i-2,j+2) - 2*qp(:,i+1,j) - qp(:,i+1,j-1) +2*qp(:,i+1,j-2)
     .- qp(:,i+1,j+1) + 2*qp(:,i+1,j+2) - 4*qp(:,i+2,j)
     .-2*qp(:,i+2,j-1) + 4*qp(:,i+2,j-2) - 2*qp(:,i+2,j+1)
     .+ 4*qp(:,i+2,j+2))/(70.*hx*hy**2)


      qxxy(:,i,j) = (2*qp(:,i,j-1) + 4*qp(:,i,j-2) - 2*qp(:,i,j+1)
     .- 4*qp(:,i,j+2) +qp(:,i-1,j-1) + 2*qp(:,i-1,j-2) - qp(:,i-1,j+1)
     . - 2*qp(:,i-1,j+2) - 2*qp(:,i-2,j-1) - 4*qp(:,i-2,j-2)
     . + 2*qp(:,i-2,j+1) + 4*qp(:,i-2,j+2)+ qp(:,i+1,j-1)
     . + 2*qp(:,i+1,j-2) - qp(:,i+1,j+1)
     .- 2*qp(:,i+1,j+2) -2*qp(:,i+2,j-1) - 4*qp(:,i+2,j-2)
     . + 2*qp(:,i+2,j+1) + 4*qp(:,i+2,j+2))/(70.*hx**2*hy)


      qxxx(:,i,j) = (2*qp(:,i-1,j) + 2*qp(:,i-1,j-1) + 2*qp(:,i-1,j-2)
     . + 2*qp(:,i-1,j+1) + 2*qp(:,i-1,j+2) - qp(:,i-2,j) - qp(:,i-2,j-1)
     . - qp(:,i-2,j-2) - qp(:,i-2,j+1) - qp(:,i-2,j+2) - 2*qp(:,i+1,j)
     . - 2*qp(:,i+1,j-1) - 2*qp(:,i+1,j-2) - 2*qp(:,i+1,j+1)
     . - 2*qp(:,i+1,j+2) + qp(:,i+2,j) + qp(:,i+2,j-1) + qp(:,i+2,j-2)
     . + qp(:,i+2,j+1) +  qp(:,i+2,j+2))/(10.*hx**3)

      qyy(:,i,j) = (-350*qp(:,i,j) + 17*qp(:,i,j-1) + 68*qp(:,i,j-2)
     . + 17*qp(:,i,j+1) + 68*qp(:,i,j+2) - 10*qp(:,i-1,j)
     . + 7*qp(:,i-1,j-1) + 58*qp(:,i-1,j-2)  + 7*qp(:,i-1,j+1)
     . + 58*qp(:,i-1,j+2) - 40*qp(:,i-2,j) -  23*qp(:,i-2,j-1)
     . + 28*qp(:,i-2,j-2) - 23*qp(:,i-2,j+1) + 28*qp(:,i-2,j+2)
     . - 10*qp(:,i+1,j) + 7*qp(:,i+1,j-1) + 58*qp(:,i+1,j-2)
     . + 7*qp(:,i+1,j+1) + 58*qp(:,i+1,j+2) -  40*qp(:,i+2,j)
     . - 23*qp(:,i+2,j-1) + 28*qp(:,i+2,j-2) - 23*qp(:,i+2,j+1)
     . + 28*qp(:,i+2,j+2))/(945.*hy**2)

      qxy(:,i,j) = (qp(:,i-1,j-1) + 2*qp(:,i-1,j-2) - qp(:,i-1,j+1)
     . - 2*qp(:,i-1,j+2) + 2*qp(:,i-2,j-1) + 4*qp(:,i-2,j-2)
     . - 2*qp(:,i-2,j+1) - 4*qp(:,i-2,j+2)- qp(:,i+1,j-1)
     . - 2*qp(:,i+1,j-2) + qp(:,i+1,j+1) + 2*qp(:,i+1,j+2)
     . - 2*qp(:,i+2,j-1) - 4*qp(:,i+2,j-2) + 2*qp(:,i+2,j+1)
     . + 4*qp(:,i+2,j+2))/(100.*hx*hy)

      qxx(:,i,j) = (-350*qp(:,i,j) - 10*qp(:,i,j-1) - 40*qp(:,i,j-2)
     . - 10*qp(:,i,j+1) - 40*qp(:,i,j+2) + 17*qp(:,i-1,j)
     . + 7*qp(:,i-1,j-1) - 23*qp(:,i-1,j-2) + 7*qp(:,i-1,j+1)
     . - 23*qp(:,i-1,j+2) + 68*qp(:,i-2,j) + 58*qp(:,i-2,j-1)
     . + 28*qp(:,i-2,j-2) + 58*qp(:,i-2,j+1) + 28*qp(:,i-2,j+2)
     . + 17*qp(:,i+1,j) + 7*qp(:,i+1,j-1) - 23*qp(:,i+1,j-2)
     . + 7*qp(:,i+1,j+1) - 23*qp(:,i+1,j+2) + 68*qp(:,i+2,j)
     . + 58*qp(:,i+2,j-1) + 28*qp(:,i+2,j-2) + 58*qp(:,i+2,j+1)
     . + 28*qp(:,i+2,j+2))/(945.*hx**2)

         qy(:,i,j)  = (-288*qp(:,i,j-1) - 65*qp(:,i,j-2)
     . + 288*qp(:,i,j+1) + 65*qp(:,i,j+2) - 263*qp(:,i-1,j-1)
     .  - 15*qp(:,i-1,j-2) + 263*qp(:,i-1,j+1) + 15*qp(:,i-1,j+2)
     . - 188*qp(:,i-2,j-1) + 135*qp(:,i-2,j-2) + 188*qp(:,i-2,j+1)
     . - 135*qp(:,i-2,j+2) - 263*qp(:,i+1,j-1) - 15*qp(:,i+1,j-2)
     . + 263*qp(:,i+1,j+1) + 15*qp(:,i+1,j+2) - 188*qp(:,i+2,j-1)
     .  + 135*qp(:,i+2,j-2) + 188*qp(:,i+2,j+1)
     .  - 135*qp(:,i+2,j+2))/(1680.*hy)

         qx(:,i,j)  = (-288*qp(:,i-1,j) - 263*qp(:,i-1,j-1)
     . - 188*qp(:,i-1,j-2) - 263*qp(:,i-1,j+1) - 188*qp(:,i-1,j+2)
     . - 65*qp(:,i-2,j) - 15*qp(:,i-2,j-1) + 135*qp(:,i-2,j-2)
     . - 15*qp(:,i-2,j+1) + 135*qp(:,i-2,j+2) + 288*qp(:,i+1,j)
     . + 263*qp(:,i+1,j-1) + 188*qp(:,i+1,j-2) + 263*qp(:,i+1,j+1)
     . + 188*qp(:,i+1,j+2) + 65*qp(:,i+2,j) + 15*qp(:,i+2,j-1)
     . - 135*qp(:,i+2,j-2) + 15*qp(:,i+2,j+1)
     . - 135*qp(:,i+2,j+2))/(1680.*hx)

       end do
       end do

       return
       end
