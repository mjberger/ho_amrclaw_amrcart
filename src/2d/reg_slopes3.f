c
c------------------------------------------------------------
c
       subroutine reg_slopes3(qavg,q,qx,qy,qxx,qxy,qyy,mitot,mjtot,irr,
     &                        lstgrd,lwidth,hx,hy,xlow,ylow,mptr,
     &                        nvar,istage)

       use amr_module
       implicit double precision(a-h,o-z)

       dimension q(nvar,mitot,mjtot),qx(nvar,mitot,mjtot),
     &           qy(nvar,mitot,mjtot),irr(mitot,mjtot)
       dimension qavg(nvar,mitot,mjtot)
       dimension qxx(nvar,mitot,mjtot), qyy(nvar,mitot,mjtot)
       dimension qxy(nvar,mitot,mjtot)
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
c  now set for quadratic reconstruction. Limiting done elsewhere
c  initialized to zero in method (in case no slopes at all)
c
c  this version converts to pointwise and then does pointwise derivatives
c  formulas should be high order quadratics

       q = qavg ! this will become pointwise vals
       hx2 = 2.d0*hx
       hy2 = 2.d0*hy
       hxsq = hx*hx
       hysq = hy*hy
c
c      compute second derivatives first since needed to convert
c      note that division by h^2 postponed til next loop
       do j = 2, mjtot-1
       do i = 2, mitot-1
          qxx(:,i,j) = q(:,i+1,j)-2.d0*q(:,i,j)+q(:,i-1,j)
          qyy(:,i,j) = q(:,i,j+1)-2.d0*q(:,i,j)+q(:,i,j-1)
          ! and while were at it, one more
          qxy(:,i,j) =  ((q(:,i+1,j+1) - q(:,i-1,j+1))
     &                 - (q(:,i+1,j-1) - q(:,i-1,j-1)))/(hx2*hy2)
       end do
       end do


       do j = 2, mjtot-1
       do i = 2, mitot-1
          q(:,i,j)  = q(:,i,j) -(qxx(:,i,j)+qyy(:,i,j))/24.d0
          qxx(:,i,j) = qxx(:,i,j)/hxsq
          qyy(:,i,j) = qyy(:,i,j)/hysq
       end do
       end do

       ! finally do gradients
       do j = 2, mjtot-1
       do i = 2, mitot-1
          qx(:,i,j) = (q(:,i+1,j)-q(:,i-1,j))/hx2
          qy(:,i,j) = (q(:,i,j+1)-q(:,i,j-1))/hy2
       end do
       end do

c
       return
       end
