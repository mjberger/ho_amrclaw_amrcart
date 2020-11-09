c
c------------------------------------------------------------
c
       subroutine reg_slopes(q,qx,qy,qxx,qxy,qyy,mitot,mjtot,irr,lstgrd,
     &                   lwidth,hx,hy,xlow,ylow,mptr,nvar)

       use amr_module
       implicit double precision(a-h,o-z)

       dimension q(nvar,mitot,mjtot),qx(nvar,mitot,mjtot),
     &           qy(nvar,mitot,mjtot),irr(mitot,mjtot)
       dimension qxx(nvar,mitot,mjtot), qyy(nvar,mitot,mjtot)
       dimension qxy(nvar,mitot,mjtot)
       logical  regular, quad, nolimiter
       include "cuserdt.i"
       common /order2/ ssw, quad, nolimiter

c      regular(i,j) = ((i.gt. lwidth).and.(i.le.mitot-lwidth).and.
c    &                 (j.gt. lwidth).and.(j.le.mjtot-lwidth))

c
c      # ssw = slope switch (1. for slopes, 0 for donor cell 0 slopes)
c      # now set in amrcart
c
c   ssw = 0 for ctu (no slopes), 1 for muscl (second order).
c   compute slopes using muscl limiter
c   q contains either conserved or primitive variables
c
c  now set for quadratic reconstruction. Limiting done elsewhere
c  initialized in method (in case no slopes at all)

       hx2 = 2.d0*hx
       hy2 = 2.d0*hy
       hxsq = hx*hx
       hysq = hy*hy
c
       if (ssw .eq. 0.d0) go to 99
c
       do 12 j = 1, mjtot
       do 10 i = 2, mitot-1
       !!  for now do for all if (irr(i,j) .ne. lstgrd) go to 10
       do 11 m = 1, nvar
          qx(m,i,j)  = (q(m,i+1,j) - q(m,i-1,j))/hx2
          qxx(m,i,j) = (q(m,i+1,j)-2.d0*q(m,i,j)+q(m,i-1,j))/hxsq
 11    continue
 10    continue
 12    continue
c
       do 22 j = 2, mjtot-1
       do 21 i = 1, mitot
       !! for now do for all if (irr(i,j) .ne. lstgrd) go to 21
       do 20 m = 1, nvar
          qy(m,i,j) = (q(m,i,j+1) - q(m,i,j-1))/hy2
          qyy(m,i,j) = (q(m,i,j+1)-2.d0*q(m,i,j)+q(m,i,j-1))/hysq
          ! next line relies on qxy=0 at first and last rows/cols, so
          ! qxy only exists in interior
          if (i .gt. 1 .and. i .lt. mitot) then
             qxy(m,i,j) = (qx(m,i,j+1) - qx(m,i,j-1))/hy2
          endif
 20    continue
 21    continue
 22    continue
c
 99    continue
       return
       end
