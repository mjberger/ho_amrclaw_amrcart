c
c-----------------------------------------------------------------------
c
       subroutine estdt(val,irr,mitot,mjtot,nvar,dx,dy,dtgrid,
     1                  lwidth,aux,naux)
c
       use amr_module, only : rinfinity
       implicit double precision (a-h, o-z)
       dimension val(nvar,mitot,mjtot), irr(mitot,mjtot)
       include "cuserdt.i"

c

       dtgrid = rinfinity ! initialize in case all solid cells
       speedmax =  0.d0
       do 20 j = lwidth+1, mjtot-lwidth
       do 10 i = lwidth+1, mitot-lwidth
          if (irr(i,j) .ne. -1) then
             rho  = val(1,i,j)
             u   = val(2,i,j)/rho
             v   = val(3,i,j)/rho
             etot  = val(4,i,j)
             velsq  = (u*u+v*v)
             p    = gamma1*(etot-.5d0*rho*velsq)
             c    = dsqrt(gamma*p/rho)
             speed = dsqrt(velsq) + c
             speedmax = max(speedmax,speed)
          endif
 10    continue
 20    continue
c     
!!effh = (dx*dy)/(2.d0*dx + 2.d0*dy)
          effh = min(dx,dy)
          if (speedmax .gt. 0.d0) then
             dtgrid = cflcart * effh / speedmax
          endif
c
c      write(*,*) "setting dt = ",dtgrid

       return
       end
