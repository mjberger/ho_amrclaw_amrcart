c
c-----------------------------------------------------------------------
c
       subroutine estdt(val,maxip1,maxjp1,nvar,dx,dy,dt,irr,
     1                  mitot,mjtot,lwidth)
c
       implicit double precision (a-h, o-z)
       dimension val(maxip1,maxjp1,nvar), irr(mitot,mjtot)
       common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                 ismp,gradThreshold
       common/steadydata/steadydiff

!       dt = 1.d+20
!
!       do 20 j = 2, maxjp1-1
!       do 10 i = 1, maxip1-1
!         if (irr(i-1+lwidth,j-1+lwidth) .ne. -1) then
!            if(iprob .ne. 25) then
!             rho  = val(i,j,1)
!             ux   = dabs(val(i,j,2))
!             uy   = dabs(val(i,j,3))
!             eng  = val(i,j,4)
!             usq  = (ux**2 + uy**2 ) /rho
!             eint = eng - .5*usq
!             p    = eint*(gamma-1.d0)
!             c    = dsqrt(gamma*p/rho)
!             dt1  = dx/(c + ux/rho)
!             dt2  = dy/(c + uy/rho)
!             dt   = dmin1(dt,dt1,dt2)
!             else
!             dt = dx*dy / (2*dx + 2*dy)
!             endif
!         endif
! 10    continue
! 20    continue
!c
!       dt = cfl * dt

       steadydiff = 1.d10
       dt = 0.d0

!      dt    = 0.147114317E-01
       return
       end
