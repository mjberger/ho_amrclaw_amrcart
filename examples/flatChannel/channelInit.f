c
c-----------------------------------------------------------------------
c
      subroutine channelInit(xin,yin,state,time)

      implicit double precision(a-h,o-z)
      dimension state(4)
      include "cuserdt.i"
c
c     # For the point (x,y), compute the true solution (or initial
c     #  conditions) at this point 
c     # returns conserved variables
c

       !u = 0.d0
       !v = 0.d0
       u = 1.d0
       !v = 0.1d0
       v = 0.d0
       p = 1.d0

       x = xin - u*time
       y = yin - v*time

       !rho = 0.5d0
       !rho = y - 0.1d0*x + 0.5d0
       rho = x**2
       !rho = y**2  + 3.d0*x*y + 0.5*x + 1.d0
       !rho = y**3 + 2.d0*x*y + 0.5d0
       !rho = sin(x) + cos(y)

       state(1) = rho
       state(2) = rho*u
       state(3) = rho*v
       velsq = u*u + v*v
       state(4) = 0.5d0*rho*velsq + p/gamma1

       return
       end
