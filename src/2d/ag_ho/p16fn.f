c
c-----------------------------------------------------------------------
c
      subroutine p16fn(xcen,ycen,q)
      implicit double precision(a-h,o-z)
      dimension q(4)
c
c     # For the point (x,y), compute the true solution at this point 
c     # for linear channel, tangential flow only, const pressure
c     # so this is really a linear advection test

       rho = ycen -.1*xcen + .5
       u   = .1
       v   = .01
       pr  = 1.

       q(1) = rho
       q(2) = u
       q(3) = v
       q(4) = pr

       return
       end
