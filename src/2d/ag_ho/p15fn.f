c
c-----------------------------------------------------------------------
c
      subroutine p15fn(xcen,ycen,rho,u,v,p,poly)
      implicit double precision(a-h,o-z)
      dimension poly(10,2)
c
c     # For the point (x,y), compute the true solution at this point 
c
       rhobase = 1.d0
       velbase = 2.d0
       pbase   = 1.d0
       alpha   = .11d0
       alpha   = .3333333333d0
       rmagu   = 1.d0/dsqrt(1.d0+alpha*alpha)
       rmagv   = alpha/dsqrt(1.d0+alpha*alpha)

       rho = rhobase + (ycen - alpha*xcen)**2
       u = velbase*rmagu
       v = velbase*rmagv
       p = pbase

       return
       end
