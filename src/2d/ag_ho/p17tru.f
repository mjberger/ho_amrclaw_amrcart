c
c-----------------------------------------------------------------------
c
c
      subroutine p17tru(x,y,rho,u,v,pr)
      implicit double precision(a-h,o-z)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
c
c     # For the point (x,y), compute the true solution at this point 
c
       rhoinf = 1.d0
       pinf = 1.d0
       rm = .05d0
       vinf = rm * dsqrt(gamma*pinf/rhoinf)
       vinf2 = vinf * vinf
       einf = pinf / rhoinf / (gamma - 1.d0)
       ho = gamma * einf + vinf2 / 2.d0
c
c   a = cylinder radius 
c   centerx,centery =  the center of the cyl.
c      a = .1923538406d0
       a =  1.000567d0
       centerx = 2.d0
       centery = 2.d0

       a2 = a * a
       xx = x - centerx
       yy = y - centery

	y2 = yy * yy
	r2 = xx**2 + yy**2
	r4 = r2 * r2

	u = vinf * (1.d0 + a2 / r2 * ( 2.d0 * y2 / r2 - 1.d0) )
	v = -2.d0 * vinf * a2 * xx * yy / r4
	vel2 = u*u + v*v

	pr = pinf + .5d0 * rhoinf * (vinf2 - vel2)
	rho = rhoinf

      return
      end
