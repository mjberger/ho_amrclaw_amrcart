c
c-----------------------------------------------------------------------
c
      subroutine p15tru(xcen,ycen,rho,u,v,p,poly)
      implicit double precision(a-h,o-z)
      dimension poly(10,2)
      dimension q(4,3), qcent(5,4), area(5), celltot(4)
c
c     # For the centroid (x,y), compute the true cell avg.  To
c      compute cell average (to 3rd order). decompose
c      polygon into triangles, with centroid always
c      a vertex.  use approx. that cell avg. is
c      average of the function values on the 3 sides

       rhobase = 1.d0
       velbase = 2.d0
       pbase   = 1.d0
       alpha   = .11d0
       alpha   = .3333333333d0
       rmagu   = 1.d0/dsqrt(1.d0+alpha*alpha)
       rmagv   = alpha/dsqrt(1.d0+alpha*alpha)

       do 10 kside = 2, 7
c	 if (poly(kside,1) .eq. -1) then
	 if (poly(kside,1) .eq. -11) then
	    kend = kside - 2
	    go to 99
	 endif

	  xst = poly(kside-1,1)
	  yst = poly(kside-1,2)

	  xend = poly(kside,1)
	  yend = poly(kside,2)

          xavg = (xst+xend)/2.d0
          yavg = (yst+yend)/2.d0
          q(1,1) = rhobase + (yavg - alpha*xavg)**2
	  q(2,1) = velbase*rmagu
	  q(3,1) = velbase*rmagv
	  q(4,1) = pbase

          xavg = (xend+xcen)/2.d0
          yavg = (yend+ycen)/2.d0
 	  q(1,2) = rhobase + (yavg - alpha*xavg)**2
	  q(2,2) = velbase*rmagu
	  q(3,2) = velbase*rmagv
	  q(4,2) = pbase

          xavg = (xst+xcen)/2.d0
          yavg = (yst+ycen)/2.d0
 	  q(1,3) = rhobase + (yavg - alpha*xavg)**2
	  q(2,3) = velbase*rmagu
	  q(3,3) = velbase*rmagv
	  q(4,3) = pbase
c
c         compute area of triangle
          area(kside-1) = .5d0*((yst+ycen)*(xcen-xst) +
     &         (ycen+yend)*(xend-xcen)+(yend+yst)*(xst-xend))
	  do 11 m = 1, 4
	     qcent(kside-1,m) = (q(m,1)+q(m,2)+q(m,3))/3.d0
 11       continue
 10       continue

 99       continue

	  do 13 m = 1, 4
	     totar   = 0.d0
	     celltot(m)  = 0.d0
	     do 12 k = 1, kend
	       celltot(m) = area(k)*qcent(k,m) + celltot(m)
	       totar   = totar + area(k)
 12          continue
 13       continue

	  rho = celltot(1) / totar
	  u   = celltot(2) / totar
	  v   = celltot(3) / totar
	  p  = celltot(4) / totar

       return
       end
