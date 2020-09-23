c
c-----------------------------------------------------------------------
c
      subroutine p20tru(xcen,ycen,rho,poly,kirr,nvar,time)
      implicit double precision(a-h,o-z)
      dimension poly(10,2)
      dimension q(1,3), qcent(20,nvar), area(20), celltot(1)
c
c     # For the centroid (x,y), compute the true cell avg.  To
c      compute cell average (to 3rd order). decompose
c      polygon into triangles, with centroid always
c      a vertex.  use approx. that cell avg. is
c      average of the function values on the 3 sides

       do 10 kside = 2, 20
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
	  call p20fn(xavg,yavg,q(1,1),time)

          xavg = (xend+xcen)/2.d0
          yavg = (yend+ycen)/2.d0
	  call p20fn(xavg,yavg,q(1,2),time)

          xavg = (xst+xcen)/2.d0
          yavg = (yst+ycen)/2.d0
	  call p20fn(xavg,yavg,q(1,3),time)
c
c         compute area of triangle
          area(kside-1) = .5d0*((yst+ycen)*(xcen-xst) +
     &         (ycen+yend)*(xend-xcen)+(yend+yst)*(xst-xend))
	  do 11 m = 1,nvar 
	     qcent(kside-1,m) = (q(m,1)+q(m,2)+q(m,3))/3.d0
 11       continue
 10       continue

 99       continue

	  do 13 m = 1, nvar
	     totar   = 0.d0
	     celltot(m)  = 0.d0
	     do 12 k = 1, kend
	       celltot(m) = area(k)*qcent(k,m) + celltot(m)
	       totar   = totar + area(k)
 12          continue
 13       continue

	  rho = celltot(1) / totar

       return
       end
