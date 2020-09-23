c
c -------------------------------------------------------------
c
      subroutine outvar(rect,maxip1,maxjp1,nvar,mptr,irr,
     1                  mitot,mjtot,lwidth,lstgrd,
     2                  dx,dy,xlow,ylow,rectsm)
c
      implicit double precision (a-h,o-z)
      include "cirr.i"
      dimension rect(maxip1,maxjp1,nvar),irr(mitot,mjtot)
      real * 4  rectsm(mitot,mjtot,nvar)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      logical  eout
      data      spval/-1.e10/, eout/.true./
c
c  only output max - 1 rows and cols, since with cell centered
c  variables there is one extra cell outside the grid.
c
      maxi = maxip1 - 1
      maxj = maxjp1 - 1
      xlowb = xlow - lwidth*dx
      ylowb = ylow - lwidth*dy
      write(3,100) mptr
 100  format('*SOLN     ',i10,' is the grid - all variables')
      do 15 i = 2, maxi
      do 15 j = 2, maxj
      if (iprob .eq. 19 .or. iprob .eq. 15) then
	 iadj = i - 1 + lwidth
	 jadj = j - 1 + lwidth
	 kirr = irr(iadj,jadj)
	 if (kirr .eq. lstgrd) then
	    call makep(poly(1,1,kirr),iadj,jadj,xlowb,ylowb,dx,dy)
	    xcen = xlow + (dfloat(i)-1.5)* dx
	    ycen = ylow + (dfloat(j)-1.5)* dy
	 else if (kirr .ne. -1) then
	    xcen = xcirr(kirr)
	    ycen = ycirr(kirr)
	 endif
	 if (kirr .ne. -1) then
	   if (iprob .eq. 19) then
 	      call p19tru(xcen,ycen,rho,u,v,p,poly(1,1,kirr),kirr)
	    elseif (iprob .eq. 15) then
 	      call p15tru(xcen,ycen,rho,u,v,p,poly(1,1,kirr))
	   endif
	 endif
      endif
      do 16 ivar = 1, nvar
         if (irr(i+lwidth-1,j+lwidth-1) .eq. -1) then
c        if (irr(i+lwidth-1,j+lwidth-1) .ne. lstgrd) then
	    rectsm(i+lwidth-1,j+lwidth-1,ivar) = spval
         else
	    rectsm(i+lwidth-1,j+lwidth-1,ivar) = rect(i,j,ivar)
          endif
 16    continue
        if (eout .and. (iprob .eq. 19 .or. iprob .eq. 15)
     &     .and. kirr .ne. -1) then
 	  pcomp = .4*(rect(i,j,4)-.5*(rect(i,j,2)**2+rect(i,j,3)**2)/
     .           rect(i,j,1))
       	     rectsm(i+lwidth-1,j+lwidth-1,4) =  pcomp - p
             rectsm(i+lwidth-1,j+lwidth-1,1) = 
     &	          rect(i,j,1) - rho
             rectsm(i+lwidth-1,j+lwidth-1,2) = 
     &	          rect(i,j,2) - rho*u
             rectsm(i+lwidth-1,j+lwidth-1,3) = 
     &	          rect(i,j,3) - rho*v
        endif
 15    continue
c
      do 20 ivar = 1, nvar
         write(3,101) ((rectsm(i,j,ivar),i=lwidth+1,mitot-lwidth),
     .                                   j=lwidth+1,mjtot-lwidth)
 101     format(5e20.10)
 20   continue
c
      return
      end
