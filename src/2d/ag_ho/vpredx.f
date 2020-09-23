c
c -------------------------------------------------------------
c     
       subroutine vpredx(q,qx,ur,ul,mitot,mjtot,dt,dx,j,msize)
       implicit double precision(a-h,o-z)
       common   /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
       dimension q(mitot,mjtot,4)
       dimension qx(mitot,mjtot,4),ur(msize,4),ul(msize,4)
c
c predictor step for variables at midtime at cell edges for muscl
c working on row number j
c
       hdtdx = .5d0*dt/dx
	  do 20 i = 2, mitot-2
	     rho = q(i,j,1)
	     u   = q(i,j,2)
	     v   = q(i,j,3)
	     pr  = q(i,j,4) 
	     cs2 = gamma*pr / rho
	     cs  = dsqrt(cs2)
	     bt1 = hdtdx*(qx(i,j,4)/cs-rho*qx(i,j,2))
	     bt2 = hdtdx*cs*(qx(i,j,1)-qx(i,j,4)/cs2)
	     bt3 = hdtdx*cs*qx(i,j,3)
	     if (u-cs .le. 0.d0) bt1 = 0.d0
	     if (u .le. 0.d0) then
		bt2 = 0.d0
		bt3 = 0.d0
	     endif
	     off = .5d0 - dmax1(u+cs,0.d0)*hdtdx
c
	     ul(i+1,1) = bt1 + bt2 + rho + off*qx(i,j,1)
	     ul(i+1,2) = -bt1*cs/rho + u + off*qx(i,j,2)
	     ul(i+1,3) = bt3 + v + off*qx(i,j,3)
	     ul(i+1,4) = bt1*cs2 + pr + off*qx(i,j,4)
c
	     rho = q(i+1,j,1)
	     u   = q(i+1,j,2)
	     v   = q(i+1,j,3)
             pr =  q(i+1,j,4) 
	     cs2 = gamma*pr / rho
	     cs  = dsqrt(cs2)
	     bt2 = hdtdx*cs*(qx(i+1,j,4)/cs2-qx(i+1,j,1))
	     bt3 = -1.d0*hdtdx*cs*qx(i+1,j,3)
	     bt4 = -1.d0*hdtdx*(qx(i+1,j,2)*rho+qx(i+1,j,4)/cs)
	     if (u .gt. 0.d0) then
		bt2 = 0.d0
		bt3 = 0.d0
	     endif
	     if (u+cs .ge. 0.d0) bt4 = 0.d0
	     off = .5d0 + dmin1(u-cs,0.d0)*hdtdx
c
	     ur(i+1,1) = bt2 + bt4 + rho - off*qx(i+1,j,1)
	     ur(i+1,2) = bt4*cs/rho + u - off*qx(i+1,j,2)
	     ur(i+1,3) = bt3 + v - off*qx(i+1,j,3)
	     ur(i+1,4) = bt4*cs2 + pr - off*qx(i+1,j,4)
 20       continue
 40    continue
c
       return
       end
