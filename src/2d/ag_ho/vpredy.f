c
c -------------------------------------------------------------
c 
       subroutine vpredy(q,qy,ut,ub,mitot,mjtot,dt,dy,i,msize)
       implicit double precision(a-h,o-z)
       common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
       dimension q(mitot,mjtot,4)
       dimension qy(mitot,mjtot,4),ub(msize,4),ut(msize,4)
c
c predictor step for variables at midtime at cell edges for muscl
c        working on ith column
c
       hdtdy = .5d0*dt/dy
       do 80 j = 2, mjtot-2
	     rho = q(i,j,1)
	     u   = q(i,j,2)
	     v   = q(i,j,3)
	     pr  = q(i,j,4) 
	     cs2 = gamma*pr/rho
	     cs  = dsqrt(cs2)
	     bt1 = hdtdy*(qy(i,j,4)/cs-rho*qy(i,j,3))
	     bt2 = hdtdy*cs*(qy(i,j,1)-qy(i,j,4)/cs2)
	     bt3 = -1.d0*hdtdy*cs*qy(i,j,2)
	     if (v-cs .le. 0.d0) bt1 = 0.d0
	     if (v .le. 0.d0) then
		bt2 = 0.d0
		bt3 = 0.d0
	     endif
	     off = .5d0 - dmax1(v+cs,0.d0)*hdtdy
c
             ub(j+1,1) = bt1 + bt2 + rho + off*qy(i,j,1)
	     ub(j+1,2) = -bt3 + u + off*qy(i,j,2)
	     ub(j+1,3) = -bt1*cs/rho + v + off*qy(i,j,3)
	     ub(j+1,4) = bt1*cs2 + pr + off*qy(i,j,4)
c
	     rho = q(i,j+1,1)
	     u   = q(i,j+1,2)
	     v   = q(i,j+1,3)
	     pr =  q(i,j+1,4) 
	     cs2 = gamma*pr/rho
	     cs  = dsqrt(cs2)
	     bt2 =  hdtdy*cs*(qy(i,j+1,4)/cs2-qy(i,j+1,1))
	     bt3 =  hdtdy*cs*qy(i,j+1,2)
	     bt4 =  -1.d0*hdtdy*(qy(i,j+1,3)*rho+qy(i,j+1,4)/cs)
	     if (v .ge. 0.d0) then
		bt2 = 0.d0
		bt3 = 0.d0
	     endif
	     if (v+cs .ge. 0.d0) bt4 = 0.d0
	     off = .5d0 + dmin1(v-cs,0.d0)*hdtdy
c
	     ut(j+1,1) = bt2 + bt4 + rho - off*qy(i,j+1,1)
	     ut(j+1,2) = -bt3 + u - off*qy(i,j+1,2)
	     ut(j+1,3) = bt4*cs/rho + v - off*qy(i,j+1,3)
	     ut(j+1,4) = bt4*cs2 + pr - off*qy(i,j+1,4)
 80    continue
c
       return
       end
