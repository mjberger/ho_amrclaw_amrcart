c
c ---------------------------------------------------------------------
c
       subroutine qslopes(qp,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,
     &                    hx,hy,xlow,ylow,qxx,qxy,qyy,mptr,nvar)
      implicit double precision(a-h,o-z)
      dimension qp(mitot,mjtot,nvar),qx(mitot,mjtot,nvar),
     &          qy(mitot,mjtot,nvar),irr(mitot,mjtot)
      dimension qxx(mitot,mjtot,nvar),qxy(mitot,mjtot,nvar)
      dimension qyy(mitot,mjtot,nvar)
      dimension recon(nvar)
      common /order2/ ssw, quad, nolimiter
      include "cirr.i"
      common   /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold



      return
      end
