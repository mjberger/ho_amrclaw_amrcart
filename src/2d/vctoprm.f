c
c -------------------------------------------------------------------
c
       subroutine vctoprm(q,qp,mitot,mjtot,nvar)
c
c      convert from cell averages in conserved variables to
c      pointwise values in primitive variables.  note for future
c      that may be overwriting conserved variables
c
       use amr_module
       implicit double precision(a-h,o-z)
       include "./quadrature.i"

       dimension q(nvar,mitot,mjtot),qp(nvar,mitot,mjtot)

       include "cuserdt.i"
c
c      given a pointwise value in conserved variables, convert to prim
c      for ex. if already have slopes through conserved vars that preserve
c      cell avgs, evaluate them at selected points
c      ### might be overwriting q with qp
c
       ! hxsq = hx*hx
       ! hysq = hy*hy

       do  j = 1, mjtot
       do  i = 1, mitot
          ! in 3rd order case, pointwise val differs by laplacian from cell average
           qp(1,i,j) = q(1,i,j)
c          qp(:,i,j) = q(:,i,j) - hxsq/24.d0*qxx(:,i,j) 
c    &                         - hysq/24.d0*qyy(:,i,j)
          qp(2,i,j) = q(2,i,j)/q(1,i,j)
          qp(3,i,j) = q(3,i,j)/q(1,i,j)
          velsq = qp(2,i,j)**2 + qp(3,i,j)**2
          qp(4,i,j) = gamma1*(q(4,i,j)-.5d0*qp(1,i,j)*velsq)

       end do
       end do
c
       return
       end
