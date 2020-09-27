c
c
c     =====================================================
       subroutine qinit(meqn,mbc,mx,my,xcorn,ycorn,
     &                   dx,dy,q,maux,aux,lstgrd,irr)
c     =====================================================
c
c     # Set initial conditions for q.
c     # rotated channel problem
c
       use amr_module
       implicit double precision (a-h,o-z)
       include "quadrature.i"

       dimension q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
       dimension aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
       integer irr(1-mbc:mx+mbc, 1-mbc:my+mbc)
       dimension fakestate(4)
       data fakestate/0.5d0,0.d0,0.d0,2.5d0/

       xlow = xcorn - mbc*dx  ! makep needs this with 
       ylow = ycorn - mbc*dy  ! ghost cells
       q = 0.d0  !initialize for accumulation when computing integral
c
c      fill ghost cells too so can more easily plot before
c      time stepping starts 
       do 21 j = 1-mbc, my+mbc
       do 20 i = 1-mbc, mx+mbc

          kirr = irr(i,j)
          ! account for difference in declaration from 1..mitot or 1-mbc..mx
          call getCellAvgState(meqn,xlow,ylow,
     &                         dx,dy,q(:,i,j),lstgrd,kirr,i+mbc,j+mbc)

  20   continue
  21   continue

       return
       end
