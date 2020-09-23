

c
c
c ---------------------------------------------------------------
c
      subroutine makepout(pout,i,j,xlow,ylow,dx,dy)
      implicit double precision (a-h,o-z)
      dimension pout(20,2)
c
c     # create a list of vertices of the poutgon corresponding to grid
c     # cell (i,j).
c
      pout(1,1) = xlow + (i-1)*dx
      pout(1,2) = ylow + (j-1)*dy
      pout(2,1) = xlow + (i-1)*dx
      pout(2,2) = ylow + j*dy
      pout(3,1) = xlow + i*dx
      pout(3,2) = ylow + j*dy
      pout(4,1) = xlow + i*dx
      pout(4,2) = ylow + (j-1)*dy
      pout(5,1) = xlow + (i-1)*dx
      pout(5,2) = ylow + (j-1)*dy
c     pout(6,1) = -1
      pout(6,1) = -11
      return
      end
