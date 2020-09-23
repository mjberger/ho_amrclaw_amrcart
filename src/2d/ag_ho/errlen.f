c
c ------------------------------------------------------------------
c
      subroutine errlen(f,g,irr,mitot,mjtot,dt,dx,dy,
     *                   lstgrd,mptr)
      implicit double precision (a-h,o-z)
      dimension f(mitot,mjtot,4),g(mitot,mjtot,4),irr(mitot,mjtot)
c
c     # irregular side modifications to fluxes for length only
c
      include "cirr.i"
c
c     # initialize fluxes at boundary to zero and set done array to false.
c     # done(k) is set to true when kth boundary cell is handled and is
c     # used to insure that interfaces between boundary cells
c     # are only handled once.
c
c
c     # modify fluxes on sides of irregular cells to reflect the fact that
c     # the length of each side may be less than the full dx or dy.
c     # (Recall that the fluxes f and g are the actual fluxes multiplied by
c     # the length of the side).
c
c     fix f first
      do 10 i = 2, mitot
      do 10 j = 2, mjtot
       if (irr(i,j) .eq. lstgrd .or. irr(i-1,j) .eq. lstgrd) then
	 f(i,j,1) = f(i,j,1)*dy
	 f(i,j,2) = f(i,j,2)*dy
	 f(i,j,3) = f(i,j,3)*dy
	 f(i,j,4) = f(i,j,4)*dy
       elseif (irr(i,j) .eq. -1 .or. irr(i-1,j) .eq. -1) then
	 f(i,j,1) = 0.d0
	 f(i,j,2) = 0.d0
	 f(i,j,3) = 0.d0
	 f(i,j,4) = 0.d0
c      else both of them are irregular. will recalc. in irreg3e
       endif
 10    continue
      do 20 i = 2, mitot
      do 20 j = 2, mjtot
       if (irr(i,j) .eq. lstgrd .or. irr(i,j-1) .eq. lstgrd) then
	 g(i,j,1) = g(i,j,1)*dx
	 g(i,j,2) = g(i,j,2)*dx
	 g(i,j,3) = g(i,j,3)*dx
	 g(i,j,4) = g(i,j,4)*dx
       elseif (irr(i,j) .eq. -1 .or. irr(i,j-1) .eq. -1) then
	 g(i,j,1) = 0.d0
	 g(i,j,2) = 0.d0
	 g(i,j,3) = 0.d0
	 g(i,j,4) = 0.d0
c      else both of them are irregular. will recalc. in irreg3e
       endif
 20    continue
c
	return
	end
