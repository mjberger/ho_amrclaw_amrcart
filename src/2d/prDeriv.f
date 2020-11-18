c
c ----------------------------------------------------------
c
      subroutine prDeriv(qd,irr,mitot,mjtot,nvar,lstgrd)

      implicit real*8 (a-h, o-z)
      dimension qd(nvar,mitot,mjtot), irr(mitot,mjtot)

      do j = 1, mjtot
      do i = 1, mitot
         k = irr(i,j)
         if (k .eq. -1) cycle
         if (k .eq. lstgrd) then
            write(21,900) i,j,k,(qd(mm,i,j),mm=1,nvar)
 900        format(1x,3i5,4e15.7)
         else
            write(21,901) i,j,k,(qd(mm,i,j),mm=1,nvar)
 901        format('*',3i5,4e15.7)
         endif
      end do
      end do

      return
      end
