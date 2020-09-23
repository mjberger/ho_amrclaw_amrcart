c
c    
      subroutine readcart(val,mitot,mjtot)
      implicit double precision(a-h,o-z)
      dimension val(mitot,mjtot,4)
      dimension q(4)

10    continue
      read(8,100,end=99) i,j,(q(k),k=1,4)
 100  format(23x,i3,3x,i3,4x,4e16.7)
c     write(6,*)i,j,(q(k),k=1,4)
c     q is in primitive vars
      val(i,j,1) = q(1)
      val(i,j,2) = q(1)*q(2)
      val(i,j,3) = q(1)*q(3)
      val(i,j,4) = q(4)/.4 + .5d0*q(1)*(q(2)*q(2)+q(3)*q(3))


      go to 10 

99    return
      end
