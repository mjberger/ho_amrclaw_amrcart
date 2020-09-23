! cholesky_d.f -*-f90-*-
! Using Cholesky decomposition, cholesky_d.f solve a linear equation Ax=b,
! where A is a n by n positive definite real symmetric matrix, x and b are
! real*8 vectors length n.
!
! Time-stamp: <2015-06-25 18:05:47 takeshi>
! Author: Takeshi NISHIMATSU
! Licence: GPLv3
!
! [1] A = G tG, where G is a lower triangular matrix and tG is transpose of G.
! [2] Solve  Gy=b with  forward elimination
! [3] Solve tGx=y with backward elimination
!
! Reference: Taketomo MITSUI: Solvers for linear equations [in Japanese]
!            http://www2.math.human.nagoya-u.ac.jp/~mitsui/syllabi/sis/info_math4_chap2.pdf
!
! Comment:   This Cholesky decomposition is used in src/elastic.F and
!            src/optimize-inho-strain.F of feram http://loto.sourceforge.net/feram/ .
!!
      subroutine solve(nn,n, G, b)
          implicit none
          integer,    intent(in)    :: n,nn
          real*8,     intent(in)    :: G(nn,nn)
          real*8,     intent(out)   :: b(nn)
          real*8 :: tmp
          integer    :: i,j

      ! [2]
      do i = 1, n
         tmp = 0.0d0
         do j = 1, i-1
            tmp = tmp + G(i,j)*b(j)
         end do
         b(i) = (b(i)-tmp)/G(i,i)
      end do

      ! [3]
      do i = n, 1, -1
         tmp = 0.0d0
         do j = i+1, n
            tmp = tmp + b(j)*G(j,i)
         end do
         b(i) = (b(i)-tmp)/G(i,i)
      end do


      end subroutine
