c
c ---------------------------------------------------------------
c
      subroutine setquadrature()

      implicit double precision (a-h,o-z)
      include "quadrature.i"
      common /userdt/ cflcart,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     &        ismp,gradThreshold,pwconst,ghost_ccg,limitTile,lpChoice
      common /order2/ ssw, quad, nolimiter


       !if (ismp .eq. 5 .and. .not. quad) then
       !  write(6,*)" cannot have ismp 5 without quadratic recon."
       !  write(6,*) "resetting quad to true"
       !  quad = .true.
       !endif
c
c 1 : linear reconstruction, first order slopes on cut cells
c second order slopes in the volume
c least squares linear reconstruction on cut cells (3x3 tile)
c 
c -1 : linear reconstruction, second order slopes on cut cells
c second order slopes in the volume
c least squares quadratic reconstruction on cut cells, chop off the quadratic degrees of freedom (3x3 tile)
c 
c -10 : linear reconstruction, "pointwise" second order slopes
c second order slopes in the volume
c least squares "pointwise" quadratic reconstruction on cut cells, chop off the quadratic degrees of freedom (3x3 tile)
c 
c 2 : quadratic reconstruction
c standard finite difference reconstruction on the volume
c least squares quadratic reconstruction on cut cells (3x3 tile)
c 
c -2 : quadratic reconstruction
c cubic finite difference reconstruction on the volume, chop off the cubic degrees of freedom (we need to talk about this)
c least squares cubic reconstruction on cut cells, chop off the cubic degrees of freedom (5x5 tile)
c 
c 3 : cubic reconstruction
c finite difference reconstruction on the volume (the formulas are really really complicated, we need to talk about these)
c least squares cubic reconstruction on cut cells  (5x5 tile)
c c 
      savessw = ssw
      ssw = 2 

      if (ssw .eq. 0) then
          ! 2*1 - 1 = 1
          rline = (/ 0.d0, 0.d0, 0.d0, 0.d0 /)
          wline = (/ 1.d0, 0.d0, 0.d0, 0.d0 /)
          nlinequad = 1

          ! degree 1 quadrature
          rtri = (/ 1.d0/3.d0, 0.d0, 0.d0, 0.d0/)
          stri = (/ 1.d0/3.d0, 0.d0, 0.d0, 0.d0/)
          wtri = (/ 1.d0     , 0.d0, 0.d0, 0.d0/)
          ntriquad = 1

          ! degree 1 quadrature
          rquad = (/ 0.d0, 0.d0, 0.d0, 0.d0 /)
          squad = (/ 0.d0, 0.d0, 0.d0, 0.d0 /)
          wquad = (/ 1.d0, 0.d0, 0.d0, 0.d0 /)
          nquadquad = 1
      elseif (ssw .eq. 1) then
          ! 2*1 - 1 = 1
          rline = (/ 0.d0, 0.d0, 0.d0, 0.d0 /)
          wline = (/ 1.d0, 0.d0, 0.d0, 0.d0 /)
          nlinequad = 1

          ! degree 1 quadrature
          rtri = (/ 1.d0/3.d0, 0.d0, 0.d0, 0.d0/)
          stri = (/ 1.d0/3.d0, 0.d0, 0.d0, 0.d0/)
          wtri = (/ 1.d0     , 0.d0, 0.d0, 0.d0/)
          ntriquad = 1

          ! degree 1 quadrature
          rquad = (/ 0.d0, 0.d0, 0.d0, 0.d0 /)
          squad = (/ 0.d0, 0.d0, 0.d0, 0.d0 /)
          wquad = (/ 1.d0, 0.d0, 0.d0, 0.d0 /)
          nquadquad = 1
      elseif (ssw .eq. -10) then
          ! 2*1 - 1 = 1
          rline = (/ 0.d0, 0.d0, 0.d0, 0.d0 /)
          wline = (/ 1.d0, 0.d0, 0.d0, 0.d0 /)
          nlinequad = 1

          ! degree 1 quadrature
          rtri = (/ 1.d0/3.d0, 0.d0, 0.d0, 0.d0/)
          stri = (/ 1.d0/3.d0, 0.d0, 0.d0, 0.d0/)
          wtri = (/ 1.d0     , 0.d0, 0.d0, 0.d0/)
          ntriquad = 1

          ! degree 1 quadrature
          rquad = (/ 0.d0, 0.d0, 0.d0, 0.d0 /)
          squad = (/ 0.d0, 0.d0, 0.d0, 0.d0 /)
          wquad = (/ 1.d0, 0.d0, 0.d0, 0.d0 /)
          nquadquad = 1
      elseif (ssw .eq. -1) then ! quadratic slopes
          ! 2*1 - 1 = 1
          rline = (/ 0.d0, 0.d0, 0.d0, 0.d0 /)
          wline = (/ 1.d0, 0.d0, 0.d0, 0.d0 /)
          nlinequad = 1

          ! degree 2 quadrature
          rtri(:) = (/ 1.d0/6.d0, 2.d0/3.d0,1.d0/6.d0, 0.d0 /)
          stri(:) = (/ 1.d0/6.d0, 1.d0/6.d0,2.d0/3.d0, 0.d0 /)
          wtri(:) = (/ 1.d0/3.d0, 1.d0/3.d0,1.d0/3.d0, 0.d0 /)
          ntriquad = 3


          ! degree 3 quadrature
          rquad = (/-dsqrt(1.d0/3.d0),  dsqrt(1.d0/3.d0),
     .         dsqrt(1.d0/3.d0), -dsqrt(1.d0/3.d0)/)
          squad = (/ -dsqrt(1.d0/3.d0), -dsqrt(1.d0/3.d0),
     .          dsqrt(1.d0/3.d0),  dsqrt(1.d0/3.d0) /)
          wquad = (/ 1.d0/4.d0, 1.d0/4.d0, 1.d0/4.d0, 1.d0/4.d0 /)
          nquadquad = 4
      elseif(ssw .eq. 2) then
          ! 2*2 - 1 = 3 degree
          rline(:) =(/-dsqrt(3.d0)/3.d0,dsqrt(3.d0)/3.d0, 0.d0, 0.d0 /)
          wline(:) = (/ 0.5d0, 0.5d0, 0.d0, 0.d0 /)
          nlinequad = 2

          ! degree 2 quadrature
          rtri(1:3) = (/ 1.d0/6.d0, 2.d0/3.d0,1.d0/6.d0 /)
          stri(1:3) = (/ 1.d0/6.d0, 1.d0/6.d0,2.d0/3.d0 /)
          wtri(1:3) = (/ 1.d0/3.d0, 1.d0/3.d0,1.d0/3.d0 /)
          ntriquad = 3

          ! degree 3 quadrature
          rquad = (/-dsqrt(1.d0/3.d0),  dsqrt(1.d0/3.d0),
     .         dsqrt(1.d0/3.d0), -dsqrt(1.d0/3.d0)/)
          squad = (/ -dsqrt(1.d0/3.d0), -dsqrt(1.d0/3.d0),
     .          dsqrt(1.d0/3.d0),  dsqrt(1.d0/3.d0) /)
          wquad = (/ 1.d0/4.d0, 1.d0/4.d0, 1.d0/4.d0, 1.d0/4.d0 /)
          nquadquad = 4

          !  quadratic reconstruction
          !+ quadratic mapping
          !+ linear jacobian = degree 5 quadrature
          rtri_ho(1:7) = (/ 0.333333333333333d0,0.797426985353087d0,
     .                    0.101286507323456d0,0.101286507323456d0,
     .                    0.470142064105115d0,0.059715871789770d0,
     .                    0.470142064105115d0 /)

          stri_ho(1:7) = (/0.333333333333333d0, 0.101286507323456d0,
     .                   0.797426985353087d0, 0.101286507323456d0,
     .                   0.059715871789770d0, 0.470142064105115d0,
     .                   0.470142064105115d0/)

      wtri_ho(1:7)=(/0.225000000000000d0/2.d0, 0.125939180544827d0/2.d0,
     .               0.125939180544827d0/2.d0, 0.125939180544827d0/2.d0,
     .               0.132394152788506d0/2.d0, 0.132394152788506d0/2.d0,
     .               0.132394152788506d0/2.d0/)
       ntriquad_ho = 7



       ! quadratic reconstruction
       !+quadratic mapping
       !+linear jacobian = degree 5
       rline_ho(1:3) = (/ dsqrt(3.d0/5.d0), 0.d0,-dsqrt(3.d0/5.d0) /)
       wline_ho(1:3) = (/ 5.d0/18.d0, 8.d0/18.d0,5.d0/18.d0        /)
       nlinequad_ho = 3

      elseif(ssw .eq. -2) then
          ! 2*2 - 1 = 3 degree
          rline(:) =(/-dsqrt(3.d0)/3.d0,dsqrt(3.d0)/3.d0, 0.d0, 0.d0 /)
          wline(:) = (/ 0.5d0, 0.5d0, 0.d0, 0.d0 /)
          nlinequad = 2

          ! degree 3 quadrature
          rtri = (/ 1.d0/3.d0, 3.d0/5.d0, 1.d0/5.d0, 1.d0/5.d0 /)
          stri = (/ 1.d0/3.d0, 1.d0/5.d0, 3.d0/5.d0, 1.d0/5.d0 /)
          wtri = (/-9.d0/16.d0,25.d0/48.d0,25.d0/48.d0,25.d0/48.d0/)
          ntriquad = 4


          ! degree 3 quadrature
          rquad = (/-dsqrt(1.d0/3.d0),  dsqrt(1.d0/3.d0),
     .         dsqrt(1.d0/3.d0), -dsqrt(1.d0/3.d0)/)
          squad = (/ -dsqrt(1.d0/3.d0), -dsqrt(1.d0/3.d0),
     .          dsqrt(1.d0/3.d0),  dsqrt(1.d0/3.d0) /)
          wquad = (/ 1.d0/4.d0, 1.d0/4.d0, 1.d0/4.d0, 1.d0/4.d0 /)
          nquadquad = 4


          ! cubic reconstruction + quadratic mapping + linear jacobian = degree 7 quadrature
          rtri_ho(1:13) = (/0.333333333333333d0,0.479308067841920d0,
     .                   0.260345966079040d0,0.260345966079040d0,
     .                   0.869739794195568d0,0.065130102902216d0,
     .                   0.065130102902216d0,0.048690315425316d0,
     .                   0.312865496004874d0,0.638444188569810d0,
     .                   0.048690315425316d0,0.312865496004874d0,
     .                   0.638444188569810d0 /)
          stri_ho(1:13) = (/0.333333333333333d0,0.260345966079040d0,
     .                   0.479308067841920d0,0.260345966079040d0,
     .                   0.065130102902216d0,0.869739794195568d0,
     .                   0.065130102902216d0,0.312865496004874d0,
     .                   0.048690315425316d0,0.048690315425316d0,
     .                   0.638444188569810d0,0.638444188569810d0,
     .                   0.312865496004874d0 /)

      wtri_ho(1:13) = (/-0.149570044467682d0/2.d0,
     .       0.175615257433208d0/2.d0,
     .       0.175615257433208d0/2.d0, 0.175615257433208d0/2.d0,
     .       0.053347235608838d0/2.d0, 0.053347235608838d0/2.d0,
     .       0.053347235608838d0/2.d0, 0.077113760890257d0/2.d0,
     .       0.077113760890257d0/2.d0, 0.077113760890257d0/2.d0,
     .       0.077113760890257d0/2.d0, 0.077113760890257d0/2.d0,
     .       0.077113760890257d0/2.d0 /)
          ntriquad_ho = 13

      ! cubic reconstruction + quadratic mapping + linear jacobian = degree 7 quadrature
          rline_ho(1:4) =(/ -dsqrt((3.d0+2.d0*dsqrt(6.d0/5.d0))/7.d0),
     .                    -dsqrt((3.d0-2.d0*dsqrt(6.d0/5.d0))/7.d0),
     .                     dsqrt((3.d0-2.d0*dsqrt(6.d0/5.d0))/7.d0),
     .                     dsqrt((3.d0+2.d0*dsqrt(6.d0/5.d0))/7.d0) /)
          wline_ho(1:4) = (/ (18.d0-dsqrt(30.d0))/36.d0/2.d0,
     .                      (18.d0+dsqrt(30.d0))/36.d0/2.d0,
     .                      (18.d0+dsqrt(30.d0))/36.d0/2.d0,
     .                      (18.d0-dsqrt(30.d0))/36.d0/2.d0 /)
          nlinequad_ho = 4

      elseif(ssw .eq. 3 .or. ssw .eq. -3) then
          ! 2*2 - 1 = 3 degree
          rline(:) =(/-dsqrt(3.d0)/3.d0,dsqrt(3.d0)/3.d0, 0.d0, 0.d0 /)
          wline(:) = (/ 0.5d0, 0.5d0, 0.d0, 0.d0 /)
          nlinequad = 2

          ! degree 3 quadrature
          rtri = (/ 1.d0/3.d0, 3.d0/5.d0, 1.d0/5.d0, 1.d0/5.d0 /)
          stri = (/ 1.d0/3.d0, 1.d0/5.d0, 3.d0/5.d0, 1.d0/5.d0 /)
          wtri = (/-9.d0/16.d0,25.d0/48.d0,25.d0/48.d0,25.d0/48.d0/)
          ntriquad = 4


          ! degree 3 quadrature
          rquad = (/-dsqrt(1.d0/3.d0),  dsqrt(1.d0/3.d0),
     .         dsqrt(1.d0/3.d0), -dsqrt(1.d0/3.d0)/)
          squad = (/ -dsqrt(1.d0/3.d0), -dsqrt(1.d0/3.d0),
     .          dsqrt(1.d0/3.d0),  dsqrt(1.d0/3.d0) /)
          wquad = (/ 1.d0/4.d0, 1.d0/4.d0, 1.d0/4.d0, 1.d0/4.d0 /)
          nquadquad = 4


          ! cubic reconstruction + cubic mapping + quadratic jacobian = degree 12 quadrature
          rtri_ho(:) = (/2.3565220452389998d-02, 4.8821738977380502d-01,
     .    4.8821738977380502d-01, 1.2055121541107899d-01,
     .    4.3972439229445998d-01, 4.3972439229445998d-01,
     .    4.5757922997576800d-01, 2.7121038501211597d-01,
     .    2.7121038501211597d-01, 7.4484770891682806d-01,
     .    1.2757614554158600d-01, 1.2757614554158600d-01,
     .    9.5736529909357904d-01, 2.1317350453210000d-02,
     .    2.1317350453210000d-02, 1.1534349453469800d-01,
     .    2.7571326968551402d-01, 6.0894323577978804d-01,
     .    2.7571326968551402d-01, 6.0894323577978804d-01,
     .    1.1534349453469800d-01, 2.2838332222257000d-02,
     .    2.8132558098993998d-01, 6.9583608678780295d-01,
     .    2.8132558098993998d-01, 6.9583608678780295d-01,
     .    2.2838332222257000d-02, 2.5734050548330001d-02,
     .    1.1625191590759699d-01, 8.5801403354407302d-01,
     .    1.1625191590759699d-01, 8.5801403354407302d-01,
     .    2.5734050548330001d-02 /)
      stri_ho(:) = (/4.8821738977380502d-01, 4.8821738977380502d-01,
     . 2.3565220452389998d-02, 4.3972439229445998d-01,
     . 4.3972439229445998d-01, 1.2055121541107899d-01,
     . 2.7121038501211597d-01, 2.7121038501211597d-01,
     . 4.5757922997576800d-01, 1.2757614554158600d-01,
     . 1.2757614554158600d-01, 7.4484770891682806d-01,
     . 2.1317350453210000d-02, 2.1317350453210000d-02,
     . 9.5736529909357904d-01, 2.7571326968551402d-01,
     . 6.0894323577978804d-01, 1.1534349453469800d-01,
     . 1.1534349453469800d-01, 2.7571326968551402d-01,
     . 6.0894323577978804d-01, 2.8132558098993998d-01,
     . 6.9583608678780295d-01, 2.2838332222257000d-02,
     . 2.2838332222257000d-02, 2.8132558098993998d-01,
     . 6.9583608678780295d-01, 1.1625191590759699d-01,
     . 8.5801403354407302d-01, 2.5734050548330001d-02,
     . 2.5734050548330001d-02, 1.1625191590759699d-01,
     . 8.5801403354407302d-01/)
      wtri_ho(:) = (/2.5731066440454999d-02/2.d0,
     .2.5731066440454999d-02/2.d0,
     . 2.5731066440454999d-02/2.d0, 4.3692544538038003d-02/2.d0,
     . 4.3692544538038003d-02/2.d0, 4.3692544538038003d-02/2.d0,
     . 6.2858224217885006d-02/2.d0, 6.2858224217885006d-02/2.d0,
     . 6.2858224217885006d-02/2.d0, 3.4796112930709000d-02/2.d0,
     . 3.4796112930709000d-02/2.d0, 3.4796112930709000d-02/2.d0,
     . 6.1662610515590003d-03/2.d0, 6.1662610515590003d-03/2.d0,
     . 6.1662610515590003d-03/2.d0, 4.0371557766381003d-02/2.d0,
     . 4.0371557766381003d-02/2.d0, 4.0371557766381003d-02/2.d0,
     . 4.0371557766381003d-02/2.d0, 4.0371557766381003d-02/2.d0,
     . 4.0371557766381003d-02/2.d0, 2.2356773202303001d-02/2.d0,
     . 2.2356773202303001d-02/2.d0, 2.2356773202303001d-02/2.d0,
     . 2.2356773202303001d-02/2.d0, 2.2356773202303001d-02/2.d0,
     . 2.2356773202303001d-02/2.d0, 1.7316231108659000d-02/2.d0,
     . 1.7316231108659000d-02/2.d0, 1.7316231108659000d-02/2.d0,
     . 1.7316231108659000d-02/2.d0, 1.7316231108659000d-02/2.d0,
     . 1.7316231108659000d-02/2.d0 /)
          ntriquad_ho = 33

      ! cubic reconstruction + cubic mapping + quadratic jacobian = degree 11 quadrature
          rline_ho(:) =(/ -9.324695142031521d-01,-6.612093864662646d-01,
     .    -2.386191860831969d-01, 2.386191860831969d-01,
     .    6.612093864662646d-01, 9.324695142031521d-01 /)
       wline_ho(:) = (/ 1.713244923791705d-01/2.d0,
     . 3.607615730481386d-01/2.d0,
     . 4.679139345726913d-01/2.d0,4.679139345726913d-01/2.d0,
     . 3.607615730481386d-01/2.d0,1.713244923791705d-01/2.d0  /)
       nlinequad_ho = 6


      endif


      ssw = savessw


      return
      end
