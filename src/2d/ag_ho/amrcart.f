c
c ----------------------------------------------------------------
c
       program amrcart
c
c main program
c
c program solves the hyperbolic 2-d equation:
c
c              u  +  f(u)    + g(u)   = 0
c               t         x        y
c
c  using rectangular mesh refinement.
c  no rotated rectangles are used in this version.
c  code assumes lower left domain at (0,0).
c
      implicit double precision (a-h,o-z)
      dimension          work(20)
      character * 8      pltfile, infile, outfile
      character * 12     rstfile
      logical            vtime,rest,steady,quad,nolimiter
      logical            graf
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead
      include "cirr.i"
      include "calloc.i"
      include "quadrature.i"
      common  /space/   lfree(150,2), lenf, idimf
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold, ilts,ibc,imodal
      common /center/ sx,sy
      common /stats/   evol, rvol, rvoll(maxlv), lentot, lenmax
      common /cloops/  xloops(10),yloops(10),nloops
      common /order2/   ssw, quad, nolimiter
      common /runparam/ ftime,outtime,nstopc, iprint
      common/eqncount/ieqncount
      common/steadydata/steadydiff, isteadysim
      integer oldmode

c
c 
        lwidth = 4  ! number ghost cells
c
        infile  = 'incart'
        outfile = 'outcart'
        pltfile = 'plot.3'
        rstfile = 'restart.data'

        iprint = -10
c
c     open(3,file=pltfile,status='unknown',form='formatted')
      open(5,file=infile,status='old',form='formatted')
      open(13,file='fort.13',status='unknown',form='formatted')
c
c     read(5,*) iout,nstopc,mxnest,iorder,iousr,nvar
      read(5,*) iout
      read(5,*) nstopc
      read(5,*) mxnest
      read(5,*) iorder
      read(5,*) ismp   
      read(5,*) iouser
      read(5,*) nvar
c     read(5,*) kcheck,ibuffx,ibuffy,ismp
      read(5,*) kcheck
      read(5,*) ibuffx
      read(5,*) ibuffy
c     read(5,*) ssw,quad,nolimiter
      read(5,*) ssw
      read(5,*) quad
      read(5,*) nolimiter
      read(5,*) (intrat(i),i=1,mxnest)

c     read(5,*) tol,cut,cdist,hxposs(1),hyposs(1),possk(1),cfl
      read(5,*) tol
      read(5,*) cut 
      cdist = 1.5d0
      !!read(5,*) hxposs(1)
      !!read(5,*) hyposs(1)
      read(5,*) ncells_x
      read(5,*) ncells_y
      read(5,*) possk(1)
      read(5,*) cfl
c     read(5,106) vtime,graf,rest,steady
      read(5,106) vtime
      read(5,106) graf
      read(5,106) rest
      read(5,106) steady
 106  format(4l1)
c     read(5,*) gamma, xprob, yprob, iprob
      gamma = 1.4d0
      read(5,*) xprob
      read(5,*) yprob
      hxposs(1) = xprob/dfloat(ncells_x)
      hyposs(1) = yprob/dfloat(ncells_y)
      read(5,*) iprob
      read(5,*) nloops
      do i = 1,  nloops
        read(5,*)xloops(i),yloops(i)
      end do

      read(5,*) outtime
      read(5,*) ftime
      read(5,*) ieqncount

      read(5,*) ilts
      read(5,*) ihob
      read(5,*) ibc

      read(5,*) sx
      read(5,*) sy

      read(5,*) imodal

      if(nstopc .eq. -2) then
      isteadysim = 1
      else
      isteadysim = 0
      endif


       if (ismp .eq. 5 .and. .not. quad) then
         write(6,*)" cannot have ismp 5 without quadratic recon."
         write(6,*) "resetting quad to true"
         quad = .true.
       endif

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
                                                                                                                1.d0/6.d0, 1.d0/6.d0,2.d0/3.d0, 0.d0 /)
      wtri_ho(1:7)=(/0.225000000000000d0/2.d0, 0.125939180544827d0/2.d0,
     .0.125939180544827d0/2.d0, 0.125939180544827d0/2.d0,
     .0.132394152788506d0/2.d0, 0.132394152788506d0/2.d0,
     .0.132394152788506d0/2.d0/)
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
     .                   0.312865496004874d0 /)                                                                                                                1.d0/6.d0, 1.d0/6.d0,2.d0/3.d0, 0.d0 /)
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
     . 8.5801403354407302d-01/)                                                                                                                1.d0/6.d0, 1.d0/6.d0,2.d0/3.d0, 0.d0 /)
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




       bzonex = dfloat(ibuffx)
       bzoney = dfloat(ibuffy)
c
c      ## if ar(k)/ar(lstgrd) less than gradThreshold dont use gradients in that cell
       gradThreshold = 1.d-7  ! other possibility is hx*hy. 

       if (rest) go to 10
          open(6,file=outfile,status='unknown',form='formatted')
          lentot = 0
          lenmax = 0
          rvol   = 0.0d0
          gamma1 = gamma - 1
          do 8 i   = 1, mxnest
 8           rvoll(i) = 0.0d0
          evol   = 0.0d0
          call   stst1
          call   domain (nvar,lfix,vtime,quad)
          call   setgrd (nvar,cut,cdist,lfix,quad)
          time = 0.0d0
          nstart = 0
          go to 20
c
 10   continue
      open(6,file=outfile,status='old',access='append',form='formatted')
      open(9,file=rstfile,status='old',form='unformatted')
      rewind 9
      call restrt(nsteps,time,nvar,lfix)
      nstart = nsteps
c
c  print out program parameters for this run
c
 20   write(6,107)tol,iorder,ismp,kcheck,ibuffx,ibuffy,cut,
     1            cdist,mxnest,cfl,gamma,ssw,quad,
     2            xprob, yprob,lfix,(intrat(i),i=1,mxnest)
107   format(/
     *       ' mesh refinement parameters:',//,
     *       ' error tol            ',e12.5,/,
     *       ' order of integrator     ',i9,/,
     *       ' small cell problem      ',i9,/,
     *       ' (1=wave prop,2=flux rd,3=boxes)'   ,/,
     *       ' error checking interval ',i9,/,
     *       ' buffer zone size (x)    ',i9,/,
     *       ' buffer zone size (y)    ',i9,/,
     *       ' volume ratio cutoff  ',e12.5,/,
     *       ' cluster separation   ',e12.5,/,
     *       ' max. refinement level   ',i9,/,
     *       ' cfl # (if var. delt) ',e12.5,/,
     *       ' gamma                ',e12.5,/,
     *       ' ssw                  ',e12.5,/,
     *       ' quad (fit at cut cells) ',l9,/,
     *       ' xprob (right bndry)  ',e12.5,/,
     *       ' yprob ( top  bndry)  ',e12.5,/,
c    *       ' alpha (flat plate)   ',e12.5,/,
c    *       ' Re (viscous flow)    ',e12.5,/,
     *       ' # fixed ref. levels     ',i9,/,
     *       ' refinement ratios       ',11i5,//)
c
c
      call outtre (mstart,.false.,nvar)
      nplot = 0
C       call dumptec(1, lfine,nvar,nplot,time)


!      call conck(1,nvar,time)


c     call valout(0,1,lfine,time,nvar)
c   
c     # print initial data on outcart:
c     if (.not. graf) call outtre(mstart,.true.,nvar)
c

      dx   = rnode(9, 1)
      dy   = rnode(10,1)
      xlow = rnode(1,1) - lwidth*dx
      ylow = rnode(2,1) - lwidth*dy

      call tick(nvar,iout,nstart,nstopc,cut,cdist,vtime,
     1          work,time,iousr,steady,lfix,quad,nplot)

!      call outgeom(1, node(5,1), node(6,1), alloc(node(14,1)),
!     .             xlow, ylow, dx, dy)
c      call outtre(mstart,graf,nvar)
c     call dumptec(1, lfine,nvar,nplot,time)
c
c  checkpoint everything for possible future restart
c
c     call check(nstopc,time,lfix)
      if (iprob .eq. 19 .or. iprob .eq. 20) then
         call errdrive(1,1,time,nvar,iprob)
      endif
c
c report on statistics
c
      write(6,*)
      write(6,*)
      write(6,*) 'current space usage = ',lentot
      write(6,*) 'maximum space usage = ',lenmax
      write(6,*) 'number of cells advanced in time = ',rvol
      do 60 level = 1,mxnest
 60     write(6,*) 'number of cells advanced in time on level ', level,
     & ' = ', rvoll(level)
 
      write(6,*) 'number of cells advanced for error estimation   = ', 
     &   evol
      if (evol+rvol .gt. 0.) then
         ratmet = rvol / (evol+rvol) * 100.0d0
      else 
         ratmet = 0.0d0
      endif
      write(6,*) 'percentage of cells advanced in time  = ', ratmet
c
      write(6,108)
108   format(//,' ------  end of mr integration --------  ')
c
      stop
      end
