      real*8 rline,wline,rtri,stri,wtri,rquad,squad,wquad
      real*8 rline_ho,rtri_ho,stri_ho,wtri_ho,wline_ho
      integer nlinequad,ntriquad,nquadquad,ntriquad_ho,nlinequad_ho

      common /linequadrature/ rline(4), wline(4), nlinequad
      common /triquadrature/ rtri(4),stri(4), wtri(4), ntriquad
      common /quadquadrature/ rquad(4),squad(4),wquad(4),nquadquad
      common /hotriquadrature/rtri_ho(33),stri_ho(33),wtri_ho(33),ntriquad_ho
      common /holinequadrature/ rline_ho(6), wline_ho(6), nlinequad_ho
