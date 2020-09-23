c      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c      :::::   data structure  for irregular cell info.
c      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
 
      parameter (irrsize=9533)
	  parameter ( ijsize = 5200 )
      common  /cirr/ poly(10,2,irrsize),ar(-1:irrsize),
     * points(24,2,irrsize),wt(24,-1:irrsize),
     * xcirr(irrsize),ycirr(irrsize),ix(irrsize),
     * iy(irrsize), nxtirr(irrsize),ipad,
     * volMerge(irrsize),
     * xcentMerge(irrsize),ycentMerge(irrsize),
     * iidx(irrsize,10), jidx(irrsize,10),
     * ncount(irrsize), numhoods(ijsize,ijsize),
     * ipad2,
     * dmergeshifts(7, irrsize),
     * dcubicshifts(4, irrsize),
     * mioff(ijsize,ijsize),mjoff(ijsize,ijsize),
     * qmshift(3, irrsize),
     * bdry(4,2,irrsize), nbdry,
     * ipad3,ar_ho(-1:irrsize),
     * xcirr_ho(irrsize),ycirr_ho(irrsize),
     * ihob, ipad4,poly_ho(10,2,irrsize),
     * ipad5, ipad7,dcubicshifts_ho(4, irrsize),
     * volMerge_ho(irrsize),
     * xcentMerge_ho(irrsize),ycentMerge_ho(irrsize),
     * qmshifts_ho(7, irrsize)

