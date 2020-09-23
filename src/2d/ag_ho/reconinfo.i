c      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c      :::::   data structure  for reconstruction info.
c      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
 
	  parameter ( iijjsize = 5200 )
      common  /reconinfo/ inuf(iijjsize,iijjsize),
     * iir(iijjsize,iijjsize), jjr(iijjsize,iijjsize)
