c
c ------------------------------------------------------------
c
       subroutine ziplst(listbc,maxsp)
 
      implicit double precision (a-h,o-z)
       dimension listbc(5,maxsp)
c
c  initialize list space to 0, since use 0 terminator if number of
c  cells to be fixed is less than maxsp
c
       do 10 isp = 1, maxsp
          listbc(1,isp) = 0
 10    continue
c
       return
       end
