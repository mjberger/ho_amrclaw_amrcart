
c
c -------------------------------------------------------------------
c
       subroutine preppoly(p, rlink, pout)
       implicit double precision (a-h, o-z)
       dimension pout(20,2), p(10,2),rlink(9,2)

c
c      merge the 2 polygon descriptions into one big one -> pout.
c
       index = 0
       do 10 i = 1, 10
	  if (p(i+1,1) .eq. -11) then
c            move on
             go to 20
	  else
	    index = index + 1
	    pout(index,1) = p(i,1)
	    pout(index,2) = p(i,2)
	  endif
 10     continue
 20    continue

       do 30 i = 2, 10
	  if (rlink(i,1) .eq. -11) then
	     if (i .eq. 2) then
c              no links in rlink (i.e. reg. cell). close the loop
	       index = index + 1
	       pout(index,1) = p(1,1)
	       pout(index,2) = p(1,2)
	     endif
	     index = index + 1
	     pout(index,1) = -11
	     go to 40
	  else
	     index = index + 1
	     pout(index,1) = rlink(i,1)
	     pout(index,2) = rlink(i,2)
	  endif
 30    continue

 40    continue

       return
       end
