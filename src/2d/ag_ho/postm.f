c
c ----------------------------------------------------------
c
       subroutine postm(valbig,val,nvar,maxip1,maxjp1,mitot,mjtot,
     1                  lwidth)
c
      implicit double precision (a-h,o-z)
       dimension valbig(mitot,mjtot,nvar),val(maxip1,maxjp1,nvar)
c
c the method program has returned new values in valbig, which
c included space for lwidth dummy points all around. copy the
c interior solution back into val array which the
c rest of the amr program uses.
c
       do 10 i = lwidth+1, mitot-lwidth
       do 10 j = lwidth+1, mjtot-lwidth
       do 10 ivar = 1, nvar
	  val(i-lwidth+1,j-lwidth+1,ivar) = valbig(i,j,ivar)
 10    continue
       return
       end
