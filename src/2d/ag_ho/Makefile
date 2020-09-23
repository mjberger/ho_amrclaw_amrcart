
#
#		Makefile for the amr code:
#
#	To compile and link, type:    	make xcart
#
#	To compile a file.f type:	make file.o
#
#	To combine the source type:	make cat
#
#
# sun f77 compiler  -C is for index bounds checking
# DEBUGGER -g -fbounds-check -ffpe-trap=invalid,overflow,zero
#	$(HBOX_SRC)/outgeom.f 
FC = gfortran
FFLAGS = -g -fbounds-check -ffpe-trap=invalid,overflow,zero -Wuninitialized -gno-strict-dwarf
#FFLAGS = -O3 -mcmodel=large
#FFLAGS = -freal-8-real-16 -g -fbounds-check -ffpe-trap=invalid,overflow,zero -Wuninitialized -mcmodel=large
#FFLAGS = -O3 -mcmodel=large -freal-8-real-16

LFLAGS =  $(FFLAGS)
F77    = $(FC) 
LINK	= $(F77)

HBOX_SRC = .
#EX = channel
EX = supersonic
#    $(EX)/vrmc_colella.f\
#    $(EX)/vrm.f\

FOBJS = $(FSOURCES:.f=.o)

FSOURCES   = \
    $(EX)/vrm.f\
    $(EX)/physbd.f  \
    $(EX)/pphysbdlin.f\
    $(EX)/irregflux_gauss.f \
    $(HBOX_SRC)/qr_solve.f90 \
    $(HBOX_SRC)/limiter.f \
    $(HBOX_SRC)/limiter_srd.f \
	$(HBOX_SRC)/SRD_cellMerge_ho.f \
    $(HBOX_SRC)/merge_shifts_ho.f \
    $(HBOX_SRC)/outgeom.f \
    $(HBOX_SRC)/cholesky.f \
    $(HBOX_SRC)/solve.f \
    $(HBOX_SRC)/outgmsh.f \
    $(HBOX_SRC)/method.f \
	$(HBOX_SRC)/eval_dU_gauss.f \
    $(HBOX_SRC)/getXface_gauss.f \
    $(HBOX_SRC)/getYface_gauss.f \
    $(HBOX_SRC)/qslopes.f \
	$(HBOX_SRC)/slopes.f \
	$(HBOX_SRC)/SRD_cellMerge.f \
    $(HBOX_SRC)/advanc.f  \
    $(HBOX_SRC)/prem.f  \
    $(HBOX_SRC)/postm.f  \
    $(HBOX_SRC)/igetsp.f \
    $(HBOX_SRC)/reclam.f  \
    $(HBOX_SRC)/birect.f  \
    $(HBOX_SRC)/bisect.f  \
    $(HBOX_SRC)/bound.f  \
    $(HBOX_SRC)/bdcrse.f  \
    $(HBOX_SRC)/check.f \
    $(HBOX_SRC)/cntain.f  \
    $(HBOX_SRC)/colate.f  \
    $(HBOX_SRC)/errest.f  \
    $(HBOX_SRC)/errf1.f  \
    $(HBOX_SRC)/errf2.f  \
    $(HBOX_SRC)/setir2.f  \
    $(HBOX_SRC)/fill.f  \
    $(HBOX_SRC)/filval.f  \
    $(HBOX_SRC)/filpatch.f  \
    $(HBOX_SRC)/filpc.f  \
    $(HBOX_SRC)/flglvl.f \
    $(HBOX_SRC)/flgbuf.f  \
    $(HBOX_SRC)/fluxad.f  \
    $(HBOX_SRC)/fluxsv.f  \
    $(HBOX_SRC)/getsrc.f  \
    $(HBOX_SRC)/ginit.f  \
    $(HBOX_SRC)/grdfit.f  \
    $(HBOX_SRC)/intfil.f  \
    $(HBOX_SRC)/join.f  \
    $(HBOX_SRC)/merge.f  \
    $(HBOX_SRC)/moment.f  \
    $(HBOX_SRC)/nestck.f  \
    $(HBOX_SRC)/inbase.f  \
    $(HBOX_SRC)/ovrlap.f  \
    $(HBOX_SRC)/phys.f \
    $(HBOX_SRC)/tirr.f  \
    $(HBOX_SRC)/prepc.f  \
    $(HBOX_SRC)/prepf.f  \
    $(HBOX_SRC)/projec.f \
    $(HBOX_SRC)/signs.f  \
    $(HBOX_SRC)/findcut.f  \
    $(HBOX_SRC)/smartbis.f  \
    $(HBOX_SRC)/putnod.f  \
    $(HBOX_SRC)/putsp.f  \
    $(HBOX_SRC)/regrid.f  \
    $(HBOX_SRC)/restrt.f  \
    $(HBOX_SRC)/setgrd.f  \
    $(HBOX_SRC)/sethk.f  \
    $(HBOX_SRC)/setirr.f  \
    $(HBOX_SRC)/setuse.f  \
    $(HBOX_SRC)/stst1.f  \
    $(HBOX_SRC)/tick.f  \
    $(HBOX_SRC)/conck.f \
    $(HBOX_SRC)/unpack.f  \
    $(HBOX_SRC)/trimbd.f  \
    $(HBOX_SRC)/trimgr.f  \
    $(HBOX_SRC)/update.f  \
    $(HBOX_SRC)/nodget.f \
    $(HBOX_SRC)/lstget.f  \
    $(HBOX_SRC)/lstput.f  \
    $(HBOX_SRC)/upbnd.f  \
    $(HBOX_SRC)/ziplst.f  \
    $(HBOX_SRC)/basic.f   \
    $(HBOX_SRC)/drawrg.f \
    $(HBOX_SRC)/outirr.f  \
    $(HBOX_SRC)/outbig.f  \
    $(HBOX_SRC)/outval.f \
    $(HBOX_SRC)/valout.f  \
    $(HBOX_SRC)/errdrive.f  \
    $(HBOX_SRC)/errout20.f \
    $(HBOX_SRC)/outvar.f  \
    $(HBOX_SRC)/errout19.f  \
    $(HBOX_SRC)/outmsh.f  \
    $(HBOX_SRC)/outlev.f  \
    $(HBOX_SRC)/outtre.f  \
    $(HBOX_SRC)/amrcart.f  \
    $(HBOX_SRC)/domain.f \
    $(HBOX_SRC)/irrlen.f   \
    $(HBOX_SRC)/makep.f  \
    $(HBOX_SRC)/polyin.f   \
    $(HBOX_SRC)/fbody.f   \
    $(HBOX_SRC)/setir21.f  \
    $(HBOX_SRC)/errlen.f \
    $(HBOX_SRC)/vctoprm.f  \
    $(HBOX_SRC)/prm.f \
    $(HBOX_SRC)/vpredx.f \
    $(HBOX_SRC)/vpredy.f \
    $(HBOX_SRC)/vprmtoc.f \
    $(HBOX_SRC)/rmnv.f \
    $(HBOX_SRC)/makebx.f \
    $(HBOX_SRC)/centrd.f \
    $(HBOX_SRC)/linfcn.f \
    $(HBOX_SRC)/p18.f  \
    $(HBOX_SRC)/p18f.f \
    $(HBOX_SRC)/p18mu.f \
    $(HBOX_SRC)/p18th.f \
    $(HBOX_SRC)/p18yw.f \
    $(HBOX_SRC)/p18tru.f \
    $(HBOX_SRC)/p18xw.f \
    $(HBOX_SRC)/p18fxw.f \
    $(HBOX_SRC)/zero2.f \
    $(HBOX_SRC)/p17tru.f \
    $(HBOX_SRC)/p19tru.f \
    $(HBOX_SRC)/p20tru.f \
    $(HBOX_SRC)/p19fn.f \
    $(HBOX_SRC)/p20fn.f \
    $(HBOX_SRC)/p16fn.f \
    $(HBOX_SRC)/p19old.f \
    $(HBOX_SRC)/p15fn.f \
    $(HBOX_SRC)/p15tru.f \
    $(HBOX_SRC)/weights.f \
    $(HBOX_SRC)/transfer.f \
    $(HBOX_SRC)/preppoly.f \
    $(HBOX_SRC)/makepout.f \
    $(HBOX_SRC)/qsltemp.f \
    $(HBOX_SRC)/qsintfil.f \
    $(HBOX_SRC)/tirrfil.f  \
    $(HBOX_SRC)/addtemp.f \
    $(HBOX_SRC)/addneighbors.f \
    $(HBOX_SRC)/addMoreNeighbors.f \
    $(HBOX_SRC)/dumptec.f \
    $(HBOX_SRC)/outtec.f \
    $(HBOX_SRC)/outresid.f \
    $(HBOX_SRC)/estdt.f \
    $(HBOX_SRC)/mkbdry.f \
    $(HBOX_SRC)/fs.f \
    $(HBOX_SRC)/zeroin.f \
    $(HBOX_SRC)/isig.f \
    $(HBOX_SRC)/finit.f  \
    $(HBOX_SRC)/getirrlen.f \
    $(HBOX_SRC)/cellsClose.f \
    $(HBOX_SRC)/merge_shifts.f \
    $(HBOX_SRC)/qmslopes.f \

all: $(FOBJS) Makefile
	$(LINK) $(LFLAGS) $(FOBJS) -o $(EX)/xcart

cat : $(FSOURCES)
	cat $(FSOURCES) > cart.f

clean:
	\rm -f  *.o \
    \rm -f  $(EX)/*.o
new:
	make clean
	make
