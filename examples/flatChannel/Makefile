
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = amrclaw                  # Clawpack package to use
EXE = xamr                          # Executable to create
SETRUN_FILE = setrun.py             # File containing function to make data
OUTDIR = .                          # Directory for output
#OUTDIR = _output                    # Directory for output
SETPLOT_FILE = setplot.py           # File containing function to set plots
PLOTDIR = _plots                    # Directory for plots

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
#FC = ifort
FC = gfortran-8
#FC = ifort
#FFLAGS ?= -O -C -CB -CU -fpe0 -ftrapuv 
#FFLAGS ?= -O2 -C -CB -CU -fpe0 -ftrapuv -qopenmp  -mcmodel=large
FFLAGS ?= -g -fbounds-check -gno-strict-dwarf  -ffpe-trap=invalid,overflow,zero 

# ---------------------------------
# package sources for this program:
# ---------------------------------

#AMRLIB = $(CLAW)/amrclaw/src/2d
#AMRLIB = /Users/berger/ho_amrclaw-amrcart/src/2d
AMRLIB = /Users/berger/ho_amrclaw_amrcart/src/2d
include $(AMRLIB)/Makefile.amr_2d

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source 
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = \

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------

MODULES = \

SOURCES = \
  setprob.f \
  qinit.f \
  fbody.f \
  bc2amr.f90 \
  channelInit.f \
  $(AMRLIB)/fixed_vector_lim_allineq.f \
  $(AMRLIB)/fixed2_Simplex_allineq.f \
  $(AMRLIB)/vrmc_llf.f \
  $(AMRLIB)/mymethod.f \
  $(AMRLIB)/SRD_cellMerge.f \
  outtec.f \
  dumptec.f \

#  $(AMRLIB)/lp_qslopesWithGhostAvgQuadratic.f \
#  $(AMRLIB)/lp_qslopesPtQuad2.f \
#  $(AMRLIB)/SRD_normal2.f \
#  $(AMRLIB)/both_vector_lim_allineq.f \
#  $(AMRLIB)/lp_qslopesWithGhostQuadratic.f \

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

