
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
OUTDIR = _output                    # Directory for output
SETPLOT_FILE = setplot.py           # File containing function to set plots
PLOTDIR = _plots                    # Directory for plots

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
3FFLAGS ?= -g -C -CB -CU -fpe0 -ftrapuv 
#FFLAGS ?=  -C -CB -CU -fpe0 -ftrapuv -mcmodel=large  -g -qopenmp
#FFLAGS ?=  -C -CB -CU -fpe0 -ftrapuv -O1 -qopenmp
FFLAGS ?= -fopenmp -g -fbounds-check -ffpe-trap=invalid,overflow,zero 

# ---------------------------------
# package sources for this program:
# ---------------------------------

#AMRLIB = $(CLAW)/amrclaw/src/2d
AMRLIB = /Users/berger/amrclaw-amrcart/src/2d
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
  dumpBndry.f \
  outtec.f \
  bc2amr.f90 \
  $(AMRLIB)/vrm_llf.f \
  $(AMRLIB)/vrmc_llf.f \
  $(AMRLIB)/revdivmusclmethod2.f \
  $(AMRLIB)/fixed2_Simplex_allineq.f \
  $(AMRLIB)/fixed_vector_lim_allineq.f \

#  $(AMRLIB)/mymethod.f \
#  $(AMRLIB)/lp_qslopesPtQuad2.f \
#  $(AMRLIB)/qslopesWithGhostAvgQuadratic.f \
#  $(AMRLIB)/slopes.f \
#  $(AMRLIB)/SRD_quadratic.f \
#  $(AMRLIB)/SRD_normal2.f \
#  $(AMRLIB)/SRD_normalWithGhostQuadratic.f \
#  $(AMRLIB)/revdivmusclmethod2.f \
#  $(AMRLIB)/lp_qslopes.f \
#  $(AMRLIB)/both_vector_lim_allineq.f \

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

