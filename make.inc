
#=====================================================
# For Linux with intel version 11/12 on 64bit machines
#=====================================================
F90    = ifort
MPIF90 = mpifort
FCOPTS = -O2 -g -traceback -assume buffered_io -fpp
LDOPTS = -O2 -g -traceback -assume buffered_io -fpp

#========================================================
# Intel mkl libraries. Set LIBPATH if not in default path
#========================================================

LIBDIR = /appl/compiler/intel2013_sp1/composer_xe_2013_sp1.2.144/mkl/lib/
LIBS   =  -L$(LIBDIR) -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lpthread


