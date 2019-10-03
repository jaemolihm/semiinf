#=====================================================
# Preprocessing flags
#=====================================================
# Parallel execution:
#       __MPI         MPI parallelization
#       __MPI_MODULE      Compile with "USE mpi" instead of "include 'mpif/h'"
DFLAGS = -fpp -D__MPI


#=====================================================
# For Linux with intel version 11/12 on 64bit machines
#=====================================================
F90    = ifort
MPIF90 = mpifort
FCOPTS = $(DFLAGS) -O2 -g -traceback -assume buffered_io
LDOPTS = $(DFLAGS) -O2 -g -traceback -assume buffered_io

#========================================================
# Intel mkl libraries. Set LIBPATH if not in default path
#========================================================

LIBDIR = /appl/compiler/intel2013_sp1/composer_xe_2013_sp1.2.144/mkl/lib/
LIBS   =  -L$(LIBDIR) -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lpthread


