MPIF90=mpifort
#FFLAG=-g -check all -fpe0 -warn -traceback -debug extended
FFLAG=-g 
FFLAG=-O2 -g -traceback

LIBDIR=/appl/compiler/intel2013_sp1/composer_xe_2013_sp1.2.144/mkl
LIBS=-L$(LIBDIR) -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lpthread

OBJECTS=constants.o  hamiltonian.o  iter_bulk.o  iter_slab.o  main.o  parameters.o  postprocess_green.o

default: semiinf

semiinf: constants.o  hamiltonian.o  iter_bulk.o  iter_slab.o  main.o  parameters.o  postprocess_green.o
	$(MPIF90) ${FFLAG} ${LIBS} -o semiinf.x $^

constants.o:
	$(MPIF90) ${FFLAG} -c constants.f90

parameters.o: constants.o
	$(MPIF90) ${FFLAG} -c parameters.f90

hamiltonian.o: constants.o parameters.o
	$(MPIF90) ${FFLAG} -c hamiltonian.f90

postprocess_green.o: constants.o parameters.o hamiltonian.o
	$(MPIF90) ${FFLAG} -c postprocess_green.f90

iter_bulk.o: constants.o parameters.o hamiltonian.o postprocess_green.o
	$(MPIF90) ${FFLAG} -c iter_bulk.f90

iter_slab.o: constants.o parameters.o hamiltonian.o postprocess_green.o
	$(MPIF90) ${FFLAG} -c iter_slab.f90

main.o: constants.o parameters.o hamiltonian.o postprocess_green.o iter_bulk.o iter_slab.o
	$(MPIF90) ${FFLAG} -c main.f90

clean:
	rm semiinf.x *.o
