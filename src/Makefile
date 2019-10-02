# Makefile for semiinf.x

include ../make.inc

OBJECTS=comms.o  hamiltonian.o  iter_bulk.o  iter_slab.o  main.o  parameters.o  postprocess_green.o

default: all

all: semiinf

semiinf: comms.o  hamiltonian.o  iter_bulk.o  iter_slab.o  main.o  parameters.o  postprocess_green.o
	$(MPIF90) ${FCOPTS} ${LIBS} -o semiinf.x $^

comms.o:
	$(MPIF90) ${FCOPTS} -c comms.f90

parameters.o: comms.o
	$(MPIF90) ${FCOPTS} -c parameters.f90

hamiltonian.o: comms.o parameters.o
	$(MPIF90) ${FCOPTS} -c hamiltonian.f90

postprocess_green.o: comms.o parameters.o hamiltonian.o
	$(MPIF90) ${FCOPTS} -c postprocess_green.f90

iter_bulk.o: comms.o parameters.o hamiltonian.o postprocess_green.o
	$(MPIF90) ${FCOPTS} -c iter_bulk.f90

iter_slab.o: comms.o parameters.o hamiltonian.o postprocess_green.o
	$(MPIF90) ${FCOPTS} -c iter_slab.f90

main.o: comms.o parameters.o hamiltonian.o postprocess_green.o iter_bulk.o iter_slab.o
	$(MPIF90) ${FCOPTS} -c main.f90
	- ( cd ../; ln -fs src/semiinf.x . ; )

clean:
	- rm -f semiinf.x ../semiinf.x *.o *.mod 
