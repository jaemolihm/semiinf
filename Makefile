# Makefile for semiinf
sinclude make.inc

default: all

all: semiinf

semiinf:     
	( cd src ; $(MAKE) all || exit 1 )

clean : 
	( cd src ; $(MAKE) clean )
