#####################################################
# Hans A. Winther (2018) (hans.a.winther@gmail.com)
#####################################################

SHELL := /bin/bash

# Set compiler (MPI: mpicc-openmpi-gcc6)
CC = gcc-mp-6 

# GSL library
LIB = -L/opt/local/lib/ -lgsl -lm
INC = -I/opt/local/include/

# Options (-DUSE_MPI for MPI version)
OPTIONS = -DUSE_OMP -DWEIGHTS
C = -O3 -fopenmp $(OPTIONS)

# Executable
EXEC = CUTER

TARGETS := $(EXEC)
all: $(TARGETS)

# OBJECT FILES
OBJS = main.o gsl_wrapper.o binning.o correlators.o cosmo.o galaxycat.o grid.o

# DEPENDENCIES
main.o		    : gsl_wrapper.c binning.c cosmo.c galaxycat.c grid.c
correlators.o : galaxycat.c grid.c binning.c
cosmo.o       : gsl_wrapper.c
grid.o        : galaxycat.c

# HEADERS
HEADERS = cosmo.h gsl_wrapper.h binning.h correlators.h grid.h galaxycat.h

$(EXEC): $(OBJS)
	${CC} -o $@ $^ $C $(INC) $(LIB)

%.o: %.c $(HEADERS)
	${CC} -c -o $@ $< $C $(INC)

clean:
	rm -rf $(TARGETS) *.o $(EXEC)

