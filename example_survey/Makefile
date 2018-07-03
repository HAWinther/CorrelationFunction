#=================================================
# Hans A. Winther (2018) (hans.a.winther@gmail.com)
#=================================================

SHELL := /bin/bash

#====================================
# Executable
#====================================
EXEC = CUTER

#====================================
# Set options, uncomment to use it
#====================================
#USE_OMP     = 1
#USE_WEIGHTS = 1
#USE_MPI     = 1
#BRUTEFORCE  = 1

#====================================
# Set compiler and add defines
#====================================
CC = g++-mp-6
OPTIONS = 
ifdef USE_MPI
CC = mpicc-openmpi-gcc6
OPTIONS += -DUSE_MPI
EXEC = CUTER_MPI
endif
ifdef USE_OMP
OPTIONS += -DUSE_OMP -fopenmp
endif
ifdef USE_WEIGHTS
OPTIONS += -DWEIGHTS
endif
ifdef BRUTEFORCE
OPTIONS += -DBRUTEFORCE
endif
C = -O3 $(OPTIONS)

#====================================
# GSL library
#====================================
LIB = -L/opt/local/lib/ -lgsl -lm
INC = -I/opt/local/include/

#====================================
# OBJECT FILES
#====================================
OBJS = ../src/gsl_spline_wrapper.o ../src/cuter_library.o main.o

$(EXEC): $(OBJS)
	${CC} -o $@ $^ $C $(INC) $(LIB)

%.o: %.c
	${CC} -c -o $@ $< $C $(INC)

clean:
	rm -rf $(TARGETS) *.o ../src/*.o $(EXEC)

