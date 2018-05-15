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
OBJS = main.o

$(EXEC): $(OBJS)
	${CC} -o $@ $^ $C $(INC) $(LIB)

%.o: %.c
	${CC} -c -o $@ $< $C $(INC)

clean:
	rm -rf $(TARGETS) *.o $(EXEC)

