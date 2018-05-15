#####################################################
# Hans A. Winther (2018) (hans.a.winther@gmail.com)
#####################################################

SHELL := /bin/bash

# Set compiler
CC = gcc-mp-6

# GSL library
L = -L/opt/local/lib/ -lgsl -lgslcblas -lm
I = -I/opt/local/include/gsl/

C = -O3 -fopenmp $(OPTIONS)

EXEC = CUTER

TARGETS := $(EXEC)
all: $(TARGETS)

# OBJECT FILES
OBJS = main.o

# DEPENDENCIES
main.o		:

# HEADERS
HEADERS = 

CUTER: $(OBJS)
	${CC} -o $@ $^ $C $I $L $(SFML)

%.o: %.cpp $(HEADERS)
	${CC} -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o $(EXEC)

