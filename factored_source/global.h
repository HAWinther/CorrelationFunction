#ifndef _GLOBAL_H
#define _GLOBAL_H

#ifdef USE_OMP
#include <omp.h>
#elif defined(USE_MPI)
#include <mpi.h>
#endif

extern int mpi_rank;
extern int mpi_size;

#define MIN(x,y) (x > y ? y : x)
#define MAX(x,y) (x > y ? x : y)
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define SPEED_OF_LIGHT_IN_KM_PER_SEC 299792.458

extern double r_of_z(double z);

#endif
