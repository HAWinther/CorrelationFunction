#include <stdio.h>
#include <stdlib.h>
#if defined(USE_OMP)
#include <omp.h>
#elif defined(USE_MPI)
#include <mpi.h>
#endif
#include "../src/cuter_library.h"

#define TRUE  1
#define FALSE 0

int main(int argc, char **argv){
  int rank = 0;
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  int i;
 
  //=================================
  // Set the parameters
  //=================================
  int nbins   = 10;
  double rmin = 1.0;
  double rmax = 100.0;
  double box  = 400.0;
  char filename[100] = "periodic_mock.dat";
  
  //=================================
  // Set options and call the library
  //=================================
  CUTER_set_verbose(TRUE);
  CUTER_set_bintype_lin(); // CUTER_set_bintype_log(); 
  BinnedCorrelationFunction *bcf = CUTER_correlation_function_periodic(filename, nbins, rmin, rmax, box);
  if(bcf == NULL) exit(1);

  //=================================
  // Output correlation function
  //=================================
  if(rank == 0){
    printf("\nThe binned correlation function: \n");
    for(i = 0; i < nbins; i++){
      printf("r: %5.2lf  DD(r): %10i  CorrFunc(r): %5.2lf  PoissonErr %5.2e\n", bcf->r[i], (int) bcf->DD[i], bcf->corr_func[i], bcf->err_corr[i]);
    }
  }

  //=================================
  // Clean up
  //=================================
  free_binned_correlation_function(bcf);

#ifdef USE_MPI
  MPI_Finalize();
#endif
}
