#include <stdio.h>
#include <stdlib.h>
#if defined(USE_OMP)
#include <omp.h>
#elif defined(USE_MPI)
#include <mpi.h>
#endif
#include "cuter_library.h"

int main(int argc, char **argv){
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif
  int i;
 
  //=================================
  // Set the parameters
  //=================================
  double rmax = 100.0;
  double box  = 400.0;
  int nbins   = 10;
  char filename[100] = "HOD.dat";
  
  //=================================
  // Call the library
  //=================================
  cute_library_verbose = 1;
  BinnedCorrelationFunction *bcf = CUTER_correlation_function_periodic(filename, nbins, rmax, box);
 
  //=================================
  // Output correlation function
  //=================================
  printf("\nThe binned correlation function: \n");
  for(i = 0; i < nbins; i++){
    printf("r: %5.2lf  DD(r): %10i  CorrFunc(r): %5.2lf  PoissonErr %5.2e\n", bcf->r[i], (int) bcf->DD[i], bcf->corr_func[i], bcf->err_corr[i]);
  }

  //=================================
  // Clean up
  //=================================
  free_binned_correlation_function(bcf);
#ifdef USE_MPI
  MPI_Finalize();
#endif
}
