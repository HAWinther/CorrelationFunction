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
  double rmin    = 0.0;
  double rmax    = 192.0;
  double OmegaM  = 0.3;
  int nbins      = 32;
  char filename_galaxy[100] = "galaxy_mock.dat";
  char filename_random[100] = "random.dat";
  
  //=================================
  // Call the library
  //=================================
  CUTER_set_verbose(TRUE);
  // CUTER_set_bintype_log(); rmin = 1.0; // Switch to logarithminc binning
  BinnedCorrelationFunction *bcf = CUTER_correlation_function(filename_galaxy, filename_random, 
       FILEFORMAT_RA_DEC_Z_POSITIONS, nbins, rmin, rmax, OmegaM);
  if(bcf == NULL) exit(1);

  //=================================
  // Output correlation function
  //=================================
  if(rank == 0){
    printf("\nThe binned correlation function: \n");
    for(i = 0; i < nbins; i++){
      printf("r: %6.2lf  DD(r): %10i RR(r): %10i DR(r): %10i  CorrFunc(r): %6.2lf  PoissonErr %6.2e\n", bcf->r[i], (int) bcf->DD[i], (int) bcf->RR[i], (int) bcf->DR[i], bcf->corr_func[i], bcf->err_corr[i]);
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
