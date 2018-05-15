#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "grid.h"
#include "galaxycat.h"
#include "correlators.h"
#include "binning.h"
#include "gsl_wrapper.h"
#include "cosmo.h"

int mpi_rank = 0, mpi_size = 1;

int main(int argc, char **argv){

#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  //====================================================
  // Read parameters from stdin
  //====================================================
  if(argc < 7){
    if(mpi_rank == 0)
      printf("Run as ./corr [filename_mock] [filename_random] [filename_output] [nbins] [rmax] [OmegaM]\n");
    exit(1);
  }

  char *filename_galaxies = argv[1];
  char *filename_random   = argv[2];
  char *filename_output   = argv[3];
  int nbins               = atoi(argv[4]);
  double rmax             = atof(argv[5]);
  double OmegaM           = atof(argv[6]);
  double box;

  if(mpi_rank == 0){
    printf("\n====================================\n");
    printf("Parameters:\n");
    printf("====================================\n");
    printf("Galaxy Catalog: [%s]\n", filename_galaxies);
    printf("Random Catalog: [%s]\n", filename_random);
    printf("Rmax:           [%5.2lf] Mpc/h\n", rmax);
    printf("nbins:          [%i]\n", nbins);
    printf("OmegaM:         [%0.3lf]\n", OmegaM);
#ifdef WEIGHTS
    printf("We are including the weights\n");
#else  
    printf("We are NOT including the weights\n");
#endif
#ifdef BRUTEFORCE
    printf("Will compute pair counts BRUTE-FORCE\n");
#endif
#ifdef USE_OMP
#pragma omp parallel
    {
      if(omp_get_thread_num() == 0) 
        printf("We are using OpenMP nthreads = %i\n", omp_get_num_threads());
    }
#elif defined(USE_MPI)
    if(mpi_rank == 0)
      printf("We are using MPI nsize = %i\n", mpi_size);
#else
    printf("We are NOT using OpenMP or MPI\n");
#endif
    printf("====================================\n");
  }

  //====================================================
  // Create a spline of r(z)
  //====================================================
  global_spline_rofz = create_rofz_spline(OmegaM); 

  //====================================================
  // Count number of galaxies and randoms in the files
  //====================================================
  int ngalaxies = count_lines_in_file(filename_galaxies);
  int nrandom   = count_lines_in_file(filename_random);

  //====================================================
  // Read galaxies from file and store in catalog
  //====================================================
  GalaxyCatalog *galaxy_cat = read_galaxies_from_file(filename_galaxies, ngalaxies);
  GalaxyCatalog *random_cat = read_galaxies_from_file(filename_random, nrandom);

  //====================================================
  // Now we have positions on the sky. Compute the box
  // we need to encompas them and shift positions to [0,box]^3
  //====================================================
  compute_boxsize_shift_positions(galaxy_cat, random_cat, &box);

  //====================================================
  // Make some grids
  //====================================================
  Grid *galaxy_grid = create_grid(ngalaxies, rmax, box);
  Grid *random_grid = create_grid(nrandom, rmax, box);

  //====================================================
  // Allocate and assign galaxies to the cells
  //====================================================
  add_galaxies_to_cells(galaxy_grid, galaxy_cat);
  add_galaxies_to_cells(random_grid, random_cat);

  //====================================================
  // Allocate histograms for DD, DR and RR  
  //====================================================
  PairCountBinning *DD = create_binning(nbins, rmax);
  PairCountBinning *DR = create_binning(nbins, rmax);
  PairCountBinning *RR = create_binning(nbins, rmax);

  // Total number of (weighted) pairs used to normalize corr function
  DD->norm = 0.5*(pow2(galaxy_cat->sum_w) - galaxy_cat->sum_w2);
  RR->norm = 0.5*(pow2(random_cat->sum_w) - random_cat->sum_w2);
  DR->norm = galaxy_cat->sum_w * random_cat->sum_w;

  //====================================================
  // Compute correlation function
  //====================================================
#ifdef BRUTEFORCE
  brute_force_pair_counting(galaxy_cat, DD);
  brute_force_cross_pair_counting(random_cat, galaxy_cat, DR);
  brute_force_pair_counting(random_cat, RR);
#else
  grid_pair_counting(galaxy_grid, DD);
  grid_cross_pair_counting(galaxy_grid, random_grid, DR);
  grid_pair_counting(random_grid, RR);
#endif
  compute_correlation_function(DD, DR, RR, filename_output);

  //====================================================
  // Free up memory allocated above
  //====================================================
  free_grid(galaxy_grid);
  free_grid(random_grid);
  free_cat(galaxy_cat);
  free_cat(random_cat);
  free_binning(DD);
  free_binning(DR);
  free_binning(RR);
  Free_GSL_Spline(global_spline_rofz);

#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}

