#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#ifdef OMP
#include <omp.h>
#endif
#define MIN(x,y) (x > y ? y : x)
#define MAX(x,y) (x > y ? x : y)
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define SPEED_OF_LIGHT_IN_KM_PER_SEC 299792.458

//====================================================
// 
// Very simple code to compute the radial monopole correlation
// function for data from a galaxy survey / mocks. Using
// grids to speed it up.
// Same (or a bit faster) speed as CUTE with OpenMP
//
// Input: galaxyfile, randomfile, nbins, rmax and OmegaM
// Output: [xi(r), err(r), DD(r), DR(r), RR(r)]

// Input is assumed to have format: [RA, DEC, z, (weight)]
// with RA/DEC in degrees. Rmax in Mpc/h
//
// Defines:
// WEIGHTS : Use weights in galaxy mock
// OMP : Parallelize using OpenMP
// BRUTEFORCE : Brute-force pair-counts
// OUTPUT_PAIRS_TO_SCREEN : print paircount to screen
// as we go along
//
// Written by Hans Winther (2018)
//
//====================================================

//====================================================
// A single galaxy
//====================================================
typedef struct Galaxy{
  double x[3];
#ifdef WEIGHTS
  double w;
#endif
} Galaxy;

//====================================================
// A galaxy catalog containing a list of galaxies and info
//====================================================
typedef struct GalaxyCatalog{
  int ngalaxies;    // Number of galaxies
  Galaxy *galaxies; // List of galaxies
  double sum_w;     // Sum of weight
  double sum_w2;    // Sum of weights^2
  int allocated;    // Is galaxies allocated or not?
} GalaxyCatalog;

//====================================================
// This is a single grid-cell
//====================================================
typedef struct Cell{
  int np;           // Number of galaxies in the cell
  Galaxy *galaxy;   // Array of galaxies
} Cell;

//====================================================
// A grid containing a list of cells and general info
//====================================================
typedef struct Grid{
  Cell *cells;                // List of cells
  double cell_size;           // Size of cells in Mpc/h
  int ngrid;                  // Number of gridcells per dimension
  int ngalaxies;              // Number of galaxies in cell
  int max_ix, max_iy, max_iz; // No galaxies in cells (i,j,k) where i>max_ix or j>max_iy etc. 
  int allocated;              // Is cells allocated or not?
} Grid;

//====================================================
// Container for a binning
//====================================================
typedef struct PairCountBinning{
  int nbins;             // Number of linear bins between r=0 and r=RMAX
  double rmax;           // The RMAX we bin up to
  double norm;           // Normalization of paircount: (weighted) total number of pairs
  double *paircount;     // The paircounts
  int allocated;         // Is paircount allocated or not?
} PairCountBinning;

//====================================================
// GSL spline
//====================================================
typedef struct GSL_Spline {
  gsl_spline *spline;
  gsl_interp_accel *xacc;
  double xmin, xmax;
  int allocated;
} GSL_Spline;

//====================================================
// Global r(z) spline
//====================================================
GSL_Spline *global_spline_rofz;

//====================================================
// Functions defined below
//====================================================
Grid *create_grid(int ngalaxies, double rmax, double box);
GalaxyCatalog *read_galaxies_from_file(char *filename, int npart);
PairCountBinning *create_binning(int nbins, double rmax);
GSL_Spline *create_rofz_spline(double OmegaM);
int  ode_rofz(double z, const double r[], double drdz[], void *params);
int  count_lines_in_file(char *filename);
void free_binning(PairCountBinning *pc);
void add_galaxies_to_cells(Grid *grid, GalaxyCatalog *cat);
void brute_force_pair_counting(GalaxyCatalog *cat, PairCountBinning *pc);
void brute_force_cross_pair_counting(GalaxyCatalog *cat, GalaxyCatalog *cat2, PairCountBinning *pc);
void grid_pair_counting(Grid *grid, PairCountBinning *pc);
void grid_cross_pair_counting(Grid *grid, Grid *grid2, PairCountBinning *pc);
void compute_correlation_function(PairCountBinning *DD, PairCountBinning *DR, PairCountBinning *RR, char *filename);
void outputGalaxies(GalaxyCatalog *cat, char *filename);
void compute_boxsize_shift_positions(GalaxyCatalog *cat, GalaxyCatalog *cat2, double *box);
void Create_GSL_Spline(GSL_Spline *splinecontainer, double *x, double *y, int nx);
void Free_GSL_Spline(GSL_Spline *splinecontainer);
void free_cat(GalaxyCatalog *cat);
void free_grid(Grid *grid);
double Lookup_GSL_Spline(GSL_Spline *splinecontainer, double x);
double r_of_z(double z);

//====================================================
// All starts with main...
//====================================================

int main(int argc, char **argv){

  //====================================================
  // Read parameters from stdin
  //====================================================
  if(argc < 7){
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
#ifdef OMP
#pragma omp parallel
  {
    if(omp_get_thread_num() == 0) 
      printf("We are using OpenMP nthreads = %i\n", omp_get_num_threads());
  }
#else
  printf("We are NOT using OpenMP\n");
#endif
  printf("====================================\n");
  
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

  return 0;
}

//====================================================
// Spline of r(z) = Int_0^z dz/H(z)
//====================================================
double r_of_z(double z){
  return Lookup_GSL_Spline(global_spline_rofz, z);
}

//====================================================
// Brute-force pair counting
//====================================================
void brute_force_pair_counting(GalaxyCatalog *cat, PairCountBinning *pc){
  // Fetch data from galaxy catalog
  Galaxy *allgalaxies = cat->galaxies;
  int ngalaxies       = cat->ngalaxies;
  
  // Fetch data from binning
  int nbins   = pc->nbins;
  double rmax = pc->rmax;
  double *XX  = pc->paircount;
 
  // Other variables
  double pairs_dist = 0.0, pairs_dist2 = 0.0;
  double nbins_over_rmax = nbins / rmax, rmax2 = rmax * rmax;
  int i;
  
  // Initialize OpenMP
  int nthreads = 1, id = 0;
#ifdef OMP 
#pragma omp parallel
  {
    if(omp_get_thread_num() == 0) nthreads = omp_get_num_threads();
  }
#endif

  // Initialize binning
  double **XX_threads = malloc(sizeof(double *)*nthreads);
  for(id = 0; id < nthreads; id++){
    XX_threads[id] = malloc(sizeof(double)*nbins);
    for(i = 0; i < nbins; i++){
      XX_threads[id][i] = 0.0;
      XX[i] = 0.0;
    }
  }

  printf("\n====================================\n");
  printf("Brute-force pair counting:\n");
  printf("====================================\n");
  printf("Using n = %i threads\n", nthreads);

  // Double loop over all pairs of galaxies
  clock_t start = clock(), diff;
#ifdef OMP
#pragma omp parallel for private(id) reduction(+:pairs_dist2) reduction(+:pairs_dist) schedule(dynamic)
#endif
  for(i = 0; i < ngalaxies; i++){
#ifdef OMP
    id = omp_get_thread_num();
#else
    id = 0;
#endif
    Galaxy *p1 = &allgalaxies[i];
    int j;
    for(j = i+1; j < ngalaxies; j++){
      Galaxy *p2 = &allgalaxies[j];
     
      // Distance between galaxies
      double dist2 = pow2(p1->x[0] - p2->x[0]);
      dist2 += pow2(p1->x[1] - p2->x[1]);
      dist2 += pow2(p1->x[2] - p2->x[2]);

      // Add to bin
      if(dist2 < rmax2){
        double dist = sqrt(dist2);

        // The index in the binning
        int ibin = (int) (dist * nbins_over_rmax);

#ifdef WEIGHTS
        XX_threads[id][ibin] += p1->w * p2->w;
#else
        XX_threads[id][ibin] += 1.0;
#endif
        // Total number of pairs we have computed distances for
        pairs_dist += 1.0;
      }

      // Total number of pairs we have computed square distances for
      pairs_dist2 += 1.0;
    }
  }

  // Sum up results from different threads
  for(id = 0; id < nthreads; id++){
    for(i = 0; i < nbins; i++){
      XX[i] += XX_threads[id][i];
    }
    free(XX_threads[id]);
  }
  free(XX_threads);

  // Show timing and how many dist^2 and square roots we had to compute
  diff = clock() - start;
  printf("We have computed dist^2 for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist2, 
      (double)ngalaxies * (double)ngalaxies/2.0, pairs_dist2/((double)ngalaxies*(double) ngalaxies/2.0) * 100.0);
  printf("We have computed sqrt(dist^2) for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist, 
      (double)ngalaxies * (double)ngalaxies/2.0, pairs_dist/((double)ngalaxies*(double) ngalaxies/2.0) * 100.0);
  printf("This took: [%8.3lf] sec\n", (double)(1000 * diff / CLOCKS_PER_SEC)/1000.0);

#ifdef OUTPUT_PAIRS_TO_SCREEN
  // Output pair counts
  printf("r          XX(r):\n");
  for(i = 0; i < nbins; i++){
    double r = rmax / (double) nbins * (i+0.5);
    double xi = XX[i];
    printf("%lf   %lf\n", r, xi);    
  }
#endif
}

//====================================================
// Brute-force cross pair counting
//====================================================
void brute_force_cross_pair_counting(GalaxyCatalog *cat, GalaxyCatalog *cat2, PairCountBinning *pc){
  // Fetch data from galaxy catalog
  Galaxy *allgalaxies = cat->galaxies;
  int ngalaxies       = cat->ngalaxies;
  
  // Fetch data from galaxy catalog2
  Galaxy *allgalaxies2 = cat2->galaxies;
  int ngalaxies2       = cat2->ngalaxies;
  
  // Fetch data from binning
  int nbins   = pc->nbins;
  double rmax = pc->rmax;
  double *XY  = pc->paircount;
  
  // Other variables
  double pairs_dist = 0.0, pairs_dist2 = 0.0;
  double nbins_over_rmax = nbins / rmax, rmax2 = rmax * rmax;
  int i;
  
  // Initialize OpenMP
  int nthreads = 1, id;
#ifdef OMP 
#pragma omp parallel
  {
    if(omp_get_thread_num() == 0) nthreads = omp_get_num_threads();
  }
#endif

  // Initialize binning
  double **XY_threads = malloc(sizeof(double *)*nthreads);
  for(id = 0; id < nthreads; id++){
    XY_threads[id] = malloc(sizeof(double)*nbins);
    for(i = 0; i < nbins; i++){
      XY_threads[id][i] = 0.0;
      XY[i] = 0.0;
    }
  }

  printf("\n====================================\n");
  printf("Brute-force cross pair counting:\n");
  printf("====================================\n");
  printf("Using n = %i threads\n", nthreads);

  // Double loop over all pairs of galaxies
  clock_t start = clock(), diff;
#ifdef OMP
#pragma omp parallel for private(id) reduction(+:pairs_dist2) reduction(+:pairs_dist) schedule(dynamic)
#endif
  for(i = 0; i < ngalaxies; i++){
#ifdef OMP
    id = omp_get_thread_num();
#else
    id = 0;
#endif
    Galaxy *p1 = &allgalaxies[i];
    int j;
    for(j = 0; j < ngalaxies2; j++){
      Galaxy *p2 = &allgalaxies2[j];
     
      // Distance between galaxies
      double dist2 = pow2(p1->x[0] - p2->x[0]);
      dist2 += pow2(p1->x[1] - p2->x[1]);
      dist2 += pow2(p1->x[2] - p2->x[2]);

      // Add to bin
      if(dist2 < rmax2){
        double dist = sqrt(dist2);

        // The index in the binning
        int ibin = (int) (dist * nbins_over_rmax);

#ifdef WEIGHTS
        XY_threads[id][ibin] += p1->w * p2->w;
#else
        XY_threads[id][ibin] += 1.0;
#endif
        // Total number of pairs we have computed distances for
        pairs_dist += 1.0;
      }

      // Total number of pairs we have computed square distances for
      pairs_dist2 += 1.0;
    }
  }

  // Sum up results from different threads
  for(id = 0; id < nthreads; id++){
    for(i = 0; i < nbins; i++){
      XY[i] += XY_threads[id][i];
    }
    free(XY_threads[id]);
  }
  free(XY_threads);

  // Show timing and how many dist^2 and square roots we had to compute
  diff = clock() - start;
  printf("We have computed dist^2 for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist2, 
      (double)ngalaxies * (double)ngalaxies2, pairs_dist2/((double)ngalaxies*(double) ngalaxies2) * 100.0);
  printf("We have computed sqrt(dist^2) for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist, 
      (double)ngalaxies * (double)ngalaxies2, pairs_dist/((double)ngalaxies*(double) ngalaxies2) * 100.0);
  printf("This took: [%8.3lf] sec\n", (double)(1000 * diff / CLOCKS_PER_SEC)/1000.0);

#ifdef OUTPUT_PAIRS_TO_SCREEN
  // Output cross pair counts
  printf("r          XY(r):\n");
  for(i = 0; i < nbins; i++){
    double r = rmax / (double) nbins * (i+0.5);
    double xi = XY[i];
    printf("%lf   %lf\n", r, xi);    
  }
#endif
}

//====================================================
// Pair counting using grid to speed it up
//====================================================
void grid_pair_counting(Grid *grid, PairCountBinning *pc){
  // Fetch data from grid
  Cell *cells   = grid->cells;
  int ngrid     = grid->ngrid;
  int ngalaxies = grid->ngalaxies;
  int max_ix    = grid->max_ix;
  int max_iy    = grid->max_iy;
  int max_iz    = grid->max_iz;
  double cell_size = grid->cell_size;

  // Fetch data from binning
  int nbins   = pc->nbins;
  double rmax = pc->rmax;
  double *XX  = pc->paircount;

  // Other variables
  double pairs_dist2 = 0.0, pairs_dist = 0.0;
  double nbins_over_rmax = nbins / rmax, rmax2 = rmax * rmax;
  int i;

  // Initialize OpenMP
  int nthreads = 1, id = 0;
#ifdef OMP 
#pragma omp parallel
  {
    if(omp_get_thread_num() == 0) nthreads = omp_get_num_threads();
  }
#endif

  // Initialize binning
  double **XX_threads = malloc(sizeof(double *)*nthreads);
  for(id = 0; id < nthreads; id++){
    XX_threads[id] = malloc(sizeof(double)*nbins);
    for(i = 0; i < nbins; i++){
      XX_threads[id][i] = 0.0;
      XX[i] = 0.0;
    }
  }

  printf("\n====================================\n");
  printf("Pair counting using grid:\n");
  printf("====================================\n");
  printf("Using n = %i threads\n", nthreads);

  // How many cells in each direction we must search
  int delta_ncells = (int)(ceil(rmax / cell_size)) + 1;
  printf("We will go left and right: [%i] cells. The corresponding distance: +/-[%lf] Mpc/h\n", 
      delta_ncells,  delta_ncells * cell_size);

  //==========================================================
  // Loop over all the cells
  // The loops only go up to max_ix etc. since we wan to skip 
  // looping over cells that we know are empty
  //==========================================================

  int ix0, num_processed = 0;
  clock_t start = clock(), diff;
#ifdef OMP 
#pragma omp parallel for private(id) reduction(+:pairs_dist2) reduction(+:pairs_dist) schedule(dynamic)
#endif
  for(ix0 = 0; ix0 <= max_ix; ix0++){
#ifdef OMP 
    id = omp_get_thread_num();
#else
    id = 0;
#endif
    int iy0;
    for(iy0 = 0; iy0 <= max_iy; iy0++){
      int iz0;
      for(iz0 = 0; iz0 <= max_iz; iz0++){

        // Index of current cell
        int index = (ix0*ngrid + iy0)*ngrid + iz0;

        // Pointer to current cell
        Cell *curcell = &cells[index];

        // Number of galaxies in current cell
        int np_cell = curcell->np;

        // Loop over all galaxies in current cell
        int ipart_cell;
        for(ipart_cell = 0; ipart_cell < np_cell; ipart_cell++){

          // Current particle
          Galaxy *curpart_cell = &curcell->galaxy[ipart_cell];

          // We now want to loop over nearby cells by looking at cube of cells around current cell
          int ix_right = ix0 + delta_ncells <= max_ix  ? ix0 + delta_ncells : max_ix;
          int iy_right = iy0 + delta_ncells <= max_iy  ? iy0 + delta_ncells : max_iy;
          int iz_right = iz0 + delta_ncells <= max_iz  ? iz0 + delta_ncells : max_iz;
          int ix_left  = ix0 - delta_ncells >= 0       ? ix0 - delta_ncells : 0;
          int iy_left  = iy0 - delta_ncells >= 0       ? iy0 - delta_ncells : 0;
          int iz_left  = iz0 - delta_ncells >= 0       ? iz0 - delta_ncells : 0;

          // Loop over neightbor cells
          int ix, iy, iz;
          for(ix = ix_left; ix <= ix_right; ix++){
            for(iy = iy_left; iy <= iy_right; iy++){
              for(iz = iz_left; iz <= iz_right; iz++){

                // Avoid double counting so we skip cells that have been correlated with this one before
                if(ix < ix0) continue;
                if(ix == ix0 && iy < iy0) continue;
                if(ix == ix0 && iy == iy0 && iz < iz0) continue;

                // Index of neighboring cell
                int index_neighbor_cell = (ix*ngrid + iy)*ngrid + iz;
                
                // Pointer to neighboring cell
                Cell *neighborcell = &cells[index_neighbor_cell];

                // Number of galaxies in neighboring cell
                int npart_neighbor_cell = neighborcell->np;

                // Careful: if the nbor cell is the same as the current cell then 
                // we will overcount if we do all particles so only correlate with partices we haven't touched yet 
                int istart_nbor_cell = 0;
                if(curcell == neighborcell) istart_nbor_cell = ipart_cell + 1;

                // Loop over galaxies in neighbor cells
                int ipart_neighbor_cell;
                for(ipart_neighbor_cell = istart_nbor_cell; ipart_neighbor_cell < npart_neighbor_cell; ipart_neighbor_cell++){
                  
                  // Galaxy in neighboring cell
                  Galaxy *curpart_neighbor_cell = &neighborcell->galaxy[ipart_neighbor_cell];

                  // ==================================================================
                  // We now count up the pair [curpart_cell] x [curpart_neighbor_cell]
                  // ==================================================================

                  // The distance between the two galaxies
                  double dist2 = pow2(curpart_cell->x[0] - curpart_neighbor_cell->x[0]);
                  dist2 += pow2(curpart_cell->x[1] - curpart_neighbor_cell->x[1]);
                  dist2 += pow2(curpart_cell->x[2] - curpart_neighbor_cell->x[2]);

                  // Add to bin
                  if(dist2 < rmax2){
                    double dist = sqrt(dist2);

                    // The index in the binning
                    int ibin = (int) (dist * nbins_over_rmax);

#ifdef WEIGHTS
                    XX_threads[id][ibin] += curpart_cell->w * curpart_neighbor_cell->w;
#else
                    XX_threads[id][ibin] += 1.0;
#endif

                    // Total number of pairs we have computed distances for
                    pairs_dist += 1.0;
                  }
                  
                  // Total number of pairs we have computed square distances for
                  pairs_dist2 += 1.0;
                }
              }
            }
          }
        }
      }
    }

    // Show progress...
#pragma omp critical
    {
      printf("Processed (%4i / %4i)   (ThreadID: %3i)\n", num_processed, max_ix, id);
      num_processed++;
    }
  }

  // Sum up results from different threads
  for(id = 0; id < nthreads; id++){
    for(i = 0; i < nbins; i++){
      XX[i] += XX_threads[id][i];
    }
    free(XX_threads[id]);
  }
  free(XX_threads);

  // Show timing and how many dist^2 and square roots we had to compute
  diff = clock() - start;
  printf("We have computed dist^2 for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist2, 
      (double)ngalaxies * (double)ngalaxies/2.0, pairs_dist2/((double)ngalaxies*(double) ngalaxies/2.0) * 100.0);
  printf("We have computed sqrt(dist^2) for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist, 
      (double)ngalaxies * (double)ngalaxies/2.0, pairs_dist/((double)ngalaxies*(double) ngalaxies/2.0) * 100.0);
  printf("This took: [%8.3lf] sec\n", (double)(1000 * diff / CLOCKS_PER_SEC)/1000.0);

#ifdef OUTPUT_PAIRS_TO_SCREEN
  // Output pair counts
  printf("r          XX(r):\n");
  for(i = 0; i < nbins; i++){
    double r = rmax / (double) nbins * (i+0.5);
    double xi = XX[i];
    printf("%lf   %lf\n", r, xi);    
  }
#endif
}

//====================================================
// Read galaxies from file. Format: RA DEC z
//====================================================
GalaxyCatalog *read_galaxies_from_file(char *filename, int ngalaxies){
  int i;

  printf("\n====================================\n"); 
  printf("Reading from file [%s]\n", filename); 
  printf("====================================\n"); 
  printf("Galaxy file has %d galaxies\n",ngalaxies);

  // Allocate particle array
  GalaxyCatalog *cat = malloc(sizeof(GalaxyCatalog));
  cat->galaxies = malloc(sizeof(Galaxy)*ngalaxies);
  Galaxy *allgalaxies = cat->galaxies;
  cat->ngalaxies = ngalaxies;
  cat->allocated = 1;
    
  double sum_w = 0.0, sum_w2 = 0.0;

  // Read the data
  FILE *fp = fopen(filename, "r");
  for(i = 0; i < ngalaxies; i++){
    char line[1024];
    double Pos[3], RA, DEC, z, w = 1.0;

    // Read galaxy
    fgets(line, sizeof(line),fp);
#ifdef WEIGHTS
    sscanf(line, "%lf %lf %lf %lf", &RA, &DEC, &z, &w);
#else
    sscanf(line, "%lf %lf %lf", &RA, &DEC, &z);
#endif

    // Convert to positions in Mpc/h
    double costheta = cos(2.0 * M_PI / 360.0 * (90.0 - DEC));
    double sintheta = sqrt(1.0 - costheta*costheta);
    double phi      = RA * 2.0 * M_PI / 360.0;
    double r        = r_of_z(z);

    // Set positions
    Pos[0] = r * sintheta * cos(phi);
    Pos[1] = r * sintheta * sin(phi);
    Pos[2] = r * costheta;

    // Store galaxies
    allgalaxies[i].x[0] = Pos[0];
    allgalaxies[i].x[1] = Pos[1];
    allgalaxies[i].x[2] = Pos[2];
#ifdef WEIGHTS
    allgalaxies[i].w = w;
#endif  
    sum_w  += w;
    sum_w2 += w*w;
  }
  fclose(fp);
  
  // The mean weight and RMS
  printf("Mean weight: %lf  RMS: %lf\n", sum_w/(double)ngalaxies, sqrt(sum_w2/(double)ngalaxies));
  cat->sum_w = sum_w;
  cat->sum_w2 = sum_w2;

  return cat;
}

//====================================================
// Compute the boxsize we need to encompas all galaxies
// in both catalogs and shift the positions so that 
// they are inside [0, BOX]^3
//====================================================
void compute_boxsize_shift_positions(GalaxyCatalog *cat, GalaxyCatalog *cat2, double *box){
  int ngalaxies  = cat->ngalaxies;
  Galaxy *galaxies = cat->galaxies;
  int ngalaxies2 = cat2->ngalaxies;
  Galaxy *galaxies2 = cat2->galaxies;
  
  // Compute max and min position
  double max_x = -1e100, min_x = 1e100;
  double max_y = -1e100, min_y = 1e100;
  double max_z = -1e100, min_z = 1e100;
 
  int i;
  for(i = 0; i < ngalaxies; i++){
    double *Pos = &galaxies[i].x[0];
    if(Pos[0] > max_x) max_x = Pos[0];
    if(Pos[1] > max_y) max_y = Pos[1];
    if(Pos[2] > max_z) max_z = Pos[2];
    if(Pos[0] < min_x) min_x = Pos[0];
    if(Pos[1] < min_y) min_y = Pos[1];
    if(Pos[2] < min_z) min_z = Pos[2];
  }
  for(i = 0; i < ngalaxies2; i++){
    double *Pos = &galaxies2[i].x[0];
    if(Pos[0] > max_x) max_x = Pos[0];
    if(Pos[1] > max_y) max_y = Pos[1];
    if(Pos[2] > max_z) max_z = Pos[2];
    if(Pos[0] < min_x) min_x = Pos[0];
    if(Pos[1] < min_y) min_y = Pos[1];
    if(Pos[2] < min_z) min_z = Pos[2];
  }

  // Shift positions
  for(i = 0; i < ngalaxies; i++){
    double *Pos = &galaxies[i].x[0];
    Pos[0] -= min_x;
    Pos[1] -= min_y;
    Pos[2] -= min_z;
  }
  for(i = 0; i < ngalaxies2; i++){
    double *Pos = &galaxies2[i].x[0];
    Pos[0] -= min_x;
    Pos[1] -= min_y;
    Pos[2] -= min_z;
  }

  // New max_x,y,z values
  max_x = max_x-min_x;
  max_y = max_y-min_y;
  max_z = max_z-min_z;

  // The min/max positions (separations)
  printf("\n====================================\n");
  printf("Shifting particles and computing boxsize:\n");
  printf("====================================\n");
  printf("Min/Max X position: 0.0 -> %5.2lf Mpc/h\n", max_x);
  printf("Min/Max Y position: 0.0 -> %5.2lf Mpc/h\n", max_y);
  printf("Min/Max Z position: 0.0 -> %5.2lf Mpc/h\n", max_z);
  
  // The largest displacement in either direction
  *box = 1.01 * MAX(max_x, MAX(max_y, max_z));
  printf("The boxsize we will use is %5.2lf Mpc/h\n", *box);
}

//====================================================
// Function to create a new grid
// Estimates the gridsize from rmax, box and nparticles
// and initalizes all th cells
//====================================================
Grid *create_grid(int ngalaxies, double rmax, double box){
  
  //====================================================
  // Estimate optimal size of grid
  //====================================================
  int ngrid1 = (int)(8.0*box/rmax);                  // 8 cells to get to rmax
  int ngrid2 = (int)(pow(0.5*ngalaxies,0.3333333));  // 2 galaxies per cells on average
  int ngrid  = (int) ceil( MIN( ngrid1, ngrid2 ) );

  //====================================================
  // Make our grid. This is an array of N^3 cells
  //====================================================
  int NcellTotal  = pow3(ngrid);
  Grid *grid      = malloc(sizeof(Grid));
  grid->cells     = malloc(sizeof(Cell) * NcellTotal);
  Cell *cells     = grid->cells;
  grid->ngrid     = ngrid;
  grid->ngalaxies = ngalaxies;
  grid->cell_size = box/(double) ngrid;
  grid->allocated = 1;

  printf("\n====================================\n"); 
  printf("Creating new grid\n"); 
  printf("====================================\n"); 
  printf("ngrid = %i [min(%i,%i)] CellSize: %lf Mpc/h\n", ngrid, ngrid1, ngrid2, box/(double) ngrid);
  printf("Total number of cells: [%i] Galaxies that will be added to grid: [%i]\n", NcellTotal, ngalaxies);

  //====================================================
  // Initialize the number of galaxies in all cells
  //====================================================
  int i;
  for(i = 0; i < NcellTotal; i++){
    cells[i].np = 0;
  }

  return grid;
}

//====================================================
// Count the number of lines in a file
//====================================================
int count_lines_in_file(char *filename){
  int n = 0;
  FILE *fp = fopen(filename, "r");
  while(!feof(fp)){
    char ch = fgetc(fp);
    if(ch == '\n') n++;
  }
  return n;
}

//====================================================
// Add a list of galaxies to a grid. We first
// count how many are in each cell, allocate and then add
//====================================================
void add_galaxies_to_cells(Grid *grid, GalaxyCatalog *cat){
  // Fetch data from grid
  Cell *cells = grid->cells;
  int ngrid   = grid->ngrid;
  double cell_size = grid->cell_size;
  
  // Fetch data from galaxy catalog
  Galaxy *allgalaxies = cat->galaxies;
  int ngalaxies       = cat->ngalaxies;
  
  // Other variables
  int NcellTotal = pow3(ngrid), i;
  int max_ix = 0, max_iy = 0, max_iz = 0;

  //====================================================
  // Count how many particles there are in each cell
  //====================================================
  for(i = 0; i < ngalaxies; i++){
    // Position of current galaxy
    double *Pos = &allgalaxies[i].x[0];

    // Determine the cell-coordinate the galaxies belongs to
    int ix = (int)(Pos[0]/cell_size);
    int iy = (int)(Pos[1]/cell_size);
    int iz = (int)(Pos[2]/cell_size);

    // Check for errors
    if(ix >= ngrid || iy >= ngrid || iz >= ngrid){
      printf("Error in add_galaxies_to_grid (%i,%i,%i)  ngrid = %i\n", ix, iy, iz, ngrid);
      exit(1);
    }

    // Find the maximum values (so we can skip looping through them later)
    if(ix > max_ix) max_ix = ix;
    if(iy > max_iy) max_iy = iy;
    if(iz > max_iz) max_iz = iz;

    // Index of cell particle belongs to in grid
    int index = (ix*ngrid + iy)*ngrid + iz;

    // Increase the number of galaxies we have in the cell
    cells[index].np += 1; 
  }

  // Max grid coordinates where we have galaxies
  max_ix = MIN(max_ix,ngrid-1);
  max_iy = MIN(max_iy,ngrid-1);
  max_iz = MIN(max_iz,ngrid-1);
  grid->max_ix = max_ix;
  grid->max_iy = max_iy;
  grid->max_iz = max_iz;

  //====================================================
  // Now that we know how many galaxies are in each cell
  // Allocate the particle array in each cell
  //====================================================
  int nempty = 0;
  for(i = 0; i < NcellTotal; i++){
    grid->cells[i].galaxy = malloc(sizeof(Galaxy) * cells[i].np);
    if(cells[i].np == 0) nempty += 1;
  }

  printf("\n====================================\n");
  printf("Adding galaxies to grid\n");
  printf("====================================\n");
  printf("There are [%i] empty cells (%lf %%) in the grid\n", nempty, (double)nempty / (double) NcellTotal * 100.0);

  //====================================================
  // Go through galaxies and add them to the grid cells
  // where they belong
  //====================================================

  // Book-keeping array when adding particles (tells how many we have added in each cell so far)
  int *ncurrent = malloc(sizeof(int)*NcellTotal);
  for(i = 0; i < NcellTotal; i++) ncurrent[i] = 0;

  for(i = 0; i < ngalaxies; i++){

    // Current particle position
    double *Pos = &allgalaxies[i].x[0];

    // Determine the cell-coordinate the galaxies belongs to
    int ix = (int)(Pos[0]/cell_size);
    int iy = (int)(Pos[1]/cell_size);
    int iz = (int)(Pos[2]/cell_size);

    // Index of cell particle belongs to in grid
    int index = (ix*ngrid + iy)*ngrid + iz;

    // Assign particle to the array in the given cell
    int partnumber = ncurrent[index];
    cells[index].galaxy[partnumber].x[0] = Pos[0];
    cells[index].galaxy[partnumber].x[1] = Pos[1];
    cells[index].galaxy[partnumber].x[2] = Pos[2];
#ifdef WEIGHTS
    cells[index].galaxy[partnumber].w = allgalaxies[i].w;
#endif

    // Increase book-keeping variable 
    ncurrent[index]++;
  }
  free(ncurrent);
}

//====================================================
// Cross pair counts using grid to speed it up
// Cross seems to be faster if we loop over the coarsest 
// grid first so call in order (galaxy_grid, random_grid)
//====================================================
void grid_cross_pair_counting(Grid *grid, Grid *grid2, PairCountBinning *pc){
  // Fetch data from the grid
  Cell *cells    = grid->cells;
  int ngrid      = grid->ngrid;
  int ngalaxies  = grid->ngalaxies;
  int max_ix     = grid->max_ix;
  int max_iy     = grid->max_iy;
  int max_iz     = grid->max_iz;

  // Fetch data from the grid2
  Cell *cells2   = grid2->cells;
  int ngrid2     = grid2->ngrid;
  int ngalaxies2 = grid2->ngalaxies;
  int max_ix2    = grid2->max_ix;
  int max_iy2    = grid2->max_iy;
  int max_iz2    = grid2->max_iz;
  double cell_size2 = grid2->cell_size;

  // Fetch data from the binning
  int nbins      = pc->nbins;
  double rmax    = pc->rmax;
  double *XY     = pc->paircount;

  // Other variables
  double pairs_dist2 = 0.0, pairs_dist = 0.0;
  double nbins_over_rmax = nbins / rmax, rmax2 = rmax * rmax;
  int i;

  // Initialize OpenMP
  int nthreads = 1, id = 0;
#ifdef OMP 
#pragma omp parallel
  {
    if(omp_get_thread_num() == 0) nthreads = omp_get_num_threads();
  }
#endif

  // Initialize binning
  double **XY_threads = malloc(sizeof(double *)*nthreads);
  for(id = 0; id < nthreads; id++){
    XY_threads[id] = malloc(sizeof(double)*nbins);
    for(i = 0; i < nbins; i++){
      XY_threads[id][i] = 0.0;
      XY[i] = 0.0;
    }
  }

  printf("\n====================================\n");
  printf("Cross pair counts using grid:\n");
  printf("====================================\n");
  printf("Using n = %i threads\n", nthreads);

  // How many cells in each direction we must search in the second grid
  int delta_ncells2 = (int)(ceil(rmax / cell_size2)) + 2;
  printf("We will go left and right: [%i] cells. The corresponding distance: +/-[%lf] Mpc/h\n", 
      delta_ncells2,  delta_ncells2 * cell_size2);

  //==========================================================
  // Loop over all the cells in grid1
  // The loops only go up to max_ix etc. since we wan to skip 
  // looping over cells that we know are empty
  //==========================================================

  int ix0, num_processed = 0;
  clock_t start = clock(), diff;
#ifdef OMP 
#pragma omp parallel for private(id) reduction(+:pairs_dist2) reduction(+:pairs_dist) schedule(dynamic)
#endif
  for(ix0 = 0; ix0 <= max_ix; ix0++){
#ifdef OMP 
    id = omp_get_thread_num();
#else
    id = 0;
#endif
    int iy0;
    for(iy0 = 0; iy0 <= max_iy; iy0++){
      int iz0;
      for(iz0 = 0; iz0 <= max_iz; iz0++){

        // Index of current cell
        int index = (ix0*ngrid + iy0)*ngrid + iz0;

        // Pointer to current cell
        Cell *curcell = &cells[index];

        // Number of galaxies in current cell
        int np_cell = curcell->np;

        // Loop over all galaxies in current cell
        int ipart_cell;
        for(ipart_cell = 0; ipart_cell < np_cell; ipart_cell++){

          // Current particle
          Galaxy *curpart_cell = &curcell->galaxy[ipart_cell];

          //========================================
          // Now we find the index of the second grid 
          // this grid corresponds to
          //========================================
          int ix_grid2 = (int)(ix0 * ngrid2/(double)ngrid);
          int iy_grid2 = (int)(iy0 * ngrid2/(double)ngrid);
          int iz_grid2 = (int)(iz0 * ngrid2/(double)ngrid);

          // We now want to loop over nearby cells by looking at cube of cells around current cell
          int ix2_right = ix_grid2 + delta_ncells2 <= max_ix2  ? ix_grid2 + delta_ncells2 : max_ix2;
          int iy2_right = iy_grid2 + delta_ncells2 <= max_iy2  ? iy_grid2 + delta_ncells2 : max_iy2;
          int iz2_right = iz_grid2 + delta_ncells2 <= max_iz2  ? iz_grid2 + delta_ncells2 : max_iz2;
          int ix2_left  = ix_grid2 - delta_ncells2 >= 0        ? ix_grid2 - delta_ncells2 : 0;
          int iy2_left  = iy_grid2 - delta_ncells2 >= 0        ? iy_grid2 - delta_ncells2 : 0;
          int iz2_left  = iz_grid2 - delta_ncells2 >= 0        ? iz_grid2 - delta_ncells2 : 0;

          // Loop over neightbor cells
          int ix2, iy2, iz2;
          for(ix2 = ix2_left; ix2 <= ix2_right; ix2++){
            for(iy2 = iy2_left; iy2 <= iy2_right; iy2++){
              for(iz2 = iz2_left; iz2 <= iz2_right; iz2++){

                // Index of neighboring cell
                int index_neighbor_cell = (ix2*ngrid2 + iy2)*ngrid2 + iz2;
                
                // Pointer to neighboring cell
                Cell *neighborcell = &cells2[index_neighbor_cell];

                // Number of galaxies in neighboring cell
                int npart_neighbor_cell = neighborcell->np;

                // Loop over galaxies in neighbor cells
                int ipart_neighbor_cell;
                for(ipart_neighbor_cell = 0; ipart_neighbor_cell < npart_neighbor_cell; ipart_neighbor_cell++){
                  
                  // Galaxy in neighboring cell
                  Galaxy *curpart_neighbor_cell = &neighborcell->galaxy[ipart_neighbor_cell];

                  // ==================================================================
                  // We now count up the pair [curpart_cell] x [curpart_neighbor_cell]
                  // ==================================================================

                  // The distance between the two galaxies
                  double dist2 = pow2(curpart_cell->x[0] - curpart_neighbor_cell->x[0]);
                  dist2 += pow2(curpart_cell->x[1] - curpart_neighbor_cell->x[1]);
                  dist2 += pow2(curpart_cell->x[2] - curpart_neighbor_cell->x[2]);

                  // Add to bin
                  if(dist2 < rmax2){
                    double dist = sqrt(dist2);

                    // The index in the binning
                    int ibin = (int) (dist * nbins_over_rmax);

#ifdef WEIGHTS
                    XY_threads[id][ibin] += curpart_cell->w * curpart_neighbor_cell->w;
#else
                    XY_threads[id][ibin] += 1.0;
#endif

                    // Total number of pairs we have computed distances for
                    pairs_dist += 1.0;
                  }
                  
                  // Total number of pairs we have computed square distances for
                  pairs_dist2 += 1.0;
                }
              }
            }
          }
        }
      }
    }
    
    // Show progress...
#pragma omp critical
    {
      printf("Processed (%4i / %4i)   (ThreadID: %3i)\n", num_processed, max_ix, id);
      num_processed++;
    }
  }

  // Sum up results from different threads
  for(id = 0; id < nthreads; id++){
    for(i = 0; i < nbins; i++){
      XY[i] += XY_threads[id][i];
    }
    free(XY_threads[id]);
  }
  free(XY_threads);

  // Show timing and how many dist^2 and square roots we had to compute
  diff = clock() - start;
  printf("We have computed dist^2 for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist2, 
      (double)ngalaxies * (double)ngalaxies2, pairs_dist2/((double)ngalaxies*(double) ngalaxies2) * 100.0);
  printf("We have computed sqrt(dist^2) for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist, 
      (double)ngalaxies * (double)ngalaxies2, pairs_dist/((double)ngalaxies*(double) ngalaxies2) * 100.0);
  printf("This took: [%8.3lf] sec\n", (double)(1000 * diff / CLOCKS_PER_SEC)/1000.0);

#ifdef OUTPUT_PAIRS_TO_SCREEN
  // Output cross pair counts
  printf("r          XY(r):\n");
  for(i = 0; i < nbins; i++){
    double r = rmax / (double) nbins * (i+0.5);
    double xi = XY[i];
    printf("%lf   %lf\n", r, xi);    
  }
#endif
}

//====================================================
// Compute correlation function using LS
//====================================================
void compute_correlation_function(PairCountBinning *DD, PairCountBinning *DR, PairCountBinning *RR, char *filename){
  // Fetch data from binning
  int nbins      = DD->nbins;
  double rmax    = DD->rmax;
  double norm_DD = DD->norm;
  double norm_DR = DR->norm;
  double norm_RR = RR->norm;
  int i;

  printf("\n====================================\n");
  printf("Correlation function using LS estimator:\n");
  printf("====================================\n");
  printf("Outputfile [%s] has format [r  xi  err_xi  DD  DR  RR]\n", filename);
  
  FILE *fp = fopen(filename, "w");
  for(i = 0; i < nbins; i++){
    // Center of bin
    double r = rmax / (double) nbins * (i+0.5);
    
    // Compute correlation function with Poisson errors
    double corr = 0.0, err_corr = 0.0;
    if(DD->paircount[i] != 0.0){
      
      double one_over_sqrtDD = 1.0/sqrt(DD->paircount[i]);
      double one_over_sqrtDR = 1.0/sqrt(DR->paircount[i]);
      double one_over_sqrtRR = 1.0/sqrt(RR->paircount[i]);

      double normed_DD = DD->paircount[i] / norm_DD;
      double normed_DR = DR->paircount[i] / norm_DR;
      double normed_RR = RR->paircount[i] / norm_RR;
      
      corr = (normed_DD - 2.0*normed_DR + normed_RR) / normed_RR;
      err_corr = (1.0 + corr) * (one_over_sqrtDD + one_over_sqrtDR + one_over_sqrtRR);
    }

    // Output to file
#ifdef WEIGHTS
    fprintf(fp,"%le  %le  %le  %le  %le  %le\n", r, corr, err_corr, DD->paircount[i], DR->paircount[i], RR->paircount[i]);
#else
    fprintf(fp,"%le  %le  %le  %d   %d   %d\n",  r, corr, err_corr, (int)DD->paircount[i], (int)DR->paircount[i], (int)RR->paircount[i]);
#endif
  }
  fclose(fp);
}

//====================================================
// Allocate memory for a binning
//====================================================
PairCountBinning *create_binning(int nbins, double rmax){
  PairCountBinning *pc = malloc(sizeof(PairCountBinning));
  pc->paircount = malloc(sizeof(double) * nbins);
  pc->nbins = nbins;
  pc->rmax  = rmax;
  pc->norm  = 1.0;
  pc->allocated = 1;
  return pc;
}

//====================================================
// Free memory associated with a binning
//====================================================
void free_binning(PairCountBinning *pc){
  if(pc->allocated){
    free(pc->paircount);
    free(pc);
  }
}

//====================================================
// Free the memory associated with a galaxy catalog
//====================================================
void free_cat(GalaxyCatalog *cat){
  if(cat->allocated){
    free(cat->galaxies);
    free(cat);
  }
}

//====================================================
// Free the memory associated with a grid
//====================================================
void free_grid(Grid *grid){
  if(grid->allocated){
    int ngrid = grid->ngrid;
    int NcellTotal = pow3(ngrid), i;
    for(i = 0; i < NcellTotal; i++){
      free(grid->cells[i].galaxy);
    }
    free(grid->cells);
    free(grid);
  }
}

//====================================================
// Output galaxy catalog in physical coordinates
//====================================================
void outputGalaxies(GalaxyCatalog *cat, char *filename){
  int ngalaxies = cat->ngalaxies;
  Galaxy *galaxies = cat->galaxies;
  FILE *fp = fopen(filename, "w");

  int i;
  for(i = 0; i < ngalaxies; i++){
    Galaxy *curgalaxy = &galaxies[i];
#ifdef WEIGHTS
    fprintf(fp, "%lf  %lf  %lf  %lf\n", curgalaxy->x[0], curgalaxy->x[1], curgalaxy->x[2], curgalaxy->w);
#else
    fprintf(fp, "%lf  %lf  %lf\n", curgalaxy->x[0], curgalaxy->x[1], curgalaxy->x[2]);
#endif
  }
  fclose(fp);
}

//====================================================
// The ODE dr/dz = 1/H(z) in units of Mpc/h for LCDM
// Hubble_Length = c/H0
//====================================================
int ode_rofz(double z, const double r[], double drdz[], void *params){
  const double Hubble_Length_in_Mpch = SPEED_OF_LIGHT_IN_KM_PER_SEC / 100.0;
  double OmegaM = *(double *) params;
  
  drdz[0] =  Hubble_Length_in_Mpch / sqrt(OmegaM*pow3(1.0+z) + 1.0 - OmegaM);
  return GSL_SUCCESS;
}

//====================================================
// Integrate the ODE for r(z) and return a spline of the result
//====================================================
GSL_Spline* create_rofz_spline(double OmegaM){
  const int npts      = 1000;
  const double zmax   = 5.0;
  const double deltaz = zmax/(double) (npts-1);

  // Set up ODE system
  gsl_odeiv2_system sys_rofz = {ode_rofz, NULL, 1, &OmegaM};
  gsl_odeiv2_driver *ode = gsl_odeiv2_driver_alloc_y_new (&sys_rofz, gsl_odeiv2_step_rk2, 1e-10, 1e-10, 0.0);

  double *z_arr = malloc(sizeof(double) * npts);
  double *rofz_arr = malloc(sizeof(double) * npts);

  // Initial conditions r(z=0) = 0
  double ode_z = z_arr[0] = 0.0;
  double rofz_now[1] = {0.0};
  rofz_arr[0] = 0.0;

  // Integration over redshift
  int i;
  for(i = 1; i < npts; i++){
    double znow = i * deltaz;

    // Integrate one step
    int status = gsl_odeiv2_driver_apply(ode, &ode_z, znow, rofz_now);
    if(status != GSL_SUCCESS){
      printf("Error in integrating z = %f  r = %f\n", znow, rofz_now[0]);
      exit(1);
    }

    // Store values
    z_arr[i] = znow;
    rofz_arr[i] = rofz_now[0];
  }

  // Spline up the results
  GSL_Spline *rofz_spline = malloc(sizeof(GSL_Spline));
  Create_GSL_Spline(rofz_spline, z_arr, rofz_arr, npts);

  // Free memory
  free(z_arr);
  free(rofz_arr);

  return rofz_spline;
}

//====================================================
// Create a GSL spline
//====================================================
void Create_GSL_Spline(GSL_Spline *splinecontainer, double *x, double *y, int nx){
  splinecontainer->xacc   = gsl_interp_accel_alloc();
  splinecontainer->xmin   = x[0];
  splinecontainer->xmax   = x[nx-1];
  splinecontainer->spline = gsl_spline_alloc(gsl_interp_cspline, nx);
  gsl_spline_init(splinecontainer->spline, x, y, nx);
  splinecontainer->allocated = 1;
}

//====================================================
// Free memory of a GSL spline
//====================================================
void Free_GSL_Spline(GSL_Spline *splinecontainer){
  if(splinecontainer->allocated){
    gsl_interp_accel_free(splinecontainer->xacc);
    gsl_spline_free(splinecontainer->spline);
    splinecontainer->allocated = 0;
    free(splinecontainer);
  }
}

//====================================================
// Lookup a value from a GSL spline
//====================================================
double Lookup_GSL_Spline(GSL_Spline *splinecontainer, double x){
  double xx = x;
  if(x < splinecontainer->xmin) xx = splinecontainer->xmin;
  if(x > splinecontainer->xmax) xx = splinecontainer->xmax;
  return gsl_spline_eval(splinecontainer->spline, xx, splinecontainer->xacc);
}

