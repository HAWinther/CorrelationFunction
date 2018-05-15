#include "correlators.h"

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
#ifdef USE_OMP 
#pragma omp parallel
  {
    if(omp_get_thread_num() == 0) nthreads = omp_get_num_threads();
  }
#elif defined(USE_MPI)
  nthreads = mpi_size;
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

  if(mpi_rank == 0){
    printf("\n====================================\n");
    printf("Brute-force pair counting:\n");
    printf("====================================\n");
    printf("Using n = %i threads\n", nthreads);
  }

  // Double loop over all pairs of galaxies
  clock_t start = clock(), diff;
  int istart = 0, iend = ngalaxies;
#ifdef USE_OMP
#pragma omp parallel for private(id) reduction(+:pairs_dist2) reduction(+:pairs_dist) schedule(dynamic)
#elif defined(USE_MPI)
  int i_per_task = ngalaxies / nthreads;
  istart = i_per_task * mpi_rank;
  iend   = MIN( i_per_task * (mpi_rank + 1), ngalaxies);
#endif
  for(i = istart; i < iend; i++){
#ifdef USE_OMP
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

#ifdef USE_MPI
  // Gather results from all CPUs
  double *recv = malloc(sizeof(double)*nbins*nthreads);
  MPI_Allgather(XX_threads[mpi_rank], nbins, MPI_DOUBLE, recv, nbins, MPI_DOUBLE, MPI_COMM_WORLD);
  for(id = 0; id < nthreads; id++)
    for(i = 0; i < nbins; i++)
      XX_threads[id][i] = recv[id * nbins + i];
  free(recv);
#endif

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
  if(mpi_rank == 0){
    printf("We have computed dist^2 for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist2, 
        (double)ngalaxies * (double)ngalaxies/2.0, pairs_dist2/((double)ngalaxies*(double) ngalaxies/2.0) * 100.0);
    printf("We have computed sqrt(dist^2) for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist, 
        (double)ngalaxies * (double)ngalaxies/2.0, pairs_dist/((double)ngalaxies*(double) ngalaxies/2.0) * 100.0);
    printf("This took: [%8.3lf] sec\n", (double)(1000 * diff / CLOCKS_PER_SEC)/1000.0);
  }

#ifdef OUTPUT_PAIRS_TO_SCREEN
  // Output pair counts
  if(mpi_rank == 0)
    printf("r          XX(r):\n");
  for(i = 0; i < nbins; i++){
    double r = rmax / (double) nbins * (i+0.5);
    double xi = XX[i];
    if(mpi_rank == 0)
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
  int nthreads = 1, id = 0;
#ifdef USE_OMP 
#pragma omp parallel
  {
    if(omp_get_thread_num() == 0) nthreads = omp_get_num_threads();
  }
#elif defined(USE_MPI)
  nthreads = mpi_size;
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

  if(mpi_rank == 0){
    printf("\n====================================\n");
    printf("Brute-force cross pair counting:\n");
    printf("====================================\n");
    printf("Using n = %i threads\n", nthreads);
  }

  // Double loop over all pairs of galaxies
  clock_t start = clock(), diff;
  int istart = 0, iend = ngalaxies;
#ifdef USE_OMP
#pragma omp parallel for private(id) reduction(+:pairs_dist2) reduction(+:pairs_dist) schedule(dynamic)
#elif defined(USE_MPI)
  int i_per_task = ngalaxies / nthreads;
  istart = i_per_task * mpi_rank;
  iend   = MIN( i_per_task * (mpi_rank + 1), ngalaxies);
#endif
  for(i = istart; i < iend; i++){
#ifdef USE_OMP
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

#ifdef USE_MPI
  // Gather results from all CPUs
  double *recv = malloc(sizeof(double)*nbins*nthreads);
  MPI_Allgather(XY_threads[mpi_rank], nbins, MPI_DOUBLE, recv, nbins, MPI_DOUBLE, MPI_COMM_WORLD);
  for(id = 0; id < nthreads; id++)
    for(i = 0; i < nbins; i++)
      XY_threads[id][i] = recv[id * nbins + i];
  free(recv);
#endif

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
  if(mpi_rank == 0){
    printf("We have computed dist^2 for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist2, 
        (double)ngalaxies * (double)ngalaxies2, pairs_dist2/((double)ngalaxies*(double) ngalaxies2) * 100.0);
    printf("We have computed sqrt(dist^2) for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist, 
        (double)ngalaxies * (double)ngalaxies2, pairs_dist/((double)ngalaxies*(double) ngalaxies2) * 100.0);
    printf("This took: [%8.3lf] sec\n", (double)(1000 * diff / CLOCKS_PER_SEC)/1000.0);
  }

#ifdef OUTPUT_PAIRS_TO_SCREEN
  // Output cross pair counts
  if(mpi_rank == 0)
    printf("r          XY(r):\n");
  for(i = 0; i < nbins; i++){
    double r = rmax / (double) nbins * (i+0.5);
    double xi = XY[i];
    if(mpi_rank == 0)
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
#ifdef USE_OMP 
#pragma omp parallel
  {
    if(omp_get_thread_num() == 0) nthreads = omp_get_num_threads();
  }
#elif defined(USE_MPI)
  nthreads = mpi_size;
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

  if(mpi_rank == 0){
    printf("\n====================================\n");
    printf("Pair counting using grid:\n");
    printf("====================================\n");
    printf("Using n = %i threads\n", nthreads);
  }

  // How many cells in each direction we must search
  int delta_ncells = (int)(ceil(rmax / cell_size)) + 1;
  if(mpi_rank == 0)
    printf("We will go left and right: [%i] cells. The corresponding distance: +/-[%lf] Mpc/h\n", 
        delta_ncells,  delta_ncells * cell_size);

  //==========================================================
  // Loop over all the cells
  // The loops only go up to max_ix etc. since we wan to skip 
  // looping over cells that we know are empty
  //==========================================================

  int ix0, num_processed = 0;
  clock_t start = clock(), diff;
  int istart = 0, iend = max_ix + 1;
#ifdef USE_OMP 
#pragma omp parallel for private(id) reduction(+:pairs_dist2) reduction(+:pairs_dist) schedule(dynamic)
#elif defined(USE_MPI)
  int i_per_task = (max_ix + 1) / nthreads;
  istart = i_per_task * mpi_rank;
  iend   = i_per_task * (mpi_rank + 1);
  if(mpi_rank == mpi_size - 1) iend = max_ix + 1;
#endif
  for(ix0 = istart; ix0 < iend; ix0++){
#ifdef USE_OMP 
    id = omp_get_thread_num();
#elif defined(USE_MPI)
    id = mpi_rank;
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
#ifdef USE_OMP
#pragma omp critical
    {
      printf("Processed (%4i / %4i)   (ThreadID: %3i)\n", num_processed, max_ix, id);
      num_processed++;
    }
#endif
  }

#ifdef USE_MPI
  // Gather results from all CPUs
  double *recv = malloc(sizeof(double)*nbins*nthreads);
  MPI_Allgather(XX_threads[mpi_rank], nbins, MPI_DOUBLE, recv, nbins, MPI_DOUBLE, MPI_COMM_WORLD);
  for(id = 0; id < nthreads; id++)
    for(i = 0; i < nbins; i++)
      XX_threads[id][i] = recv[id * nbins + i];
  free(recv);
#endif

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
  if(mpi_rank == 0){
    printf("We have computed dist^2 for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist2, 
        (double)ngalaxies * (double)ngalaxies/2.0, pairs_dist2/((double)ngalaxies*(double) ngalaxies/2.0) * 100.0);
    printf("We have computed sqrt(dist^2) for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist, 
        (double)ngalaxies * (double)ngalaxies/2.0, pairs_dist/((double)ngalaxies*(double) ngalaxies/2.0) * 100.0);
    printf("This took: [%8.3lf] sec\n", (double)(1000 * diff / CLOCKS_PER_SEC)/1000.0);
  }

#ifdef OUTPUT_PAIRS_TO_SCREEN
  // Output pair counts
  if(mpi_rank == 0)
    printf("r          XX(r):\n");
  for(i = 0; i < nbins; i++){
    double r = rmax / (double) nbins * (i+0.5);
    double xi = XX[i];
    if(mpi_rank == 0)
      printf("%lf   %lf\n", r, xi);    
  }
#endif
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
#ifdef USE_OMP 
#pragma omp parallel
  {
    if(omp_get_thread_num() == 0) nthreads = omp_get_num_threads();
  }
#elif defined(USE_MPI)
  nthreads = mpi_size;
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

  if(mpi_rank == 0){
    printf("\n====================================\n");
    printf("Cross pair counts using grid:\n");
    printf("====================================\n");
    printf("Using n = %i threads\n", nthreads);
  }

  // How many cells in each direction we must search in the second grid
  int delta_ncells2 = (int)(ceil(rmax / cell_size2)) + 2;
  if(mpi_rank == 0)
    printf("We will go left and right: [%i] cells. The corresponding distance: +/-[%lf] Mpc/h\n", 
        delta_ncells2,  delta_ncells2 * cell_size2);

  //==========================================================
  // Loop over all the cells in grid1
  // The loops only go up to max_ix etc. since we wan to skip 
  // looping over cells that we know are empty
  //==========================================================

  int ix0, num_processed = 0;
  clock_t start = clock(), diff;
  int istart = 0, iend = max_ix + 1;
#ifdef USE_OMP 
#pragma omp parallel for private(id) reduction(+:pairs_dist2) reduction(+:pairs_dist) schedule(dynamic)
#elif defined(USE_MPI)
  int i_per_task = (max_ix + 1) / nthreads;
  istart = i_per_task * mpi_rank;
  iend   = i_per_task * (mpi_rank + 1);
  if(mpi_rank == mpi_size - 1) iend = max_ix + 1;
#endif
  for(ix0 = istart; ix0 < iend; ix0++){
#ifdef USE_OMP 
    id = omp_get_thread_num();
#elif defined(USE_MPI)
    id = mpi_rank;
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
#ifdef USE_OMP
#pragma omp critical
    {
      printf("Processed (%4i / %4i)   (ThreadID: %3i)\n", num_processed, max_ix, id);
      num_processed++;
    }
#endif
  }

#ifdef USE_MPI
  // Gather results from all CPUs
  double *recv = malloc(sizeof(double)*nbins*nthreads);
  MPI_Allgather(XY_threads[mpi_rank], nbins, MPI_DOUBLE, recv, nbins, MPI_DOUBLE, MPI_COMM_WORLD);
  for(id = 0; id < nthreads; id++)
    for(i = 0; i < nbins; i++)
      XY_threads[id][i] = recv[id * nbins + i];
  free(recv);
#endif

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
  if(mpi_rank == 0){
    printf("We have computed dist^2 for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist2, 
        (double)ngalaxies * (double)ngalaxies2, pairs_dist2/((double)ngalaxies*(double) ngalaxies2) * 100.0);
    printf("We have computed sqrt(dist^2) for [%.0lf] pairs out of [%.0lf]. That is %lf%%\n", pairs_dist, 
        (double)ngalaxies * (double)ngalaxies2, pairs_dist/((double)ngalaxies*(double) ngalaxies2) * 100.0);
    printf("This took: [%8.3lf] sec\n", (double)(1000 * diff / CLOCKS_PER_SEC)/1000.0);
  }

#ifdef OUTPUT_PAIRS_TO_SCREEN
  // Output cross pair counts
  if(mpi_rank == 0)
    printf("r          XY(r):\n");
  for(i = 0; i < nbins; i++){
    double r = rmax / (double) nbins * (i+0.5);
    double xi = XY[i];
    if(mpi_rank == 0)
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

  if(mpi_rank == 0){
    printf("\n====================================\n");
    printf("Correlation function using LS estimator:\n");
    printf("====================================\n");
    printf("Outputfile [%s] has format [r  xi  err_xi  DD  DR  RR]\n", filename);
  }

  FILE *fp; 
  if(mpi_rank == 0)
    fp = fopen(filename, "w");
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
    if(mpi_rank == 0){
#ifdef WEIGHTS
      fprintf(fp,"%le  %le  %le  %le  %le  %le\n", r, corr, err_corr, DD->paircount[i], DR->paircount[i], RR->paircount[i]);
#else
      fprintf(fp,"%le  %le  %le  %d   %d   %d\n",  r, corr, err_corr, (int)DD->paircount[i], (int)DR->paircount[i], (int)RR->paircount[i]);
#endif
    }
  }
  if(mpi_rank == 0)
    fclose(fp);
}

