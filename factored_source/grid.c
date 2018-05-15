#include "grid.h"

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

  if(mpi_rank == 0){
    printf("\n====================================\n"); 
    printf("Creating new grid\n"); 
    printf("====================================\n"); 
    printf("ngrid = %i [min(%i,%i)] CellSize: %lf Mpc/h\n", ngrid, ngrid1, ngrid2, box/(double) ngrid);
    printf("Total number of cells: [%i] Galaxies that will be added to grid: [%i]\n", NcellTotal, ngalaxies);
  }

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

  if(mpi_rank == 0){
    printf("\n====================================\n");
    printf("Adding galaxies to grid\n");
    printf("====================================\n");
    printf("There are [%i] empty cells (%lf %%) in the grid\n", nempty, (double)nempty / (double) NcellTotal * 100.0);
  }

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
