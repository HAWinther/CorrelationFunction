#ifndef _GRID_H
#define _GRID_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "galaxycat.h"

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

Grid *create_grid(int ngalaxies, double rmax, double box);
void add_galaxies_to_cells(Grid *grid, GalaxyCatalog *cat);
void free_grid(Grid *grid);

#endif
