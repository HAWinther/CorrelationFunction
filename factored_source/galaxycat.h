#ifndef _GALAXYCAT_H
#define _GALAXYCAT_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"

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
// Functions defined below
//====================================================
GalaxyCatalog *read_galaxies_from_file(char *filename, int npart);
int  count_lines_in_file(char *filename);
void outputGalaxies(GalaxyCatalog *cat, char *filename);
void compute_boxsize_shift_positions(GalaxyCatalog *cat, GalaxyCatalog *cat2, double *box);
void free_cat(GalaxyCatalog *cat);

#endif
