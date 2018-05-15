#ifndef _CORRELATORS_H
#define _CORRELATORS_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "global.h"
#include "galaxycat.h"
#include "grid.h"
#include "binning.h"

void brute_force_pair_counting(GalaxyCatalog *cat, PairCountBinning *pc);
void brute_force_cross_pair_counting(GalaxyCatalog *cat, GalaxyCatalog *cat2, PairCountBinning *pc);
void grid_pair_counting(Grid *grid, PairCountBinning *pc);
void grid_cross_pair_counting(Grid *grid, Grid *grid2, PairCountBinning *pc);
void compute_correlation_function(PairCountBinning *DD, PairCountBinning *DR, PairCountBinning *RR, char *filename);

#endif
