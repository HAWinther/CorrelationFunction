#ifndef _BINNING_H
#define _BINNING_H
#include <stdlib.h>

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

PairCountBinning *create_binning(int nbins, double rmax);
void free_binning(PairCountBinning *pc);

#endif
