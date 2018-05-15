#include "binning.h"

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

