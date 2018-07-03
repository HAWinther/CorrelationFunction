#ifndef _GALAXY_HEADER
#define _GALAXY_HEADER

//====================================================
// A single galaxy
//====================================================
typedef struct Galaxy{
  double x[3];
#ifdef VELOCITY
  double v[3];
#endif
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

#endif
