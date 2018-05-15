#include "gsl_wrapper.h"

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
