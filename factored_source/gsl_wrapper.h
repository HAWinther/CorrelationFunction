#ifndef _MYGSL_HEADER
#define _MYGSL_HEADER
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

//====================================================
// GSL spline
//====================================================
typedef struct GSL_Spline {
  gsl_spline *spline;
  gsl_interp_accel *xacc;
  double xmin, xmax;
  int allocated;
} GSL_Spline;

void Create_GSL_Spline(GSL_Spline *splinecontainer, double *x, double *y, int nx);
void Free_GSL_Spline(GSL_Spline *splinecontainer);
double Lookup_GSL_Spline(GSL_Spline *splinecontainer, double x);

#endif
