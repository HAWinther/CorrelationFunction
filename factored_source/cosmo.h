#ifndef _COSMO_H
#define _COSMO_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gsl_wrapper.h"
#include "global.h"

GSL_Spline *global_spline_rofz;
GSL_Spline* create_rofz_spline(double OmegaM);
int ode_rofz(double z, const double r[], double drdz[], void *params);
double r_of_z(double z);

#endif
