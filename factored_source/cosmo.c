#include "cosmo.h"

//====================================================
// Spline of r(z) = Int_0^z dz/H(z)
//====================================================
double r_of_z(double z){
  return Lookup_GSL_Spline(global_spline_rofz, z);
}

//====================================================
// The ODE dr/dz = 1/H(z) in units of Mpc/h for LCDM
// Hubble_Length = c/H0
//====================================================
int ode_rofz(double z, const double r[], double drdz[], void *params){
  const double Hubble_Length_in_Mpch = SPEED_OF_LIGHT_IN_KM_PER_SEC / 100.0;
  double OmegaM = *(double *) params;
  
  drdz[0] =  Hubble_Length_in_Mpch / sqrt(OmegaM*pow3(1.0+z) + 1.0 - OmegaM);
  return GSL_SUCCESS;
}

//====================================================
// Integrate the ODE for r(z) and return a spline of the result
//====================================================
GSL_Spline* create_rofz_spline(double OmegaM){
  const int npts      = 1000;
  const double zmax   = 5.0;
  const double deltaz = zmax/(double) (npts-1);

  // Set up ODE system
  gsl_odeiv2_system sys_rofz = {ode_rofz, NULL, 1, &OmegaM};
  gsl_odeiv2_driver *ode = gsl_odeiv2_driver_alloc_y_new (&sys_rofz, gsl_odeiv2_step_rk2, 1e-10, 1e-10, 0.0);

  double *z_arr = malloc(sizeof(double) * npts);
  double *rofz_arr = malloc(sizeof(double) * npts);

  // Initial conditions r(z=0) = 0
  double ode_z = z_arr[0] = 0.0;
  double rofz_now[1] = {0.0};
  rofz_arr[0] = 0.0;

  // Integration over redshift
  int i;
  for(i = 1; i < npts; i++){
    double znow = i * deltaz;

    // Integrate one step
    int status = gsl_odeiv2_driver_apply(ode, &ode_z, znow, rofz_now);
    if(status != GSL_SUCCESS){
      printf("Error in integrating z = %f  r = %f\n", znow, rofz_now[0]);
      exit(1);
    }

    // Store values
    z_arr[i] = znow;
    rofz_arr[i] = rofz_now[0];
  }

  // Spline up the results
  GSL_Spline *rofz_spline = malloc(sizeof(GSL_Spline));
  Create_GSL_Spline(rofz_spline, z_arr, rofz_arr, npts);

  // Free memory
  free(z_arr);
  free(rofz_arr);

  return rofz_spline;
}

