#include "wavetypes.h"

//! Calc by recurrence all Spherical Rkl's at one x,y,z point.
//AK!! lmax is bad, need to change such that the loops run accordingly, i.e., < lmax, not <=lmax
void SphericalRkl(double *rklmax, double x, double y, double z, double kv, int lmax)
{
  double rv = sqrt(x*x+y*y+z*z);
  double kr = kv*rv ;
  for (int j=0; j<=lmax; j++)
    rklmax[j] = gsl_sf_bessel_jl(j,kr) ;
}

// Calc all Coulomb waves at one x,y,z point.
void CoulombRkl(double *rklmax, double x, double y, double z, double kv, int lmax, double radfn)
{
  double rv = sqrt(x*x+y*y+z*z);
  double ro = kv*rv;
  double z_ion = radfn;
  double z_el = -1.0;
  double eta = z_ion*z_el/kv;
  double lminc=0;
  double fc_array[lmax+1];
  double F_exponent;
  gsl_sf_coulomb_wave_F_array(lminc, lmax, eta, ro, fc_array, &F_exponent) ;
  for (int j=0; j<=lmax; j++) {
    if (ro>0)
      rklmax[j] = fc_array[j]/ro ;
    else 
// if ro=0, above expression gives infinity. This is necessary, but probably is an approximation.
// Need to find more elegant solution.
     rklmax[j] = 0 ;
 }
}

//! Calc all radial part functions Rkl's at one x,y,z point.
void Rkl(double radfn, double *rklmax, double x, double y, double z, double kv, int lmax)
{
  if (radfn==0)
      SphericalRkl(rklmax,x,y,z,kv,lmax);
  else
      CoulombRkl(rklmax,x,y,z,kv,lmax,radfn);
}


