#include "tools.h"
#include "complexno.h"
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_sf_coulomb.h"
#include <math.h>

//Colection of functions for calculating the radial part for different types of waves

//Calc Spherical Rkl at x,y,z.
void SphericalRkl(double *rklmax, double x, double y, double z, double kv, int lmax);


// Calc all Coulomb waves at one x,y,z point.
void CoulombRkl(double *rklmax, double x, double y, double z, double kv, int lmax);

//! Calc all radial part functions Rkl's at one x,y,z point.
void Rkl(double radfn, double *rklmax, double x, double y, double z, double kv, int lmax);



