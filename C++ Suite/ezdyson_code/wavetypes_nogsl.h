#include "tools.h"
#include "complexno.h"
#include <math.h>

//Colection of functions for calculating the radial part for different types of waves

//!Calc Spherical Rk0 at x,y,z.
//double Rk0(double kv, double kr);
//!Calc Spherical Rk1 at x,y,z.
//double Rk1(double kv, double kr);

//!Calc Spherical Rk0 at x,y,z, 2k*jl(kr)
inline double Rk0(double kr)
{
  double rk0 = sin(kr)/kr;
  return rk0;
}

//!Calc Spherical Rk1 at x,y,z.
inline double Rk1(double kr)
{
  double rk1 = (sin(kr)/pow(kr,2)-cos(kr)/kr);
  return rk1;
}

//!Calc Spherical Rkl for kr<<l.
double LimRkl(double kr, int l);
//!Calc l max for which to use Rkl recurrence
//int CalcLimL(double kr);
//! Calc by recurrence all Spherical Rkl's at one x,y,z point.
void SphericalRkl(double *rklmax, double x, double y, double z, double kv, int lmax);

/*! Calc all Coulomb waves at one x,y,z point.
Adapted from web PL/1 code: http://www.plasmaphysics.org.uk/programs/coulombc.htm  */
void CoulombRkl(double *rklmax, double x, double y, double z, double kv, int lmax, double radfn);
//! Procedure for calc Coulomb waves (see PL/1 code)
double DQG16(double *QG, double ro, double eta, double XL, double XU, int l1);
//! Procedure for calc of Coulomb waves (see PL/1 code)
double FCT(double x, double ro, double eta, int l1);
//! Calc Gamma function for Coulomb wave prefactor
Complex Gamma(double *QG, double KVG, double XDG, const Complex& LE);
//! Procedure for calc of complex Gamma function (see PL/1 code)
Complex DQG16G(double *QG, const Complex& Z, double XL, double XU);
//! Procedure for calc of complex Gamma funcion (see PL/1 code)
inline Complex FCTG(double x, const Complex& Z)
{
  return Complex(CPower(x,Z-1.0))*exp(-x);
}


//! Calc all radial part functions Rkl's at one x,y,z point.
void Rkl(double radfn, double *rklmax, double x, double y, double z, double kv, int lmax);



