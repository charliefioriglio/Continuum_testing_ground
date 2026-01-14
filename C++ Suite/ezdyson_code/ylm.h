#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

//! Spherical harmonics, Ylm part only, for the PAD calculation
double ThetaYlm(double theta, int lv, int mv);
//! Legendre polynomilas
double Pl(double theta, int lv);
