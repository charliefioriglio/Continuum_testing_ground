#ifndef MATH_SPECIAL_H
#define MATH_SPECIAL_H

#include <complex>
#include <vector>

namespace MathSpecial {

    // Spherical Bessel function j_l(x)
    double SphericalBesselJ(int l, double x);

    // Spherical Harmonic Y_lm(theta, phi)
    // Returns complex value Y_lm
    std::complex<double> SphericalHarmonicY(int l, int m, double theta, double phi);

}

#endif
