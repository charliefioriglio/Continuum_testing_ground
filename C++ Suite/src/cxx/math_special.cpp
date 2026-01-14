#include "math_special.h"
#include <cmath>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace MathSpecial {

    double CylBesselJ(double nu, double x) {
        if (x < 0) return 0.0; // Assume x>=0
        if (x < 1e-10) {
            return (nu < 1e-10) ? 1.0 : 0.0;
        }

        // Series expansion for small x
        // Use threshold x < 15.0 or x < sqrt(nu+1)*something?
        // Series converges well when x < nu.
        // Asymptotic works when x >> nu.
        
        // Let's use x < 20.0 for series?
        // If x is large, terms get large and cancellation error occurs.
        
        if (x < 12.0) {
             double x2 = x * 0.5;
             double term = std::pow(x2, nu) / std::tgamma(nu + 1);
             double sum = term;
             double x2_sq = - x2 * x2; // - (x/2)^2
             
             for (int k = 1; k < 100; ++k) {
                 term *= x2_sq / (k * (nu + k));
                 sum += term;
                 if (std::abs(term) < 1e-15 * std::abs(sum)) break;
             }
             return sum;
        } else {
             // Asymptotic expansion (Hanckel's)
             // J_nu(x) ~ sqrt(2/pi/x) * [ P cos(chi) - Q sin(chi) ]
             // chi = x - (nu/2 + 1/4)pi
             // P ~ 1 - ...
             // Q ~ ...
             // For simplicity, 0th order asymptotic:
             // sqrt(2/pi/x) * cos(x - nu*pi/2 - pi/4)
             
             // Let's try to add a few correction terms if possible, but 0th order is decent for x > 12 usually?
             // Actually x>12 might have only 1e-3 accuracy with 0th order.
             // We need better.
             
             // P(nu, x) = 1 - (4nu^2-1)(4nu^2-9)/(2! (8x)^2) + ...
             // Q(nu, x) = (4nu^2-1)/(1! 8x) - ...
             
             double mu = 4.0 * nu * nu;
             double term8x = 8.0 * x;
             
             double P = 1.0 - (mu - 1.0)*(mu - 9.0) / (2.0 * term8x * term8x);
             double Q = (mu - 1.0) / term8x - (mu - 1.0)*(mu - 9.0)*(mu - 25.0) / (6.0 * std::pow(term8x, 3));
             
             double chi = x - (nu * 0.5 + 0.25) * M_PI;
             double sqrt_factor = std::sqrt(2.0 / (M_PI * x));
             
             return sqrt_factor * (P * std::cos(chi) - Q * std::sin(chi));
        }
    }

    double SphericalBesselJ(int l, double x) {
        if (std::abs(x) < 1e-8) {
            return (l == 0) ? 1.0 : 0.0;
        }

        // If l < x, forward recursion is stable enough.
        // If l > x, use backward recursion (Miller's Method).
        if (l < x) {
             double j_minus_1 = std::sin(x) / x;
             double j_0 = std::sin(x) / (x * x) - std::cos(x) / x;
             
             if (l==0) return std::sin(x)/x;
             if (l==1) return j_0;

             double j_curr = j_0; // j1
             double j_prev = std::sin(x)/x; // j0
             
             for (int i = 2; i <= l; ++i) {
                 double factor = (2.0 * (i - 1) + 1.0) / x;
                 double temp = factor * j_curr - j_prev;
                 j_prev = j_curr;
                 j_curr = temp;
             }
             return j_curr;
        } else {
            // Backward Recursion
            // Start high
            int l_start = l + 20 + int(x); // Heuristic
            
            double j_next = 0.0;
            double j_curr = 1.0e-30;
            
            double computed_j_l = 0.0;
            
            for (int k = l_start; k >= 1; --k) {
                 double factor = (2.0 * k + 1.0) / x;
                 double j_prev = factor * j_curr - j_next;
                 
                 if (k == l) computed_j_l = j_curr; 
                 if ((k-1) == l) computed_j_l = j_prev; 

                 j_next = j_curr;
                 j_curr = j_prev;
            }
            // Now j_curr is j0.
            double true_j0 = std::sin(x) / x;
            double scale = true_j0 / j_curr;
            
            return computed_j_l * scale;
        }
    }

    // Adapted from reference sph.C
    std::complex<double> SphericalHarmonicY(int l, int m, double theta, double phi) {
        double x = std::cos(theta);
        int abs_m = std::abs(m);
        if (abs_m > l) return 0.0;
        
        // Compute Associated Legendre Polynomial P_l^m(x)
        double P_lm = 0.0;
        
        // 1. P_m^m
        double sin_theta = std::sin(theta);
        double cur = 1.0;
        // (2m-1)!!
        for (int i = 1; i <= abs_m; i++) cur *= (2*i - 1);
        if (abs_m % 2 != 0) cur = -cur;
        cur *= std::pow(sin_theta, abs_m);
        
        if (l == abs_m) {
            P_lm = cur;
        } else {
            // 2. P_{m+1}^m
            double P_m_m = cur;
            double P_mp1_m = x * (2 * abs_m + 1) * P_m_m;
            if (l == abs_m + 1) {
                P_lm = P_mp1_m;
            } else {
                // 3. P_l^m recursion
                double P_prev = P_mp1_m;
                double P_prev2 = P_m_m;
                double P_curr = 0.0;
                for (int ll = abs_m + 2; ll <= l; ll++) {
                    P_curr = ((2.0 * ll - 1.0) * x * P_prev - (ll + abs_m - 1.0) * P_prev2) / (ll - abs_m);
                    P_prev2 = P_prev;
                    P_prev = P_curr;
                }
                P_lm = P_curr;
            }
        }
        
        // Normalization
        double num = 1.0;
        double den = 1.0;
        for (int k = 1; k <= (l - abs_m); k++) num *= k;
        for (int k = 1; k <= (l + abs_m); k++) den *= k;
        
        double norm = std::sqrt( ((2.0*l + 1.0) / (4.0 * M_PI)) * (num / den) );
        double res = norm * P_lm;
        
        std::complex<double> phase = std::exp(std::complex<double>(0, m * phi));
        
        // Handle negative m: Y_{l,-m} = (-1)^m Y_{l,m}^*
        if (m < 0) {
            // (-1)^m for negative m is (-1)^(-|m|) = (-1)^|m|
            if (abs_m % 2 != 0) res = -res;
        }

        return res * phase;
    }

}
