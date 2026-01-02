#include "math_special.h"
#include <cmath>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace MathSpecial {

    double SphericalBesselJ(int l, double x) {
        if (std::abs(x) < 1e-8) {
            return (l == 0) ? 1.0 : 0.0;
        }

        // If l < x, forward recursion is stable enough.
        // If l > x, use backward recursion (Miller's Method).
        if (l < x) {
             double j_minus_1 = std::sin(x) / x;
             double j_0 = std::sin(x) / (x * x) - std::cos(x) / x;
             
             if (l==0) return j_minus_1; // Wait, j_0 is usually l=0. Code above calls sin/x j_0? 
             // Logic: j0 = sin(x)/x. j1 = sin(x)/x^2 - cos(x)/x.
             if (l==0) return std::sin(x)/x;
             if (l==1) return std::sin(x)/(x*x) - std::cos(x)/x;

             double j_curr = std::sin(x)/(x*x) - std::cos(x)/x; // j1
             double j_prev = std::sin(x)/x; // j0
             
             for (int i = 2; i <= l; ++i) {
                 double j_next = (2.0 * i - 1.0) / x * j_curr - j_prev; // Using ratio relation
                 // Relation: j_{l-1} + j_{l+1} = (2l+1)/x * j_l
                 // j_{l+1} = (2l+1)/x * j_l - j_{l-1}
                 // Here i is taking place of "l" in formula?
                 // If loop starts l=2. we compute j2 from j1, j0.
                 // (2*1 + 1)/x * j1 - j0.  (Using l=1 index in formula).
                 // My loop uses `i`. Formula needs `l` (the index of `curr`).
                 // Formula: j_{i} = (2(i-1)+1)/x * j_{i-1} - j_{i-2}
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
            // Arbitrary start
            double j_next = 0.0;
            double j_curr = 1.0e-30;
            
            double scale_factor = 0.0;
            bool found = false;
            double target_val = 0.0;
            
            for (int k = l_start; k >= 0; --k) {
                // j_{k-1} = (2k+1)/x * j_k - j_{k+1}
                double factor = (2.0 * k + 1.0) / x;
                double j_prev = factor * j_curr - j_next;
                
                if (k == l) target_val = j_curr; // Capture the unscaled value for target l
                
                // Shift
                j_next = j_curr;
                j_curr = j_prev;
                
                // Normalize using j0
                if (k == 0) {
                    // j_curr is now j(-1) ?? No.
                    // Loop ends after k=0 iteration.
                    // In k=0 iteration:
                    //   Compute j_-1 from j0, j1? No.
                    //   We want to compute down to j0.
                    // Loop condition k>=1?
                }
            }
            // Let's refine loop:
            // We know j_{l+1}, j_{l+2}... 
            // Relation: j_{k} = (2k+3)/x * j_{k+1} - j_{k+2}?? 
            // No. j_{k-1} + j_{k+1} = (2k+1)/x * j_k
            // j_{k-1} = (2k+1)/x * j_k - j_{k+1}
            // If we are at k, we compute k-1.
            
            j_next = 0.0;      // j_{start+1} effectively
            j_curr = 1.0e-30;  // j_{start}
            
            double computed_j_l = 0.0;
            
            for (int k = l_start; k >= 1; --k) {
                 double factor = (2.0 * k + 1.0) / x;
                 double j_prev = factor * j_curr - j_next;
                 
                 if (k == l) computed_j_l = j_curr; // Store unscaled j_l
                 if ((k-1) == l) computed_j_l = j_prev; // Store if l was k-1

                 j_next = j_curr;
                 j_curr = j_prev;
            }
            // Now j_curr is j0.
            // True j0 = sin(x)/x
            double true_j0 = std::sin(x) / x;
            double scale = true_j0 / j_curr;
            
            return computed_j_l * scale;
        }
    }

    // Adapted from reference sph.C
    // Hardcoded coefficients for up to l=4 to save space, generic for higher if needed?
    // The reference `sph.C` has up to l=10. I will implement up to l=5 completely.
    // Normalized spherical harmonics.
    
    std::complex<double> SphericalHarmonicY(int l, int m, double theta, double phi) {
        // Condon-Shortley phase convention is usually included in std definition.
        // Y_lm = N * P_l^m(cos theta) * e^{im phi}
        // m can be negative.
        // P_l^{-m} = (-1)^m (l-m)!/(l+m)! P_l^m
        
        // We will calculate P_l^|m|(x) where x = cos(theta)
        // using standard recursion.
        
        double x = std::cos(theta);
        int abs_m = std::abs(m);
        if (abs_m > l) return 0.0;
        
        // Compute Associated Legendre Polynomial P_l^m(x)
        // Re-implement simplified recursion for P_l^m
        double P_lm = 0.0;
        
        // 1. P_m^m
        double sin_theta = std::sin(theta);
        double cur = 1.0;
        // (2m-1)!! = 1 * 3 * ... * (2m-1)
        for (int i = 1; i <= abs_m; i++) cur *= (2*i - 1);
        // * (-1)^m * sin^m
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
        // N = sqrt( (2l+1)/4pi * (l-m)!/(l+m)! )
        double num = 1.0;
        double den = 1.0;
        
        for (int k = 1; k <= (l - abs_m); k++) num *= k;
        for (int k = 1; k <= (l + abs_m); k++) den *= k;
        
        double norm = std::sqrt( ((2.0*l + 1.0) / (4.0 * M_PI)) * (num / den) );
        
        double res = norm * P_lm;
        
        // Phase e^{im phi}
        std::complex<double> phase = std::exp(std::complex<double>(0, m * phi));
        
        // Handle negative m: Y_{l,-m} = (-1)^m Y_{l,m}^*
        if (m < 0) {
            if (abs_m % 2 != 0) res = -res;
        }

        return res * phase;
    }

}
