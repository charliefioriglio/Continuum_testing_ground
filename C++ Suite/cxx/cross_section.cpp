#include "cross_section.h"
#include "tools.h"
#include "math_special.h"
#include "continuum.h"
#include "tools.h" // For Constants if needed
#include <cmath>
#include <complex>
#include <iostream>
// #include <omp.h>

// Constants
constexpr double HARTREE_TO_EV = 27.211386;
constexpr double EV_TO_HARTREE = 1.0 / HARTREE_TO_EV;
constexpr double C_SPEED_AU = 137.035999; 

// Plane wave expansion:
// psi_{klm} = Y_{lm}(r_hat) * j_l(kr)
// C_{klm} = i^l * integral( phi_dyson(r) * r_alpha * psi_{klm}(r) )
//
// In calculate.py, "planewave_expansion" returns Y * R.
// The integral calculation sums: conj(psi) * dipole * DO.
// Wait, eq 4: C_{klm} = i^l * int( phi * r_alpha * psi ).
// calculate.py line 124: integrand = conj(psi) * dipole * DO.
// This implies psi in Eq 4 is actually the conjugate already? 
// No, standard overlap <psi| dipole | phi>.
// If psi is the final state (continuum), then it's <psi| r | phi> = integral( psi* r phi ).
// So calculate.py is correct: conj(psi).
// Note: integral is over volume.
// Prefactor i^l is NOT in calculate.py's integral loop, but |C|^2 makes phase irrelevant?
// "Instead of |C|^2, we used the more accurate (C^L)* C^R..." = |C|^2 for same orbital.
// Since i^l is a phase, |i^l|^2 = 1. So it drops out for Total Cross Section.
// It assumes real Dyson orbital?
// The C++ `cklm.C` uses `phasefactor`.
// I will implement <psi| r | phi> and square the magnitude.

std::vector<double> CrossSectionCalculator::ComputeTotalCrossSection(
    const Dyson& dyson_L,
    const Dyson& dyson_R,
    const UniformGrid& grid,
    const std::vector<double>& photon_energies_ev,
    double ionization_energy_ev,
    int l_max
) {
    std::vector<double> results;
    results.reserve(photon_energies_ev.size());
    
    // Pre-compute Dyson values on grid to avoid re-evaluation in inner loops
    std::vector<double> phi_L_vals;
    std::vector<double> phi_R_vals;
    phi_L_vals.reserve(grid.nx * grid.ny * grid.nz);
    phi_R_vals.reserve(grid.nx * grid.ny * grid.nz);
    
    // Fill grids
    double x0 = grid.xmin; double y0 = grid.ymin; double z0 = grid.zmin;
    double step = grid.dx; // Assume uniform
    
    // Loop once to fill
    for (int ix = 0; ix < grid.nx; ++ix) {
        for (int iy = 0; iy < grid.ny; ++iy) {
            for (int iz = 0; iz < grid.nz; ++iz) {
                 double x = x0 + ix * step;
                 double y = y0 + iy * step;
                 double z = z0 + iz * step;
                 phi_L_vals.push_back(dyson_L.evaluate(x, y, z));
                 phi_R_vals.push_back(dyson_R.evaluate(x, y, z));
            }
        }
    }

    // Optimization: Loop Grid (Outer) -> Energy (Inner)
    // Avoids re-evaluating Dyson and Ylm for every energy.
    // Dyson(r) and Ylm(theta,phi) are Energy-independent.
    // Bessel(kr) depends on k(E).
    
    int num_energies = photon_energies_ev.size();
    
    // Store partial sums for each energy: total_sum[E]
    std::vector<double> energy_sums(num_energies, 0.0);
    
    // Pre-calculate k values
    std::vector<double> k_values;
    k_values.reserve(num_energies);
    for(double E_ev : photon_energies_ev) {
        double E_au = E_ev * EV_TO_HARTREE;
        double IE_au = ionization_energy_ev * EV_TO_HARTREE;
        double eKE = E_au - IE_au;
        if(eKE > 0) {
            k_values.push_back(std::sqrt(2.0 * eKE));
        } else {
            k_values.push_back(0.0);
        }
    }
    
    // Iterate Grid points once
    int nx = grid.nx; int ny = grid.ny; int nz = grid.nz;
    double dV = step * step * step;
    
    // To accumulate "sum_alpha" correctly over grid, we need separate accumulators per (Energy, alpha) ?
    // The integral is I = Sum( phi * r_alpha * conj(Y) * jl ).
    // Cross Section ~ |I|^2.
    // So we must compute the full Integral I for each Energy FIRST.
    // We cannot just sum |term|^2 over grid.
    // So we need accumulators for the Integers (Complex) for each E, l, m, alpha.
    
    // Storage: vector of Complex Accumulators.
    // Size: Energies * (L_max+1) * (2L+1) * 3 (alphas) * 2 (L/R) ??
    // That's manageable. 
    // l_max ~ 3. (L+1)^2 ~ 16. Alpha=3. Total ~ 48.
    // Energies ~ 20. Total ~ 1000 accumulators per L/R.
    
    struct PartialWaveAccumulator {
        std::complex<double> val_L;
        std::complex<double> val_R;
    };
    
    // Indexing: [energy_idx][l][m_idx][alpha]
    // Flattened or nested vectors.
    // Let's use a flat vector for efficiency.
    // Map (l, m) -> lm_index = l^2 + (m+l) ?
    // Standard ordering: l=0(m=0), l=1(m=-1,0,1)...
    // Total LM count = (l_max + 1)^2.
    
    int num_lm = (l_max + 1) * (l_max + 1);
    int num_alpha = 3;
    
    // [energy][lm][alpha]
    std::vector<std::vector<std::vector<PartialWaveAccumulator>>> partial_sums(
        num_energies, 
        std::vector<std::vector<PartialWaveAccumulator>>(
            num_lm,
            std::vector<PartialWaveAccumulator>(num_alpha, {0.0, 0.0})
        )
    );
    
    // Loop Grid
    int idx = 0;
    for (int ix = 0; ix < nx; ++ix) {
        double x = x0 + ix * step;
        for (int iy = 0; iy < ny; ++iy) {
             double y = y0 + iy * step;
             for (int iz = 0; iz < nz; ++iz) {
                double z = z0 + iz * step;
                
                double phi_L = phi_L_vals[idx];
                double phi_R = phi_R_vals[idx];
                idx++;
                
                if (std::abs(phi_L) < 1e-15 && std::abs(phi_R) < 1e-15) continue;
                
                double r_sq = x*x + y*y + z*z;
                if (r_sq < 1e-18) continue;
                double r = std::sqrt(r_sq);
                
                double theta = std::acos(z/r);
                double phi_ang = std::atan2(y, x);
                
                double dip_x = x;
                double dip_y = y;
                double dip_z = z;
                
                // For each LM
                int lm_idx = 0;
                for(int l=0; l<=l_max; ++l) {
                    
                    // Precompute Ylm for this grid point (reused across energies)
                    // But we need Ylm loop inside here.
                    
                    for(int m=-l; m<=l; ++m) {
                         std::complex<double> Ylm_conj = std::conj(MathSpecial::SphericalHarmonicY(l, m, theta, phi_ang));
                         
                         // Precompute term common to all energies: phi * conj(Y) * r_alpha
                         // For each alpha
                         std::complex<double> termL_base = phi_L * Ylm_conj;
                         std::complex<double> termR_base = phi_R * Ylm_conj;
                         
                         // Optimization: Inner Loop over Energies
                         for(int e=0; e<num_energies; ++e) {
                             if(k_values[e] <= 0) continue;
                             
                             double k = k_values[e];
                             double jl = MathSpecial::SphericalBesselJ(l, k*r); // Most expensive part potentially?
                             // Optimization: jl depends only on (l, k*r).
                             
                             double factor = jl * dV; // Include dV here
                             
                             // Integrate for each alpha
                             // Alpha 0 (x)
                             partial_sums[e][lm_idx][0].val_L += termL_base * dip_x * factor;
                             partial_sums[e][lm_idx][0].val_R += termR_base * dip_x * factor;
                             
                             // Alpha 1 (y)
                             partial_sums[e][lm_idx][1].val_L += termL_base * dip_y * factor;
                             partial_sums[e][lm_idx][1].val_R += termR_base * dip_y * factor;
                             
                             // Alpha 2 (z)
                             partial_sums[e][lm_idx][2].val_L += termL_base * dip_z * factor;
                             partial_sums[e][lm_idx][2].val_R += termR_base * dip_z * factor;
                         }
                         
                         lm_idx++;
                    }
                }
             }
        }
    }
    
    // Final Assembly
    for(int e=0; e<num_energies; ++e) {
        if(k_values[e] <= 0) {
            results.push_back(0.0);
            continue;
        }
        
        double E_ph_au = photon_energies_ev[e] * EV_TO_HARTREE;
        double k = k_values[e];
        double prefactor = (8.0 * M_PI * k * E_ph_au) / C_SPEED_AU;
        
        double total_sum = 0.0;
        
        // Sum over LM and Alpha
        for(int lm=0; lm<num_lm; ++lm) {
            double sum_alpha = 0.0;
            for(int alpha=0; alpha<3; ++alpha) {
                std::complex<double> amp_L = partial_sums[e][lm][alpha].val_L;
                std::complex<double> amp_R = partial_sums[e][lm][alpha].val_R;
                
                sum_alpha += (std::conj(amp_L) * amp_R).real();
            }
            total_sum += sum_alpha;
        }
        
        total_sum /= 3.0; // Average polarization
        double norms = dyson_L.qchem_norm * dyson_R.qchem_norm;
        double sigma = prefactor * total_sum * norms * 2.0;
        
        results.push_back(sigma);
    }

    return results;
}


