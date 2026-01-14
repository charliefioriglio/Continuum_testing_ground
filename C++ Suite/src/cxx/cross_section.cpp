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
    // Global accumulator (merged results from threads)
    std::vector<std::vector<std::vector<PartialWaveAccumulator>>> final_partial_sums(
        num_energies, 
        std::vector<std::vector<PartialWaveAccumulator>>(
            num_lm,
            std::vector<PartialWaveAccumulator>(num_alpha, {0.0, 0.0})
        )
    );
    
    // Parallel Region
    #pragma omp parallel
    {
        // Thread-local Accumulator
        auto thread_partial_sums = final_partial_sums; // Copy structure (zeros)
        // Reset just in case copy not zero
        for(auto& e_vec : thread_partial_sums)
            for(auto& lm_vec : e_vec)
                for(auto& a_val : lm_vec) a_val = {0.0, 0.0};
        
        #pragma omp for
        for (int ix = 0; ix < nx; ++ix) {
            double x = x0 + ix * step;
            for (int iy = 0; iy < ny; ++iy) {
                 double y = y0 + iy * step;
                 for (int iz = 0; iz < nz; ++iz) {
                    double z = z0 + iz * step;
                    
                    int idx = ix * (ny * nz) + iy * nz + iz;
                    double phi_L = phi_L_vals[idx];
                    double phi_R = phi_R_vals[idx];
                    
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
                        for(int m=-l; m<=l; ++m) {
                             std::complex<double> Ylm_conj = std::conj(MathSpecial::SphericalHarmonicY(l, m, theta, phi_ang));
                             
                             std::complex<double> termL_base = phi_L * Ylm_conj;
                             std::complex<double> termR_base = phi_R * Ylm_conj;
                             
                             for(int e=0; e<num_energies; ++e) {
                                 if(k_values[e] <= 0) continue;
                                 
                                 double k = k_values[e];
                                 double jl = MathSpecial::SphericalBesselJ(l, k*r);
                                 double factor = jl * dV;
                                 
                                 // Alpha 0 (x)
                                 thread_partial_sums[e][lm_idx][0].val_L += termL_base * dip_x * factor;
                                 thread_partial_sums[e][lm_idx][0].val_R += termR_base * dip_x * factor;
                                 
                                 // Alpha 1 (y)
                                 thread_partial_sums[e][lm_idx][1].val_L += termL_base * dip_y * factor;
                                 thread_partial_sums[e][lm_idx][1].val_R += termR_base * dip_y * factor;
                                 
                                 // Alpha 2 (z)
                                 thread_partial_sums[e][lm_idx][2].val_L += termL_base * dip_z * factor;
                                 thread_partial_sums[e][lm_idx][2].val_R += termR_base * dip_z * factor;
                             }
                             lm_idx++;
                        }
                    }
                 }
            }
        } // End Grid Loop (Parallel)
        
        // Merge thread results
        #pragma omp critical
        {
            for(int e=0; e<num_energies; ++e) {
                for(int lm=0; lm<num_lm; ++lm) {
                    for(int a=0; a<num_alpha; ++a) {
                        final_partial_sums[e][lm][a].val_L += thread_partial_sums[e][lm][a].val_L;
                        final_partial_sums[e][lm][a].val_R += thread_partial_sums[e][lm][a].val_R;
                    }
                }
            }
        }
    } // End Parallel
    
    // Final Assembly (Use merged results)
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
                std::complex<double> amp_L = final_partial_sums[e][lm][alpha].val_L;
                std::complex<double> amp_R = final_partial_sums[e][lm][alpha].val_R;
                
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


#include "point_dipole.h"

std::vector<double> CrossSectionCalculator::ComputePointDipoleCrossSection(
    const Dyson& dyson_L,
    const Dyson& dyson_R,
    const UniformGrid& grid,
    const std::vector<double>& photon_energies_ev,
    double ionization_energy_ev,
    int l_max,
    double dipole_magnitude
) {
    std::vector<double> results;
    results.reserve(photon_energies_ev.size());
    
    // Fill Grid Caches
    std::vector<double> phi_L_vals;
    std::vector<double> phi_R_vals;
    phi_L_vals.reserve(grid.nx * grid.ny * grid.nz);
    phi_R_vals.reserve(grid.nx * grid.ny * grid.nz);
    
    double x0 = grid.xmin; double y0 = grid.ymin; double z0 = grid.zmin;
    double step = grid.dx;
    
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
    
    std::vector<double> k_values;
    for(double E_ev : photon_energies_ev) {
        double E_au = E_ev * EV_TO_HARTREE;
        double IE_au = ionization_energy_ev * EV_TO_HARTREE;
        double eKE = E_au - IE_au;
        if(eKE > 0) k_values.push_back(std::sqrt(2.0 * eKE));
        else k_values.push_back(0.0);
    }
    
    PointDipole pd(dipole_magnitude);
    
    int nx = grid.nx; int ny = grid.ny; int nz = grid.nz;
    double dV = step * step * step;

    // Loop Energies
    for(size_t e=0; e < photon_energies_ev.size(); ++e) {
        if (k_values[e] <= 0) {
            results.push_back(0.0);
            continue;
        }
        
        double k = k_values[e];
        double total_sum = 0.0;
        
        // Sum over m and N
        for (int lam = -l_max; lam <= l_max; ++lam) {
            Eigensystem sys = pd.GetEigensystem(lam, l_max);
            
            int n_modes = sys.l_vals.size();
            for (int N = 0; N < n_modes; ++N) {
                double eigval = sys.eigvals[N];
                if (eigval < -0.25) continue;
                
                double L_eff = 0.5 * (-1.0 + std::sqrt(1.0 + 4.0 * eigval));
                
                // Compute Moment vector I_alpha for this mode (L and R)
                // I_alpha = Sum_l c_{lN}^* Integral( phi * r_alpha * Ylm* * j_{Leff}(kr) )
                
                std::complex<double> I_L[3] = {0.0, 0.0, 0.0};
                std::complex<double> I_R[3] = {0.0, 0.0, 0.0};
                
                // Pre-calculate integrals for each l component: J_{l,alpha}
                // Actually easier to Loop Grid -> Accumulate contributions to I_alpha directly
                // accumulating separately for each l component is standard, 
                // but here c_{lN} mixes them.
                // Optim: Loop Grid -> Compute Ylm* -> Add to accumulators for each l.
                
                std::vector<std::vector<std::complex<double>>> partial_integrals_L(n_modes, std::vector<std::complex<double>>(3, 0.0));
                std::vector<std::vector<std::complex<double>>> partial_integrals_R(n_modes, std::vector<std::complex<double>>(3, 0.0));
                
                // BUT we are inside Energy Loop -> Grid loop is expensive.
                // We should invert: Loop Grid -> Loop Energies?
                // But eigenmodes/Leff depend on nothing, but Bessel depends on k.
                // We can stick to Energy -> Grid if performance allows. The user asked for 7 D values * 2 Orbitals = 14 runs.
                // 14 runs is fine.
                // Speed up: Precompute Ylm on grid? Done in PWE implicitly?
                // Let's just do OpenMP the grid loop here.
                
                // We need to compute I_L/R for this SPECIFIC mode N.
                // Or compute for ALL modes N at once to save grid pass?
                // Yes, do all modes at once.
            }
        }
        
        // Correct approach: Energy -> Grid
        // Accumulate moments for ALL m, N. 
        // Need storage: moments[lam_idx][N][alpha]
        // lam goes -l_max to l_max.
        // N goes 0 to n_modes(lam).
        
        struct Moment { std::complex<double> val_L; std::complex<double> val_R; };
        // Flattened storage? Map is slow.
        // Let's use `vector<vector<Moment[3]>>` where outer is m, inner is N.
        
        std::vector<std::vector<std::vector<Moment>>> mode_moments; 
        // mode_moments[lam_offset][N][alpha]
        mode_moments.resize(2*l_max + 1);
        
        // Setup systems and resize
        for (int lam = -l_max; lam <= l_max; ++lam) {
            Eigensystem sys = pd.GetEigensystem(lam, l_max);
            mode_moments[lam + l_max].resize(sys.l_vals.size(), std::vector<Moment>(3, {0.0, 0.0}));
        }
        
        #pragma omp parallel
        {
            auto local_moments = mode_moments;
            // Zero out
            for(auto& v1 : local_moments) for(auto& v2 : v1) for(auto& m : v2) m = {0.0, 0.0};
            
            #pragma omp for
            for (int ix = 0; ix < nx; ++ix) {
                for (int iy = 0; iy < ny; ++iy) {
                    for (int iz = 0; iz < nz; ++iz) {
                        int idx = ix * ny * nz + iy * nz + iz;
                        double phi_L = phi_L_vals[idx];
                        double phi_R = phi_R_vals[idx];
                        
                        if (std::abs(phi_L) < 1e-15 && std::abs(phi_R) < 1e-15) continue;
                        
                        double x = x0 + ix * step;
                        double y = y0 + iy * step;
                        double z = z0 + iz * step;
                        double r_sq = x*x + y*y + z*z;
                        if (r_sq < 1e-18) continue;
                        double r = std::sqrt(r_sq);
                        double theta = std::acos(z/r);
                        double phi_ang = std::atan2(y, x);
                        
                        double dip[3] = {x, y, z};
                        
                        // Per m
                        for (int lam = -l_max; lam <= l_max; ++lam) {
                            Eigensystem sys = pd.GetEigensystem(lam, l_max);
                            int n_modes = sys.l_vals.size();
                            
                            // We need Sum_l c_{lN}* Y_{l,lam}*
                            // Optimization: Calculate Y_{l,lam}* for all l first.
                            std::vector<std::complex<double>> Y_vals(n_modes);
                            for(int i=0; i<n_modes; ++i) {
                                Y_vals[i] = std::conj(MathSpecial::SphericalHarmonicY(sys.l_vals[i], lam, theta, phi_ang));
                            }
                            
                            for (int N = 0; N < n_modes; ++N) {
                                double eigval = sys.eigvals[N];
                                if (eigval < -0.25) continue;
                                double L_eff = 0.5 * (-1.0 + std::sqrt(1.0 + 4.0 * eigval));
                                
                                // Radial part
                                double bessel = MathSpecial::CylBesselJ(L_eff + 0.5, k*r);
                                double radial = std::sqrt(M_PI/(2.0*k*r)) * bessel; 
                                if (k*r < 1e-10) { // Limit case
                                   radial = (L_eff < 0.1) ? 1.0 : 0.0;
                                }
                                
                                // Angular part: Omega* = Sum_l c_{lN}* Y* = (Sum c_{lN} Y)*
                                // Wait, omega = Sum c Y. omega* = Sum c* Y*. (Coeffs are real).
                                std::complex<double> omega_conj = 0.0;
                                for(int i=0; i<n_modes; ++i) {
                                    omega_conj += sys.eigvecs[i][N] * Y_vals[i];
                                }
                                
                                std::complex<double> basis_val = omega_conj * radial * dV;
                                
                                // Accumulate
                                for(int alpha=0; alpha<3; ++alpha) {
                                    local_moments[lam+l_max][N][alpha].val_L += phi_L * dip[alpha] * basis_val;
                                    local_moments[lam+l_max][N][alpha].val_R += phi_R * dip[alpha] * basis_val;
                                }
                            }
                        }
                    } 
                }
            } // End Grid
            
            #pragma omp critical
            {
                for(int m=0; m < (int)local_moments.size(); ++m) {
                    for(int N=0; N < (int)local_moments[m].size(); ++N) {
                         for(int a=0; a<3; ++a) {
                             mode_moments[m][N][a].val_L += local_moments[m][N][a].val_L;
                             mode_moments[m][N][a].val_R += local_moments[m][N][a].val_R;
                         }
                    }
                }
            }
        } // End Parallel
        
        // Sum Contribution
        double prefactor = (8.0 * M_PI * k * (photon_energies_ev[e] * EV_TO_HARTREE)) / C_SPEED_AU;
        
        for(int lam = -l_max; lam <= l_max; ++lam) {
             for(int N = 0; N < (int)mode_moments[lam+l_max].size(); ++N) {
                 double sum_alpha = 0.0;
                 for(int a=0; a<3; ++a) {
                     std::complex<double> mL = mode_moments[lam+l_max][N][a].val_L;
                     std::complex<double> mR = mode_moments[lam+l_max][N][a].val_R;
                     sum_alpha += std::real(std::conj(mL) * mR);
                 }
                 total_sum += sum_alpha;
             }
        }
        
        results.push_back(total_sum * prefactor / 3.0 * dyson_L.qchem_norm * dyson_R.qchem_norm * 2.0);
    }
    
    return results;
}
