#include "beta.h"
#include "continuum.h"
#include "rotation.h"
#include "tools.h"
#include <cmath>
#include <complex>
#include <iostream>

std::complex<double> BetaCalculator::ComputeNumericalMatrixElement(
    const Dyson& dyson,
    const UniformGrid& grid,
    const double* k_vec,
    const double* pol_vec
) {
    std::complex<double> integral(0.0, 0.0);
    double dV = grid.dx * grid.dy * grid.dz;
    
    double k_mag = std::sqrt(k_vec[0]*k_vec[0] + k_vec[1]*k_vec[1] + k_vec[2]*k_vec[2]);
    int l_max = 10 + int(k_mag * 15.0); 
    
    for (int ix = 0; ix < grid.nx; ++ix) {
        double x = grid.xmin + ix * grid.dx;
        for (int iy = 0; iy < grid.ny; ++iy) {
            double y = grid.ymin + iy * grid.dy;
            for (int iz = 0; iz < grid.nz; ++iz) {
                double z = grid.zmin + iz * grid.dz;
                
                double dyson_val = dyson.evaluate(x, y, z);
                if (std::abs(dyson_val) < 1e-12) continue; 
                
                double r_dot_eps = x * pol_vec[0] + y * pol_vec[1] + z * pol_vec[2];
                double r_vec[3] = {x, y, z};
                std::complex<double> psi_val = Continuum::EvaluatePlaneWaveExpansion(k_vec, r_vec, l_max);
                
                integral += std::conj(psi_val) * r_dot_eps * dyson_val;
            }
        }
    }
    
    return integral * dV;
}

std::vector<BetaResult> BetaCalculator::CalculateBeta(
    const Dyson& dyson_L,
    const Dyson& dyson_R,
    const UniformGrid& grid,
    const AngleGrid& angle_grid,
    const std::vector<double>& photoelectron_energies_ev 
) {
    std::vector<BetaResult> results;
    const double HARTREE_EV = 27.211386;
    
    double pol_lab[3] = {0.0, 0.0, 1.0};      
    double k_par_lab[3] = {0.0, 0.0, 1.0};    
    double k_perp1_lab[3] = {1.0, 0.0, 0.0};  
    double k_perp2_lab[3] = {0.0, 1.0, 0.0};  
    
    for (double E_eV : photoelectron_energies_ev) {
        double E_au = E_eV / HARTREE_EV;
        double k_mag = std::sqrt(2.0 * E_au);
        
        double sum_sigma_par = 0.0;
        double sum_sigma_perp = 0.0;
        
        for (const auto& orient : angle_grid.points) {
            RotationMatrix R;
            R.SetFromEuler(orient.alpha, orient.beta, orient.gamma);
            RotationMatrix RT = R.Transpose();
            
            double pol_mol[3] = {pol_lab[0], pol_lab[1], pol_lab[2]};
            RT.Apply(pol_mol[0], pol_mol[1], pol_mol[2]);
            
            auto rotate_k = [&](const double* k_lab_vec, double k_mag) {
                double k_mol[3] = {k_lab_vec[0]*k_mag, k_lab_vec[1]*k_mag, k_lab_vec[2]*k_mag};
                RT.Apply(k_mol[0], k_mol[1], k_mol[2]);
                return std::vector<double>{k_mol[0], k_mol[1], k_mol[2]};
            };

            auto k_par_mol = rotate_k(k_par_lab, k_mag);
            auto k_perp1_mol = rotate_k(k_perp1_lab, k_mag);
            auto k_perp2_mol = rotate_k(k_perp2_lab, k_mag);
            
            std::complex<double> M_par_L = ComputeNumericalMatrixElement(dyson_L, grid, k_par_mol.data(), pol_mol);
            std::complex<double> M_par_R = ComputeNumericalMatrixElement(dyson_R, grid, k_par_mol.data(), pol_mol);
            
            std::complex<double> A_par_L = M_par_L;
            std::complex<double> A_par_R = std::conj(M_par_R); 
            double sigma_par_orient = std::real(A_par_L * A_par_R);
            
            std::complex<double> M_perp1_L = ComputeNumericalMatrixElement(dyson_L, grid, k_perp1_mol.data(), pol_mol);
            std::complex<double> M_perp1_R = ComputeNumericalMatrixElement(dyson_R, grid, k_perp1_mol.data(), pol_mol);
            double sigma_perp1_orient = std::real(M_perp1_L * std::conj(M_perp1_R));
            
            std::complex<double> M_perp2_L = ComputeNumericalMatrixElement(dyson_L, grid, k_perp2_mol.data(), pol_mol);
            std::complex<double> M_perp2_R = ComputeNumericalMatrixElement(dyson_R, grid, k_perp2_mol.data(), pol_mol);
            double sigma_perp2_orient = std::real(M_perp2_L * std::conj(M_perp2_R));
            
            double sigma_perp_orient = 0.5 * (sigma_perp1_orient + sigma_perp2_orient);
            
            sum_sigma_par += sigma_par_orient * orient.weight;
            sum_sigma_perp += sigma_perp_orient * orient.weight;
        }
        
        double norms = dyson_L.qchem_norm * dyson_R.qchem_norm;
        double sigma_par = sum_sigma_par * norms;
        double sigma_perp = sum_sigma_perp * norms;
        
        double beta = 2.0 * (sigma_par - sigma_perp) / (sigma_par + 2.0 * sigma_perp);
        
        results.push_back({E_eV, sigma_par, sigma_perp, beta});
    }
    
    return results;
}
