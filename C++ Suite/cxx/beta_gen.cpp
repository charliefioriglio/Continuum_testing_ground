#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "tools.h"
#include "molecule.h"
#include "dyson.h"
#include "num_eikr.h"
#include "grid.h"
#include "angle_grid.h"
#include "beta.h"

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " input_file output_file [--beta]" << std::endl;
        // Output file arg is kept to match dyson_gen signature but maybe unused or used for beta.csv
        return 1;
    }

    std::string input_file = argv[1];
    std::string output_file = argv[2]; // Unused for Beta, but arguments must match python call structure?
    // Actually, python doesn't call this yet. I will call it manually.
    
    std::ifstream in(input_file);
    if (!in) {
        std::cerr << "Error: Could not open input file " << input_file << std::endl;
        return 1;
    }

    Molecule mol;
    int n_atoms;
    if (!(in >> n_atoms)) return 1;

    for (int i = 0; i < n_atoms; ++i) {
        std::string sym;
        int idx;
        double x, y, z;
        in >> sym >> idx >> x >> y >> z;
        mol.add_atom(sym, idx, x, y, z);
    }
    
    int n_shells;
    in >> n_shells;
    for (int i = 0; i < n_shells; ++i) {
        int atom_idx, l, n_prim;
        bool is_pure;
        in >> atom_idx >> l >> is_pure >> n_prim;
        
        std::vector<double> exps(n_prim);
        std::vector<double> coeffs(n_prim);
        for (int j = 0; j < n_prim; ++j) {
            in >> exps[j] >> coeffs[j];
        }
        mol.add_shell_to_atom(atom_idx, l, is_pure, exps, coeffs);
    }
    
    int num_dyson_orbs;
    if (!(in >> num_dyson_orbs)) num_dyson_orbs = 1; 

    std::vector<Dyson> dysons;
    for (int d = 0; d < num_dyson_orbs; ++d) {
        int n_coeffs;
        double norm_val = 1.0;
        in >> n_coeffs >> norm_val; 
        std::vector<double> coeffs(n_coeffs);
        for (int i = 0; i < n_coeffs; ++i) in >> coeffs[i];
        
        std::string label = (d==0) ? "Left" : "Right";
        Dyson d_obj(&mol, coeffs, label);
        d_obj.qchem_norm = norm_val;
        dysons.push_back(d_obj);
    }
    
    double x0, x1, y0, y1, z0, z1, step;
    in >> x0 >> x1 >> y0 >> y1 >> z0 >> z1 >> step;
    
    // Renormalize
    for(auto& do_obj : dysons) {
        do_obj.renormalize(x0, x1, y0, y1, z0, z1, step);
    }
    
    // Make Grid for integration
    UniformGrid grid(x0, x1, y0, y1, z0, z1, step);

    // Default Beta Energies
    std::vector<double> beta_energies = {
        0.0001, 0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5
    };
    
    int n_points = 150; // Default (Hardcoded)
    bool use_hardcoded = true;

    // Parse args
    for(int i=3; i<argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--energies") {
            beta_energies.clear();
            int j = i + 1;
            while(j < argc) {
                std::string val = argv[j];
                if (val.substr(0, 2) == "--") break; 
                try {
                    beta_energies.push_back(std::stod(val));
                } catch(...) {
                    break; 
                }
                j++;
            }
            i = j - 1;
        } else if (arg == "--points") {
            if (i + 1 < argc) {
                try {
                    n_points = std::stoi(argv[i+1]);
                    use_hardcoded = false;
                    i++;
                } catch (...) {}
            }
        }
    }
    
    if (beta_energies.empty()) {
        std::cerr << "Error: No energies specified." << std::endl;
        return 1;
    }

    AngleGrid angle_grid;
    if (use_hardcoded && n_points == 150) {
        std::cerr << "Generating Hardcoded Angle Grid (150 pts)..." << std::endl;
        angle_grid.GenerateHardcoded();
    } else {
        std::cerr << "Generating Repulsion Angle Grid (" << n_points << " pts)..." << std::endl;
        angle_grid.GenerateRepulsion(n_points);
    }
    
    const Dyson& L_orig = dysons[0];
    const Dyson& R_orig = (dysons.size() > 1) ? dysons[1] : dysons[0];

    // Centering Logic (User Request)
    // 1. Calculate Centroid using the bounding box of the grid
    std::cout << "Calculating Dyson centroid..." << std::endl;
    // Use the grid parameters from input for bounding box
    Dyson::Vector3 centroid = L_orig.get_centroid(grid.xmin, grid.xmax, grid.ymin, grid.ymax, grid.zmin, grid.zmax, grid.dx);
    
    std::cout << "Dyson Centroid: (" << centroid.x << ", " << centroid.y << ", " << centroid.z << ")" << std::endl;
    std::cout << "Shifting molecule to center Dyson orbital at (0,0,0)..." << std::endl;
    
    // 2. Shift Molecule
    mol.shift_geometry(-centroid.x, -centroid.y, -centroid.z);
    
    // 3. Update Dyson Basis Caches
    for(auto& d : dysons) {
        d.update_geometry();
    }
    
    // Get updated references
    const Dyson& L = dysons[0];
    const Dyson& R = (dysons.size() > 1) ? dysons[1] : dysons[0];

    // Use NumEikr for optimized calculation (ezDyson logic)
    std::cout << "Calculating Beta parameters using NumEikr (ezDyson algorithm)..." << std::endl;
    NumEikr num_eikr;
    // Note: Dyson objects L and R passed separately
    num_eikr.compute(grid, L, R, angle_grid, beta_energies);
    
    // Write Results
    std::ofstream beta_file(output_file);
    beta_file << "eKE,SigmaPar,SigmaPerp,Beta\n";
    
    for(size_t i=0; i<beta_energies.size(); ++i) {
        double eKE = beta_energies[i];
        double par = num_eikr.get_sigma_par(i);
        double perp = num_eikr.get_sigma_perp(i);
        
        // Beta formula: 2(Par - Perp) / (Par + 2*Perp)
        double denom = par + 2.0 * perp;
        double beta = 0.0;
        if (std::abs(denom) > 1e-14) {
            beta = 2.0 * (par - perp) / denom;
        }
        
        // Apply scaling for absolute cross sections (optional, beta is ratio)
        // ezDyson scale: norm * kwave * ene * dyson_norm
        // We output raw Par/Perp proportional values here, consistent with user needs for Beta.
        // If absolute XS needed, apply scale.
        // par *= scale; perp *= scale;
        
        beta_file << eKE << "," << par << "," << perp << "," << beta << "\n";
    }
    beta_file.close();
    
    std::cout << "Done. Written to " << output_file << std::endl;
    
    return 0;
}
