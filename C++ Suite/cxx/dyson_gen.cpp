#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "tools.h"
#include "molecule.h"
#include "dyson.h"
#include "grid.h"
#include "cross_section.h"
#include "cross_section.h"

// Input format expected from Python script:
// ATOM_COUNT
// SYMBOL INDEX X Y Z
// ...
// SHELL_COUNT
// ATOM_IDX L IS_PURE N_PRIM
// EXP COEF
// ...
// DYSON_COEFF_COUNT
// C_0 C_1 ...
// GRID_PARAMS
// X0 X1 Y0 Y1 Z0 Z1 STEP

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " input_file output_file" << std::endl;
        return 1;
    }

    std::string input_file = argv[1];
    std::string output_file = argv[2];
    
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
    if (!(in >> num_dyson_orbs)) num_dyson_orbs = 1; // Backward compatibility fallback (unlikely to hit if I control input)
    // Actually, old format had "n_dyson" (count of coeffs) where I now put "num_dyson_orbs".
    // To safe-guard, the old format had large int (e.g. 500) for n_coeffs. "num_dyson_orbs" will be 1 or 2.
    // If input is "500", we might fail. 
    // I will assume I always update python writer.
    
    std::vector<Dyson> dysons;
    
    // I need to handle the case where the file might just have the n_coeffs directly (old format).
    // But since I'm rewriting python, let's just assume new format 100%.
    
    for (int d = 0; d < num_dyson_orbs; ++d) {
        int n_coeffs;
        double norm_val = 1.0;
        in >> n_coeffs >> norm_val; // Expecting "N NORM"
        // If old input logic, fallback? 
        // If just N, stream parsing might fail or read next token. 
        // Python writer guarantees we send norm now.
        
        std::vector<double> coeffs(n_coeffs);
        for (int i = 0; i < n_coeffs; ++i) in >> coeffs[i];
        
        // Label
        std::string label = (d==0) ? "Left/Primary" : "Right/Secondary";
        Dyson d_obj(&mol, coeffs, label);
        d_obj.qchem_norm = norm_val; // Store new property
        dysons.push_back(d_obj);
    }
    
    double x0, x1, y0, y1, z0, z1, step;
    in >> x0 >> x1 >> y0 >> y1 >> z0 >> z1 >> step;
    
    // Centering Logic (Identical to beta_gen.cpp)
    // 1. Calculate Centroid using the bounding box of the grid
    std::cout << "Calculating Dyson centroid..." << std::endl;
    const Dyson& L_orig = dysons[0];
    Dyson::Vector3 centroid = L_orig.get_centroid(x0, x1, y0, y1, z0, z1, step);
    
    std::cout << "Dyson Centroid: (" << centroid.x << ", " << centroid.y << ", " << centroid.z << ")" << std::endl;
    std::cout << "Shifting molecule and grid to center Dyson orbital at (0,0,0)..." << std::endl;
    
    // 2. Shift Molecule
    mol.shift_geometry(-centroid.x, -centroid.y, -centroid.z);
    
    // 3. Update Dyson Basis Caches
    for(auto& d : dysons) {
        d.update_geometry();
    }
    
    // 4. Shift Grid Limits to match new centered frame
    // If the orbital moves by -C, the grid window must also move by -C to stay centered on it
    x0 -= centroid.x; x1 -= centroid.x;
    y0 -= centroid.y; y1 -= centroid.y;
    z0 -= centroid.z; z1 -= centroid.z;
    
    // Renormalize all
    for(auto& do_obj : dysons) {
        do_obj.renormalize(x0, x1, y0, y1, z0, z1, step);
    }
    
    // Evaluate grid (Primary/Left orbital)
    UniformGrid grid(x0, x1, y0, y1, z0, z1, step);
    std::vector<double> data;
    data.reserve(grid.nx * grid.ny * grid.nz);
    
    const Dyson& primary_dyson = dysons[0];

    for (int ix = 0; ix < grid.nx; ++ix) {
        for (int iy = 0; iy < grid.ny; ++iy) {
            for (int iz = 0; iz < grid.nz; ++iz) {
                double x = grid.xmin + ix * grid.dx;
                double y = grid.ymin + iy * grid.dy;
                double z = grid.zmin + iz * grid.dz;
                data.push_back(primary_dyson.evaluate(x, y, z));
            }
        }
    }
    
    grid.write_binary(output_file, data);
    std::cout << "Wrote " << data.size() << " points to " << output_file << std::endl;
    // Print norms
    // Print norms
    std::cout << "Norm(L) [Grid]: " << dysons[0].normalization_factor 
              << ", [QChem]: " << dysons[0].qchem_norm << std::endl;
    if (dysons.size() > 1) {
        std::cout << "Norm(R) [Grid]: " << dysons[1].normalization_factor 
                  << ", [QChem]: " << dysons[1].qchem_norm << std::endl;
    }

    // XS Params
    double ie_ev, e_min, e_max;
    int l_max, n_pts;
    if (in >> ie_ev >> l_max >> n_pts >> e_min >> e_max) {
        std::cerr << "Computing Cross Sections..." << std::endl;
        std::cerr << "IE: " << ie_ev << " eV, L_max: " << l_max << std::endl;
        std::cerr << "Energy Grid: " << e_min << " - " << e_max << " eV (" << n_pts << " points)" << std::endl;
        
        std::vector<double> energies;
        if (n_pts > 1) {
            double step = (e_max - e_min) / (n_pts - 1);
            for(int i=0; i<n_pts; ++i) energies.push_back(e_min + i*step);
        } else {
            energies.push_back(e_min);
        }
        
        const Dyson& L = dysons[0];
        const Dyson& R = (dysons.size() > 1) ? dysons[1] : dysons[0];
        
        std::vector<double> xs = CrossSectionCalculator::ComputeTotalCrossSection(
            L, R, grid, energies, ie_ev, l_max
        );
        
        // Write XS to file "cross_section.txt"
        std::ofstream xs_file("cross_section.txt");
        xs_file << "# Energy(eV) CrossSection(au)\n";
        for(size_t i=0; i<energies.size(); ++i) {
            xs_file << energies[i] << " " << xs[i] << "\n";
        }
        xs_file.close();
        std::cerr << "Cross sections written to cross_section.txt" << std::endl;
    }
    
    return 0;
}
