#ifndef BETA_H
#define BETA_H

#include <vector>
#include "dyson.h"
#include "grid.h"
#include "angle_grid.h"

struct BetaResult {
    double energy;
    double sigma_par;
    double sigma_perp;
    double beta;
};

class BetaCalculator {
public:
    static std::vector<BetaResult> CalculateBeta(
        const Dyson& dyson_L,
        const Dyson& dyson_R,
        const UniformGrid& grid,
        const AngleGrid& angle_grid,
        const std::vector<double>& photoelectron_energies_ev
    );

private:
    static std::complex<double> ComputeNumericalMatrixElement(
        const Dyson& dyson,
        const UniformGrid& grid,
        const double* k_vec,
        const double* pol_vec
    );
};

#endif
