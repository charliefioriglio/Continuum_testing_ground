#ifndef CROSS_SECTION_H
#define CROSS_SECTION_H

#include "dyson.h"
#include "grid.h"
#include <vector>
#include <string>

// Computes total cross section for a given Dyson orbital and energy range
class CrossSectionCalculator {
public:
    // Calculates sigma for a list of PHOTON energies (eV)
    static std::vector<double> ComputeTotalCrossSection(
        const Dyson& dyson_L,
        const Dyson& dyson_R,
        const UniformGrid& grid,
        const std::vector<double>& photon_energies_ev,
        double ionization_energy_ev,
        int l_max
    );

private:

};

#endif
