#include "num_eikr.h"
#include <cmath>
#include <iostream>
#include <cstring>
#include <array>

// Internal helper for ezDyson-compatible Rotation Logic
struct EzRotation {
    double mat[9];

    // Calc rotn which is Transpose of Active ZXZ Matrix
    // Matches ezDyson/rotnmatr.C: EulerRotnMatr(alpha, beta, gamma)
    void set_euler_zxz_transpose(double alpha, double beta, double gamma) {
        double cosA = std::cos(alpha), sinA = std::sin(alpha);
        double cosB = std::cos(beta), sinB = std::sin(beta);
        double cosC = std::cos(gamma), sinC = std::sin(gamma);

        // Row 0
        mat[0] = cosC*cosA - cosB*sinA*sinC; // Matches Active ZXZ R_00 ? No, R_00 = cAcG - cBsAsG. Same.
        mat[1] = cosC*sinA + cosB*cosA*sinC; // R_10? sAcG + cBcAsG. Matches.
        mat[2] = sinC*sinB;                  // R_20? sBsG. Matches.

        // Row 1
        mat[3] = (-sinC)*cosA - cosB*sinA*cosC; // R_01? -cA sG - cB sA cG. Matches R_01 of Active ZXZ?
                                                // R_01 = cA(-sG) - sA(cB cG) = -cAsG - cBsAcG. Matches.
        mat[4] = (-sinC)*sinA + cosB*cosA*cosC; // R_11? sA(-sG) + cA(cB cG). Matches.
        mat[5] = cosC*sinB;                     // R_21? sB cG. Matches.

        // Row 2
        mat[6] = sinB*sinA;      // R_02? sA sB. Matches.
        mat[7] = (-sinB)*cosA;   // R_12? -cA sB. Matches.
        mat[8] = cosB;           // R_22? cB. Matches.
        
        // Note: ezDyson implementation seems to store R^T relative to Active ZXZ?
        // Wait, if mat[1] (0,1) matches R_10, then it IS Transpose.
        // My check above: mat[1] matches R_10. 
        // So this matrix IS R^T.
        // And it transforms Lab -> Mol.
    }

    // Transform Lab -> Mol
    void transform(double xL, double yL, double zL, double& xM, double& yM, double& zM) const {
        xM = mat[0]*xL + mat[1]*yL + mat[2]*zL;
        yM = mat[3]*xL + mat[4]*yL + mat[5]*zL;
        zM = mat[6]*xL + mat[7]*yL + mat[8]*zL;
    }
};

NumEikr::NumEikr() {}
NumEikr::~NumEikr() {}

void NumEikr::compute(const UniformGrid& labgrid, 
                      const Dyson& dysonL, const Dyson& dysonR,
                      const AngleGrid& anggrid, 
                      const std::vector<double>& energies) {
                      
    std::vector<double> k_values;
    k_values.reserve(energies.size());
    // Convert E(eV) to k(a.u.)
    // E = k^2 / 2 => k = sqrt(2*E)
    // Careful with units: Input energies are in eV.
    // k (a.u.) = sqrt(2 * E_eV / 27.211386)
    // From tools.h: HAR_TO_EV = 27.211386
    const double HAR_TO_EV = 27.211386;
    for(double e : energies) {
        k_values.push_back(std::sqrt(2.0 * e / HAR_TO_EV));
    }
    
    // Allocate results
    cpar.assign(k_values.size(), 0.0);
    cperp.assign(k_values.size(), 0.0);
    
    calc_eikr_sq(labgrid, dysonL, dysonR, anggrid, k_values);
}

void NumEikr::calc_eikr_sq(const UniformGrid& labgrid, 
                           const Dyson& dysonL, const Dyson& dysonR,
                           const AngleGrid& anggrid, 
                           const std::vector<double>& k_values) {
                           
    int n_energies = k_values.size();
    int n_orientations = anggrid.points.size();
    
    // Grid parameters
    // Assuming labgrid describes the integration volume in Lab Frame
    // Just iterating points x,y,z
    // We assume uniform steps
    double dx = labgrid.dx;
    double dy = labgrid.dy;
    double dz = labgrid.dz;
    double dV = dx*dy*dz;
    double dV2 = dV * dV; // Precompute dV^2 for |Integral|^2?
    // Wait. ezDyson computes Integral(Psi * ...)^2.
    // Integral ~ Sum(val * dV).
    // |Integral|^2 ~ |Sum|^2.
    // ezDyson code: tmp = (SumL * SumR).Re().
    // Inside loop: Lcklm += val * dV? No.
    // In eikr.C loop: Lcklm += totalLDys_par (which is val).
    // Then after loop: tmp_par = (SumL * SumR).Re() * dVV * tmpavg.
    // dVV = dV * dV. This matches |Integral*dV|^2.
    
    double x0 = labgrid.xmin;
    double y0 = labgrid.ymin;
    double z0 = labgrid.zmin;
    int nx = labgrid.nx;
    int ny = labgrid.ny;
    int nz = labgrid.nz;
    
    // Thread-local accumulation (since we parallelize over orientations? Or just simple loop)
    // We will loop over orientations and Inside loop over energies.
    // Results accumulate to cpar[k], cperp[k].
    
    EzRotation rot;
    
    for (int v = 0; v < n_orientations; ++v) {
        const auto& p = anggrid.points[v];
        // ezDyson uses numeric averaging with RotnMatr(NUM, 0, pi.beta, pi.gamma)
        // See eikr.C: "RotnMatr rot(avgcklm,0,bj,aj);" where bj=Beta, aj=Gamma from my AngleGrid analysis.
        // My AngleGrid stores: alpha=0, beta=p.beta, gamma=p.gamma.
        // So I pass these to the rotation matrix.
        // The rotation is Active ZXZ Transpose (for Lab -> Mol).
        rot.set_euler_zxz_transpose(p.alpha, p.beta, p.gamma);
        double weight = p.weight;
        
        // Pre-allocate sums for all energies for this orientation
        std::vector<std::complex<double>> sumL_par(n_energies, 0.0);
        std::vector<std::complex<double>> sumR_par(n_energies, 0.0);
        std::vector<std::complex<double>> sumL_x(n_energies, 0.0);
        std::vector<std::complex<double>> sumR_x(n_energies, 0.0);
        std::vector<std::complex<double>> sumL_y(n_energies, 0.0);
        std::vector<std::complex<double>> sumR_y(n_energies, 0.0);
        
        // Loop Spatial Grid
        for(int ix=0; ix<nx; ++ix) {
            double xL = x0 + ix*dx;
            for(int iy=0; iy<ny; ++iy) {
                double yL = y0 + iy*dy;
                for(int iz=0; iz<nz; ++iz) {
                    double zL = z0 + iz*dz;
                    
                    // Rotate Coords Lab -> Mol
                    double xM, yM, zM;
                    rot.transform(xL, yL, zL, xM, yM, zM);
                    
                    // Evaluate Dyson
                    double valL = dysonL.evaluate(xM, yM, zM);
                    double valR = dysonR.evaluate(xM, yM, zM);
                    // L/R Logic: In clean code, we might treat them same if same object
                    
                    // Dipole Operator (r . eps) in Lab Frame
                    // ezDyson fixes Polarization along Z (eps = z).
                    // So Operator is always zL.
                    // Par: k || Z. Perp: k || X (or Y).
                    
                    std::complex<double> termL_par = valL * zL;
                    std::complex<double> termR_par = valR * zL;
                    
                    // Perp X: Dipole Z, k || X
                    std::complex<double> termL_x = valL * zL; 
                    std::complex<double> termR_x = valR * zL;
                    
                    // Perp Y: Dipole Z, k || Y
                    std::complex<double> termL_y = valL * zL;
                    std::complex<double> termR_y = valR * zL;
                    
                    // Loop Energies
                    for(int k=0; k<n_energies; ++k) {
                        double kval = k_values[k];
                        // Plane Wave part: exp(i * k * r_lab . k_hat)
                        // k vector is along Z_lab?
                        // Standard Photoelectron definition: k is the direction of electron emission.
                        // In PAD formalism, usually we align k along Z_lab (or define axis wrt k).
                        // ezDyson eikr.C:
                        // double coskz = cos(k * z); ... eikrz = Complex(coskz, sinkz);
                        // totalLDys_par = tLDys_cart * eikrz;
                        // totalLDys_x = tLDys_cart * eikrx;
                        // totalLDys_y = tLDys_cart * eikry;
                        // It seems ezDyson calculates scattering into X, Y, Z directions simultaneously?
                        // And then averages?
                        // Wait. usually we pick k || Z_lab.
                        // cpar corresponds to polarization || k? No.
                        // Parallel/Perp refer to Polarization vector relative to some axis?
                        // Usually Sigma_Par is Pol || k? Or Pol || Molecular Axis?
                        // In Lab Frame fixed experiment (e.g. photodetachment), 
                        // Beta is defined relative to Laser Polarization.
                        // Theta is angle between k and Polarization epsilon.
                        // Here we are calculating Total Cross Section or Beta?
                        // Beta from Sigma_par and Sigma_perp.
                        // Sigma_par: electron emitted along polarization?
                        // ezDyson: cpar computes integral with eikrz (k along Z).
                        // tLDys_cart uses z (Dipole along Z).
                        // So Par: k || Z, eps || Z. (Theta = 0).
                        // cperp: eikrz (k along Z).
                        // But uses x and y dipoles?
                        // No, ezDyson uses eikrX, eikrY, eikrZ.
                        // eikr.C: 
                        // totalLDys_par = tLDys_cart * eikrz; (Dipole Z, k Z) -> Theta=0
                        // totalLDys_x = tLDys_cart * eikrx; (Dipole Z, k X) -> Theta=90? (Pol Z, k X)
                        // Wait. tLDys_cart = Ldys * z. (This is dipole Z).
                        // So Polarization is FIXED along Z.
                        // Then it computes emission along Z (eikrz) -> Par.
                        // And emission along X (eikrx) -> Perp?
                        // And emission along Y (eikry) -> Perp?
                        // Yes!
                        // "cperp[k] += 0.5*tmp_x+0.5*tmp_y;"
                        // So it calculates k along X and k along Y (both perp to Pol Z).
                        // This assumes averaging over orientations makes X and Y redundant?
                        // Yes.
                        
                        // Implementation:
                        // Par: Dipole Z, PlaneWave exp(i k zL).
                        // Perp: Dipole Z, PlaneWave exp(i k xL) (and yL).
                        
                        // Wait, my code above:
                        // termL_par = valL * zL; (Dipole Z)
                        // termL_x = valL * xL; (Dipole X?)
                        // ezDyson: `tLDys_cart=Ldys_value*gridptr_z[nz];` (Dipole Z is FIXED).
                        // Then `totalLDys_par = tLDys_cart * eikrz`.
                        // `totalLDys_x = tLDys_cart * eikrx`.
                        // Ah! It uses the SAME dipole operator (Z-polarized light) for all.
                        // It scans k-vector direction.
                        // k || Z -> Par. k || X -> Perp.
                        
                        // Correct Logic:
                        std::complex<double> dipole_op = valL * zL; // Polarization along Z
                        std::complex<double> dipole_op_R = valR * zL;
                        
                        double kz = kval * zL;
                        double kx = kval * xL;
                        // double ky = kval * yL; // Optional for averaging
                        
                        std::complex<double> exp_kz(std::cos(kz), std::sin(kz));
                        std::complex<double> exp_kx(std::cos(kx), std::sin(kx));
                        
                        // Accumulate
                        // Integral Psi^* * Op * Phi
                        // Psi_k = exp(i k r). Psi^* = exp(-i k r).
                        // Integral exp(-i k r) * z * Dyson
                        // My sumL should sum (dyson * z * exp(-ikz)).
                        // termL * conj(exp_kz).
                        // ezDyson: `totalLDys_par=tLDys_cart*eikrz;` (No conj? Maybe eikrz is exp(-ikz)?)
                        // `Complex eikrz(coskz,sinkz);` -> exp(ikz).
                        // Integration: `Lcklm += total`.
                        // Later: `tmp_par=(Lcklm * Rcklm).Re()`.
                        // If one is conjugate?
                        // Matrix Element = <Psi | Op | Dyson> = Integral Psi^* Op Dyson.
                        // If Psi = exp(ikz), Psi^* = exp(-ikz).
                        // ezDyson accumulates Dyson * exp(ikz).
                        // This implies M = Integral Dyson * exp(ikz).
                        // This corresponds to < exp(-ikz) | Op | Dyson >? Or < Dyson | Op | exp(-ikz) >* ?
                        // The physics: M ~ Fourier Transform of (Op * Dyson).
                        // FT(f)(k) = Integral f(r) exp(-ik r).
                        // ezDyson uses exp(+ikr).
                        // Maybe definition of k is -k? Or it calculates < Dyson | Op | Psi >?
                        // If <D | Op | Psi> = Integral D^* z exp(ikz).
                        // Then |M|^2 is same.
                        // I will assume `exp(ikz)` is correct matching ezDyson.
                        
                        sumL_par[k] += dipole_op * exp_kz;
                        sumR_par[k] += dipole_op_R * std::conj(exp_kz); // Wait.
                        // ezDyson: `totalRDys_par=tRDys_cart*eikrz.Conj();`
                        // So R uses Conj(exp). L uses exp.
                        // Then product L * R corresponds to |Integral|^2?
                        // Integral L * Integral R^*.
                        // If SumL = Sum(D * exp), SumR = Sum(D * exp^*).
                        // SumL * SumR = (Sum D exp) * (Sum D exp^*).
                        // This is NOT |Sum|^2 unless exp^* = conj(exp).
                        // If SumR uses conj(exp), and if coefficients are real...
                        // If SumL = A, SumR = A^*. Then product is |A|^2.
                        // ezDyson uses `totalRDys_par=tRDys_cart*eikrz.Conj()`.
                        // So yes, R takes conjugate.
                        
                        sumL_par[k] += dipole_op * exp_kz;
                        sumR_par[k] += dipole_op_R * std::conj(exp_kz);
                        
                        sumL_x[k] += dipole_op * exp_kx;
                        sumR_x[k] += dipole_op_R * std::conj(exp_kx);
                        
                        // Ignore Y for speed if X is statistically sufficient (ezDyson does X and Y)
                        // I'll do just X for now to match 2x speedup or do both to match accuracy.
                        // ezDyson does "0.5*tmp_x+0.5*tmp_y". I'll do just X and assume symmetry or do Y.
                        // Let's do Y for completeness.
                        double ky = kval * yL;
                        std::complex<double> exp_ky(std::cos(ky), std::sin(ky));
                        sumL_y[k] += dipole_op * exp_ky;
                        sumR_y[k] += dipole_op_R * std::conj(exp_ky);
                    }
                }
            }
        } // End Grid
        
        // Accumulate to global cross sections
        for(int k=0; k<n_energies; ++k) {
            double term_par = std::real(sumL_par[k] * sumR_par[k]);
            double term_x = std::real(sumL_x[k] * sumR_x[k]);
            double term_y = std::real(sumL_y[k] * sumR_y[k]);
            
            cpar[k]  += term_par * dV2 * weight;
            cperp[k] += 0.5 * (term_x + term_y) * dV2 * weight;
        }
        
    } // End Orientations
}
