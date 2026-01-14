#ifndef _orbital_h
#define _orbital_h

#include "gauss.h"
#include "aobasis.h"
#include "xyzgrid.h"

//!Coeffs of any orbital(left-right) in AO basis.
class Orbital
{
  //! Number of AO basis functions.
  int nbasis;
  //! Restr or unrestr calc in QChem?
  bool ifrestr;
  //! N-th ionization transition in QChem.
  int transno;
  //! Coefficients of orbital in AO basis(Left and Right Dyson).
  double *Lcoeff;   
  double *Rcoeff;
  //! Keeps track of small values of AO coefficients
  short *zeromap;
  //! left and right Dyson MO norms 
  double Lnorm;
  double Rnorm;

  bool IfAlloc() const { return nbasis; }
  void Alloc(int nb);
  void Free();
  void Set();
  void AllocAndCopyData(int nb, bool ifr, int tran, double ln, double rn, double *Lc, double *Rc, short *zm);
  void SetZeroMap(double thresh=0.0001); 
 public:
	
  Orbital();
  Orbital(const Orbital& other);
  Orbital(int nb, bool ifr, int tran, double rn, double ln, double *Lc, double *Rc);
  //Make a Dyson orbital from other by retaining only the part on a given atom 
  Orbital(const Orbital& other, int atom);
  ~Orbital();
  

  Orbital& operator=(const Orbital& other);

  //! Print Orbital coeffs in AO basis.
  void Print(FILE *fout=stdout) const;
  //! Add up all basis fns with proper coefs to get orbital value at x,y,z.
  void CalcDyson_xyz(double &Ldys_value, double &Rdys_value, double x, double y, double z) const;
  //! Calc left Dyson orbital norm and center
  double GetOrbitalNormAndCenter(double *center, const XYZGrid& grid) const;
  //!  Calculate overlap with another orbital
  void GetOverlap(const Orbital& other, const XYZGrid& grid, double& s_right_right, double& s_left_left) const;

  double GetLDysonNorm() const { return Lnorm; }
  double GetRDysonNorm() const { return Rnorm; }
  //! This is actually norm^2.
  double GetLRNorm() const { return Lnorm*Rnorm; }
  //! Get coeffs of L,R  Dyson MO in AO basis
  double GetLDysAO(int n) const { return Lcoeff[n]; }
  double GetRDysAO(int n) const { return Rcoeff[n]; }
};

#endif
