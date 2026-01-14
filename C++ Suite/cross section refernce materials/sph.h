#ifndef _sph_h
#define _sph_h

#include "complexno.h"


//! Class for fast calculation of spherical functions, i^l Ylm*eimf, normalized to 1.
class SPH
{
  //! Max angular momentum implemented
  int lmax; 
  //precomputed coefs
  double **thetacoeff;
  
  //! Private constructor
  SPH(int lmax=11);
 
public:
  ~SPH();

  Complex GetValue(double x, double y, double z, int lv, int mv) const;
  int LMax() const { return lmax; }
  
  //! Access function
  friend SPH& theSPH();
};

#endif
