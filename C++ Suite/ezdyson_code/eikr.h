#ifndef _eikr_h
#define _eikr_h

#include "complexno.h"
#include "xyzgrid.h"
#include "klmgrid.h"
#include "orbital.h"
#include "tools.h"
#include "rotnmatr.h"
#include "anglegrid.h"

//! Numerical Eikr intergrals
class NumEikr
{
  int nkv;
  double *cpar;
  double *cperp;

  bool IfAlloc() const { return nkv; }
  void Alloc(int nk);
  void Free();
  void Set();

  NumEikr(const NumEikr& other);
  NumEikr& operator=(const NumEikr& other);
 
//! Compute cross sections with numerical REPULSION averaging
  void CalcEikrSq(XYZGrid& labgrid, Orbital& dysorb,
                  AngleGrid& anggrid,KLMPoints& klmgrid);

public:
  NumEikr();
  
  NumEikr(XYZGrid& labgrid, Orbital& dysorb, AngleGrid& anggrid, KLMPoints& klmgrid);

  void IniNumEikr(XYZGrid& labgrid, Orbital& dysorb,
                    AngleGrid& anggrid, KLMPoints& klmgrid);

  ~NumEikr();
  
  int GetNKV() const { return nkv; }
  //! Coeffs of k-th wavenumber and lm-th spherical harmonic, e.g., l=1, m=-1 has lm = 2.
  double GetCPar(int k) const { return cpar[k]; }
  double GetCPerp(int k) const { return cperp[k]; }
  void PrintDysonScanZ(XYZGrid& labgrid, Orbital& dysorb, KLMPoints& klmgrid);

};
#endif

