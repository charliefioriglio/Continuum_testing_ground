#ifndef _cklmcoeff_h
#define _cklmcoeff_h

#include "complexno.h"
#include "xyzgrid.h"
#include "klmgrid.h"
#include "radialfns.h"
#include "orbital.h"
#include "tools.h"


//! Cklm coefficients
class CklmCoeff
{
  int lmax;
  int nkv;
  //! Number of Cklm coefficients: nkv*nlmv.
  int ncklm;
  //! Matrix of all square coeffs, ||Cklm||^2
  double **cklmsq;
  //! Matrix of cross terms C*klmCkl'm. Dimensions [ncklm][lmax+1] 
  double ***crosscklm;
  double **sumcklm;

  bool IfAlloc() const { return ncklm; }
  void Alloc(int nc, int lmx, int nk);
  void Free();
  void Set();

  CklmCoeff(const CklmCoeff& other);
  CklmCoeff& operator=(const CklmCoeff& other);
 
  // Calc x,y,z componenets of Cklm for subsequent averaging
  void CalcCklmSq(XYZGrid& labgrid, Orbital& dysorb,
                  RadialFunction& allrkl,KLMPoints& klmgrid);

  // Calc Cklm for a particular molecular orientation and laser polarization
  void CalcCklmSq(XYZGrid& labgrid, Orbital& dysorb,
                  double *RIonz, RadialFunction& allrkl,KLMPoints& klmgrid);

  
public:
  CklmCoeff();
  
  CklmCoeff(XYZGrid& labgrid, Orbital& dysorb, RadialFunction& allrkl, KLMPoints& klmgrid);

  CklmCoeff(XYZGrid& labgrid, Orbital& dysorb, double *RIonz, RadialFunction& allrkl, KLMPoints& klmgrid);

  void IniCklmCoeff(XYZGrid& labgrid, Orbital& dysorb,
                    RadialFunction& allrkl, KLMPoints& klmgrid);

  void IniCklmCoeff(XYZGrid& labgrid, Orbital& dysorb,
                    double *RIonz, RadialFunction& allrkl, KLMPoints& klmgrid);

  ~CklmCoeff();
  
  //! Print all avg square Cklm coeffs, ||Cklm^2||.
  void PrintCklmSq(KLMPoints& klmgrid, FILE *fout=stdout);
  void PrintAvgCklmSq(KLMPoints& klmgrid, FILE *fout=stdout);

  int GetNCklm() const { return ncklm; }
  int GetNKV() const { return nkv; }
  int GetLMax() const { return lmax; }
  //! Coeffs of k-th wavenumber and lm-th spherical harmonic, e.g., l=1, m=-1 has lm = 2.
  double GetCklmSq(int v, int k, int lm) const { return cklmsq[v][k*(lmax+1)*(lmax+1)+lm]; }
  double GetCrossCklm(int v, int k, int lm, int l2) const { return crosscklm[v][k*(lmax+1)*(lmax+1)+lm][l2]; } 
  double GetSumCklm(int v, int k) const { return sumcklm[v][k]; }
  void PrintDysonScanZ(XYZGrid& labgrid, Orbital& dysorb, RadialFunction& allrkl, KLMPoints& klmgrid);
};
#endif
