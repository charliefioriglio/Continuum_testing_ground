#ifndef _cklmcoeff_h
#define _cklmcoeff_h

#include "complexno.h"
#include "xyzgrid.h"
#include "klmgrid.h"
#include "orbital.h"
#include "tools.h"
#include "wavetypes.h"

//! Cklm coefficients
class CklmCoeff
{
  int lmax;
  int nkv;
  //! Number of Cklm coefficients: nkv*nlmv.
  int ncklm;
  //! Matrix of all square coeffs, ||Cklm||^2
  double *cklmsq;
  //! Matrix of cross terms C*klmCkl'm. Dimensions [ncklm][lmax+1] 
  double **crosscklm;
  double *sumcklm;
    
  bool IfAlloc() const { return ncklm; }
  void Alloc(int nc, int lmx, int nk);
  void Free();
  void Set();

  CklmCoeff(const CklmCoeff& other);
  CklmCoeff& operator=(const CklmCoeff& other);
 
  // Calc x,y,z componenets of Cklm for subsequent averaging
  void CalcCklmSq(XYZGrid& labgrid, Orbital& dysorb,
                  double radfn,KLMPoints& klmgrid);

  // Calc Cklm for a particular molecular orientation and laser polarization
  void CalcCklmSq(XYZGrid& labgrid, Orbital& dysorb,
                  double *RIonz, double radfn,KLMPoints& klmgrid);
  // Calc overlap of Dyson orbital with plane wave for a particular molecular orientation 
  void CalcCklmSq_Overlap(XYZGrid& labgrid, Orbital& dysorb,
			  double radfn,KLMPoints& klmgrid);


public:
  CklmCoeff();
  
  CklmCoeff(XYZGrid& labgrid, Orbital& dysorb, double radfn, KLMPoints& klmgrid, bool do_overlap=false);

  CklmCoeff(XYZGrid& labgrid, Orbital& dysorb, double *RIonz, double radfn, KLMPoints& klmgrid);

  void IniCklmCoeff(XYZGrid& labgrid, Orbital& dysorb,
                    double radfn, KLMPoints& klmgrid);

  void IniCklmCoeff(XYZGrid& labgrid, Orbital& dysorb,
                    double *RIonz, double radfn, KLMPoints& klmgrid);

  void IniCklmCoeff_Overlap(XYZGrid& labgrid, Orbital& dysorb,
			    double radfn, KLMPoints& klmgrid);


  ~CklmCoeff();
  
  //! Print all avg square Cklm coeffs, ||Cklm^2||.
  void PrintCklmSq(KLMPoints& klmgrid, FILE *fout=stdout);
  void PrintAvgCklmSq(KLMPoints& klmgrid, FILE *fout=stdout);

  int GetNCklm() const { return ncklm; }
  int GetNKV() const { return nkv; }
  int GetLMax() const { return lmax; }
  //! Coeffs of k-th wavenumber and lm-th spherical harmonic, e.g., l=1, m=-1 has lm = 2.
  double GetCklmSq(int k, int lm) const { return cklmsq[k*(lmax+1)*(lmax+1)+lm]; }
  double GetCrossCklm(int k, int lm, int l2) const { return crosscklm[k*(lmax+1)*(lmax+1)+lm][l2]; } 
  double GetSumCklm(int k) const { return sumcklm[k]; }
  void PrintDysonScanZ(XYZGrid& labgrid, Orbital& dysorb, KLMPoints& klmgrid);

};
#endif

