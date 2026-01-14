#ifndef _radialfunction_h
#define _radialfunction_h

#include "tools.h"
#include "klmgrid.h"
#include "xyzgrid.h"
//#include "gsl/gsl_sf_bessel.h"
//#include "gsl/gsl_sf_coulomb.h"
#include "readwrite.h"
#include <math.h>
#include "wavetypes.h"


//!Radial part of the spherical waves, Rkl
class RadialFunction
{
  //! Type of Rkl functions (Spherical, Coulomb waves).
  RadFn radfn;
  //! Number of Rkl functions.
  int nrkl;
  //! Number of k values.
  int nkv;
  //! l max for which to calc Rkl.
  int lmax;
  //! Number of (x,y,z) grid points.
  int nxyz;
  //! Matrix of all Rkl's at all x,y,z grid points
  double **rkl;

  bool IfAlloc() const { return nrkl*nxyz; }
  void Alloc(int nk, int lmx, int nxyzpts);
  void Free();
  void Set();
  void IniRadialFunction(const char* xmlFileName, const XYZGrid& labgrid, const KLMPoints& klmgrid);
  //void IniRadialFunction(RadFn rfn, int nk, int lmx, int nxyzpts, double** const rkl, double** const rklij);
  //! Calc all Rkl values (nrkl*nxyz) at x,y,z grid points.
  void CalcAllRkl(const XYZGrid& labgrid, const KLMPoints& klmgrid);

  //! Calc overlaps between different l, same k radial functions.
  void CalcRklOverlaps(const KLMPoints& klmgrid);
  //double GetRklOverlap(int k, int l1, int l2) const { return rklij[k*(lmax+1)+l1][l2]; }


  RadialFunction(const RadialFunction& other);
  RadialFunction& operator=(const RadialFunction& other);

 public:
	
  RadialFunction();
  RadialFunction(const char* xmlFileName, const XYZGrid& labgrid, const KLMPoints& klmgrid);
  ~RadialFunction();
  
  //! Print Rkl matrix, for all k,l and x,y,z values.
  void Print(const char *filename) const;
  //! Calc norm matrix RklRkl'
  void NormMatrix(double rmax, int nrv, const KLMPoints& klmgrid);

  int NXYZ() const { return nxyz; }
  int NKV() const { return nkv; }
  int NRkl() const { return nrkl; }
  int LMax() const { return lmax; }
  inline double& GetRkl_xyz(int kl, int i) const { return rkl[kl][i]; }
  
  //double* GetRklPtr(int k, int l) { return rkl[k*(lmax+1)+l]; }; 
};

#endif
