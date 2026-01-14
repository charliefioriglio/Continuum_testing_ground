#ifndef _AOBasis_h
#define _AOBasis_h

#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "tools.h"
#include "gauss.h"
#include "xyzgrid.h"

#include <armadillo>

//!AO basis information
class AOBasis
{
  //! number of AOs
  int nbasis;
  int natoms;
  //! Array [nbasis] of contracted gaussians
  Gauss *basis;
  double *xgeom;    
  double *ygeom;    
  double *zgeom;    
  //! Number of basis functions per atom
  int *gtos_atom;
  int purecart[4];      
  bool is_lowdin_ready; 
  arma::mat X_mat;
  arma::mat X_min_mat;
  //Do not really need it, only for debug
  arma::mat S_mat;
  
  bool IfAlloc() const { return nbasis*natoms; }
  void Alloc(int nb, int na);
  void Free();

  AOBasis() : nbasis(0), natoms(0) {} 
  //  virtual AOBasis& operator=(const AOBasis& other) = 0;

 public:

  ~AOBasis() ;
  
  void IniAOBasis(const char* xmlFileName);
  void PrintAOBasis(FILE *fout=stdout);

  int NAtoms() const { return natoms; }
  int NBasis() const { return nbasis; }
  Gauss& GetGauss(int i) const { return basis[i]; }
  double CalcGaussValue(int ng, int nat, double x, double y, double z) const
    {
      //check(!(ng<nbasis && nat < natoms),"AOBasis::GetGaussValue() : invalid ng/na");
      return basis[ng].GaussValue(x,y,z,xgeom[nat],ygeom[nat],zgeom[nat]);
    }

  int GetGTOsAtom(int i) const { return gtos_atom[i]; }
  //Shift the basis in the new coordinate system 
  void ShiftGeom(double *newcenter);
  //Shift origin to the specified atom
  void ShiftGeomToAtom(int at);
  //Get position of atom 
  void GetGeom(int at, double *coord) const;
  //Calculate overlap
  void CalcOverlapMatrix(const XYZGrid& grid, arma::mat& Smatrix) const;
  //Calculate Lowdin's forward and backward transformations
  void CalcLowdinTransf(const XYZGrid& grid);
  
  arma::mat& GetLowdinTransf_X() { 

    //std::cout << "is_lowdin_ready=" << is_lowdin_ready << std::endl;
    check(!is_lowdin_ready,"AOBasis::GetLowdinTransf_X(): Not initialized");
    return X_mat;
  }
  arma::mat& GetLowdinTransf_Xm() {

    check(!is_lowdin_ready,"AOBasis::GetLowdinTransf_Xm(): Not initialized");
      return X_min_mat;
    }
  arma::mat& GetOverlap() {

    check(!is_lowdin_ready,"AOBasis::GetOverlap(): Not initialized");
      return S_mat;
  }
  
  friend AOBasis& TheAOBasis();
};

extern  AOBasis& TheAOBasis();


#endif
