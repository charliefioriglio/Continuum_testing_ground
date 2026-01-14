#ifndef _rotnmatr_h
#define _rotnmatr_h

#include "tools.h"
#include "anglegrid.h"

typedef enum {NOAVG,AVG,NUM} MolAvg;


//! Class that contains all matrices for mol/lab frame conversion
class RotnMatr
{
  //! Array of rotn matrices(3x3) for each angle grid point: [nabc][9]
  double rotn[XYZ*XYZ];

  //! Calc rotation matrix for cyllindrical rotation (one angle: alpha)..
//  void CalcRotnMatr(MolAvg moav, double alpha);
  //! Calc rotation matrix (Euler transformation) for alpha, beta, gamma Euler angles.
  void CalcRotnMatr(double aplha, double beta, double gamma);

 public:
  RotnMatr();
  RotnMatr(const RotnMatr& other);
  RotnMatr(MolAvg moav, double alpha);
  RotnMatr(MolAvg moav, double alpha, double beta, double gamma);
  
  RotnMatr& operator=(const RotnMatr& other);

  double GetRotnMatrElem(int i) const { return rotn[i]; }
  double& operator[](int i) { return rotn[i];}

  //double* GetPtrToRotn() { return rotn; }
  void Print(const char *filename) const;

  //! Calc rotn 3x3 Euler rotation matrix for a a certain molec orientation: alpha,beta,gamma.
  RotnMatr& EulerRotnMatr(double alpha, double beta, double gamma);
  //! Calc rotn 3x3 rotation matrix around x by alpha.
  RotnMatr& XRotnMatr(double alpha);
  //! Calc rotn 3x3 rotation matrix around y by alpha.
  RotnMatr& YRotnMatr(double alpha);
  //! Calc rotn 3x3 rotation matrix around z by alpha.
  RotnMatr& ZRotnMatr(double alpha);
  //! Calc determinant of 3x3 rotation matrix.
  double DetRotnMatr() const;
  //! Get inverse of rotation matrix
  RotnMatr GetInvMatr() const;
  //! Convert x coord from molec to lab frame or the reverse 
  double XLabMol( double xL, double yL, double zL) const;
  //! Convert y coord from molec to lab frame or the reverse 
  double YLabMol(double xL, double yL, double zL) const;
  //! Convert z coord from molec to lab frame or the reverse 
  double ZLabMol(double xL, double yL, double zL) const;
  //! Convert x,y,z from molec to labe frame or the reverse.
  inline void GenLabMol(double &xM, double &yM, double &zM, double xL, double yL, double zL) const
    {
      xM = rotn[0]*xL+rotn[1]*yL+rotn[2]*zL;
      yM = rotn[3]*xL+rotn[4]*yL+rotn[5]*zL;
      zM = rotn[6]*xL+rotn[7]*yL+rotn[8]*zL;
    }
  
  //! Convert x,y,z from molec to labe frame or the reverse, for cyllindrical avg around X.
  inline void CylXLabMol(double &xM, double &yM, double &zM, double xL, double yL, double zL) const
    {
      xM = xL;
      yM = rotn[4]*yL+rotn[5]*zL;
      zM = rotn[7]*yL+rotn[8]*zL;
    }

  //! Convert x,y,z from molec to labe frame or the reverse, for cyllindrical avg around Y.
  inline void CylYLabMol(double &xM, double &yM, double &zM, double xL, double yL, double zL) const
    {
      xM = rotn[0]*xL+rotn[2]*zL;
      yM = yL;
      zM = rotn[6]*xL+rotn[8]*zL;
    }

  //! Convert x,y,z from molec to labe frame or the reverse, for cyllindrical avg around Z.
  inline void CylZLabMol(double &xM, double &yM, double &zM, double xL, double yL, double zL) const
    {
       xM = rotn[0]*xL+rotn[1]*yL;
       yM = rotn[3]*xL+rotn[4]*yL;
       zM = zL;
    }
    
  //! Choose which function to use for cylindrical rotn of X,Y,Z axes.
  void CylLabMol(MolAvg molavg, double &xM, double &yM, double &zM, double xL, double yL, double zL) const;

  //! Convert dipole moment in lab frame 
  void DipMomInLab(double *L_M, double *M) const;
    
};

/*! Form matrix of molec orientation distr before ionization, after polarized light excitation.
  Averaging function for all (alpha,beta,gamma) orientations(nang).*/
void RempiAvgFn(double *avgfn, double *Mex, double *Rex, RotnMatr *inva, int nang);
/*! Excitation probability for a molecule with orientation corresponding to dipole mom in lab frame,L_Mex[x,y,z],
  by excitation laser with polarization Rex[x,y,z]. */
double RempiProb(double *L_Mex, double *Rex);
double CylAvgfunction();
double AvgFunction(int molavg, double alpha, double beta, double gamma);
#endif
