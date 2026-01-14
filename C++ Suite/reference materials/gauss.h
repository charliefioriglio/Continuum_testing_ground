#ifndef _gauss_h
#define _gauss_h
#include "tools.h"


//!Contracted gaussian.
class Gauss
{
  //! number of contracted gaussians
  int ncontr;
  //! Total angular momentum
  int angmom;
  //! x, y, z powers (define l,m)
  int ia, ja, ka;  
  //! 2 in case we have pure d gauss
  int dx2y2, dz2; 
  //! 3 in case we have pure f gauss
  int puref;
  //! which pure f? 1-4,6,7; fxyz is fa=5
  int fa;
  //! 4 in case we have pure g gauss
  int pureg;
  //! which pure g? 1-9.
  int ga;
  //! exponents
  double *alpha;    
  //! coefficients
  double *coeff;   
  //! Norms of individual uncontracted gaussians, not sure if we need them
  double *norm;      
  //! precomputed coeff*norm for each uncontracted GTO
  double *coeffnorm;  

  //! Checks if it is allocated
  bool IfAlloc() const { return ncontr; }
  void Alloc(int nc);
  void Free();

 public:
  Gauss() : ncontr(0) {};
  Gauss(const Gauss& other);
  ~Gauss() ;
  
  void IniGauss(int angmom, int nc, int i, int j, int k, int dxy, int dz, int puref, int f, int pureg, int g, double *a, double *c, double sc);
  Gauss& operator=(const Gauss& other);

  //! Calc norm of contracted gaussian.
  double TotalNorm() const;
  void GaussNorm() const;
  //! Precompute coeff*norm for each uncontracted GTO.
  void GaussCoeffNorm() const;

  //! Calc value of gaussian at x,y,z.
  double GaussValue(double x, double y, double z, double x0, double y0, double z0) const;

  //! Angular momentum
  int AngMom() const { return angmom; }
  //! Checks for consistency of angular momentum
  void CheckAngMom() const;
  //! No. contracted gaussians.
  int NContr() const { return ncontr; }
  //! Exponents matrix. 
  double Alpha(int j) const { return alpha[j]; }
  //! Contraction coeffs matrix.
  double Coeff(int j) const { return coeff[j]; }
  //! Renormalized contraction coeffs matrix.
  double CoeffNorm(int j) const { return coeffnorm[j]; }
  //! Norm of contracted gaussian.
  double Norm(int j) const {return norm[j]; }
  //Powers of x,y,z for cartesian gaussians.
  int Ia() const { return ia; }
  int Ja() const { return ja; }
  int Ka() const { return ka; }
  //!Pure d fns? Yes(2), No(0).
  int Dx2y2() const { return dx2y2; }
  int Dz2() const { return dz2; }
  int PureF() const { return puref; }
  int Fa() const { return fa; }
};


//!Calculates norm of uncontracted cart gaussians.
double CartNorm(int ib, int jb, int kb, double alph);
double PureNorm(int lpure, double alph);



#endif
