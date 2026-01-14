#ifndef _anglegrid_h
#define _anglegrid_h

#include "tools.h"

//! The dimension of Euler angles space
#define ABG 3
typedef enum {A,B,G} Angles;

//! general alpha, beta, gamma Euler angles grid
class AngleGrid
{
  //! Number grid points (a,b,c).
  int nabc;
  //! Number of points for rotation around each axis.
//  int npoints[ABG];
//  double abcmin[ABG];
//  double abcmax[ABG];
  //! The Euler angles points 
  double *thegrid[ABG];
  //! Increments (deltas) for alpha, beta, gamma
//  double dabc[ABG];     

  bool IfAlloc() const { return nabc; }
  void Alloc();
  void Free();
  void Set();

  //!Calc da,db,dc for the angles grid.
//  void CalcDaDbDc();
  //!Generate alpha, beta, gamma values for *thegrid[ABG].
  void CalcAngleGrid();

  AngleGrid& operator=(const AngleGrid& other);
  AngleGrid(const AngleGrid& other);

 public:
  AngleGrid();
//  AngleGrid(const int *const nabc, const double *min, const double *max);
  ~AngleGrid();

  //! Print grid info
//  void PrintGridInfo(FILE *outfile=stdout) const;
  //!Print to file the grid as thegrid[ABG][nabc].
//  void Print(const char *filename) const;

  //! No. alpha, beta, gamma grid points.
  int NABC() const { return nabc; }
  //! Returns grid dimensions, e.g., NPoints(0)=nalpha
//  int NPoints(Angles abc) const { return npoints[abc]; }
  //! Returns point from the grid, e.g., GetPoint(0,i) returns alpha[i]
  double GetPoint(Angles abc, int i) const { return thegrid[abc][i]; }
  //! Returns value of alpha for cyllindrical averaging.
//  double GetPoint(int i) const { return thegrid[A][i]; }
  //! Returns a point from the grid by absolute address
//  void GetPoint(int absaddr, double& a, double &b, double &c) const;
  //!Spacing for alpha, beta, gamma.
//  double DABC(Angles abc) const { return dabc[abc]; }

  //!Calc angle points for isotropic averaging.
  //friend void CalcIsoGrid(const AngleGrid&, double** isogrid, double* dbeta);
};

//! Grid for isotropic averaging, it is a bit tricky...
/*
class ISOGrid
{
  int nabc;
  double *isogrid[ABG];
  double *dbeta;
  
 public:
  ISOGrid(const AngleGrid& anggrid);
  ~ISOGrid();
  
  inline int NABC() const { return nabc; }
  inline double GetDBeta(int pnt) const { return dbeta[pnt]; }
  inline double GetPoint(Angles Ang, int pnt) const { return isogrid[Ang][pnt]; }
};
*/

//!Pointer to Euler angles grid (molec orientations in lab frame) in input.
//long PtrToAnglesGrid(FILE *input);

//!Read npoints[ABG] from input file & Euler angles range
//void ReadAllAngleGridInfo(int* npts, double* abcmin, double* abcmax, const char* xmlFileName);

#endif




