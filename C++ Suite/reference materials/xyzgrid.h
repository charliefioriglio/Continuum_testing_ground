#ifndef _xyzgrid_h
#define _xyzgrid_h

#include "tools.h"
#include "aik_xml_parser.h"

//general x,y,z grid in a.u.
class XYZGrid
{
  //! Number grid points (x,y,z).
  int nxyz;
  //! Number of points along X, Y, Z
  int npoints[XYZ];
  double xyzmin[XYZ];
  double xyzmax[XYZ];
  //! The grid points for XYZ
  double *thegrid[XYZ];
  //! Increments (deltas) for x,y,z
  double dxyz[XYZ];     

  bool IfAlloc() const { return nxyz; }
  void Alloc(const int *const npts);

  void Free();
  void Set();
  //!Calc dx,dy,dz for the x,y,z grid.
  void CalcDxDyDz();
  //!Generate x,y and z values for *thegrid[XYZ].
  void CalcXYZGrid();

  //! To make sure noone uses it!
  XYZGrid& operator=(const XYZGrid& other);
 
public:

  XYZGrid();
  // XYZGrid(int *npts, const double* min, const double* max);
  XYZGrid(const char* xmlFileName);
  XYZGrid(const XYZGrid& other); 
  ~XYZGrid();

  //  void IniXYZGrid(FILE *input);

  //! Shift the center of the XYZGrid by xyzcenter[i],[j],[k].
  void ShiftXYZGrid(double *xyzcenter);
  //! Print the grid info to std output.
  void PrintGridInfo(FILE *fout=stdout) const;
  //!Print to file the grid as thegrid[XYZ][nxyz].
  void Print(const char *filename) const;

  //! No. xyz grid points.
  int NXYZ() const { return nxyz; }
  //! Returns grid dimensions, e.g., NPoints(0)=nx
  int NPoints(Carts xyz) const { return npoints[xyz]; }
  //! Returns point from the grid, e.g., GetPoint(0,i) returns x[i]
  double GetPoint(Carts xyz, int i) const { return thegrid[xyz][i]; }
  /*! Returns a point from the grid by absolute address. 
      axis1,2,3 give the ordering of grid pts by axis, e.g., X,Y,Z; Z,X,Y; Y,X,Z; ...
      pt1,2,3 are the points on each axis, e.g., correspondingly: (x,y,z); (z,x,y); (y,x,z); ... */
  void GetPoint(int absaddr, double &pt1, double &pt2, double &pt3, Carts axis1, Carts axis2, Carts axis3) const;
  //!Spacing along x, y or z.
  double DXYZ(Carts xyz) const { return dxyz[xyz]; }
  double *GetGridPtr(Carts xyz) const { return thegrid[xyz]; }
};

#endif




