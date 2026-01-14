#include "xyzgrid.h"


void XYZGrid::Alloc(const int *const npts)
{
  check((!(npts[X]>0) && !(npts[Y]>0) && !(npts[Z]>0)),"XYZGrid::Alloc() : Invalid nxyz\n");
  Free();

  double gsize=0.0;

  memcpy(npoints, npts, sizeof(int)*XYZ);
  nxyz = npoints[X]*npoints[Y]*npoints[Z];
  for (int i=0; i<XYZ; i++)
    {
      thegrid[i] = new double[npoints[i]];
      gsize+=sizeof(double)*npoints[i];
      check(!(thegrid[i]>NULL), "XYZGrid::Alloc() : Allocation failed");
    }
  
  printf("XYZGrid::Alloc:  %lf MB allocated\n",gsize/MEGABITE);
}

void XYZGrid::Free()
{
  if(IfAlloc())
    {
      for (int i=0; i<XYZ; i++)
	if(NULL!=thegrid[i])     
	  {
	    delete [] thegrid[i]; 
	    thegrid[i]=NULL;
	  }
      memset(npoints,0,XYZ*sizeof(int));
      nxyz=0;
    }
}

void XYZGrid::Set()
{  
  nxyz=0;
  memset(npoints,0,XYZ*sizeof(int));
  for (int i=0; i<XYZ; i++)
    thegrid[i]=NULL;
}

XYZGrid::XYZGrid() 
{
  Set();
}

XYZGrid::XYZGrid(const XYZGrid& other) 
{
  Set();
  Alloc(other.npoints);
  memcpy(xyzmin,other.xyzmin,sizeof(double)*XYZ);
  memcpy(xyzmax,other.xyzmax,sizeof(double)*XYZ);
  memcpy(dxyz,other.dxyz,sizeof(double)*XYZ);
  for(int i=0; i<XYZ; i++)
    memcpy(thegrid[i],other.thegrid[i],sizeof(double)*npoints[i]);
}


XYZGrid::XYZGrid(const char* xmlFileName)
{
  std::ifstream xml_file(xmlFileName); 

  xml_node node_labgrid("lab_xyz_grid",xml_file);
  xml_node node_axis(node_labgrid,"axis",0);

  Set();

  int npts[XYZ];
  for (int i=0; i<XYZ; i++)
    {
      npts[i]=node_axis.read_int_value("n_points");
      //Read and convert to a.u.:
      xyzmin[i]=node_axis.read_double_value("min") * ANGS_TO_BOHR;
      xyzmax[i]=node_axis.read_double_value("max") * ANGS_TO_BOHR;
    }
  Alloc(npts);

  CalcDxDyDz();
  CalcXYZGrid();
  
  xml_file.close();
}

XYZGrid::~XYZGrid() 
{ 
  Free();
}

//!Calc dx,dy,dz for the x,y,z grid.
void XYZGrid::CalcDxDyDz()
{
  dxyz[X] = (xyzmax[X]-xyzmin[X])/(npoints[X]-1);
  dxyz[Y] = (xyzmax[Y]-xyzmin[Y])/(npoints[Y]-1);
  dxyz[Z] = (xyzmax[Z]-xyzmin[Z])/(npoints[Z]-1);
}

//!Generate x,y and z values for *thegrid[XYZ].
void XYZGrid::CalcXYZGrid()
{
  for (int i=0; i<npoints[X]; i++)
      thegrid[X][i] = xyzmin[X]+i*dxyz[X];
  for (int j=0; j<npoints[Y]; j++)
      thegrid[Y][j] = xyzmin[Y]+j*dxyz[Y];
  for (int k=0; k<npoints[Z]; k++)
      thegrid[Z][k] = xyzmin[Z]+k*dxyz[Z];
}

//! Shift the center of the XYZGrid by xyzcenter[i],[j],[k].
void XYZGrid::ShiftXYZGrid(double *xyzcenter)
{
  for (int i=0; i<npoints[X]; i++)
    thegrid[X][i] -= xyzcenter[X];
  for (int i=0; i<npoints[Y]; i++)
    thegrid[Y][i] -= xyzcenter[Y];
  for (int i=0; i<npoints[Z]; i++)
    thegrid[Z][i] -= xyzcenter[Z];
}

//! Print the grid info to std output.
void XYZGrid::PrintGridInfo(FILE *fout) const
{
  fprintf(fout,"\nPrinting numerical xyz information:\n");
  fprintf(fout,"x %5d pts   %10.6lf to %10.6lf a.u.\n",npoints[X],xyzmin[X],xyzmax[X]);
  fprintf(fout,"y %5d pts   %10.6lf to %10.6lf a.u.\n",npoints[Y],xyzmin[Y],xyzmax[Y]);
  fprintf(fout,"z %5d pts   %10.6lf to %10.6lf a.u.\n",npoints[Z],xyzmin[Z],xyzmax[Z]);
  fprintf(fout,"dx dy dz:   %10.6lf  %10.6lf %10.6lf a.u.\n",dxyz[X],dxyz[Y],dxyz[Z]);
  fflush(stdout);
}

//!Print to file the grid as thegrid[XYZ][nxyz].
void XYZGrid::Print(const char *filename) const
{
  FILE *outfile = fopen(filename,"w");
  fprintf(outfile,"\n\n         i   X               Y               Z  (BOHR)\n");

  for (int n=0,i=0; i<npoints[X]; i++)
     for (int j=0; j<npoints[Y]; j++)
        for (int k=0; k<npoints[Z]; k++)
          {
            fprintf(outfile,"%10d  %+13lE  %+13lE  %+13lE\n",n,thegrid[X][i],thegrid[Y][j],thegrid[Z][k]);
            fflush(outfile);
            n++;
          }
  fclose(outfile);
}


void XYZGrid::GetPoint(int absaddr, double &pt1, double &pt2, double &pt3, Carts axis1, Carts axis2, Carts axis3) const
{
  //Absolute address on the grid is addr=nz + ny*nzmax + nx*nymax*nzmax
  int ia, ka, ja;
  int n0 = npoints[axis1];
  int n1 = npoints[axis2];
  int n2 = npoints[axis3];
  int n12 = n1*n2;
  //!Calc the address of x, y and z.
  ia = absaddr/n12;
  ja = (absaddr-ia*n12)/n2;
  ka = absaddr-ia*n12-ja*n2;
  pt1 = GetPoint(axis1,ia);
  pt2 = GetPoint(axis2,ja);
  pt3 = GetPoint(axis3,ka);
}


