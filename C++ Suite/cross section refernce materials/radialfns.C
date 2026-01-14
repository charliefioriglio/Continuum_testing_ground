#include "radialfns.h"

void RadialFunction::Alloc(int nk, int lmx, int nxyzpts)
{
  check(!(nk*nxyzpts>0) || !(lmx>=0),"RadialFunction::Alloc() : Invalid nkv,nxyz or lmax\n");

  Free();
  double gsize=0.0;
  
  nkv=nk;
  lmax=lmx;
  nrkl=nk*(lmx+1);
  nxyz=nxyzpts;
  rkl = new double*[nrkl];
  gsize+=sizeof(double)*nrkl;

  
  check(!(rkl>NULL), "RadialFunction::Alloc() : Allocation failed");

  for (int i=0; i<nrkl; i++)
    {
      rkl[i] = new double[nxyz];
      gsize+=sizeof(double)*nxyz;
      check(!(rkl[i]>NULL), "RadialFunction::Alloc() : Allocation failed");
    }

  printf("\nRadialFunction::Alloc:  %lf MB allocated\n",gsize/MEGABITE);

}

void RadialFunction::Free()
{
  if (IfAlloc())
   {
     for (int i=0; i<nrkl; i++)
         if (NULL!=rkl[i])
          {
            delete [] rkl[i];
            rkl[i]=NULL;
          }
     if (NULL!=rkl)
      {
        delete [] rkl;
        rkl=NULL;
      }

     nkv=0;
     nrkl=0;
     nxyz=0;
   }
}

void RadialFunction::Set()
{
  nkv=0;
  nrkl=0;
  nxyz=0;
  for (int i=0; i<nrkl; i++)
      rkl[i] = NULL;
  rkl = NULL;

}

RadialFunction::RadialFunction()
{
  Set();
}

RadialFunction::RadialFunction(const char* xmlFileName, const XYZGrid& labgrid, const KLMPoints& klmgrid)
{
  Set();
  IniRadialFunction(xmlFileName,labgrid,klmgrid);
}

RadialFunction::~RadialFunction() 
{ 
  Free();
}

void RadialFunction::IniRadialFunction(const char* xmlFileName, const XYZGrid& labgrid, const KLMPoints& klmgrid)
{
  int nk = klmgrid.NKV();
  int lmx = klmgrid.LMax();
  int nxyzpts = labgrid.NXYZ();
  Alloc(nk,lmx,nxyzpts);

  simpleXMLparser xmlF;
  xmlF.assignFile(xmlFileName); 
  xmlF.reset().node("root").node("free_electron");
  std::string tmpStr;
  tmpStr=xmlF.value("radial_functions");
  if (tmpStr=="spherical")
    radfn=SphR;
  else if (tmpStr=="coulomb")
    radfn=CoulR;
  else if (tmpStr=="coulombzerok")
    radfn=CoulR0l;
  else if (tmpStr=="library")
    radfn=LibR;
  else
    xmlF.exitOnFormatError(true);

  std::cout << "Radial function = " << radfn <<'\n';

  CalcAllRkl(labgrid,klmgrid);
  //CalcRklOverlaps(klmgrid);
}

//! Calc all Rkl values (nrkl*nxyz) at x,y,z grid points.
void RadialFunction::CalcAllRkl(const XYZGrid& labgrid, const KLMPoints& klmgrid)
{
  double xi, yi, zi;
  int lmax = klmgrid.LMax();
  double *rklmax = new double[lmax+1];  
  check(!(rklmax>NULL), "RadialFunction::CalcAllRkl() : Allocation failed");

  for (int i=0; i<nxyz; i++)
    {
      labgrid.GetPoint(i,xi,yi,zi,X,Y,Z);
      for (int k=0,no=0; k<nkv; k++)
        {
          Rkl(radfn,rklmax,xi,yi,zi,klmgrid.GetKV(k),lmax);
          for (int l=0; l<=lmax; l++,no++)
	    //rkl[k*(lmax+1)+l][i] = rklmax[l];
	    rkl[no][i] = rklmax[l];
        }
    }
  delete [] rklmax;
}

//! Print Rkl matrix, for all k,l and x,y,z values.
void RadialFunction::Print(const char *filename) const
{
  FILE *outfile = fopen(filename,"w");
  fprintf(outfile,"\nRkl values\n");
  fprintf(outfile,"         i   ");
  for (int k=0; k<nkv; k++)
      for (int l=0; l<=lmax; l++)
          fprintf(outfile,"R%d_%d           ",k,l);
  fprintf(outfile,"\n");

  for (int i=0; i<nxyz; i++)
    {
      fprintf(outfile,"%10d",i);
      for (int k=0; k<nkv; k++)
          for (int l=0; l<=lmax; l++)
              fprintf(outfile,"  %+13lE",rkl[k*(lmax+1)+l][i]);
      fprintf(outfile,"\n");
    }
  fclose(outfile);
}  

//! Calc norm matrix RklRkl'
void RadialFunction::NormMatrix(double rmax, int nrv, const KLMPoints& klmgrid)
{
  double dr = rmax*ANGS_TO_BOHR/(nrv-1);
  double ri;
  int kl, kl1;
  
  double *rklmax = new double[lmax+1];
  double **rkl_r = new double*[nrkl];
  for (int i=0; i<nrkl; i++)
     rkl_r[i] = new double[nrv];
  for (int i=0; i<nrkl; i++)
   for (int j=0; j<nrv; j++)
    rkl_r[i][j] = 0.0;

  double **rklnorm = new double*[nrkl];
  for (int i=0; i<nrkl; i++)
     rklnorm[i] = new double[lmax+1];
  for (int i=0; i<nrkl; i++)
   for (int j=0; j<=lmax; j++)
      rklnorm[i][j] = 0.0;

  ri = 0.0;
  for (int i=0; i<nrv; i++)
   { 
    for (int k=0,no=0; k<nkv; k++)
      {
        Rkl(radfn,rklmax,ri,0.0,0.0,klmgrid.GetKV(k),lmax);
        for (int l=0; l<=lmax; l++,no++)
           rkl_r[no][i] = rklmax[l];
      }
     ri += dr;

   }

  delete [] rklmax;

  ri = 0.0;
  for (int i=0; i<nrv; i++)
   {
     for (int k=0; k<nkv; k++)
      for (int l=0; l<=lmax; l++)
       for (int l1=0; l1<=lmax; l1++)
        {
          kl = k*(lmax+1) + l;
          kl1 = k*(lmax+1) +l1;
          rklnorm[kl][l1] += rkl_r[kl][i]*rkl_r[kl1][i]*ri*ri*dr;
        }
     ri += dr;
   }

  for (int i=0; i<nrkl; i++)
     delete [] rkl_r[i];
  delete [] rkl_r;

  for (int k=0; k<nkv; k++)
   {
     fprintf(stderr,"k=%d\n",k);
     for (int l=0; l<=lmax; l++)
      {
        fprintf(stderr," l=%d",l);
        kl = k*(lmax+1) + l;
        for (int l1=0; l1<=lmax; l1++)
           fprintf(stderr,"  %lE",rklnorm[kl][l1]);
        fprintf(stderr,"\n");
      }
   }

  for (int i=0; i<nrkl; i++)
     delete [] rklnorm[i];
  delete [] rklnorm;
}
