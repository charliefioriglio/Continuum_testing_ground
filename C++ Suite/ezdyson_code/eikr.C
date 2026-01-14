#include "eikr.h"
#include "math.h"
#include "omp.h"


void NumEikr::Alloc(int nk)
{
  check(!(nk>0),"NumEikr::Alloc() : Invalid nkv\n");

  Free();

  double gsize=0.0;

  nkv = nk;

  cpar = new double[nkv];
  gsize+=sizeof(double)*nkv;
  check(!cpar, "NumEikr::Alloc() : Allocation failed");
  memset(cpar,0,sizeof(double)*nkv);

  cperp = new double[nkv];
  gsize+=sizeof(double)*nkv;
  check(!cperp, "NumEikr::Alloc() : Allocation failed");
  memset(cperp,0,sizeof(double)*nkv);

  printf("\nNumEikr::Alloc:  %lf MB allocated\n",gsize/MEGABITE);
}

void NumEikr::Free()
{
  if (IfAlloc())
   {
      if (NULL!=cpar)
      {
        delete [] cpar;
        cpar=NULL;
      }
      if (NULL!=cperp)
      {
        delete [] cperp;
        cperp=NULL;
      }
   }
  nkv=0;
}

void NumEikr::Set()
{
  nkv=0;
  cpar = NULL;
  cperp = NULL;
}


NumEikr::NumEikr()
{
  Set();
}

//! Compute cross sections with numerical REPULSION averaging
NumEikr::NumEikr(XYZGrid& labgrid, Orbital& dysorb,
                     AngleGrid& anggrid, KLMPoints& klmgrid)
{
  Set();
  IniNumEikr(labgrid,dysorb,anggrid,klmgrid);
}

NumEikr::~NumEikr()
{
  Free();
}

//! Compute cross sections with numerical REPULSION averaging
void NumEikr::IniNumEikr(XYZGrid& labgrid, Orbital& dysorb,
                             AngleGrid& anggrid, KLMPoints& klmgrid)
{
  int nk = klmgrid.NKV();
  Alloc(nk);

  PrintDysonScanZ(labgrid,dysorb,klmgrid);

  CalcEikrSq(labgrid,dysorb,anggrid,klmgrid);
}

//! Compute cross sections with numerical REPULSION averaging
void NumEikr::CalcEikrSq(XYZGrid& labgrid, Orbital& dysorb,
                           AngleGrid& anggrid, KLMPoints& klmgrid)
{
  int nxyz = labgrid.NXYZ();
  int nptx = labgrid.NPoints(X);
  int npty = labgrid.NPoints(Y);
  int nptz = labgrid.NPoints(Z);
  double *gridptr_x=labgrid.GetGridPtr(X);
  double *gridptr_y=labgrid.GetGridPtr(Y);
  double *gridptr_z=labgrid.GetGridPtr(Z);
  double dV = labgrid.DXYZ(X)*labgrid.DXYZ(Y)*labgrid.DXYZ(Z);
  double Ldys_value, Rdys_value;

  int nabc=anggrid.NABC();

  Complex  **Lcklm_par[nabc], **Rcklm_par[nabc];
  Complex  **Lcklm_x[nabc], **Rcklm_x[nabc], **Lcklm_y[nabc], **Rcklm_y[nabc];

  int nthreads=1;
  nthreads=omp_get_max_threads();
  printf("Executing CalcEikrSq on %d threads\n",nthreads); fflush(stdout);

  for (int v=0; v<nabc; v++)
    {
      Lcklm_par[v] = new Complex *[nkv];
      Rcklm_par[v] = new Complex *[nkv];
      Lcklm_x[v] = new Complex *[nkv];
      Rcklm_x[v] = new Complex *[nkv];
      Lcklm_y[v] = new Complex *[nkv];
      Rcklm_y[v] = new Complex *[nkv];
      for(int i=0; i<nkv; i++)
        {
          Lcklm_par[v][i]=new Complex[nthreads];
          Rcklm_par[v][i]=new Complex[nthreads];
          Lcklm_x[v][i]=new Complex[nthreads];
          Rcklm_x[v][i]=new Complex[nthreads];
          Lcklm_y[v][i]=new Complex[nthreads];
          Rcklm_y[v][i]=new Complex[nthreads];
          memset(Rcklm_par[v][i],0,nthreads*sizeof(Complex));
          memset(Lcklm_par[v][i],0,nthreads*sizeof(Complex));
          memset(Rcklm_x[v][i],0,nthreads*sizeof(Complex));
          memset(Lcklm_x[v][i],0,nthreads*sizeof(Complex));
          memset(Rcklm_y[v][i],0,nthreads*sizeof(Complex));
          memset(Lcklm_y[v][i],0,nthreads*sizeof(Complex));
        }
    }

#pragma omp parallel for
  for (int v=0; v<nabc; v++) {
    int thread_id=omp_get_thread_num();
    double aj, bj, wj, tmpavg;
    double xM, yM, zM;
    aj = anggrid.GetPoint(A,v);
    bj = anggrid.GetPoint(B,v);
    wj = anggrid.GetPoint(G,v);
    MolAvg avgcklm=NUM;
    RotnMatr rot(avgcklm,0,bj,aj);
    tmpavg = wj;
    for(int nx=0; nx<nptx; nx++)
     {
      for(int ny=0; ny<npty; ny++)
        {
          for(int nz=0; nz<nptz; nz++)
            {
              double Ldys_value, Rdys_value;
              unsigned long long int gridaddr = nx*npty*nptz + ny*nptz + nz;
              rot.GenLabMol(xM,yM,zM,gridptr_x[nx],gridptr_y[ny],gridptr_z[nz]);
              dysorb.CalcDyson_xyz(Ldys_value,Rdys_value,xM,yM,zM);
              Complex tLDys_cart, tRDys_cart, totalLDys_par, totalRDys_par;
              Complex totalLDys_x, totalRDys_x, totalLDys_y, totalRDys_y; 
              tLDys_cart=Ldys_value*gridptr_z[nz];
              tRDys_cart=Rdys_value*gridptr_z[nz];
               for (int k=0; k<nkv; k++)
                {
                  double coskz = cos(klmgrid.GetKV(k)*gridptr_z[nz]);
                  double sinkz = sin(klmgrid.GetKV(k)*gridptr_z[nz]);
                  double coskx = cos(klmgrid.GetKV(k)*gridptr_x[nx]);
                  double sinkx = sin(klmgrid.GetKV(k)*gridptr_x[nx]);
                  double cosky = cos(klmgrid.GetKV(k)*gridptr_y[ny]);
                  double sinky = sin(klmgrid.GetKV(k)*gridptr_y[ny]);
                  Complex eikrx(coskx,sinkx);
                  Complex eikry(cosky,sinky);
                  Complex eikrz(coskz,sinkz);
                  totalLDys_par=tLDys_cart*eikrz;
                  totalRDys_par=tRDys_cart*eikrz.Conj();
                  totalLDys_x=tLDys_cart*eikrx;
                  totalRDys_x=tRDys_cart*eikrx.Conj();
                  totalLDys_y=tLDys_cart*eikry;
                  totalRDys_y=tRDys_cart*eikry.Conj();
                  Lcklm_par[v][k][thread_id] += totalLDys_par;
                  Rcklm_par[v][k][thread_id] += totalRDys_par;
                  Lcklm_x[v][k][thread_id] += totalLDys_x;
                  Rcklm_x[v][k][thread_id] += totalRDys_x;
                  Lcklm_y[v][k][thread_id] += totalLDys_y;
                  Rcklm_y[v][k][thread_id] += totalRDys_y;
                }
            }
        }
    }

   for (int k=0; k<nkv; k++)
     for(int i=1; i<nthreads; i++)
       {
         Lcklm_par[v][k][0] += Lcklm_par[v][k][i];
         Rcklm_par[v][k][0] += Rcklm_par[v][k][i];
         Lcklm_x[v][k][0] += Lcklm_x[v][k][i];
         Rcklm_x[v][k][0] += Rcklm_x[v][k][i];
         Lcklm_y[v][k][0] += Lcklm_y[v][k][i];
         Rcklm_y[v][k][0] += Rcklm_y[v][k][i];
       }

    double dVV=dV*dV;

  for (int k=0; k<nkv; k++)
        {
          double tmp_par=(Lcklm_par[v][k][0]*Rcklm_par[v][k][0]).Re();
          double tmp_x=(Lcklm_x[v][k][0]*Rcklm_x[v][k][0]).Re();
          double tmp_y=(Lcklm_y[v][k][0]*Rcklm_y[v][k][0]).Re();
          tmp_par*=dVV*tmpavg;
          tmp_x*=dVV*tmpavg;
          tmp_y*=dVV*tmpavg;
          cpar[k] += tmp_par;
          cperp[k] += 0.5*tmp_x+0.5*tmp_y;
        }
  }

  for (int v=0; v<nabc; v++)
    {
       for(int i=0; i<nkv; i++)
       {
         delete [] Rcklm_par[v][i];
         delete [] Lcklm_par[v][i];
         delete [] Rcklm_x[v][i];
         delete [] Lcklm_x[v][i];
         delete [] Rcklm_y[v][i];
         delete [] Lcklm_y[v][i];
       }
    }
}

void NumEikr::PrintDysonScanZ(XYZGrid& labgrid, Orbital& dysorb, KLMPoints& klmgrid)
{
  FILE *file_orbsR=fopen("orbitalR.dat","w");
  FILE *file_orbsL=fopen("orbitalL.dat","w");
  double *gridptr_x=labgrid.GetGridPtr(X);
  double *gridptr_y=labgrid.GetGridPtr(Y);
  double *gridptr_z=labgrid.GetGridPtr(Z);
  double xyz[XYZ], Ldys_value, Rdys_value, rkl;
  
  fprintf(file_orbsR, "X(a.u.) RDys_val RD*X Y(a.u.) RDys_val RD*Y Z(a.u.) RDys_val RD*Z \n");
  fprintf(file_orbsL, "X(a.u.) LDys_val LD*X Y(a.u.) LDys_val LD*Y Z(a.u.) LDys_val LD*Z \n");
 
  for(int i=0; i<labgrid.NPoints(Z); i++)
    {
      xyz[X]=gridptr_x[i];
      xyz[Y]=gridptr_y[i];
      xyz[Z]=gridptr_z[i];
      dysorb.CalcDyson_xyz(Ldys_value,Rdys_value,xyz[X],0.0,0.0);
      fprintf(file_orbsR,"%lE %lE %lE ",xyz[X], Rdys_value, Rdys_value*xyz[X]);
      fprintf(file_orbsL,"%lE %lE %lE ",xyz[X], Ldys_value, Ldys_value*xyz[X]);

      dysorb.CalcDyson_xyz(Ldys_value,Rdys_value,0.0,xyz[Y],0.0);
      fprintf(file_orbsR,"%lE %lE %lE ",xyz[Y], Rdys_value, Rdys_value*xyz[Y]);
      fprintf(file_orbsL,"%lE %lE %lE ",xyz[Y], Ldys_value, Ldys_value*xyz[Y]);

      dysorb.CalcDyson_xyz(Ldys_value,Rdys_value,0.0,0.0,xyz[Z]);
      fprintf(file_orbsR,"%lE %lE %lE\n",xyz[Z], Rdys_value, Rdys_value*xyz[Z]);
      fprintf(file_orbsL,"%lE %lE %lE\n",xyz[Z], Ldys_value, Ldys_value*xyz[Z]);
    }   
  
  fclose(file_orbsR);
  fclose(file_orbsL);

}


