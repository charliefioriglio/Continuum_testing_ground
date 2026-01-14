#include "cklm.h"
#include "sph.h"
#include "math.h"
#include "omp.h"

void CklmCoeff::Alloc(int nc, int lmx, int nk)
{
  check(!(nc>0),"CklmCoeff::Alloc() : Invalid ncklm\n");

  Free();

  double gsize=0.0;

  lmax = lmx;
  ncklm=nc;
  nkv = nk;

  cklmsq = new double*[XYZ];
  gsize+=sizeof(double *)*XYZ;
  check(!cklmsq, "CklmCoeff::Alloc() : Allocation failed");
  for (int v=0; v<XYZ; v++)
    {
      cklmsq[v] = new double[ncklm];
      gsize+=sizeof(double)*ncklm;
      check(!cklmsq[v], "CklmCoeff::Alloc() : Allocation failed");
      memset(cklmsq[v],0,sizeof(double)*ncklm);
    }

  crosscklm = new double **[XYZ];
  gsize+=sizeof(double **)*XYZ;
  check(!crosscklm, "CklmCoeff::Alloc() : Allocation failed");
  for (int v=0; v<XYZ; v++)
    {
      crosscklm[v] = new double *[ncklm];
      gsize+=sizeof(double *)*ncklm;
      check(!crosscklm[v], "CklmCoeff::Alloc() : Allocation failed");
      for (int i=0; i<ncklm; i++)
        {
          crosscklm[v][i] = new double[lmax+1];
          gsize+=sizeof(double)*(lmax+1);
          check(!crosscklm[v][i], "CklmCoeff::Alloc() : Allocation failed");
          memset(crosscklm[v][i],0,sizeof(double)*(lmax+1));
        }
     }
  sumcklm = new double*[XYZ];
  gsize+=sizeof(double)*XYZ;
  check(!sumcklm, "CklmCoeff::Alloc() : Allocation failed");
  for (int v=0; v<XYZ; v++)
    {
      sumcklm[v] = new double[nkv];
      gsize+=sizeof(double)*nkv;
      check(!sumcklm[v], "CklmCoeff::Alloc() : Allocation failed");
      memset(sumcklm[v],0,sizeof(double)*nkv);
    }

  printf("\nCklmCoeff::Alloc:  %lf MB allocated\n",gsize/MEGABITE);
}

void CklmCoeff::Free()
{
  if (IfAlloc())
   {
     for (int v=0; v<XYZ; v++)
         if (NULL!=cklmsq[v])
        {
          delete [] cklmsq[v];
          cklmsq[v]=NULL;
        }
      if (NULL!=cklmsq)
      {
        delete [] cklmsq;
        cklmsq=NULL;
      }
     for (int v=0; v<XYZ; v++)
       {
         for (int i=0; i<ncklm; i++)
             if (NULL!=crosscklm[v][i])
            {
              delete [] crosscklm[v][i];
              crosscklm[v][i]=NULL;
            }
          if (NULL!=crosscklm[v])
         {
           delete [] crosscklm[v];
           crosscklm[v]=NULL;
         }
       }
       if (NULL!=crosscklm)
      {
        delete [] crosscklm;
        crosscklm=NULL;
      }
     for (int v=0; v<XYZ; v++)
         if (NULL!=sumcklm[v])
        {
          delete [] sumcklm[v];
          sumcklm[v]=NULL;
        }
      if (NULL!=sumcklm)
      {
        delete [] sumcklm;
        sumcklm=NULL;
      }
   }
  ncklm=0;
  lmax=0;
  nkv=0;
}

void CklmCoeff::Set()
{
  lmax =0;
  ncklm=0;
  nkv=0;
  cklmsq = NULL;
  crosscklm = NULL;
  sumcklm = NULL;
}


CklmCoeff::CklmCoeff()
{
  Set();
}

//! Compute Cklm averaged over all orientations/polarizations
CklmCoeff::CklmCoeff(XYZGrid& labgrid, Orbital& dysorb,
                     RadialFunction& allrkl, KLMPoints& klmgrid)
{
  Set();
  IniCklmCoeff(labgrid,dysorb,allrkl,klmgrid);
}

//! Compute Cklm for a particular molecular orientation and laser polarization
CklmCoeff::CklmCoeff(XYZGrid& labgrid, Orbital& dysorb,
                     double *RIonz, RadialFunction& allrkl, KLMPoints& klmgrid)
{
  Set();
  IniCklmCoeff(labgrid,dysorb,RIonz,allrkl,klmgrid);
}

CklmCoeff::~CklmCoeff()
{
  Free();
}

//! Compute Cklm averaged over all orientations/polarizations
void CklmCoeff::IniCklmCoeff(XYZGrid& labgrid, Orbital& dysorb,
                          RadialFunction& allrkl, KLMPoints& klmgrid)
{
  int lmx = klmgrid.LMax();
  int nk = klmgrid.NKV();
  int nc = nk*klmgrid.NLMV();
  Alloc(nc,lmx,nk);
  
  //Print Dyson and Rkl scans 
  PrintDysonScanZ(labgrid,dysorb,allrkl,klmgrid);
  
  CalcCklmSq(labgrid,dysorb,allrkl,klmgrid);

  PrintAvgCklmSq(klmgrid);
}

//! Compute Cklm for a particular molecular orientation and laser polarization
void CklmCoeff::IniCklmCoeff(XYZGrid& labgrid, Orbital& dysorb,
                             double *RIonz, RadialFunction& allrkl,
                             KLMPoints& klmgrid)
{
  int lmx = klmgrid.LMax();
  int nk = klmgrid.NKV();
  int nc = nk*klmgrid.NLMV();
  Alloc(nc,lmx,nk);

  //Print Dyson and Rkl scans 
  PrintDysonScanZ(labgrid,dysorb,allrkl,klmgrid);

  CalcCklmSq(labgrid,dysorb,RIonz,allrkl,klmgrid);

  PrintCklmSq(klmgrid);
}

/*! Compute cklm on numerical grid. 
  Here we compute 3 components (x,y,z) to later obtain a final averaged Cklm. */
void CklmCoeff::CalcCklmSq(XYZGrid& labgrid, Orbital& dysorb,
                           RadialFunction& allrkl, KLMPoints& klmgrid)
{
  extern SPH& theSPH();
  static SPH *sph=&(theSPH());

  int nlmv = klmgrid.NLMV();
  int nxyz = labgrid.NXYZ();
  int nptx = labgrid.NPoints(X);
  int npty = labgrid.NPoints(Y);
  int nptz = labgrid.NPoints(Z);
  double *gridptr_x=labgrid.GetGridPtr(X);
  double *gridptr_y=labgrid.GetGridPtr(Y);
  double *gridptr_z=labgrid.GetGridPtr(Z);
  double dV = labgrid.DXYZ(X)*labgrid.DXYZ(Y)*labgrid.DXYZ(Z);

  Complex  **Lcklm[XYZ], **Rcklm[XYZ];

  int nthreads=1;
  nthreads=omp_get_max_threads();
  printf("Executing CalcCklmSq on %d threads\n",nthreads); fflush(stdout);

  for (int v=0; v<XYZ; v++)
    {
      Lcklm[v] = new Complex *[ncklm];
      Rcklm[v] = new Complex *[ncklm];
      for(int i=0; i<ncklm; i++)
	{
	  Lcklm[v][i]=new Complex[nthreads];
	  Rcklm[v][i]=new Complex[nthreads];
	  memset(Rcklm[v][i],0,nthreads*sizeof(Complex));
	  memset(Lcklm[v][i],0,nthreads*sizeof(Complex));
	}  
    }

  //printf("Alloc done\n"); fflush(stdout);

#pragma omp parallel for
  for(int nx=0; nx<nptx; nx++)
    {
      int thread_id=omp_get_thread_num();
      for(int ny=0; ny<npty; ny++)
	{
	  for(int nz=0; nz<nptz; nz++) 
	    {
              double Ldys_valueX, Rdys_valueX, Ldys_valueY, Rdys_valueY, Ldys_valueZ, Rdys_valueZ;
              unsigned long long int gridaddr = nx*npty*nptz + ny*nptz + nz;
              dysorb.CalcDyson_xyz(Ldys_valueX,Rdys_valueX,gridptr_x[nz],gridptr_y[ny],gridptr_z[nx]);
              dysorb.CalcDyson_xyz(Ldys_valueY,Rdys_valueY,gridptr_x[nx],gridptr_y[nz],gridptr_z[ny]);
              dysorb.CalcDyson_xyz(Ldys_valueZ,Rdys_valueZ,gridptr_x[nx],gridptr_y[ny],gridptr_z[nz]);

              Complex tLDys_cart[XYZ], tRDys_cart[XYZ], totalLDys_cart[XYZ], totalRDys_cart[XYZ];
              tLDys_cart[X]=Ldys_valueX*gridptr_z[nz];
              tLDys_cart[Y]=Ldys_valueY*gridptr_z[nz];
              tLDys_cart[Z]=Ldys_valueZ*gridptr_z[nz];
              tRDys_cart[X]=Rdys_valueX*gridptr_z[nz];
              tRDys_cart[Y]=Rdys_valueY*gridptr_z[nz];
              tRDys_cart[Z]=Rdys_valueZ*gridptr_z[nz];

	      for (int k=0, no=0, kl=0; k<nkv; k++)
		for(int l=0; l<=lmax; l++,kl++)
		  {
		    double rkl(allrkl.GetRkl_xyz(kl,gridaddr));

		    for (int m=-l; m<=l; m++, no++)
		      {
			Complex ylm=sph->GetValue(gridptr_x[nx],gridptr_y[ny],gridptr_z[nz],l,m);
			Complex ylm_conj(ylm.Conj());
			//Complex tLdys_value(Ldys_value*rkl*ylm),tRdys_value(Rdys_value*rkl*ylm.Conj());
			//i^l factor is included in SphHarm
			for (int v=0; v<XYZ; v++)
			  {
			    totalLDys_cart[v]=tLDys_cart[v]*ylm*rkl;
			    totalRDys_cart[v]=tRDys_cart[v]*ylm_conj*rkl;
			    Lcklm[v][no][thread_id] += totalLDys_cart[v];
			    Rcklm[v][no][thread_id] += totalRDys_cart[v];
			  }
		      }
		  }
	    }
	}
    }
  
  //printf("Loop done\n"); fflush(stdout);

  for (int k=0,no=0; k<nkv; k++)
    for (int lm=0; lm<nlmv; lm++,no++)
      for (int v=0; v<XYZ; v++)
	for(int i=1; i<nthreads; i++)
	  {
	    //Add everything to the Lcklm[v][no][0]
	    Lcklm[v][no][0] += Lcklm[v][no][i];
	    Rcklm[v][no][0] += Rcklm[v][no][i];
	  }

  double dVV=dV*dV;

  for (int k=0,no=0; k<nkv; k++)
    for (int lm=0; lm<nlmv; lm++,no++)
      for (int v=0; v<XYZ; v++)
	{
	  
	  double tmp=(Lcklm[v][no][0]*Rcklm[v][no][0]).Re();
	  tmp*=dVV;
	  cklmsq[v][no] += tmp;
	}
  
  for (int k=0,no=0; k<nkv; k++)
    {
      int klm = k*(lmax+1)*(lmax+1);
      for (int l=0; l<=lmax; l++)
        for (int m=-l; m<=l; m++, no++)
          for (int l2=abs(m); l2<=lmax; l2++)
            {
              int no2 = klm + l2*(l2+1) + m;
              for (int v=0; v<XYZ; v++)
		{
		  double tmp=(Lcklm[v][no][0]*Rcklm[v][no2][0]).Re();
		  tmp*=dVV;
		  crosscklm[v][no][l2] += tmp;
		  //crosscklm[v][no][l2] += (Lcklm[v][no]*Rcklm[v][no2]).Re();
		}
            }
    }
  
  for (int k=0,no=0; k<nkv; k++)
    for (int l=0; l<=lmax; l++)
      for (int m=-l; m<=l; m++, no++)
        for (int v=0; v<XYZ; v++)
	  sumcklm[v][k] += cklmsq[v][no];
  
  //printf("Dealloc\n"); fflush(stdout);
  for (int v=0; v<XYZ; v++)
    {
      for(int i=0; i<ncklm; i++)
	{
	  delete [] Rcklm[v][i];
	  delete [] Lcklm[v][i];
	}
    }
  //printf("Dealloc-done\n"); fflush(stdout);
}

//! Compute Cklm for a particular molecular orientation and laser polarization
void CklmCoeff::CalcCklmSq(XYZGrid& labgrid, Orbital& dysorb,
                           double *RIonz, RadialFunction& allrkl, KLMPoints& klmgrid)
{
  extern SPH& theSPH() ;
  static SPH *sph=&(theSPH());

  int nlmv = klmgrid.NLMV();
  int nxyz = labgrid.NXYZ();
  int nptx = labgrid.NPoints(X);
  int npty = labgrid.NPoints(Y);
  int nptz = labgrid.NPoints(Z);
  double *gridptr_x=labgrid.GetGridPtr(X);
  double *gridptr_y=labgrid.GetGridPtr(Y);
  double *gridptr_z=labgrid.GetGridPtr(Z);
  double dV = labgrid.DXYZ(X)*labgrid.DXYZ(Y)*labgrid.DXYZ(Z);
  double xyzL[XYZ], Ldys_value, Rdys_value;
  double unitvect[XYZ];

  Complex  **Lcklm[XYZ], **Rcklm[XYZ];

   int nthreads=1;
   nthreads=omp_get_max_threads();
   printf("Executing CalcCklmSq on %d threads\n",nthreads); fflush(stdout);

  for (int v=0; v<XYZ; v++)
    {
       Lcklm[v] = new Complex *[ncklm];
       Rcklm[v] = new Complex *[ncklm];
       for(int i=0; i<ncklm; i++)
       {
         Lcklm[v][i]=new Complex[nthreads];
         Rcklm[v][i]=new Complex[nthreads];
         memset(Rcklm[v][i],0,nthreads*sizeof(Complex));
         memset(Lcklm[v][i],0,nthreads*sizeof(Complex));
       }  
     }

  unitvect[X]=RIonz[X]/sqrt(pow(RIonz[X],2)+pow(RIonz[Y],2)+pow(RIonz[Z],2));
  unitvect[Y]=RIonz[Y]/sqrt(pow(RIonz[X],2)+pow(RIonz[Y],2)+pow(RIonz[Z],2));
  unitvect[Z]=RIonz[Z]/sqrt(pow(RIonz[X],2)+pow(RIonz[Y],2)+pow(RIonz[Z],2));

// SG: Maybe add an "if unitvec not equal to 0" somewhere below? could save 1/3 or 2/3 of the time required.

#pragma omp parallel for
  for(int nx=0; nx<nptx; nx++)
    {
      int thread_id=omp_get_thread_num();
      for(int ny=0; ny<npty; ny++)
        {
          for(int nz=0; nz<nptz; nz++)
            {
              double Ldys_value, Rdys_value;
              unsigned long long int gridaddr = nx*npty*nptz + ny*nptz + nz;
              dysorb.CalcDyson_xyz(Ldys_value,Rdys_value,gridptr_x[nx],gridptr_y[ny],gridptr_z[nz]);

              Complex tLDys_cart[XYZ], tRDys_cart[XYZ], totalLDys_cart[XYZ], totalRDys_cart[XYZ];
              tLDys_cart[X]=Ldys_value*gridptr_x[nx];
              tLDys_cart[Y]=Ldys_value*gridptr_y[ny];
              tLDys_cart[Z]=Ldys_value*gridptr_z[nz];
              tRDys_cart[X]=Rdys_value*gridptr_x[nx];
              tRDys_cart[Y]=Rdys_value*gridptr_y[ny];
              tRDys_cart[Z]=Rdys_value*gridptr_z[nz];

              for (int k=0, no=0, kl=0; k<nkv; k++)
                for(int l=0; l<=lmax; l++,kl++)
                  {
                    double rkl(allrkl.GetRkl_xyz(kl,gridaddr));

                    for (int m=-l; m<=l; m++, no++)
                      {
                        Complex ylm=sph->GetValue(gridptr_x[nx],gridptr_y[ny],gridptr_z[nz],l,m);
                        Complex ylm_conj(ylm.Conj());
                        //i^l factor is included in SphHarm

                        for (int v=0; v<XYZ; v++)
                          {
                            totalLDys_cart[v]=tLDys_cart[v]*ylm*rkl*unitvect[v];
                            totalRDys_cart[v]=tRDys_cart[v]*ylm_conj*rkl*unitvect[v];
                            Lcklm[v][no][thread_id] += totalLDys_cart[v];
                            Rcklm[v][no][thread_id] += totalRDys_cart[v];
                          }
                      }
                  }
            }
        }
    }

   for (int k=0,no=0; k<nkv; k++)
     for (int lm=0; lm<nlmv; lm++,no++)
       for (int v=0; v<XYZ; v++)
       for(int i=1; i<nthreads; i++)
         {
           //Add everything to the Lcklm[v][no][0]
           Lcklm[v][no][0] += Lcklm[v][no][i];
           Rcklm[v][no][0] += Rcklm[v][no][i];
         }

    double dVV=dV*dV;

  for (int k=0,no=0; k<nkv; k++)
    for (int lm=0; lm<nlmv; lm++,no++)
      for (int v=0; v<XYZ; v++)
        {
          double tmp=(Lcklm[v][no][0]*Rcklm[v][no][0]).Re();
          tmp*=dVV;
          cklmsq[v][no] += tmp;
        }

  for (int k=0,no=0; k<nkv; k++)
    {
      int klm = k*(lmax+1)*(lmax+1);
      for (int l=0; l<=lmax; l++)
        for (int m=-l; m<=l; m++, no++)
          for (int l2=abs(m); l2<=lmax; l2++)
            {

              int no2 = klm + l2*(l2+1) + m;
              for (int v=0; v<XYZ; v++)
               {
                 double tmp=(Lcklm[v][no][0]*Rcklm[v][no2][0]).Re();
                 tmp*=dVV;
                 crosscklm[v][no][l2] += tmp;
                 //crosscklm[v][no][l2] += (Lcklm[v][no]*Rcklm[v][no2]).Re();
               }
            }
    }

  for (int k=0,no=0; k<nkv; k++)
    for (int l=0; l<=lmax; l++)
      for (int m=-l; m<=l; m++, no++)
        for (int v=0; v<XYZ; v++)
          sumcklm[v][k] += cklmsq[v][no];

  for (int v=0; v<XYZ; v++)
    {
       for(int i=0; i<ncklm; i++)
       {
         delete [] Rcklm[v][i];
         delete [] Lcklm[v][i];
       }
    }
}



//! Print all square Cklm coeffs, ||Cklm^2||, for averaged molecular orientations.
void CklmCoeff::PrintAvgCklmSq(KLMPoints& klmgrid, FILE *fout)
{
  fprintf(fout,"\n||Cklm^2|| printout\n");
  fprintf(fout,"k,a.u.         l      m      Re(Cklm[x]^2)     Re(Cklm[y]^2)     Re(Cklm[z]^2)     Re(|Cklm|^2)\n");

  int nlmv = klmgrid.NLMV();
  int no;
  double eng, kau;
  for(int k=0,no=0; k<nkv; k++)
    for (int lm=0; lm<nlmv; lm++,no++)
      {
      fprintf(fout,"%.3lE  %5d  %5d",klmgrid.GetKV(k),klmgrid.GetLV(lm),klmgrid.GetMV(lm));
      fprintf(fout," %17.6lE %17.6lE %17.6lE %17.6lE\n",cklmsq[0][no],cklmsq[1][no],cklmsq[2][no],
               (cklmsq[0][no]+cklmsq[1][no]+cklmsq[2][no])/3);
      }
  fprintf(fout,"\n  E_k,eV     Sum(|Cklm[x]|^2)  Sum(|Cklm[y]|^2)  Sum(|Cklm[z]|^2)  Sum(|(Cklm)|^2)\n");
  for (int k=0; k<nkv; k++)
    {
      kau = klmgrid.GetKV(k);
      eng = 0.5*kau*kau*HAR_TO_EV;
      fprintf(fout,"%7.3lf",eng);
      fprintf(fout," %17.6lE %17.6lE %17.6lE %17.6lE\n",sumcklm[0][k],sumcklm[1][k],sumcklm[2][k],
             (sumcklm[0][k]+sumcklm[1][k]+sumcklm[2][k])/3);
    }
}

//! Print all square Cklm coeffs, ||Cklm^2||.
void CklmCoeff::PrintCklmSq(KLMPoints& klmgrid, FILE *fout)
{
  fprintf(fout,"\n||Cklm^2|| printout\n");
  fprintf(fout,"k,a.u.         l      m      Re(Cklm[x]^2)     Re(Cklm[y]^2)     Re(Cklm[z]^2)     Re(|Cklm|^2)\n");

  int nlmv = klmgrid.NLMV();
  int no;
  double eng, kau;
  for(int k=0,no=0; k<nkv; k++)
    for (int lm=0; lm<nlmv; lm++,no++)
      {
      fprintf(fout,"%.3lE  %5d  %5d",klmgrid.GetKV(k),klmgrid.GetLV(lm),klmgrid.GetMV(lm));
      fprintf(fout," %17.6lE %17.6lE %17.6lE %17.6lE\n",cklmsq[0][no],cklmsq[1][no],cklmsq[2][no],
               (cklmsq[0][no]+cklmsq[1][no]+cklmsq[2][no]));
      }
  fprintf(fout,"\n  E_k,eV     Sum(|Cklm[x]|^2)  Sum(|Cklm[y]|^2)  Sum(|Cklm[z]|^2)  Sum(|(Cklm)|^2)\n");
  for (int k=0; k<nkv; k++)
    {
      kau = klmgrid.GetKV(k);
      eng = 0.5*kau*kau*HAR_TO_EV;
      fprintf(fout,"%7.3lf",eng);
      fprintf(fout," %17.6lE %17.6lE %17.6lE %17.6lE\n",sumcklm[0][k],sumcklm[1][k],sumcklm[2][k],
             (sumcklm[0][k]+sumcklm[1][k]+sumcklm[2][k]));
    }

//SG: disabled printing of cross-terms
  
//  fprintf(fout,"\nAverage Cross-terms\n");
//  fprintf(fout,"   L     M     L' Re(CklmCkl'm),a.u.\n");
//  for (int k=0,no=0; k<nkv; k++)
//    for (int l=0; l<=lmax; l++)
//      for (int m=-l; m<=l; m++, no++)
//       {
//         fprintf(fout,"%4d  %+4d  ",l,m);
//         for (int l2=0; l2<=lmax; l2++)
//           {
//           fprintf(fout,"%4d",l2);
//           for (int v=0; v<XYZ; v++)
//             fprintf(fout,"  %+.3lE",crosscklm[v][no][l2]);
//           }
//         fprintf(fout,"\n");
//       }
//  fflush(fout);
}


void CklmCoeff::PrintDysonScanZ(XYZGrid& labgrid, Orbital& dysorb, RadialFunction& allrkl, KLMPoints& klmgrid)
{
  FILE *file_orbsR=fopen("orbitalR.dat","w");
  FILE *file_orbsL=fopen("orbitalL.dat","w");
//  FILE *file_rkl=fopen("rkl.dat","w");
  double *gridptr_x=labgrid.GetGridPtr(X);
  double *gridptr_y=labgrid.GetGridPtr(Y);
  double *gridptr_z=labgrid.GetGridPtr(Z);
  double xyz[XYZ], Ldys_value, Rdys_value, rkl;
//  int midpointx=labgrid.NPoints(X)/2;
//  int midpointy=labgrid.NPoints(Y)/2;
//  int midpointz=labgrid.NPoints(Z)/2;
  
  fprintf(file_orbsR, "X,A RDys_val RD*X Y,A RDys_val RD*Y Z,A RDys_val RD*Z \n");
  fprintf(file_orbsL, "X,A LDys_val LD*X Y,A LDys_val LD*Y Z,A LDys_val LD*Z \n");
 
  for(int i=0; i<labgrid.NPoints(Z); i++)
    {
      xyz[X]=gridptr_x[i];
      xyz[Y]=gridptr_y[i];
      xyz[Z]=gridptr_z[i];
      dysorb.CalcDyson_xyz(Ldys_value,Rdys_value,xyz[X],0.0,0.0);
      fprintf(file_orbsR,"%lE %lE %lE ",xyz[X]*BOHR_TO_ANGS, Rdys_value, Rdys_value*xyz[X]);
      fprintf(file_orbsL,"%lE %lE %lE ",xyz[X]*BOHR_TO_ANGS, Ldys_value, Ldys_value*xyz[X]);

      dysorb.CalcDyson_xyz(Ldys_value,Rdys_value,0.0,xyz[Y],0.0);
      fprintf(file_orbsR,"%lE %lE %lE ",xyz[Y]*BOHR_TO_ANGS, Rdys_value, Rdys_value*xyz[Y]);
      fprintf(file_orbsL,"%lE %lE %lE ",xyz[Y]*BOHR_TO_ANGS, Ldys_value, Ldys_value*xyz[Y]);

      dysorb.CalcDyson_xyz(Ldys_value,Rdys_value,0.0,0.0,xyz[Z]);
      fprintf(file_orbsR,"%lE %lE %lE\n",xyz[Z]*BOHR_TO_ANGS, Rdys_value, Rdys_value*xyz[Z]);
      fprintf(file_orbsL,"%lE %lE %lE\n",xyz[Z]*BOHR_TO_ANGS, Ldys_value, Ldys_value*xyz[Z]);
    }   
  
  fclose(file_orbsR);
  fclose(file_orbsL);

}
