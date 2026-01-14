#include "cklm.h"
#include "sph.h"
#include "math.h"
#include "omp.h"
#include "clebsh_gordan_coeff.h"

void CklmCoeff::Alloc(int nc, int lmx, int nk)
{
  check(!(nc>0),"CklmCoeff::Alloc() : Invalid ncklm\n");

  Free();

  double gsize=0.0;

  lmax = lmx;
  ncklm=nc;
  nkv = nk;

    cklmsq = new double[ncklm];
    gsize+=sizeof(double)*ncklm;
    check(!cklmsq, "CklmCoeff::Alloc() : Allocation failed");
    memset(cklmsq,0,sizeof(double)*ncklm);

    crosscklm = new double *[ncklm];
    gsize+=sizeof(double *)*ncklm;
    check(!crosscklm, "CklmCoeff::Alloc() : Allocation failed");
    for (int i=0; i<ncklm; i++)
      {
        crosscklm[i] = new double[(lmax+1)*(lmax+1)];
        gsize+=sizeof(double)*(lmax+1)*(lmax+1);
        check(!crosscklm[i], "CklmCoeff::Alloc() : Allocation failed");
        memset(crosscklm[i],0,sizeof(double)*(lmax+1)*(lmax+1));
      }
    sumcklm = new double[nkv];
    gsize+=sizeof(double)*nkv;
    check(!sumcklm, "CklmCoeff::Alloc() : Allocation failed");
    memset(sumcklm,0,sizeof(double)*nkv);

  printf("\nCklmCoeff::Alloc:  %lf MB allocated\n",gsize/MEGABITE);
}

void CklmCoeff::Free()
{
  if (IfAlloc())
   {
      if (NULL!=cklmsq)
      {
        delete [] cklmsq;
        cklmsq=NULL;
      }
      for (int i=0; i<ncklm; i++)
          if (NULL!=crosscklm[i])
         {
           delete [] crosscklm[i];
           crosscklm[i]=NULL;
         }
       if (NULL!=crosscklm)
      {
        delete [] crosscklm;
        crosscklm=NULL;
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
                     double radfn, KLMPoints& klmgrid, bool do_overlap)
{
  Set();
  if(!do_overlap)
    IniCklmCoeff(labgrid,dysorb,radfn,klmgrid);
  else
    IniCklmCoeff_Overlap(labgrid,dysorb,radfn,klmgrid);
}

//! Compute Cklm for a particular molecular orientation and laser polarization
CklmCoeff::CklmCoeff(XYZGrid& labgrid, Orbital& dysorb,
                     double *RIonz, double radfn, KLMPoints& klmgrid)
{
  Set();
  IniCklmCoeff(labgrid,dysorb,RIonz,radfn,klmgrid);
}

CklmCoeff::~CklmCoeff()
{
  Free();
}

//! Compute Cklm averaged over all orientations/polarizations
void CklmCoeff::IniCklmCoeff(XYZGrid& labgrid, Orbital& dysorb,
                          double radfn, KLMPoints& klmgrid)
{
  int lmx = klmgrid.LMax();
  int nk = klmgrid.NKV();
  int nc = nk*klmgrid.NLMV();
  Alloc(nc,lmx,nk);
  
  //Print Dyson and Rkl scans 
  PrintDysonScanZ(labgrid,dysorb,klmgrid);
  
  CalcCklmSq(labgrid,dysorb,radfn,klmgrid);

  PrintAvgCklmSq(klmgrid);
}

//! Compute Cklm for a particular molecular orientation and laser polarization
void CklmCoeff::IniCklmCoeff(XYZGrid& labgrid, Orbital& dysorb,
                             double *RIonz, double radfn,
                             KLMPoints& klmgrid)
{
  int lmx = klmgrid.LMax();
  int nk = klmgrid.NKV();
  int nc = nk*klmgrid.NLMV();
  Alloc(nc,lmx,nk);

  //Print Dyson and Rkl scans 
  PrintDysonScanZ(labgrid,dysorb,klmgrid);

  CalcCklmSq(labgrid,dysorb,RIonz,radfn,klmgrid);

  PrintCklmSq(klmgrid);
}

//! Compute Cklm for a particular molecular orientation and laser polarization
void CklmCoeff::IniCklmCoeff_Overlap(XYZGrid& labgrid, Orbital& dysorb,
				     double radfn,
				     KLMPoints& klmgrid)
{
  int lmx = klmgrid.LMax();
  int nk = klmgrid.NKV();
  int nc = nk*klmgrid.NLMV();
  Alloc(nc,lmx,nk);

  //Print Dyson and Rkl scans 
  //PrintDysonScanZ(labgrid,dysorb,klmgrid);
  
  CalcCklmSq_Overlap(labgrid,dysorb,radfn,klmgrid);

  PrintCklmSq(klmgrid);
}

/*! Compute cklm on numerical grid. 
  Here we compute 3 components (x,y,z) to later obtain a final averaged Cklm. */
void CklmCoeff::CalcCklmSq(XYZGrid& labgrid, Orbital& dysorb,
                           double radfn, KLMPoints& klmgrid)
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

#pragma omp parallel for
  for(int nx=0; nx<nptx; nx++)
    {
      int thread_id=omp_get_thread_num();
      for(int ny=0; ny<npty; ny++)
	{
	  for(int nz=0; nz<nptz; nz++) 
	    {
              double Ldys_value, Rdys_value ;
              unsigned long long int gridaddr = nx*npty*nptz + ny*nptz + nz;
              dysorb.CalcDyson_xyz(Ldys_value,Rdys_value,gridptr_x[nx],gridptr_y[ny],gridptr_z[nz]);

              Complex tLDys_cart, tRDys_cart, totalLDys_cart[XYZ], totalRDys_cart[XYZ];
              tLDys_cart=Ldys_value*sqrt(pow(gridptr_x[nx],2)+pow(gridptr_y[ny],2)+pow(gridptr_z[nz],2)) ;
              tRDys_cart=Rdys_value*sqrt(pow(gridptr_x[nx],2)+pow(gridptr_y[ny],2)+pow(gridptr_z[nz],2)) ;

               for (int k=0, no=0, kl=0; k<nkv; k++)
                {
                double *rklmax = new double[lmax+1];
                Rkl(radfn,rklmax,gridptr_x[nx],gridptr_y[ny],gridptr_z[nz],klmgrid.GetKV(k),lmax);
                for(int l=0; l<=lmax; l++,kl++)
                  {
                    double rkl = rklmax[l];
                    int sign=l%4;
                    Complex phasefactor= (!sign) ? Complex(1,0.0) :
                      ((1==sign)? Complex(0.0, -1) :
                       ((2==sign)? Complex(-1,0.0) : Complex(0.0,1)));
                    for (int m=-l; m<=l; m++, no++)
                      {
                        Complex ylm=sph->GetValue(gridptr_x[nx],gridptr_y[ny],gridptr_z[nz],l,m);
                        ylm = ylm*phasefactor;
                        Complex ylm_conj(ylm.Conj());
                        for (int v=0; v<XYZ; v++)
                          {
                            Complex ylmtwo=sph->GetValue(gridptr_x[nx],gridptr_y[ny],gridptr_z[nz],1,v-1);
                            Complex ylmtwo_conj(ylmtwo.Conj());
                            totalLDys_cart[v]=tLDys_cart*ylm*ylmtwo*rkl;
                            totalRDys_cart[v]=tRDys_cart*ylm_conj*ylmtwo_conj*rkl;
                            Lcklm[v][no][thread_id] += totalLDys_cart[v];
                            Rcklm[v][no][thread_id] += totalRDys_cart[v];
                                }
                            }
                        }
                  delete [] rklmax;
               }
            }
        }
    }

 
  for (int k=0,no=0; k<nkv; k++)
    for (int lm=0; lm<nlmv; lm++,no++)
      for (int v=0; v<XYZ; v++)
	for(int i=1; i<nthreads; i++)
	  {
	    Lcklm[v][no][0] += Lcklm[v][no][i];
	    Rcklm[v][no][0] += Rcklm[v][no][i];
	  }

  double dVV=dV*dV;

  clebsh_gordan_coeff cgc;

  for (int k=0,no=0; k<nkv; k++)
      for (int l=0; l<=lmax; l++)
        for (int m1=-l; m1<=l; m1++, no++)
           for (int m21=-l; m21<=l; m21++)
              for (int m22=-l; m22<=l; m22++)
                 for (int v1=0; v1<XYZ; v1++)
                    for (int v2=0; v2<XYZ; v2++)
                       for (int ltot=(abs(l-1)); ltot<=(l+1); ltot++) {
                          if ((m21+v1)==(m22+v2)){
                             int no1 = k*(lmax+1)*(lmax+1) + l*(l+1) + m21;
                             int no2 = k*(lmax+1)*(lmax+1) + l*(l+1) + m22;
                             double tmp=(Lcklm[v1][no1][0]*Rcklm[v2][no2][0]).Re();
                             tmp*=dVV;
                             double cgc1 = cgc.cgc[l][(ltot-abs(l-1))][m1+l][1] ;
                             double cgc21 = cgc.cgc[l][(ltot-abs(l-1))][m21+l][v1]  ;
                             double cgc22 = cgc.cgc[l][(ltot-abs(l-1))][m22+l][v2]  ;
                             tmp*=cgc1*cgc1*cgc21*cgc22;
                             tmp/=2*ltot+1;
                             cklmsq[no] += tmp ;
                           }
                       }

  for (int k=0,no=0; k<nkv; k++)
      for (int l=0; l<=lmax; l++)
        for (int m11=-l; m11<=l; m11++, no++)
          for (int l2=0; l2<=lmax; l2++)
           for (int m12=-l2; m12<=l2; m12++)
             for (int m21=-l; m21<=l; m21++)
                for (int m22=-l2; m22<=l2; m22++)
                   for (int v1=0; v1<XYZ; v1++)
                      for (int v2=0; v2<XYZ; v2++)
                         for (int ltot=(abs(l-1)); ltot<=(l+1); ltot++) {
                           for (int ltot2=(abs(l2-1)); ltot2<=(l2+1); ltot2++) 
                            if ((m21+v1-1)==(m22+v2-1) && m12==m11 && ltot==ltot2){
                               int no1 = k*(lmax+1)*(lmax+1) + l*(l+1) + m21;
                               int no2 = k*(lmax+1)*(lmax+1) + l2*(l2+1) + m22;
                               double tmp=(Lcklm[v1][no1][0]*Rcklm[v2][no2][0]).Re();
                               tmp*=dVV;
                               double cgc11 = cgc.cgc[l][(ltot-abs(l-1))][m11+l][1] ;
                               double cgc12 = cgc.cgc[l2][(ltot2-abs(l2-1))][m12+l2][1] ;
                               double cgc21 = cgc.cgc[l][(ltot-abs(l-1))][m21+l][v1] ;
                               double cgc22 = cgc.cgc[l2][(ltot2-abs(l2-1))][m22+l2][v2] ;
                               tmp*=cgc11*cgc12*cgc21*cgc22;
                               tmp/=2*ltot+1;
                               crosscklm[no][l2] += tmp;
                             }
                         }
 
  for (int k=0,no=0; k<nkv; k++)
    for (int l=0; l<=lmax; l++)
      for (int m=-l; m<=l; m++, no++)
	sumcklm[k] += cklmsq[no];
  
  for (int v=0; v<XYZ; v++)
    {
      for(int i=0; i<ncklm; i++)
	{
	  delete [] Rcklm[v][i];
	  delete [] Lcklm[v][i];
	}
    }

}

//! Compute Cklm for a particular molecular orientation and laser polarization
void CklmCoeff::CalcCklmSq(XYZGrid& labgrid, Orbital& dysorb,
                           double *RIonz, double radfn, KLMPoints& klmgrid)
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
                {
                double *rklmax = new double[lmax+1];
                Rkl(radfn,rklmax,gridptr_x[nx],gridptr_y[ny],gridptr_z[nz],klmgrid.GetKV(k),lmax);
                for(int l=0; l<=lmax; l++,kl++)
                  {
                    double rkl = rklmax[l];
                    int sign=l%4;
                    Complex phasefactor= (!sign) ? Complex(1,0.0) :
                      ((1==sign)? Complex(0.0, -1) :
                       ((2==sign)? Complex(-1,0.0) : Complex(0.0,1)));
                    for (int m=-l; m<=l; m++, no++)
                      {
                        Complex ylm=sph->GetValue(gridptr_x[nx],gridptr_y[ny],gridptr_z[nz],l,m);
                        ylm = ylm*phasefactor;
                        Complex ylm_conj(ylm.Conj());
                        for (int v=0; v<XYZ; v++)
                          {
                            totalLDys_cart[v]=tLDys_cart[v]*ylm*rkl*unitvect[v];
                            totalRDys_cart[v]=tRDys_cart[v]*ylm_conj*rkl*unitvect[v];
                            Lcklm[v][no][thread_id] += totalLDys_cart[v];
                            Rcklm[v][no][thread_id] += totalRDys_cart[v];
                          }
                      }
                  }
                  delete [] rklmax;
               }
            }
        }
    }

   for (int k=0,no=0; k<nkv; k++)
     for (int lm=0; lm<nlmv; lm++,no++)
       for (int v=0; v<XYZ; v++)
       for(int i=1; i<nthreads; i++)
         {
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
          cklmsq[no] += tmp*3./(4.*M_PI);
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
                 crosscklm[no][l2] += tmp*3./(4.*M_PI);
               }
            }
    }

  for (int k=0,no=0; k<nkv; k++)
    for (int l=0; l<=lmax; l++)
      for (int m=-l; m<=l; m++, no++)
          sumcklm[k] += cklmsq[no];

  for (int v=0; v<XYZ; v++)
    {
       for(int i=0; i<ncklm; i++)
       {
         delete [] Rcklm[v][i];
         delete [] Lcklm[v][i];
       }
    }
}

//! Compute overlap of DO and PW for a particular molecular orientation
void CklmCoeff::CalcCklmSq_Overlap(XYZGrid& labgrid, Orbital& dysorb,
				   double radfn, KLMPoints& klmgrid)
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

  Complex  **Lcklm[XYZ], **Rcklm[XYZ];

   int nthreads=1;
   nthreads=omp_get_max_threads();
   printf("Executing CalcCklmSq_Overlap on %d threads\n",nthreads); fflush(stdout);

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

  //May be able to optimize code below
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
                {
                double *rklmax = new double[lmax+1];
                Rkl(radfn,rklmax,gridptr_x[nx],gridptr_y[ny],gridptr_z[nz],klmgrid.GetKV(k),lmax);
                for(int l=0; l<=lmax; l++,kl++)
                  {
                    double rkl = rklmax[l];
                    int sign=l%4;
                    Complex phasefactor= (!sign) ? Complex(1,0.0) :
                      ((1==sign)? Complex(0.0, -1) :
                       ((2==sign)? Complex(-1,0.0) : Complex(0.0,1)));
                    for (int m=-l; m<=l; m++, no++)
                      {
                        Complex ylm=sph->GetValue(gridptr_x[nx],gridptr_y[ny],gridptr_z[nz],l,m);
                        ylm = ylm*phasefactor;
                        Complex ylm_conj(ylm.Conj());
                        for (int v=0; v<XYZ; v++)
                          {
                            totalLDys_cart[v]=tLDys_cart[v]*ylm*rkl;
                            totalRDys_cart[v]=tRDys_cart[v]*ylm_conj*rkl;
                            Lcklm[v][no][thread_id] += totalLDys_cart[v];
                            Rcklm[v][no][thread_id] += totalRDys_cart[v];
                          }
                      }
                  }
                  delete [] rklmax;
               }
            }
        }
    }

   for (int k=0,no=0; k<nkv; k++)
     for (int lm=0; lm<nlmv; lm++,no++)
       for (int v=0; v<XYZ; v++)
       for(int i=1; i<nthreads; i++)
         {
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
          cklmsq[no] += tmp*3./(4.*M_PI);
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
                 crosscklm[no][l2] += tmp*3./(4.*M_PI);
               }
            }
    }

  for (int k=0,no=0; k<nkv; k++)
    for (int l=0; l<=lmax; l++)
      for (int m=-l; m<=l; m++, no++)
          sumcklm[k] += cklmsq[no];

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
  fprintf(fout,"k,a.u.         l      m     Re(|Cklm|^2)\n");

  int nlmv = klmgrid.NLMV();
  int no;
  double eng, kau;
  for(int k=0,no=0; k<nkv; k++)
    for (int lm=0; lm<nlmv; lm++,no++)
      {
      fprintf(fout,"%.3lE  %5d  %5d",klmgrid.GetKV(k),klmgrid.GetLV(lm),klmgrid.GetMV(lm));
      fprintf(fout," %17.6lE\n",(cklmsq[no]));
      }
  fprintf(fout,"\n  E_k,eV     Sum(|(Cklm)|^2)\n");
  for (int k=0; k<nkv; k++)
    {
      kau = klmgrid.GetKV(k);
      eng = 0.5*kau*kau*HAR_TO_EV;
      fprintf(fout,"%7.3lf",eng);
      fprintf(fout," %17.6lE\n",(sumcklm[k]));
    }
}

//! Print all square Cklm coeffs, ||Cklm^2||.
void CklmCoeff::PrintCklmSq(KLMPoints& klmgrid, FILE *fout)
{
  fprintf(fout,"\n||Cklm^2|| printout\n");
  fprintf(fout,"k,a.u.         l      m      Re(|Cklm|^2)\n");

  int nlmv = klmgrid.NLMV();
  int no;
  double eng, kau;
  for(int k=0,no=0; k<nkv; k++)
    for (int lm=0; lm<nlmv; lm++,no++)
      {
      fprintf(fout,"%.3lE  %5d  %5d",klmgrid.GetKV(k),klmgrid.GetLV(lm),klmgrid.GetMV(lm));
      fprintf(fout," %17.6lE\n",(cklmsq[no]));
      }
  fprintf(fout,"\n  E_k,eV     Sum(|(Cklm)|^2)\n");
  for (int k=0; k<nkv; k++)
    {
      kau = klmgrid.GetKV(k);
      eng = 0.5*kau*kau*HAR_TO_EV;
      fprintf(fout,"%7.3lf",eng);
      fprintf(fout," %17.6lE\n",(sumcklm[k]));
    }

#if 0  
//SG: disabled printing of cross-terms
  
  fprintf(fout,"\nAverage Cross-terms\n");
  fprintf(fout,"   L     M     L' Re(CklmCkl'm),a.u.\n");
  for (int k=0,no=0; k<nkv; k++)
    for (int l=0; l<=lmax; l++)
      for (int m=-l; m<=l; m++, no++)
	{
	  fprintf(fout,"%4d  %+4d  ",l,m);
        for (int l2=0; l2<=lmax; l2++)
	  {
	    fprintf(fout,"%4d",l2);
	    //for (int v=0; v<XYZ; v++)
	      fprintf(fout,"  %+.3lE",crosscklm[no][l2]);
	  }
	fprintf(fout,"\n");
       }
  fflush(fout);
#endif  
}


void CklmCoeff::PrintDysonScanZ(XYZGrid& labgrid, Orbital& dysorb, KLMPoints& klmgrid)
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


