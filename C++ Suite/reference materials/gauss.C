#include "gauss.h"

#define MYPOW(x,x2,n) ((n==0) ? 1: ((n==1) ? x : ((n==2) ? x2 :  ((n==3) ? (x2*x) :  ((n==4) ? (x2*x2) : pow(x,n))))))

void Gauss::Alloc(int nc)
{
  check(!(nc>0), "Gauss::Alloc() : Invalid ncontr\n");

  Free();

  ncontr=nc;
  alpha=new double[ncontr];
  coeff=new double[ncontr];
  norm=new double[ncontr];
  coeffnorm=new double[ncontr];
}

void Gauss::Free()
{
  if(IfAlloc())
   {
     delete [] alpha;
     delete [] coeff;
     delete [] norm;
     delete [] coeffnorm;
     ncontr=0;
   }
}

Gauss::Gauss(const Gauss& other) : ncontr(0) 
{
  IniGauss(other.angmom,other.ncontr,other.ia,other.ja,other.ka,other.dx2y2,other.dz2,other.puref,other.fa,
           other.pureg,other.ga,other.alpha,other.coeff,1.0);
}

void Gauss::CheckAngMom() const
{
  int ijk=ia+ja+ka;
  if(!dx2y2 && !puref && !pureg)
    check(!ijk==angmom,"CheckAngMom() : angular momentum does not match\n");
  else
    {
      if(dx2y2)
	check(2!=angmom,"CheckAngMom() : Should be d-function\n");
      if(puref)
	check(3!=angmom,"CheckAngMom() : Should be f-function\n");
      if(pureg)
        check(4!=angmom,"CheckAngMom() : Should be g-function\n");
    }
}

Gauss::~Gauss() 
{ 
  Free();
}

void Gauss::IniGauss(int angm, int nc, int i, int j, int k, int dxy, int dz, int purf, int f, int purg, int g, double *a, double *c, double sc)
{
  check(IfAlloc(),"Gauss::IniGauss() : Already initialized\n");
  Alloc(nc);
  angmom=angm;
  memcpy(alpha,a,sizeof(double)*ncontr);
  for(int n=0; n<ncontr; n++)
    alpha[i]*=sc;
  memcpy(coeff,c,sizeof(double)*ncontr);
  ia = i; ja = j; ka = k; dx2y2 = dxy; dz2 = dz; puref=purf; fa=f; pureg=purg; ga=g;
  GaussNorm();
  GaussCoeffNorm();
}

Gauss& Gauss::operator=(const Gauss& other)
{
  if (this!=&other)
    
  IniGauss(other.angmom,other.ncontr,other.ia,other.ja,other.ka,other.dx2y2,other.dz2,other.puref,other.fa,
           other.pureg,other.ga,other.alpha,other.coeff,1.0);
  return *this; 
}

double Gauss::TotalNorm() const
{
  double totalnorm = 0.0;
  
  for (int i=0; i<ncontr; i++)
    for (int j=i; j<ncontr; j++)
      {
	double tot_coeff=(j==i)? 1.0 : 2.0;
	tot_coeff*=coeff[i]*coeff[j]*norm[i]*norm[j];
	double alpha_ij=0.5*(alpha[i]+alpha[j]);

	if (2==dx2y2)
	  totalnorm += tot_coeff*PureNorm(2,alpha_ij);
	else 
	  if (2==dz2)
	    totalnorm += tot_coeff*PureNorm(2,alpha_ij);
	  else 
	    if (3==puref)
	      {
                  totalnorm += tot_coeff*PureNorm(3,alpha_ij);
	      }
             else
               if (4==pureg)
                {
                    totalnorm += tot_coeff*PureNorm(4,alpha_ij);
              }
	    else //S, P, and all  pure carts
	      totalnorm += tot_coeff*CartNorm(ia,ja,ka,alpha_ij);
      }
  return totalnorm;
}

void Gauss::GaussNorm() const
{
  for (int i=0; i<ncontr; i++)
    {
      if (2==dx2y2)
	norm[i] = 1./sqrt(PureNorm(2,alpha[i]));
      else 
	if (2==dz2)
	  norm[i] = 1./sqrt(PureNorm(2,alpha[i]));
	else 
	  if (3==puref)
	    {
		norm[i] = 1./sqrt(PureNorm(3,alpha[i]));
	    }  
          else
            if (4==pureg)
              {
                  norm[i] = 1./sqrt(PureNorm(4,alpha[i]));
            }
	  else
	    norm[i] = 1./sqrt(CartNorm(ia,ja,ka,alpha[i]));
    }
}

//Precompute coeff*norm for each uncontracted GTO.
void Gauss::GaussCoeffNorm() const
{
  //Calc total norm of the contracted gaussian, everything included
  double total=sqrt(TotalNorm());
 
  for (int i=0; i<ncontr; i++)
    coeffnorm[i]= coeff[i]*norm[i]/total;
}


//Calc value of gaussian at x,y,z.
double Gauss::GaussValue(double x, double y, double z, double x0, double y0, double z0) const
{
  double Xat=x-x0,Yat=y-y0,Zat=z-z0;
  double cart_gto=0.0;
  double Xat2=Xat*Xat;
  double Yat2=Yat*Yat;
  double Zat2=Zat*Zat;
  double rsqscaled = Xat2+Yat2+Zat2;

  if(2==dx2y2)
    {
      double rx2y2 = sqrt(3)/2* (Xat2-Yat2); //pow(Xat,2)-pow(Yat,2);
      for (int i=0; i<ncontr; i++)
	cart_gto += coeffnorm[i]*rx2y2*exp(-alpha[i]*rsqscaled);
    }
  else 
    if (2==dz2)
      {
	double rz2 = 0.5*(2*Zat2 - Xat2 - Yat2); //2*pow(Zat,2)-pow(Xat,2)-pow(Yat,2);
	for (int i=0; i<ncontr; i++)
	  cart_gto += coeffnorm[i]*rz2*exp(-alpha[i]*rsqscaled);
      }
    else if (3==puref)
      {
	double rf=0.0;
	switch (fa)
	  {
	  case 1:
            rf = sqrt(5)/sqrt(8)*(Yat*(3*Xat2-Yat2)); //Yat*(pow(Zat,2)-pow(Xat,2));
	    break;
	  case 3:
            rf = sqrt(3)/sqrt(8)*(4*Yat*Zat2-Yat*Xat2-Yat*Yat2); //2*pow(Yat,3)-3*pow(Xat,2)*Yat-3*pow(Zat,2)*Yat;
	    break;
	  case 4:
            rf = 0.5*(2*Zat2*Zat-3*Xat2*Zat-3*Yat2*Zat); //2*pow(Zat,3)-3*pow(Xat,2)*Zat-3*pow(Yat,2)*Zat;
	    break;
          case 5:
            rf = sqrt(3)/sqrt(8)*(4*Xat*Zat2-Xat*Xat2-Xat*Yat2); //2*pow(Xat,3)-3*pow(Yat,2)*Xat-3*pow(Zat,2)*Xat;
            break;
	  case 6:
            rf = sqrt(15)/2*(Zat*(Xat2-Yat2)); //Zat*(pow(Xat,2)-pow(Yat,2));
	    break;
	  case 7:
            rf = sqrt(5)/sqrt(8)*(Xat*(Xat2-3*Yat2)); //Xat*(pow(Zat,2)-pow(Yat,2));
	    break;
	  }    
	for (int i=0; i<ncontr; i++)
	  cart_gto += coeffnorm[i]*rf*exp(-alpha[i]*rsqscaled);
      }
    else if (4==pureg)
      {
        double rg=0.0;
        switch (ga)
          {
          case 1:
            rg = sqrt(35)/sqrt(4)*(Xat2*Xat*Yat-Xat*Yat2*Yat);
            break;
          case 2:
            rg = sqrt(35)/sqrt(8)*(3*Xat2*Yat*Zat-Yat2*Yat*Zat);
            break;
          case 3:
            rg = sqrt(5)/sqrt(4)*(6*Xat*Yat*Zat2-Xat2*Xat*Yat-Xat*Yat2*Yat);
            break;
          case 4:
            rg = sqrt(5)/sqrt(8)*(4*Xat*Zat2*Zat-3*Xat2*Xat*Zat-3*Xat*Yat2*Zat);
            break;
          case 5:
            rg = 1./sqrt(64)*(8*Zat2*Zat2+3*Xat2*Xat2+3*Yat2*Yat2-24*Xat2*Zat2-24*Yat2*Zat2+6*Xat2*Yat2);
            break;
          case 6:
            rg = sqrt(5)/sqrt(8)*(4*Yat*Zat2*Zat-3*Xat2*Yat*Zat-3*Yat2*Yat*Zat);
            break;
          case 7:
            rg = sqrt(5)/sqrt(16)*(6*Xat2*Zat2-Xat2*Xat2-6*Yat2*Zat2+Yat2*Yat2);
            break;
          case 8:
            rg = sqrt(35)/sqrt(8)*(Xat2*Xat*Zat-3*Xat*Yat2*Zat);
            break;
          case 9:
            rg = sqrt(35)/sqrt(64)*(Xat2*Xat2+Yat2*Yat2-6*Xat2*Yat2);
            break;
          }
        for (int i=0; i<ncontr; i++)
          cart_gto += coeffnorm[i]*rg*exp(-alpha[i]*rsqscaled);
      }
    else
      {
	//pow(Xat,ia)*pow(Yat,ja)*pow(Zat,ka);
	double rxyz = MYPOW(Xat,Xat2,ia)*MYPOW(Yat,Yat2,ja)*MYPOW(Zat,Zat2,ka);
	for (int i=0; i<ncontr; i++)   
	  cart_gto += coeffnorm[i]*rxyz*exp(-alpha[i]*rsqscaled);
      }
  return cart_gto;
}


double CartNorm(int ib, int jb, int kb, double alph)
{
  double alphasc = 2.*alph;
  double gnorm = pow(M_PI/alphasc,1.5)*dfac(2*ib-1)*dfac(2*jb-1)*dfac(2*kb-1);
  gnorm /= pow(2.*alphasc,ib+jb+kb);
  return gnorm;
}

double PureNorm(int lpure, double alph)
{
  double alphasc = 2.*alph;
  double gnorm = pow(M_PI/alphasc,1.5)*dfac(2*lpure-1);
  gnorm /= pow(2.*alphasc,lpure);
  return gnorm;
}


