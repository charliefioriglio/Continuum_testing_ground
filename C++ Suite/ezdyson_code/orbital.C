#include "orbital.h"


void Orbital::Alloc(int nb)
{
  check(!(nb>0),"Orbital::Alloc() : Invalid nbasis\n");
  
  Free();
  
  nbasis=nb;
  Lcoeff=new double[nbasis];
  Rcoeff=new double[nbasis];
  zeromap=new short[nbasis];
}

void Orbital::Free()
{
  if(IfAlloc())
   {
     if (NULL!=Lcoeff)
      {
        delete [] Lcoeff;
        Lcoeff=NULL;
      }
     if (NULL!=Rcoeff)
      {  
        delete [] Rcoeff;
        Rcoeff=NULL;
      }
     if (NULL!=zeromap)
       {  
	 delete [] zeromap;
	 zeromap=NULL;
      }
     nbasis=0;
   }
}

void Orbital::Set()
{
  nbasis=0;
  Lcoeff=NULL;
  Rcoeff=NULL;
  zeromap=NULL;
}


void Orbital::AllocAndCopyData(int nb, bool ifr, int tran, double ln, double rn, double *Lc, double *Rc, short *zm)
{
  Alloc(nb);
  ifrestr = ifr;
  transno = tran;
  Lnorm = ln;
  Rnorm = rn;
  memcpy(Lcoeff,Lc,sizeof(double)*nbasis);
  memcpy(Rcoeff,Rc,sizeof(double)*nbasis);
  memcpy(zeromap,zm,sizeof(short)*nbasis);
}

void Orbital::SetZeroMap(double thresh)
{
  for(int i=0; i<nbasis; i++)
    zeromap[i]=(fabs(Lcoeff[i])<thresh && fabs(Rcoeff[i])<thresh) ? 0 : 1;
}

Orbital::Orbital() 
{
  Set();
}

Orbital::Orbital(const Orbital& other)
{
  Set();
  AllocAndCopyData(other.nbasis,other.ifrestr,other.transno,other.Lnorm,other.Rnorm,other.Lcoeff,other.Rcoeff,other.zeromap);
}

Orbital::Orbital(int nb, bool ifr, int tran, double ln, double rn, double *Lc, double *Rc)
{
  Set();
  Alloc(nb);
  ifrestr = ifr;
  transno = tran;
  Lnorm = ln;
  Rnorm = rn;
  memcpy(Lcoeff,Lc,sizeof(double)*nbasis);
  memcpy(Rcoeff,Rc,sizeof(double)*nbasis);
  //AllocAndCopyData(nb,ifr,tran,Lc,Rc);
  SetZeroMap();
}


#if 0
//Testing algorithm #2, described in the Dyson_mc.tex: Norms are conserved, but the results are crazy
Orbital::Orbital(const Orbital& other, int atom)
{
  Set();
  Alloc(other.nbasis);
  ifrestr = other.ifrestr;
  transno = other.transno;
  
  memset(Lcoeff,0.0,sizeof(double)*nbasis);
  memset(Rcoeff,0.0,sizeof(double)*nbasis);
  Lnorm=0.0;
  Rnorm=0.0;
  
  AOBasis *pbasis=&(TheAOBasis());
  int natoms=pbasis->NAtoms();
  check(!(atom<natoms), "Orbital::Orbital(const Orbital& other, int atom) : Invalid atom="+atom);

  arma::mat X=pbasis->GetLowdinTransf_X();
  arma::mat Xm=pbasis->GetLowdinTransf_Xm();

  arma::mat right_orb(other.Rcoeff,nbasis,1,true,false);
  arma::mat left_orb(other.Lcoeff,nbasis,1,true,false);
  
  std::cout <<  std::setw(8) << std::setprecision(5) << std::fixed;

#if 0
  //Check overlap
  arma::mat S=pbasis->GetOverlap();
  arma::mat tmp=trans(right_orb)*S*right_orb;
  tmp.print("Proper right norm^2: Rt_orbxSxR_orb");
  tmp=trans(left_orb)*S*left_orb;
  tmp.print("Proper left norm^2: LtorbxSxL_orb");
#endif
  
  arma::mat symm_right_orb=Xm*right_orb;
  arma::mat symm_left_orb=Xm*left_orb; 
  
  std::cout << "Right norm in AO=" << norm(right_orb) << " norm in MO=" << norm(symm_right_orb) << std::endl;
  std::cout << "Left norm in AO=" << norm(left_orb) << " norm in MO=" << norm(symm_left_orb) << std::endl;
  
  //right_orb.raw_print("Right orb");
  //symm_right_orb.raw_print("Symm Right orb");
  //end debug printing
  
  //Split the symmetric orbital according to raw atomic contributions
  for(int i=0, nb=0; i<natoms; i++) {
    if(!(i==atom))
      for(int j=0; j<pbasis->GetGTOsAtom(i); j++) {
	symm_left_orb[nb+j]=0.0;
	symm_right_orb[nb+j]=0.0;
      }
    nb+=pbasis->GetGTOsAtom(i);
  }

  Lnorm=norm(symm_left_orb);
  Rnorm=norm(symm_right_orb);
  
  //std::cout << "Symm atomic Rnorm=" << Rnorm <<  " and  Lnorm=" <<  Lnorm << std::endl;

  //Normalize them to one
  symm_left_orb/=Lnorm;
  symm_right_orb/=Rnorm;

  //Now transform them back
  arma::mat final_right_orb=X*symm_right_orb;
  arma::mat final_left_orb=X*symm_left_orb; 

#if 0  
  //Check overlaps:
  tmp=trans(final_right_orb)*S*final_right_orb;
  tmp.print("Proper right norm^2: Rt_orbxSxR_orb");
  tmp=trans(final_left_orb)*S*final_left_orb;
  tmp.print("Proper left norm^2: LtorbxSxL_orb");
#endif
  
  //Now copy them back to Rcoeff and Lcoeff
  memcpy(Rcoeff,final_right_orb.memptr(),sizeof(double)*nbasis);
  memcpy(Lcoeff,final_left_orb.memptr(),sizeof(double)*nbasis);
  
  //These are proper norms now
  Lnorm = other.Lnorm*Lnorm;
  Rnorm = other.Rnorm*Rnorm;

  SetZeroMap();

  //Print();
  //exit(1);
}
#endif


#if 1 //Another algorithm using stricly localized orbitals 
//This code can be cleaned to remove uneccesary debug calculations
Orbital::Orbital(const Orbital& other, int atom)
{
  Set();
  Alloc(other.nbasis);
  ifrestr = other.ifrestr;
  transno = other.transno;
  
  memset(Lcoeff,0.0,sizeof(double)*nbasis);
  memset(Rcoeff,0.0,sizeof(double)*nbasis);
  Lnorm=0.0;
  Rnorm=0.0;
  
  AOBasis *pbasis=&(TheAOBasis());
  int natoms=pbasis->NAtoms();
  check(!(atom<natoms), "Orbital::Orbital(const Orbital& other, int atom) : Invalid atom="+atom);

  arma::mat X=pbasis->GetLowdinTransf_X();
  arma::mat Xm=pbasis->GetLowdinTransf_Xm();
  arma::mat S=pbasis->GetOverlap();
  std::cout <<  std::setw(8) << std::setprecision(5) << std::fixed;

  //Make coloumns
  arma::mat right_orb(other.Rcoeff,nbasis,1,true,false);
  arma::mat left_orb(other.Lcoeff,nbasis,1,true,false);
  
  //Check overlap
  arma::mat tmp=trans(right_orb)*S*right_orb;
  tmp.print("Proper right norm^2: Rt_orbxSxR_orb");
  tmp=trans(left_orb)*S*left_orb;
  tmp.print("Proper left norm^2: LtorbxSxL_orb");
  
  arma::mat symm_right_orb=Xm*right_orb;
  arma::mat symm_left_orb=Xm*left_orb; 
  
  std::cout << "Right norm in AO=" << norm(right_orb) << " norm in MO=" << norm(symm_right_orb) << std::endl;
  std::cout << "Left norm in AO=" << norm(left_orb) << " norm in MO=" << norm(symm_left_orb) << std::endl;
  
  //Split the orbital according to raw atomic contributions
  for(int i=0, nb=0; i<natoms; i++) {
    //std::cout << "i="<<i << " nbas[i]=" <<  pbasis->GetGTOsAtom(i) << " offset=" << nb <<std::endl;
    if(i==atom)
      for(int j=0; j<pbasis->GetGTOsAtom(i); j++) {
	Lcoeff[nb+j]=other.Lcoeff[nb+j];
	Rcoeff[nb+j]=other.Rcoeff[nb+j];
      }
    nb+=pbasis->GetGTOsAtom(i);
  }

  //Now orbital is localized, but it is not normalzied
  arma::mat atomic_right_orb(Rcoeff,nbasis,1,true,false);
  arma::mat atomic_left_orb(Lcoeff,nbasis,1,true,false);
  
  Rnorm=norm(atomic_right_orb);
  Lnorm=norm(atomic_left_orb);
  
  std::cout << "Norms in AOs:  Rnorm=" << Rnorm <<  " and  Lnorm=" <<  Lnorm << std::endl;
  
  arma::mat atomic_symm_right_orb=Xm*atomic_right_orb;
  arma::mat atomic_symm_left_orb=Xm*atomic_left_orb; 
  Rnorm=norm(atomic_symm_right_orb);
  Lnorm=norm(atomic_symm_left_orb);
  //The same as ao.t()xSxao 
  std::cout << "Proper norms in MO:  Rnorm=" << Rnorm <<  " and Lnorm=" <<  Lnorm << std::endl;

  //Now let us normalize  local left and right orbitals
  for(int i=0; i<nbasis; i++) {
    Rcoeff[i]/=Rnorm;
    Lcoeff[i]/=Lnorm;
  }
  //They are now properly normalized to one.
  
#if 0 //uncomment this for norm-preserving version
  //All is good now, but the norms are not additive, because we are missing overlap contribution.
  //Lets compute it:
  double roverlap = dot(symm_right_orb,atomic_symm_right_orb);
  double loverlap = dot(symm_left_orb,atomic_symm_left_orb);
  std::cout << "Overlap between localized and original DO:  Right=" << roverlap <<  " and Left=" <<  loverlap << std::endl;
  //Use these overlaps instead of proper norms -- this way the norms of localized orbitals add up
  std::cout << "Use overlap to fix the norm... The local orbitals are properly normalized and total norm is conserved." << std::endl; 
  Rnorm=sqrt(roverlap);
  Lnorm=sqrt(loverlap);
#endif  
  
  //These are proper norms now
  Lnorm = other.Lnorm*Lnorm;
  Rnorm = other.Rnorm*Rnorm;
  
  SetZeroMap();

  //Print();
  //exit(1);
}
#endif


Orbital::~Orbital() 
{ 
  Free();
}

Orbital& Orbital::operator=(const Orbital& other)
{
  if(this!=&other)
    AllocAndCopyData(other.nbasis,other.ifrestr,other.transno,other.Lnorm,other.Rnorm,other.Lcoeff,other.Rcoeff,other.zeromap);
  return *this; 
}


//! Print Orbital coeffs in AO basis.
void Orbital::Print(FILE *fout) const
{
  fprintf(fout,"Left norm (sqrt|dyson_l x dyson_l|)=%lf\n",Lnorm);
  fprintf(fout,"Right norm (sqrt|dyson_r x dyson_r|)=%lf\n",Rnorm);
  fprintf(fout,"LeftxRight=%lf\n",Lnorm*Rnorm);

  fprintf(fout,"\nDysonMOs in AO basis\n");
  fprintf(fout,"  AO     [LEFT]              [RIGHT]        ZERO_MAP\n");
  for (int i=0; i<nbasis; i++)
    fprintf(fout,"%3d %15.8lf     %15.8lf      %d\n",i+1,Lcoeff[i],Rcoeff[i],zeromap[i]);
  fprintf(fout,"-----------------------------------------------------------------------------------\n");
  fflush(fout);
}

//! Add up all basis fns with proper coefs to get orbital value at x,y,z.
void Orbital::CalcDyson_xyz(double &Ldys_value, double &Rdys_value, double x, double y, double z) const
{
  static AOBasis *pbasis=&(TheAOBasis());
  static int nbasis = pbasis->NBasis();

  double onegto;

  //Initialize to zero.
  Ldys_value = 0.0;
  Rdys_value = 0.0;

  for(int n=0,k=0,i=0; i<nbasis; i++,k++)
    {
      if(k==pbasis->GetGTOsAtom(n))
       {
         n++;
         k=0;
       }
     
      if(zeromap[i]) //Do calcs only if non-zero coefficient
	{
	  onegto = pbasis->CalcGaussValue(i,n,x,y,z);
	  Ldys_value += Lcoeff[i]*onegto;
	  Rdys_value += Rcoeff[i]*onegto;
	}
    }

  //Slater test
/*    double rsq = x*x+y*y+z*z;
    Ldys_value = 0.23*z*(2.57*exp(-20.96*rsq)+2.41*exp(-4.8*rsq)+1.865*exp(-1.46*rsq));
    Ldys_value +=0.38*z*0.575*exp(-0.483*rsq)+0.43*z*0.128*exp(-0.146*rsq);
    Ldys_value +=0.22*z*0.02856*exp(-0.0483*rsq)+0.05*z*0.00637*exp(-0.01319*rsq);

    Ldys_value = z*exp(-1.2*BOHR_TO_ANGS*sqrt(rsq));
    Rdys_value = Ldys_value;  */
}

double Orbital::GetOrbitalNormAndCenter(double *xyzcenter, const XYZGrid& grid) const
{
  double ldys,rdys,ldyson_sq,rdyson_sq;
  double lnorm=0.0,rnorm=0.0,lrnorm=0.0, dV=grid.DXYZ(X)*grid.DXYZ(Y)*grid.DXYZ(Z);
  double xi, yi, zi;
  int i;

  memset(xyzcenter,0,sizeof(double)*XYZ);

  for(i=0; i<grid.NXYZ(); i++)
    {
      grid.GetPoint(i,xi,yi,zi,X,Y,Z);
      CalcDyson_xyz(ldys,rdys,xi,yi,zi);
      ldyson_sq=ldys*ldys;
      rdyson_sq=rdys*rdys;
      lnorm+=ldyson_sq;
      rnorm+=rdyson_sq;
      lrnorm+=ldys*rdys;
      
      xyzcenter[X] += xi*ldyson_sq;
      xyzcenter[Y] += yi*ldyson_sq;
      xyzcenter[Z] += zi*ldyson_sq;
    }
  
  for (i=0; i<XYZ; i++)
    xyzcenter[i]/=lnorm;

  lnorm*=dV;
  rnorm*=dV;
  lrnorm*=dV;
  lnorm=sqrt(lnorm);
  rnorm=sqrt(rnorm);
  lrnorm=sqrt(lrnorm);

  fprintf(stdout,"\nNorm of the left Dyson orbital integrated on the grid: %15.8lf\n", lnorm);
  fprintf(stdout,"Norm of the right Dyson orbital integrated on the grid: %15.8lf\n", rnorm);
  fprintf(stdout,"Left-Right Norm of the Dyson orbital integrated on the grid: %15.8lf\n", lrnorm);
  fprintf(stdout,"Xavg Yavg Zavg:  %lf  %lf  %lf a.u.\n",xyzcenter[X],xyzcenter[Y],xyzcenter[Z]);
  fprintf(stdout,"-----------------------------------------------------------------------------------\n");

  //Now compute X^2 Y^2 and Z^2
  double second_mom[XYZ];
  second_mom[X]=second_mom[Y]=second_mom[Z]=0.0;

  for(i=0; i<grid.NXYZ(); i++)
    {
      grid.GetPoint(i,xi,yi,zi,X,Y,Z);
      CalcDyson_xyz(ldys,rdys,xi,yi,zi);
      ldyson_sq=ldys*ldys;
      lnorm+=ldyson_sq;

      xi-=xyzcenter[X];
      yi-=xyzcenter[Y];
      zi-=xyzcenter[Z];

      second_mom[X] += xi*xi*ldyson_sq;
      second_mom[Y] += yi*yi*ldyson_sq;
      second_mom[Z] += zi*zi*ldyson_sq;
    }

  double size=0.0;

  for (i=0; i<XYZ; i++)
    {
      second_mom[i]/=lnorm;
      second_mom[i]*=(BOHR_TO_ANGS*BOHR_TO_ANGS);
      size+=second_mom[i];
    }

  fprintf(stdout,"<X^2> <Y^2> <Z^2> <R^2>:  %lf  %lf  %lf %lf Angs\n",
	  second_mom[X],second_mom[Y],second_mom[Z],size);
  fprintf(stdout,"-----------------------------------------------------------------------------------\n");
  fflush(stdout);
  return lnorm;
}

//!  Calculate overlap with another orbital
void Orbital::GetOverlap(const Orbital& other, const XYZGrid& grid, double& s_rr, double& s_ll) const {
  
  double ldys,rdys;
  double other_ldys, other_rdys;
  double dV=grid.DXYZ(X)*grid.DXYZ(Y)*grid.DXYZ(Z);
  double xi, yi, zi;
  int i;

  s_rr=s_ll = 0.0;
  
  for(i=0; i<grid.NXYZ(); i++)
    {
      grid.GetPoint(i,xi,yi,zi,X,Y,Z);
      CalcDyson_xyz(ldys,rdys,xi,yi,zi);
      other.CalcDyson_xyz(other_ldys,other_rdys,xi,yi,zi);
      s_rr+=rdys*other_rdys;
      s_ll+=ldys*other_ldys;
     }
  
  s_rr*=dV;
  s_ll*=dV;

  fprintf(stdout,"\nOverlap of the two right Dyson orbitals integrated on the grid: %15.8lf\n", s_rr);
  fprintf(stdout,"\nOverlap of the two right Dyson orbitals integrated on the grid: %15.8lf\n", s_ll);
  fflush(stdout);
}


