#include "rotnmatr.h"

RotnMatr::RotnMatr()
{
  memset(rotn,0,XYZ*XYZ*sizeof(double));
}

RotnMatr::RotnMatr(const RotnMatr& other) 
{
  memcpy(rotn,other.rotn,XYZ*XYZ*sizeof(double));
}

//RotnMatr::RotnMatr(MolAvg moav, double alpha) 
//{
//  CalcRotnMatr(moav,alpha);
//}

RotnMatr::RotnMatr(MolAvg moav, double alpha, double beta, double gamma) 
{
  CalcRotnMatr(alpha,beta,gamma);
}

RotnMatr& RotnMatr::operator=(const RotnMatr& other)
{
  if(this!=&other)
      memcpy(rotn,other.rotn,XYZ*XYZ*sizeof(double));
  return *this; 
}

//! Calc rotation matrix for cyllindrical rotation (one angle: alpha).
/*void RotnMatr::CalcRotnMatr(MolAvg molavg, double alpha)
{
  switch (molavg)
    {
    case CYLX:
      XRotnMatr(alpha);
      break;
    case CYLY:
      YRotnMatr(alpha);
      break;
    case CYLZ:
      ZRotnMatr(alpha);
      break;
    default:
      perror("RotnMatr::CalcRotnMatr() : invalid molavg");
      break;
    }
}
*/
//! Calc rotation matrix (Euler transformation) for alpha, beta, gamma Euler angles.
void RotnMatr::CalcRotnMatr(double alpha, double beta, double gamma)
{
  EulerRotnMatr(alpha,beta,gamma);
}


void RotnMatr::Print(const char *filename) const 
{
  FILE *outfile = fopen(filename,"a");
  fprintf(outfile,"Rotational matrix\n");
  int addr=0;
  for (int i=0; i<XYZ; i++)
    {
      for (int j=0; j<XYZ; j++)
	fprintf(outfile,"%8.2lf",rotn[addr++]);
      fprintf(outfile,"\n");
    }
  fprintf(outfile,"\n");
  fclose(outfile);
}


//!Calc rotn 3x3 Euler rotation matrix for a a certain molec orientation: alpha,beta,gamma.
RotnMatr&  RotnMatr::EulerRotnMatr(double alpha, double beta, double gamma)
{
/*  double c1=cos(alpha), s1=sin(alpha);
  double c2=cos(beta), s2=sin(beta);
  double c3=cos(gamma), s3=sin(gamma);

  rotn[0] = c2*c3;
  rotn[1] = s1*s3-c1*c3*s2;
  rotn[2] = c3*s1*s2+c1*s3;

  rotn[3] = s2;
  rotn[4] = c1*c2;
  rotn[5] = -c2*s1;

  rotn[6] = -c2*c3;
  rotn[7] = c3*s1+c1*s2*s3;
  rotn[8] = c1*c3-s1*s2*s3;      */


  double cosA=cos(alpha), sinA=sin(alpha);
  double cosB=cos(beta), sinB=sin(beta);
  double cosC=cos(gamma), sinC=sin(gamma);

  rotn[0] = cosC*cosA-cosB*sinA*sinC;
  rotn[1] = cosC*sinA+cosB*cosA*sinC;
  rotn[2] = sinC*sinB;

  rotn[3] = (-sinC)*cosA-cosB*sinA*cosC;
  rotn[4] = (-sinC)*sinA+cosB*cosA*cosC;
  rotn[5] = cosC*sinB;

  rotn[6] = sinB*sinA;
  rotn[7] = (-sinB)*cosA;
  rotn[8] = cosB;         

  return *this;
}

 
//!Calc rotn 3x3 rotation matrix around x by alpha.
RotnMatr&  RotnMatr::XRotnMatr(double alpha)
{
  rotn[0] = 1.0;
  rotn[1] = 0.0;
  rotn[2] = 0.0;

  rotn[3] = 0.0;
  rotn[4] = cos(alpha);
  rotn[5] = sin(alpha);
  
  rotn[6] = 0.0;
  rotn[7] = -sin(alpha);
  rotn[8] = cos(alpha);
  
  return *this;
}

//! Calc rotn 3x3 rotation matrix around y by alpha.
RotnMatr& RotnMatr::YRotnMatr(double alpha)
{
  rotn[0] = cos(alpha);
  rotn[1] = 0.0;
  rotn[2] = sin(alpha);

  rotn[3] = 0.0;
  rotn[4] = 1.0;
  rotn[5] = 0.0;

  rotn[6] = -sin(alpha);
  rotn[7] = 0.0;
  rotn[8] = cos(alpha);

  return *this;
}

//! Calc rotn 3x3 rotation matrix around z by alpha.
RotnMatr& RotnMatr::ZRotnMatr(double alpha)
{
  rotn[0] = cos(alpha);
  rotn[1] = sin(alpha);
  rotn[2] = 0.0;

  rotn[3] = -sin(alpha);
  rotn[4] = cos(alpha);
  rotn[5] = 0.0;

  rotn[6] = 0.0;
  rotn[7] = 0.0;
  rotn[8] = 1.0;
  
  return *this;
}


//! Calc determinant of 3x3 rotation matrix.
double RotnMatr::DetRotnMatr() const 
{
  double detA;
  detA = rotn[0]*rotn[4]*rotn[8]+rotn[3]*rotn[7]*rotn[2]+rotn[6]*rotn[1]*rotn[5];
  detA -= rotn[2]*rotn[4]*rotn[6]+rotn[5]*rotn[7]*rotn[0]+rotn[8]*rotn[1]*rotn[3];
  
  return detA;
}

//! Calc inverse of rotation matrix a[9].
RotnMatr RotnMatr::GetInvMatr() const
{
  double detA = DetRotnMatr();
  detA = 1/detA;
  
  RotnMatr inva;

  inva[0] = detA*(rotn[4]*rotn[8]-rotn[5]*rotn[7]);
  inva[1] = detA*(rotn[2]*rotn[7]-rotn[1]*rotn[8]);
  inva[2] = detA*(rotn[1]*rotn[5]-rotn[2]*rotn[4]);

  inva[3] = detA*(rotn[5]*rotn[6]-rotn[3]*rotn[8]);
  inva[4] = detA*(rotn[0]*rotn[8]-rotn[2]*rotn[6]);
  inva[5] = detA*(rotn[2]*rotn[3]-rotn[0]*rotn[5]);

  inva[6] = detA*(rotn[3]*rotn[7]-rotn[4]*rotn[6]);
  inva[7] = detA*(rotn[1]*rotn[6]-rotn[0]*rotn[7]);
  inva[8] = detA*(rotn[0]*rotn[4]-rotn[1]*rotn[3]);
  
  return inva;
}


//! Convert x coord from molec to lab frame or the reverse for (inv) rotn matrix (inv)a[9].
double RotnMatr::XLabMol(double xL, double yL, double zL) const
{
  double x = rotn[0]*xL+rotn[1]*yL+rotn[2]*zL;
  return x;
}

//! Convert y coord from molec to lab frame or the reverse for (inv) rotn matrix (inv)a[9].
double RotnMatr::YLabMol(double xL, double yL, double zL) const
{
  double y = rotn[3]*xL+rotn[4]*yL+rotn[5]*zL;
  return y;
}

//! Convert z coord from molec to lab frame or the reverse for (inv) rotn matrix (inv)a[9].
double RotnMatr::ZLabMol(double xL, double yL, double zL) const
{
  double z = rotn[6]*xL+rotn[7]*yL+rotn[8]*zL;
  return z;
}


/*
//! Convert x,y,z from molec to labe frame or the reverse.
void RotnMatr::GenLabMol(double &xM, double &yM, double &zM, double xL, double yL, double zL) const
{
  xM = rotn[0]*xL+rotn[1]*yL+rotn[2]*zL;
  yM = rotn[3]*xL+rotn[4]*yL+rotn[5]*zL;
  zM = rotn[6]*xL+rotn[7]*yL+rotn[8]*zL;
}

//! Convert x,y,z from molec to labe frame or the reverse, for cyllindrical avg around X.
void RotnMatr::CylXLabMol(double &xM, double &yM, double &zM, double xL, double yL, double zL) const
{
  xM = xL;
  yM = rotn[4]*yL+rotn[5]*zL;
  zM = rotn[7]*yL+rotn[8]*zL;
}  

//! Convert x,y,z from molec to labe frame or the reverse, for cyllindrical avg around Y.
void RotnMatr::CylYLabMol(double &xM, double &yM, double &zM, double xL, double yL, double zL) const
{
  xM = rotn[0]*xL+rotn[2]*zL;
  yM = yL;
  zM = rotn[6]*xL+rotn[8]*zL;
}

//! Convert x,y,z from molec to labe frame or the reverse, for cyllindrical avg around Z.
void RotnMatr::CylZLabMol(double &xM, double &yM, double &zM, double xL, double yL, double zL) const
{
  xM = rotn[0]*xL+rotn[1]*yL;
  yM = rotn[3]*xL+rotn[4]*yL;
  zM = zL;
}
*/

//! Choose which function to use for cylindrical rotn of X,Y,Z axes.
/*void RotnMatr::CylLabMol(MolAvg molavg, double &xM, double &yM, double &zM, double xL, double yL, double zL) const
{
  switch(molavg)
    {
    case CYLX:
      CylXLabMol(xM,yM,zM,xL,yL,zL);
      break;
    case CYLY:
      CylYLabMol(xM,yM,zM,xL,yL,zL);
      break;
    case CYLZ:
      CylZLabMol(xM,yM,zM,xL,yL,zL);
      break;
    default:
      perror("RontMatr::CylLabMol(): invalid molavg");
      break;
    }
}
*/

//! Convert dipole moment in lab frame for a certain rotn matrix inva[9].
void RotnMatr::DipMomInLab(double *L_M, double *M) const 
{
  L_M[0] = XLabMol(M[0],M[1],M[2]);
  L_M[1] = YLabMol(M[0],M[1],M[2]);
  L_M[2] = ZLabMol(M[0],M[1],M[2]);
  L_M[3] = XLabMol(M[3],M[4],M[5]);
  L_M[4] = YLabMol(M[3],M[4],M[5]);
  L_M[5] = ZLabMol(M[3],M[4],M[5]);
}


/*! Form matrix of molec orientation distr before ionization, after polarized light excitation.
  Averaging function for all (alpha,beta,gamma) orientations(nang).*/
void RempiAvgFn(double *avgfn, double *Mex, double *Rex, RotnMatr *inva, int nang)
{
  double L_Mex[6];
  FILE *outfile = fopen("avgfn.chk","w");
  fprintf(outfile,"Rempi Avg Fn\n");
  fprintf(outfile,"i avgfn   LMex_x LMex_y LMex_z   LMex_x LMex_y LMex_z\n");
  for (int i=0; i<nang; i++)
    {
      inva[i].DipMomInLab(L_Mex,Mex);
      avgfn[i] = RempiProb(L_Mex,Rex);
      fprintf(outfile,"%10d   %+lE      %+lE   %+lE   %+lE",i,avgfn[i],L_Mex[0],L_Mex[1],L_Mex[2]);
      fprintf(outfile,"      %+lE   %+lE   %+lE\n",L_Mex[3],L_Mex[4],L_Mex[5]);
    }
  fclose(outfile);
}


/* ! Excitation probability for a molecule with orientation corresponding to dipole mom in lab frame,L_Mex[x,y,z],
  by excitation laser with polarization Rex[x,y,z]. */
double RempiProb(double *L_Mex, double *Rex)
{
  double prob_abc = L_Mex[0]*Rex[0]*L_Mex[3]*Rex[0];
  prob_abc += L_Mex[1]*Rex[1]*L_Mex[4]*Rex[1];
  prob_abc += L_Mex[2]*Rex[2]*L_Mex[5]*Rex[2];
  prob_abc /= 4*pow(M_PI,3);

  return prob_abc;
}

double CylAvgfunction()
{
  return 0.5/M_PI;
}
  
/* 
 1/4pi comes from the unit solid angle over a sphere. We average the value of coeffs over alpha and
beta that span a whole sphere, so we have to divide by the surface of the sphere. The area element is
actually d(alpha)*d(cos(beta)), and that is why we calc dbeta in isogrid, it is more complicated.
The averaging over gamma involves a cyllindrical averaging at each point on the sphere in case the
object does not have cyllindrical symmetry. For this averaging we should divide by 2pi for each gamma
point, and this is added in isogrid.DBeta() too. For the points where we do not consider gamma
averaging (alpha 0, 180) there is no 1/2pi factor. This is the only way to get the same number for 1s
orbital with isotropic or no averaging.
*/
double AvgFunction(int molavg, double alpha, double beta, double gamma)
{
  double weightavg;
  switch (molavg)
  {
    case NUM:
      weightavg = 0.25/M_PI;
      break;
/*    case REMPI:
      weightavg = 0.0;
      fprintf(stderr,"avg function: REMPI not implemented yet\n");
      break;
*/
  }
  return weightavg;
}
