#include "wavetypes.h"

//!Calc Spherical Rkl for kr<<l.
double LimRkl(double kr, int l)
{
  int num = 2*l+1;
  double rkl = pow(kr,double(l))/dfac(num);
  return rkl;
}


//! Calc by recurrence all Spherical Rkl's at one x,y,z point.
//AK!! lmax is bad, need to change such that the loops run accordingly, i.e., < lmax, not <=lmax
// SG: I think lmax is fine?

void SphericalRkl(double *rklmax, double x, double y, double z, double kv, int lmax)
{
  double rv = sqrt(x*x+y*y+z*z);
  double kr = kv*rv;
  
  //max l for which to use recurrence; for l > lim, use LimRkl function
  int lim;
  if (kr > 1.4)
     lim = 10;
  else if (kr > 1.2)
     lim = 9;
  else if (kr > 8.4)
     lim = 8;
  else if (kr > 0.55)
     lim = 7;
  else if (kr > 0.285)
     lim = 6;
  else if (kr > 0.14)
     lim = 5;
  else if (kr > 0.06)
     lim =4;
  else if (kr > 0.013)
     lim =3;
  else
     lim =2;
  lim = ((lim <= lmax) ? lim : lmax);
  
  static double thresh=1E-6;

  if(kr > thresh)
    {
      rklmax[0] = Rk0(kr);
      rklmax[1] = Rk1(kr);
      for (int j=2; j<=lim; j++)
	rklmax[j] = (2*j-1)*rklmax[j-1]/kr-rklmax[j-2];

      for (int j=lim+1; j<=lmax; j++)
	rklmax[j] = LimRkl(kr,j);
    }
  else
     for (int j=0; j<=lmax; j++)
         rklmax[j] = LimRkl(kr,j);
}

/*! Calc all Coulomb waves at one x,y,z point.
Adapted from web PL/1 code: http://www.plasmaphysics.org.uk/programs/coulombc.htm  */
void CoulombRkl(double *rklmax, double x, double y, double z, double kv, int lmax, double radfn)
{
  double rv = sqrt(x*x+y*y+z*z);           //calc r, ro, eta

  static int if1sttime=1;
  static double QG[16];              //weights (i think) for integration
  if(if1sttime)
    {
      QG[0]=9.894009349916499325961541735E-01;  QG[1]=2.71524594117540948517805725E-02;
      QG[2]=9.445750230732325760779884155E-01;  QG[3]=6.22535239386478928628438370E-02;
      QG[4]=8.656312023878317438804678977E-01;  QG[5]=9.51585116824927848099251076E-02;
      QG[6]=7.554044083550030338951011948E-01;  QG[7]=1.246289712555338720524762822E-01;
      QG[8]=6.178762444026437484466717640E-01;  QG[9]=1.495959888165767320815017305E-01; 
      QG[10]=4.580167776572273863424194430E-01; QG[11]=1.691565193950025381893120790E-01;
      QG[12]=2.816035507792589132304605015E-01; QG[13]=1.826034150449235888667636680E-01; 
      QG[14]=9.50125098376374401853193354E-02;  QG[15]=1.894506104550684962853967232E-01;
      
      for (int i=0; i<16; i++)
	QG[i] *= 0.5;

      if1sttime=0;
    }

  double ro = kv*rv;
  double z_ion = radfn;
  double z_el = -1;
  double eta = z_ion*z_el/kv;

  int il, l, l1;
  double FLN, FLL=0.0, FLU=0.0;
  double ULIMO = 1E28*kv;
  for (il=lmax; il>=lmax-1; il--)            //calc Rkl for two highest l first
   {
     l = il;
     l1=l+1;                 
     Complex LE(double(l1),eta);                                                           //argument of Gamma function   
     double KVG = 1E-7, XDG = 5.0;
     Complex GLE = Gamma(QG,KVG,XDG,LE);  
     double GAE = exp(-M_PI*eta/2)/(pow(2.0,double(l))*sqrt(GLE.Norm()));                  //calc prefactor

     // Calc integral    
     double XD00 = 0.5;
     double XDL = XD00/sqrt(l1);                          
     double XL, XU=0.0, XD0 = XDL, XD, FL = 0.0, QY = 0.0, QYL, QYU = 0.0, QYM;
     int jj = 1, di=1;
     double KV = 1E-7, CONVGC = 1.0;

     while (CONVGC >= KV)                               //Numerical integration step at one k, r, l
      {
        QYL = QY;
        XL = XU;  
        XD = MIN(XD0,M_PI/(ro/pow(cosh(XL),2.0)+fabs(2*eta)));
        if (0.0==XL) 
           XD *= 0.5;
        XU = XL + XD;
        QY = DQG16(QG,ro,eta,XL,XU,l1);
        QYU += fabs(QY);
        FL += QY;
        QYM = fabs(FL)/jj;

        if (QYM > 0.0)
           CONVGC = (fabs(QY)+fabs(QYL))/(2*QYM);
        jj++;
      }
     if (1==di)
      {
        di = -1;
        CONVGC = 1.0;
        FL -= QY;
        jj--;
        XU -= XD;
        QYU -= fabs(QY);
        XD0 = XDL/2;

        //Repeat last integration step with half stepwidth
        while (CONVGC >= KV)                             
         {
           QYL = QY;
           XL = XU;
           XD = MIN(XD0,0.5*M_PI/(ro/pow(cosh(XL),2.0)+fabs(2*eta)));
           if (0==XL)
              XD *= 0.5;
           XU = XL + XD;
           QY = DQG16(QG,ro,eta,XL,XU,l1);
           QYU += fabs(QY);
           FL += QY;
           QYM = fabs(FL)/jj;
           if (QYM > 0.0)
              CONVGC = (fabs(QY)+fabs(QYL))/(2*QYM);
           jj++;
         }
      }
 
     double CFL = GAE*pow(ro,double(l1));
     FL *= CFL;
     if ((QYU/2*CFL > ULIMO) || ((QYU/2*CFL/FL > 1E28)))
        fprintf(stderr,"CoulWave(l=%d): Desired accuracy cannot be reached\n",il);
     if (lmax==il) 
        FLU = FL;
     else
        FLL = FL;
     rklmax[il] = FL/ro;
   }

    FLN = FLL;
    FLL = FLU;
    double FLN1, FLN2;
    for (il=lmax-1; il>=1; il--)
     {
       l = il;
       l1 = l+1;
       FLU = FLL;
       FLL = FLN;
       FLN1 = (2*l+1)*(eta+l*l1/ro)*FLL;
       FLN2 = l*sqrt(pow(double(l1),2.0)+pow(eta,2.0))*FLU;
       if ((fabs(FLN1)+fabs(FLN2))/fabs(FLN1-FLN2) > ULIMO)
          fprintf(stderr,"CoulWave(l=%d): Desired accuracy cannot be reached\n",il-1);
       FLN = (FLN1-FLN2)/(l1*sqrt(pow(double(l),2.0)+pow(eta,2.0)));
       rklmax[il-1] = FLN/ro;
     }
}

//! Procedure for calc Coulomb waves (see PL/1 code)
double DQG16(double *QG, double ro, double eta, double XL, double XU, int l1)
{
  double QY;
  double A = 0.5*(XU+XL);
  double B = XU-XL;
  double LY =0.0;
  double C;
  for (int i=0; i<15; i+=2)
   {
     C = QG[i]*B;
     LY += QG[i+1]*(FCT(A+C,ro,eta,l1)+FCT(A-C,ro,eta,l1));
   }
  QY = LY*B;

  return QY;
}

//! Procedure for calc of Coulomb waves (see PL/1 code)
double FCT(double x, double ro, double eta, int l1)
{
  double y;
  double chx = cosh(x);
  double chx2 = pow(chx,2.0);
  y = cos(ro*sqrt(chx2-1)/chx-2*eta*x)/pow(chx2,double(l1));

  return y;
}

//! Calc Gamma function for Coulomb wave prefactor.
Complex Gamma(double *QG, double KVG, double XDG, const Complex& Z)
{
  Complex G;
  double GR, GI;
  double RZ = Z.GetRe(), IZ = fabs(Z.GetIm());
  double EIZ = (IZ > 0.0) ? exp(0.5*M_PI/IZ) : 1E10/KVG;

  double XL0 = MAX(RZ-1.0,1.0);
  double XL, XU = XL0;
  XDG = XDG/(floor(RZ/32.0)+1);

  Complex QY;
  double QYR, QYI, QYMR, QYMI;
  double CONVGCR, CONVGCI, CONVGC = 1.0;
  int jj = 1, di = 1;
  while (CONVGC >= KVG)
   {
     if (1==di)
      {
        XL = XU;
        XU = XL + MIN(XL*EIZ-XL,XDG);
      }
     else
      {
        XU = XL;
        XL = MAX(XU/EIZ,XU-XDG);
      }
     QY = DQG16G(QG,Z,XL,XU);
     QYR = QY.GetRe();
     QYI = QY.GetIm();
     G += QY;
     GR = G.GetRe(); 
     GI = G.GetIm(); 
     QYMR = fabs(GR)/jj;
     QYMI = fabs(GI)/jj;
     if (QYMR > 0.0)
        CONVGCR = fabs(QYR/QYMR);
     if (QYMI > 0.0)
        CONVGCI = fabs(QYI/QYMI);
     CONVGC = MAX(CONVGCR,CONVGCI);
     jj++;
   }

  if (1==di)
   {
     CONVGC = 1.0;
     di = -1;
     XL = XL0;
     while (CONVGC >= KVG) //additional test step after convergence is achieved
      {
        if (1==di)
         {
           XL = XU;
           XU = XL + MIN(XL*EIZ-XL,XDG);
         }
        else
         {
           XU = XL;
           XL = MAX(XU/EIZ,XU-XDG);
         }
        QY = DQG16G(QG,Z,XL,XU);
        QYR = QY.GetRe();
        QYI = QY.GetIm();
        G += QY;
        GR = G.GetRe(); 
        GI = G.GetIm(); 
        QYMR = fabs(GR)/jj;
        QYMI = fabs(GI)/jj;
        if (QYMR > 0.0)
           CONVGCR = fabs(QYR/QYMR);
        if (QYMI > 0.0)
           CONVGCI = fabs(QYI/QYMI);
        CONVGC = MAX(CONVGCR,CONVGCI);
        jj++;
      }
   }

  return G;
}

//! Procedure for calc of complex Gamma function (see PL/1 code)
Complex DQG16G(double *QG, const Complex& Z, double XL, double XU)
{
  Complex LY;
  double A = 0.5*(XU+XL);
  double B = XU - XL;
  double C;
 
  for (int i=0; i<15; i+=2)
   {
     C = QG[i]*B;
     LY += QG[i+1]*(FCTG(A+C,Z)+FCTG(A-C,Z));
   }
  
  return LY*B;
}

//! Calc all radial part functions Rkl's at one x,y,z point.

void Rkl(double radfn, double *rklmax, double x, double y, double z, double kv, int lmax)
{
  if (radfn==0)
      SphericalRkl(rklmax,x,y,z,kv,lmax);
  else if (x == 0 && y==0 && z==0)
// if r=0, the coulomb wave crashes. This is to circumvent this from happening, but it may be an approximation
// Need to find a more elegant solution
      for (int j=0; j<=lmax; j++)
        rklmax[j] = 0;
  else if (kv < 0.085)
// Currently, code is numerically unstable and inefficient at low k. Can either use gsl library or fix code..
      for (int j=0; j<=lmax; j++)
        rklmax[j] = 0;
  else
      CoulombRkl(rklmax,x,y,z,kv,lmax,radfn);
}


