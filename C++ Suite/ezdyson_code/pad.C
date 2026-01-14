#include "pad.h"
#include "ylm.h"

PAD::PAD(int nth, int nk) : ntheta(nth), nkv(nk) 
{
  pad=new double[nkv*ntheta];
  crosspad=new double[nkv*ntheta];
  totalpad=new double[nkv*ntheta];  
  memset(pad,0,nkv*ntheta*sizeof(double));
  memset(crosspad,0,nkv*ntheta*sizeof(double));
  memset(totalpad,0,nkv*ntheta*sizeof(double));
}

PAD::~PAD()
{
  delete [] pad;
  delete [] crosspad;
  delete [] totalpad;
}


//! Caution: this is unscaled PAD
// SG: PAD calculated for averaged Cklm
//AIK: CalcPad can be called several times, to aggregare contributions from multiple centers.
//Use coeffitient (scale_coeff, which is norm of DO*energy-independent constants) to scale 
void PAD::CalcPad(CklmCoeff& allCklm, double scale_coeff) const
{
  printf("\nCalculating PADs and cross sections\nntheta=%d nkv=%d scale_coeff=%lf\n", ntheta,nkv, scale_coeff); fflush(stdout);
  //pad and crosspad are aggregating the data, but totalpad is not
  
  double delta=DeltaTheta(), theta, tmp, tmp2, coeff;
  int lmax=allCklm.GetLMax();

  for (int k=0; k<nkv; k++) //LOOP BY k
    for (int t=0; t<ntheta; t++) //Loop by theta
      {
	theta = t*delta;
	//Loops by i=(l,m) 
	for (int i=0, l=0; l<=lmax; l++)
	  for (int m=-l; m<=l; m++, i++)
	    {
	      tmp = ThetaYlm(theta,l,m);
	      for (int l2=0; l2<=lmax; l2++)
		{
		  if (l2==l) //diagonal contributons
		    {
                      coeff = allCklm.GetCklmSq(k,i)*scale_coeff;
		      pad[k*ntheta+t] += coeff*tmp*tmp;//*costsq;
		    }
		  else //cross-terms
		    {
		      tmp2 = ThetaYlm(theta,l2,m);
                      coeff = allCklm.GetCrossCklm(k,i,l2)*scale_coeff;
		      crosspad[k*ntheta+t] += coeff*tmp*tmp2;
		    }
		}
	    }
	totalpad[k*ntheta+t] = (pad[k*ntheta+t] + crosspad[k*ntheta+t]);
      }
  printf("PAD: done\n");
}


#if 0  //This one without i^l factor  and with cross-cklm terms for coherences
void PAD::CalcMCPadWithCoherences(KLMPoints& klmgrid, CklmCoeff& allCklm1, CklmCoeff& allCklm2,
				  double *center1, double *center2,
				  double scale_coeff1, double scale_coeff2, double *xsec) const
{
  printf("\nCalculating MC PADs and cross sections with coherences\nntheta=%d nkv=%d scale_coeff1=%lf scale_coeff1=%lf\n", ntheta,nkv,
	 scale_coeff1, scale_coeff2); fflush(stdout);
  //pad and crosspad are aggregating the data, but totalpad is not
  
  double delta=DeltaTheta(), theta, tmp, tmp2, coeff;
  int lmax=allCklm1.GetLMax();

  double rab[XYZ];
  double dist;
  for(int i=0; i<XYZ; i++) {
      rab[i]=center1[i]-center2[i];
      dist+=rab[i]*rab[i];
    }
  dist=sqrt(dist);

  std::cout << "Distance between center1 and center2: " << dist << " x=" << rab[X] << " y=" << rab[Y] << " z=" << rab[Z] << std::endl;
  std::complex<double> ci(0.,1.);
  
  for (int k=0; k<nkv; k++) //LOOP BY k
    {
      double kwave=klmgrid.GetKV(k);
      std::complex<double> ci(0.,1.); //i
      
      for (int t=0; t<ntheta; t++) //Loop by theta
	{
	  theta = t*delta;
	  
	  //Only distance along Z matters:
	  std::complex<double> coherence(0.0, rab[Z]*cos(theta)*kwave);
	  double coh=std::real(exp(coherence));
	  //std::cout << "Theta=" << theta << " k=" << k << 
	  //" Coherence factor=" << coherence << " real(exp(coherence))=" << coh << std::endl;

	  //Loops by i=(l,m) 
	  for (int i=0, l=0; l<=lmax; l++)
	    for (int m=-l; m<=l; m++, i++)
	      {
		tmp = ThetaYlm(theta,l,m);
		for (int l2=0; l2<=lmax; l2++)
		  {
		    if (l2==l) //diagonal contributons
		      {
			//Need fabs to avoid NANs when numbers are too small!
			double prod1=sqrt(fabs(allCklm1.GetCklmSq(k,i)*scale_coeff1));
			double prod2=sqrt(fabs(allCklm2.GetCklmSq(k,i)*scale_coeff2));			
			pad[k*ntheta+t] += coh*prod1*prod2*tmp*tmp;
			//xsec[k]+=coh*prod1*prod2;
			xsec[k]+=coh*prod1*prod2*sin(theta)*tmp*tmp;
			//std::cout << "P1="<<prod1 <<" P2="<< prod2 << " pad=" <<   pad[k*ntheta+t] << std::endl;
		      }
		    else //cross-terms
		      {
			tmp2 = ThetaYlm(theta,l2,m);
			//Assume l on the bra, l' on the ket: phase=(i^l)^*(i^l2)... I am not sure what I am doing here
			double term1=allCklm1.GetCrossCklm(k,i,l2)*scale_coeff1;
			double term2=allCklm2.GetCrossCklm(k,i,l2)*scale_coeff2;
			double prod1=(term1>=0.0) ? sqrt(fabs(term1)) : -sqrt(fabs(term1));
			double prod2=(term2>=0.0) ? sqrt(fabs(term2)) : -sqrt(fabs(term2));
#if 1			
			crosspad[k*ntheta+t] += coh*prod1*prod2*tmp*tmp2; // no phase
			//with phase
#else  
			std::complex<double> phase_bra=std::conj(pow(ci,l)), phase_ket=pow(ci,l2);
			crosspad[k*ntheta+t] += std::real(phase_bra*phase_ket)*coh*prod1*prod2*tmp*tmp2;
			//This seems to give correct answer for simple MC test for total xsec, but beta become different
			//at larger Ekin. Need to look into it some more.
#endif			
			//std::cout << "P1="<<prod1 <<" P2="<< prod2 << " xpad=" <<   crosspad[k*ntheta+t] << std::endl;
		      }
		  }
	      }
	  totalpad[k*ntheta+t] = (pad[k*ntheta+t] + crosspad[k*ntheta+t]);
	}
      //this is equivalent to integrating over theta
      xsec[k]=xsec[k]*delta*2.0*M_PI;
    }
  printf("PAD: done\n");
}


#else  //Experimental
void PAD::CalcMCPadWithCoherences(KLMPoints& klmgrid, CklmCoeff& allCklm1, CklmCoeff& allCklm2,
				  double *center1, double *center2,
				  double scale_coeff1, double scale_coeff2, double *xsec) const
{
  printf("\nCalculating MC PADs and cross sections with coherences\nntheta=%d nkv=%d scale_coeff1=%lf scale_coeff1=%lf\n", ntheta,nkv,
	 scale_coeff1, scale_coeff2); fflush(stdout);
  //pad and crosspad are aggregating the data, but totalpad is not
  
  double delta=DeltaTheta(), theta, tmp, tmp2, coeff;
  int lmax=allCklm1.GetLMax();

  double rab[XYZ];
  double dist;
  for(int i=0; i<XYZ; i++) {
      rab[i]=center1[i]-center2[i];
      dist+=rab[i]*rab[i];
    }
  dist=sqrt(dist);

  std::cout << "Distance between center1 and center2: " << dist << " x=" << rab[X] << " y=" << rab[Y] << " z=" << rab[Z] << std::endl;
  std::complex<double> ci(0.,1.);
  
  for (int k=0; k<nkv; k++) //LOOP BY k
    {
      double kwave=klmgrid.GetKV(k);
      std::complex<double> ci(0.,1.); //i
      
      for (int t=0; t<ntheta; t++) //Loop by theta
	{
	  theta = t*delta;
	  
	  //Only distance along Z matters:
	  std::complex<double> coherence(0.0, rab[Z]*cos(theta)*kwave);
	  double coh=std::real(exp(coherence));
	  //std::cout << "Theta=" << theta << " k=" << k << 
	  //" Coherence factor=" << coherence << " real(exp(coherence))=" << coh << std::endl;

	  //Loops by i=(l,m) 
	  for (int i=0, l=0; l<=lmax; l++)
	    for (int m=-l; m<=l; m++, i++)
	      {
		tmp = ThetaYlm(theta,l,m);
		for (int l2=0; l2<=lmax; l2++)
		  {
		    if (l2==l) //diagonal contributons
		      {
			//Need fabs to avoid NANs when numbers are too small!
			double prod1=sqrt(fabs(allCklm1.GetCklmSq(k,i)*scale_coeff1));
			double prod2=sqrt(fabs(allCklm2.GetCklmSq(k,i)*scale_coeff2));			
			pad[k*ntheta+t] += coh*prod1*prod2*tmp*tmp;
			//xsec[k]+=coh*prod1*prod2;
			xsec[k]+=coh*prod1*prod2*sin(theta)*tmp*tmp;
			//std::cout << "P1="<<prod1 <<" P2="<< prod2 << " pad=" <<   pad[k*ntheta+t] << std::endl;
		      }
		    else //cross-terms
		      {
			tmp2 = ThetaYlm(theta,l2,m);
			//Assume l on the bra, l' on the ket: phase=(i^l)^*(i^l2)... I am not sure what I am doing here
			double term1=allCklm1.GetCrossCklm(k,i,l2)*scale_coeff1;
			double term2=allCklm2.GetCrossCklm(k,i,l2)*scale_coeff2;
			double prod1=(term1>=0.0) ? sqrt(fabs(term1)) : -sqrt(fabs(term1));
			double prod2=(term2>=0.0) ? sqrt(fabs(term2)) : -sqrt(fabs(term2));
			if (dist == 0.0 ) //add cross pad only to the same-senter
			  crosspad[k*ntheta+t] += coh*prod1*prod2*tmp*tmp2; // no phase
			//with phase
			/* do nothing for different center 
			   std::complex<double> phase_bra=std::conj(pow(ci,l)), phase_ket=pow(ci,l2);
			   crosspad[k*ntheta+t] += std::real(phase_bra*phase_ket)*coh*prod1*prod2*tmp*tmp2;  
			*/
		      }
		  }
	      }
	  totalpad[k*ntheta+t] = (pad[k*ntheta+t] + crosspad[k*ntheta+t]);
	}
      //this is equivalent to integrating over theta
      xsec[k]=xsec[k]*delta*2.0*M_PI;
    }
  printf("PAD: done\n");
}
#endif

void PAD::Print(const KLMPoints& klmgrid) const
{
  //AIK: need to rewrite appropriately
  FILE *outfile = fopen("pad.dat","w");
  double delta=DeltaTheta(), theta;
  
  for (int k=0,no=0; k<nkv; k++) //LOOP BY k
    {
      double kwave=klmgrid.GetKV(k);
      double ene=kwave*kwave/2.*HAR_TO_EV;
      fprintf(outfile,"Energy=%lf eV\n",ene); 
      fprintf(outfile,"theta,rad   sq_contrib(lm)   cross_contrib(lml’m)   total_PAD\n");
      for (int t=0; t<ntheta; t++,no++) //Loop by theta
	{
	  theta = t*delta;
	  fprintf(outfile," %13.6lf %17.6lE %17.6lE %17.6lE\n",theta,pad[no],crosspad[no],totalpad[no]);
	}
      fprintf(outfile,"__________________________________\n");
    }
  fclose(outfile);
}

//! Calculates total xsec by numeric integration, need it for consistency check
double PAD::TotalXSec(int kv) const
{
  double xsec=0.0, delta=DeltaTheta();

  int t=0;
  xsec+=totalpad[kv*ntheta+t]*sin(t*delta)*delta/2.;   
  
  for(t=1; t<ntheta-1; t++)
    xsec+=totalpad[kv*ntheta+t]*sin(t*delta)*delta; 
  
  xsec+=totalpad[kv*ntheta+t]*sin(t*delta)*delta/2.; 
  
  return xsec*2.*M_PI; //f integration 2PI;
}

//! returns parallel x-sec, i.e. theta=0
double PAD::XSec_par(int kv) const
{
  return totalpad[kv*ntheta+0]; 
}

//! returns perpendicular x-sec, i.e. theta=pi/2
double PAD::XSec_perp(int kv) const
{
  double delta=DeltaTheta();
  int midpoint=ntheta/2;
  double theta=midpoint*delta;

  //printf("PAD::XSec_Perp Ntheta=%d midpoint=%d value/(pi/2)=%lf",ntheta, midpoint,theta/(M_PI/2.));
  
  //AIK: not sure about 2PI
  if(fabs(theta/(M_PI/2.)-1.0)<0.00001)
    return totalpad[kv*ntheta+midpoint];
  else
    {
      printf("PAD::XSec_Perp() : Warning: The midpoint is not pi/2\n");
      return 0.5*(totalpad[kv*ntheta+midpoint]+totalpad[kv*ntheta+midpoint+1]);
    }
}
