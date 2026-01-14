#include "ylm.h"
#include "tools.h"

// SG: Note that this is NOT where ylm is computed for |k>
// Pl is not used anywhere else in the program
// ThetaYlm is used only to compute PADs
// Ylms computed for |k> are in sph.C

double Pl(double theta, int lv)
{
  double pl=1.0;
  
  if(0!=lv)
    {
      double cost = cos(theta);
  
      switch (lv)
	{
	case 1:
	  pl=cost;
	  break;
	case 2:
	  pl=0.5*(3.0*cost*cost-1.0);
	  break;
	case 3:
	  pl=0.5*(5.0*pow(cost,3)-3.0*cost);
	  break;
	case 4:
	  pl=1./8.*(35.*pow(cost,4)-30.*pow(cost,2)+ 3.);
	  break;
	case 5:
	  pl=1./8.*(63.*pow(cost,5)-70.*pow(cost,3)+15.*cost);
	  break;
	default:
	  pl=0.0; //AIK! Place warning
	  break;
	}
    }
  return pl;
}

double ThetaYlm(double theta, int lv, int mv)
{
  double epsint = sin(theta);
  double cost = cos(theta);


  double ylm;
  switch (lv)
    {
    case 0:
      if (mv==0)
        ylm = 0.5*sqrt(1.0/pi);
      break;
    case 1:
	if (mv==0)
	  ylm = 0.5*sqrt(3.0/pi)*cost;
	if (mv==1 || mv==-1)
	  ylm = double(mv)*(-1.0)*0.5*sqrt(1.5/pi)*epsint;
	break;
    case 2:
      if (mv==0)
        ylm = 0.25*sqrt(5.0/pi)*(3.0*pow(cost,2.0)-1.0);
      if (mv==1 || mv==-1)
        ylm = double(mv)*(-1.0)*0.5*sqrt(7.5/pi)*cost*epsint;
      if (mv==2 || mv==-2)
        ylm = 0.25*sqrt(7.5/pi)*pow(epsint,2);
      break;
    case 3:
      if (mv==0)
        ylm = 0.25*sqrt(7.0/pi)*(5.0*pow(cost,3.0)-3.0*cost);
      if (mv==1 || mv==-1)
        ylm = double(mv)*(-1.0)*0.125*sqrt(21.0/pi)*(5.0*pow(cost,2.0)-1.0)*epsint;
      if (mv ==2 || mv==-2)
        ylm = 0.25*sqrt(52.5/pi)*cost*pow(epsint,2);
      if (mv==3 || mv==-3)
        ylm = (double(mv)/3.0)*(-1.0)*0.125*sqrt(35.0/pi)*pow(epsint,3);
    break;
    case 4:
      if (mv==0)
        ylm = (3.0/16.0)*sqrt(1.0/pi)*(35.0*pow(cost,4.0)-30.0*pow(cost,2.0)+3.0);
      if (mv==1 || mv==-1)
        ylm = double(mv)*(-1.0)*(3.0/8.0)*sqrt(5.0/pi)*(7.0*pow(cost,3.0)-3.0*cost)*epsint;
      if (mv==2 || mv==-2)
        ylm = (3.0/8.0)*sqrt(2.5/pi)*(7.0*pow(cost,2.0)-1.0)*pow(epsint,2);
      if (mv==3 || mv==-3)
        ylm = (double(mv)/3.0)*(-1.0)*(3.0/8.0)*sqrt(35.0/pi)*cost*pow(epsint,3);
      if (mv==4 || mv==-4)
        ylm = (3.0/16.0)*sqrt(17.5/pi)*pow(epsint,4);
    break;
    case 5:
      if (mv==0)
        ylm = (1.0/16.0)*sqrt(11.0/pi)*(63.0*pow(cost,5.0)-70.0*pow(cost,3.0)+15.0*cost);
      if (mv==1 || mv==-1)
        ylm = double(mv)*(-1.0)*(1.0/16.0)*sqrt(82.5/pi)*(21.0*pow(cost,4.0)-14.0*pow(cost,2.0)+1.0)*epsint;
      if (mv==2 || mv==-2)
        ylm = 0.125*sqrt(577.5/pi)*(3.0*pow(cost,3.0)-cost)*pow(epsint,2);
      if (mv==3 || mv==-3)
        ylm = (double(mv)/3.0)*(-1.0)*(1.0/32.0)*sqrt(385.0/pi)*(9.0*pow(cost,2.0)-1.0)*pow(epsint,3);
      if (mv==4 || mv==-4)
        ylm = (3.0/16.0)*sqrt(192.5/pi)*cost*pow(epsint,4);
      if (mv==5 || mv==-5)
        ylm = (double(mv)/5.0)*(-1.0)*(3.0/32.0)*sqrt(77.0/pi)*pow(epsint,5);
    break;
    case 6:
      if (mv==0)
        ylm = (1.0/32.0)*sqrt(13.0/pi)*(231.0*pow(cost,6.0)-315.0*pow(cost,4.0)+105.0*pow(cost,2.0)-5.0);
      if (mv==1 || mv==-1)
        ylm = double(mv)*(-1.0)*(1.0/16.0)*sqrt(136.5/pi)*(33.0*pow(cost,5.0)-30.0*pow(cost,3.0)+5.0*cost)*epsint;
      if (mv==2 || mv==-2)
        ylm = (1.0/64.0)*sqrt(1365.0/pi)*(33.0*pow(cost,4.0)-18.0*pow(cost,2.0)+1.0)*pow(epsint,2);
      if (mv==3 || mv==-3)
        ylm = (double(mv)/3.0)*(-1.0)*(1.0/32.0)*sqrt(1365.0/pi)*(11.0*pow(cost,3.0)-3.0*cost)*pow(epsint,3);
      if (mv==4 || mv==-4)
        ylm = (3.0/32.0)*sqrt(45.5/pi)*(11.0*pow(cost,2.0)-1.0)*pow(epsint,4);
      if (mv==5 || mv==-5)
        ylm = (double(mv)/5.0)*(-1.0)*(3.0/32.0)*sqrt(1001.0/pi)*cost*pow(epsint,5);
      if (mv==6 || mv==-6)
        ylm = (1.0/64.0)*sqrt(3003.0/pi)*pow(epsint,6);
    break;
    case 7:
      if (mv==0)
        ylm = (1.0/32.0)*sqrt(15.0/pi)*(429.0*pow(cost,7.0)-693.0*pow(cost,5.0)+315.0*pow(cost,3.0)-35.0*cost);
      if (mv==1 || mv==-1)
        ylm = double(mv)*(-1.0)*(1.0/64.0)*sqrt(52.5/pi)*(429.0*pow(cost,6.0)-495.0*pow(cost,4.0)+135.0*pow(cost,2.0)-5.0)*epsint;
      if (mv==2 || mv==-2)
        ylm = (3.0/64.0)*sqrt(35.0/pi)*(143.0*pow(cost,5.0)-110.0*pow(cost,3.0)+15.0*cost)*pow(epsint,2);
      if (mv==3 || mv==-3)
        ylm = (double(mv)/3.0)*(-1.0)*(3.0/64.0)*sqrt(17.5/pi)*(143.0*pow(cost,4.0)-66.0*pow(cost,2.0)+3.0)*pow(epsint,3);
      if (mv==4 || mv==-4)
        ylm = (3.0/32.0)*sqrt(192.5/pi)*(13.0*pow(cost,3.0)-3.0*cost)*pow(epsint,4);
      if (mv==5 || mv==-5)
        ylm = (double(mv)/5.0)*(-1.0)*(3.0/64.0)*sqrt(192.5/pi)*(13.0*pow(cost,2.0)-1.0)*pow(epsint,5);
      if (mv==6 || mv==-6)
        ylm = (3.0/64.0)*sqrt(5005.0/pi)*cost*pow(epsint,6);
      if (mv==7 || mv==-7)
        ylm = (double(mv)/7.0)*(-1.0)*(3.0/64.0)*sqrt(357.5/pi)*pow(epsint,7);
      break;
    case 8:
      if (mv==0)
        ylm = (1.0/256.0)*sqrt(17.0/pi)*(6435.0*pow(cost,8.0)-12012.0*pow(cost,6.0)+6930.0*pow(cost,4.0)-1260.0*pow(cost,2.0)+35.0);
      if (mv==1 || mv==-1)
        ylm = double(mv)*(-1.0)*(3.0/64.0)*sqrt(8.5/pi)*(715.0*pow(cost,7.0)-1001.0*pow(cost,5.0)+385.0*pow(cost,3.0)-35.0*cost)*epsint;
      if (mv==2 || mv==-2)
        ylm = (3.0/128.0)*sqrt(595.0/pi)*(143.0*pow(cost,6.0)-143.0*pow(cost,4.0)+33.0*pow(cost,2.0)-1.0)*pow(epsint,2);
      if (mv==3 || mv==-3)
        ylm = (double(mv)/3.0)*(-1.0)*(1.0/64.0)*sqrt(9817.5/pi)*(39.0*pow(cost,5.0)-26.0*pow(cost,3.0)+3.0*cost)*pow(epsint,3);
      if (mv==4 || mv==-4)
        ylm = (3.0/128.0)*sqrt(654.5/pi)*(65.0*pow(cost,4.0)-26.0*pow(cost,2.0)+1.0)*pow(epsint,4);
      if (mv==5 || mv==-5)
        ylm = (double(mv)/5.0)*(-1.0)*(3.0/64.0)*sqrt(8508.5/pi)*(5.0*pow(cost,3.0)-cost)*pow(epsint,5);
      if (mv==6 || mv==-6)
        ylm = (1.0/128.0)*sqrt(7293.0/pi)*(15.0*pow(cost,2.0)-1.0)*pow(epsint,6);
      if (mv==7 || mv==-7)
        ylm = (double(mv)/7.0)*(-1.0)*(3.0/64.0)*sqrt(6077.5/pi)*cost*pow(epsint,7);
      if (mv==8 || mv==-8)
        ylm = (3.0/256.0)*sqrt(6077.5/pi)*pow(epsint,8);
      break;
    case 9:
      if (mv==0)
        ylm = (1.0/256.0)*sqrt(19.0/pi)*(12155.0*pow(cost,9.0)-25740.0*pow(cost,7.0)+18018.0*pow(cost,5.0)-4620.0*pow(cost,3.0)+315.0*cost);
      if (mv==1 || mv==-1)
        ylm = double(mv)*(-1.0)*(3.0/256.0)*sqrt(47.5/pi)*(2431.0*pow(cost,8.0)-4004.0*pow(cost,6.0)+2002.0*pow(cost,4.0)-308.0*pow(cost,2.0)+7.0)*epsint;
      if (mv==2 || mv==-2)
        ylm = (3.0/128.0)*sqrt(1045.0/pi)*(221.0*pow(cost,7.0)-273.0*pow(cost,5.0)+91.0*pow(cost,3.0)-7.0*cost)*pow(epsint,2);
      if (mv==3 || mv==-3)
        ylm = (double(mv)/3.0)*(-1.0)*(1.0/256.0)*sqrt(21945.0/pi)*(221.0*pow(cost,6.0)-195.0*pow(cost,4.0)+39.0*pow(cost,2.0)-1.0)*pow(epsint,3);
      if (mv==4 || mv==-4)
        ylm = (3.0/128.0)*sqrt(47547.5/pi)*(17.0*pow(cost,5.0)-10.0*pow(cost,3.0)+cost)*pow(epsint,4);
      if (mv==5 || mv==-5)
        ylm = (double(mv)/5.0)*(-1.0)*(3.0/256.0)*sqrt(2717.0/pi)*(85.0*pow(cost,4.0)-30.0*pow(cost,2.0)+1.0)*pow(epsint,5);
      if (mv==6 || mv==-6)
        ylm = (1.0/128.0)*sqrt(40755.0/pi)*(17.0*pow(cost,3.0)-3.0*cost)*pow(epsint,6);
      if (mv==7 || mv==-7)
        ylm = (double(mv)/7.0)*(-1.0)*(3.0/512.0)*sqrt(13585.0/pi)*(17.0*pow(cost,2.0)-1.0)*pow(epsint,7);
      if (mv==8 || mv==-8)
        ylm = (3.0/256.0)*sqrt(115472.5/pi)*cost*pow(epsint,8);
      if (mv==9 || mv==-9)
        ylm = (double(mv)/9.0)*(-1.0)*(1.0/512.0)*sqrt(230945.0/pi)*pow(epsint,9);
    break;
    case 10:
      if (mv==0)
        ylm = (1.0/512.0)*sqrt(21.0/pi)*(46189.0*pow(cost,10.0)-109395.0*pow(cost,8.0)+90090.0*pow(cost,6.0)-30030.0*pow(cost,4.0)+3465.0*pow(cost,2.0)-63.0);
      if (mv==1 || mv==-1)
        ylm = double(mv)*(-1.0)*(1.0/256.0)*sqrt(577.5/pi)*(4199.0*pow(cost,9.0)-7956.0*pow(cost,7.0)+4914.0*pow(cost,5.0)-1092.0*pow(cost,3.0)+63.0*cost)*epsint;
      if (mv==2 || mv==-2)
        ylm = (3.0/512.0)*sqrt(192.5/pi)*(4199.0*pow(cost,8.0)-6188.0*pow(cost,6.0)+2730.0*pow(cost,4.0)-364.0*pow(cost,2.0)+7.0)*pow(epsint,2);
      if (mv==3 || mv==-3)
        ylm = (double(mv)/3.0)*(-1.0)*(3.0/256.0)*sqrt(5005.0/pi)*(323.0*pow(cost,7.0)-357.0*pow(cost,5.0)+105.0*pow(cost,3.0)-7.0*cost)*pow(epsint,3);
      if (mv==4 || mv==-4)
        ylm = (3.0/256.0)*sqrt(2502.5/pi)*(323.0*pow(cost,6.0)-255.0*pow(cost,4.0)+45.0*pow(cost,2.0)-1)*pow(epsint,4);
      if (mv==5 || mv==-5)
        ylm = (double(mv)/5.0)*(-1.0)*(3.0/256.0)*sqrt(1001.0/pi)*(323.0*pow(cost,5.0)-170.0*pow(cost,3.0)+15.0*cost)*pow(epsint,5);
      if (mv==6 || mv==-6)
        ylm = (3.0/1024.0)*sqrt(5005.0/pi)*(323.0*pow(cost,4.0)-102.0*pow(cost,2.0)+3.0)*pow(epsint,6);
      if (mv==7 || mv==-7)
        ylm = (double(mv)/7.0)*(-1.0)*(3.0/512.0)*sqrt(85085.0/pi)*(19.0*pow(cost,3.0)-3.0*cost)*pow(epsint,7);
      if (mv==8 || mv==-8)
        ylm = (1.0/512.0)*sqrt(127627.5/pi)*(19.0*pow(cost,2.0)-1.0)*pow(epsint,8);
      if (mv==9 || mv==-9)
        ylm = (double(mv)/9.0)*(-1.0)*(1.0/512.0)*sqrt(4849845.0/pi)*cost*pow(epsint,9);
      if (mv==10 || mv==-10)
        ylm = (1.0/1024.0)*sqrt(969969.0/pi)*pow(epsint,10);
      break;
    }

  return ylm;
}






