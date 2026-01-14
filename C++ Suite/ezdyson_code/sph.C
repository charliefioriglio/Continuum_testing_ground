#include "sph.h"


SPH::SPH(int lm)
{
  lmax=lm;
  thetacoeff=new double *[lmax];
  for(int lv=0; lv<lmax; lv++)
    thetacoeff[lv]=new double[lv+1];

 
  thetacoeff[0][0]= 0.5*sqrt(1.0/pi);
  thetacoeff[1][0]= 0.5*sqrt(3.0/pi);
  thetacoeff[1][1]=(-1.0)*0.5*sqrt(1.5/pi);
  

  thetacoeff[2][0]=0.25*sqrt(5.0/pi);
  thetacoeff[2][1]=(-1.0)*0.5*sqrt(7.5/pi);
  thetacoeff[2][2]=0.25*sqrt(7.5/pi);

  thetacoeff[3][0]=0.25*sqrt(7.0/pi);
  thetacoeff[3][1]=(-1.0)*0.125*sqrt(21.0/pi);
  thetacoeff[3][2]=0.25*sqrt(52.5/pi);
  thetacoeff[3][3]=(-1.0)*0.125*sqrt(35.0/pi);

  thetacoeff[4][0]=(3.0/16.0)*sqrt(1.0/pi);
  thetacoeff[4][1]=(-1.0)*(3.0/8.0)*sqrt(5.0/pi);
  thetacoeff[4][2]=(3.0/8.0)*sqrt(2.5/pi);
  thetacoeff[4][3]=(-1.0)*(3.0/8.0)*sqrt(35.0/pi);
  thetacoeff[4][4]=(3.0/16.0)*sqrt(17.5/pi);

  thetacoeff[5][0]=(1.0/16.0)*sqrt(11.0/pi);
  thetacoeff[5][1]=(-1.0)*(1.0/16.0)*sqrt(82.5/pi);
  thetacoeff[5][2]=0.125*sqrt(577.5/pi);
  thetacoeff[5][3]=(-1.0)*(1.0/32.0)*sqrt(385.0/pi);
  thetacoeff[5][4]=(3.0/16.0)*sqrt(192.5/pi);
  thetacoeff[5][5]=(-1.0)*(3.0/32.0)*sqrt(77.0/pi);

  thetacoeff[6][0]=(1.0/32.0)*sqrt(13.0/pi);
  thetacoeff[6][1]=(-1.0)*(1.0/16.0)*sqrt(136.5/pi);
  thetacoeff[6][2]=(1.0/64.0)*sqrt(1365.0/pi);
  thetacoeff[6][3]=(-1.0)*(1.0/32.0)*sqrt(1365.0/pi);
  thetacoeff[6][4]=(3.0/32.0)*sqrt(45.5/pi);
  thetacoeff[6][5]=(-1.0)*(3.0/32.0)*sqrt(1001.0/pi);
  thetacoeff[6][6]=(1.0/64.0)*sqrt(3003.0/pi);

  thetacoeff[7][0]=(1.0/32.0)*sqrt(15.0/pi);
  thetacoeff[7][1]=(-1.0)*(1.0/64.0)*sqrt(52.5/pi);
  thetacoeff[7][2]=(3.0/64.0)*sqrt(35.0/pi);
  thetacoeff[7][3]=(-1.0)*(3.0/64.0)*sqrt(17.5/pi);
  thetacoeff[7][4]=(3.0/32.0)*sqrt(192.5/pi);
  thetacoeff[7][5]=(-1.0)*(3.0/64.0)*sqrt(192.5/pi);
  thetacoeff[7][6]=(3.0/64.0)*sqrt(5005.0/pi);
  thetacoeff[7][7]=(-1.0)*(3.0/64.0)*sqrt(357.5/pi);

  thetacoeff[8][0]=(1.0/256.0)*sqrt(17.0/pi);
  thetacoeff[8][1]=(-1.0)*(3.0/64.0)*sqrt(8.5/pi);
  thetacoeff[8][2]=(3.0/128.0)*sqrt(595.0/pi);
  thetacoeff[8][3]=(-1.0)*(1.0/64.0)*sqrt(9817.5/pi);
  thetacoeff[8][4]=(3.0/128.0)*sqrt(654.5/pi);
  thetacoeff[8][5]=(-1.0)*(3.0/64.0)*sqrt(8508.5/pi);
  thetacoeff[8][6]=(1.0/128.0)*sqrt(7293.0/pi);
  thetacoeff[8][7]=(-1.0)*(3.0/64.0)*sqrt(6077.5/pi);
  thetacoeff[8][8]=(3.0/256.0)*sqrt(6077.5/pi);

  thetacoeff[9][0]=(1.0/256.0)*sqrt(19.0/pi);
  thetacoeff[9][1]=(-1.0)*(3.0/256.0)*sqrt(47.5/pi);
  thetacoeff[9][2]=(3.0/128.0)*sqrt(1045.0/pi);
  thetacoeff[9][3]=(-1.0)*(1.0/256.0)*sqrt(21945.0/pi);
  thetacoeff[9][4]=(3.0/128.0)*sqrt(47547.5/pi);
  thetacoeff[9][5]=(-1.0)*(3.0/256.0)*sqrt(2717.0/pi);
  thetacoeff[9][6]=(1.0/128.0)*sqrt(40755.0/pi);
  thetacoeff[9][7]=(-1.0)*(3.0/512.0)*sqrt(13585.0/pi);
  thetacoeff[9][8]=(3.0/256.0)*sqrt(115472.5/pi);
  thetacoeff[9][9]=(-1.0)*(1.0/512.0)*sqrt(230945.0/pi);


  thetacoeff[10][0]=(1.0/512.0)*sqrt(21.0/pi);
  thetacoeff[10][1]=(-1.0)*(1.0/256.0)*sqrt(577.5/pi);
  thetacoeff[10][2]=(3.0/512.0)*sqrt(192.5/pi);
  thetacoeff[10][3]=-1.0*(3.0/256.0)*sqrt(5005.0/pi);
  thetacoeff[10][4]=(3.0/256.0)*sqrt(2502.5/pi);
  thetacoeff[10][5]=-1.0*(3.0/256.0)*sqrt(1001.0/pi);
  thetacoeff[10][6]= (3.0/1024.0)*sqrt(5005.0/pi);
  thetacoeff[10][7]=-1.0*(3.0/512.0)*sqrt(85085.0/pi);
  thetacoeff[10][8]=(1.0/512.0)*sqrt(127627.5/pi);
  thetacoeff[10][9]= -1.0*(1.0/512.0)*sqrt(4849845.0/pi);
  thetacoeff[10][10]=(1.0/1024.0)*sqrt(969969.0/pi);

}

SPH::~SPH()
{
  for(int lv=0; lv<lmax; lv++)
    delete [] thetacoeff[lv];
  
  delete [] thetacoeff;
}


Complex SPH::GetValue(double x, double y, double z, int lv, int mv) const
{
  double rv = sqrt(x*x+y*y+z*z);
  double re=0., im=0., cost=1.;
  
  if(fabs(rv)>0.000001)
    {
      cost = z/rv;
// SG: re and im are used later to compute sin(theta) * e^(m* i * phi)
      re = x/rv;
      im = (mv<0) ? -y/rv : y/rv;
    }
 
  //theta polynomial 
  double tpol=1.0; //if mv==lv tpol=1.0
  int abs_mv=(mv>=0) ? mv : -mv;
  
  /*  if |mv|==l, tpol=1, if |mv==l-1, tpol=cost,
      for all the rest use the following switch: */
  
  if(abs_mv!=lv) 
    {
      if(abs_mv==lv-1)
	tpol=cost;
      else
	{
	  switch (lv)
	    {
	    case 2:
	      tpol = 3.0*pow(cost,2)-1.0;
	      break;
	    case 3:
	      tpol = (!abs_mv) ?  5.0*pow(cost,3)-3.0*cost :  5.0*pow(cost,2)-1.0;
	      break;
	    case 4:
	      switch (abs_mv)
		{
		case 0:
		  tpol = 35.0*pow(cost,4)-30.0*pow(cost,2)+3.0;
		  break;
		case 1:
		  tpol = 7.0*pow(cost,3)-3.0*cost;
		  break;
		case 2:
		  tpol = 7.0*pow(cost,2)-1.0;
		  break;
		}
	      break;
	    case 5:
	      switch (abs_mv)
		{
		case 0:
		  tpol=63.0*pow(cost,5)-70.0*pow(cost,3)+15.0*cost;
		  break;
		case 1:
		  tpol = 21*pow(cost,4)-14.0*pow(cost,2)+1.0;
		  break;
		case 2:
		  tpol= 3.0*pow(cost,3)-cost;
		  break;
		case 3:
		  tpol = 9.0*pow(cost,2)-1.0;
		  break;
		}
	      break;
	    case 6:
	      switch (abs_mv)
		{
		case 0:
		  tpol = 231.0*pow(cost,6)-315.0*pow(cost,4)+105.0*pow(cost,2)-5.0;
		  break;
		case 1:
		  tpol = 33.0*pow(cost,5)-30.0*pow(cost,3)+5.0*cost;
		  break;
		case 2:
		  tpol = 33.0*pow(cost,4)-18.0*pow(cost,2)+1.0;
		  break;
		case 3:
		  tpol = 11.0*pow(cost,3)-3.0*cost;
		  break;
		case 4:
		  tpol = 11.0*pow(cost,2)-1.0;
		  break;
		}
	      break;
	    case 7:
	      switch (abs_mv)
		{
		case 0:
		  tpol=429.0*pow(cost,7)-693.0*pow(cost,5)+315.0*pow(cost,3)-35.0*cost;
		  break;
		case 1:
		  tpol = 429.0*pow(cost,6)-495.0*pow(cost,4)+135.0*pow(cost,2)-5.0;
		  break;
		case 2:
		  tpol = 143.0*pow(cost,5)-110.0*pow(cost,3)+15.0*cost;
		  break;
		case 3:
		  tpol = 143.0*pow(cost,4)-66.0*pow(cost,2)+3.0;
		  break;
		case 4:
		  tpol = 13.0*pow(cost,3)-3.0*cost;
		  break;
		case 5:
		  tpol = 13.0*pow(cost,2)-1.0;
		  break;
		}
	      break;
	    case 8:
	       switch (abs_mv)
		{
		case 0: 
		  tpol=6435.0*pow(cost,8)-12012.0*pow(cost,6)+6930.0*pow(cost,4)-1260.0*pow(cost,2)+35.0;
		  break;
		case 1:
		  tpol = 715.0*pow(cost,7)-1001.0*pow(cost,5)+385.0*pow(cost,3)-35.0*cost;
		  break;
		case 2:
		  tpol = 143.0*pow(cost,6)-143.0*pow(cost,4)+33.0*pow(cost,2)-1.0;
		  break;
		case 3:
		 tpol = 39.0*pow(cost,5)-26.0*pow(cost,3)+3.0*cost;
		 break;
		case 4:
		  tpol = 65.0*pow(cost,4)-26.0*pow(cost,2)+1.0;
		  break;
		case 5:
		  tpol = 5.0*pow(cost,3)-cost;
		  break;
		case 6:
		 tpol = 15.0*pow(cost,2)-1.0;
		 break;
		}
	       break;
	    case 9:
	      switch (abs_mv)
		{
		case 0: 
		  tpol = 12155.0*pow(cost,9)-25740.0*pow(cost,7)+18018.0*pow(cost,5)-4620.0*pow(cost,3)+315.0*cost;
		  break;
		  case 1:
		    tpol = 2431.0*pow(cost,8)-4004.0*pow(cost,6)+2002.0*pow(cost,4)-308.0*pow(cost,2)+7.0;
		    break;
		case 2:
		  tpol = 221.0*pow(cost,7)-273.0*pow(cost,5)+91.0*pow(cost,3)-7.0*cost;
		  break;
		case 3:
		    tpol = 221.0*pow(cost,6)-195.0*pow(cost,4)+39.0*pow(cost,2)-1.0;
		    break;
		case 4:
		  tpol = 17*pow(cost,5)-10.0*pow(cost,3)+cost;
		  break;
		case 5:
		  tpol = 85.0*pow(cost,4)-30.0*pow(cost,2)+1.0;
		  break;
		case 6:
		  tpol = 17.0*pow(cost,3)-3.0*cost;
		  break;
		case 7:
		  tpol = 17.0*pow(cost,2)-1.0;
		  break;
		}
	      break;
	    case 10:
	      switch (abs_mv)
		{
		case 0: 
		  tpol = 46189.0*pow(cost,10)-109395.0*pow(cost,8)+90090.0*pow(cost,6)-30030.0*pow(cost,4)+3465.0*pow(cost,2)-63.0;
		  break;
		case 1:
		  tpol = 4199.0*pow(cost,9)-7956.0*pow(cost,7)+4914.0*pow(cost,5)-1092.0*pow(cost,3)+63.0*cost;
		  break;
		case 2:
		tpol = 4199.0*pow(cost,8)-6188.0*pow(cost,6)+2730.0*pow(cost,4)-364.0*pow(cost,2)+7.0;
		break;
		case 3:
		  tpol = 323.0*pow(cost,7)-357.0*pow(cost,5)+105.0*pow(cost,3)-7.0*cost;
		  break;
		case 4:
		  tpol = 323.0*pow(cost,6)-255.0*pow(cost,4)+45.0*pow(cost,2)-1;
		 break;
		case 5:
		  tpol = 323.0*pow(cost,5)-170.0*pow(cost,3)+15.0*cost;
		  break;
		case 6:
		  tpol = 323.0*pow(cost,4)-102.0*pow(cost,2)+3.0;
		  break;
		case 7:
		  tpol = 19.0*pow(cost,3)-3.0*cost;
		  break;
		case 8:
		  tpol = 19.0*pow(cost,2)-1.0;
		  break;
		}
	      break;
	    default:
	      printf("l=%d m=%d\n", lv,mv); fflush(stdout);
	      error("SPH::GetValue(): Invalid l");
	      break;
	    }
	}

    }
  
  tpol*=thetacoeff[lv][abs_mv];
 

#if 0 
  //now take care of i^l term
  int sign=lv%4; //Phases are 1, i, -1, -i 
  Complex ylm= (!sign) ? Complex(tpol,0.0) : 
    ((1==sign)? Complex(0.0,tpol) : 
     ((2==sign)? Complex(-tpol,0.0) : Complex(0.0,-tpol)));
#endif 

  //now take care of i^-l term
//  int sign=lv%4; 
  //Phases i^l are 1, i, -1, -i; 
  //i^-l are: 1, -i, -1, i  
//  Complex ylm= (!sign) ? Complex(tpol,0.0) : 
//    ((1==sign)? Complex(0.0, -tpol) : 
//     ((2==sign)? Complex(-tpol,0.0) : Complex(0.0,tpol)));

  Complex ylm = tpol;
  if(abs_mv) // sin and e^imf term
    {
// SG: epsint is e^ (m*i*phi) * sin(theta)
      Complex epsint(re,im);
      ylm*= (abs_mv==1) ? epsint : CPower(epsint,abs_mv);
      
      if( mv<0 && mv%2 ) //odd mv && mv <0
	ylm*=-1.0;
    }
  



  return ylm;
}


SPH& theSPH()
{
  static SPH thesph;
  return thesph;
}


