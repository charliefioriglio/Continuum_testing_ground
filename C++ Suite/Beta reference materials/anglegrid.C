#include "anglegrid.h"

void AngleGrid::Alloc()
{
  Free();
  nabc = 150;
  for (int i=0; i<ABG ; i++)
   {
     thegrid[i] = new double [nabc];
     memset(thegrid[i],0,sizeof(double)*nabc);
   }
}

void AngleGrid::Free()
{
  if (IfAlloc())
   {
     for (int i=0; i<ABG; i++)
       if (NULL!=thegrid[i])
	 {
	   delete [] thegrid[i]; 
	   thegrid[i]=NULL;
	 }
     nabc=0;
   }
}

void AngleGrid::Set()
{
  nabc=0;
}

//AngleGrid::AngleGrid(const AngleGrid& other)
//{
//  Set();
//}

AngleGrid::AngleGrid()
{
  Set();
  Alloc();
  CalcAngleGrid();
}

AngleGrid::~AngleGrid() 
{ 
  Free();
}

//!Calc da,db,dc for the alpha, beta, gamma grid.
/*void AngleGrid::CalcDaDbDc()
{
  double range;
  for (int i=0; i<ABG; i++)
    {
      if (npoints[i]>0)
	{
	  range = abcmax[i]-abcmin[i];
	  if (0==i || 2==i)
            dabc[i] = M_PI*range/npoints[i];
	  if (1==i) //AIK: i do not get this, but it gives correct iso average
            dabc[i] = M_PI*range/(npoints[i]-1);
	}
      else
	dabc[i] = 0.0;
    }
}
*/
//!Generate x,y and z values for *thegrid[ABG].
void AngleGrid::CalcAngleGrid()
{
thegrid[A][0] = 1.07393;
thegrid[B][0] = 2.37487;
thegrid[G][0] = 0.00649606;
thegrid[A][1] = 3.31052;
thegrid[B][1] = 2.01802;
thegrid[G][1] = 0.00677567;
thegrid[A][2] = 1.74962;
thegrid[B][2] = 0.19298;
thegrid[G][2] = 0.00640311; 
thegrid[A][3] = 1.64604;
thegrid[B][3] = 1.42732;
thegrid[G][3] = 0.00666473; 
thegrid[A][4] = 6.10750;
thegrid[B][4] = 0.94148;
thegrid[G][4] = 0.00669922; 
thegrid[A][5] = 1.75961;
thegrid[B][5] = 2.87768;
thegrid[G][5] = 0.00663394; 
thegrid[A][6] = 5.53027;
thegrid[B][6] = 0.43346;
thegrid[G][6] = 0.00674343;
thegrid[A][7] = 5.52842;
thegrid[B][7] = 2.23837;
thegrid[G][7] = 0.00667127; 
thegrid[A][8] = 2.28479;
thegrid[B][8] = 2.38175;
thegrid[G][8] = 0.00668776;
thegrid[A][9] = 3.58752;
thegrid[B][9] = 1.86742;
thegrid[G][9] = 0.00676721; 
thegrid[A][10] = 1.66523;
thegrid[B][10] = 0.48430;
thegrid[G][10] = 0.00665191; 
thegrid[A][11] = 0.84226;
thegrid[B][11] = 1.47576;
thegrid[G][11] = 0.00663109; 
thegrid[A][12] = 3.29536;
thegrid[B][12] = 2.33447;
thegrid[G][12] = 0.00663458; 
thegrid[A][13] = 5.57673;
thegrid[B][13] = 1.61340;
thegrid[G][13] = 0.00674962; 
thegrid[A][14] = 2.47174;
thegrid[B][14] = 1.85509;
thegrid[G][14] = 0.00667875; 
thegrid[A][15] = 0.92423;
thegrid[B][15] = 2.01268;
thegrid[G][15] = 0.00676506; 
thegrid[A][16] = 3.93291;
thegrid[B][16] = 0.75689;
thegrid[G][16] = 0.00663678; 
thegrid[A][17] = 3.31288;
thegrid[B][17] = 1.06096;
thegrid[G][17] = 0.00664742; 
thegrid[A][18] = 1.90607;
thegrid[B][18] = 1.28682;
thegrid[G][18] = 0.00663493; 
thegrid[A][19] = 5.84085;
thegrid[B][19] = 2.45608;
thegrid[G][19] = 0.00677988; 
thegrid[A][20] = 0.38988;
thegrid[B][20] = 1.22613;
thegrid[G][20] = 0.00666949; 
thegrid[A][21] = 3.04104;
thegrid[B][21] = 0.85503;
thegrid[G][21] = 0.00678447; 
thegrid[A][22] = 3.73432;
thegrid[B][22] = 2.48750;
thegrid[G][22] = 0.00678442; 
thegrid[A][23] = 0.06801;
thegrid[B][23] = 2.59037;
thegrid[G][23] = 0.00674655; 
thegrid[A][24] = 4.53341;
thegrid[B][24] = 1.12758;
thegrid[G][24] = 0.00660888; 
thegrid[A][25] = 3.95200;
thegrid[B][25] = 1.08387;
thegrid[G][25] = 0.00662420; 
thegrid[A][26] = 4.40989;
thegrid[B][26] = 1.71823;
thegrid[G][26] = 0.00663710; 
thegrid[A][27] = 3.00258;
thegrid[B][27] = 2.89856;
thegrid[G][27] = 0.00638865; 
thegrid[A][28] = 0.72053;
thegrid[B][28] = 1.19880;
thegrid[G][28] = 0.00674299; 
thegrid[A][29] = 3.05150;
thegrid[B][29] = 1.82581;
thegrid[G][29] = 0.00667288; 
thegrid[A][30] = 3.43130;
thegrid[B][30] = 0.72273;
thegrid[G][30] = 0.00663664; 
thegrid[A][31] = 4.27813;
thegrid[B][31] = 2.03217;
thegrid[G][31] = 0.00663669; 
thegrid[A][32] = 0.33053;
thegrid[B][32] = 2.31736;
thegrid[G][32] = 0.00669460; 
thegrid[A][33] = 0.21824;
thegrid[B][33] = 1.47451;
thegrid[G][33] = 0.00674184; 
thegrid[A][34] = 1.26162;
thegrid[B][34] = 2.64602;
thegrid[G][34] = 0.00666331; 
thegrid[A][35] = 2.18568;
thegrid[B][35] = 1.75559;
thegrid[G][35] = 0.00668894; 
thegrid[A][36] = 4.36603;
thegrid[B][36] = 0.58855;
thegrid[G][36] = 0.00668811; 
thegrid[A][37] = 0.72968;
thegrid[B][37] = 1.75632;
thegrid[G][37] = 0.00678161; 
thegrid[A][38] = 0.09293;
thegrid[B][38] = 1.74660;
thegrid[G][38] = 0.00666844; 
thegrid[A][39] = 1.39264;
thegrid[B][39] = 1.57420;
thegrid[G][39] = 0.00649848; 
thegrid[A][40] = 0.61877;
thegrid[B][40] = 0.64023;
thegrid[G][40] = 0.00672334; 
thegrid[A][41] = 3.83971;
thegrid[B][41] = 1.67455;
thegrid[G][41] = 0.00663329; 
thegrid[A][42] = 1.33503;
thegrid[B][42] = 1.85446;
thegrid[G][42] = 0.00657182; 
thegrid[A][43] = 2.72987;
thegrid[B][43] = 1.02968;
thegrid[G][43] = 0.00678480; 
thegrid[A][44] = 5.64007;
thegrid[B][44] = 1.28404;
thegrid[G][44] = 0.00678987; 
thegrid[A][45] = 0.64928;
thegrid[B][45] = 2.56314;
thegrid[G][45] = 0.00668689; 
thegrid[A][46] = 4.61242;
thegrid[B][46] = 1.94129;
thegrid[G][46] = 0.00668108; 
thegrid[A][47] = 3.37637;
thegrid[B][47] = 0.21312;
thegrid[G][47] = 0.00665181; 
thegrid[A][48] = 5.84322;
thegrid[B][48] = 0.13739;
thegrid[G][48] = 0.00659894; 
thegrid[A][49] = 3.88503;
thegrid[B][49] = 1.37820;
thegrid[G][49] = 0.00675097; 
thegrid[A][50] = 3.30267;
thegrid[B][50] = 1.37309;
thegrid[G][50] = 0.00676891; 
thegrid[A][51] = 2.80485;
thegrid[B][51] = 2.42800;
thegrid[G][51] = 0.00677254; 
thegrid[A][52] = 5.96564;
thegrid[B][52] = 1.21872;
thegrid[G][52] = 0.00678260; 
thegrid[A][53] = 3.63200;
thegrid[B][53] = 2.18767;
thegrid[G][53] = 0.00672571; 
thegrid[A][54] = 3.21888;
thegrid[B][54] = 2.62084;
thegrid[G][54] = 0.00664667; 
thegrid[A][55] = 5.15715;
thegrid[B][55] = 2.30207;
thegrid[G][55] = 0.00640289; 
thegrid[A][56] = 6.16719;
thegrid[B][56] = 2.00292;
thegrid[G][56] = 0.00656302; 
thegrid[A][57] = 5.35457;
thegrid[B][57] = 1.40952;
thegrid[G][57] = 0.00671165; 
thegrid[A][58] = 0.54376;
thegrid[B][58] = 1.48686;
thegrid[G][58] = 0.00670240; 
thegrid[A][59] = 4.22616;
thegrid[B][59] = 1.22480;
thegrid[G][59] = 0.00655946; 
thegrid[A][60] = 4.85400;
thegrid[B][60] = 2.13311;
thegrid[G][60] = 0.00665147; 
thegrid[A][61] = 3.31674;
thegrid[B][61] = 1.69066;
thegrid[G][61] = 0.00671663; 
thegrid[A][62] = 1.37092;
thegrid[B][62] = 1.25267;
thegrid[G][62] = 0.00668897; 
thegrid[A][63] = 4.28857;
thegrid[B][63] = 2.52738;
thegrid[G][63] = 0.00670018; 
thegrid[A][64] = 1.13002;
thegrid[B][64] = 1.45655;
thegrid[G][64] = 0.00646668; 
thegrid[A][65] = 0.71738;
thegrid[B][65] = 2.25519;
thegrid[G][65] = 0.00646788; 
thegrid[A][66] = 0.54974;
thegrid[B][66] = 2.00568;
thegrid[G][66] = 0.00663172; 
thegrid[A][67] = 0.55092;
thegrid[B][67] = 2.85707;
thegrid[G][67] = 0.00672750; 
thegrid[A][68] = 1.24132;
thegrid[B][68] = 2.12000;
thegrid[G][68] = 0.00657072; 
thegrid[A][69] = 3.64290;
thegrid[B][69] = 0.96483;
thegrid[G][69] = 0.00637692; 
thegrid[A][70] = 1.64918;
thegrid[B][70] = 1.75524;
thegrid[G][70] = 0.00672365; 
thegrid[A][71] = 0.21685;
thegrid[B][71] = 2.03425;
thegrid[G][71] = 0.00674388; 
thegrid[A][72] = 3.08618;
thegrid[B][72] = 0.51224;
thegrid[G][72] = 0.00667000; 
thegrid[A][73] = 0.92048;
thegrid[B][73] = 0.89061;
thegrid[G][73] = 0.00666926; 
thegrid[A][74] = 2.12227;
thegrid[B][74] = 0.67305;
thegrid[G][74] = 0.00675168; 
thegrid[A][75] = 2.48303;
thegrid[B][75] = 1.24734;
thegrid[G][75] = 0.00664912; 
thegrid[A][76] = 5.11643;
thegrid[B][76] = 0.87917;
thegrid[G][76] = 0.00663779; 
thegrid[A][77] = 4.98922;
thegrid[B][77] = 0.58117;
thegrid[G][77] = 0.00677741; 
thegrid[A][78] = 5.85691;
thegrid[B][78] = 1.50352;
thegrid[G][78] = 0.00670043; 
thegrid[A][79] = 5.20647;
thegrid[B][79] = 2.01029;
thegrid[G][79] = 0.00660064; 
thegrid[A][80] = 0.55287;
thegrid[B][80] = 0.93519;
thegrid[G][80] = 0.00656175; 
thegrid[A][81] = 4.78971;
thegrid[B][81] = 1.36289;
thegrid[G][81] = 0.00666901; 
thegrid[A][82] = 0.20077;
thegrid[B][82] = 0.97133;
thegrid[G][82] = 0.00638311; 
thegrid[A][83] = 1.62743;
thegrid[B][83] = 0.77734;
thegrid[G][83] = 0.00668824; 
thegrid[A][84] = 1.16900;
thegrid[B][84] = 0.64195;
thegrid[G][84] = 0.00678113; 
thegrid[A][85] = 5.35471;
thegrid[B][85] = 2.56591;
thegrid[G][85] = 0.00665192; 
thegrid[A][86] = 3.02272;
thegrid[B][86] = 1.19311;
thegrid[G][86] = 0.00672800; 
thegrid[A][87] = 6.06776;
thegrid[B][87] = 1.72994;
thegrid[G][87] = 0.00638353; 
thegrid[A][88] = 2.18861;
thegrid[B][88] = 1.16124;
thegrid[G][88] = 0.00639049; 
thegrid[A][89] = 4.82467;
thegrid[B][89] = 2.76802;
thegrid[G][89] = 0.00675165; 
thegrid[A][90] = 5.11437;
thegrid[B][90] = 1.23525;
thegrid[G][90] = 0.00663724; 
thegrid[A][91] = 2.76769;
thegrid[B][91] = 1.94278;
thegrid[G][91] = 0.00639167; 
thegrid[A][92] = 1.87105;
thegrid[B][92] = 2.56288;
thegrid[G][92] = 0.00676403; 
thegrid[A][93] = 4.61828;
thegrid[B][93] = 0.31374;
thegrid[G][93] = 0.00675114; 
thegrid[A][94] = 4.49672;
thegrid[B][94] = 2.25460;
thegrid[G][94] = 0.00666951; 
thegrid[A][95] = 3.57552;
thegrid[B][95] = 1.54648;
thegrid[G][95] = 0.00662358; 
thegrid[A][96] = 2.38265;
thegrid[B][96] = 0.91810;
thegrid[G][96] = 0.00666779; 
thegrid[A][97] = 2.76245;
thegrid[B][97] = 1.35613;
thegrid[G][97] = 0.00663807; 
thegrid[A][98] = 4.72005;
thegrid[B][98] = 1.64956;
thegrid[G][98] = 0.00668906; 
thegrid[A][99] = 4.04847;
thegrid[B][99] = 2.27038;
thegrid[G][99] = 0.00678518; 
thegrid[A][100] = 4.67308;
thegrid[B][100] = 0.80996;
thegrid[G][100] = 0.00666891; 
thegrid[A][101] = 0.03871;
thegrid[B][101] = 0.41045;
thegrid[G][101] = 0.00666294; 
thegrid[A][102] = 5.47292;
thegrid[B][102] = 0.74431;
thegrid[G][102] = 0.00671233; 
thegrid[A][103] = 0.39803;
thegrid[B][103] = 1.74535;
thegrid[G][103] = 0.00670227; 
thegrid[A][104] = 5.04133;
thegrid[B][104] = 1.53611;
thegrid[G][104] = 0.00677799; 
thegrid[A][105] = 1.91134;
thegrid[B][105] = 1.60435;
thegrid[G][105] = 0.00676458; 
thegrid[A][106] = 2.47181;
thegrid[B][106] = 1.55181;
thegrid[G][106] = 0.00677420; 
thegrid[A][107] = 2.18850;
thegrid[B][107] = 1.45794;
thegrid[G][107] = 0.00666742; 
thegrid[A][108] = 0.02788;
thegrid[B][108] = 1.21354;
thegrid[G][108] = 0.00659765; 
thegrid[A][109] = 1.64398;
thegrid[B][109] = 1.09534;
thegrid[G][109] = 0.00672893; 
thegrid[A][110] = 2.99010;
thegrid[B][110] = 2.15902;
thegrid[G][110] = 0.00665960; 
thegrid[A][111] = 5.39785;
thegrid[B][111] = 1.07567;
thegrid[G][111] = 0.00677759; 
thegrid[A][112] = 2.22545;
thegrid[B][112] = 2.06928;
thegrid[G][112] = 0.00676543; 
thegrid[A][113] = 0.13728;
thegrid[B][113] = 0.69204;
thegrid[G][113] = 0.00657994; 
thegrid[A][114] = 4.95468;
thegrid[B][114] = 1.83188;
thegrid[G][114] = 0.00675244; 
thegrid[A][115] = 3.95088;
thegrid[B][115] = 2.77746;
thegrid[G][115] = 0.00666695; 
thegrid[A][116] = 4.48075;
thegrid[B][116] = 1.42929;
thegrid[G][116] = 0.00672720; 
thegrid[A][117] = 3.89107;
thegrid[B][117] = 1.99852;
thegrid[G][117] = 0.00664683; 
thegrid[A][118] = 5.28079;
thegrid[B][118] = 1.71301;
thegrid[G][118] = 0.00674483; 
thegrid[A][119] = 1.04117;
thegrid[B][119] = 1.17178;
thegrid[G][119] = 0.00669447; 
thegrid[A][120] = 1.52149;
thegrid[B][120] = 2.35523;
thegrid[G][120] = 0.00672329; 
thegrid[A][121] = 6.18079;
thegrid[B][121] = 1.46674;
thegrid[G][121] = 0.00659747; 
thegrid[A][122] = 5.73604;
thegrid[B][122] = 0.97673;
thegrid[G][122] = 0.00678974; 
thegrid[A][123] = 1.91169;
thegrid[B][123] = 2.24095;
thegrid[G][123] = 0.00678624; 
thegrid[A][124] = 5.92459;
thegrid[B][124] = 0.67027;
thegrid[G][124] = 0.00674990; 
thegrid[A][125] = 5.51662;
thegrid[B][125] = 1.93475;
thegrid[G][125] = 0.00666400; 
thegrid[A][126] = 3.78850;
thegrid[B][126] = 0.47169;
thegrid[G][126] = 0.00668276; 
thegrid[A][127] = 2.46795;
thegrid[B][127] = 2.66075;
thegrid[G][127] = 0.00666590; 
thegrid[A][128] = 4.85331;
thegrid[B][128] = 1.07081;
thegrid[G][128] = 0.00637544; 
thegrid[A][129] = 2.43693;
thegrid[B][129] = 0.42167;
thegrid[G][129] = 0.00664386; 
thegrid[A][130] = 0.84941;
thegrid[B][130] = 0.36992;
thegrid[G][130] = 0.00667255; 
thegrid[A][131] = 4.11191;
thegrid[B][131] = 1.78975;
thegrid[G][131] = 0.00637753; 
thegrid[A][132] = 1.05014;
thegrid[B][132] = 1.74413;
thegrid[G][132] = 0.00676400; 
thegrid[A][133] = 3.60250;
thegrid[B][133] = 1.25456;
thegrid[G][133] = 0.00663239; 
thegrid[A][134] = 2.66381;
thegrid[B][134] = 0.70579;
thegrid[G][134] = 0.00669928; 
thegrid[A][135] = 5.14617;
thegrid[B][135] = 3.06199;
thegrid[G][135] = 0.00664242; 
thegrid[A][136] = 1.60611;
thegrid[B][136] = 2.05766;
thegrid[G][136] = 0.00676106; 
thegrid[A][137] = 5.85919;
thegrid[B][137] = 2.12962;
thegrid[G][137] = 0.00672395; 
thegrid[A][138] = 6.20983;
thegrid[B][138] = 2.29422;
thegrid[G][138] = 0.00666835; 
thegrid[A][139] = 3.03865;
thegrid[B][139] = 1.51664;
thegrid[G][139] = 0.00677843; 
thegrid[A][140] = 5.79858;
thegrid[B][140] = 1.83453;
thegrid[G][140] = 0.00658002; 
thegrid[A][141] = 4.28145;
thegrid[B][141] = 0.90821;
thegrid[G][141] = 0.00672617; 
thegrid[A][142] = 2.75753;
thegrid[B][142] = 1.65087;
thegrid[G][142] = 0.00666076; 
thegrid[A][143] = 4.16086;
thegrid[B][143] = 1.50476;
thegrid[G][143] = 0.00662597; 
thegrid[A][144] = 5.81163;
thegrid[B][144] = 2.78176;
thegrid[G][144] = 0.00668821; 
thegrid[A][145] = 4.79705;
thegrid[B][145] = 2.47069;
thegrid[G][145] = 0.00664327; 
thegrid[A][146] = 1.95620;
thegrid[B][146] = 0.94743;
thegrid[G][146] = 0.00664268; 
thegrid[A][147] = 2.56970;
thegrid[B][147] = 2.17897;
thegrid[G][147] = 0.00667836; 
thegrid[A][148] = 1.91953;
thegrid[B][148] = 1.92705;
thegrid[G][148] = 0.00678659; 
thegrid[A][149] = 1.30330;
thegrid[B][149] = 0.94949;
thegrid[G][149] = 0.00674897; 
}

/*
ISOGrid::ISOGrid(const AngleGrid& anggrid)
{
  nabc = 150;

  for (int i=0; i<nabc; i++)
    {
      isogrid[A][n] = anggrid.GetPoint(A,i);
      isogrid[B][n] = anggrid.GetPoint(B,0);
      isogrid[G][n] = 0.0;
      dbeta[n] = 0.5*fabs(cos(anggrid.GetPoint(B,1)) - cos(anggrid.GetPoint(B,0)));
      n++;
    }   
 
  for (int i=0; i<na; i++)
    for (int j=1; j<nb-1; j++)
      for (int k=0; k<nc; k++)
	{
	  isogrid[A][n] = anggrid.GetPoint(A,i);
	  isogrid[B][n] = anggrid.GetPoint(B,j);
	  isogrid[G][n] = anggrid.GetPoint(C,k);
	  dbeta[n] = anggrid.DABC(C)/(2*M_PI)*0.5*fabs(cos(anggrid.GetPoint(B,j+1)) - cos(anggrid.GetPoint(B,j-1)));
	  n++;
	}
  
  for (int i=0; i<na; i++)
    {
      isogrid[A][n] = anggrid.GetPoint(A,i);
      isogrid[B][n] = anggrid.GetPoint(B,nb-1);
      isogrid[G][n] = 0.0;
      dbeta[n] = 0.5*fabs(cos(anggrid.GetPoint(B,nb-1)) - cos(anggrid.GetPoint(B,nb-2)));
      n++;
    }
}

ISOGrid::~ISOGrid()
{
  delete [] dbeta;
  
  for (int v=0; v<ABG; v++)
    delete [] isogrid[v];
}

//! Print the grid info to std output.
void AngleGrid::PrintGridInfo(FILE *outfile) const
{
  fprintf(outfile,"\nAngleGrid\n");
  fprintf(outfile,"alpha %5d pts   %9.7lf to %9.7lf PI\n",npoints[A],abcmin[A],abcmax[A]);
  fprintf(outfile,"beta  %5d pts   %9.7lf to %9.7lf PI\n",npoints[B],abcmin[B],abcmax[B]);
  fprintf(outfile,"gamma %5d pts   %9.7lf to %9.7lf PI\n",npoints[G],abcmin[G],abcmax[G]);
  fprintf(outfile,"dAlpha dBeta dGamma:   %9.7lf %9.7lf %9.7lf PI\n",dabc[A],dabc[B],dabc[G]);
  fflush(outfile);
}

//!Print to file the grid as thegrid[ABG][nabc].
void AngleGrid::Print(const char *filename) const
{
  FILE *outfile = fopen(filename,"w");
  if (npoints[B]*npoints[G]>0)
   {
     fprintf(outfile,"\n\n         i   Alpha          Beta           Gamma\n");
     for (int n=0,i=0; i<npoints[A]; i++)
        for (int j=0; j<npoints[B]; j++)
          for (int k=0; k<npoints[G]; k++)
            {
              fprintf(outfile,"%10d  %9.7lf  %9.7lf  %9.7lf\n",n,thegrid[A][i]/M_PI,thegrid[B][j]/M_PI,thegrid[G][k]/M_PI);
              n++;
            }
   }
  else
   {
     fprintf(outfile,"\n\n         i   Alpha\n");
     for (int i=0; i<npoints[A]; i++)
         fprintf(outfile,"%10d  %10.8lf\n",i,thegrid[A][i]);
   }        
  fclose(outfile);
}
*/

/*
 void AngleGrid::GetPoint(int absaddr, double& a, double &b, double &c) const
{
  int ia, ka, ja;
  //!Calc the address of alpha, beta, gamma
  int nbc = npoints[B]*npoints[G];
  ia = absaddr/nbc;
  ja = (absaddr-ia*nbc)/npoints[G];
  ka = absaddr-ia*nbc-ja*npoints[G];
  a = GetPoint(A,ia);
  b = GetPoint(B,ja);
  c = GetPoint(C,ka);
}
*/

//!Pointer to Euler angles grid (molec orientations in lab frame) in input.
/*
long PtrToAnglesGrid(FILE *input)
{
  char str6[] = "$AnglesGrid";
  char *ptr = NULL, string[256];

  rewind(input);
  while (ptr == NULL && !feof(input))
      {
        fgets(string,255,input);
        ptr = strstr(string,str6);
      }
  long anglegridptr=ftell(input);

  return anglegridptr;
}
*/
/*
 void ReadAllAngleGridInfo(int* npts, double* abcmin, double* abcmax, const char* xmlFileName)
{
  simpleXMLparser xmlF;
  xmlF.assignFile(xmlFileName); 

  xmlF.reset().node("root").node("averaging").node("angles_grid");

  for (int i=0; i<ABG; i++)
    {
      xmlF.node("angle", i+1);
      npts[i]=xmlF.getIntValue("n_points");
      abcmin[i]=xmlF.getDoubleValue("min");
      abcmax[i]=xmlF.getDoubleValue("max");
      xmlF.stepBack();
    }
}
*/

