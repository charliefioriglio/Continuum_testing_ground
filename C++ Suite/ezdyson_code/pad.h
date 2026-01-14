//! Class for PAD 
#include "cklm.h"
#include "readwrite.h"

class PAD
{
  int ntheta;
  int nkv;
  double *pad;
  double *crosspad;
  double *totalpad;
  
  double DeltaTheta() const { return M_PI/(ntheta-1); }

 public:
  PAD(int nth, int nkv);
  ~PAD();
  void CalcPad(CklmCoeff& allCklm, double scale_coeff=1.0) const;
  void CalcMCPadWithCoherences(KLMPoints& klmgrid, CklmCoeff& allCklm1, CklmCoeff& allCklm2,
			       double *center1, double *center2,
			       double scale_coeff1,double scale_coeff2, double *xsec) const;
  void Print(const KLMPoints& klmgrid) const;
  double TotalXSec(int kv) const;
  //! returns parallel x-sec, i.e. theta=0
  double XSec_par(int kv) const; 
  //! returns perpendicular x-sec, i.e. theta=pi/2
  double XSec_perp(int kv) const; 
};
