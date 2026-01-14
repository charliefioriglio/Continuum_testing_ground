#ifndef _klmpoints_h
#define _klmpoints_h

#include "sph.h"
#include "tools.h"
#include "aik_xml_parser.h"

//! Collection of k,l,m points.
class KLMPoints
{
  //Number of k values.
  int nkv;
  //! min and max values for k.
  double kmin;
  double kmax;
  //! grid of wavevector points, k.
  double *kv;
  //! Number of (l,m) values.
  int nlmv;
  //! max ang mom l considered.
  int lmax;
  //! Arrays of l,m values as funtion of the absolute address
  int *lv;
  int *mv;

  bool IfAlloc() const { return nkv*nlmv; }
  void Alloc(int nk, int nlm);
  void Free();
  void Set();
  //void AllocAndCopyData(int nk,double kmn,double kmx,int nlm,int lmx,double *k,int *l,int *m);
  
  //!Generate k grid: kv[nkv].
  void CalcKGrid();
  //!From l and m matrices. lv matrix does not have unique values:
  //! same l for each m value,e.g. lv[0 1 1 1], mv[0 -1 0 1]
  void CalcLMPoints();

  KLMPoints(const KLMPoints& other);
  KLMPoints& operator=(const KLMPoints& other);

 public:
  KLMPoints();
  //  KLMPoints(int nkv,double kmn,double kmx,int nlm, int lmx);
  KLMPoints(const char* xmlFileName);
  ~KLMPoints();


  //! Print k,l,m grid info.
  void PrintGridInfo(FILE *outf=stdout) const;

  //! Get k and l value for the n-th spherical wave, by absolute address.
  void GetKLPoint(int absaddr, double &kpt, double &lpt);

  inline int NKV() const { return nkv; } 
  inline int NLMV() const { return nlmv; }
  inline int LMax() const { return lmax; }
  inline double GetKV(int i) const { return kv[i]; }
  //! Returns l correspomdig to point i
  inline int GetLV(int i) const { return lv[i]; }
  //! Returns m correspomdig to point i
  inline int GetMV(int i) const { return mv[i]; }
  int GetKLPoint(int i) const;
};


#endif
