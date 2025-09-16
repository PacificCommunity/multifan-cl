/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

dvector operator * (const dvector& _v,const banded_lower_triangular_dmatrix& M)
{
  ADUNCONST(dvector,v)
  int vmin=v.indexmin();
  int vmax=v.indexmax();
  int mmin=M.indexmin();
  int mmax=M.indexmax();
  int bw=M.bandwidth();
  dvector tmp(vmin,vmax);
  if (vmin != mmin || vmax != mmax)
  {
    cerr << "size mismatch in" << endl <<
      "dvector operator * (const dvector& v,"
      "const banded_lower_triangular_dmatrix& M)" << endl;
  }
  for (int i=mmin;i<=mmax;i++)
  {
    MY_DOUBLE_TYPE ssum=0.0;
    for (int j=i;j<=min(i+bw-1,mmax);j++)
    {
      ssum+=v(j)*M(j,i);
    }
    tmp(i)=ssum;
  }
  return tmp;
}

/*
dmatrix make_dmatrix(const banded_lower_triangular_dmatrix& M)
{
  int mmin=M.indexmin();
  int mmax=M.indexmax();
  int bw=M.bandwidth();
  dmatrix tmp(mmin,mmax,mmin,mmax);
  tmp.initialize();
  for (int i=mmin;i<=mmax;i++)
  {
    MY_DOUBLE_TYPE ssum=0.0;
    for (int j=i;j<=min(i+bw-1,mmax);j++)
    {
      tmp(j,i)=M(j,i);
    }
  }
  return tmp;
}
*/
