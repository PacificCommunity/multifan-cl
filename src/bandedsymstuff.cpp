/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include <admodel.h>

dvector solve(const banded_symmetric_dmatrix& _S,const dmatrix& u,
  const dmatrix& v,const dvector& x)
{
  ADUNCONST(banded_symmetric_dmatrix,S)
  int mmin=S.indexmin();
  int mmax=S.indexmax();
  if (mmin != u(u.indexmin()).indexmin() ||
    mmax != u(u.indexmin()).indexmax() ||
    mmin != v(u.indexmin()).indexmin() ||
    mmax != v(u.indexmin()).indexmax() ||
    x.indexmin() != mmin ||
    x.indexmax() != mmax )
  {
    cerr << "shape error in "
    "solve(const symmetric_tridiagonal_dmatrix& S,const dmatrix& u," << endl;
    ad_exit(1);
  }
  int imin=u.indexmin();
  int imax=u.indexmax();
  dmatrix z(imin,imax,mmin,mmax);
  dmatrix w(imin,imax,mmin,mmax);
  dvector lambda(imin,imax);
  for (int i=imin;i<=imax;i++)
  {
    z(i) = solve(S,u(i));
    w(i) = solve(S,v(i));
    for (int j=imin;j<i;j++)
    {
      z(i) -= (w(j)*u(i))/(1+lambda(j))*z(j);
      w(i) -= (z(j)*v(i))/(1+lambda(j))*w(j);
    }
    lambda(i)=v(i)*z(i);
  }

  dvector y = solve(S,x);
  for (int i=imin;i<=imax;i++)
  {
    y-=(w(i)*x)/(1+lambda(i))*z(i);
  }
  return y;
}
