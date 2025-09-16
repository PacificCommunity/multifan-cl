/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include <admodel.h>
#include "phi_stuff.h"
#include "tridiagonal_dmatrix.h"
  dmatrix make_dmatrix(const symmetric_tridiagonal_dmatrix& _S)
  {
    ADUNCONST(symmetric_tridiagonal_dmatrix,S)
    int mmin=S.indexmin();
    int mmax=S.indexmax();
    dmatrix tmp(mmin,mmax,mmin,mmax);
    tmp.initialize();
    tmp(mmin,mmin)=S(mmin,mmin);
    tmp(mmin,mmin+1)=S(mmin,mmin+1);
    tmp(mmin+1,mmin)=S(mmin+1,mmin);
    for (int i=mmin+1;i<mmax;i++)
    {
      tmp(i+1,i)=S(i+1,i);
      tmp(i,i)=S(i,i);
      tmp(i,i+1)=S(i,i+1);
    }
    tmp(mmax,mmax-1)=S(mmax,mmax-1);
    tmp(mmax,mmax)=S(mmax,mmax);
    return tmp;
  }
    

dvector operator * (const symmetric_tridiagonal_dmatrix & _S,
  const dvector & _v)
{
  ADUNCONST(symmetric_tridiagonal_dmatrix,S)
  ADUNCONST(dvector,v)
   if (v.indexmin() != S.indexmin() ||
       v.indexmax() != S.indexmax() )
   {
     cerr << "incompatible shape in " 
       "dvector operator * (const symmetric_tridiagonal_dmatrix & S,"
       "const symmetric_tridiagonal_dmatrix & Sinv);" << endl;
     ad_exit(1);
   }
   int mmin=v.indexmin();
   int mmax=v.indexmax();
   dvector b=S.get_diag();
   dvector a=S.get_subdiag();
   dvector tmp(mmin,mmax);
   tmp(mmin)=b(1)*v(1)+a(2)*v(2);
   for (int i=mmin+1;i<mmax;i++)
   {
     tmp(i)=a(i)*v(i-1)+b(i)*v(i)+a(i+1)*v(i+1);
   }
   tmp(mmax)=a(mmax)*v(mmax-1)+b(mmax)*v(mmax);
   return tmp;
}

dvector solve(const symmetric_tridiagonal_dmatrix& _S,const dmatrix& u,
  const dmatrix& v,const dvector& x)
{
  ADUNCONST(symmetric_tridiagonal_dmatrix,S)
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
MY_DOUBLE_TYPE ln_det(const symmetric_tridiagonal_dmatrix& _S,const dmatrix& u,
  const dmatrix& v,const int& _ierr)
{
  ADUNCONST(symmetric_tridiagonal_dmatrix,S)
  ADUNCONST(int,ierr)
  int mmin=S.indexmin();
  int mmax=S.indexmax();
  if (mmin != u(u.indexmin()).indexmin() ||
    mmax != u(u.indexmin()).indexmax() ||
    mmin != v(u.indexmin()).indexmin() ||
    mmax != v(u.indexmin()).indexmax() )
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

  MY_DOUBLE_TYPE ld=2.0*ln_det_choleski(S,ierr);
  for (int i=imin;i<=imax;i++)
  {
    MY_DOUBLE_TYPE tmp=1.0+v(i)*z(i);
    if (tmp<=0.0)
    {
      cerr << "error in "
     "solve(const symmetric_tridiagonal_dmatrix& S,const dmatrix& u," << endl;
      ierr=1;
      ad_exit(1);
    }
    ld+=log(tmp);
  }
  return ld;
}

