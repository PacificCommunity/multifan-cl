/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include <admodel.h>
#include "tridiagonal_dmatrix.h"
dvector operator *(const dvector& r,const symmetric_tridiagonal_dmatrix& S)
{
  return S*r;
}

dvector solve(symmetric_tridiagonal_dmatrix& S,const dvector& r)
{
  int mmin=S.indexmin();
  int mmax=S.indexmax();
  dvector b=S.get_diag();
  dvector a=S.get_subdiag();

  dvector  bet(mmin,mmax);
  dvector  gam(mmin,mmax);
  dvector  u(mmin,mmax);
  dvector  v(mmin,mmax);

  if (S.get_diag(mmin) == 0.0) cerr <<"Error 2 in solve"<< endl;
  bet(1)=b(mmin);
  u(mmin)=r(mmin)/(b(mmin));
  for (int j=mmin+1;j<=mmax;j++) 
  {
    gam[j]=a(j)/bet(j-1);
    bet(j)=b(j)-a(j)*gam(j);
    if (bet(j) == 0.0) cerr <<"Error 2 in solve"<< endl;
    u[j]=(r[j]-a(j)*u(j-1))/bet(j);
  }
  v=u;
  for (int j=(mmax-1);j>=mmin;j--)
    v[j] -= gam[j+1]*v[j+1];
  return v;
}
MY_DOUBLE_TYPE  ln_det_choleski(const symmetric_tridiagonal_dmatrix& _M,const int& _ierr)
{
  ADUNCONST(int,ierr)
  ADUNCONST(symmetric_tridiagonal_dmatrix,M)
  int mmin=M.indexmin();
  int mmax=M.indexmax();
  
  symmetric_tridiagonal_dmatrix L(mmin,mmax);
#ifndef SAFE_INITIALIZE
    L.initialize();
#endif

  int i,j,k;
  MY_DOUBLE_TYPE tmp;
  if (M(mmin,mmin)<=0)
  {
    if (ierr==0)
      cerr << "Error matrix not positive definite in choleski_decomp"
        <<endl;
    ierr=1;
    return 0;
  }
  L(mmin,mmin)=sqrt(M(mmin,mmin));
  for (i=mmin;i<=mmin+1;i++)
  {
    L(i,mmin)=M(i,mmin)/L(mmin,mmin);
  }

  for (i=mmin+1;i<=mmax;i++)
  {
    if (i>2)
    {	
      tmp=M(i,i-1);
      L(i,i-1)=tmp/L(i-1,i-1);
    }
    tmp=M(i,i);
    if (i>1)	
      tmp-=L(i,i-1)*L(i,i-1);
    if (tmp<=0)
    {
      if (ierr==0)
        cerr << "Error matrix not positive definite in choleski_decomp"
          <<endl;
      ierr=1;
      return 0;
    }
    L(i,i)=sqrt(tmp);
  }
  MY_DOUBLE_TYPE ln_det=0.0;
  for (int i=mmin;i<=mmax;i++)
  {
    ln_det+=log(L(i,i));
  }
  return ln_det;
}
