/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
#define TINY 1.e-16;
void ludcmp(const banded_symmetric_dmatrix& _a,const ivector& _indx,
  const MY_DOUBLE_TYPE& _d)
{
  int i=0;
  int imax=0;
  int j=0;
  int k=0;
  int n=0;
  MY_DOUBLE_TYPE& d=(MY_DOUBLE_TYPE&)_d;
  dmatrix& a=(dmatrix&)_a;
  ivector& indx=(ivector&)_indx;

  n=a.colsize();
  int lb=a.colmin();
  int ub=a.colmax();

  MY_DOUBLE_TYPE big,dum,sum,temp;

  dvector vv(lb,ub);


  d=1.0;

  for (i=lb;i<=ub;i++)
  {
    big=0.0;
    for (j=lb;j<=ub;j++)
    {
      temp=fabs(a[i][j]);
      if (temp > big)
      {
        big=temp;
      }
    }
    if (big == 0.0) 
    {
      cerr << "Error in matrix inverse -- matrix singular in inv(dmatrix)\n";
    }
    vv[i]=1.0/big;
  }



  for (j=lb;j<=ub;j++)
  {
    for (i=lb;i<j;i++) 
    {
      sum=a[i][j];
      for (k=lb;k<i;k++)
      {
        sum = sum - a[i][k]*a[k][j];
      }
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<=ub;i++) 
    {
      sum=a[i][j];
      for (k=lb;k<j;k++)
      {
        sum = sum - a[i][k]*a[k][j];
      }
      a[i][j]=sum;
      dum=vv[i]*fabs(sum);
      if ( dum >= big)
      {
        big=dum;
        imax=i;
      }
    }
    if (j != imax)
    {
      for (k=lb;k<=ub;k++)
      {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      d = -d;
      vv[imax]=vv[j];
    }
    indx[j]=imax;

    if (a[j][j] == 0.0)
    {
      a[j][j]=TINY;
    }

    if (j != n)
    {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=ub;i++)
      {
        a[i][j] = a[i][j] * dum;
      }
    }
  }
}
