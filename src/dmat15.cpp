/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#define HOME_VERSION
#include "all.hpp"
//#include "newmult.hpp"

#ifdef __TURBOC__
#pragma hdrstop

#endif

#ifdef __ZTC__
#include <iostream.hpp>
#endif

#ifdef __SUN__

#endif
#ifdef __NDPX__

#endif
// **************************************************************
// **************************************************************
// set to 1 to test exception handling for choleski decomp failure
int test_adexception_flag=0;
// **************************************************************
// **************************************************************

extern "C" void adfloat_except(int k);
/*
{
  throw myexception();
  //throw ad_float_exception(1);
}
*/

dmatrix choleski_decomp(_CONST dmatrix& MM)
{
  // kludge to deal with constantness
  dmatrix & M= * (dmatrix *) &MM;
  if (M.colsize() != M.rowsize())
  {
    cerr << "Error in chol_decomp. Matrix not square" << endl;
    ad_exit(1);
  }
  int rowsave=M.rowmin();
  int colsave=M.colmin();
  M.rowshift(1);
  M.colshift(1);
  int n=M.rowmax();

  dmatrix L(1,n,1,n);
#ifndef SAFE_INITIALIZE
    L.initialize();
#endif

  int i,j,k;
  MY_DOUBLE_TYPE tmp;
    if (M(1,1)<=0)
    {
      cerr << "Error matrix not positive definite in choleski_decomp"
        <<endl;
      ad_exit(1);
    }
  L(1,1)=sqrt(M(1,1));
  for (i=2;i<=n;i++)
  {
    L(i,1)=M(i,1)/L(1,1);
  }

  for (i=2;i<=n;i++)
  {
    for (j=2;j<=i-1;j++)
    {
      tmp=M(i,j);
      for (k=1;k<=j-1;k++)
      {
        tmp-=L(i,k)*L(j,k);
      }
      L(i,j)=tmp/L(j,j);
    }
    tmp=M(i,i);
    for (k=1;k<=i-1;k++)
    {
      tmp-=L(i,k)*L(i,k);
    }
  // TEST
     
    if (i==25 && test_adexception_flag) 
       tmp=-1.0;
    if (tmp<=0)
    {
      cerr << "Error matrix not positive definite in choleski_decomp"
        <<endl;

       throw myexception();

      ad_exit(1);
    }
    L(i,i)=sqrt(tmp);
  }
  L.rowshift(rowsave);
  L.colshift(colsave);

  return L;
}
static dmatrix onerror(dmatrix& L,int & ierror)
{
  ierror=1;
  return L;
}

dmatrix choleski_decomp_error(_CONST dmatrix& MM,int& ierror)
{
  ierror=0;
  // kludge to deal with constantness
  dmatrix & M= * (dmatrix *) &MM;
  if (M.colsize() != M.rowsize())
  {
    cerr << "Error in chol_decomp. Matrix not square" << endl;
    ad_exit(1);
  }
  int rowsave=M.rowmin();
  int colsave=M.colmin();
  M.rowshift(1);
  M.colshift(1);
  int n=M.rowmax();

  dmatrix L(1,n,1,n);
#ifndef SAFE_INITIALIZE
    L.initialize();
#endif

  int i,j,k;
  MY_DOUBLE_TYPE tmp;
    if (M(1,1)<=0)
    {
      return onerror(L,ierror);
    }
  L(1,1)=sqrt(M(1,1));
  for (i=2;i<=n;i++)
  {
    L(i,1)=M(i,1)/L(1,1);
  }

  for (i=2;i<=n;i++)
  {
    for (j=2;j<=i-1;j++)
    {
      tmp=M(i,j);
      for (k=1;k<=j-1;k++)
      {
        tmp-=L(i,k)*L(j,k);
      }
      L(i,j)=tmp/L(j,j);
    }
    tmp=M(i,i);
    for (k=1;k<=i-1;k++)
    {
      tmp-=L(i,k)*L(i,k);
    }
    if (tmp<=0)
    {
      return onerror(L,ierror);
    }
    L(i,i)=sqrt(tmp);
  }
  L.rowshift(rowsave);
  L.colshift(colsave);

  return L;
}

dmatrix choleski_decomp_neghess_error(_CONST dmatrix& MM,int& ierror)
{
  ierror=0;
  // kludge to deal with constantness
  dmatrix & M= * (dmatrix *) &MM;
  if (M.colsize() != M.rowsize())
  {
    cerr << "Error in chol_decomp. Matrix not square" << endl;
    ad_exit(1);
  }
  int rowsave=M.rowmin();
  int colsave=M.colmin();
  M.rowshift(1);
  M.colshift(1);
  int n=M.rowmax();

  dmatrix L(1,n,1,n);
#ifndef SAFE_INITIALIZE
    L.initialize();
#endif

  int i,j,k;
  MY_DOUBLE_TYPE tmp;
    if (M(1,1)>=0)
    {
      return onerror(L,ierror);
    }
  L(1,1)=sqrt(-M(1,1));
  for (i=2;i<=n;i++)
  {
    L(i,1)=-M(i,1)/L(1,1);
  }

  for (i=2;i<=n;i++)
  {
    for (j=2;j<=i-1;j++)
    {
      tmp=-M(i,j);
      for (k=1;k<=j-1;k++)
      {
        tmp-=L(i,k)*L(j,k);
      }
      L(i,j)=tmp/L(j,j);
    }
    tmp=-M(i,i);
    for (k=1;k<=i-1;k++)
    {
      tmp-=L(i,k)*L(i,k);
    }
    if (tmp<=0)
    {
      return onerror(L,ierror);
    }
    L(i,i)=sqrt(tmp);
  }
  L.rowshift(rowsave);
  L.colshift(colsave);

  return L;
}

#undef HOME_VERSION
