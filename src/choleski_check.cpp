/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
/**
 * \file
 * Description not yet available.
 */
#include "fvar.hpp"

#ifdef __TURBOC__
#pragma hdrstop
#include <iostream.h>
#endif

#ifdef __ZTC__
#include <iostream.hpp>
#endif

#ifdef __SUN__
#include <iostream.h>
#endif
#ifdef __NDPX__
#include <iostream.h>
#endif

/**
 * Description not yet available.
 * \param
 */
int choleski_check(const dmatrix& MM,MY_DOUBLE_TYPE bound)
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

  int retval=0;
  int i,j,k;
  MY_DOUBLE_TYPE tmp;
    if (M(1,1)<=bound)
    {
      M(1,1)=bound;
      retval=1;
      return retval;
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
    if (tmp<=bound)
    {
      retval=1;
      return retval;
    }
    L(i,i)=sqrt(tmp);
  }
  return retval;
}

