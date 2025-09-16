/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#define HOME_VERSION
#include "fvar.hpp"
#ifdef __TURBOC__
  #pragma hdrstop
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

MY_DOUBLE_TYPE dmin(double,MY_DOUBLE_TYPE);
MY_DOUBLE_TYPE dmax(double, MY_DOUBLE_TYPE);

  dvariable boundp( dvariable xx, MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, prevariable& fpen,
    _CONST MY_DOUBLE_TYPE& s)
{
  dvariable t,y,x;
  x=xx/s;

  t=fmin + (fmax-fmin)*(sin(x*1.570795)+1.)/2.;

  if (x < -.9999)
  {
    fpen+=(x+0.9999)*(x+0.9999);
  }

  if (x > 0.9999)
  {
    fpen+=(x-0.9999)*(x-0.9999);
  }

  if (x < -1.)
  {
    fpen+=1000.*(x+1.)*(x+1.);
  }

  if (x > 1.)
  {
    fpen+=1000.*(x-1.)*(x-1.);
  }

  return(t);
}

MY_DOUBLE_TYPE boundp( MY_DOUBLE_TYPE xx, MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, MY_DOUBLE_TYPE& fpen,
  _CONST MY_DOUBLE_TYPE& s)
{
  MY_DOUBLE_TYPE t,x;
  x=xx/s;
  t=fmin + (fmax-fmin)*(sin(x*1.570795)+1)/2;

  if (x < -.9999)
  {
    fpen+=(x+0.9999)*(x+0.9999);
  }

  if (x > 0.9999)
  {
    fpen+=(x-0.9999)*(x-0.9999);
  }

  if (x < -1)
  {
    fpen+=1000*(x+1)*(x+1);
  }

  if (x > 1)
  {
    fpen+=1000*(x-1)*(x-1);
  }

  return(t);
}

MY_DOUBLE_TYPE boundpin(MY_DOUBLE_TYPE x, MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax,_CONST MY_DOUBLE_TYPE& s)
{
  MY_DOUBLE_TYPE tinv;

  if (x < fmin)
  {
    printf("variable out of bounds in boundpin: variable = %lg", x);
    printf("; min = %lg", fmin);
    printf("; max = %lg\n", fmax);

    x=dmin(fmin+.001,fmin+.01*(fmax-fmin));
  }

  if (x > fmax)
  {
    printf("variable out of bounds in boundpin: variable = %lg", x);
    printf("; min = %lg", fmin);
    printf("; max = %lg\n", fmax);

    x=dmax(fmax-.001,fmax-.01*(fmax-fmin));
  }

  tinv=asin(2.*(x-fmin)/(fmax-fmin)-1.)/1.570795;
  return(s*tinv);
}

#undef HOME_VERSION
