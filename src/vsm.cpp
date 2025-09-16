/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

dvariable sfabs1(const prevariable & x)
{
  const MY_DOUBLE_TYPE eps=1.e-7;
  return sqrt(eps+square(x));
}

dvariable smax1(const prevariable & x,const prevariable & y)
{
  return 0.5*(x+y+0.1*sfabs1((x-y)*10.));
}

dvariable smax1(const dvar_vector& v)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  dvariable tmp=smax1(v(mmin),v(mmin+1));
  for (int i=mmin+2;i<=mmax-1; i++)
  {
    tmp=smax1(tmp,v(i));
  } 
  return tmp;
}
