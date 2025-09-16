/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

  
#define HOME_VERSION
#include "fvar.hpp"


void set_value(const dvar_vector& _x,const dvector& y,const dvar_vector& v, const int& _ii)
{
  dvar_vector& x=(dvar_vector&) _x;
  int& ii=(int&) _ii;
  int min=x.indexmin();
  int max=x.indexmax();
  for (int i=min;i<=max;i++)
  {
    if (y(i)) x(i)=v(ii++);
  }
}

void set_value(const dvar_vector& _x,const dvector& y,const dvar_vector& v, const int& _ii, 
  CGNU_DOUBLE fmin,CGNU_DOUBLE fmax,const dvariable& fpen)
{
  dvar_vector& x=(dvar_vector&) _x;
  int& ii=(int&) _ii;
  int min=x.indexmin();
  int max=x.indexmax();
  for (int i=min;i<=max;i++)
  {
    if (y(i)) x(i)=boundp(v(ii++),fmin,fmax,fpen);
  }
}


void set_value_inv(const dvector& x,const dvector& y ,const dvector& _v, const int& _ii)
{
  dvector& v=(dvector&) _v;
  int& ii=(int&) _ii;
  int min=x.indexmin();
  int max=x.indexmax();
  for (int i=min;i<=max;i++)
  {
    if (y(i)) v(ii++)=x(i);
  }
}

void set_value_inv(const dvector& x,const dvector& y,const dvector& _v, const int& _ii, 
  CGNU_DOUBLE fmin,CGNU_DOUBLE fmax)
{
  dvector& v=(dvector&) _v;
  int& ii=(int&) _ii;
  int min=x.indexmin();
  int max=x.indexmax();
  for (int i=min;i<=max;i++)
  {
    if (y(i)) v(ii++)=boundpin(x(i),fmin,fmax);
  }
}

#undef HOME_VERSION
