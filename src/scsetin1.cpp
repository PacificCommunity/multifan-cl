/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE



#define HOME_VERSION
#include "fvar.hpp"
#include "scbd.hpp"


/*
void set_value_inv(const dvar_vector& x,const dvector& y,const dvector& _v,const int& _ii,
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s)
{
  dvar_vector& v= (dvar_vector&) _v;
  int& ii= (int&) _ii;
  int min=x.indexmin();
  int max=x.indexmax();
  for (int i=min;i<=max;i++)
  {
    if (y(i)) v(ii++)=boundpin(x(i),fmin,fmax,s);
  }
}
*/

#undef HOME_VERSION
