/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


  

#include "all.hpp"
#include "scbd.hpp"

void set_value_partial(const dvar_vector& _x,const dvar_vector& v, const int& _ii, int n,
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, const dvariable& fpen, MY_DOUBLE_TYPE s)
{
  dvar_vector& x=(dvar_vector&) _x;
  int& ii=(int&) _ii;
  int min=x.indexmin();
  int max=min+n-1;
  #ifdef SAFE_ARRAYS
    if (max >x.indexmax())
    {
      cerr << "index out of range in set_value_patial(const dvar_vector&, ... "
           << endl;
    }
  #endif
  for (int i=min;i<=max;i++)
  {
    x(i)=boundp(v(ii++)/s,fmin,fmax,fpen);
  }
}

void set_value_partial(const dvar_matrix _x,const dvar_vector& v,
  const int& _ii, int n,MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, const dvariable& fpen,
  MY_DOUBLE_TYPE s)
{
  ADUNCONST(dvar_matrix,x)
  ADUNCONST(int,ii)
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    set_value_partial(x(i),v,ii,n,fmin,fmax,fpen,s);
  }
}

void set_value_partial(const dvar_vector& _x,const dvar_vector& v, const int& _ii, int n,
  MY_DOUBLE_TYPE s)
{
  dvar_vector& x=(dvar_vector&) _x;
  int& ii=(int&) _ii;
  int min=x.indexmin();
  int max=min+n-1;
  #ifdef SAFE_ARRAYS
    if (max >x.indexmax())
    {
      cerr << "index out of range in set_value_patial(const dvar_vector&, ... "
           << endl;
    }
  #endif
  for (int i=min;i<=max;i++)
  {
    x(i)=v(ii++)/s;
  }
}


#undef HOME_VERSION
