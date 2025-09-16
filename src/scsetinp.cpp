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


void set_value_inv_partial(const dvar_vector& x,const dvector& _v, const int& _ii, int n,
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, MY_DOUBLE_TYPE s)
{
  dvector& v=(dvector&)_v;
  int& ii=(int&)_ii;
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
    v(ii++)=s*boundpin(x(i),fmin,fmax);
  }
}


void set_value_inv_partial(const dvector& x,const dvector& _v, const int& _ii, int n,
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, MY_DOUBLE_TYPE s)
{
  dvector& v=(dvector&)_v;
  int& ii=(int&)_ii;
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
    v(ii++)=s*boundpin(x(i),fmin,fmax);
  }
}

void set_value_inv_partial(const dvar_matrix& x,const dvector& _v, const int& _ii, int n,
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, MY_DOUBLE_TYPE s)
{
  dvector& v=(dvector&)_v;
  int& ii=(int&)_ii;
  int min=x.indexmin();
  int max=x.indexmax();
  for (int i=min;i<=max;i++)
  {
    set_value_inv_partial(x(i),v,ii,n,fmin,fmax,s);
    //set_value_inverse_partial(x(i),v,ii,n,fmin,fmax,s);
  }
}

/*
void set_value_inv_partial(const dvar_vector& x,const dvector& _v, const int& _ii, int n,
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax)
{
  dvector& v=(dvector&)_v;
  int& ii=(int&)_ii;
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
    v(ii++)=boundpin(x(i),fmin,fmax);
  }
}
*/

void set_value_inv_partial(const dvar_vector& x,const dvector& _v, const int& _ii, int n,
  MY_DOUBLE_TYPE s)
{
  dvector& v=(dvector&)_v;
  int& ii=(int&)_ii;
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
    v(ii++)=s*value(x(i));
  }
}


#undef HOME_VERSION
