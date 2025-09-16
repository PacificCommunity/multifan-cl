/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


  

#include "all.hpp"
#include "scbd.hpp"


void set_value(const dvar_vector& _x,const dvector& y,const dvar_vector& v, const int& _ii, 
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, const dvariable& fpen,MY_DOUBLE_TYPE s)
{
  dvar_vector& x=(dvar_vector&) _x;
  int& ii=(int&) _ii;
  int min=x.indexmin();
  int max=x.indexmax();
  for (int i=min;i<=max;i++)
  {
    if (y(i)) x(i)=boundp(v(ii++),fmin,fmax,fpen,s);
  }
}

void set_value(const dvar_vector& _x,const dvar_vector& v,int flag,
  const int& _ii,MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, const dvariable& fpen,MY_DOUBLE_TYPE s)
{
  dvar_vector& x=(dvar_vector&) _x;
  int& ii=(int&) _ii;
  int min=x.indexmin();
  int max=x.indexmax();
  if (flag)
  {
    for (int i=min;i<=max;i++)
    {
      x(i)=boundp(v(ii++),fmin,fmax,fpen,s);
    } 
  }
}

void set_value(const dvar_vector& _x,const dvar_vector& v,
  const int& _ii,MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax, const dvariable& fpen,MY_DOUBLE_TYPE s)
{
  dvar_vector& x=(dvar_vector&) _x;
  int& ii=(int&) _ii;
  int min=x.indexmin();
  int max=x.indexmax();
  for (int i=min;i<=max;i++)
  {
    x(i)=boundp(v(ii++),fmin,fmax,fpen,s);
  } 
}

void set_value_inv(const dvar_vector& _x,const dvar_vector& v,int flag,
  const int& _ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s)
{
  dvar_vector& x=(dvar_vector&) _x;
  int& ii=(int&) _ii;
  int min=x.indexmin();
  int max=x.indexmax();
  if (flag)
  {
    for (int i=min;i<=max;i++)
    {
      x(i)=boundpin(v(ii++),fmin,fmax,s);
    } 
  }
}

void set_value_inv(const dvar_vector& _x,const dvar_vector& v,
  const int& _ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s)
{
  dvar_vector& x=(dvar_vector&) _x;
  int& ii=(int&) _ii;
  int min=x.indexmin();
  int max=x.indexmax();
  for (int i=min;i<=max;i++)
  {
    x(i)=boundpin(v(ii++),fmin,fmax,s);
  } 
}



#undef HOME_VERSION
