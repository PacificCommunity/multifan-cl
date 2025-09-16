/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"

par_ofstream& operator << (par_ofstream& pof,int x)
{
  (ofstream&)pof << x; 
  return pof;
}
par_ofstream& operator << (par_ofstream& pof,const char * s)
{
  (ofstream&)pof << s; 
  return pof;
}
 par_ofstream& operator << (par_ofstream& pof,const prevariable& v)
{
  if (v < 1.e-60 && v > -1.e-60) 
  {
    (ofstream&)pof << " 0 "; 
  }
  else
  {
    (ofstream&)pof << v; 
  }
  return pof;
}
 par_ofstream& operator << (par_ofstream& pof,const dvar_vector& v)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    pof << " " << v[i];
  }
  pof << endl;
  return pof;
}
 par_ofstream& operator << (par_ofstream& pof,const dvar_matrix& v)
{
  int mmin=v.rowmin();
  int mmax=v.rowmax();
  for (int i=mmin;i<=mmax;i++)
  {
    if (allocated(v[i]))
      pof << v[i];
  }
  return pof;
}
 par_ofstream& operator << (par_ofstream& pof,const dvar3_array& v)
{
  if (allocated(v))
  {
    int mmin=v.slicemin();
    int mmax=v.slicemax();
    for (int i=mmin;i<=mmax;i++)
    {
      pof << v[i];
    }
  }
  return pof;
}
 par_ofstream& operator << (par_ofstream& pof,const dvar4_array& v)
{
  int mmin=v.hslicemin();
  int mmax=v.hslicemax();
  for (int i=mmin;i<=mmax;i++)
  {
    pof << v[i];
  }
  return pof;
}
 par_ofstream& operator << (par_ofstream& pof,MY_DOUBLE_TYPE v)
{
  if (v < 1.e-60 && v > -1.e-60) 
  {
    (ofstream&)pof << " 0 "; 
  }
  else
  {
    (ofstream&)pof << v << " "; 
  }
  return pof;
}
 par_ofstream& operator << (par_ofstream& pof,const dvector& v)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    pof << " " << v[i];
  }
  pof << endl;
  return pof;
}
 par_ofstream& operator << (par_ofstream& pof,const ivector& v)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    pof << " " << v[i];
  }
  pof << endl;
  return pof;
}
 par_ofstream& operator << (par_ofstream& pof,const dmatrix& v)
{
  int mmin=v.rowmin();
  int mmax=v.rowmax();
  for (int i=mmin;i<=mmax;i++)
  {
    pof << v[i];
  }
  return pof;
}
 par_ofstream& operator << (par_ofstream& pof,const d3_array& v)
{
  int mmin=v.slicemin();
  int mmax=v.slicemax();
  for (int i=mmin;i<=mmax;i++)
  {
    pof << v[i];
  }
  return pof;
}
 par_ofstream& operator << (par_ofstream& pof,const d4_array& v)
{
  int mmin=v.hslicemin();
  int mmax=v.hslicemax();
  for (int i=mmin;i<=mmax;i++)
  {
    pof << v[i];
  }
  return pof;
}
