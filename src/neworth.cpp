/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#define USE_DD_NOT
#include <fvar.hpp>
//#include <qd/fpu.h>
#undef qdmatrix;
//qd_real sum(const qdvector&);

//void test_fpu(void);

dmatrix orthpoly_constant_begin_end(int n,int deg,int begin_degree,
  int nconst_begin,int end_degree,int nconst_end)
{
  unsigned int old_cw;
  //fpu_fix_start(&old_cw);
  unsigned int old_cw2;
  //fpu_fix_start(&old_cw2);
  MY_DOUBLE_TYPE qx;
  qx=100;
  int i; int j; int ia; int is; int ik;
  dmatrix ocoff(0,deg,1,n);
  ocoff.initialize();

  MY_DOUBLE_TYPE ssum;
  ocoff(0)=sqrt(double(n));
  if (nconst_begin==0) nconst_begin=1;
  if (nconst_end==0) nconst_end=1;
  if (nconst_begin>n-1)
  {
    cerr << "nconst_begin too large in orthpoly_constant_begin"
         << endl;
    ad_exit(1);
  }
  if (deg>n-nconst_begin-nconst_end+1)
  {
    cerr << "deg too large in orthpoly_constant_begin"
         << endl;
    ad_exit(1);
  }
  if (nconst_begin>n-nconst_end)
  {
    cerr << "This sucks" << endl;
    ad_exit(1);
  }
  
  MY_DOUBLE_TYPE pi=3.14159265358979323844;
  int icount=1;
  for (is=1; is<=deg; is++)
  {
#if !defined(NO_MY_DOUBLE_TYPE)
    ocoff(is,1)=pow(-1.0L,double(is));
#else
    ocoff(is,1)=pow(-1.0,double(is));
#endif
  }
  for (is=1; is<=deg; is++)
  {
#if !defined(NO_MY_DOUBLE_TYPE)
    ocoff(is,1)=pow(-1.0L,double(is));
#else
    ocoff(is,1)=pow(-1.0,double(is));
#endif
    int begin_lag=0;
    for (j=2; j<=n; j++)
    {
      int jj=j;
      if (j<=nconst_begin  && is>=begin_degree)
      {
        jj=1;
        begin_lag++;
      }
      else if (j<=n-nconst_end+1 || is<end_degree)
      {
        jj=j-begin_lag;
      }
      else if (j>n-nconst_end+1 && is>=end_degree)
      {
        jj=n-nconst_end+1-begin_lag;
      }
      MY_DOUBLE_TYPE maxj=n-nconst_end+1-begin_lag;
      MY_DOUBLE_TYPE qdis=is;
      MY_DOUBLE_TYPE qdjj=jj;
      MY_DOUBLE_TYPE qdone=1.0;
      MY_DOUBLE_TYPE qdtwo=2.0;
      icount++;
      ocoff(is,j)=pow(-qdone+(qdjj-qdone)/(maxj-qdone)*qdtwo,is);
    }
  }
  
  for (is=0; is<=deg; is++) 
  {
    ocoff(is)/=norm(ocoff(is));
    // for extra robustness
    for (ik=0; ik<is; ik++)
    {
      ssum=ocoff(is)*ocoff(ik);
      ocoff(is)-=ssum*ocoff(ik);
    }
    ocoff(is)/=norm(ocoff(is));
    for (ik=is+1; ik<=deg; ik++)
    {
      ssum=ocoff(is)*ocoff(ik);
      ocoff(ik)-=ssum*ocoff(is);
    }
  }
  
#if defined(USE_DD)
  dmatrix ttmp=make_dmatrix(ocoff);
#else
  dmatrix ttmp=ocoff;
#endif
  return trans(ttmp);
}

