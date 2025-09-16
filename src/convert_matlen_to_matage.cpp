/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

dvar_vector age_at_length_calc(int nlint,MY_DOUBLE_TYPE shlen,MY_DOUBLE_TYPE filen,
  const dvar_vector& vb_coff,int nage,const ivector& pf,int mult)
{
  dvar_vector agedist(1,mult*nlint);
  agedist.initialize();
  MY_DOUBLE_TYPE minv=1.0/mult;
  if (pf(226)==0)
  {
    dvariable rho=exp(-vb_coff(3));
    dvariable d=1.0/(vb_coff(2)-vb_coff(1));
    dvariable x = vb_coff(1)*d;
    dvariable y = 1.0-pow(rho,nage-1);

    MY_DOUBLE_TYPE v=shlen-0.5*filen*minv;
    //v=0.999*value(vb_coff(2))-filen;
    cout << "Linf = " 
#if !defined(NO_MY_DOUBLE_TYPE)
         << vb_coff(1)+(vb_coff(2)-vb_coff(1))/(1-pow(rho,nage-1.0L)) << endl;
#else
         << vb_coff(1)+(vb_coff(2)-vb_coff(1))/(1-pow(rho,nage-1.0)) << endl;
#endif
    for (int ill=1;ill<=mult*nlint;ill++)
    {
      //int il=(ill-1)/mult+1;
      v+=filen*minv;
      //dvariable tmpx= (v-vb_coff(1))/(vb_coff(2)-vb_coff(1));
      dvariable tmp= v*d-x;
      //dvariable tmp1x= 1.-tmpx*(1.0-pow(rho,nage-1));
      dvariable tmp1= 1.-tmp*y;
      //cout << tmp1 << " " << tmp1x << endl;
      if (tmp1<=0.0) tmp1=1.e-20;
      dvariable age= 1.-log(tmp1)/vb_coff(3);
      if (age<1) age=1;
      if (age>nage) age=nage;
      agedist(ill)=age;
    }
    return agedist;
  }
  else
  {
    dvariable T=0.0;
    if (pf(226)==1)
      T=-exp(vb_coff(4));
    else
      T=exp(vb_coff(4));
    dvariable c1=pow(vb_coff(1),-1.0/T);
    dvariable cN=pow(vb_coff(2),-1.0/T);
    dvariable rho=exp(-vb_coff(3));
    MY_DOUBLE_TYPE v=shlen-0.5*minv*filen;
    dvariable y = 1.0-pow(rho,nage-1);
    for (int ill=1;ill<=mult*nlint;ill++)
    {
      v+=minv*filen;
      dvariable pv=pow(v,-1.0/T);
      dvariable tmp= (pv-c1)/(cN-c1);
      dvariable tmp1= 1.-tmp*y;
      if (tmp1<=0.0) tmp1=1.e-20;
      dvariable age= 1.-log(tmp1)/vb_coff(3);
      if (age<1) age=1;
      if (age>nage) age=nage;
      agedist(ill)=age;
    }
    return agedist;
  }
}

dvar_vector maturity_length_to_age(dvar_vector& alc,
  dvector& cpmature_at_length,int nage,int nlint,int mult)
{
  dvar_vector mat_at_age(1,nage);
  dvar_vector prop_at_age(1,nage);
  mat_at_age.initialize();
  prop_at_age.initialize();
  for (int ill=1;ill<=mult*nlint;ill++)
  {
    int il=(ill-1)/mult+1;
    MY_DOUBLE_TYPE cage=value(alc(ill));

    int jj=int(cage);

    dvariable sf=daves_kludge1(alc(ill));
    dvariable sf1=1.0-sf;
    mat_at_age(jj)+=sf1*cpmature_at_length(il);
    prop_at_age(jj)+=sf1;
    if (jj<nage)
    {
      mat_at_age(jj+1)+=sf*cpmature_at_length(il);
      prop_at_age(jj+1)+=sf;
    }
    else
    {
      mat_at_age(jj)+=sf*cpmature_at_length(il);
      prop_at_age(jj)+=sf;
    }
  }  
  mat_at_age=elem_div(mat_at_age,1.e-10+prop_at_age);
   
  return mat_at_age;
}
