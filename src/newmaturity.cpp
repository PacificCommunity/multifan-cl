/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
dvector setm11(const dvector & v)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  dvector tmp=v/(v(mmax)-v(mmin))*2;
  tmp=tmp-tmp(mmin)-1.0;
  return tmp;
}

dvar_vector setm11(const dvar_vector & v)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  dvar_vector tmp=v/(v(mmax)-v(mmin))*2;
  tmp=tmp-tmp(mmin)-1.0;
  return tmp;
}


dvar_vector max(const dvar_vector & v,MY_DOUBLE_TYPE d)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  dvar_vector tmp(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    if (v(i)>d)
     tmp(i)=v(i);
    else
     tmp(i)=d;
  }
  return tmp;
}

dvar_vector min(const dvar_vector & v,MY_DOUBLE_TYPE d)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  dvar_vector tmp(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    if (v(i)<d)
     tmp(i)=v(i);
    else
     tmp(i)=d;
  }
  return tmp;
}

dvar_vector dvar_len_fish_stock_history::maturity_length_to_age_simple_spline
  (int isp)
{
//  vcubic_spline_function csf(fmid,cpmature_at_length);
  MY_DOUBLE_TYPE eps=1.e-6;
  dvar_vector ml(1,nage);
  dvar_vector tmp(1,nage);  //NMD_19jan2022
  tmp.initialize();
  dvariable fpen=0.0;
  dvar_vector vbc;
//  dvector cpml;   NMD_14Dec2022
  dvector cpml(cpmature_at_length.indexmin(),cpmature_at_length.indexmax());    //NMD_19Sep2018, 14Dec2022
  cpml.initialize();
  if (isp==1)
  {
    vbc=vb_coff;
    cpml=cpmature_at_length;    //NMD_19Sep2018
//    cpml+=1.e-8;       //NMD_14Dec2022
  }
  else
  {
    vbc=pmsd->vb_coff(isp);
    cpml=value(pmsd->cpmature_at_length(isp));    //NMD_19Sep2018
//    cpml+=1.e-8;       //NMD_14Dec2022
  }
  vcubic_spline_function csf(fmid,cpml);    //NMD_19Sep2018

  dvariable rho=exp(-vbc(3));
  dvariable xn1inv=1.0/(1.-pow(rho,nage-1));
  if (parest_flags(226)==0)
  {
    dvariable diff=vbc(2)-vbc(1);
    for (int j=1;j<=nage;j++)
    {
      ml(j)=vbc(1)+diff*(1.-pow(rho,j-1))*xn1inv;
    }
  }
  else
  {
    dvariable T;
    if (parest_flags(226)==1)
      T=exp(vbc(4));   // change sign to make simpler
    else
      T=-exp(vbc(4));
    dvariable c1=pow(vbc(1),1.0/T);
    dvariable cN=pow(vbc(2),1.0/T);
    dvariable diff=cN-c1;
    for (int j=1;j<=nage;j++)
    {
      dvariable tt=c1+diff*(1.-pow(rho,j-1))*xn1inv;
      ml(j)=pow(tt,T);
    }
  }

//  dvar_vector matage=csf(ml);
  dvar_vector matage(1,nage);
  matage.initialize();
  tmp=csf(ml);  //NMD_19jan2022
  for (int j=1;j<=nage;j++)
  {
    tmp(j)=posfun(tmp(j),0.001,fpen);
    if (tmp(j)<0.0)
    {
      cout << "Error: negative maturity at age from spline function" << endl;
      cout <<  " j = " << j << "  " << tmp(j) << endl;
      cout << "tmp(j)" << endl;
      cout << tmp(j) << endl;
      ad_exit(1);
    }
    matage(j)=tmp(j);
  }  
  dvar_vector tmp3= matage/smax1(matage);
  tmp3=-posfun(1.0-tmp3,0.001,fpen)+1.0;
  if (min(tmp3)<0.0 )
  {
    cout << " min case matage " << matage <<  endl;
    cout << " smax1(matage) " << smax1(matage) <<  endl;
    cout << " tmp3 " << tmp3 <<  endl;
    ad_exit(1);
  }
  if (max(tmp3)>1.0 )
  {
    cout << " max case matage " << matage <<  endl;
    cout << " smax1(matage) " << smax1(matage) <<  endl;
    cout << " tmp3 " << tmp3 <<  endl;
    ad_exit(1);
  }
  matage=tmp3;
  return matage;
}
dvar_vector dvar_len_fish_stock_history::maturity_length_to_age_weighted_spline
  (int isp)
{
//  vcubic_spline_function csf(fmid,cpmature_at_length);  //NMD_19Sep2018
  dvar_vector ml(1,nage);
  dvar_vector sigma(1,nage);
  dvar_vector vbc;
  dvar_vector vrc;
//  dvector cpml;    //NMD_19Sep2018
  dvector cpml(cpmature_at_length.indexmin(),cpmature_at_length.indexmax());    //NMD_19Sep2018, 14Dec2022
  cpml.initialize();
  if (isp==1)
  {
    vbc=vb_coff;
    vrc=var_coff;
    cpml=cpmature_at_length;    //NMD_19Sep2018
//    cpml+=1.e-8;       //NMD_14Dec2022
  }
  else
  {
    vbc=pmsd->vb_coff(isp);
    vrc=pmsd->var_coff(isp);
    cpml=value(pmsd->cpmature_at_length(isp));    //NMD_19Sep2018
//    cpml+=1.e-8;       //NMD_14Dec2022
  }
//  cpml+=1.e-08;       //NMD_14Dec2022
  vcubic_spline_function csf(fmid,cpml);    //NMD_19Sep2018
  dvariable rho=exp(-vbc(3));
  dvariable temp2=1.-pow(rho,nage-1);
  dvariable xn1inv=1.0/temp2;
  
  if (parest_flags(226)==0)
  {
    dvariable diff=vbc(2)-vbc(1);
    for (int j=1;j<=nage;j++)
    {
      dvariable tmp=(1.-pow(rho,j-1))*xn1inv;
      sigma(j)=vrc(1)*exp(vrc(2)*(-1.0+2.0*tmp));
      ml(j)=vbc(1)+diff*tmp;
    }
  }
  else
  {
    dvariable T;
    if (parest_flags(226)==1)
      T=exp(vbc(4));   // change sign to make simpler
    else
      T=-exp(vbc(4));
    dvariable c1=pow(vbc(1),1.0/T);
    dvariable cN=pow(vbc(2),1.0/T);
    dvariable diff=cN-c1;
    for (int j=1;j<=nage;j++)
    {
      dvariable tmp=(1.-pow(rho,j-1))*xn1inv;
      dvariable tt=c1+diff*tmp;
      ml(j)=pow(tt,T);
    }
    dvar_vector scaled_ml=setm11(ml);   //scale the ml between -1 and 1
    sigma=vrc(1)*exp(vrc(2)*scaled_ml);
  }
#if !defined(NO_MY_DOUBLE_TYPE)
  const MY_DOUBLE_TYPE v[]={-2.0,-1.5,-1.0L,-0.5,0.0,0.5,1.0L,1.5,2.0};
#else
  const MY_DOUBLE_TYPE v[]={-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0};
#endif
  const MY_DOUBLE_TYPE cw0=1;
#if !defined(NO_MY_DOUBLE_TYPE)
  const MY_DOUBLE_TYPE cw1=exp(-.5*.5*.5L); // 0.5
#else
  const MY_DOUBLE_TYPE cw1=exp(-.5*.5*.5); // 0.5
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
  const MY_DOUBLE_TYPE cw2=exp(-.5L);  // 1.0
#else
  const MY_DOUBLE_TYPE cw2=exp(-.5);  // 1.0
#endif
  const MY_DOUBLE_TYPE cw3=exp(-1.125); // 1.5
  const MY_DOUBLE_TYPE cw4=exp(-2.0);  // 2.0
  const MY_DOUBLE_TYPE cwsum=cw0+2*cw1+2*cw2+2*cw3+2*cw4;
  const MY_DOUBLE_TYPE w0=cw0/cwsum;
  const MY_DOUBLE_TYPE w1=cw1/cwsum;
  const MY_DOUBLE_TYPE w2=cw2/cwsum;
  const MY_DOUBLE_TYPE w3=cw3/cwsum;
  const MY_DOUBLE_TYPE w4=cw4/cwsum;
  const MY_DOUBLE_TYPE w[]={w4,w3,w2,w1,w0,w1,w2,w3,w4};
  const MY_DOUBLE_TYPE * pw=w+4;
  const MY_DOUBLE_TYPE * pv=v+4;
  dvar_matrix tmp(-4,4,1,nage);
  dvar_vector matage(1,nage);
  matage.initialize();
  MY_DOUBLE_TYPE lmin=fmid(1);
  MY_DOUBLE_TYPE lmax=fmid(nlint);
  cout << lmin << " " << lmax << endl;
  dvar_vector tmp1;
  //dvar_vector tmpx1;
  dvar_vector tmp2;
  //dvar_vector tmpx2;
  //dvar_vector tmpx;
  dvariable fpen=0.0;
  for (int i=-4;i<=4;i++)
  {
    //dvar_vector tmp1=max(ml+pv[i]*sigma,lmin);
    //dvar_vector tmp2=min(tmp1,lmax);
    tmp1=max(ml+pv[i]*sigma,lmin);
    tmp2=min(tmp1,lmax);
    tmp(i)=csf(tmp2);
    tmp(i)=posfun(tmp(i),0.001,fpen);
    if (min(tmp(i))<0.0)
    {
      cout <<  " i = " << i << "  " << min(tmp(i)) << endl;
      cout << "tmp(i)" << endl;
      cout << tmp(i) << endl;
      cout << "tmp2" << endl;
      cout << tmp2 << endl;
      cout << tmp1 << endl;
      ad_exit(1);
    }
    matage+=pw[i]*tmp(i);
  }
  dvar_vector tmp3= matage/smax1(matage);
  tmp3=-posfun(1.0-tmp3,0.001,fpen)+1.0;
  if (min(tmp3)<0.0 )
  {
    cout << " min case matage " << matage <<  endl;
    cout << " smax1(matage) " << smax1(matage) <<  endl;
    cout << " tmp3 " << tmp3 <<  endl;
    ad_exit(1);
  }
  if (max(tmp3)>1.0 )
  {
    cout << " max case matage " << matage <<  endl;
    cout << " smax1(matage) " << smax1(matage) <<  endl;
    cout << " tmp3 " << tmp3 <<  endl;
    ad_exit(1);
  }
  //matage/=smax1(matage);
  matage=tmp3;
  return matage;
}


/*
dvar_vector smooth_lengths(dvar_vector& len,dvector & cpmature_at_length)
{
  int mmin=len.indexmin();
  int mmax=len.indexmax();
  
    if (age<1)
    {
      age=1;
    }
    else if (age>nage)
    {
      age=nage;
    }
    else if (age<1.1)
    {
      dvariable u=age-1.0;
      u=  u*u*( 20. - 100. *u);
      age= 1.+u;
    }
    else if (age>nage-.1)
    {
      dvariable u=nage-age;
      u=  u*u*( 20. - 100. *u);
      age=nage-u;
    }
    return age;
  }
*/

