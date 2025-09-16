/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

dvariable len_multiplier_calcs (dvar_len_fish_stock_history& fsh,
  int ir,int ip,int fi,int pp,int gp,int sz,const dvector& average_sample_size,
  const dvector& OP)
{
  dvariable len_tmp=0.0;
  dvariable cvN;
  dvar_vector fp20=fsh.fish_pars(20);
  MY_DOUBLE_TYPE ss=fsh.len_sample_size(ir,ip,fi);
  dvariable vN=exp(fsh.fish_pars(14,pp));
  dvariable size_multiplier=
     pow(ss/average_sample_size(gp),fp20(pp));
  dvariable xx=vN*size_multiplier;
  dvariable fpen=0.0;
  xx=posfun(xx-0.04,.01,fpen)+0.04;
  MY_DOUBLE_TYPE ub=1000.;
  if (fsh.parest_flags(336))
    ub=fsh.parest_flags(336);
  cvN=-posfun(-xx+ub,.01,fpen)+ub;
  dvariable ln=log(cvN);
  const int nslots=sz;
  dvector pp1=fsh.m_estimator_coffs(nslots);
  int ioff=1;
  int nss=static_cast<int>(pp1(1+ioff));
  ioff++;
  dvector sample_sizes(1,nss);
  sample_sizes=pp1(1+ioff,nss+ioff).shift(1);
  ioff+=nss;
  MY_DOUBLE_TYPE ftoff=pp1(1+ioff);
  ioff++;
  dvector aa(1,7);
  aa=pp1(1+ioff,7+ioff).shift(1);
  MY_DOUBLE_TYPE mvv=mean(log(sample_sizes));
  const MY_DOUBLE_TYPE l1000=log(1000.);
  dvariable t1=(ln-l1000)/l1000;
  dvariable tmp=exp((aa(6)+aa(7)*t1)*ln);
  // KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
  len_tmp-=
    gammln(exp(ln)+tmp+exp(aa(1)+t1*(aa(2)+t1*(aa(3)+t1*(aa(4)+t1*aa(5))))));
    // before it had no1 added
    //gammln(tmp+exp(aa(1)+t1*(aa(2)+t1*(aa(3)+t1*(aa(4)+t1*aa(5))))));
  // KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
  dvar_vector sobs=cvN*OP;
  for (int ii=1;ii<=sobs.indexmax();ii++)
  {
    if (value(sobs(ii))>0.0)
    len_tmp+=gammln(sobs(ii)+ftoff);
  }
  return len_tmp;
}
dvariable wght_multiplier_calcs (dvar_len_fish_stock_history& fsh,
  int ir,int ip,int fi,int pp,int gp,int sz,const dvector& average_sample_size,
  const dvector& OP)
{
  dvariable wght_tmp=0.0;
  dvariable cvN;
  dvar_vector fp21=fsh.fish_pars(21);
  MY_DOUBLE_TYPE ss=fsh.wght_sample_size(ir,ip,fi);
  dvariable vN=exp(fsh.fish_pars(15,pp));
  dvariable size_multiplier=
     pow(ss/average_sample_size(gp),fp21(pp));
  dvariable xx=vN*size_multiplier;
  dvariable fpen=0.0;
  xx=posfun(xx-0.04,.01,fpen)+0.04;
  MY_DOUBLE_TYPE ub=1000.;
  if (fsh.parest_flags(337))
    ub=fsh.parest_flags(337);
  cvN=-posfun(-xx+ub,.01,fpen)+ub;
  dvariable ln=log(cvN);
  const int nslots=sz;
  dvector pp1=fsh.m_estimator_coffs(nslots);
  int ioff=1;
  int nss=static_cast<int>(pp1(1+ioff));
  ioff++;
  dvector sample_sizes(1,nss);
  sample_sizes=pp1(1+ioff,nss+ioff).shift(1);
  ioff+=nss;
  MY_DOUBLE_TYPE ftoff=pp1(1+ioff);
  ioff++;
  dvector aa(1,7);
  aa=pp1(1+ioff,7+ioff).shift(1);
  MY_DOUBLE_TYPE mvv=mean(log(sample_sizes));
  const MY_DOUBLE_TYPE l1000=log(1000.);
  dvariable t1=(ln-l1000)/l1000;
  dvariable tmp=exp((aa(6)+aa(7)*t1)*ln);
  // KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
  wght_tmp-=
    gammln(exp(ln)+tmp+exp(aa(1)+t1*(aa(2)+t1*(aa(3)+t1*(aa(4)+t1*aa(5))))));
    // before it had no1 added
    //gammln(tmp+exp(aa(1)+t1*(aa(2)+t1*(aa(3)+t1*(aa(4)+t1*aa(5))))));
  // KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
  dvar_vector sobs=cvN*OP;
  for (int ii=1;ii<=sobs.indexmax();ii++)
  {
    if (value(sobs(ii))>0.0)
    wght_tmp+=gammln(sobs(ii)+ftoff);
  }
  return wght_tmp;
}
