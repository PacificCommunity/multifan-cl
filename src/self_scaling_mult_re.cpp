/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

dvariable len_self_scaling_multinomial_re(dvar_len_fish_stock_history& fsh)
{
 
  dvariable rho=fsh.length_rho;
  dvariable psi=fsh.length_psi;
  dvariable var=exp(fsh.log_length_variance);
  int nlint=fsh.nlint;
  
  dvar_matrix vS(1,nlint,1,nlint);

  // correlation for random effects
  dvar_vector rp(0,nlint-1);
  rp(0)=1.0;
  switch (fsh.parest_flags(297))
  {
  case 0:
    rp(1)=rho;  //LN1  LN2
    break;
  case 1:   // LN3m
    rp(1)=rho+psi/(1.0+square(rho+psi)/(1.0-square(rho)));
    break;
  case 2:   // LN3
    {
      dvariable phi2=-1+(2.0-sfabs(rho))*psi;
      rp(1)=rho/(1.0-phi2);
    }
    break;
  default:
    cerr << "illegal value for fsh.parest_flags(297)" << endl;
    ad_exit(1);
  }
  switch (fsh.parest_flags(297))
  {
  case 0:
  case 1:   // LN3m
    for (int i=2;i<nlint;i++)
    {
      rp(i)=rho*rp(i-1);
    }
    break;
  case 2:   // LN3
    {
      dvariable phi2=-1+(2.0-sfabs(rho))*psi;
      for (int i=2;i<nlint;i++)
      {
        rp(i)=rho*rp(i-1)+phi2*rp(i-2);
      }
    }
    break;
  default:
    cerr << "illegal value for fsh.parest_flags(297)" << endl;
    ad_exit(1);
  }
  // multiply by variance
  rp*=var;
  for (int i=1;i<=nlint;i++)
  {
    vS(i,i)=1.0;
    for (int j=1;j<i;j++)
    {
      vS(i,j)=rp(i-j);
      vS(j,i)=rp(i-j);
    }
  }
  dvar_matrix ch=choleski_decomp(vS);
  dvar_matrix vSinv=inv(vS);
  dvariable ld=0.5*ln_det(vS);
  
  fsh.lSinv=value(vSinv);
  const MY_DOUBLE_TYPE a=0.75;
  const MY_DOUBLE_TYPE beta1=0.99500;
  const MY_DOUBLE_TYPE beta2=0.9990;
  const int n=nlint;
  MY_DOUBLE_TYPE offset=(n-1)/10.;
  if (n>20) offset+=(n-20)/10.;
  if (n>40) offset+=(n-40)/20.;
  if (n>60) offset+=(n-60)/20.;

  dvar_vector vN=exp(fsh.fish_pars(14));
  dvar_vector avn;
  if (n<50)
    avn=pow(vN,a);
  else if (n<75)
    avn=pow(vN,0.9*a);
  else
    avn=pow(vN,0.65*a);

#if !defined(NO_MY_DOUBLE_TYPE)
  dvar_vector alpha=beta1+(beta2-beta1)*elem_div(avn-1.0L,avn);
#else
  dvar_vector alpha=beta1+(beta2-beta1)*elem_div(avn-1.0,avn);
#endif
  dvar_vector svn=pow(vN,alpha);


  dvariable len_tmp=0.;
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    int ntimes;
    if (fsh.parest_flags(142)==0)
    {
      ntimes=fsh.num_fish_periods(ir);
    }
    else 
    {
      ntimes=min(fsh.parest_flags(142),fsh.num_fish_periods(ir));
    }
    funnel_dvariable ftmp;
    //dvariable ftmp;
    dvariable end_tmp=0.0;
    int blocksize=40;
    //for (int ipp=1;ipp<=10;ipp+=blocksize)
    for (int ipp=1;ipp<=ntimes;ipp+=blocksize)
    {
      dvariable loc_tmp=0.0;
      ad_begin_funnel();
      cout << " ipp= " << ipp << endl;
      for (int ip=ipp;ip<=ipp+blocksize-1;ip++)
      {
        if (ip>ntimes) break;
        for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
        {
          int pp=fsh.parent(ir,ip,fi);
          dvector OP=fsh.len_freq(ir,ip,fi);
          dvar_vector P=fsh.tprob(ir,ip,fi);
          if (fsh.len_sample_size(ir,ip,fi)>0)
          {
            loc_tmp+=ld;
            loc_tmp-=gammln(svn(pp)+1.0+offset);
            dvar_vector sobs=vN(pp)*OP;
            for (int ii=1;ii<=n;ii++)
            {
              if (value(sobs(ii))>0.0)
              loc_tmp+=gammln(sobs(ii));
            }
            ad_exit(1);
           /*
            loc_tmp+=lognormal_multinomial_log_likelihood(vN(pp),OP,P,
              vSinv,fsh.lSinv);
           */
          }
        }
      }
      ftmp=loc_tmp;
      end_tmp=ftmp;
      len_tmp+=end_tmp;
    }
  }
  return len_tmp;
}



dvariable wt_self_scaling_multinomial_re(dvar_len_fish_stock_history& fsh)
{
  dvariable rho=fsh.weight_rho;
  dvariable psi=fsh.weight_psi;
  dvariable var=exp(fsh.log_weight_variance);
  int nwint=fsh.nwint;
  
  dvar_matrix vS(1,nwint,1,nwint);

  // correlation for random effects
  dvar_vector rp(0,nwint-1);
  rp(0)=1.0;
  switch (fsh.parest_flags(297))
  {
  case 0:
    rp(1)=rho;  //LN1  LN2
    break;
  case 1:   // LN3m
    rp(1)=rho+psi/(1.0+square(rho+psi)/(1.0-square(rho)));
    break;
  case 2:   // LN3
    {
      dvariable phi2=-1+(2.0-sfabs(rho))*psi;
      rp(1)=rho/(1.0-phi2);
    }
    break;
  default:
    cerr << "illegal value for fsh.parest_flags(297)" << endl;
    ad_exit(1);
  }
  switch (fsh.parest_flags(297))
  {
  case 0:
  case 1:   // LN3m
    for (int i=2;i<nwint;i++)
    {
      rp(i)=rho*rp(i-1);
    }
    break;
  case 2:   // LN3
    {
      dvariable phi2=-1+(2.0-sfabs(rho))*psi;
      for (int i=2;i<nwint;i++)
      {
        rp(i)=rho*rp(i-1)+phi2*rp(i-2);
      }
    }
    break;
  default:
    cerr << "illegal value for fsh.parest_flags(297)" << endl;
    ad_exit(1);
  }
  // multiply by variance
  rp*=var;
  for (int i=1;i<=nwint;i++)
  {
    vS(i,i)=1.0;
    for (int j=1;j<i;j++)
    {
      vS(i,j)=rp(i-j);
      vS(j,i)=rp(i-j);
    }
  }
  dvar_matrix ch=choleski_decomp(vS);
  dvar_matrix vSinv=inv(vS);
  dvariable ld=0.5*ln_det(vS);
  fsh.wtSinv=value(vSinv);
  const MY_DOUBLE_TYPE a=0.75;
  const MY_DOUBLE_TYPE beta1=0.99500;
  const MY_DOUBLE_TYPE beta2=0.9990;
  const int n=nwint;
  MY_DOUBLE_TYPE offset=(n-1)/10.;
  if (n>20) offset+=(n-20)/10.;
  if (n>40) offset+=(n-40)/20.;
  if (n>60) offset+=(n-60)/20.;

  dvar_vector vN=exp(fsh.fish_pars(15));
  dvar_vector avn;
  if (n<50)
    avn=pow(vN,a);
  else if (n<75)
    avn=pow(vN,0.9*a);
  else
    avn=pow(vN,0.65*a);

#if !defined(NO_MY_DOUBLE_TYPE)
  dvar_vector alpha=beta1+(beta2-beta1)*elem_div(avn-1.0L,avn);
#else
  dvar_vector alpha=beta1+(beta2-beta1)*elem_div(avn-1.0,avn);
#endif
  dvar_vector svn=pow(vN,alpha);


  dvariable len_tmp=0.;
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    int ntimes;
    if (fsh.parest_flags(142)==0)
    {
      ntimes=fsh.num_fish_periods(ir);
    }
    else 
    {
      ntimes=min(fsh.parest_flags(142),fsh.num_fish_periods(ir));
    }
    funnel_dvariable ftmp;
    //dvariable ftmp;
    dvariable end_tmp=0.0;
    int blocksize=40;
    //for (int ipp=1;ipp<=10;ipp+=blocksize)
    for (int ipp=1;ipp<=ntimes;ipp+=blocksize)
    {
      dvariable loc_tmp=0.0;
      ad_begin_funnel();
      cout << " ipp= " << ipp << endl;
      for (int ip=ipp;ip<=ipp+blocksize-1;ip++)
      {
        if (ip>ntimes) break;
        for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
        {
          int pp=fsh.parent(ir,ip,fi);
          dvector OP=fsh.wght_freq(ir,ip,fi);
          dvar_vector P=fsh.wtprob(ir,ip,fi);
          if (fsh.wght_sample_size(ir,ip,fi)>0)
          {
            loc_tmp+=ld;
            loc_tmp-=gammln(svn(pp)+1.0+offset);
            dvar_vector sobs=vN(pp)*OP;
            for (int ii=1;ii<=n;ii++)
            {
              if (value(sobs(ii))>0.0)
              loc_tmp+=gammln(sobs(ii));
            }
            ad_exit(1);
           /*
            loc_tmp+=lognormal_multinomial_log_likelihood(vN(pp),OP,P,
              vSinv,fsh.wtSinv);
           */
          }
        }
      }
      ftmp=loc_tmp;
      end_tmp=ftmp;
      len_tmp+=end_tmp;
    }
  }
  return len_tmp;
}
