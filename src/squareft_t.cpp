/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE



#define HOME_VERSION
#include "all.hpp"

#include "variable.hpp"
#include "newmprot.hpp"


dvariable square_fit_t(dvar_len_fish_stock_history& fsh, d3_array& total_freq)
{
  dvariable var=exp(fsh.log_length_variance);
  dvariable tot_tmp=0.;
  dvariable tmp_sigma=0.;
  int nlint =fsh.nlint;
  MY_DOUBLE_TYPE eps=1.0/fsh.nlint;
  MY_DOUBLE_TYPE avp=1.0/fsh.nlint;
  if (fsh.parest_flags(193))
    eps=fsh.parest_flags(193)/(100.*nlint);
  fsh.ppstf->lencontrib=0.;
  dvariable v=exp(fsh.log_length_dof);
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

    for (int ip=1;ip<=ntimes;ip++)
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
      {
        if (fsh.len_sample_size(ir,ip,fi)>0)
        {
          dvariable len_tmp=0.;
          MY_DOUBLE_TYPE tau2=1.0/total_freq(ir,ip,fi);
          dvar_vector es=eps+elem_prod(1.0-fsh.tprob(ir,ip,fi),
                           fsh.tprob(ir,ip,fi)/avp);
          dvar_vector esquigle=pow(es,fsh.length_tot_exp);
          dvar_vector varQ=esquigle*(var*tau2);
          dvar_vector sinv=1.0/sqrt(varQ);
          dvariable ln_det_choleski_inv=sum(log(sinv));

          int pp=nlint;
          const MY_DOUBLE_TYPE lppi2=0.5*pp*log(3.1415926535);
          if (!fsh.pmsd || fsh.pmsd->fisn(ir,ip,fi)==1)
          {
            dvar_vector diff=fsh.len_freq(ir,ip,fi)-fsh.tprob(ir,ip,fi);
            dvar_vector e=elem_prod(sinv,diff);
            len_tmp = -gammln(0.5*(v+pp)) + gammln(0.5*v)
                + lppi2 +(0.5*pp)*log(v)
                +0.5*(v+pp)*log(1.0+e*e/v)
                - ln_det_choleski_inv;
          }
          else
          {
            if (ir==fsh.pmsd->reg_in_catch(ir,ip,fi,1))
            {
              int mmin=fsh.tprob(ir,ip,fi).indexmin();
              int mmax=fsh.tprob(ir,ip,fi).indexmax();
              dvar_vector tmp(mmin,mmax);
              tmp=fsh.tprob(ir,ip,fi);
              int n=fsh.pmsd->fisn(ir,ip,fi);
              for (int i=2;i<=n;i++)
              {
                int rr=fsh.pmsd->reg_in_catch(ir,ip,fi,i);
                tmp+=fsh.tprob(rr,ip,fi);
              }
              tmp/=double(n);

              dvar_vector diff=fsh.len_freq(ir,ip,fi)-tmp;
              dvar_vector e=elem_prod(sinv,diff);

              len_tmp = -gammln(0.5*(v+pp)) + gammln(0.5*v)
                  + lppi2 +(0.5*pp)*log(v)
                  +0.5*(v+pp)*log(1.0+e*e/v)
                  - ln_det_choleski_inv;
            }
          }

          //contrib by fishery ....PK jun26-07
          fsh.ppstf->lencontrib(fsh.parent(ir,ip,fi)) += value(len_tmp);
          tot_tmp += len_tmp;
        }
      }
    }
  }
  return tot_tmp;
}

dvariable square_fit_t_wght(dvar_len_fish_stock_history& fsh, d3_array& total_wght)
{
  dvariable var=exp(fsh.log_weight_variance);
  dvariable tot_tmp=0.;
  dvariable tmp_sigma=0.;
  int nwint =fsh.nwint;
  MY_DOUBLE_TYPE eps=0.1/nwint;
  MY_DOUBLE_TYPE avp=1.0/nwint;
  if (fsh.parest_flags(193))
    eps=fsh.parest_flags(193)/(100.*nwint);
  fsh.ppstf->wghtcontrib=0.;
  dvariable v=exp(fsh.log_weight_dof);
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

    for (int ip=1;ip<=ntimes;ip++)
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
      {
        if (fsh.wght_sample_size(ir,ip,fi)>0)
        {
          dvariable wght_tmp=0.;
          MY_DOUBLE_TYPE tau2=1.0/total_wght(ir,ip,fi);
          dvar_vector es=eps+elem_prod(1.0-fsh.wtprob(ir,ip,fi),
                           fsh.wtprob(ir,ip,fi)/avp);
          dvar_vector esquigle=pow(es,fsh.weight_tot_exp);
          dvar_vector varQ=esquigle*(var*tau2);
          dvar_vector sinv=1.0/sqrt(varQ);
          dvariable ln_det_choleski_inv=sum(log(sinv));

          int pp=nwint;
          const MY_DOUBLE_TYPE lppi2=0.5*pp*log(3.1415926535);
          if (!fsh.pmsd || fsh.pmsd->fisn(ir,ip,fi)==1)
          {
            dvar_vector diff=fsh.wght_freq(ir,ip,fi)-fsh.wtprob(ir,ip,fi);
            dvar_vector e=elem_prod(sinv,diff);
            wght_tmp = -gammln(0.5*(v+pp)) + gammln(0.5*v)
                + lppi2 +(0.5*pp)*log(v)
                +0.5*(v+pp)*log(1.0+e*e/v)
                - ln_det_choleski_inv;
          }
          else
          {
            if (ir==fsh.pmsd->reg_in_catch(ir,ip,fi,1))
            {
              int mmin=fsh.wtprob(ir,ip,fi).indexmin();
              int mmax=fsh.wtprob(ir,ip,fi).indexmax();
              dvar_vector tmp(mmin,mmax);
              tmp=fsh.wtprob(ir,ip,fi);
              int n=fsh.pmsd->fisn(ir,ip,fi);
              for (int i=2;i<=n;i++)
              {
                int rr=fsh.pmsd->reg_in_catch(ir,ip,fi,i);
                tmp+=fsh.wtprob(rr,ip,fi);
              }
              tmp/=double(n);

              dvar_vector diff=fsh.wght_freq(ir,ip,fi)-tmp;
              dvar_vector e=elem_prod(sinv,diff);

              wght_tmp = -gammln(0.5*(v+pp)) + gammln(0.5*v)
                  + lppi2 +(0.5*pp)*log(v)
                  +0.5*(v+pp)*log(1.0+e*e/v)
                  - ln_det_choleski_inv;
            }
          }

          //contrib by fishery ....PK jun26-07
          fsh.ppstf->wghtcontrib(fsh.parent(ir,ip,fi)) += value(wght_tmp);
          tot_tmp += wght_tmp;
        }
      }
    }
  }
  return tot_tmp;
}

#undef HOME_VERSION
