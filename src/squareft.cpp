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

dvariable objective_function(dvar_fish_stock_history& fsh)
{
  dvariable xxxx=0.0;

  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    xxxx+=norm2(log(.001+fsh.obs_catch(ir))-fsh.catch(ir)); // catch is 
                                              // already "logged"
    xxxx+=norm2(log(fsh.tot_catch(ir))-log(fsh.obs_tot_catch(ir)));
  }
  return xxxx;
}

int done_it_already=0;

// Modified min. chi^2


dvariable square_fit0(dvar_len_fish_stock_history& fsh, d3_array& total_freq)
{
  dvariable len_tmp=0.;
  dvariable tmp_sigma=0.;
  int nlint =fsh.nlint;
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
          if (fsh.parest_flags(311)==0)
          {
            len_tmp-=total_freq(ir,ip,fi)*sum(log(exp(-square(elem_div(
              (fsh.len_freq(ir,ip,fi)-fsh.tprob(ir,ip,fi)),
              sqr(.01+fsh.len_freq(ir,ip,fi)))))+.001));
            if (fsh.parest_flags(161)==0)
            {
              len_tmp+=0.5*sum(log(.01+fsh.len_freq(ir,ip,fi)));
            }
          }
          else
          {
            dvar_vector & tp =fsh.tprob(ir,ip,fi);
            len_tmp-=total_freq(ir,ip,fi)*sum(log(exp(-square(elem_div(
              (fsh.tc_len_freq(ir,ip,fi)-fsh.tc_tprob(ir,ip,fi)),
              sqr(.01+fsh.tc_len_freq(ir,ip,fi)))))+.001));
            if (fsh.parest_flags(161)==0)
            {
              len_tmp+=0.5*sum(log(.01+fsh.tc_len_freq(ir,ip,fi)));
            }
          }
        }
      }
    }
  }
  return len_tmp;
}

dvariable square_fit0a(dvar_len_fish_stock_history& fsh, d3_array& total_freq)
{
  dvariable len_tmp=0.;
  dvariable tmp_sigma=0.;
  int nlint =fsh.nlint;
  MY_DOUBLE_TYPE wt=0.1/nlint;
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
          dvar_vector r2=0.5*elem_div(
            square(fsh.len_freq(ir,ip,fi)-fsh.tprob(ir,ip,fi)),
              wt+fsh.len_freq(ir,ip,fi));

          len_tmp-=total_freq(ir,ip,fi)*
            sum(log(exp(-r2) + 0.001));

          if (fsh.parest_flags(161)==0)
          {
            len_tmp+=0.5*sum(log(wt+fsh.len_freq(ir,ip,fi)));
          }
        }
      }
    }
  }
  return len_tmp;
}



// put int the factor of 2 and use proper distribution

dvariable square_fit(dvar_len_fish_stock_history& fsh, d3_array& total_freq)
{
  dvariable len_tmp=0.;
  dvariable tmp_sigma=0.;
  int nlint =fsh.nlint;
  MY_DOUBLE_TYPE wt=0.001/nlint;
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
          dvar_vector r2=0.5*elem_div(
            square(fsh.len_freq(ir,ip,fi)-fsh.tprob(ir,ip,fi)),
              wt+fsh.tprob(ir,ip,fi));

          len_tmp-=total_freq(ir,ip,fi)*
            sum(log(exp(-r2) + 0.03/(1.0+r2)));

          if (fsh.parest_flags(161)==0)
          {
            len_tmp+=0.5*sum(log(wt+fsh.tprob(ir,ip,fi)));
          }
        }
      }
    }
  }
  return len_tmp;
}

// LH function as in YFT paper

dvariable square_fita(dvar_len_fish_stock_history& fsh, d3_array& total_freq)
{
  dvariable tot_tmp=0.;
  dvariable tmp_sigma=0.;
  int nlint =fsh.nlint;
  //double wt=0.1/nlint;
  MY_DOUBLE_TYPE wt=1.0/nlint;
  fsh.ppstf->lencontrib=0.;
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
          dvar_vector esquigle=elem_prod(1.0-fsh.len_freq(ir,ip,fi),
                           fsh.len_freq(ir,ip,fi));
          dvar_vector varQ=esquigle+wt;
          if (!fsh.pmsd || fsh.pmsd->fisn(ir,ip,fi)==1)
          {
            dvar_vector r2=elem_div(
              square(fsh.len_freq(ir,ip,fi)-fsh.tprob(ir,ip,fi)),
              2*varQ*tau2);
            len_tmp -= sum(log(exp(-r2) + 0.001));
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
              dvar_vector r2=elem_div(
                square(fsh.len_freq(ir,ip,fi)-tmp),
                2*varQ*tau2);
              len_tmp -= sum(log(exp(-r2) + 0.001));
            }
          }

            

          if (fsh.parest_flags(161)==0)
          {
            len_tmp+=0.5*sum(log(6.283186*varQ));
            len_tmp+=nlint*log(sqrt(tau2));
          }
          //contrib by fishery ....PK jun26-07
          fsh.ppstf->lencontrib(fsh.parent(ir,ip,fi)) += value(len_tmp);
          fsh.ppstf->lencontrib_by_realization(ir,ip,fi) = value(len_tmp);
          tot_tmp += len_tmp;
        }
      }
    }
  }
  return tot_tmp;
}

dvariable square_fita_wght(dvar_len_fish_stock_history& fsh,
  d3_array& total_wght)
{
  dvariable tot_tmp=0.;
  dvariable tmp_sigma=0.;
  int nwint =fsh.nwint;
  fsh.ppstf->wghtcontrib=0.;
  MY_DOUBLE_TYPE wt=1.0/nwint;
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
          dvar_vector esquigle=elem_prod(1.0-fsh.wght_freq(ir,ip,fi),
                           fsh.wght_freq(ir,ip,fi));
          dvar_vector varQ=esquigle+wt;
          if (!fsh.pmsd || fsh.pmsd->fisn(ir,ip,fi)==1)
          {
            dvar_vector r2=elem_div(
              square(fsh.wght_freq(ir,ip,fi)-fsh.wtprob(ir,ip,fi)),
              2*varQ*tau2);
            wght_tmp-= sum(log(exp(-r2) + 0.001));
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
              dvar_vector r2=elem_div(
                square(fsh.wght_freq(ir,ip,fi)-tmp),
                2*varQ*tau2);
              wght_tmp -= sum(log(exp(-r2) + 0.001));
            }
          }

          if (fsh.parest_flags(161)==0)
          {
            wght_tmp+=0.5*sum(log(6.283186*varQ));
            wght_tmp+=nwint*log(sqrt(tau2));
          }
          //contrib by fishery ....PK jun26-07
          fsh.ppstf->wghtcontrib(fsh.parent(ir,ip,fi)) += value(wght_tmp);
          fsh.ppstf->wghtcontrib_by_realization(ir,ip,fi) = value(wght_tmp);
          tot_tmp += wght_tmp;
        }
      }
    }
  }
  return tot_tmp;
}
#undef HOME_VERSION
