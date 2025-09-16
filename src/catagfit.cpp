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

dvariable catch_at_age_fit(dvar_len_fish_stock_history& fsh,int print_switch)
{
  dvariable f=0.0;
  if (fsh.age_nage==fsh.nage)
  {
    for (int ir=1;ir<=fsh.num_regions;ir++)
    {
      for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)
      {
        for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
        {                                        
          if (fsh.age_sample_size(ir,ip,fi)>0)
          {
            dvar_vector tmp_catch=exp(fsh.catch(ir,ip,fi));
            tmp_catch/=sum(tmp_catch);
            f+=100.0*norm2(tmp_catch-fsh.age_freq(ir,ip,fi));
          }
        }
      }
    }
  }
  else if (fsh.age_nage<fsh.nage)
  {
    for (int ir=1;ir<=fsh.num_regions;ir++)
    {
      for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)
      {
        for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
        {                                        
          if (fsh.age_sample_size(ir,ip,fi)>0)
          {
            dvar_vector tmp_catch=exp(fsh.catch(ir,ip,fi));
            tmp_catch/=sum(tmp_catch);
            tmp_catch(fsh.age_nage)+=sum(tmp_catch(fsh.age_nage+1,fsh.nage));
            f+=100.0*norm2(tmp_catch(1,fsh.age_nage)-fsh.age_freq(ir,ip,fi));
          }
        }
      }
    }
  }
  else if (fsh.age_nage>fsh.nage)
  {
    dvector tmp_age_freq(1,fsh.nage);
    for (int ir=1;ir<=fsh.num_regions;ir++)
    {
      for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)
      {
        for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
        {                                        
          if (fsh.age_sample_size(ir,ip,fi)>0)
          {
            dvector& freq=fsh.age_freq(ir,ip,fi);
            dvar_vector tmp_catch=exp(fsh.catch(ir,ip,fi));
            tmp_catch/=sum(tmp_catch);
            tmp_age_freq=freq(1,fsh.nage);
            tmp_age_freq(fsh.nage)
              +=sum(fsh.age_freq(fsh.nage+1,fsh.age_nage));
            f+=100.0*norm2(tmp_catch-tmp_age_freq);
          }
        }
      }
    }
  }
  return f;
}



dvariable catch_at_age_fit_dirichlet_multinomial(dvar_len_fish_stock_history& fsh,int print_switch)
{
  dvariable f=0.0;
  dvar_vector A=mfexp(fsh.fish_pars(13));
  if (fsh.age_nage==fsh.nage)
  {
    for (int ir=1;ir<=fsh.num_regions;ir++)
    {
      for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)
      {
        for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
        {                                        
          if (fsh.age_sample_size(ir,ip,fi)>0)
          {
            dvar_vector tmp_catch=exp(fsh.catch(ir,ip,fi));
            tmp_catch/=sum(tmp_catch);
            tmp_catch+=1.e-8;
            //f+=100.0*norm2(tmp_catch-fsh.age_freq(ir,ip,fi));

            int fishery=fsh.parent(ir,ip,fi);
  
            dvar_vector alpha=A(fishery)*tmp_catch;
  
            f-=gammln(A(fishery))
             + sum(gammln(alpha+fsh.age_freq(ir,ip,fi)))
             -gammln(A(fishery)+fsh.age_sample_size(ir,ip,fi))
             -sum(gammln(alpha));

          }
        }
      }
    }
  }
  else if (fsh.age_nage<fsh.nage)
  {
    cerr << "not implemented yet" << endl;
    ad_exit(1);
    for (int ir=1;ir<=fsh.num_regions;ir++)
    {
      for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)
      {
        for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
        {                                        
          if (fsh.age_sample_size(ir,ip,fi)>0)
          {
            dvar_vector tmp_catch=exp(fsh.catch(ir,ip,fi));
            tmp_catch/=sum(tmp_catch);
            tmp_catch(fsh.age_nage)+=sum(tmp_catch(fsh.age_nage+1,fsh.nage));
            f+=100.0*norm2(tmp_catch(1,fsh.age_nage)-fsh.age_freq(ir,ip,fi));
          }
        }
      }
    }
  }
  else if (fsh.age_nage>fsh.nage)
  {
    cerr << "not implemented yet" << endl;
    ad_exit(1);
    dvector tmp_age_freq(1,fsh.nage);
    for (int ir=1;ir<=fsh.num_regions;ir++)
    {
      for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)
      {
        for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
        {                                        
          if (fsh.age_sample_size(ir,ip,fi)>0)
          {
            dvector& freq=fsh.age_freq(ir,ip,fi);
            dvar_vector tmp_catch=exp(fsh.catch(ir,ip,fi));
            tmp_catch/=sum(tmp_catch);
            tmp_age_freq=freq(1,fsh.nage);
            tmp_age_freq(fsh.nage)
              +=sum(fsh.age_freq(fsh.nage+1,fsh.age_nage));
            f+=100.0*norm2(tmp_catch-tmp_age_freq);
          }
        }
      }
    }
  }
  return f;
}

dvariable catch_at_age_fit_multinomial_logit(dvar_len_fish_stock_history& fsh,int print_switch)
{
  dvariable f=0.0;
  dvar_vector A=mfexp(fsh.fish_pars(13));
  if (fsh.age_nage==fsh.nage)
  {
    for (int ir=1;ir<=fsh.num_regions;ir++)
    {
      for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)
      {
        for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
        {                                        
          if (fsh.age_sample_size(ir,ip,fi)>0)
          {
            dvar_vector& logcatch=fsh.catch(ir,ip,fi);
            dvar_vector logit=logcatch(1,fsh.nage-1)-logcatch(fsh.nage);
            int fishery=fsh.parent(ir,ip,fi);
          }
        }
      }
    }
  }
  else if (fsh.age_nage<fsh.nage)
  {
    cerr << "not implemented yet" << endl;
    ad_exit(1);
    for (int ir=1;ir<=fsh.num_regions;ir++)
    {
      for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)
      {
        for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
        {                                        
          if (fsh.age_sample_size(ir,ip,fi)>0)
          {
            dvar_vector tmp_catch=exp(fsh.catch(ir,ip,fi));
            tmp_catch/=sum(tmp_catch);
            tmp_catch(fsh.age_nage)+=sum(tmp_catch(fsh.age_nage+1,fsh.nage));
            f+=100.0*norm2(tmp_catch(1,fsh.age_nage)-fsh.age_freq(ir,ip,fi));
          }
        }
      }
    }
  }
  else if (fsh.age_nage>fsh.nage)
  {
    cerr << "not implemented yet" << endl;
    ad_exit(1);
    dvector tmp_age_freq(1,fsh.nage);
    for (int ir=1;ir<=fsh.num_regions;ir++)
    {
      for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)
      {
        for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
        {                                        
          if (fsh.age_sample_size(ir,ip,fi)>0)
          {
            dvector& freq=fsh.age_freq(ir,ip,fi);
            dvar_vector tmp_catch=exp(fsh.catch(ir,ip,fi));
            tmp_catch/=sum(tmp_catch);
            tmp_age_freq=freq(1,fsh.nage);
            tmp_age_freq(fsh.nage)
              +=sum(fsh.age_freq(fsh.nage+1,fsh.age_nage));
            f+=100.0*norm2(tmp_catch-tmp_age_freq);
          }
        }
      }
    }
  }
  return f;
}

//this was the old catach at age fit code
/*
dvariable catch_at_age_fit(dvar_len_fish_stock_history& fsh,int print_switch)
{
  d3_array total_freq=total_freq_calc(fsh,print_switch);
  dvar4_array tmp_catch=exp(fsh.catch);
  d4_array tmp_obs_catch(1,fsh.num_regions,1,fsh.num_fish_periods,
    1,fsh.num_fish_incidents,1,fsh.nage);
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)  // Loop over fishing periods
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        for (int j=1;j<=fsh.nage;j++)
        {
          tmp_obs_catch(ir,ip,fi,j)=fsh.len_freq(ir,ip,fi,j);
        }
        tmp_catch(ir,ip,fi)=tmp_catch(ir,ip,fi)/sum(tmp_catch(ir,ip,fi));
        tmp_obs_catch(ir,ip,fi)=tmp_obs_catch(ir,ip,fi)/
          sum(tmp_obs_catch(ir,ip,fi));
      }
    } 
  }
  dvariable len_tmp=0.;
  for (ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)  // Loop over fishing periods
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        len_tmp-=total_freq(ir,ip,fi)*sum(log(exp(-square(elem_div(
          (tmp_obs_catch(ir,ip,fi)-tmp_catch(ir,ip,fi)),
        sqr(.01+tmp_catch(ir,ip,fi)))))+.0001));
      }
    }
  }
  cout << "leaving catch-age fit" << endl;
  return len_tmp;
}
*/

#undef HOME_VERSION
