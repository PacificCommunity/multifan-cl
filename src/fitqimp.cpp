/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

/*
dvariable dvar_fish_stock_history::fit_implicit_q(void)
{
  int i;
  dvariable fpen=0.0;
  dvar_matrix loge_by_fishery(1,num_fisheries,1,num_real_fish_times);
  dvar_matrix log_q_implicit_by_fishery(1,num_fisheries,1,num_real_fish_times);
  for (i=1;i<=num_fisheries;i++)
  {
    loge_by_fishery(i)
      =log(1.e-15+effort_by_fishery(i)(1,num_real_fish_times(i)));
    loge_by_fishery(i)-=mean(loge_by_fishery(i));
  }
  for (i=1;i<=num_fisheries;i++)
  {
    for (int nt=1;nt<=num_real_fish_times(i);nt++)
    {
      int rr=realization_region(i,nt);
      int rp=realization_period(i,nt);
      int ri=realization_incident(i,nt);
      log_q_implicit_by_fishery(i,nt)=log(1.e-15+fm_level(rr,rp,ri));
    }
    log_q_implicit_by_fishery(i)
      -= mean(log_q_implicit_by_fishery(i));
  }
  MY_DOUBLE_TYPE penwt=age_flags(160);
  for (i=1;i<=num_fisheries;i++)
  {
    //fpen+=penwt*norm2(loge_by_fishery(i)-log_q_implicit_by_fishery(i));
    fpen-=
      log(.01+exp(-penwt
        *norm2(loge_by_fishery(i)-log_q_implicit_by_fishery(i))));
  }
  return fpen;
}
*/  
 
dvariable dvar_fish_stock_history::fit_implicit_q(void)
{
  int i;
  dvariable tmppen=0.0;
  dvar_matrix loge_by_fishery(1,num_fisheries,1,num_real_fish_times);
  dvar_matrix log_q_implicit_by_fishery(1,num_fisheries,1,num_real_fish_times);
  for (i=1;i<=num_fisheries;i++)
  {
    for (int nt=1;nt<=num_real_fish_times(i);nt++)
    {
      int rr=realization_region(i,nt);
      int rp=realization_period(i,nt);
      int ri=realization_incident(i,nt);
      log_q_implicit_by_fishery(i,nt)=log(1.e-15+fm_level(rr,rp,ri));
    }
  }
  ivector ff60=column(fish_flags,60);
  if (sum(ff60)==0)   // no grouping for catchability
  {
    for (i=1;i<=num_fisheries;i++)
    {
      loge_by_fishery(i)
        =log(1.e-15+effort_by_fishery(i)(1,num_real_fish_times(i)));
      loge_by_fishery(i)-=mean(loge_by_fishery(i));
    }

    for (i=1;i<=num_fisheries;i++)
    {
      log_q_implicit_by_fishery(i)
        -= mean(log_q_implicit_by_fishery(i));
    }
  }
  else
  {
    int ngroups=gfish_index.indexmax();
    dvar_vector grouped_e_means(1,ngroups);
    dvar_vector grouped_q_means(1,ngroups);
    dvector numbers(1,ngroups);
    grouped_q_means.initialize();
    grouped_e_means.initialize();
    numbers.initialize();
    for (i=1;i<=num_fisheries;i++)
    {
      loge_by_fishery(i)
        =log(1.e-15+effort_by_fishery(i)(1,num_real_fish_times(i)));
      grouped_e_means(ff60(i))+=sum(loge_by_fishery(i));
      numbers(ff60(i))+=num_real_fish_times(i);
    }

    for (i=1;i<=num_fisheries;i++)
    {
      grouped_q_means(ff60(i))+=
        sum(log_q_implicit_by_fishery(i));
    }
    for (i=1;i<=ngroups;i++)
    {
      grouped_q_means(i)/=numbers(i);
      grouped_e_means(i)/=numbers(i);
    }
    for (i=1;i<=num_fisheries;i++)
    {
      loge_by_fishery(i)
        -= grouped_e_means(ff60(i));
    }

    for (i=1;i<=num_fisheries;i++)
    {
      log_q_implicit_by_fishery(i)
        -= grouped_q_means(ff60(i));
    }
  }
  for (int fi=1;fi<=num_fisheries;fi++)
  { 
    if (fish_flags(fi,65)==0)
    {
      int nrft=num_real_fish_times(fi);
      if (fish_flags(fi,13)==0)
      {
        MY_DOUBLE_TYPE eff_wt=10.;
        dvar_vector tmp_square=
          square(loge_by_fishery(fi)-log_q_implicit_by_fishery(fi));
        //dvar_vector tmp_square=square(effort_dev_coffs(fi)(1,nrft));
        dvariable eff_dev_pen=-sum(log(1.e-10+exp(-eff_wt*tmp_square)+
          .05*exp(-eff_wt/5.*tmp_square)));
        tmppen+=eff_dev_pen;
      }
      else if (fish_flags(fi,13)>0)
      {
        MY_DOUBLE_TYPE eff_wt=fish_flags(fi,13);
        dvar_vector tmp_square=
          square(loge_by_fishery(fi)-log_q_implicit_by_fishery(fi));
        //dvar_vector tmp_square=square(effort_dev_coffs(fi)(1,nrft));
        dvariable eff_dev_pen=-sum(log(1.e-10+exp(-eff_wt*tmp_square)+
          .05*exp(-eff_wt/5.*tmp_square)));
        tmppen+=eff_dev_pen;
      }
      else if (fish_flags(fi,13)<0)
      {
        dvar_vector eff_wt=-fish_flags(fi,13)*
          sqrt(.01+effort_by_fishery(fi)(1,nrft));
        dvar_vector tmp_square=
          square(loge_by_fishery(fi)-log_q_implicit_by_fishery(fi));
        //dvar_vector tmp_square=square(effort_dev_coffs(fi)(1,nrft));
        dvariable eff_dev_pen=
          -sum(log(1.e-10+exp(elem_prod(-eff_wt,tmp_square))+
          .05*exp(elem_prod(-eff_wt/5,tmp_square))));
        tmppen+=eff_dev_pen;
      }
    }
    else
    {
      int nrft=num_real_fish_times(fi);
      if (fish_flags(fi,13)==0)
      {
        MY_DOUBLE_TYPE eff_wt=10.;
        dvar_vector tmp_square=
          square(loge_by_fishery(fi)-log_q_implicit_by_fishery(fi));
        //dvar_vector tmp_square=square(effort_dev_coffs(fi)(1,nrft));
        dvariable eff_dev_pen=-sum(log(1.e-10+exp(-eff_wt*tmp_square)+
          .01/(1.0+eff_wt/3.*tmp_square)));
        tmppen+=eff_dev_pen;
      }
      else if (fish_flags(fi,13)>0)
      {
        MY_DOUBLE_TYPE eff_wt=fish_flags(fi,13);
        dvar_vector tmp_square=
          square(loge_by_fishery(fi)-log_q_implicit_by_fishery(fi));
        //dvar_vector tmp_square=square(effort_dev_coffs(fi)(1,nrft));
        dvariable eff_dev_pen=-sum(log(1.e-10+exp(-eff_wt*tmp_square)+
          .01/(1.0+eff_wt/3.*tmp_square)));
        tmppen+=eff_dev_pen;
      }
      else if (fish_flags(fi,13)<0)
      {
        dvar_vector eff_wt=-fish_flags(fi,13)*
          sqrt(.01+effort_by_fishery(fi)(1,nrft));
        dvar_vector tmp_square=
          square(loge_by_fishery(fi)-log_q_implicit_by_fishery(fi));
        //dvar_vector tmp_square=square(effort_dev_coffs(fi)(1,nrft));
        dvariable eff_dev_pen=
           -sum(log(1.e-10+exp(elem_prod(-eff_wt,tmp_square))+
          .01/(1.0+elem_prod(eff_wt/5,tmp_square))));
        tmppen+=eff_dev_pen;
      }
    }
  }
  return tmppen;
}

