/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

void dvar_len_fish_stock_history::
  do_popes_approximation(int ir,int ip,dvariable& ffpen,dvariable& ffpen1)
{
  MY_DOUBLE_TYPE pwght1=100.0;
  MY_DOUBLE_TYPE pwght2=100.0;
  dvar_vector M=e_nat_mort_by_period(ir,ip);
  dvar_vector wt=mean_weight(ir,ip,1);
  int nfi=num_fish_incidents(ir,ip);
  dvar_matrix sel(1,nfi,1,nage);
  dvar_vector lambda(1,nfi);
  dvar_matrix F(1,nfi,1,nage);
  ivector wtflag=column(data_fish_flags,1);
  dvar_matrix * pfm=0;
  //dvar_matrix * pN=0;
  dvar_vector * ptm=0;
  dvar_vector * psurv=0;
  dvar_vector * pnumfish=0;
  dvar3_array * uuu1=0;
  dvar_vector * uuu=0;
  if (af170q0==0)
  {
    pnumfish=&(num_fish(ir,ip));
    uuu=&num_fish(ir,ip);
    uuu1=&num_fish;
    ptm=&tot_mort(ir,ip);
    //pN=&N;
    pfm=&fish_mort(ir,ip);
    psurv=&(survival)(ir,ip);
    cout << pnumfish << " " << uuu << " " << &((*uuu1)(ir,ip)) << endl;
  }
  else
  {
    pnumfish=&num_fish_q0(ir,ip);
    ptm=&tot_mort_q0(ir,ip);
    //pN=&N_q0;
    pfm=&fish_mort_q0(ir,ip);
    psurv=&(survival_q0)(ir,ip);
  }
  const dvar_vector& enf=mfexp(*pnumfish);
  dvar_vector NM2=elem_prod(enf,-0.5*M);

  for (int fi=1;fi<=nfi;fi++)
  {
    sel(fi)=mfexp(incident_sel(ir,ip,fi));
    if (missing_catch_for_incident_flag(ir,ip,fi)==0)
    {
      if (wtflag(fi)==0)
      {
        lambda(fi)=obs_tot_catch(ir,ip,fi)/(sel(fi)*NM2);
      }
      else
      {
        lambda(fi)=obs_tot_catch(ir,ip,fi)/(wt*elem_prod(sel(fi),NM2));
      }
      F(fi)=lambda(fi)*sel(fi);
    }
    else
    {
      int i=parent(ir,ip,fi);
      dvector dv=fml_designvec(ir,ip,fi);
      ivector num_ifmlrp=implicit_fml_bounds(3);
      dvar_vector ests=implicit_fm_level_regression_pars(i)(1,num_ifmlrp(i));
      // reorganize log effort by fishery so that this is simple
      MY_DOUBLE_TYPE leff=log_effort_by_fishery(i,fish_times(ir,ip,fi));
      dvariable tmp2 = dv * inv(fml_R(i)) * ests;
      fm_level(ir,ip,fi)=tmp2+leff;
      lambda(fi)=exp(fm_level(ir,ip,fi));
      F(fi)=lambda(fi)*sel(fi);
    }
  }
  //for (int fi=1;fi<=nfi;fi++)
  //{
    // impose bounds constraints
    dvar_vector omega(1,nage);
    dvar_vector Z=colsum(F);
    for (int j=1;j<=nage;j++)
    { 
      if (Z(j) < pamin1)
      {
        (*psurv)(j)=1-Z(j);
      }
      else if (Z(j) < pamin2)
      {
        (*psurv)(j)=1-Z(j);
        ffpen+=pwght1*square(Z(j)-pamin1);
      }
      else
      {
        ffpen+=pwght1*square(Z(j)-pamin1);
        dvariable ptmp=0.0;
        dvariable tmp=pamin2-posfun(pamin2-Z(j),0.02,ptmp);
        ffpen1+=pwght2*ptmp;
        omega(j)=tmp/Z(j);
        Z(j)=tmp;
        (*psurv)(j)=1-Z(j);
        for (int fi=1;fi<=nfi;fi++)
        {
          F(fi,j)*=omega(j);
        }
      }
    } 
  //}
}
  /* 
  for (fi=1;fi<=nfi;fi++)
  {
    // get selectivity and fish mort 
    if (missing_catch_for_incident_flag(ir,ip,fi)==0)
    {
      sel(++ii)=mfexp(incident_sel(ir,ip,fi));
      otc(ii)=obs_tot_catch(ir,ip,fi);
    }
    else
    {
      int i=parent(ir,ip,fi);
      dvector dv=fml_designvec(ir,ip,fi);
      ivector num_ifmlrp=implicit_fml_bounds(3);
      dvar_vector ests=implicit_fm_level_regression_pars(i)(1,num_ifmlrp(i));
      // reorganize log effort by fishery so that this is simple
      MY_DOUBLE_TYPE leff=log_effort_by_fishery(i,fish_times(ir,ip,fi));
      dvariable tmp2 = dv * inv(fml_R(i)) * ests;
      dvariable log_fmlevel= tmp2+leff;
      fm_level(ir,ip,fi)=log_fmlevel;
      dvar_vector logfm=log_fmlevel+incident_sel(ir,ip,fi);
      cout << "logfm" << endl << logfm << endl << " exp(logfm)"
           << endl   << exp(logfm) << endl;
      other_mort+= exp(log_fmlevel+incident_sel(ir,ip,fi));
    }
   */   
