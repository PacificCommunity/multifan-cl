/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"
#ifdef __MSC__
  MY_DOUBLE_TYPE mfexp(const MY_DOUBLE_TYPE& );
  dvariable mfexp(const prevariable& );
  dvar_vector mfexp(const dvar_vector& );
  dvector mfexp(const dvector& );
#else
  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
  dvariable mfexp(prevariable& );
  //dvar_vector mfexp(dvar_vector& );
  //dvector mfexp(dvector& );
#endif

  extern dvar_vector * psv;
  extern dvar_len_fish_stock_history * pcfsh;


void dvar_fish_stock_history::xcatch_equations_calc(dvar_vector& sv)
{
  pcfsh->set_global_variance();
  set_global_vars_flag=1;
  pcfsh->mean_lengths_by_year_calc();
  pcfsh->mean_weights_by_year_calc();
  //greport("beginning catch_equations_calc");
  tmprecr.initialize();
  tmpinitpop.initialize();
  MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
  // calculate the totalfishing mortality and survival rates for each
  // fishing period in each region
  //cout <<  num_fish_periods << endl;
  int ir;

  //if (age_flags(69))
  if (num_regions>1 && !age_flags(114))
  {
    setup_diffusion();
  }

  get_initial_population(sv,0,0);
  tot_mort.initialize();
  fish_mort.initialize();

  int current_year=1;
  ivector rip(1,num_regions);
  rip=1;
  int finished_flag=1;
  int break_flag=0;
  ivector tmp_mp(1,num_regions);
  ivector tmp_ip(1,num_regions);
  for (ir=1;ir<=num_regions;ir++)
  {
    num_fish(ir,1)=N(ir,1);
  }

  dvar_vector lbio1=log(1.e-20+pcfsh->calculate_the_current_biomass(rip));
  dvar_vector lcurrent_rbio(1,num_regions);
  lcurrent_rbio.initialize();
   
  do
  {
    finished_flag=1;
    for (int ir=1;ir<=num_regions;ir++)
    {
      break_flag=0;
      int& ip=rip(ir);
      do
      {
        // need to calculate fish_mort(ir,ip) and tot_mort(ir,ip)
        // here
        get_fishing_and_total_mortality_for_this_period(ir,ip,lcurrent_rbio);
        do_fish_mort_intermediate_calcs(ir,ip);

        if (ip>=num_fish_periods(ir)) break;
        finished_flag=0;

        if (year(ir,ip+1)==year(ir,ip))
        {
          num_fish(ir,ip+1)=num_fish(ir,ip)-tot_mort(ir,ip);
        }
        else
        {
          // age the fish
          num_fish(ir,ip+1,1)=N(ir,year(ir,ip+1),1);

          if (nage>2)
          --num_fish(ir,ip+1)(2,nage-1)=
            num_fish(ir,ip)(1,nage-2)-tot_mort(ir,ip)(1,nage-2);

          num_fish(ir,ip+1,nage)=
            log(1.e-10 + mfexp(num_fish(ir,ip,nage-1)-tot_mort(ir,ip,nage-1))
              + mfexp(num_fish(ir,ip,nage)-tot_mort(ir,ip,nage)) );
        }
        if (move_flags(ir,ip))
        {
          tmp_ip(ir)=ip;
          tmp_mp(ir)=move_index(ir,ip);
          break_flag=1;
        } 
        ip++;
      }
      while (!break_flag); 
    }
    // move the fish
    {
      if (num_regions>1)
      {
        check_sanity(tmp_mp);
        dvar_matrix tmp=fast_diffusion_calcs(nage,num_regions,num_fish,
          Dad(tmp_mp(1)),rip);
        for (int ir=1;ir<=num_regions;ir++)
        {
          num_fish(ir,rip(ir))=log(5.e-10+tmp(ir));
        }
      }
      dvar_vector current_biomass=pcfsh->calculate_the_current_biomass(rip);
      lcurrent_rbio=log(current_biomass)-lbio1;
    }
  } // need to decide when to quit
  while (!finished_flag);
  
  for (ir=1;ir<=num_regions;ir++)
  {
    N(ir,year(ir,1))=num_fish(ir,1);
    for (int ip=2;ip<=num_fish_periods(ir);ip++)  
    {
      if (year(ir,ip)>year(ir,ip-1))
        N(ir,year(ir,ip))=num_fish(ir,ip);
    }
  }
  
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)  
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {        
        catch(ir,ip,fi)=fish_mort_calcs(ir,ip,fi)+num_fish(ir,ip);
      }
    }
  }
  
 
 /*
  cout << "num_fish" << endl;
  cout << num_fish << endl;
  cout << "catch" << endl;
  cout << catch << endl;
  exit(1);
 */
 
}

dvar_vector dvar_len_fish_stock_history::calculate_the_current_biomass
  (ivector& rip)
{
  int ir;
  dvar_vector vvar=global_vars;
  //double lwc=len_wt_coff;
  dvar_vector current_biomass_by_region(1,num_regions);
  for (ir=1;ir<=num_regions;ir++)
  {
    if (!age_flags(112))
    {
      current_biomass_by_region(ir)
         =exp(num_fish(ir,rip(ir)))*mean_weight_yr(ir,1);
    }
    else
    {
      current_biomass_by_region(ir)=sum(exp(num_fish(ir,rip(ir))));
    }
  }
  return current_biomass_by_region;
}

void dvar_fish_stock_history::
  get_fishing_and_total_mortality_for_this_period(int ir,int ip,
  dvar_vector& lcurrent_rbio)
{
  dvar_vector tmp_bio=fish_pars(8);
    
  dvar_vector& tm=tot_mort(ir,ip);
  tm.initialize();
  for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
  {                                         // incidents for this period
    int i=parent(ir,ip,fi);
    int rr=realization_region(i,1);
    int rp=realization_period(i,1);
    int ri=realization_incident(i,1);
    dvar_vector& fm=fish_mort(ir,ip,fi);
    dvar_vector& fs=incident_sel(rr,rp,ri);
    MY_DOUBLE_TYPE eff=effort(ir,ip,fi);
    const prevariable& cat=catchability(ir,ip,fi);
    fm=fs+tmp_bio(i)*lcurrent_rbio(ir)+eff+cat;
    if (age_flags(34)>0)
    {
      fm+=effort_devs(ir,ip,fi);
    }
    tm+=mfexp(fm);
  }
  tm+=mfexp(get_nat_mort_region(ir)(year(ir,ip))+fraction(ir,ip));
  survival(ir,ip)=mfexp(-tm);
}
void dvar_fish_stock_history::
  get_fishing_and_total_mortality_for_this_period_q0(int ir,int ip,
  dvar_vector& lcurrent_rbio)
{
  dvar_vector tmp_bio=fish_pars(8);
  dvar_vector& tm=tot_mort_q0(ir,ip);
  tm.initialize();
  for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
  {                                         // incidents for this period
    int i=parent(ir,ip,fi);
    int rr=realization_region(i,1);
    int rp=realization_period(i,1);
    int ri=realization_incident(i,1);
    dvar_vector& fm=fish_mort_q0(ir,ip,fi);
    dvar_vector& fs=incident_sel(rr,rp,ri);
    MY_DOUBLE_TYPE eff=effort(ir,ip,fi);
    const prevariable& cat=catchability_q0(ir,ip,fi);
    fm=fs+tmp_bio(i)*lcurrent_rbio(ir)+eff+cat;
    if (age_flags(34)>0)
    {
      fm+=effort_devs(ir,ip,fi);
    }
    tm+=mfexp(fm);
  }
  tm+=mfexp(get_nat_mort_region(ir)(year(ir,ip))+fraction(ir,ip));
  survival_q0(ir,ip)=mfexp(-tm);
}
