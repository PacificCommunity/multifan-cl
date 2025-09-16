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


void  print_FF(dvar_fish_stock_history& fsh);
void  calc_FF_mean(dvar_fish_stock_history& fsh);


void dvar_len_fish_stock_history::catch_equations_calc_implicit_mc(dvar_vector& sv,
  dvariable& ffpen)
{
  dvariable fpen1=0.0;
   missing_catch_counter.initialize();
  //greport("beginning catch_equations_calc");
  tmprecr.initialize();
  tmpinitpop.initialize();
  MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
  tot_mort.initialize();
  // calculate the totalfishing mortality and survival rates for each
  // fishing period in each region

  if (age_flags(69) && !age_flags(114))
  {
      setup_diffusion();
  }
  get_initial_population(sv,0,0);

  int current_recruitment_period=1;
  ivector rip(1,num_regions);
  rip=1;
  do
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      int& ip=rip(ir);
      if (recruitment_period(ir,ip)==current_recruitment_period)
      {
        num_fish(ir,ip)=N(ir,current_recruitment_period);
	do_newton_raphson(ir,ip,ffpen,fm_level(ir,ip),
          missing_catch_counter,fm_level_devs,fpen1);
        do
        {
          if (ip>=num_fish_periods(ir)) break;
          if (recruitment_period(ir,ip+1)==current_recruitment_period)
          {
            num_fish(ir,ip+1)=num_fish(ir,ip)-tot_mort(ir,ip);
            ip++;
	    do_newton_raphson(ir,ip,ffpen,fm_level(ir,ip),
              missing_catch_counter,fm_level_devs,fpen1);
          }
          else
          {
	    // does this assume a fishery in each recruitment period
            num_fish(ir,ip+1,1)=N(ir,current_recruitment_period+1,1);

            --num_fish(ir,ip+1)(2,nage)=
              num_fish(ir,ip)(1,nage-1)-tot_mort(ir,ip)(1,nage-1);

            num_fish(ir,ip+1,nage)=
              log(1.e-10 + mfexp(num_fish(ir,ip,nage-1)-tot_mort(ir,ip,nage-1))
                + mfexp(num_fish(ir,ip,nage)-tot_mort(ir,ip,nage)) );
  
            if (current_recruitment_period <num_recruitment_periods)
            {
              N(ir,current_recruitment_period+1)(2,nage)
	        =num_fish(ir,ip+1)(2,nage);
            }

            ip++;
            break;
          }
        }
        while (1);
      }
      else   // there were no fisheries in this region for this year
      {
        if (current_recruitment_period <num_recruitment_periods)
        {
          if (!pmsd)
          {
            for (int j=1;j<nage;j++)      // Loop over age classes
            {
              N(ir,current_recruitment_period+1,j+1)=
  	        N(ir,current_recruitment_period,j)-
                nat_mort(current_recruitment_period,j);
            }
          }
          else
          {
            for (int j=1;j<nage;j++)      // Loop over age classes
            {
              N(ir,current_recruitment_period+1,j+1)=
  	        N(ir,current_recruitment_period,j)-
                get_nat_mort_region(ir)(current_recruitment_period,j);
            }
          }
        }
      }
    }
    if (current_recruitment_period <num_recruitment_periods)
    {
      // Changed af(57) to af(53) J.H. 27 Nov 01
      if (age_flags(53))
      {
        if ( !((current_recruitment_period)%age_flags(53)) )
        {
          if (num_regions>1) do_the_diffusion(current_recruitment_period+1,
            sv,N);
        }
      }
      else
      {
        if (num_regions>1) do_the_diffusion(current_recruitment_period+1,sv,N);
      }
    }
    current_recruitment_period++;
  }
  while (current_recruitment_period <=num_recruitment_periods);
  //greport("B catch_equations_calc");

  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)  
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {        
        catch(ir,ip,fi)=fish_mort(ir,ip,fi)-log(1.e-10+tot_mort(ir,ip))+
          log(one_plus-survival(ir,ip))+num_fish(ir,ip);
      }
    }
  }
  //greport("leaving catch_equations_calc");
  print_FF(*this);
  calc_FF_mean_partial(*this);
}


