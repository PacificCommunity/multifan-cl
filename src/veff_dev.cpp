/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

void dvar_fish_stock_history::effort_devs_calc()
{
  effort_devs.initialize();
  if  (projection_sim_index==0 || parest_flags(234)==0)
  {
    for (int i=1;i<=num_fisheries;i++)
    {
      if (fish_flags(i,4)>0)
      {
        for (int j=1;j<=num_fish_times(i);j++)
        {
          effort_devs(realization_region(i,j),
  		  realization_period(i,j),realization_incident(i,j))
  					      =effort_dev_coffs(i,j);
        }
      }
    }
  }
  else
  { 
    int num_sims=age_flags(20);
    int pf234=parest_flags(234);
    int pf235=parest_flags(235);
    if (pf234>0)
    {
      if (!allocated(simulated_effort_devs))
      {
        ivector j1=num_real_fish_times+1;
        imatrix i1=imatrix(1,num_sims,j1);
        imatrix i2=imatrix(1,num_sims,num_fish_times);
        simulated_effort_devs.allocate(1,num_sims,1,num_fisheries,
          i1,i2);
        int iseed=1673;
        if (parest_flags(236)>0)
         iseed=parest_flags(236);
        random_number_generator rng(iseed);
        simulated_effort_devs.fill_randn(rng);
        MY_DOUBLE_TYPE sd=0.1;
        if (pf235>0) sd=pf235/100.;
        int i;
        for (i=1;i<=num_sims;i++)
          simulated_effort_devs(i)*=sd;
      }
    }
    for (int i=1;i<=num_fisheries;i++)
    {
      if (fish_flags(i,4)>0)
      {
        for (int j=1;j<=num_real_fish_times(i);j++)
        {
          effort_devs(realization_region(i,j),
  		  realization_period(i,j),realization_incident(i,j))
  					      =effort_dev_coffs(i,j);
        }
        for (int j=num_real_fish_times(i)+1;j<=num_fish_times(i);j++)
        {
          effort_devs(realization_region(i,j),
            realization_period(i,j),realization_incident(i,j))
               = simulated_effort_devs(projection_sim_index,i,j);
        }
      }
    }
  }
}
