/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


#include "all.hpp"

d3_array dvar_len_fish_stock_history::report_catch_by_year(void)
{
  int i,j,k;
  int ny=max(really_true_year);
  d3_array catch_by_fishery(1,num_fisheries,1,ny,1,nage);
  catch_by_fishery.initialize();
  for (i=1;i<=num_regions;i++)
  {
    for (j=1;j<=num_fish_periods(i);j++)
    {
      for (k=1;k<=num_fish_incidents(i,j);k++)
      {
        catch_by_fishery(parent(i,j,k),really_true_year(i,j))
          +=exp(value(catch(i,j,k))); 
      }
    }
  }
  return catch_by_fishery;
}

d3_array dvar_len_fish_stock_history::report_estimated_catch_by_year(void)
{
  int i,j,k;
  int ny=max(really_true_year);
  d3_array estimated_catch_by_fishery(1,num_fisheries,1,ny,1,nage);
  estimated_catch_by_fishery.initialize();
  for (i=1;i<=num_regions;i++)
  {
    for (j=1;j<=num_fish_periods(i);j++)
    {
      for (k=1;k<=num_fish_incidents(i,j);k++)
      {
        if (age_sample_size(i,j,k)>0)
          estimated_catch_by_fishery(parent(i,j,k),really_true_year(i,j))
            +=age_freq(i,j,k); 
      }
    }
  }
  return estimated_catch_by_fishery;
}
