/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

void dvar_fish_stock_history::un_normalize_year(void)
{
  really_true_month.initialize();
  really_true_year.initialize();
  if (direction_flag==-1)
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
      {
        //fra[i].month+=(13-month_1);
        int tm=true_month(ir,ip);
        if (tm<=13-month_1)
        {
          really_true_month(ir,ip)=tm+month_1-1;
          really_true_year(ir,ip)=true_year(ir,ip);
        }
        else
        {
          really_true_month(ir,ip)=tm+month_1-13;
          really_true_year(ir,ip)=true_year(ir,ip)+1;
        }
      }
    }
  }	
  else if (direction_flag==1)
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
      {
        //fra[i].month+=(13-month_1);
        int tm=true_month(ir,ip);
        if (tm<=13-month_1)
        {
          really_true_month(ir,ip)=tm+month_1-1;
          really_true_year(ir,ip)=true_year(ir,ip)-1;
        }
        else
        {
          really_true_month(ir,ip)=tm+month_1-13;
          really_true_year(ir,ip)=true_year(ir,ip);
        }
      }
    }
  }
  else if (month_1==1)
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
      {
        really_true_month(ir,ip)=true_month(ir,ip);
        really_true_year(ir,ip)=true_year(ir,ip);
      }
    }
  }
}
