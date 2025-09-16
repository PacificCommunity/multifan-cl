/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

void dvar_len_fish_stock_history::set_some_flags(ivector& dataswitch,
  int min_year)
{
  if (!(!dataswitch))
  { 
    if (dataswitch.indexmax()>4)
    {
      if (age_flags(57)==0) 
      {
        if (dataswitch(8)) 
        {
          age_flags(57)=dataswitch(8);
        }
        else
        {
          age_flags(57)=1;
        }
      }
      if (dataswitch(5)) 
      {
        if (dataswitch(5)!=min_year)
        {
          cerr << " Now dataswitch(5) must equal " << min_year << endl;
          cerr << " you have " << dataswitch(5) << endl;
          ad_exit(1);
        }
        year1=min_year;
      }
      else
      {
        year1=min_year;
      }
    }
    else
    {
      year1=min_year;
    }
  }
  if (sum(data_fish_flags(2)))
  {
    have_projection_periods_flag=1;
    for (int fi=1;fi<=num_fisheries;fi++)
    {
      int year=data_fish_flags(2,fi);
      int month=data_fish_flags(3,fi);
      kludge_the_date(year,month);
      projection_year(fi)=year;
      projection_month(fi)=month;
    }
  }
  else
  {
    have_projection_periods_flag=0;
  }
}
