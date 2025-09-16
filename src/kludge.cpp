/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"


void dvar_fish_stock_history::kludge_the_date(int & year,int & month)
{
  year-=year1-1;
  if (month_1 !=1)
  {
    if (direction_flag==-1)
    {
      if (month>=month_1)
      {
        month-=(month_1-1);
      }
      else
      {
        year-=1;
        month+=12;
        month-=(month_1-1);
      }
    }
    else if (direction_flag==1)
    {
      month+=(13-month_1);
      if (month>12)
      {
        year+=1;
        month-=12;
      }
    }
    if (month_factor!=0 && month_factor!=1)
    {         
      date_struc newdate=get_new_time(year,
        month,1,month_factor,first_time);
      year=newdate.year;
      month=newdate.month;
    }
  }
}
