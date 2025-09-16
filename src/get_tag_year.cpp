/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

void normalize_tag_dates(
  ivector & tag_year,
  ivector & tag_month,
  ivector & true_tag_year,
  int month_1,
  int year1,
  int direction_flag,
  int month_factor,
  int first_time,
  int it)
{
  date_struc new_date;  
  true_tag_year(it)=tag_year(it);
  tag_year(it)-=year1-1;
  if (month_1 !=1)
  {
    if (direction_flag==-1)
    {
      if (tag_month(it) >= month_1)
      {
        tag_month(it)-=(month_1-1);
      }
      else
      {
        tag_year(it)-=1;
        tag_month(it)+=12;
        tag_month(it)-=(month_1-1);
      }
    }  
    else if (direction_flag==1)
    {
      tag_month(it)+=(13-month_1);
      if (tag_month(it)>12)
      {
        tag_year(it)+=1;
        tag_month(it)-=12;
      }
    }
  }
  if (month_factor!=0 && month_factor!=1)
  {         
    date_struc newdate=get_new_time(tag_year(it),
      tag_month(it),1,month_factor,first_time);
    tag_year(it)=newdate.year;
    tag_month(it)=newdate.month;
  }
}
