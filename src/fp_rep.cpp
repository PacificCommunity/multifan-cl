/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
int dvar_fish_stock_history::have_fish_periods(ivector& rip)
{
  int flag=0;
  int ir;
  for (ir=1;ir<=num_regions;ir++)
  {
    if (rip(ir)<=num_fish_periods(ir))
    {
      flag=1;
      break;
    }
  }
  return flag;
}

ivector dvar_fish_stock_history::get_min_time_indices(ivector& rip)
{
  ivector dates(1,num_regions);
  int ir;
  for (ir=1;ir<=num_regions;ir++)
  {
    if (rip(ir)<=num_fish_periods(ir))
    {
      dates(ir)= 48*(year(ir,rip(ir))-1)+4*(month(ir,rip(ir))-1)+
        week(ir,rip(ir));
    }
    else
    {
      dates(ir)=1000000;
    }
  }
  int min_date=min(dates);
  int ii=0;
  for (ir=1;ir<=num_regions;ir++)
  {
    if (dates(ir)==min_date) ii++;
  }
  if (ii==0)
  {
    cerr << "Error" << endl;
    ad_exit(1);
  }
  ivector tmp(1,ii);
  ii=0;
  for (ir=1;ir<=num_regions;ir++)
  {
    if (dates(ir)==min_date) 
    {
      ii++;
      tmp(ii)=ir;
    }
  }
  return tmp;
}
  
void dvar_fish_stock_history::make_fishing_period_report(void)
{
  ofstream ofs("fishing_period.rep");
  int ir; 
  ivector rip(1,num_regions);
  rip=1;
  do
  {
    ivector indices=get_min_time_indices(rip);
  
    int in=indices(1);
  
  
    ofs  << setw(4) << year(in,rip(in)) << " " 
         << setw(2) << month(in,rip(in)) << " " 
         << setw(2) << week(in,rip(in)) << " ";
  
    int ii=1;
    for (ir=1;ir<=num_regions;ir++)
    {
      if (ii <= indices.indexmax() && ir==indices(ii))
      {
        if (move_index(ir,rip(ir))>0)
        {
          ofs << "  M "
              << setw(3) << rip(ir) << ":" 
              << setw(2)<< num_fish_incidents(ir,rip(ir)) << "   ";
        }
        else
        {
          ofs << "  x "
              << setw(3) << rip(ir) << ":" 
              << setw(2)<< num_fish_incidents(ir,rip(ir)) << "   ";
        }
        rip(ir)++;
        ii++;
      }
      else
      {
        ofs << "             ";
      }
    }  
    ofs << endl;
  }
  while(have_fish_periods(rip));
  
}     

