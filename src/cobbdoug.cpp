/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#define HOME_VERSION
#include "all.hpp"

#ifndef __GNUC__
//#include <mf_menu.h>
#endif

void dvar_fish_stock_history::
  effort_multiplier_for_cobb_douglas_production_function(void)
{
  //cout << "A" << endl;
  int i;
  int num_global_fishing_periods=global_fishing_periods(1,num_fish_periods(1));
  //cout << "B" << endl;
  for (int ir=2;ir<=num_regions;ir++)
  {
    if (num_global_fishing_periods <
      global_fishing_periods(ir,num_fish_periods(ir)))
      num_global_fishing_periods=
        global_fishing_periods(ir,num_fish_periods(ir));
  }
  //cout << "C" << endl;
  ivector ff52=column(fish_flags,52);
  if (allocated(grouped_effort)) grouped_effort.deallocate();
  grouped_effort.allocate(min(ff52),max(ff52),
    1,num_global_fishing_periods);
  grouped_effort.initialize();

  //cout << "D" << endl;
  for (i=1;i<=num_fisheries;i++)
  {
    for (int nt=1;nt<=num_fish_times(i);nt++)
    {
      int rp=realization_period(i,nt);
      int ri=realization_incident(i,nt);
      int rr=realization_region(i,nt);
      int gfp=global_fishing_periods(rr,rp);
      grouped_effort(ff52(i),gfp)+=effort(rr,rp,ri);
    }
  } 

  //cout << "E" << endl;
  int imin=grouped_effort.indexmin();
  int imax=grouped_effort.indexmax();
  for (i=imin;i<=imax;i++)
  {
    grouped_effort(i)/=mean(grouped_effort(i))+1.e-20;
  }
  //cout << "F" << endl;
}
