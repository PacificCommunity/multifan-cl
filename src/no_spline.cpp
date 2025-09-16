/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
  
void block_error();
extern dvar_vector * const_vb;
extern dvar_vector * const_var;

void dvar_len_fish_stock_history::no_splines(void)
{
  for (int i=1;i<=num_fisheries;i++)
  {
    if (fish_flags(i,57) !=3 )
    {
      int rr=realization_region(i,1);
      int rp=realization_period(i,1);
      int ri=realization_incident(i,1);
      incident_sel(rr,rp,ri) = 1.;
    }
  }
}
