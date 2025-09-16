/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


#include "all.hpp"
#include "scbd.hpp"


int dvar_fish_stock_history::new_fishing_selectivity_interface_nvar(void)
{
  int nv=0;
  ivector ff19=column(fish_flags,19);
  ivector ff51=column(fish_flags,51);
  if (sum(ff51))
  {
    cerr << "grouping not implemented yet" << endl;
    ad_exit(1);
  }
  for (int fi=1;fi<=num_fisheries;fi++)
  {
    if (ff19(fi)>0)
    {
      for (int it=1;it<=num_fish_times(fi);it++)
      {
        if (ff19(fi)>sel_dev_coffs(fi,it).indexmax())    
        {
          cerr << "ff19 to large for fishery " << fi 
               << " incident " << it << endl;
          ad_exit(1);
        }
        int rr=realization_region(fi,it);
        int rp=realization_period(fi,it);
        int ri=realization_incident(fi,it);

	if(len_sample_size(rr,rp,ri)>0 || wght_sample_size(rr,rp,ri)>0)
        {
          nv+=ff19(fi);
        }
      }
    }         
  }
  return nv;
}

void dvar_fish_stock_history::new_fishing_selectivity_interface_reset
  (int& ii,dvar_vector& x,dvariable& fpen,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s)
{
  ivector ff19=column(fish_flags,19);
  for (int fi=1;fi<=num_fisheries;fi++)
  {
    if (ff19(fi)>0)
    {
      for (int it=1;it<=num_fish_times(fi);it++)
      {
        int rr=realization_region(fi,it);
        int rp=realization_period(fi,it);
        int ri=realization_incident(fi,it);

	if(len_sample_size(rr,rp,ri)>0 || wght_sample_size(rr,rp,ri)>0)
        {
          set_value_partial(sel_dev_coffs(fi,it),x,ii,ff19(fi),
            fmin,fmax,fpen,s);
        }
      }
    }         
  }
}
void dvar_fish_stock_history::new_fishing_selectivity_interface_xinit
  (int& ii,dvector& x,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s)
{
  ivector ff19=column(fish_flags,19);
  for (int fi=1;fi<=num_fisheries;fi++)
  {
    if (ff19(fi)>0)
    {
      for (int it=1;it<=num_fish_times(fi);it++)
      {
        int rr=realization_region(fi,it);
        int rp=realization_period(fi,it);
        int ri=realization_incident(fi,it);

	if(len_sample_size(rr,rp,ri)>0 || wght_sample_size(rr,rp,ri)>0)
        {
          set_value_inv_partial(value(sel_dev_coffs(fi,it)),x,ii,ff19(fi),
            fmin,fmax,s);
        }
      }
    }         
  }
}
