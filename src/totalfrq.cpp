/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE



#define HOME_VERSION
#include "all.hpp"

#include "variable.hpp"
#include "newmprot.hpp"

d3_array total_freq_calc(dvar_len_fish_stock_history& fsh,int print_switch)
{
  d3_array total_freq(1,fsh.num_regions,1,fsh.num_fish_periods,1,
    fsh.num_fish_incidents);
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)  
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        int pp1=fsh.parent(ir,ip,fi);
        total_freq(ir,ip,fi)=
          mmiinn(fsh.len_sample_size(ir,ip,fi)/fsh.effective_len_size(pp1),
          1000./fsh.effective_len_size(pp1));
      }
    }
  }
  return total_freq;
}

d3_array total_wght_calc(dvar_len_fish_stock_history& fsh,int print_switch)
{
  ofstream * pppofs = 0;
  /*
  if (print_switch)
  {
    pppofs =new ofstream("wghtscale");
  }
  */
  d3_array total_wght(1,fsh.num_regions,1,fsh.num_fish_periods,1,
    fsh.num_fish_incidents);
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)  
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        int pp1=fsh.parent(ir,ip,fi);
        total_wght(ir,ip,fi)=
          mmiinn(fsh.wght_sample_size(ir,ip,fi)/fsh.effective_weight_size(pp1),
           1000./fsh.effective_weight_size(pp1));
        if (pppofs)
        {
          *pppofs << ir << " " << ip << " " << fi << " " << pp1 << " " 
              << total_wght(ir,ip,fi) << " "  
              << fsh.wght_sample_size(ir,ip,fi) << endl;
        }
      }
    }
  }
  if (pppofs)
  {
    delete pppofs;
    pppofs=0;
  }
  return total_wght;
}

dvar4_array weights_calc(dvar_len_fish_stock_history& fsh,
  d3_array& total_freq,int print_switch)
{
  dvar4_array weights;
  weights.allocate(fsh.len_freq);
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++) 
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        if (fsh.len_sample_size(ir,ip,fi)>0)
        {
          weights(ir,ip,fi)=(.0001+
            elem_prod(fsh.tprob(ir,ip,fi),1.-fsh.tprob(ir,ip,fi)))/
	    total_freq(ir,ip,fi);
        }
        else
        {
          weights(ir,ip,fi).initialize();
        }
      }
    }
  }
  return weights;
}
#undef HOME_VERSION
