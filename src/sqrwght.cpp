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


// put int the factor of 2 and use proper distribution

dvariable square_fit_wght(dvar_len_fish_stock_history& fsh,
  d3_array& total_freq)
{
  dvariable wght_tmp=0.;
  dvariable tmp_sigma=0.;
  int nwint =fsh.nwint;
  // dave f put this back to 0.1/nwint oct 20 2014
  //double wt=0.001/nwint;
  MY_DOUBLE_TYPE wt=0.1/nwint;
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    int ntimes;
    if (fsh.parest_flags(142)==0)
    {
      ntimes=fsh.num_fish_periods(ir);
    }
    else 
    {
      ntimes=min(fsh.parest_flags(142),fsh.num_fish_periods(ir));
    }

    for (int ip=1;ip<=ntimes;ip++)
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
      {
        if (fsh.wght_sample_size(ir,ip,fi)>0)
        {
          dvar_vector r2=0.5*elem_div(
            square(fsh.wght_freq(ir,ip,fi)-fsh.wtprob(ir,ip,fi)),
              wt+fsh.wtprob(ir,ip,fi));

          wght_tmp-=total_freq(ir,ip,fi)*
            sum(log(exp(-r2) + 0.03/(1.0+r2)));

          if (fsh.parest_flags(161)==0)
          {
            wght_tmp+=0.5*sum(log(wt+fsh.wtprob(ir,ip,fi)));
          }
        }
      }
    }
  }
  return wght_tmp;
}

// LH function as in YFT paper

#undef HOME_VERSION
