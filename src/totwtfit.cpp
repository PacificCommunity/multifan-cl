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
extern dvector mean_weights_kludge;

dvariable total_weight_fit(dvar_fish_stock_history& fsh,int print_switch,
  int avg_calc_flag)
{
  dvariable xy=0.;
  dvariable tmp1=0.;
  int ntimes;
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
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
        if (fsh.obs_tot_catch(ir,ip,fi)>0)
        {
          if (avg_calc_flag)
          {
            fsh.totalcatch_by_weight(ir)+=fsh.tot_catch(ir,ip,fi);
            fsh.obstotalcatch_by_weight(ir)+=fsh.obs_tot_catch(ir,ip,fi);
            fsh.numtotalcatch_by_weight(ir)+=1;
          }
          else
          {
            tmp1+= .01*fsh.age_flags(144)*square(log(1.+
              exp(fsh.catch(ir,ip,fi))*mean_weights_kludge)
		      -log(1.+fsh.obs_tot_catch(ir,ip,fi)));
          }
        }
      }
    }
  }
  xy+=tmp1;
  cout << " mean_weights_kludge " << endl << mean_weights_kludge << endl;
  if (print_switch)
  {
    cout << " after total weight fit = " << xy << endl;
    cout << "exp(q0)= " << setprecision(10) <<exp(fsh.q0) << "  "
       << setprecision(10) << exp(fsh.q0(1)) << endl;
    cout << "total catch contribution = "<< tmp1 << endl;
  }
  return xy;
}

#undef HOME_VERSION
