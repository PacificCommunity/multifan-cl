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

dvector cbiocalc(dvar_len_fish_stock_history& fsh);
dvector cbiocalc_wt(dvar_len_fish_stock_history& fsh);
extern dvector mean_weights_kludge;


void simreport(dvar_len_fish_stock_history& fsh)
{
  {
    adstring fn=adstring("nmort") + str(fsh.age_flags(20)) 
       + adstring(".rep");
    ofstream of(fn);
    of << exp(value(fsh.nat_mort(1))) << endl;
  }

/*
  {
    adstring fn=adstring("biom") + str(fsh.age_flags(20)) 
      + adstring(".rep");
    ofstream of(fn);

    if (fsh.parest_flags(41)==0)
    {
      dvector rbio=cbiocalc(fsh);
      rbio=rbio/max(rbio);
      of << rbio << endl;
    }
    else
    {
      dvector rbio=cbiocalc_wt(fsh);
      rbio=rbio/max(rbio);
      of << rbio << endl;
    }
  }
*/

  {
    adstring fn=adstring("recr") + str(fsh.age_flags(20)) 
       + adstring(".rep");
    ofstream of(fn);
    {
      dvector relrecr=value(column(fsh.N(1),1));
      of << setprecision(4)  << relrecr << endl << endl;
    }
  }
  {
    adstring fn=adstring("growth") + str(fsh.age_flags(20)) 
       + adstring(".rep");
    ofstream of(fn);
    of << setprecision(4) << fsh.vb_coff << endl;
  }

  {
    adstring fn=adstring("catch") + str(fsh.age_flags(20)) 
      + adstring(".rep");
    ofstream of(fn);
    for (int i=1;i<=fsh.num_fisheries;i++)
    {
      for (int nt =1;nt<=fsh.num_fish_times(i);nt++)
      {
        int rr=fsh.realization_region(i,nt);  //NMD_8May2018
        int rp=fsh.realization_period(i,nt);
        int ri=fsh.realization_incident(i,nt);
        of << setprecision(4) << "  " 
           << exp(value(fsh.catchability(rr,rp,ri)));  //NMD_8May2018
      }
      of << endl << endl; 
    }
  }
}
