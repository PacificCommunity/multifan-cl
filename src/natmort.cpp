/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#define USE_NO_LAPACKE
#pragma implementation "newmult.hpp"
#pragma implementation "variable.hpp"
#include "version.h"
#include "all.hpp"



dvar4_array  dvar_fish_stock_history::
get_initial_equilibrium_survival_natmort
  (ivector & mps,ivector * pq_flag)
{
  // Now we need to determine the total morality in each
  // between movement period

  int af57=age_flags(57);
  dvar4_array tmort(1,af57,
    1,mps,1,num_regions,1,nage);  
//  dvar4_array surv(1,af57,
//    1,mps,1,nage,1,num_regions);  
  dvar4_array local_surv(1,af57,
    1,mps,1,nage,1,num_regions);  //NMD_24jan2022
  if (!allocated(csurv))
  {
    csurv.allocate(1,af57,1,mps,1,nage,1,num_regions);  
    csurv_chk.allocate(0,10);
    for (int i=0;i<=10;i++)
    {
      csurv_chk(i).allocate(1,af57,1,mps,1,nage,1,num_regions);  
    }
    csurv_chk.initialize();
  }

  tmort.initialize();  
  i3_array ny(1,num_regions,1,age_flags(57),1,imatrix(1,num_regions,mps));
  ny.initialize();
  local_surv.initialize(); //NMD_9May2022
  
  int ir;
  int newyear_flag=0;
  for (ir=1;ir<=num_regions;ir++)
  {
    int oldyr=year(ir,1);
    int mpp=1;
    int is=(oldyr-1)%age_flags(57)+1;
    {
//      if (age_flags(128)==0)  
      if (age_flags(128)==0 || pq_flag)     //NMD_6Jun2022
      {
        tmort(is,mpp,ir)+=mfexp(nat_mort(year(ir,1))+fraction(ir,1));
      }
      else
      {
        tmort(is,mpp,ir)+=age_flags(128)/100.*           //NMD_4jun2025
          mfexp(get_nat_mort_region(ir)(year(ir,1))+fraction(ir,1));
      }
    }

    for (int ip=2;ip<=num_real_fish_periods(ir);ip++)  
    {
      int yr=year(ir,ip);
      if (yr>af57) break;   //NMD_25jan2022 - break after calendar year recrs  
      int is=(yr-1)%age_flags(57)+1;
      if (yr>1)
      {
        if (yr>oldyr)  // new year and we ony do this for first year
        {
//          break;   //NMD_25jan2022 - only works for af57=1
          newyear_flag=1;   // same conditions as for totmort method
          oldyr=yr;
          mpp=1;
        }
        else if (move_flags(ir,ip))
        {
          mpp++;
        }
      }

      {
//        if (age_flags(128)==0)
        if (age_flags(128)==0 || pq_flag)     //NMD_3Jun2022
          tmort(is,mpp,ir)+=mfexp(nat_mort(year(ir,ip))
            +fraction(ir,ip));
        else
          tmort(is,mpp,ir)+=
            age_flags(128)/100.*mfexp(nat_mort(year(ir,ip)) //NMD_4Jun2025
            +fraction(ir,ip));
      }
    }
  }
  for (int is=1;is<=age_flags(57);is++)
  {
    for (int pi=1;pi<=mps(is);pi++)
    {
      if (pq_flag)
      {
//        surv(is,pi)=mfexp(-trans(tmort(is,pi))-0.0001);
        local_surv(is,pi)=mfexp(-trans(tmort(is,pi))-0.0001);
      }
      else
      {
//        surv(is,pi)=mfexp(-trans(tmort(is,pi))-0.0001);
        local_surv(is,pi)=mfexp(-trans(tmort(is,pi))-0.0001);
      }
//      csurv(is,pi)=value(surv(is,pi));
//      csurv_chk(iloop,is,pi)=value(surv(is,pi));
      csurv(is,pi)=value(local_surv(is,pi));
      csurv_chk(iloop,is,pi)=value(local_surv(is,pi));
    }
  }
//  return surv;
  return local_surv;
}
