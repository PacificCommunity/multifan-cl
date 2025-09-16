/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


#define HOME_VERSION
#include "all.hpp"
//#include "here.h"
//#include "variable.hpp"
#include "newmprot.hpp"

extern dvar_vector * psv;

dvariable good_penalties(dvar_len_fish_stock_history& fsh,
  dvariable& mean_length_constraints,int print_switch,ofstream* pof_pen)
{
  dvariable xy=0.;
  dvariable ppf_tmp=0.0;   //NMD_21nov2023

    // biomass depletion target 
    //  if (fsh.age_flags(170)>0)
    //  {
    //    if (fsh.age_flags(175)>0)
    //    {
    //      xy+=fsh.calculate_the_totalbiomass_depletion();
    //    }
    //  }
    // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
    // penalties to reduce the size of diffusion

    if (fsh.age_flags(184)==2)
    {
      xy+=fsh.movement_coffs_option_2_penalty();
    }
    else
    {
      if (fsh.age_flags(68))
      {
        if (fsh.age_flags(27)==0)
        {
          xy+=5.*norm2(fsh.diff_coffs);
        }
        else if (fsh.age_flags(27)>0)
        {
          xy+=fsh.age_flags(27)/10.0*norm2(fsh.diff_coffs);
        }
        else if (fsh.age_flags(27)<0)
        {
          xy-=fsh.age_flags(27)/10.0*norm2(fsh.diff_coffs-fsh.diff_coffs_prior);
        }
      }
    }
    ppf_tmp=value(xy);
    if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_21nov2023
    {
      fsh.ppstf->diff_coff_pen=ppf_tmp;
      ppf_tmp=0.0;
    }
    
    return xy;
}
