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
#include "variable.hpp"
#include "newmprot.hpp"

extern dvar_vector * psv;

dvariable no_penalties(dvar_len_fish_stock_history& fsh,
  dvariable& mean_length_constraints,int print_switch,ofstream* pof_pen)
{
   dvariable xy=0.0;
   dvariable ppf_tmp=0.0;   //NMD_21nov2023
   ppf_tmp=xy;
   // don't do this if we are using orthogonal polynomials
   if (fsh.parest_flags(155)==0)
   {
     if (fsh.age_flags(70))
     {
       if (!fsh.age_flags(178))
       {
         xy+=norm2(fsh.region_rec_diff_sums);
       }
       else
       {
         xy+=100.*norm2(fsh.rec_delta);
       }
       if (fsh.pmsd)
       {
         int ns=fsh.pmsd->num_species;
         for (int is=2;is<=ns;is++)
         {
           if (!fsh.age_flags(178))
           {
             xy+=norm2(fsh.pmsd->region_rec_diff_sums(is));
           }
           else
           {
             xy+=norm2(fsh.pmsd->rec_delta(is));
           }
         }
       }
     }
   }
   ppf_tmp=xy-ppf_tmp;  //NMD_21nov2023
   if (fsh.ppstf && !fsh.af170q0 && print_switch)
   {
     fsh.ppstf->reg_recr_pen=ppf_tmp;
     ppf_tmp=0.0;
   }
   
   return xy;
}
