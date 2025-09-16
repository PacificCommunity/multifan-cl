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

dvariable catch_at_length_fit(dvar_len_fish_stock_history& fsh,
  int print_switch)
{
  dvariable tmp=0.;
  if (fsh.age_flags(38)>0)
  {
    d3_array total_freq=total_freq_calc(fsh,print_switch);
    dvar4_array weights=weights_calc(fsh,total_freq,print_switch);
    tmp+=robust_freq_partial_fit(fsh.len_freq, fsh.tprob,weights,fsh);
  }
  else
  {
    d3_array total_freq=total_freq_calc(fsh,print_switch);
    if (fsh.parest_flags(140)==0)
    {
      dvariable len_tmp=0.0;
      switch (fsh.parest_flags(141))
      {
      case 0:
        len_tmp=square_fit(fsh,total_freq);
        break;
      case 1:
        len_tmp=square_fit0(fsh,total_freq);
        break;
      case 2:
        len_tmp=square_fit0a(fsh,total_freq);
        break;
      case 3:
        {
          if (fsh.parest_flags(311)==0)
          {
            len_tmp=fsh.square_fita(total_freq,fsh.len_sample_size,
              fsh.len_freq,fsh.tprob,
              fsh.ppstf->lencontrib,
              fsh.ppstf->lencontrib_by_realization,1,
              print_switch);

          }
          else
          {
            len_tmp=fsh.square_fita(total_freq,fsh.len_sample_size,
              fsh.tc_len_freq,fsh.tc_tprob,
              fsh.ppstf->lencontrib,
              fsh.ppstf->lencontrib_by_realization,1,
              print_switch);
          }
        }
        break;
      case 4:
        len_tmp=dirichlet_multinomial_fit(fsh,total_freq);
        break;
      case 5:
        len_tmp=dirichlet_multinomial_mixture_fit(fsh,total_freq);
        break;
      case 6:
        len_tmp=multinomial_fit(fsh,total_freq);
        break;
      case 7:
        // this is logistic normal or logistic_t
        if (fsh.parest_flags(299)==0)
          len_tmp=logistic_normal_fit(fsh,total_freq);
        else
        {
          if (fsh.parest_flags(311)==0)
          {
            len_tmp=logistic_normal_fit_heteroscedastic(fsh,total_freq);
          }
          else
          {
            len_tmp=tail_compressed_logistic_normal_fit_heteroscedastic
              (fsh,total_freq);
          }
        }
        break;
      case 8:
         len_tmp=square_fit_t(fsh,total_freq); // this it t dist
         break;
      case 9:
         // IIIII
         //len_tmp=xlen_self_scaling_multinomial_re_multi_rho_multi_var
          // (fsh,total_freq,print_switch);

         len_tmp=len_self_scaling_multinomial_re_multi_rho_multi_var
           (fsh,total_freq,print_switch);
         break;
      case 10:
         {
           int newswitch=1;
//           if (newswitch==0)   //NMD_9Dec2019
           if (fsh.parest_flags(138)==1)
             len_tmp= old_len_self_scaling_multinomial_nore
               (fsh,total_freq,print_switch);
           else
             len_tmp= new_len_self_scaling_multinomial_nore
               (fsh,total_freq,print_switch);
         }
         break;
      case 11:
         if (fsh.parest_flags(320))
         {
           len_tmp=len_dm_nore(fsh,total_freq,print_switch);
         }
         else
         {
           len_tmp=len_dm_nore_notc(fsh,total_freq,print_switch);
         }
         break;
      case 12:
         // IIIII
         len_tmp=new_len_self_scaling_multinomial_re_multi_rho_multi_var
           (fsh,total_freq,print_switch);
         break;
      default:
        cerr << "Illegal value for parest_flags(141)" << endl; 
        ad_exit(1);
      }
   
      tmp+=len_tmp;
      cout << "Length frequency data " << setprecision(9) << setw(19) << len_tmp << endl;//NMD
    }
  /*
    else
    {
      dvariable len_tmp=mean_fit(fsh,total_freq);
      tmp+=len_tmp;
    }
  */
  }
  return tmp;
}

#undef HOME_VERSION
