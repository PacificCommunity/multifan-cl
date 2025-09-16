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

dvariable catch_at_weight_fit(dvar_len_fish_stock_history& fsh,
  int print_switch)
{
  d3_array total_wght=total_wght_calc(fsh,print_switch);
  dvariable tmp=0.0;
  if (fsh.parest_flags(140)==0)
  {
    dvariable wght_tmp=0.0;
    int i=fsh.parest_flags(139);
    switch (fsh.parest_flags(139))
    {
    case 0:
      wght_tmp=square_fit_wght(fsh,total_wght);
      break;
    case 1:
      cout << "Error need to finish this for option " << i << endl;
      exit(1);
      //wght_tmp=square_fit0_wght(fsh,total_wght);
      break;
    case 2:
      cout << "Error need to finish this for option " << i << endl;
      exit(1);
      //wght_tmp=square_fit0a_wght(fsh,total_wght);
      break;
    case 3:
      {
        /*
         dvariable tmpx=square_fita_wght(fsh,total_wght);
         dvariable tmpy=fsh.square_fita(total_wght,fsh.wght_sample_size,
           fsh.wght_freq,fsh.wtprob);
         cout << tmpx << " " << tmpy << endl;
         ad_exit(1);
        */
        if (fsh.parest_flags(301)==0)
        {
          wght_tmp=fsh.square_fita(total_wght,fsh.wght_sample_size,
            fsh.wght_freq,fsh.wtprob,
            fsh.ppstf->wghtcontrib,
            fsh.ppstf->wghtcontrib_by_realization,2,
            print_switch);
        }
        else
        {
          wght_tmp=fsh.square_fita(total_wght,fsh.wght_sample_size,
            fsh.tc_wght_freq,fsh.tc_wtprob,
            fsh.ppstf->wghtcontrib,
            fsh.ppstf->wghtcontrib_by_realization,2,
            print_switch);
        }
        break;
      }
    case 7:
      if (fsh.parest_flags(289)==0)
        wght_tmp=logistic_normal_weight_fit(fsh,total_wght);
      else
        if (fsh.parest_flags(301)==0)
        {
          wght_tmp=logistic_normal_weight_fit_heteroscedastic(fsh,total_wght);
        }
        else
        {
          wght_tmp=tail_compressed_weight_logistic_normal_fit_heteroscedastic
            (fsh,total_wght);
        }
      break;
    case 8:
      wght_tmp=square_fit_t_wght(fsh,total_wght);
      break;
    case 9:
      wght_tmp=wght_self_scaling_multinomial_re_multi_rho_multi_var(fsh,
        total_wght,print_switch);
      break;
    case 10:
       {
         int newswitch=1;
//         if (newswitch==0)   //NMD_9Dec2019
         if (fsh.parest_flags(137)==1)
           wght_tmp=old_wght_self_scaling_multinomial_nore(fsh,total_wght,
             print_switch);
         else
           wght_tmp=new_wght_self_scaling_multinomial_nore(fsh,total_wght,
             print_switch);
       }
       break;

    case 11:
      if (fsh.parest_flags(330))
      {
        wght_tmp=wght_dm_nore(fsh,total_wght,print_switch);
      }
      else
      {
        wght_tmp=wght_dm_nore_notc(fsh,total_wght,print_switch);
      }
      break;
    case 12:
       // IIIII
       wght_tmp=new_wght_self_scaling_multinomial_re_multi_rho_multi_var
         (fsh,total_wght,print_switch);
       break;
    default:
      cerr << "Illegal value for parest_flags(139)" << endl; 
      ad_exit(1);
    }
    tmp+=wght_tmp;
    cout << "Weight frequency data " << setprecision(9) << setw(19) << wght_tmp << endl;//NMD
  }
  return tmp;
}

#undef HOME_VERSION
