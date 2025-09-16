/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#define HOME_VERSION
#include "all.hpp"

dvariable robust_freq_partial_fit(d4_array & obs, dvar4_array& pred,
  dvar4_array& weights,dvar_len_fish_stock_history& fsh)
{
  RETURN_ARRAYS_INCREMENT(); //Need this statement because the function
                             //returns a variable type
  dvariable v_hat=0.;
  MY_DOUBLE_TYPE width=3.0;
  MY_DOUBLE_TYPE pcon=0.05;
  MY_DOUBLE_TYPE width2=width*width;
  MY_DOUBLE_TYPE a=0.7;
  MY_DOUBLE_TYPE a2;

  dvar4_array diff;     // These are the residuals
  diff.allocate(obs);
  dvar4_array diff2;     // These are the residuals
  diff2.allocate(obs);
  dvar4_array wdiff2;     // These are the residuals
  wdiff2.allocate(obs);
  a2=a*a;

  int tot=0;
  int ir;
  for (ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++) 
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        if (fsh.len_sample_size(ir,ip,fi)>0)
        {
          diff(ir,ip,fi) = obs(ir,ip,fi)-pred(ir,ip,fi);//These are the resids
          diff2(ir,ip,fi) = square(diff(ir,ip,fi)); // These are the squared 
                                                    // residuals
          wdiff2(ir,ip,fi) = elem_div(diff2(ir,ip,fi),weights(ir,ip,fi));
                        //These are the weighted  squared residuals
          v_hat += sum(wdiff2(ir,ip,fi));
          tot+=size_count(wdiff2(ir,ip,fi));
        }
      }
    }
  }
  v_hat=v_hat/double(tot);
  MY_DOUBLE_TYPE b=2.*pcon/(width*sqrt(3.14159));  // This is the weight for the
                                           // "robustifying" term
  dvariable log_likelihood =0.;
  dvariable var_part =  0.;
  for (ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++) 
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                             // incidents for this period
        if (fsh.len_sample_size(ir,ip,fi)>0)
        {
          log_likelihood += -sum(log((1.-pcon)*exp(wdiff2(ir,ip,fi)/
             (-2*a2*v_hat))
            + b/(1.+square(wdiff2(ir,ip,fi)/(width2*a2*v_hat)))));
          var_part +=
            0.5*(sum(log(weights(ir,ip,fi))))+
           (MY_DOUBLE_TYPE)( size_count(diff(ir,ip,fi)))*log(a2*v_hat);
        }
      }
    }
  }  
  log_likelihood += var_part;
  RETURN_ARRAYS_DECREMENT(); // Need this to decrement the stack increment
  return(log_likelihood);
}

#undef HOME_VERSION

