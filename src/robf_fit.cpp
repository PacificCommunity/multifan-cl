/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#define HOME_VERSION
#include <fvar.hpp>

dvariable robust_freq_fit(d3_array & obs, dvar3_array& pred,
  dvar3_array& weights)
{
  RETURN_ARRAYS_INCREMENT(); //Need this statement because the function
			     //returns a variable type
  dvariable v_hat;
  MY_DOUBLE_TYPE width=3.0;
  MY_DOUBLE_TYPE pcon=0.05;
  MY_DOUBLE_TYPE width2=width*width;
  MY_DOUBLE_TYPE a=0.7;
  MY_DOUBLE_TYPE a2;

  a2=a*a;
  dvar3_array diff = obs-pred;     // These are the residuals
  dvar3_array diff2 = square(diff); // These are the squared residuals
  dvar3_array wdiff2 = elem_div(diff2,weights); // These are the weighted
						// squared residuals
  v_hat = mean(wdiff2)+1.e-80; // add 1.e-80 so that a perfect fit wont't
			      // produce log(0).
  MY_DOUBLE_TYPE b=2.*pcon/(width*sqrt(3.14159));  // This is the weight for the
					   // "robustifying" term
  dvariable log_likelihood = -sum(log((1.-pcon)*exp(wdiff2/(-2*a2*v_hat))
    + b/(1.+square(wdiff2/(width2*a2*v_hat)))));

  dvariable var_part = 
    0.5*(sum(log(weights))+double(size_count(diff))*log(a2*v_hat));
  //log_likelihood += 0.5*size_count(diff)*log(a2*v_hat);
  log_likelihood += var_part;

  RETURN_ARRAYS_DECREMENT(); // Need this to decrement the stack increment
			     // caused by RETURN_ARRAYS_INCREMENT();
  return(log_likelihood);
}

dvariable robust_freq_fit(d3_array & obs, dvar3_array& pred,
  dvar3_array& weights, MY_DOUBLE_TYPE& min_sig)
{
  RETURN_ARRAYS_INCREMENT(); //Need this statement because the function
			     //returns a variable type
  dvariable v_hat;
  MY_DOUBLE_TYPE width=3.0;
  MY_DOUBLE_TYPE pcon=0.05;
  MY_DOUBLE_TYPE width2=width*width;
  MY_DOUBLE_TYPE a=0.7;
  MY_DOUBLE_TYPE a2;

  a2=a*a;
  dvar3_array diff = obs-pred;     // These are the residuals
  dvar3_array diff2 = square(diff); // These are the squared residuals
  dvar3_array wdiff2 = elem_div(diff2,weights); // These are the weighted
						// squared residuals
  v_hat = mean(wdiff2)+1.e-80; // add 1.e-80 so that a perfect fit wont't
			      // produce log(0).
  MY_DOUBLE_TYPE b=2.*pcon/(width*sqrt(3.14159));  // This is the weight for the
					   // "robustifying" term
  dvariable log_likelihood = -sum(log((1.-pcon)*exp(wdiff2/(-2*a2*v_hat))
    + b/(1.+square(wdiff2/(width2*a2*v_hat)))));

  dvariable var_part = 0.5*(sum(log(weights))+
    (MY_DOUBLE_TYPE)( size_count(diff))*log(min_sig+a2*v_hat));
  //log_likelihood += 0.5*size_count(diff)*log(a2*v_hat);
  log_likelihood += var_part;

  RETURN_ARRAYS_DECREMENT(); // Need this to decrement the stack increment
			     // caused by RETURN_ARRAYS_INCREMENT();
  return(log_likelihood);
}

#undef HOME_VERSION

