/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


#include "all.hpp"

  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
 long int diff(long int&,const long int&);


void dvar_len_fish_stock_history::cfast_pred_frequency_calc
  (dvar_len_fish_stock_history& cfsh)
{
  #ifdef __ZTC__
    dvar4_array sig = sqr(cfsh.vars);
  #else
    dvar4_array sig = sqrt(cfsh.vars);
  #endif
  tprob.initialize();
  int i;
  int j;
  int js;
  dvariable temp;
  dvariable fdiff2;
  MY_DOUBLE_TYPE prob;
  MY_DOUBLE_TYPE u;
  MY_DOUBLE_TYPE a2 = 35.91908881;
  MY_DOUBLE_TYPE a3 = 785.858644;
  int break_flag=0;
  long int ztmp=gradient_structure::totalbytes();
  long int ztmp1;
  int ir;
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        if (len_sample_size(ir,ip,fi)>0)
        {
          js=1;
          dvar_vector& tprobb=tprob(ir,ip,fi);
          dvar_vector& propp=prop(ir,ip,fi);
          const dvector& mean_len=value(cfsh.mean_length(ir,ip,fi));
          const dvector& sigg=value(sig(ir,ip,fi));
          for (i=1; i<=nlint; i++) /* L1200  */
          {
	    MY_DOUBLE_TYPE& fmidd=fmid(i);
	    for (j=js; j<=nage; j++) /* L1300  */
	    {
	      MY_DOUBLE_TYPE temp=(fmidd-mean_len(j))/sigg(j);
	      if ( fabs(temp) <= 3.03)
	      {
	        if ( fabs(temp) <= 3.00)
	        {
                  prob=exp(-temp*temp/2.e0)/sigg(j);
	        }
	        else if (temp > 0.e0)
	        {
	          u=temp-3.03;
	          prob=(u*u*(a2+a3*u))/sigg(j);
	        }
	        else if (temp < 0.e0)
	        {
	          u=temp+3.03;
	          prob=(u*u*(a2-a3*u))/sigg(j);
	        }
	        tprobb(i) += propp(j) * prob;
	      }
	      else
	      {
	        if (mean_len(j)>fmidd)
	        {
	          //goto la1301;
	          break_flag=1;
	        }
	        else
	        {
	          js=j;
	        }
	        if (break_flag ==1) break;
	      }
	      if (break_flag ==1) break;
	    }
	    break_flag=0;
          }
        }
      }
    }
  }
  ztmp1=diff(ztmp,gradient_structure::totalbytes());

  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        if (len_sample_size(ir,ip,fi)>0)
        {
	  tprob(ir,ip,fi) = tprob(ir,ip,fi) / sum(tprob(ir,ip,fi));
        }
      }
    }
  }
}
