/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"
#ifdef __MSC__
  MY_DOUBLE_TYPE mfexp(const MY_DOUBLE_TYPE& );
  dvariable mfexp(const prevariable& );
  dvar_vector mfexp(const dvar_vector& );
  dvector mfexp(const dvector& );
#else
  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
  dvariable mfexp(BOR_CONST prevariable& );
  //dvar_vector mfexp(dvar_vector& );
  //dvector mfexp(dvector& );
#endif

dvariable dvar_fish_stock_history::fit_implicit_catchability(void)
{
  dvariable f=0.0;
  dvar_vector fv(1,num_fisheries);
  int i,fi;
  for (i=1;i<=num_fisheries;i++)
  {
    dvar_vector q(1,num_fish_times(i));
    int rr=realization_region(i,1);
    for (fi=1;fi<=num_fish_times(i);fi++)
    {
      int rr=realization_region(i,fi);
      int rp=realization_period(i,fi);
      q(fi)=catchability(rr,rp,realization_incident(i,fi));
    }
  
    dvar_vector& iq = implicit_catchability(i);
    dvariable iqm=mean(iq);
    dvariable qm=mean(q);

    //dvar_vector r=iq-q-(iqm-qm);
    dvar_vector r=iq-q;
    dvar_vector r2=square(r);


    fv(i)=-sum(log(exp(-r2*50)+.05/(1+50*r2)));

    f+=fv(i);

  }
  return f;
}

