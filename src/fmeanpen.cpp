/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

dvariable dvar_len_fish_stock_history::fmeanpen()
{
  dvariable fpen=0.;

  for (int ii=1; ii<=nfmbound; ii++)
  {
    int i=ifper(ii);
    int jj=ifinc(ii);
    int rr=realization_region(i,jj);
    int rp=realization_period(i,jj);
    int ri=realization_incident(i,jj);
    int j=iageclass(ii);

    if ( mean_length(rr,rp,ri,j) < fmmin(ii) ) 
    {
      fpen+=1000.0*
        square(mean_length(rr,rp,ri,j)-fmmin(ii));
    }
    else if ( mean_length(rr,rp,ri,j)  > fmmax(ii))
    {
      fpen+=1000.0*
	square(mean_length(rr,rp,ri,j)-fmmax(ii));
    }
  }
  return fpen;
}
