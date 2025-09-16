/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE



#include "all.hpp"


void  dvar_len_fish_stock_history::big_fish_catch(void)
{
  int icut=parest_flags(180);
  dvar_vector mlbs(1,num_fisheries);
  dvar_matrix tlbs(1,num_fisheries,icut,nlint);
  int ir;
  int i;
  for (i=1;i<=num_fisheries;i++)
  {
    for (int ii=icut;ii<=nlint;ii++)
    {
      tlbs(i,ii)=lengthbsel(i,ii);
    }
  }
  for (i=1;i<=num_fisheries;i++)
  {
    tlbs(i)/=mean(lengthbsel(i));
  }
  tlbs=log(tlbs);
  const dvector& cfmid=cube(fmid(icut,nlint));
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        //dvar_vector& fs=incident_sel(ir,ip,fi);
        int i=parent(ir,ip,fi);
        dvar_vector& fs=tlbs(i);
        MY_DOUBLE_TYPE eff=effort(ir,ip,fi);
        const prevariable& cat=catchability(ir,ip,fi);
        dvar_vector fm=exp(fs+eff+cat);
        dvar_vector z=fm+exp(get_nat_mort_coff_region(ir)+fraction(ir,ip));
        const dvar_vector& ld=len_dist(ir,ip)(icut,nlint);
        dvar_vector C=
          elem_prod(elem_div(fm,z),elem_prod(1-exp(-z),ld));
        RB(ir,ip,fi)=sum(elem_prod(C,cfmid));
      }
    }
  }
}
