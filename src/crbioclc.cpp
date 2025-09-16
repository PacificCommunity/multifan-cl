/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
extern dvector mean_weights_kludge;

/*
dvector cbiocalc(dvar_len_fish_stock_history& fsh)
{
  dvector bio(1,fsh.nyears);
  bio(1)=sum(elem_prod(exp(value(fsh.num_fish(1))),
    pow(value(fsh.mean_length(1,1)),3)));
  int i=2;
  for (int ip=2;ip<=fsh.num_fish_periods;ip++)
  {
    if (fsh.year(ip)>fsh.year(ip-1))
    {
      bio(i++)=sum(elem_prod(exp(value(fsh.num_fish(ip))),
	       pow(value(fsh.mean_length(ip,1)),3)));
    }
  }
  return bio;
}

dvector cbiocalc2(dvar_len_fish_stock_history& fsh)
{
  dvector bio(1,fsh.nyears);
  dvector dot(1,fsh.nage);
  dot.fill_seqadd(1,0);
  dot(1)=0.0; 
  dvector tmp=elem_prod(exp(value(fsh.num_fish(1))),
    pow(value(fsh.mean_length(1,1)),3));
  bio(1)=dot*tmp;
  int i=2;
  for (int ip=2;ip<=fsh.num_fish_periods;ip++)
  {
    if (fsh.year(ip)>fsh.year(ip-1))
    {
      dvector tmp=elem_prod(exp(value(fsh.num_fish(ip))),
	pow(value(fsh.mean_length(ip,1)),3));
      bio(i++)=dot*tmp;
    }
  }
  return bio;
}

dvector cbiocalc_wt(dvar_len_fish_stock_history& fsh)
{
  dvector bio(1,fsh.nyears);
  bio(1)=exp(value(fsh.num_fish(1)))* mean_weights_kludge;
  int i=2;
  for (int ip=2;ip<=fsh.num_fish_periods;ip++)
  {
    if (fsh.year(ip)>fsh.year(ip-1))
    {
      bio(i++)=exp(value(fsh.num_fish(ip)))* mean_weights_kludge;
    }
  }
  return bio;
}
*/
