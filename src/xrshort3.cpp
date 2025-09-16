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
  dvariable mfexp(prevariable& );
  //dvar_vector mfexp(dvar_vector& );
  //dvector mfexp(dvector& );
#endif

  extern dvar_vector * psv;
  extern dvar_len_fish_stock_history * pcfsh;

void dvar_fish_stock_history::total_fish_mort_and_survival(int ir,int ip,
  int yr,const dvar_vector & tmp)
{
  // calculate the total fishing mortality and survival rates for 
  // fishing period ip in region ir
  //cout <<"xrshort3.cpp " << fraction(ir) << endl;
  dvar_vector& tm=tot_mort(ir,ip);
  // check if we need this!
  tm.initialize();
  dvar_matrix& fm=fish_mort(ir,ip);
  for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
  {                                         // incidents for this period
    int i=parent(ir,ip,fi);
#if !defined(NO_MY_DOUBLE_TYPE)
    fm(fi)+=(tmp(i)-1.0L)*
#else
    fm(fi)+=(tmp(i)-1.0)*
#endif
      (rel_biomass_by_region(ir,yr)-rel_biomass_by_region(ir,1));
    tm+=mfexp(fm(fi));
  }
  tm+=mfexp(get_nat_mort_region(ir)(year(ir,ip))+fraction(ir,ip));
  survival(ir,ip)=mfexp(-tm);
}

