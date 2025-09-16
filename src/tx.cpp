/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

void get_bh_alpha_and_beta(dvar_len_fish_stock_history& fsh,
  dvariable& alpha,dvariable& beta)
{
  int lag=fsh.age_flags(147);
  int numflag=fsh.age_flags(149);

  int ibh;
  int jbh;
  if(fsh.age_flags(199)>0)
  {
    ibh=fsh.last_real_year - fsh.age_flags(199)+1;
    jbh=fsh.last_real_year - fsh.age_flags(200);
  }
  else
  {
    ibh=1;
    jbh=fsh.last_real_year;
  }

  dvar_vector rec(1,fsh.last_real_year);
  dvar_vector bio(1,fsh.last_real_year);
  bio.initialize();
  rec.initialize();
  dvar_vector & pmature = fsh.pmature;
  if (!value(sum(pmature)))
  {
    for (int j=1;j<=fsh.nage;j++)
    {
      pmature=1.0;
    }
  } 
  for (int i=1;i<=fsh.last_real_year;i++)  
  {
    // numflag determines if spawning stock is in numbers or weight
    // **********************************************************
    // **********************************************************
    // If effort is turned off need to use the old numbers at age
    // which are (hopefully) stored in fsh.Nsave. we should figure out
    // a check for this so that the code wonr't silently break at
    // some point
    // DF May 20 2013
    // **********************************************************
    // **********************************************************
    if (fsh.af170q0==0)
      bio(i)+= calculate_the_biomass(i,fsh,fsh.N,numflag,pmature);
    else
      bio(i)+= calculate_the_biomass(i,fsh,fsh.Nsave,numflag,pmature);
    for (int ir=1;ir<=fsh.num_regions;ir++)
    {
      rec(i)+=fsh.exp_N(ir,i,1);
    }
  }
  dvariable b1=bio(1);
  dvariable steepness;
  dvariable phi;

  int cs=fsh.get_current_species();

  if (cs==1)
  {
    phi=fsh.phi;
    steepness=fsh.sv(29);
    beta=b1*(fsh.sv(21)+0.001);
  }
  else
  {
    phi=fsh.pmsd->phi(cs);
    steepness=fsh.pmsd->sv(cs,29);
    beta=b1*(fsh.pmsd->sv(cs,21)+0.001);
  }

  dvariable logdelta=log(4.0*steepness/(phi*(1.0-steepness)));
  dvariable loga=logdelta+log(beta);
  alpha=exp(loga);
}

dvariable get_bh_recruitment_for_projections(int i,
  dvar_len_fish_stock_history& fsh,dvariable& alpha,dvariable& beta)
{
  int numflag=fsh.age_flags(149);
  dvariable bio;
  // **********************************************************
  // **********************************************************
  // If effort is turned off need to use new number at age
  // which are (hopefully) stored in fsh.N_q0. 
  // DF May 20 2013
  // **********************************************************
  // **********************************************************
  if (fsh.af170q0==0)
    bio= calculate_the_biomass(i,fsh,fsh.N,numflag,fsh.pmature);
  else
    bio= calculate_the_biomass(i,fsh,fsh.N_q0,numflag,fsh.pmature);
  dvariable rec=alpha*bio/(beta+bio);
  if (fsh.age_flags(161))
  {
    rec*=exp(0.5*fsh.bh_variance);
  }
  return rec;
}
