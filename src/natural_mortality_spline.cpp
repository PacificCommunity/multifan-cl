/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

  void dvar_fish_stock_history::natural_mortality_splines(void)
  {
    pmsd_error();
    int nd=parest_flags(121);
    int nsame=0;
    if (parest_flags(122)>1) nsame=parest_flags(122)-1;
    dvar_vector spline_coffs=age_pars(5)(1,nd);
    dvector nodes(1,nd);
    dvector ages(1,nage-nsame);
#if !defined(NO_MY_DOUBLE_TYPE)
    nodes.fill_seqadd(0,1.0/(nd-1.0L));
#else
    nodes.fill_seqadd(0,1.0/(nd-1.0));
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
    ages.fill_seqadd(0,1.0/(nage-nsame-1.0L));
#else
    ages.fill_seqadd(0,1.0/(nage-nsame-1.0));
#endif
    vcubic_spline_function csf(nodes,spline_coffs);
    dvar_vector nm=csf(ages);
    for (int iy=1;iy<=nyears;iy++)  // Loop over fishing periods
    {
      for (int j=1;j<=nage;j++)         // Loop over age classes
      {
        nat_mort(iy)(1,nage-nsame)=nm;
      }
      if (nsame>0)
      {
        nat_mort(iy)(nage-nsame+1,nage)=nm(nage-nsame);
      }
    }
  }

  void dvar_len_fish_stock_history::natural_mortality_lorenzen(void)
  {
    dvar_vector log_scaled_length(1,nage);
    dvar_vector scaled_length(1,nage);
    dvariable rho=exp(-(vb_coff(3)));
#if !defined(NO_MY_DOUBLE_TYPE)
    dvariable t2=1.-pow(rho,nage-1.0L);
#else
    dvariable t2=1.-pow(rho,nage-1.0);
#endif
    dvariable a=vb_coff(1)/vb_coff(2);
    if (!parest_flags(173))
    {
      for (int j=1;j<=nage;j++)
      {
        dvariable t1=1.-pow(rho,j-1);
        log_scaled_length(j)=log(a+(1.0-a)*t1/t2);
      }
    }
    else
    {
      for (int j=1;j<=nage;j++)
      {
        dvariable t1=1.-pow(rho,j-1);
        scaled_length(j)=a+(1.0-a)*t1/t2;
      }
      int num=parest_flags(173);
      for (int j=1;j<num;j++)
      {
        scaled_length(j+1)+=age_pars(3,j)/vb_coff(2);
      }
      log_scaled_length=log(scaled_length);
    }
    dvar_vector nm=age_pars(5,1)+age_pars(5,2)*log_scaled_length;
    for (int iy=1;iy<=nyears;iy++)  // Loop over fishing periods
    {
      nat_mort(iy)=nm;
    }
    if (pmsd)
    {
      int ns;
      ns=pmsd->num_species;
      for(int is=2;is<=ns;is++)
      {
        int ng=get_nage(is);
	dvar_vector vbc=get_vb_coff(is);
        dvar_vector ap3=get_age_pars_species(is,3);
        dvar_vector ap5=get_age_pars_species(is,5);
        dvar_vector log_scaled_length(1,ng);
        dvar_vector scaled_length(1,ng);
        dvariable rho=exp(-(vbc(3)));
#if !defined(NO_MY_DOUBLE_TYPE)
        dvariable t2=1.-pow(rho,ng-1.0L);
#else
        dvariable t2=1.-pow(rho,ng-1.0);
#endif
        dvariable a=vbc(1)/vbc(2);
        if (!parest_flags(173))
        {
          for (int j=1;j<=ng;j++)
          {
            dvariable t1=1.-pow(rho,j-1);
            log_scaled_length(j)=log(a+(1.0-a)*t1/t2);
          }
        }
        else
        {
          for (int j=1;j<=ng;j++)
          {
            dvariable t1=1.-pow(rho,j-1);
            scaled_length(j)=a+(1.0-a)*t1/t2;
          }
          int num=parest_flags(173);
          for (int j=1;j<num;j++)
          {
            scaled_length(j+1)+=ap3(j)/vbc(2);
          }
          log_scaled_length=log(scaled_length);
        }
        dvar_vector nm=ap5(1)+ap5(2)*log_scaled_length;
        for (int iy=1;iy<=nyears;iy++)  // Loop over fishing periods
        {
          pmsd->nat_mort(is,iy)=nm;
        }
      }
    }
  }

  void dvar_len_fish_stock_history::natural_mortality_double_normal(void)
  {
    pmsd_error();
    dvar_vector nm(1,nage);
    dvariable fp9=age_pars(5,1);
    dvariable fp10=1.0+age_pars(5,2);
    dvariable efp11=exp(age_pars(5,3));
    dvariable a=age_pars(5,4);
    int mmax=nm.indexmax();
    MY_DOUBLE_TYPE jmid=0.5*(1+mmax);
    MY_DOUBLE_TYPE jmid1=1.0-jmid;
    for (int j=1;j<=mmax;j++)
    {
      MY_DOUBLE_TYPE fj=(j-jmid)/jmid1;
      if (fj<=value(fp9))
      {
        //v(j)= pow(2,-square(fj/(fp10*efp11)));
        nm(j)=pow(2,-square(fj/(fp10*efp11)));
      }
      else
      {
        nm(j)=pow(2,-square(fj*efp11/fp10));
      }
    }
    nm/=(1.e-20+sum(nm));
    nm=a+log(nm);
    
    for (int iy=1;iy<=nyears;iy++)  // Loop over fishing periods
    {
      for (int j=1;j<=nage;j++)         // Loop over age classes
      {
        nat_mort(iy)=nm;
      }
    }
  }
