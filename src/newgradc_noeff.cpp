/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"
#include "here.h"
#include "globals.h"

 // _mf_cl_globals::_mf_cl_globals() : dep_var_labels(1,MAX_DEPVAR_LABELS);

  dvariable get_dd_agedif(dvar_vector& sv,dvar_vector& rs);
  dvar_vector sbcalc(dvar_len_fish_stock_history& fsh);
  void do_average_exploitation(dvar_len_fish_stock_history& fsh,int& ii,
    int& num_grads,dvar_vector& dep_vars,ofstream& ofl);

  void depvars_bounds_check(int i,int n);

  int dep_gradients_calc_noeff(dvar_len_fish_stock_history& fsh)
  {
    //HERE
    ofstream tmpout("tmp_noeff.out");
    ofstream ofl("deplabel_noeff.tmp");
    //HERE
    int num_grads=20000;
    int ii=1;
    //HERE
    dvar_vector dep_vars(1,num_grads);
    int i,ir;
    int ny=fsh.nyears;
    dvar_matrix reg_rbio=unnormalized_reg_vbiocalc(fsh,1);  // 1 for no effort
    dvar_matrix adult_reg_rbio=unnormalized_adult_reg_vbiocalc(fsh,1);
    dvar_vector adult_rbio=adult_vbiocalc(fsh,1);   // 1 for no effort

    // Abbreviated calculation based on parest_flags(37) NMD_24july2023
    if (fsh.parest_flags(37)!=1)
    {
//  - start of abbreviation
    
    // do relative biomass gradients
      dvar_vector rbio=vbiocalc(fsh,1);   // 1 for no effort
      for (i=1;i<=ny;i++)  
      {
        dep_vars(ii) << log(rbio(i));
        ofl << dep_vars(ii++) << endl;
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("rbio_noeff(")  
           + str(i) + adstring(")") << endl;
      }

      for (i=1;i<=ny;i++)  
      {
        dep_vars(ii) << log(adult_rbio(i));
        ofl << dep_vars(ii++) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("adult_rbio_noeff(")  
           + str(i) + adstring(")") << endl;
      }

      for (ir=1;ir<=fsh.num_regions;ir++)
      {
        for (int i=1;i<=fsh.nyears;i++)  
        {
          dep_vars(ii) << log(reg_rbio(ir,i));
          ofl << dep_vars(ii++) << endl; 
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("ln_reg_bio_noeff(") +str(ir) +adstring(",") 
             + str(i) + adstring(")") << endl;
        }
      }

      for (ir=1;ir<=fsh.num_regions;ir++)
      {
        for (int i=1;i<=fsh.nyears;i++)  
        {
          dep_vars(ii) << log(adult_reg_rbio(ir,i));
          ofl << dep_vars(ii++) << endl; 
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("ln_adult_reg_bio_noeff(") +str(ir) +adstring(",") 
             + str(i) + adstring(")") << endl;
        }
      }
    }
//  - end of abbreviation

//   average adult_biomass over specified calendar years  //NMD_3May2021
    dvariable avg_unfished_adult_biomass=0.0;
    int num_for_average=10;
    if (fsh.parest_flags(60)!=0) num_for_average=fsh.parest_flags(60);
    int af57=fsh.age_flags(57);
    int yr1st=ny-((num_for_average+1)*af57)+1;
    int yrlast=ny-af57;
    int nfa_calc = num_for_average * af57;
    if (!fsh.pmsd)  //- single species case
    {
      dvariable avg_adult_biomass=0.0;
      for (i=yr1st;i<=yrlast;i++)
      {
        avg_adult_biomass+=adult_rbio(i);
      }
      avg_adult_biomass/=nfa_calc;      
      dep_vars(ii) << log(avg_adult_biomass);
      ofl << setprecision(8) << dep_vars(ii++) << endl;
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("average_adult_rbio_noeff(")
        + str(nfa_calc) + adstring("_periods)") << endl;

      // duplicate as needed for the SB_recent ratio dep_var  //NMD_21jun2022
      dep_vars(ii) << log(avg_adult_biomass);
      ofl << setprecision(8) << dep_vars(ii++) << endl;
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("average_adult_rbio_noeff(")
        + str(nfa_calc) + adstring("_periods)") << endl;  //NMD_21jun2022
    }
    else if (!sum(column(fsh.pmsd->species_flags,2))) // - multi-species
    {
      for (int is=1;is<=fsh.pmsd->num_species;is++)
      {
        int mmin,mmax;
        mmin=fsh.pmsd->region_bounds(is,1);
        mmax=fsh.pmsd->region_bounds(is,2);
        dvariable avg_adult_biomass=0.0;
        for (i=yr1st;i<=yrlast;i++)
        {
          for (int ir=mmin;ir<=mmax;ir++)
          {
            avg_adult_biomass+=adult_reg_rbio(ir,i);
          }
        }
        avg_adult_biomass/=nfa_calc;
        dep_vars(ii) << log(avg_adult_biomass);
        ofl << setprecision(8) << dep_vars(ii++) << endl;
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("average_adult_rbio_noeff(")
          + str(nfa_calc) + adstring("_periods_species_")
          + str(is) + adstring(")") << endl;
	
      // duplicate as needed for the SB_recent ratio dep_var  //NMD_21jun2022
        dep_vars(ii) << log(avg_adult_biomass);
        ofl << setprecision(8) << dep_vars(ii++) << endl;
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("average_adult_rbio_noeff(")
          + str(nfa_calc) + adstring("_periods_species_")
          + str(is) + adstring(")") << endl;
      }
    }
    else        // multi-sex
    {
      ivector sf2=column(fsh.pmsd->species_flags,2);
      for (int is=1;is<=fsh.pmsd->num_species;is++)
      {
        if (sf2(is))   // get female
        {
          int mmin,mmax;
          mmin=fsh.pmsd->region_bounds(is,1);
          mmax=fsh.pmsd->region_bounds(is,2);
          dvariable avg_adult_biomass=0.0;
          for (i=yr1st;i<=yrlast;i++)
          {
            for (int ir=mmin;ir<=mmax;ir++)
            {
              avg_adult_biomass+=adult_reg_rbio(ir,i);
            }
          }
          avg_adult_biomass/=nfa_calc;
          dep_vars(ii) << log(avg_adult_biomass);
          ofl << setprecision(8) << dep_vars(ii++) << endl;
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("average_adult_rbio_noeff(")
            + str(nfa_calc) + adstring("_periods_sex_")
            + str(is) + adstring(")") << endl;

      // duplicate as needed for the SB_recent ratio dep_var  //NMD_21jun2022
          dep_vars(ii) << log(avg_adult_biomass);
          ofl << setprecision(8) << dep_vars(ii++) << endl;
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("average_adult_rbio_noeff(")
            + str(nfa_calc) + adstring("_periods_sex_")
            + str(is) + adstring(")") << endl;

        }
      }
    }      //NMD_3May2021

    //  ***************************************************************
    int ndep=ii-1;
    mfglobals.dep_vars_values.allocate(1,ndep);
    int ij;
    for (ij=1;ij<=ndep;ij++)
    {
      mfglobals.dep_vars_values(ij) = value(dep_vars(ij));
    }
    return ndep;
  }
