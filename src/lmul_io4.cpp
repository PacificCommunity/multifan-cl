/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


#define HOME_VERSION
#include "all.hpp"
    extern int  _version_no;

    par_ofstream& operator << (par_ofstream& pof, dvar_len_fish_stock_history& fsh)
    {
      if (fsh.parest_flags(197))
      {
        fsh.parest_flags.initialize();
        fsh.parest_flags(30)=1;
      }
      //if (fsh.parest_flags(155)>0)  fsh.parest_flags(156)=1; 
      fsh.parest_flags(200)=_version_no;
      pof << "# The parest_flags" << endl;
      pof << fsh.parest_flags << endl;
      if (fsh.pmsd)
      {
        if (fsh.pmsd->num_species>1)
        {
          pof <<  "# multi-species parest flags"<<endl;
          {
            pof <<  fsh.pmsd->parest_flags <<endl;
          }
        }
      }
      pof << "# The number of age classes" << endl;
      pof << fsh.nage << endl;
      if (fsh.pmsd)
      {
        pof << "# Multi-species the number of age classes" << endl;
        pof << fsh.pmsd->nage;
      }
      pof << * (dvar_fish_stock_history *) &fsh;
      pof << "# The logistic normal parameters" << endl;
      pof << "# log_length variance    length rho     "
         " log_length_dof     length_psi   length_exp" << endl;
      pof << "   " << fsh.log_length_variance 
          << "         " << fsh.length_rho 
          << "         " << fsh.log_length_dof 
          << "         " << fsh.length_psi 
          << "         " << fsh.length_exp << endl;

      pof << "# log_weight variance    weight rho     "
         " log_weight_dof   weight_psi   weight exp " << endl;
      pof << "   " << fsh.log_weight_variance 
          << "         " << fsh.weight_rho 
          << "         " << fsh.log_weight_dof 
          << "         " << fsh.weight_psi 
          << "         " << fsh.weight_exp << endl;
      pof << "# length_tot_exp weight_tot_exp  " << endl;
      pof <<  fsh.length_tot_exp << "    " 
          << fsh.weight_tot_exp << endl;
      pof << endl << "# maturity at length" << endl;
      pof <<  fsh.cpmature_at_length;
      if (fsh.pmsd)    //NMD_aug30_21018
      {
        int numsp=fsh.pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          pof << fsh.pmsd->cpmature_at_length(i);
        }
      }    //NMD_aug30_21018
      
      pof << endl << "# The von Bertalanffy parameters" << endl;

      pof << fsh.vb_coff(1) << "  " << fsh.fmin1 << "  " << fsh.fmax1 << endl;
      pof << fsh.vb_coff(2) << "  " << fsh.fminl << "  " << fsh.fmaxl << endl;
      pof << fsh.vb_coff(3) << "  " << fsh.rhomin << "  " << fsh.rhomax << endl;
      pof << endl << "# Extra par for Richards" << endl;
      pof << fsh.vb_coff(4) << endl;

      if (fsh.pmsd)
      {
        int numsp=fsh.pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
        pof << endl << "# The von Bertalanffy parameters "
            " Species" << i << endl;
        pof << fsh.pmsd->vb_coff(i,1) << "  " << fsh.pmsd->fmin1(i) 
            << "  " << fsh.pmsd->fmax1(i) << endl;
        pof << fsh.pmsd->vb_coff(i,2) << "  " << fsh.pmsd->fminl(i) 
            << "  " << fsh.pmsd->fmaxl(i) << endl;
        pof << fsh.pmsd->vb_coff(i,3) << "  " << fsh.pmsd->rhomin(i) 
            << "  " << fsh.pmsd->rhomax(i) << endl;
        if (_version_no > 1039)
        {
          pof << endl << "# Extra par for Richards species " 
              << i  << endl;
          pof << fsh.pmsd->vb_coff(i,4) << endl;
        }
      }
    }
    //cout << "# First Length bias parameters" << endl;
    int i;
    //for (i=fsh.vb_bias.indexmin();i<=fsh.vb_bias.indexmax();i++)
    //{
    //  cout <<"lmul_io4.cpp "  << fsh.vb_bias(i) << endl;
    //}
    pof << "# First Length bias parameters" << endl;
    pof << fsh.vb_bias << endl;

    pof << "# Common first Length bias flags" << endl;
    pof << fsh.common_vb_bias << endl;
    pof << "# Common first Length bias coffs" << endl;
    pof << fsh.common_vb_bias_coffs << endl;

    if (!allocated(fsh.Ninit_standard))
    {
      fsh.Ninit_standard.allocate(1,fsh.num_regions,1,fsh.last_real_year);
      fsh.Ninit_standard.initialize();
    }
    pof << "#Recruitment standard" << endl;
    pof << fsh.num_regions << "  "  << fsh.last_real_year << endl;
    pof << fsh.Ninit_standard << endl;
    if (!allocated(fsh.Ninit_orthogonal))
    {
      fsh.Ninit_orthogonal.allocate(1,fsh.num_regions,1,fsh.last_real_year);
      fsh.Ninit_orthogonal.initialize();
    }
    pof << "#Recruitment orthogonal" << endl;
    pof << fsh.Ninit_orthogonal  << endl;
    pof << "# Seasonal growth parameters" << endl;
    pof << fsh.sv << endl;
    if (fsh.pmsd)
    {
      pof << "# Extra multi-species growth parameters" << endl;
      pof << fsh.pmsd->sv << endl;
    }
    pof << "# Cohort specific growth deviations" << endl;
    pof <<  fsh.growth_dev << endl;

    pof << "# Variance parameters" << endl;
    pof << fsh.var_coff(1) << "  " << fsh.vmin1 << "  " << fsh.vmax1 << endl;
    pof << fsh.var_coff(2) << "  " << fsh.vminl << "  " << fsh.vmaxl << endl;


      if (fsh.pmsd)
      {
        int numsp=fsh.pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          pof << fsh.pmsd->var_coff(i,1) << "  " << fsh.pmsd->vmin1(i) 
              << fsh.pmsd->vmax1(i)<<endl;
          pof << fsh.pmsd->var_coff(i,2) << "  " << fsh.pmsd->vminl(i) 
              << "  " << fsh.pmsd->vmaxl(i)<<endl;
        }
      }
      pof << "# kludged_equilib_coffs" <<endl;
      pof << fsh.kludged_equilib_coffs <<endl;
      pof << "# kludged_equilib_level_coffs" << endl;
      pof << fsh.kludged_equilib_level_coffs << endl;
      pof << "# new orthogonal coefficients" << endl;
      pof << fsh.num_new_weights << endl;
      if (allocated(fsh.new_orth_recr))
      {
        pof << fsh.new_orth_recr << endl;
      }
      if (fsh.pmsd)
      {
        int num_species=fsh.pmsd->num_species;
        for (int is=2;is<=num_species;is++)
        {
          pof << fsh.pmsd->num_new_weights(is) << endl;
          if (allocated(fsh.pmsd->new_orth_recr(is)))
          {
            pof << fsh.pmsd->new_orth_recr(is) << endl;
          }
        }
      } 

      pof << "# The number of mean constraints" << endl;
      pof << fsh.nfmbound << endl;
      for (i=1; i<=fsh.nfmbound; i++)
      {
	pof <<  fsh.ifper(i) << "  " << fsh.ifinc(i)   << "  "
	    << fsh.iageclass(i) << "  " << fsh.fmmin(i) << "  "
	    << fsh.fmmax(i) << "  " << fsh.pdown(i) << "  "
	    << fsh.pup(i) << endl;
      }
      pof << "# The diffusion coefficients" << endl;
      pof << fsh.D << endl;

      pof << "# First year in model" << endl;
      pof << fsh.year1 << endl;

      pof << "# MULTIFAN-CL compilation version number" << endl;
      pof << fsh.current_version << endl;
      
      pof << "# The grouped_catch_dev_coffs flag" << endl;
      if (allocated(fsh.grouped_catch_dev_coffs))
      {
        pof << " 1" << endl;
        pof << "# The grouped_catch_dev_coffs" << endl;
        pof << fsh.grouped_catch_dev_coffs << endl;
      }
      else
        pof << " 0" << endl;
      if (fsh.age_flags(92)>0)
        pof <<"# implict total mortality constraints penalty value is "
            << fsh.imppen  << endl;
      pof <<  "# Historical_flags" << endl;      //NMD_24Mar2015
      pof << "# The parest_flags" << endl;
      pof << fsh.parest_flags << endl;
      pof <<  "# Do multsp Historical_flags" << endl;      //NMD_24Mar2015
      if (fsh.pmsd)
      {
        if (fsh.pmsd->num_species>1)
        {
          pof <<  "# multi-species parest flags"<<endl;
          {
            pof <<  fsh.pmsd->parest_flags <<endl;
          }
        }
      }

      pof << fsh.age_flags << endl;
      if (fsh.pmsd)
      {
        if (fsh.pmsd->num_species>1)
        {
          pof <<  "# multi-species age flags"<<endl;
          {
            pof <<  fsh.pmsd->age_flags <<endl;
          }
        }
      }
      pof << "# fish flags"<<endl;
      pof << setw(4) << fsh.fish_flags << endl;
      pof << "#    End_historical_flags" << endl;      //NMD_18Dec2014

      return pof;
    }


#undef HOME_VERSION
