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

  _mf_cl_globals::_mf_cl_globals() : dep_var_labels(1,MAX_DEPVAR_LABELS)
  {
     max_depvar_labels=MAX_DEPVAR_LABELS;
  } 

  dvariable get_dd_agedif(dvar_vector& sv,dvar_vector& rs);
  dvar_vector sbcalc(dvar_len_fish_stock_history& fsh);
  void do_average_exploitation(dvar_len_fish_stock_history& fsh,int& ii,
    int& num_grads,dvar_vector& dep_vars,ofstream& ofl);

   void depvars_bounds_check(int i,int n)
   {
     if (i > n)
     {
       cerr << "Need to increase num_grads in dep_gradients_calc2"
               " current value is " << n << endl;
       exit(1);
     }
   }


  int dep_gradients_calc2(dvar_len_fish_stock_history& fsh)
  {
    //HERE
    ofstream tmpout("tmp.out");
    ofstream ofl("deplabel.tmp");
    //HERE
    int num_grads=50000;
    int ii=1;
    //HERE
    dvar_vector dep_vars(1,num_grads);
    // do relative biomass gradients
    dvar_vector rbio=vbiocalc(fsh); 
    int ny=0;
    ny=fsh.nyears;
    int i;
    int ir;

    dvar_vector adult_rbio=adult_vbiocalc(fsh); 
    dvar_matrix adult_reg_rbio=unnormalized_adult_reg_vbiocalc(fsh);

    // Abbreviated calculation based on parest_flags(37) NMD_24july2023
    if (fsh.parest_flags(37)!=1)
    {
//  - start of abbreviation

    // tag loss stuff 17 aug 2021
    // for now just do tagmort(1)
      if (fsh.parest_flags(360))
      {
        dep_vars(ii) << fsh.tagmort(1);
        ofl << dep_vars(ii++) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("tagmort(")  
           + str(1) + adstring(")") << endl;
      }
   
      for (i=1;i<=ny;i++)  
      {
        dep_vars(ii) << log(rbio(i));
        ofl << dep_vars(ii++) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("rbio(")  
           + str(i) + adstring(")") << endl;
      }

      for (i=1;i<=ny;i++)  
      {
        dep_vars(ii) << log(adult_rbio(i));
        ofl << dep_vars(ii++) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("adult_rbio(")  
           + str(i) + adstring(")") << endl;
      }

      dvar_matrix reg_rbio=unnormalized_reg_vbiocalc(fsh);
      for (ir=1;ir<=fsh.num_regions;ir++)
      {
        for (int i=1;i<=ny;i++)  
        {
          dep_vars(ii) << log(reg_rbio(ir,i));
          ofl << dep_vars(ii++) << endl; 
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("ln_reg_bio(") +str(ir) +adstring(",") 
             + str(i) + adstring(")") << endl;
        }
      }

      for (ir=1;ir<=fsh.num_regions;ir++)
      {
        for (int i=1;i<=ny;i++)  
        {
          dep_vars(ii) << log(adult_reg_rbio(ir,i));
          ofl << dep_vars(ii++) << endl; 
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("ln_adult_reg_bio(") +str(ir) +adstring(",") 
             + str(i) + adstring(")") << endl;
        }
      }
    }
//  - end of abbreviation

    int num_for_average=fsh.age_flags(57);
    if (!fsh.pmsd)  //- single species case
    {
      dvariable avg_adult_biomass=0.0;
      for (i=ny-num_for_average+1;i<=ny;i++)
      {
        avg_adult_biomass+=adult_rbio(i);
      }
      avg_adult_biomass/=num_for_average;
      dep_vars(ii) << log(avg_adult_biomass);
      ofl << setprecision(8) << dep_vars(ii++) << endl;
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("adult_rbio(")
        + str(ny) + adstring(")") << endl;
    }
    else if (!sum(column(fsh.pmsd->species_flags,2))) // - multi-species
    {
      for (int is=1;is<=fsh.pmsd->num_species;is++)
      {
        int mmin,mmax;
        mmin=fsh.pmsd->region_bounds(is,1);
        mmax=fsh.pmsd->region_bounds(is,2);
        dvariable avg_adult_biomass=0.0;
        for (i=ny-num_for_average+1;i<=ny;i++)
        {
          for (int ir=mmin;ir<=mmax;ir++)
          {
            avg_adult_biomass+=adult_reg_rbio(ir,i);
          }
        }
        avg_adult_biomass/=num_for_average;
        dep_vars(ii) << log(avg_adult_biomass);
        ofl << setprecision(8) << dep_vars(ii++) << endl;
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("adult_rbio(")
          + str(ny) + adstring("_species_") + str(is) + adstring(")") << endl;
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
          for (i=ny-num_for_average+1;i<=ny;i++)
          {
            for (int ir=mmin;ir<=mmax;ir++)
            {
              avg_adult_biomass+=adult_reg_rbio(ir,i);
            }
          }
          avg_adult_biomass/=num_for_average;
          dep_vars(ii) << log(avg_adult_biomass);
          ofl << setprecision(8) << dep_vars(ii++) << endl;
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("adult_rbio(")
            + str(ny) + adstring("_sex_") + str(is) + adstring(")") << endl;
        }
      }
    }

//   average recent adult_biomass  //NMD_22jun2022
    num_for_average=4;
    if (fsh.parest_flags(59)!=0) num_for_average=fsh.parest_flags(59);
    int af57=fsh.age_flags(57);
    int yr1st=ny-((num_for_average)*af57)+1;
    int yrlast=ny;
    int nfa_calc = num_for_average * af57;
    dvariable avg_adult_biomass=0.0;

    if (!fsh.pmsd)  //- single species case
    {
      for (i=yr1st;i<=yrlast;i++)
      {
        avg_adult_biomass+=adult_rbio(i);
      }
      avg_adult_biomass/=nfa_calc;
      dep_vars(ii) << log(avg_adult_biomass);
      ofl << setprecision(8) << dep_vars(ii++) << endl;
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("adult_rbio(")
        + adstring("recent") + adstring(")") << endl;
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
        ofl << adstring("adult_rbio(")
          + adstring("recent") + adstring("_species_") + str(is) + adstring(")") << endl;
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
          ofl << adstring("adult_rbio(")
            + adstring("recent") + adstring("_sex_") + str(is) + adstring(")") << endl;
        }
      }
    }


    // Abbreviated calculation based on parest_flags(37) NMD_24july2023
    if (fsh.parest_flags(37)!=1)
    {
//  - start of abbreviation
    

      dep_vars(ii) << fsh.sv(10);
      ofl << dep_vars(ii++) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << "sv(10)" << endl;
 
      dep_vars(ii) << fsh.sv(11);
      ofl << dep_vars(ii++) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << "sv(11)" << endl;
 
      dep_vars(ii) << fsh.sv(12);
      ofl << dep_vars(ii++) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << "sv(12)" << endl;

      if (fsh.age_flags(120))
      {
        do_average_exploitation(fsh,ii,num_grads,dep_vars,ofl);
      }
 
      dep_vars(ii) << log(rbio(ny)/max(rbio));
      ofl << dep_vars(ii++) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << "ln_rbio_maxr" << endl;
 
      dep_vars(ii) << log(rbio(ny)/rbio(1));
      ofl << dep_vars(ii++) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("log(rbio(") + str(ny) + adstring(")")
          << adstring("/rbio(") + str(1) + adstring("))") << endl;
 
      dvariable mrb=mean(rbio);
      dep_vars(ii) << log(rbio(ny)/mrb);
      ofl << dep_vars(ii++) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << "ln_ln_rbio_mean" << endl;
 
      // write index of first log bio in temporary file
      tmpout << ii << "  " ;
 
      for (i=1;i<=ny;i++)
      {
        dep_vars(ii) << log(rbio(i)/mrb);
        ofl << dep_vars(ii++) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("ln_rel_bio(") + str(i) + adstring(")") << endl;
      }

      for (i=1;i<=ny;i++)
      {
        dep_vars(ii) << log(rbio(i));
        ofl << dep_vars(ii++) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("ln_abs_bio(") + str(i) + adstring(")") << endl;
      }
      {
        dvar_vector rbio=sbcalc(fsh);
        if (fsh.age_flags(95)>0)
          dep_vars(ii) << log(rbio(ny)/mean(rbio(1,fsh.age_flags(95))));
        else
          dep_vars(ii) << log(rbio(ny)/rbio(1));
        ofl << dep_vars(ii++) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("relative_spawning_biomass") << endl;
      }
      dvar_vector rec(1,ny);
      rec.initialize();
      for (ir=1;ir<=fsh.num_regions;ir++)
      {
        for (int i=1;i<=ny;i++)
        {
          rec(i)+=exp(fsh.N(ir,i,1));
        }
      }
      dvariable avgrecr=log(mean(rec));
      rec=log(1.e-10+rec);
      for (i=1;i<=ny;i++)
      {
        dep_vars(ii) << rec(i)-avgrecr;
        ofl << dep_vars(ii++) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("ln_relative_recr(") + str(i) + adstring(")") << endl;
      }
      //cout<< "XXXXXXX newgradc.cpp:207: "<<rec<<endl;
      for (i=1;i<=ny;i++)
      {
        dep_vars(ii) << rec(i);
        ofl << dep_vars(ii++) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("ln_abs_recr(") + str(i) + adstring(")") << endl;
      }


      for (i=1;i<=fsh.num_fisheries;i++)
      {
        int nt=fsh.num_fish_times(i);
        int rr=fsh.realization_region(i,nt);
        int rp=fsh.realization_period(i,nt);
        int ri=fsh.realization_incident(i,nt);
        int rp1=fsh.realization_period(i,1);
        int ri1=fsh.realization_incident(i,1);

        dep_vars(ii) << 
        fsh.catchability(rr,fsh.realization_period(i,nt),
          fsh.realization_incident(i,nt)) 
          - (fsh.catchability(rr,fsh.realization_period(i,1),
          fsh.realization_incident(i,1)));

        ofl << dep_vars(ii++) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("lcatch-diff(") + str(i) + adstring(")") << endl;

        if (fsh.fish_flags(i,27))
        {
          const MY_DOUBLE_TYPE tpi=2.e0*3.14159e0;
          dvariable tmp=0.0;
          dvariable tmp1=0.0;;
          tmp=fsh.fish_pars(1,i)*sin(tpi*(fsh.true_month(rr,rp)/12.-fsh.fish_pars(2,i)));
          tmp1=fsh.fish_pars(1,i)*sin(tpi*(fsh.true_month(rr,rp1)/12.-fsh.fish_pars(2,i)));
          dep_vars(ii) << 
          fsh.catchability(rr,fsh.realization_period(i,nt),
            fsh.realization_incident(i,nt)) - tmp
            - (fsh.catchability(rr,fsh.realization_period(i,1),
            fsh.realization_incident(i,1))-tmp1);

          ofl << dep_vars(ii++) << endl; 
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("unseasoned_lcatch-diff(") + str(i) + adstring(")") << endl;
        } 
      }
    }
//  - end of abbreviation
    

    if (fsh.age_flags(145))
    {
      int mmin;
      int mmax;
      if (fsh.parest_flags(37)!=1)
      {
//  - start of abbreviation


//       int mmin=fsh.predicted_yield_bh.indexmin();
//       int mmax=fsh.predicted_yield_bh.indexmax();
        mmin=fsh.predicted_yield_bh.indexmin();
        mmax=fsh.predicted_yield_bh.indexmax();
        for (i=mmin;i<=mmax;i++)
        {
          MY_DOUBLE_TYPE x=fsh.predicted_yield_bh_x(i);
          dep_vars(ii) << fsh.predicted_yield_bh(i);
          ofl << dep_vars(ii++) << endl; 
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("pred_yield_for_F_mult ") 
             + str(x,0,1)  << endl;
        }
//    JH 21/03/02  Adds equilib. spawning biomass to dep_vars
        mmin=fsh.predicted_eqbio_bh.indexmin();
        mmax=fsh.predicted_eqbio_bh.indexmax();
        for (i=mmin;i<=mmax;i++)
        {
          MY_DOUBLE_TYPE x=fsh.predicted_yield_bh_x(i);
          dep_vars(ii) << fsh.predicted_eqbio_bh(i);
          ofl << dep_vars(ii++) << endl; 
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("pred_equilib_SB_for_F_mult ") 
             + str(x,0,1)  << endl;
        }
//    -------------------------------------------------------
//    JH 21/03/02  Adds equilib. total biomass to dep_vars
        mmin=fsh.predicted_eqtotbio_bh.indexmin();
        mmax=fsh.predicted_eqtotbio_bh.indexmax();
        for (i=mmin;i<=mmax;i++)
        {
          MY_DOUBLE_TYPE x=fsh.predicted_yield_bh_x(i);
          dep_vars(ii) << fsh.predicted_eqtotbio_bh(i);
          ofl << dep_vars(ii++) << endl; 
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("pred_equilib_totbio_for_F_mult ") 
             + str(x,0,1)  << endl;
        }

//    -------------------------------------------------------
//    JH 21/03/02  Adds AdultB/AdultB(msy) to dep_vars
        for (i=1;i<=ny;i++)
        {
          dep_vars(ii) << fsh.ab_ratio(i);
          ofl << dep_vars(ii++) << endl; 
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("AdultB/AdultBmsy(") + str(i) + adstring(")") << endl;
        }
// Add dep_vars on the log-scale //NMD_31may2022
        for (i=1;i<=ny;i++)
        {
          dep_vars(ii) << log(fsh.ab_ratio(i));
          ofl << dep_vars(ii++) << endl; 
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("ln_AdultB/AdultBmsy(") + str(i) + adstring(")") << endl;
        }

//    -------------------------------------------------------
//    JH 21/03/02  Adds TotalB/TotalB(msy) to dep_vars
        for (i=1;i<=ny;i++)
        {
          dep_vars(ii) << fsh.tb_ratio(i);
          ofl << dep_vars(ii++) << endl; 
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("TotalB/TotalBmsy(") + str(i) + adstring(")") << endl;
        }
// Add dep_vars on the log-scale //NMD_31may2022
        for (i=1;i<=ny;i++)
        {
          dep_vars(ii) << log(fsh.tb_ratio(i));
          ofl << dep_vars(ii++) << endl; 
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("ln_TotalB/TotalBmsy(") + str(i) + adstring(")") << endl;
        }

//    -------------------------------------------------------
//    JH 27/03/02  Adds F/F(msy) to dep_vars
        for (i=1;i<=ny;i++)
        {
          dep_vars(ii) << fsh.F_ratio(i);
          ofl << dep_vars(ii++) << endl; 
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("F/Fmsy(") + str(i) + adstring(")") << endl;
        }
// Add dep_vars on the log-scale //NMD_31may2022
        for (i=1;i<=ny;i++)
        {
          dep_vars(ii) << log(fsh.F_ratio(i));
          ofl << dep_vars(ii++) << endl; 
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("ln_F/Fmsy(") + str(i) + adstring(")") << endl;
        }
 //    -------------------------------------------------------
      }
//  - end of abbreviation

      dvariable avg_F_ratio=0.0;
      dvariable avg_SB_ratio=0.0;
      int offset=af57;
      for (i=yr1st;i<=yrlast;i++)
      {
        avg_F_ratio+=fsh.F_ratio(i-offset);
        avg_SB_ratio+=fsh.ab_ratio(i);
      }
      avg_F_ratio/=nfa_calc;
      avg_SB_ratio/=nfa_calc;

      dep_vars(ii) << avg_F_ratio;
      ofl << setprecision(8) << dep_vars(ii++) << endl;
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("average_F/Fmsy(")
        + adstring("recent") + adstring(")") << endl;

      dep_vars(ii) << avg_SB_ratio;
      ofl << setprecision(8) << dep_vars(ii++) << endl;
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("average_SB/SBmsy(")
        + adstring("recent") + adstring(")") << endl;

      dep_vars(ii) << fsh.MMSY;
      ofl << setprecision(8) << dep_vars(ii++) << endl;
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("MSY") << endl;

      if (fsh.parest_flags(37)!=1)
      {
//  - start of abbreviation
        mmin=fsh.predicted_recruitment_bh.indexmin();
        mmax=fsh.predicted_recruitment_bh.indexmax();
        for (i=mmin;i<=mmax;i++)
        {
          MY_DOUBLE_TYPE x=fsh.predicted_recruitment_bh_x(i);
          dep_vars(ii) << fsh.predicted_recruitment_bh(i);
          ofl << dep_vars(ii++) << endl; 
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("pred_rec_for_bio_level ") 
             + str(x,0,1) << endl;
        }
      }
//  - end of abbreviation
       
    }

    if (fsh.parest_flags(37)!=1)
    {
//  - start of abbreviation

      if (fsh.age_flags(150))
      {
        dvariable r=fsh.sv(22)+0.2;
        dvariable eta=fsh.sv(23)+1.0;
// John H. 26/10/2001
//      dvariable m=fsh.sv(24);
        dvariable m=fsh.sv(24)+2.0;
        dvariable k=fsh.biomass(1)*eta;
        dvariable Bmsy;
        dvariable Msy;
        get_msy_pt(Bmsy,Msy,k,m,r);

        dep_vars(ii) << r;
        ofl << dep_vars(ii++) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("pella_t r") << endl;

        dep_vars(ii) << eta;
        ofl << dep_vars(ii++) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("pella_t eta") << endl;

        dep_vars(ii) << k;
        ofl << dep_vars(ii++) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("pella_t k") << endl;

        dep_vars(ii) << m;
        ofl << dep_vars(ii++) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("pella_t m") << endl;

        dep_vars(ii) << Bmsy;
        ofl << dep_vars(ii++) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("pella_t Bmsy") << endl;

        dep_vars(ii) << Msy;
        ofl << dep_vars(ii++) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("pella_t Msy") << endl;

        int mmin=fsh.predicted_yield_pt.indexmin();
        int mmax=fsh.predicted_yield_pt.indexmax();
        for (i=mmin;i<=mmax;i++)
        {
          MY_DOUBLE_TYPE x=fsh.predicted_yield_pt_x(i);
          dep_vars(ii) << fsh.predicted_yield_pt(i);
          ofl << dep_vars(ii++) << endl; 
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("pred_yield_for_effort_level ") 
             + str(x,0,1) << endl;
        }
      }

      for (i=1;i<=fsh.nage;i++)
      {
        dvariable tmp=setup_alpha(i,fsh.sv,fsh.nage,fsh.age_flags);
        dep_vars(ii) <<  tmp;
        ofl << dep_vars(ii++) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("diffusion_rate_for_age_class(") + str(i) + adstring(")") 
            << endl;
      }

      if (fsh.parest_flags(157)>0)
      {
        dvar_vector rel_strength(1,ny);
        rel_strength.initialize();
        for (int i=1;i<=ny;i++)
        {
          for (int ir=1;ir<=fsh.num_regions;ir++)
          {
            rel_strength(i)+=exp(fsh.N(ir,i,1));
          }
        }
        rel_strength=rel_strength-mean(rel_strength);
        rel_strength/=sqrt((.001+norm2(rel_strength))/rel_strength.size());
        dvariable agedif=get_dd_agedif(fsh.sv,rel_strength);
        dep_vars(ii) << agedif;
        ofl << dep_vars(ii++) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << "dd_age_dif" << endl;
        cout << "dd_age_dif" << endl;
      }

     // new stuff for projections Df 16/02/2010
    /* 
     for (int ir=1;ir<=fsh.num_regions;ir++) 
     {
       int ip=fsh.num_real_fish_periods(ir);
       for (int j=1;j<=fsh.nage;j++)
       {
         dep_vars(ii) << exp(fsh.num_fish(ir,ip,j));
         ofl << dep_vars(ii++) << endl; 
         depvars_bounds_check(ii,num_grads);
         ofl << adstring("exp_num_fish(")  
            + str(ir) +adstring(",")
            + str(ip) +adstring(",")
            + str(j)  + adstring(")") << endl;
       }
     }
    */
    
   
    // save values of dependent variables
//  ***************************************************************
// 
//  ***************************************************************
      if (fsh.age_flags(73)>0)
      {
        int mmax=fsh.nage;
        if (fsh.age_flags(81)>0)
          mmax=fsh.nage-fsh.age_flags(81)+1;
        for (int j=1;j<=mmax;j++)
        {
          dep_vars(ii) << fsh.age_pars(2,j);
          ofl << dep_vars(ii++) << endl; 
          depvars_bounds_check(ii,num_grads);
          ofl << adstring("age_pars(")  
            + str(2) +adstring(",")
            + str(j)  + adstring(")") << endl;
        }
      }

      for (int i=1;i<=fsh.nage;i++)
      {
        dep_vars(ii) << fsh.nat_mort(1,i);
        ofl << dep_vars(ii++) << endl;
        depvars_bounds_check(ii,num_grads);
        ofl << "nat_mort(" << i << ")" << endl;
      }
    }
//  - end of abbreviation

    
    int ndep=ii-1;
    mfglobals.dep_vars_values.allocate(1,ndep);
    int ij;
    for (ij=1;ij<=ndep;ij++)
    {
      mfglobals.dep_vars_values(ij) = value(dep_vars(ij));
    }
    return ndep;
  }

void do_average_exploitation(dvar_len_fish_stock_history& fsh,int& ii,
  int& num_grads,dvar_vector& dep_vars,ofstream& ofl)
{
  dvar_matrix totnum(1,fsh.nyears,1,fsh.nage);
  dvar_matrix ssum(1,fsh.nyears,1,fsh.nage);
  dvar_matrix totcatch(1,fsh.nyears,1,fsh.nage);
  totnum.initialize();
  totcatch.initialize();
  
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++) 
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) 
      {                                 
        totcatch(fsh.year(ir,ip))+=exp(fsh.catch(ir,ip,fi));
      }
    }
  }
  int i;
  for (i=1;i<=fsh.nyears;i++)
  {
    for (int ir=1;ir<=fsh.num_regions;ir++)
    {
      totnum(i)+=exp(fsh.N(ir,i));
    }
  }
  dvar_vector avg_f1(1,fsh.nyears);
  dvar_vector avg_f2(1,fsh.nyears);
  dvar_vector num_f1(1,fsh.nyears);
  dvar_vector num_f2(1,fsh.nyears);
  dvar_matrix tf(1,fsh.nyears,1,fsh.nage);
  dvar_matrix tn(1,fsh.nyears,1,fsh.nage);
  avg_f1.initialize();
  tf.initialize();
  tn.initialize();
  avg_f2.initialize();


  for (int iy=1;iy<=fsh.nyears;iy++) 
  {
    int j;
    for (j=2;j<=5;j++)
    {
      avg_f1(iy)+=totcatch(iy,j);
      num_f1(iy)+=totnum(iy,j);
    }
    for (j=6;j<=fsh.nage;j++)
    {
      avg_f2(iy)+=totcatch(iy,j);
      num_f2(iy)+=totnum(iy,j);
    }
  }

  avg_f1=log(1.e-10+elem_div(avg_f1,num_f1));
  avg_f2=log(1.e-10+elem_div(avg_f2,num_f2));

  for (i=1;i<=fsh.nyears;i++)
  {
    dep_vars(ii) << avg_f1(i);
    ofl << dep_vars(ii++) << endl; 
    depvars_bounds_check(ii,num_grads);
    ofl << adstring("ln_avg_exp_2-5(") + str(i) + adstring(")") << endl;
  }
  for (i=1;i<=fsh.nyears;i++)
  {
    dep_vars(ii) << avg_f2(i);
    ofl << dep_vars(ii++) << endl; 
    depvars_bounds_check(ii,num_grads);
    ofl << adstring("ln_avg_exp_6-nage(") + str(i) + adstring(")") << endl;
  }
}

dvariable get_dd_agedif(dvar_vector& sv,dvar_vector& rs)
{
  dvariable rrmin=min(rs);
  dvariable rrmax=max(rs);
  dvariable tmp;
  tmp=sv(3)*rrmin;
  tmp=1./(1.+exp(-tmp));
  dvariable minage=1.9*(tmp-.5L);
  tmp=sv(3)*rrmax;
  tmp=1./(1.+exp(-tmp));
  dvariable maxage=1.9*(tmp-.5L);
  return maxage-minage;
}

/*
// void do_average_fish_mort(dvar_len_fish_stock_history& fsh,int& ii,
//   dvar_vector& dep_vars,ofstream& ofl)
// {
//   dvar3_array ssum(1,fsh.num_regions,1,fsh.nyears,1,fsh.nage);
//   dvar_matrix tf(1,fsh.nyears,1,fsh.nage);
//   dvar_matrix tn(1,fsh.nyears,1,fsh.nage);
//   dvar_vector avg_f1(1,fsh.nyears);
//   dvar_vector avg_f2(1,fsh.nyears);
//   dvar_vector num_f1(1,fsh.nyears);
//   dvar_vector num_f2(1,fsh.nyears);
//   ssum.initialize();
//   avg_f1.initialize();
//   tf.initialize();
//   tn.initialize();
//   avg_f2.initialize();
//  
//   for (int ir=1;ir<=fsh.num_regions;ir++)
//   {
//     for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++) 
//     {
//       for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) 
//       {                                 
//         ssum(ir,fsh.year(ir,ip))+=exp(fsh.fish_mort(ir,ip,fi));
//       }
//     }
//   }
//   for (int iy=1;iy<=fsh.nyears;iy++) 
//   {
//     for (int ir=1;ir<=fsh.num_regions;ir++)
//     {
//       tf(iy)+=elem_prod(ssum(ir,iy),exp(fsh.N(ir,iy)));
//       tn(iy)+=exp(fsh.N(ir,iy));
//     }
//   }
//
//   for (iy=1;iy<=fsh.nyears;iy++) 
//   {
//     for (int j=2;j<=5;j++)
//     {
//       avg_f1(iy)+=tf(iy,j);
//       num_f1(iy)+=tn(iy,j);
//     }
//     for (j=6;j<=9;j++)
//     {
//       avg_f2(iy)+=tf(iy,j);
//       num_f2(iy)+=tn(iy,j);
//     }
//   }
//
//   avg_f1=log(1.e-10+elem_div(avg_f1,num_f1));
//   avg_f2=log(1.e-10+elem_div(avg_f2,num_f2));
//
//    for (int i=1;i<=fsh.nyears;i++)
//   // for (int i=1;i<=2;i++)
//    {
//      dep_vars(ii) << avg_f1(i);
//      ofl << dep_vars(ii++) << endl; 
//      depvars_bounds_check(ii,num_grads);
//      ofl << adstring("ln_avg_f2-5(") + str(i) + adstring(")") << endl;
//    }
//    //for (i=1;i<=2;i++)
//    for (i=1;i<=fsh.nyears;i++)
//    {
//      dep_vars(ii) << avg_f2(i);
//      ofl << dep_vars(ii++) << endl; 
//      depvars_bounds_check(ii,num_grads);
//      ofl << adstring("ln_avg_f6-9(") + str(i) + adstring(")") << endl;
//    }
//  }
*/
