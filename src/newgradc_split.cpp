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

  dvariable get_dd_agedif(dvar_vector& sv,dvar_vector& rs);
  dvar_vector sbcalc(dvar_len_fish_stock_history& fsh);
  void do_average_exploitation(dvar_len_fish_stock_history& fsh,int& ii,
    int& num_grads,dvar_vector& dep_vars,ofstream& ofl);

  void depvars_bounds_check(int i,int n);
  void do_average_exploitation_split(dvar_len_fish_stock_history& fsh,int& ii,
    int& num_grads,dvar_vector& dep_vars,ofstream& ofl,const int & depcount);

  int dep_gradients_calc2_split(dvar_len_fish_stock_history& fsh)
  {
    int depcount=0;
     int pf229=fsh.parest_flags(229); 
     int pf230=fsh.parest_flags(230); 
     int offset=pf229-1;
    //HERE
    ofstream tmpout("tmp.out");
    ofstream ofl("deplabel.tmp");
    //HERE
    int num_grads=50000;
    int ii=0;
    //HERE
    dvar_vector dep_vars(1,num_grads);

    // do relative biomass gradients
    dvar_vector rbio=vbiocalc(fsh);
    int ny=0;
    if (fsh.do_fishery_projections_flag==0)
      ny=fsh.num_real_years;
    else
      ny=fsh.nyears;
    int i;
    for (i=1;i<=ny;i++)  
    {
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
      dep_vars(ii-offset) << log(rbio(i));
      ofl << dep_vars(ii-offset) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("rbio(")  
         + str(i) + adstring(")") << endl;
        }
    }
    dvar_vector adult_rbio=adult_vbiocalc(fsh);
    for (i=1;i<=ny;i++)  
    {
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
      dep_vars(ii-offset) << log(adult_rbio(i));
      ofl << dep_vars(ii-offset) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("adult_rbio(")  
         + str(i) + adstring(")") << endl;
        }
    }

    dvar_matrix reg_rbio=unnormalized_reg_vbiocalc(fsh);
    int ir;
    for (ir=1;ir<=fsh.num_regions;ir++)
    {
      for (int i=1;i<=ny;i++)  
      {
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
        dep_vars(ii-offset) << log(reg_rbio(ir,i));
        ofl << dep_vars(ii-offset) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("ln_reg_bio(") +str(ir) +adstring(",") 
           + str(i) + adstring(")") << endl;
        }
      }
    }

    dvar_matrix adult_reg_rbio=unnormalized_adult_reg_vbiocalc(fsh);
    for (ir=1;ir<=fsh.num_regions;ir++)
    {
      for (int i=1;i<=ny;i++)  
      {
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
        dep_vars(ii-offset) << log(adult_reg_rbio(ir,i));
        ofl << dep_vars(ii-offset) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("ln_adult_reg_bio(") +str(ir) +adstring(",") 
           + str(i) + adstring(")") << endl;
        }
      }
    }

    {
      dvar_vector msum(1,fsh.nage);
      msum.initialize();
      for (int i=1;i<=ny;i++)  
      {
        for (int j=1;j<=fsh.nage;j++)
        {
          msum(j)+=fsh.nat_mort(i,j);
        }
      }
      for (int j=1;j<=fsh.nage;j++)
      {
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
        dep_vars(ii-offset) << msum(j);
        ofl << dep_vars(ii-offset) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << "msum(" << str(i) << ")" << endl;
     }
      }
    }

     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
    dep_vars(ii-offset) << fsh.nat_mort(1,1);
    ofl << dep_vars(ii-offset) << endl; 
    depvars_bounds_check(ii,num_grads);
    ofl << "nat_mort" << endl;
     }
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
    dep_vars(ii-offset) << fsh.sv(10);
    ofl << dep_vars(ii-offset) << endl; 
    depvars_bounds_check(ii,num_grads);
    ofl << "sv(10)" << endl;
     }

     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
    dep_vars(ii-offset) << fsh.sv(11);
    ofl << dep_vars(ii-offset) << endl; 
    depvars_bounds_check(ii,num_grads);
    ofl << "sv(11)" << endl;
     }

     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
    dep_vars(ii-offset) << fsh.sv(12);
    ofl << dep_vars(ii-offset) << endl; 
    depvars_bounds_check(ii,num_grads);
    ofl << "sv(12)" << endl;
     }

    if (fsh.age_flags(120))
    {
      do_average_exploitation_split(fsh,ii,num_grads,dep_vars,ofl,depcount);
    }

     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
    dep_vars(ii-offset) << log(rbio(ny)/max(rbio));
    ofl << dep_vars(ii-offset) << endl; 
    depvars_bounds_check(ii,num_grads);
    ofl << "ln_rbio_maxr" << endl;
     }

     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
    dep_vars(ii-offset) << log(rbio(ny)/rbio(1));
    ofl << dep_vars(ii-offset) << endl; 
    depvars_bounds_check(ii,num_grads);
    ofl << adstring("log(rbio(") + str(ny) + adstring(")")
        << adstring("/rbio(") + str(1) + adstring("))") << endl;
        }

    //HERE
    // if (++ii>= pf229 && ii <= pf230) 
    // { 
    //   depcount++;
    //dep_vars(ii-offset) << log(-log(0.99*rbio(fsh.nyears)/max(rbio)));
   // ofl << dep_vars(ii-offset) << endl; 
   // ofl << "ln_ln_rbio_maxr" << endl;
    // }

    //HERE
    dvariable mrb=mean(rbio);
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
    dep_vars(ii-offset) << log(rbio(ny)/mrb);
    ofl << dep_vars(ii-offset) << endl; 
    depvars_bounds_check(ii,num_grads);
    ofl << "ln_ln_rbio_mean" << endl;
     }

    // write index of first log bio in temporary file
    //HERE
    tmpout << ii << "  " ;

    for (i=1;i<=ny;i++)
    {
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
      dep_vars(ii-offset) << log(rbio(i)/mrb);
      ofl << dep_vars(ii-offset) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("ln_rel_bio(") + str(i) + adstring(")") << endl;
     }
    }

    for (i=1;i<=ny;i++)
    {
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
      dep_vars(ii-offset) << log(rbio(i));
      ofl << dep_vars(ii-offset) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("ln_abs_bio(") + str(i) + adstring(")") << endl;
     }
    }
    {
      dvar_vector rbio=sbcalc(fsh);
     if (++ii>= pf229 && ii <= pf230) 
     { 
      depcount++;
      if (fsh.age_flags(95)>0)
        dep_vars(ii-offset) << log(rbio(ny)/mean(rbio(1,fsh.age_flags(95))));
      else
        dep_vars(ii-offset) << log(rbio(ny)/rbio(1));
      ofl << dep_vars(ii-offset) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("relative_spawning_biomass") << endl;
     }
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
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
        dep_vars(ii-offset) << rec(i)-avgrecr;
        ofl << dep_vars(ii-offset) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("ln_relative_recr(") + str(i) + adstring(")") << endl;
     }
      }

      for (i=1;i<=ny;i++)
      {
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
        dep_vars(ii-offset) << rec(i);
        ofl << dep_vars(ii-offset) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("ln_abs_recr(") + str(i) + adstring(")") << endl;
     }
      }


    for (i=1;i<=fsh.num_fisheries;i++)
    {
      int nt=fsh.num_fish_times(i);
      int rr=fsh.realization_region(i,nt);
      int rp=fsh.realization_period(i,nt);
      int ri=fsh.realization_incident(i,nt);
      int rp1=fsh.realization_period(i,1);
      int ri1=fsh.realization_incident(i,1);

     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
      dep_vars(ii-offset) << 
      fsh.catchability(rr,fsh.realization_period(i,nt),
        fsh.realization_incident(i,nt)) 
        - (fsh.catchability(rr,fsh.realization_period(i,1),
        fsh.realization_incident(i,1)));

      ofl << dep_vars(ii-offset) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("lcatch-diff(") + str(i) + adstring(")") << endl;
     }

      if (fsh.fish_flags(i,27))
      {
        const MY_DOUBLE_TYPE tpi=2.e0*3.14159e0;
        dvariable tmp=0.0;
        dvariable tmp1=0.0;;
        tmp=fsh.fish_pars(1,i)*sin(tpi*(fsh.true_month(rr,rp)/12.-fsh.fish_pars(2,i)));
        tmp1=fsh.fish_pars(1,i)*sin(tpi*(fsh.true_month(rr,rp1)/12.-fsh.fish_pars(2,i)));
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
        dep_vars(ii-offset) << 
        fsh.catchability(rr,fsh.realization_period(i,nt),
          fsh.realization_incident(i,nt)) - tmp
          - (fsh.catchability(rr,fsh.realization_period(i,1),
          fsh.realization_incident(i,1))-tmp1);

        ofl << dep_vars(ii-offset) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("unseasoned_lcatch-diff(") + str(i) + adstring(")") << endl;
     }
      } 
    }

    if (fsh.age_flags(145))
    {
      int mmin=fsh.predicted_yield_bh.indexmin();
      int mmax=fsh.predicted_yield_bh.indexmax();
      for (i=mmin;i<=mmax;i++)
      {
        MY_DOUBLE_TYPE x=fsh.predicted_yield_bh_x(i);
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
        dep_vars(ii-offset) << fsh.predicted_yield_bh(i);
        ofl << dep_vars(ii-offset) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("pred_yield_for_F_mult ") 
           + str(x,0,1)  << endl;
        }
      }
//    JH 21/03/02  Adds equilib. spawning biomass to dep_vars
      mmin=fsh.predicted_eqbio_bh.indexmin();
      mmax=fsh.predicted_eqbio_bh.indexmax();
      for (i=mmin;i<=mmax;i++)
      {
        MY_DOUBLE_TYPE x=fsh.predicted_yield_bh_x(i);
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
        dep_vars(ii-offset) << fsh.predicted_eqbio_bh(i);
        ofl << dep_vars(ii-offset) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("pred_equilib_SB_for_F_mult ") 
           + str(x,0,1)  << endl;
        }
      }
//    -------------------------------------------------------
//    JH 21/03/02  Adds equilib. total biomass to dep_vars
      mmin=fsh.predicted_eqtotbio_bh.indexmin();
      mmax=fsh.predicted_eqtotbio_bh.indexmax();
      for (i=mmin;i<=mmax;i++)
      {
        MY_DOUBLE_TYPE x=fsh.predicted_yield_bh_x(i);
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
        dep_vars(ii-offset) << fsh.predicted_eqtotbio_bh(i);
        ofl << dep_vars(ii-offset) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("pred_equilib_totbio_for_F_mult ") 
           + str(x,0,1)  << endl;
        }
      }
//    -------------------------------------------------------
//    JH 21/03/02  Adds AdultB/AdultB(msy) to dep_vars
      for (i=1;i<=ny;i++)
      {
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
        dep_vars(ii-offset) << fsh.ab_ratio(i);
        ofl << dep_vars(ii-offset) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("AdultB/AdultBmsy(") + str(i) + adstring(")") << endl;
     }
      }
//    -------------------------------------------------------
//    JH 21/03/02  Adds TotalB/TotalB(msy) to dep_vars
      for (i=1;i<=ny;i++)
      {
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
        dep_vars(ii-offset) << fsh.tb_ratio(i);
        ofl << dep_vars(ii-offset) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("TotalB/TotalBmsy(") + str(i) + adstring(")") << endl;
     }
      }
//    -------------------------------------------------------
//    JH 27/03/02  Adds F/F(msy) to dep_vars
      for (i=1;i<=ny;i++)
      {
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
        dep_vars(ii-offset) << fsh.F_ratio(i);
        ofl << dep_vars(ii-offset) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("F/Fmsy(") + str(i) + adstring(")") << endl;
     }
      }
//    -------------------------------------------------------

      mmin=fsh.predicted_recruitment_bh.indexmin();
      mmax=fsh.predicted_recruitment_bh.indexmax();
      for (i=mmin;i<=mmax;i++)
      {
        MY_DOUBLE_TYPE x=fsh.predicted_recruitment_bh_x(i);
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
        dep_vars(ii-offset) << fsh.predicted_recruitment_bh(i);
        ofl << dep_vars(ii-offset) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("pred_rec_for_bio_level ") 
           + str(x,0,1) << endl;
        }
      }
    }

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

     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
      dep_vars(ii-offset) << r;
      ofl << dep_vars(ii-offset) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("pella_t r") << endl;
     }

     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
      dep_vars(ii-offset) << eta;
      ofl << dep_vars(ii-offset) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("pella_t eta") << endl;
     }

     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
      dep_vars(ii-offset) << k;
      ofl << dep_vars(ii-offset) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("pella_t k") << endl;
     }

     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
      dep_vars(ii-offset) << m;
      ofl << dep_vars(ii-offset) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("pella_t m") << endl;
     }

     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
      dep_vars(ii-offset) << Bmsy;
      ofl << dep_vars(ii-offset) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("pella_t Bmsy") << endl;
     }

     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
      dep_vars(ii-offset) << Msy;
      ofl << dep_vars(ii-offset) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("pella_t Msy") << endl;
     }

      int mmin=fsh.predicted_yield_pt.indexmin();
      int mmax=fsh.predicted_yield_pt.indexmax();
      for (i=mmin;i<=mmax;i++)
      {
        MY_DOUBLE_TYPE x=fsh.predicted_yield_pt_x(i);
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
        dep_vars(ii-offset) << fsh.predicted_yield_pt(i);
        ofl << dep_vars(ii-offset) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("pred_yield_for_effort_level ") 
           + str(x,0,1) << endl;
        }
      }
    }

    for (i=1;i<=fsh.nage;i++)
    {
      dvariable tmp=setup_alpha(i,fsh.sv,fsh.nage,fsh.age_flags);
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
      dep_vars(ii-offset) <<  tmp;
      ofl << dep_vars(ii-offset) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << adstring("diffusion_rate_for_age_class(") + str(i) + adstring(")") 
          << endl;
        }
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
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
      dep_vars(ii-offset) << agedif;
      ofl << dep_vars(ii-offset) << endl; 
      depvars_bounds_check(ii,num_grads);
      ofl << "dd_age_dif" << endl;
     }
      cout << "dd_age_dif" << endl;
    }
    // save values of dependent variables
//  ***************************************************************
// */
//  ***************************************************************
    if (pf230==0)
    {
      cout << "the number of dependent variables is " << ii-1 << endl;
      ad_exit(1);
    }
    int ndep=depcount;
    //HERE
    mfglobals.dep_vars_values.allocate(1,ndep);
    //HERE
    int ij;
    for (ij=1;ij<=ndep;ij++)
    {
      mfglobals.dep_vars_values(ij) = value(dep_vars(ij));
    }
    //HERE
    return ndep;
  }


void do_average_exploitation_split(dvar_len_fish_stock_history& fsh,int& ii,
  int& num_grads,dvar_vector& dep_vars,ofstream& ofl,const int & _depcount)
{
  int & depcount = (int &)(_depcount);
  int pf229=fsh.parest_flags(229); 
  int pf230=fsh.parest_flags(230); 
  int offset=pf229-1;
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
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
    dep_vars(ii-offset) << avg_f1(i);
    ofl << dep_vars(ii-offset) << endl; 
    depvars_bounds_check(ii,num_grads);
    ofl << adstring("ln_avg_exp_2-5(") + str(i) + adstring(")") << endl;
     }
  }
  for (i=1;i<=fsh.nyears;i++)
  {
     if (++ii>= pf229 && ii <= pf230) 
     { 
       depcount++;
    dep_vars(ii-offset) << avg_f2(i);
    ofl << dep_vars(ii-offset) << endl; 
    depvars_bounds_check(ii,num_grads);
    ofl << adstring("ln_avg_exp_6-nage(") + str(i) + adstring(")") << endl;
     }
  }
}
