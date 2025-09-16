/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#define HOME_VERSION
#include "all.hpp"

  
dvariable dvar_fish_stock_history::kludged_initial_survival_equilibrium_penalty
  (ivector * pq_flag)
{
  int nya=age_flags(95); // num_years_for_average
  if (nya==0)
  {
    cerr << "Error age_flags(95) must be >0 " << endl;
    ad_exit(1);
  }
  ivector  mps = get_equilibrium_movements_per_season();
//  dvar4_array surv=get_initial_equilibrium_survival(mps,nya,pq_flag);
  dvar4_array local_surv=get_initial_equilibrium_survival(mps,nya,pq_flag);
  MY_DOUBLE_TYPE pen=10.;
  if (parest_flags(375)!=0)
  {
    pen=parest_flags(375);
  }
//  int mmin=surv.indexmin();
//  int mmax=surv.indexmax();
  int mmin=local_surv.indexmin();
  int mmax=local_surv.indexmax();
  dvariable tmp=0.0;
  for (int i=mmin;i<=mmax;i++)
  {
//    tmp+=pen*norm2(surv(i)-kludged_surv(i));
    tmp+=pen*norm2(local_surv(i)-kludged_surv(i));
  }

  if (print_implicit_effort_flag && !af170q0)  //NMD_1june2020
  {
    ofstream ofs("kludged_surv.rpt");
    ofs << " Actual and kludged survival report" << endl;
    ofs << "  Actual/Kludged Season  Mvmnt_period  Region  Survival " << endl;
    for (int is=1;is<=age_flags(57);is++)
    {
      for (int im=1;im<=mps(is);im++)
      {
//	dmatrix tmp=trans(value(surv(is,im)));
//	dmatrix tmp=trans(value(local_surv(is,im)));
	dmatrix ptmp=trans(value(local_surv(is,im))); //NMD_9may2022
	dmatrix kltmp=trans(value(kludged_surv(is,im)));	
        for (int ir=1;ir<=num_regions;ir++)
        {
          ofs << "Actual  " << is << setw(4) << im << setw(4) <<  ir << setw(6) << ptmp(ir) << endl;
	  ofs << "Kludge  " <<is <<  setw(4) << im << setw(4) <<  ir << setw(6) << kltmp(ir) << endl;
	} 
      }
    }
  }

//  cout << print_implicit_effort_flag << endl;
//  cout << af170q0 << endl;
  cout << "TTT kludged surv penalty " << tmp << endl;
  return tmp;
}

dvar4_array  dvar_fish_stock_history::get_kludged_initial_equilibrium_survival
  (ivector & mps,ivector * pq_flag)
{
  // Now we need to determine the total morality in each
  // between movement period
  int sd=parest_flags(374);
//  dvar4_array surv(1,age_flags(57),1,mps,1,nage,1,num_regions);  
  dvar4_array local_surv(1,age_flags(57),1,mps,1,nage,1,num_regions);  

  dvar4_array tmort(1,age_flags(57),
    1,mps,1,num_regions,1,nage);  
  //dvar4_array surv(1,age_flags(57),
  //  1,mps,1,nage,1,num_regions);  
  dvar3_array kludged_selmean(1,age_flags(57),
    1,mps,1,num_regions);  
  //dvar4_array kludged_equilib_coffs(1,age_flags(57),
  //  1,mps,1,num_regions,1,sd);  
  if (!allocated(kludged_equilib_coffs))
  {
    //kludged_equilib_coffs.allocate(1,age_flags(57),1,mps,
    //  1,num_regions,1,sd);  
  }
  tmort.initialize();  
  dvector x(1,sd);
#if !defined(NO_MY_DOUBLE_TYPE)
  x.fill_seqadd(0.0,1.0/(sd-1.0L));
#else
  x.fill_seqadd(0.0,1.0/(sd-1.0));
#endif
  

  kludged_selmean_square=0.0;
  for (int is=1;is<=age_flags(57);is++)
  {
    for (int im=1;im<=mps(is);im++)
    {
      for (int ir=1;ir<=num_regions;ir++)
      {
        dvar_vector tsel(1,sd);
        tsel=kludged_equilib_coffs(is,im,ir)(1,sd);
        kludged_selmean(is,im,ir)=log(mean(exp(tsel)));
        kludged_selmean_square+=square(kludged_selmean(is,im,ir));
        tsel-=kludged_selmean(is,im,ir);
        tsel+=kludged_equilib_level_coffs(ir);
        vcubic_spline_function csf(x,tsel);
        //int np=nage;
        dvector xx(1,nage);
#if !defined(NO_MY_DOUBLE_TYPE)
        xx.fill_seqadd(0,1.0/(nage-1.0L));
#else
        xx.fill_seqadd(0,1.0/(nage-1.0));
#endif
        //cout << " !!! Need to fix this -- kludged_equilib_surv.cpp 112 " << endl;
        tmort(is,im,ir)=mfexp(csf(xx));
        //for (int j=1;j<=nage;j++)
        //{
        //  tmort(is,im,ir,j)=0.4*j;
        //}
      }
    }
  }
  for (int is=1;is<=age_flags(57);is++)
  {
    for (int pi=1;pi<=mps(is);pi++)
    {
      if (pq_flag)
      {
//        surv(is,pi)=mfexp(-trans(tmort(is,pi))-0.0001);
        local_surv(is,pi)=mfexp(-trans(tmort(is,pi))-0.0001);
      }
      else
      {
//        surv(is,pi)=mfexp(-trans(tmort(is,pi))-0.0001);
        local_surv(is,pi)=mfexp(-trans(tmort(is,pi))-0.0001);
      }
    }
  }
//  return surv;
  return local_surv;
}

