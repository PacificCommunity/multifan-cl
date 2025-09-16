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

void dvar_fish_stock_history::xget_initial_age_structure_equilibrium(void)
{
  pmsd_error();
  if (age_flags(145)==0 && age_flags(171)==1)
  {
     cerr << "Can't have age_flags(171)=1 when SRR not active" << endl;
     ad_exit(1);
  }
  // get the initial age structure assuming equilibrium conditions
  // the inputs are the (log) equilibrium recruitment vector
  // and the (log) equilibrium survivial rate matrix
  dvar_matrix les=get_equilibrium_survival_rate();
  int ir;
  dvar_matrix tmpN(1,num_regions,1,nage+1);
  tmpN.initialize();
  // do equilbrium calculations for multi region situation
  dvar_vector tt(1,num_regions);
  tt.initialize();
  // calculate average initial recruitment
  if (af170q0==0)
  {
    for (ir=1;ir<=num_regions;ir++)
    {
      int ny=age_flags(95);
      for (int iy=2;iy<=ny;iy++)
      {
        tt(ir)+=mfexp(N(ir,iy,1));
      }
      tt(ir)/=(ny-1);
    }
  }
  else
  {
    for (ir=1;ir<=num_regions;ir++)
    {
      int ny=age_flags(95);
      for (int iy=2;iy<=ny;iy++)
      {
        tt(ir)+=mfexp(N_q0(ir,iy,1));
      }
      tt(ir)/=(ny-1);
    }
  }
 /*
  else if (age_flags(171)==0)
  {
    for (ir=1;ir<=num_regions;ir++)
    {
      int ny=age_flags(95);
      for (int iy=2;iy<=ny;iy++)
      {
        tt(ir)+=mfexp(N(ir,iy,1));
      }
      tt(ir)/=(ny-1);
    }
  }
  else if (age_flags(171)==1)
  {
    for (ir=1;ir<=num_regions;ir++)
    {
      int ny=age_flags(95);
      for (int iy=2;iy<=ny;iy++)
      {
        tt(ir)+=mfexp(N(ir,iy,1));
      }
      tt(ir)/=(ny-1);
    }
    tt/=sum(tt);
  }
  */

  for (ir=1;ir<=num_regions;ir++) tmpN(ir,1)=tt(ir);
  int nmp=mo.num_periods();
  dvar_matrix es=mfexp(les/double(nmp));
  for (int j=1;j<nage;j++)
  {
    for (int mp=1;mp<=nmp;mp++)
    {
      if (num_regions>1)
      {
        tt=Dad(move_map(mp),j)*tt;
      }
      for (ir=1;ir<=num_regions;ir++)
      {
        tt(ir)=es(ir,j)*tt(ir);
      }
    }
    for (ir=1;ir<=num_regions;ir++)
    {
      tmpN(ir,j+1)=tt(ir);
    }
  }
  // now for the last age class
  if (num_regions==1)
  {
    tmpN(1,nage)=tt(1)/(1.0-es(1,nage));
  }
  else
  {
    dvar_matrix meq(1,num_regions,1,num_regions);
    for (ir=1;ir<=num_regions;ir++) meq(ir,ir)=1.0;
    for (int mp=1;mp<=nmp;mp++)
    {
      meq=Dad(move_map(mp),nage)*meq;
      for (ir=1;ir<=num_regions;ir++)
      {
        for (int ir1=1;ir1<=num_regions;ir1++)
        {
          meq(ir,ir1)=es(ir1,nage)*Dad(move_map(mp),nage)(ir,ir1);
        }
      }
    }
    for (ir=1;ir<=num_regions;ir++)
    {
      meq(ir,ir)-=1;
    }
    tt=inv(meq)*tt;
    tt=-tt; // cause we used minus the matrix
    for (ir=1;ir<=num_regions;ir++)
    {
      tmpN(ir,nage)=tt(ir);
    }
  }
  tmpN=log(1.e-10+tmpN);

  if (af170q0==0)
  {
    for (ir=1;ir<=num_regions;ir++)
    {
      for (int j=1;j<=nage;j++)
      {
        N(ir,1,j)=tmpN(ir,j);
      }
    }
  }
  else
  {
    for (ir=1;ir<=num_regions;ir++)
    {
      for (int j=1;j<=nage;j++)
      {
        N_q0(ir,1,j)=tmpN(ir,j);
      }
    }
   /*
    for (ir=1;ir<=num_regions;ir++)
    {
      for (int j=1;j<=nage;j++)
      {
        N_q0(ir,1,j)=tmpN(ir,j);
        ssum+=exp(value(N_q0(ir,1,j)));
      }
    }
    dvariable loglam=0;
    dvariable lam=0;
    dvariable fpen=0.0;
    if (age_flags(171))
    {
      dvariable b0=calculate_the_biomass(1,*this,N_q0,age_flags(149),pmature);
      //dvariable b0=calculate_reproductive_biomass(N_q0,1);
#if !defined(NO_MY_DOUBLE_TYPE)
      lam=1.0+posfun((alpha*b0-beta)/b0-1.0L,.1,fpen);
#else
      lam=1.0+posfun((alpha*b0-beta)/b0-1.0,.1,fpen);
#endif
      //cout << setprecision(15) << alpha*b0/(beta+b0) << endl;
      if (lam<0)
      {
        cout << "negative equilibrium initial population"<< endl;
        cout << lam << endl;
        ad_exit(1);
      }
      dvariable loglam=log(lam);
      for (ir=1;ir<=num_regions;ir++)
      {
        for (int j=1;j<=nage;j++)
        {
          N_q0(ir,1,j)+=loglam;
        }
      }
    }
    cout << setscientific() << "XXX " << ssum*value(lam) << endl;
  */
  }
}

          
dvar_vector dvar_fish_stock_history::get_bh_recruitment
  (dvar3_array& NP,int iy)
{
  int lag=age_flags(147);
  int yr=max(1,iy-lag);
  int ir;
  dvar_vector v;
  if (iy<=lag) 
  {
    v=exp(recr(iy));
  }
  else
  {
    dvariable b0=calculate_the_biomass(yr,*this,N_q0,age_flags(149),pmature);
    //dvariable b0=calculate_reproductive_biomass(NP,yr);
    dvariable r=alpha*b0/(beta+b0);
    dvariable rtmp=r*exp(bh_recr_devs(iy-lag));
    v=epop_delta*rtmp;
  }
  v=log(v);

  if (!age_flags(71))
  {
    if (age_flags(101))
    {
      v+=sv(25)*rec_covars(iy);
    }
  }
  else
  {
    //v+=recr(iy)-recmean+rec_init_diff+region_rec_diffs(iy);
    v+=region_rec_diffs(iy);
    if (age_flags(101))
    {
      v+=sv(25)*rec_covars(iy);
    }
  }
  return exp(v);
}

dvariable dvar_fish_stock_history::get_bh_recruitment_multiplier
  (dvar3_array& NP,int iy)
{
  int lag=age_flags(147);
  int yr=max(1,iy-lag);
  int ir;
  int is=0;
  dvariable v;
  dvariable aalpha;  //NMD_12Apr2021
  dvariable bbeta;
  if (iy<=lag) 
  {
    v=1.0;
  }
  else
  {
    if (!pmsd || pmsd->current_species==1) //NMD_12Apr2021)
    {
      aalpha=alpha; bbeta=beta;
    }
    else
    {
      is=pmsd->current_species;
      aalpha=pmsd->alpha(is);
      bbeta=pmsd->beta(is);
    }
    dvariable b1=calculate_the_biomass(yr,*this,Nsave,age_flags(149),pmature);
    dvariable b0=calculate_the_biomass(yr,*this,N_q0,age_flags(149),pmature);
    dvariable r=aalpha*b0/(bbeta+b0);
    dvariable r1=aalpha*b1/(bbeta+b1);
    //f (b1<b0) 
    // cout <<"xequilib.cpp:"<<endl
    //      << b1 << "  " << b0 << endl 
    //      << r1 << "  " << r  << endl
    //      << r1 << "  " << bh_predicted_recruits(yr) << endl;
    if (!pmsd || pmsd->current_species==1) //NMD_12Apr2021
    {
      v=r/(1.e-20+bh_predicted_recruits(yr));
    }
    else
    {
      v=r/(1.e-20+pmsd->bh_predicted_recruits(is,yr));
    }
  }
 
  if (value(v)<1.00)
  {
    cout << setprecision(16) << v << " " << v-1.0 << endl;
    MY_DOUBLE_TYPE ssum=0.0;
    for (int ir=1;ir<=num_regions;ir++)
    {
      ssum+=sum(exp(value(N_q0(ir,yr))));
    }    
    cout << ssum << endl;
    cerr << "illegal value in get_bh_recruitment_multiplier"
         << endl;
    //ad_exit(1);
  }
 
  return log(v);
}
dvariable dvar_fish_stock_history::get_bh_annual_recruitment_multiplier
  (dvar3_array& NP,int iy)
{
  int lag=age_flags(147);
  int yr=max(1,iy-lag);
  int ns=age_flags(57);
  int annual_yr=(yr-1)/ns+1;
  int s1=(annual_yr-1)*ns+1;
  int ir;
  dvariable v;
  dvariable aalpha;  //NMD_12Apr2021
  dvariable bbeta;
  if (iy<=lag) 
  {
    v=1.0;
  }
  else
  {
    dvariable b0=0.0;
    int ic=0;
    for (int i=s1;i<=s1+ns-1;i++)
    {
      if (i>last_real_year)break;
      ic++;
      b0+=calculate_the_biomass(i,*this,NP,age_flags(149),pmature);
    }
    b0/=ic;
    if (!pmsd || pmsd->current_species==1) //NMD_12Apr2021)
    {
      aalpha=alpha; bbeta=beta;
    }
    else
    {
      int is=pmsd->current_species;
      aalpha=pmsd->alpha(is);
      bbeta=pmsd->beta(is);
    }

    dvariable r=aalpha*b0/(bbeta+b0);    
//    dvariable r=alpha*b0/(beta+b0);
//    v=r/(1.e-20+bh_predicted_recruits(annual_yr));
    if (!pmsd || pmsd->current_species==1) //NMD_12Apr2021
    {
      v=r/(1.e-20+bh_predicted_recruits(annual_yr));
    }
    else
    {
      int is=pmsd->current_species;
      v=r/(1.e-20+pmsd->bh_predicted_recruits(is,annual_yr));
    }
  }

  if (value(v)<1.00)
  {
    cout << setprecision(16) << v << " " << v-1.0 << endl;
    MY_DOUBLE_TYPE ssum=0.0;
    for (int ir=1;ir<=num_regions;ir++)
    {
      ssum+=sum(exp(value(NP(ir,yr))));
    }    
    cout << ssum << endl;
    cerr << "illegal value in get_bh_recruitment_multiplier"
         << endl;
    //ad_exit(1);
  }

  return log(v);
}


