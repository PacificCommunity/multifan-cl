/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

dvar_matrix dvar_fish_stock_history::catch_equations_calc_for_equilibrium
  (dvar_vector& sv, int navg,MY_DOUBLE_TYPE lambda,
    dvar_matrix& tmpcatch)
{
  pmsd_error();
  int ir;

  MY_DOUBLE_TYPE log_lambda=log(lambda);
  int current_year=nyears-navg+1;
  ivector rip(1,num_regions);
  ivector beginning_rip(1,num_regions);
  // set rip to point to the first fishing period in the n-r'th year 
  // for each region
  for (ir=1;ir<=num_regions;ir++)
  {
    int ip=num_fish_periods(ir);
    do
    { 
      ip--;
    }
    while(year(ir,ip)>=nyears-navg+1);
    beginning_rip(ir)=ip+1;
  } 
  dvar_matrix tmpfish(1,num_regions,1,nage);

  dvar3_array local_tot_mort(1,num_regions,
    beginning_rip,num_fish_periods,1,nage);

  dvar3_array local_survival(1,num_regions,
    beginning_rip,num_fish_periods,1,nage);


  local_tot_mort.initialize();
  tmpcatch.initialize();

  imatrix local_nfi(1,num_regions,beginning_rip,num_fish_periods);
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=beginning_rip(ir);ip<=num_fish_periods(ir);ip++)
    {
      local_nfi(ir,ip)=num_fish_incidents(ir,ip);
    }
  }

  dvar4_array local_fish_mort(1,num_regions,
    beginning_rip,num_fish_periods,1,local_nfi,1,nage);

  dvar4_array local_fish_mort_calcs(1,num_regions,
    beginning_rip,num_fish_periods,1,local_nfi,1,nage);

  dvar4_array local_catch(1,num_regions,
    beginning_rip,num_fish_periods,1,local_nfi,1,nage);

  
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=beginning_rip(ir);ip<=num_fish_periods(ir);ip++)
    {
      if (allocated(local_fish_mort(ir,ip)))
      {
        local_fish_mort(ir,ip)=log_lambda+fish_mort(ir,ip);
      }
      dvar_vector& tm=local_tot_mort(ir,ip);
      dvar_matrix& fm=local_fish_mort(ir,ip);
      
   
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        tm+=mfexp(fm(fi));
      }
      
      tm+=mfexp(get_nat_mort_region(ir)(year(ir,ip))+fraction(ir,ip));
    
      local_survival(ir,ip)=mfexp(-tm);
    }
  }

  const MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=beginning_rip(ir);ip<=num_fish_periods(ir);ip++)
    {
      const dvar_vector& tmp1=log(1.e-10+local_tot_mort(ir,ip));
      const dvar_vector& tmp2=log(one_plus-local_survival(ir,ip));
      const dvar_vector& tmp3=tmp2-tmp1;
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {        
        local_fish_mort_calcs(ir,ip,fi)=local_fish_mort(ir,ip,fi)+tmp3;
      }
    }
  } 



     
  int finished_flag=1;
  int break_flag=0;
  ivector tmp_mp(1,num_regions);
  ivector tmp_ip(1,num_regions);
  dvar3_array enum_fish(1,num_regions,beginning_rip,num_fish_periods,1,nage);
  enum_fish=-100;
  for (ir=1;ir<=num_regions;ir++)
  {
    enum_fish(ir,beginning_rip(ir),1)=pop_delta(ir);
  }


  for (int icount=1;icount<=20;icount++)
  {
    rip=beginning_rip;
    do
    {
      finished_flag=1;
      for (int ir=1;ir<=num_regions;ir++)
      {
        break_flag=0;
        int& ip=rip(ir);
        do
        {
          if (ip>=num_fish_periods(ir)) break;
          finished_flag=0;
          if (year(ir,ip+1)==year(ir,ip))
          {
            enum_fish(ir,ip+1)=enum_fish(ir,ip)-local_tot_mort(ir,ip);
          }
          else
          {
            // age the fish
            enum_fish(ir,ip+1,1)=pop_delta(ir);
  
            if (nage>2)
            --enum_fish(ir,ip+1)(2,nage-1)=
              enum_fish(ir,ip)(1,nage-2)-local_tot_mort(ir,ip)(1,nage-2);
  
            enum_fish(ir,ip+1,nage)=
              log(1.e-10 + mfexp(enum_fish(ir,ip,nage-1)
                -local_tot_mort(ir,ip,nage-1))
                + mfexp(enum_fish(ir,ip,nage)-local_tot_mort(ir,ip,nage)) );
  
          }
          if (move_flags(ir,ip))
          {
            tmp_ip(ir)=ip;
            tmp_mp(ir)=move_index(ir,ip);
            break_flag=1;
          } 
          ip++;
        }
        while (!break_flag); 
      }
      // move the fish
      {
        if (num_regions>1)
        {
          int ir;
          check_sanity(tmp_mp);
          dvar_matrix tmp=fast_diffusion_calcs(nage,num_regions,enum_fish,
            Dad(tmp_mp(1)),rip);
          for (ir=1;ir<=num_regions;ir++)
          {
            enum_fish(ir,rip(ir))=log(5.e-10+tmp(ir));
          }
        }
      }
    } // need to decide when to quit
    while (!finished_flag);
    //cout << "Finished " << endl;
    //greport("B catch_equations_calc");
    
    cout << "Iteration " << icount << endl;
    for (ir=1;ir<=num_regions;ir++)
    {
      tmpfish(ir)(1)=pop_delta(ir);
      tmpfish(ir)(2,nage)=++enum_fish(ir,num_fish_periods(ir))(1,nage-1);
      tmpfish(ir)(nage)=log(exp(tmpfish(ir)(nage))
        +exp(enum_fish(ir,num_fish_periods(ir))(nage)));
      enum_fish(ir)=-100;
      enum_fish(ir,beginning_rip(ir))=tmpfish(ir);
      enum_fish(ir,beginning_rip(ir),1)=pop_delta(ir);
    }
  } // for icount
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=beginning_rip(ir);ip<=num_fish_periods(ir);ip++)
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {        
        tmpcatch(ir)+=mfexp(local_fish_mort_calcs(ir,ip,fi)
          +tmpfish(ir));
      }
    }
  }
  tmpcatch/=double(navg);
  return exp(tmpfish);
}
    void xxx_aaa(dvar_vector& x){;}

dvariable yield_analysis_bh_daves_folly(dvar_len_fish_stock_history& fsh,
 ofstream * pof,dvariable alpha,dvariable beta,
 MY_DOUBLE_TYPE lwc,const dvar_vector& pmature)
{
  //cout<<"ENTERING yield_analysis_bh_daves_folly"<<endl;
  int ir;
  //int navg=4; // number of terminal years fishing mortality to use
     // in equilibrium calculations
  int nage=fsh.nage;
  int numflag=fsh.age_flags(149);
  int nyears=fsh.nyears;
  if(pof) (*pof) << "# Beverton-Holt yield analysis report" << endl;
  //const int nsteps=20000;
  int nsteps=fsh.age_flags(141);
  if(nsteps==0) nsteps=200;
  dvector tmp1(0,nsteps);
  dvar_vector tmp2(0,nsteps);
  dvar_vector tmp3(0,nsteps);  //JH 21/03/02 - equilib. adult biomass index
  dvar_vector tmp4(0,nsteps);  //JH 21/03/02 - equilib. total biomass index
  dvar_matrix tmp5(0,nsteps,1,fsh.num_regions); 
  dvar_matrix tmp6(0,nsteps,1,fsh.num_regions); 
  dvar_matrix tmp7(0,nsteps,1,fsh.num_regions); 
  dvar_matrix tmpcatch(1,fsh.num_regions,1,fsh.nage);
  tmp1.initialize();
  tmp2.initialize();
  tmp3.initialize();
  tmp4.initialize();
  tmp5.initialize();
  int i;
  int ii=1;
  int lastFmlt;
  int navg=fsh.age_flags(148);
  int nomit=fsh.age_flags(155);
  int tmult=fsh.age_flags(57);
  if (!navg)navg=tmult;
  if (!navg)navg=1;

  //n1 is the equilibrium population distribution in num of fish
  // for fish mort over last navg years
  dvar_matrix n1=fsh.get_equilibrium_age_structure(navg,nomit,1.e-10,tmpcatch);
    
  dvariable b0=0.0;  // This is the adult biomass
  dvariable tb0=0.0; // This is for total biomass
  if (numflag)
  {
    for (ir=1;ir<=fsh.num_regions;ir++)
    {
      b0+=pmature*n1(ir);
    }
    tb0=sum(n1); 
  }
  else
  {
    for (ir=1;ir<=fsh.num_regions;ir++)
    {
      b0+=pmature*(elem_prod(n1(ir),fsh.mean_weight_yr(1,1)))/1000.;
      tb0+=sum(elem_prod(n1(ir),fsh.mean_weight_yr(1,1)))/1000.; 
    }
  }
  dvariable n;
  if (fsh.age_flags(161))
  {
    n=alpha*exp(0.5*fsh.bh_variance)-beta/b0;
  }
  else
  {
    // this is equilibrium recruitment
    n=alpha-beta/b0;
  }
  tmp3(0)=b0*n;  // This scales up to the true equil. biomass
  tmp4(0)=tb0*n;

  lastFmlt = nsteps;
  for (i=1;i<=nsteps;i++)
  {
    //cout<<"IN FOLLY at ca line 65, step:"<<i<<"  out of "<<nsteps<<endl;
    //double lambda=i/100.0; //SDH 19-02-08 more effmult resolution
    //double lambda=i/10.0;
    //double lambda=i/2.0;
    MY_DOUBLE_TYPE lambda=i*fsh.age_flags(140)/100.0;

    dvariable Cn;
    dvar_vector Cn_by_region(1,fsh.num_regions);
    dvar_vector ABn_by_region(1,fsh.num_regions);
    dvar_vector TBn_by_region(1,fsh.num_regions);
    dvariable TBn;
    dvariable ABn;
    get_yield_at_effort_daves_folly(fsh,lambda,ABn,TBn,Cn,
      Cn_by_region,
      ABn_by_region,
      TBn_by_region,
      pmature,alpha,beta);
    if (fsh.age_flags(169)<0 && value(Cn)<0.0)  
    {
      lastFmlt = ii-1; // for later output of just positive part of
      break;
    }
    tmp1(ii)=lambda;        
    tmp3(ii)=ABn;
    tmp4(ii)=TBn;
    tmp2(ii)=Cn;
    tmp5(ii)=Cn_by_region;
    tmp6(ii)=ABn_by_region;
    tmp7(ii)=TBn_by_region;
    ii++;
  }
  // find the index which has the maximum yield
  int imax=1;
  {
    MY_DOUBLE_TYPE Cmax=value(tmp2(1));
    for (int i=2;i<=lastFmlt;i++)
    {
      if (tmp2(i)>Cmax)
      {
        Cmax=value(tmp2(i));
        imax=i;
      }
    }
  }  
  // JH 27/03/02 Calculates ratios of total and adult biomass (or numbers)
  //             to total and adult biomass (or numbers) at MSY
  // ---------------------------------------------------------------------
  const dvector & _tmp2=value(tmp2);
  dvector & ttmp2=(dvector &) _tmp2;
  const dvector & _tmp3=value(tmp3);
  dvector & ttmp3=(dvector &) _tmp3;
  const dvector & _tmp4=value(tmp4);
  dvector & ttmp4=(dvector &) _tmp4;
  get_biomass_ratio_by_time(fsh,numflag,pmature,ttmp2,ttmp3,
   ttmp4);
  //---------------------------------------------------------------------
  // PK 2003-06-10, penalizable Fmult at MSY by weighted avg.
  //---------------------------------------------------------------------
  dvariable pen = 0.0;
  if(fsh.age_flags(165)) 
  {
    dvariable exponent = fsh.age_flags(168) ? fsh.age_flags(168) : 50;
    dvariable divider = fsh.age_flags(169) ? fsh.age_flags(169) : 100000;
    dvariable Fmmsy;
    if (fsh.age_flags(169)>=0)  
    {
      dvar_vector wts = pow(exp(tmp2/divider),exponent);
      Fmmsy = tmp1*wts/sum(wts);
    }
    else
    {
      dvar_vector wts = pow(tmp2(0,lastFmlt)/tmp2(imax),exponent);
      Fmmsy = tmp1(0,lastFmlt)*wts/sum(wts);
    }
   
    cout << "Smoothed Fmmsy = " << Fmmsy << "\n"
         << "\"True\" Fmmsy     = " << tmp1(imax) << endl;
    dvariable wt = fsh.age_flags(166) ? float(fsh.age_flags(166))/100. : 1000.0;
    if(fsh.age_flags(167) == 0) 
    {
      // compare Fmmsy (which is Fmsy/F) to target
      cout<<"af165 "<<float(fsh.age_flags(165))/100<<endl;
      cout<<"Fmmsy/af165 "<<Fmmsy/(float(fsh.age_flags(165))/100)<<endl;
      pen += wt*square(log(Fmmsy/(float(fsh.age_flags(165))/100.0)));
      fsh.p_Fmmsy=value(Fmmsy);
      fsh.p_Fmmsy_pen=value(pen);
  
      cout << "test_msy.cpp"<<endl;
      cout << "Target Fmsy/F ratio" 
           << setfixed() << setprecision(4) << setw(8) 
           << (fsh.age_flags(165)/100.0)
           << " Actual Fmsy/F ratio " 
           << setfixed() << setprecision(4) << setw(8) 
           << Fmmsy
           << " Penalty " 
           << setfixed() << setprecision(2) << setw(8) 
           << pen
           << endl;
      cout << "MSY: "<< max(tmp2)<< endl;
    } 
    else
    {
      dvariable Cn;
      dvar_vector Cn_by_region(1,fsh.num_regions);
      
      dvariable tb1=0.0;  // This is for total biomass
      dvariable b1=0.0;  // This is for adult biomass
      get_yield_at_effort_daves_folly(fsh,Fmmsy,b1,tb1,Cn,Cn_by_region,
        pmature,alpha,beta);

      dvariable TBmsy; 
  
      dvariable n;
      if (fsh.age_flags(161))
      {
        TBmsy = tb1*(alpha*exp(0.5*fsh.bh_variance)-beta/tb1);  // This is the scaled-up total Bmsy
      }
      else
      {
        TBmsy = tb1*(alpha-beta/tb1);  // This is the scaled-up total Bmsy
      }
      // JH 27/03/02 Calculates ratios of total and adult biomass (or numbers)
      //             to total and adult biomass (or numbers) at MSY
      // ---------------------------------------------------------------------
      dvar_vector abio(1,fsh.nyears);  // Adult biomass
      dvar_vector tbio(1,fsh.nyears);  // Total biomass
      dvar_vector tcat(1,fsh.nyears);  // Total catch
      abio.initialize();
      tbio.initialize();
      for (int ir=1;ir<=fsh.num_regions;ir++)
      {
        for (int iy=1;iy<=fsh.nyears;iy++)  
        {
          // numflag determines if spawning stock is in numbers or weight
          if (numflag)
          {
            abio(iy)+=sum(elem_prod(pmature,exp(fsh.N(ir,iy))));
            tbio(iy)+=sum(exp(fsh.N(ir,iy)));
          }
          else
          {
            abio(iy)+=pmature*(elem_prod(exp(fsh.N(ir,iy)),
              fsh.mean_weight_yr(ir,iy)))/1000.;
            tbio(iy)+=sum(elem_prod(exp(fsh.N(ir,iy)),
              fsh.mean_weight_yr(ir,iy)))/1000.;
          }
        }
      }
      // calc avg. total nominal biomass over years given by af 148 & 155
      dvariable nbio = 0.0;
      for (i=nomit+1;i<=navg;i++) nbio += tbio(fsh.nyears-i+1);
      nbio /=(MY_DOUBLE_TYPE)( navg-nomit);
      dvariable BBmsy = nbio/TBmsy;
      // compare B/Bmsy to target
      pen += wt*square(log(BBmsy/(fsh.age_flags(165)/100.0)));               
      fsh.p_BBmsy=value(BBmsy);
      fsh.p_BBmsy_pen=value(pen);

      cout << "Target B/B(msy) ratio" 
           << setfixed() << setprecision(4) << setw(8) 
           << (fsh.age_flags(165)/100.0)
           << " Actual B/B(msy) ratio " 
           << setfixed() << setprecision(4) << setw(8) 
           << BBmsy
           << " Penalty " 
           << setfixed() << setprecision(2) << setw(8) 
           << pen
           << endl;
    }
  }
  int imsy=max_index(tmp2);                // gives the index of max(tmp2)
  dvariable tbio_at_msy=tmp4(imsy);        // equil. total biomass at MSY
  dvariable abio_at_msy=tmp3(imsy);        // equil. adult biomass at MSY
  dvariable MSY=tmp2(imsy)+1.e-10;                // MSY
  dvar_vector F_agg(1,fsh.nyears);         // Aggregate F by year
  dvariable Fmsy=MSY/(tbio_at_msy+1.e-10); // Aggregate F at MSY (F=catch/biomass)
  if (pof) 
  {
    fsh.yield_eq_region=value(tmp5);
    fsh.yield_eq_region_x=tmp1;
    const MY_DOUBLE_TYPE & _tmp6=value(Fmsy);
    MY_DOUBLE_TYPE & ttmp6=(MY_DOUBLE_TYPE &) _tmp6;
    const dvector & _tmp9=value(F_agg);
    dvector & ttmp9=(dvector &) _tmp9;
    (*pof) <<"# MSY by region"<<endl;
    (*pof) << tmp5(imsy) << endl;
    (*pof) << "# Total biomass by region at MSY"   <<endl;
    (*pof) << tmp7(imsy) << endl;
    (*pof) <<  "# Adult biomass by region at MSY"   <<endl;
    (*pof) << tmp6(imsy) << endl;
    (*pof) <<"# MSY"<<endl;
    (*pof) << MSY <<endl;
    (*pof) <<  "# F multiplier at MSY" <<endl;
    (*pof) <<  tmp1(imsy) <<endl;
    (*pof) <<  "# F at MSY"  <<endl;
    (*pof) <<  Fmsy  <<endl;
    (*pof) <<  "# Total biomass at MSY" <<endl;
    (*pof) <<  tbio_at_msy <<endl;
    (*pof) <<  "# Adult biomass at MSY" <<endl;
    (*pof) <<  abio_at_msy <<endl;

    print_yield_stuff(pof,fsh,lastFmlt, F_agg, 
                        F_agg, F_agg, F_agg,imsy); //last 4 args dummy
    //value(MSY),ttmp6,ttmp7,ttmp8,ttmp9,imsy,ttemp9);
  }
  return pen;
} 

dvar_matrix dvar_fish_stock_history::catch_equations_calc_for_equilibrium
  (dvar_vector& sv, int navg,const prevariable & lambda,
    dvar_matrix& tmpcatch)
{
  int ir;
  dvariable log_lambda=log(lambda);
  int current_year=nyears-navg+1;
  ivector rip(1,num_regions);
  ivector beginning_rip(1,num_regions);
  // set rip to point to the first fishing period in the n-r'th year 
  // for each region
  for (ir=1;ir<=num_regions;ir++)
  {
    int ip=num_fish_periods(ir);
    do
    { 
      ip--;
    }
    while(year(ir,ip)>=nyears-navg+1);
    beginning_rip(ir)=ip+1;
  } 
  dvar_matrix tmpfish(1,num_regions,1,nage);

  dvar3_array local_tot_mort(1,num_regions,
    beginning_rip,num_fish_periods,1,nage);

  dvar3_array local_survival(1,num_regions,
    beginning_rip,num_fish_periods,1,nage);


  local_tot_mort.initialize();
  tmpcatch.initialize();

  imatrix local_nfi(1,num_regions,beginning_rip,num_fish_periods);
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=beginning_rip(ir);ip<=num_fish_periods(ir);ip++)
    {
      local_nfi(ir,ip)=num_fish_incidents(ir,ip);
    }
  }

  dvar4_array local_fish_mort(1,num_regions,
    beginning_rip,num_fish_periods,1,local_nfi,1,nage);

  dvar4_array local_fish_mort_calcs(1,num_regions,
    beginning_rip,num_fish_periods,1,local_nfi,1,nage);

  dvar4_array local_catch(1,num_regions,
    beginning_rip,num_fish_periods,1,local_nfi,1,nage);

  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=beginning_rip(ir);ip<=num_fish_periods(ir);ip++)
    {
      local_fish_mort(ir,ip)=log_lambda+fish_mort(ir,ip);
      dvar_vector& tm=local_tot_mort(ir,ip);
      dvar_matrix& fm=local_fish_mort(ir,ip);
      
   
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        tm+=mfexp(fm(fi));
      }
      
      tm+=mfexp(get_nat_mort_region(ir)(year(ir,ip))+fraction(ir,ip));
    
      local_survival(ir,ip)=mfexp(-tm);
    }
  }

  const MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=beginning_rip(ir);ip<=num_fish_periods(ir);ip++)
    {
      const dvar_vector& tmp1=log(1.e-10+local_tot_mort(ir,ip));
      const dvar_vector& tmp2=log(one_plus-local_survival(ir,ip));
      const dvar_vector& tmp3=tmp2-tmp1;
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {        
        local_fish_mort_calcs(ir,ip,fi)=local_fish_mort(ir,ip,fi)+tmp3;
      }
    }
  } 


     
  int finished_flag=1;
  int break_flag=0;
  ivector tmp_mp(1,num_regions);
  ivector tmp_ip(1,num_regions);
  dvar3_array enum_fish(1,num_regions,beginning_rip,num_fish_periods,1,nage);
  enum_fish=-100;
  for (ir=1;ir<=num_regions;ir++)
  {
    enum_fish(ir,beginning_rip(ir),1)=pop_delta(ir);
  }


  for (int icount=1;icount<=20;icount++)
  {
    rip=beginning_rip;
    do
    {
      finished_flag=1;
      for (int ir=1;ir<=num_regions;ir++)
      {
        break_flag=0;
        int& ip=rip(ir);
        do
        {
          if (ip>=num_fish_periods(ir)) break;
          finished_flag=0;
          if (year(ir,ip+1)==year(ir,ip))
          {
            enum_fish(ir,ip+1)=enum_fish(ir,ip)-tot_mort(ir,ip);
          }
          else
          {
            // age the fish
            enum_fish(ir,ip+1,1)=pop_delta(ir);
  
            if (nage>2)
            --enum_fish(ir,ip+1)(2,nage-1)=
              enum_fish(ir,ip)(1,nage-2)-tot_mort(ir,ip)(1,nage-2);
  
            enum_fish(ir,ip+1,nage)=
              log(1.e-10 + mfexp(enum_fish(ir,ip,nage-1)-tot_mort(ir,ip,nage-1))
                + mfexp(enum_fish(ir,ip,nage)-tot_mort(ir,ip,nage)) );
  
          }
          if (move_flags(ir,ip))
          {
            tmp_ip(ir)=ip;
            tmp_mp(ir)=move_index(ir,ip);
            break_flag=1;
          } 
          ip++;
        }
        while (!break_flag); 
      }
      // move the fish
      {
        if (num_regions>1)
        {
          int ir;
          check_sanity(tmp_mp);
          dvar_matrix tmp=fast_diffusion_calcs(nage,num_regions,enum_fish,
            Dad(tmp_mp(1)),rip);
          for (ir=1;ir<=num_regions;ir++)
          {
            enum_fish(ir,rip(ir))=log(5.e-10+tmp(ir));
          }
        }
      }
    } // need to decide when to quit
    while (!finished_flag);
    //cout << "Finished " << endl;
    //greport("B catch_equations_calc");
    
    cout << "Iteration " << icount << endl;
    for (ir=1;ir<=num_regions;ir++)
    {
      tmpfish(ir)(1)=pop_delta(ir);
      tmpfish(ir)(2,nage)=++enum_fish(ir,num_fish_periods(ir))(1,nage-1);
      tmpfish(ir)(nage)=log(exp(tmpfish(ir)(nage))
        +exp(enum_fish(ir,num_fish_periods(ir))(nage)));
      enum_fish(ir)=-100;
      enum_fish(ir,beginning_rip(ir))=tmpfish(ir);
      enum_fish(ir,beginning_rip(ir),1)=pop_delta(ir);
    }
  } // for icount
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=beginning_rip(ir);ip<=num_fish_periods(ir);ip++)
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {        
        tmpcatch(ir)+=mfexp(local_fish_mort_calcs(ir,ip,fi)
          +tmpfish(ir,ip));
      }
    }
  }
  tmpcatch/=double(navg);
  return exp(tmpfish);
}

void get_yield_at_effort_daves_folly(dvar_len_fish_stock_history& fsh,
  MY_DOUBLE_TYPE lambda,prevariable & ABn,prevariable & TBn,prevariable & Cn,
  dvar_vector& Cn_by_region,
  dvar_vector& ABn_by_region,
  dvar_vector& TBn_by_region,
  const dvar_vector& pmature,prevariable& alpha,prevariable& beta)
{
  int navg=fsh.age_flags(148);
  int nomit=fsh.age_flags(155);
  int numflag=fsh.age_flags(149);
  int tmult=fsh.age_flags(57);
  if (!navg)navg=tmult;
  if (!navg)navg=1;
  int ir;
  dvar_matrix tmpcatch(1,fsh.num_regions,1,fsh.nage);
  dvar_matrix n1=fsh.get_equilibrium_age_structure(navg,nomit,lambda,
    tmpcatch);
  //dvar_matrix n1=fsh.catch_equations_calc_for_equilibrium(fsh.sv,
  // navg,lambda, tmpcatch);
    
  dvariable b1=0.0;   // Adult biomass
  dvar_vector b1_by_region(1,fsh.num_regions);   // Adult biomass
  dvar_vector tb1_by_region(1,fsh.num_regions);   // Adult biomass
  dvariable tb1=0.0;  // This is for total biomass
  b1_by_region.initialize();   // Adult biomass
  tb1_by_region.initialize();   // Adult biomass
  if (numflag)
  {
    for (ir=1;ir<=fsh.num_regions;ir++)
    {
      b1+=pmature*n1(ir);
      b1_by_region(ir)=pmature*n1(ir);
      tb1_by_region(ir)=sum(n1(ir));
    }
    tb1=sum(n1); 
  }
  else
  {
    for (ir=1;ir<=fsh.num_regions;ir++)
    {
      b1_by_region(ir)=
        pmature*(elem_prod(n1(ir),fsh.mean_weight_yr(1,1)))/1000.;
      b1+=b1_by_region(ir);
      tb1_by_region(ir)=
        sum(elem_prod(n1(ir),fsh.mean_weight_yr(1,1)))/1000.; 
      tb1+=tb1_by_region(ir);
    }
  }
  dvariable n;
  if (fsh.age_flags(161))
  {
    n=alpha*exp(0.5*fsh.bh_variance)-beta/b1;
  }
  else
  {
    n=alpha-beta/b1;
  }

  dvar_vector C1(1,fsh.nage);
  dvar_matrix C1_by_region(1,fsh.num_regions,1,fsh.nage);
  C1.initialize();
  C1_by_region.initialize();
  for (ir=1;ir<=fsh.num_regions;ir++)
  {
    C1+=tmpcatch(ir);
    C1_by_region(ir)=tmpcatch(ir);
  }

  if (numflag)
    Cn=n*sum(C1);
  else
    Cn=n*(C1*fsh.mean_weight_yr(1,1))/1000.;

  for (ir=1;ir<=fsh.num_regions;ir++)
  {
    if (numflag)
      Cn_by_region(ir)=n*sum(C1_by_region(ir));
    else
      Cn_by_region(ir)=n*(C1_by_region(ir)*fsh.mean_weight_yr(1,1))/1000.;
  }

  TBn=n*tb1;  // This is the scaled-up total biomass
  ABn=n*b1;   // and adult biomass
  
  TBn_by_region=n*tb1_by_region;
  ABn_by_region=n*b1_by_region;
}

void get_yield_at_effort_daves_folly(dvar_len_fish_stock_history& fsh,
  prevariable& lambda,prevariable & ABn,prevariable & TBn,prevariable & Cn,
  dvar_vector& Cn_by_region,
  const dvar_vector& pmature,prevariable& alpha,prevariable& beta)
{
  int navg=fsh.age_flags(148);
  int nomit=fsh.age_flags(155);
  int numflag=fsh.age_flags(149);
  int tmult=fsh.age_flags(57);
  if (!navg)navg=tmult;
  if (!navg)navg=1;
  int ir;
  dvar_matrix tmpcatch(1,fsh.num_regions,1,fsh.nage);
  dvar_matrix n1=fsh.get_equilibrium_age_structure(navg,nomit,lambda,
    tmpcatch);
    
  dvariable b1=0.0;   // Adult biomass
  dvariable tb1=0.0;  // This is for total biomass
  if (numflag)
  {
    for (ir=1;ir<=fsh.num_regions;ir++)
    {
      b1+=pmature*n1(ir);
    }
    tb1=sum(n1); 
  }
  else
  {
    for (ir=1;ir<=fsh.num_regions;ir++)
    {
      b1+=pmature*(elem_prod(n1(ir),fsh.mean_weight_yr(1,1)))/1000.;
      tb1+=sum(elem_prod(n1(ir),fsh.mean_weight_yr(1,1)))/1000.; 
    }
  }
  dvariable n;
  if (fsh.age_flags(161))
  {
    n=alpha*exp(0.5*fsh.bh_variance)-beta/b1;
  }
  else
  {
    n=alpha-beta/b1;
  }


  dvar_vector C1(1,fsh.nage);
  dvar_matrix C1_by_region(1,fsh.num_regions,1,fsh.nage);
  C1.initialize();
  C1_by_region.initialize();
  for (ir=1;ir<=fsh.num_regions;ir++)
  {
    C1+=tmpcatch(ir);
    C1_by_region(ir)=tmpcatch(ir);
  }

  if (numflag)
    Cn=n*sum(C1);
  else
    Cn=n*(C1*fsh.mean_weight_yr(1,1))/1000.;

  for (ir=1;ir<=fsh.num_regions;ir++)
  {
    if (numflag)
      Cn_by_region(ir)=n*sum(C1_by_region(ir));
    else
      Cn_by_region(ir)=n*(C1_by_region(ir)*fsh.mean_weight_yr(1,1))/1000.;
  }

  TBn=n*tb1;  // This is the scaled-up total biomass
  ABn=n*b1;   // and adult biomass
}
