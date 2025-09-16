/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
//% message to emacs: -*- mode: C++; fill-column: 80 -*-

#include "all.hpp"


extern int debug_flag;
  ostream & operator << (const ostream& _os, const print_double& _v)
  {
    ostream& os=(ostream&)(_os);
    print_double& v=(print_double&)(_v);
    os << v.to_double();
    return os;
  }

static dmatrix der_from_values(MY_DOUBLE_TYPE h,int k)
{
  int n=2*k;
  dmatrix A(0,n,0,n);
  int i;
  for (i=0;i<=n;i++)
  {
    A(i,0)=1.0;
    for (int j=1;j<=n;j++)
    {
      A(i,j)=A(i,j-1)*((i-k)*h)/j;
    }
  }
  return inv(A);
}


dvariable yield_analysis_bh_penalty(dvar_len_fish_stock_history& fsh,
          dvariable alpha,dvariable beta,const dvar_vector& pmature)
{
  //cout<<"START OF yield_analysis_bh_penalty"<<endl;
  if (debug_flag)
  {
    int n;
    cout << "HH alpha " << alpha << " beta  " << beta << endl;
    if (fsh.pmsd)
    {
      cout <<  "  species " << fsh.pmsd->current_species << endl;
    }
  }
  MY_DOUBLE_TYPE upper_bound=1.e+100;
  MY_DOUBLE_TYPE lower_bound=-1;
  dvariable lambda_save=-1.0;
  int nage=fsh.nage;
  int numflag=fsh.age_flags(149);
  int nyears=fsh.nyears;
  dvar_vector n1(1,nage);
  int sflag=0;
  if (fsh.pmsd)
  {
    sflag=sum(column(fsh.pmsd->species_flags,2));
  }

  int i;
  int navg=fsh.age_flags(148);
  int nomit=fsh.age_flags(155);
  int tmult=fsh.age_flags(57);
  if (!navg)navg=tmult;
  if (!navg)navg=1;

  if (!sflag)
  {
    fsh.get_fishing_mortality_by_age_by_year();  
  }
  else
  {
    fsh.get_fishing_mortality_by_age_by_year(1);  
    fsh.get_fishing_mortality_by_age_by_year(2);  
  }
  MY_DOUBLE_TYPE delta=1.e-3;
  MY_DOUBLE_TYPE delta2=delta*delta;
  MY_DOUBLE_TYPE targetratio=double(fsh.age_flags(165))/get_flag(fsh.age_flags(164),100.);

  int npoints=3;
  dmatrix A=der_from_values(delta,npoints);

  dvar_vector F(1,nage);

  dvar_vector FS;
  dvar_vector psmature;
  if (sflag)
  {
    FS.allocate(1,nage);
    FS.initialize();
    for (i=nomit+1;i<=navg;i++)
      FS+=fsh.pmsd->F_by_age_by_year(2,nyears-i+1);
    FS/=double(navg-nomit);
    psmature=fsh.get_pmature_species(2);
  }
  if (!fsh.pmsd || fsh.pmsd->current_species==1 )
  {
    F.initialize();
    for (i=nomit+1;i<=navg;i++)
      F+=fsh.F_by_age_by_year(nyears-i+1);
  }
  else
  {
    int cs=fsh.pmsd->current_species;
    for (i=nomit+1;i<=navg;i++)
      F+=fsh.pmsd->F_by_age_by_year(cs,nyears-i+1);
  }
  F/=double(navg-nomit);

  //cout <<"yieldbhp.cpp " << mean(F) << endl;
  dvar_vector yv(-npoints,npoints);
  dvariable Cn;
  MY_DOUBLE_TYPE target=double(fsh.age_flags(165))/100.0;
  dvariable lambda=exp(fsh.old_log_lambda); // start at target multiplier value
  if(lambda>10)lambda=10.0; //YT 2018-05-31
  dvariable fp1;
  dvariable fp2;
  dvariable Bmsy, sBmsy;
  if (!sflag)
  {
    yv(0)=get_yield_at_multiplier2(lambda,F,nyears,n1,pmature,
      alpha,beta,fsh,Bmsy,sBmsy);
  }
  else
  {
    yv(0)=get_yield_at_multiplier2(lambda,F,nyears,n1,pmature,
      alpha,beta,fsh,Bmsy,sBmsy,FS,psmature);
  }
  if (yv(0)<lambda_save && lambda_save>-0.5)
  {
    upper_bound=min(upper_bound,value(lambda_save));
  }
  if (yv(0)>lambda_save && lambda_save>-0.5)
  {
    lower_bound=max(lower_bound,value(lambda_save));
  }
    
  MY_DOUBLE_TYPE oy=value(yv(0));
  dvariable oldiff; 
  int icount=0;
  dvariable h;
  do
  {
    delta=value(lambda)-(value(lambda)-delta);
    dvariable lm1=lambda-delta;
    dvariable lm2=lambda-2.0*delta;
    dvariable lm3=lambda-3.0*delta;

    if (!sflag)
    {
      yv(-3)=get_yield_at_multiplier2(lm3,F,nyears,n1,pmature,
        alpha,beta,fsh,Bmsy,sBmsy);
    }
    else
    {
      yv(-3)=get_yield_at_multiplier2(lm3,F,nyears,n1,pmature,
        alpha,beta,fsh,Bmsy,sBmsy,FS,psmature);
    }
    if (!sflag)
    {
      yv(-2)=get_yield_at_multiplier2(lm2,F,nyears,n1,pmature,
        alpha,beta,fsh,Bmsy,sBmsy);
    }
    else
    {
      yv(-2)=get_yield_at_multiplier2(lm2,F,nyears,n1,pmature,
        alpha,beta,fsh,Bmsy,sBmsy,FS,psmature);
    }
    if (!sflag)
    {
      yv(-1)=get_yield_at_multiplier2(lm1,F,nyears,n1,pmature,
        alpha,beta,fsh,Bmsy,sBmsy);
    }
    else
    {
      yv(-1)=get_yield_at_multiplier2(lm1,F,nyears,n1,pmature,
        alpha,beta,fsh,Bmsy,sBmsy,FS,psmature);
    }

    dvariable lp1=lambda+delta;
    dvariable lp2=lambda+2.0*delta;
    dvariable lp3=lambda+3.0*delta;
    if (!sflag)
    {
      yv(3)=get_yield_at_multiplier2(lp3,F,nyears,n1,pmature,
        alpha,beta,fsh,Bmsy,sBmsy);
    }
    else
    {
      yv(3)=get_yield_at_multiplier2(lp3,F,nyears,n1,pmature,
        alpha,beta,fsh,Bmsy,sBmsy,FS,psmature);
    }
    if (!sflag)
    {
      yv(2)=get_yield_at_multiplier2(lp2,F,nyears,n1,pmature,
        alpha,beta,fsh,Bmsy,sBmsy);
    }
    else
    {
      yv(2)=get_yield_at_multiplier2(lp2,F,nyears,n1,pmature,
        alpha,beta,fsh,Bmsy,sBmsy,FS,psmature);
    }
    if (!sflag)
    {
      yv(1)=get_yield_at_multiplier2(lp1,F,nyears,n1,pmature,
        alpha,beta,fsh,Bmsy,sBmsy);
    }
    else
    {
      yv(1)=get_yield_at_multiplier2(lp1,F,nyears,n1,pmature,
        alpha,beta,fsh,Bmsy,sBmsy,FS,psmature);
    }
   
    dvar_vector tmp=A*yv.shift(0);
    yv.shift(-npoints);
   
   
    fp1=tmp(1);
    fp2=tmp(2);
 

    lambda_save=lambda;
    if (fp2< -.001)
    {
      h=-fp1/fp2;
      lambda+=h;
//      cout << " h = " << h <<  " lambda = " << lambda << endl;
    }
    else if (fp2>=0)
    {
      if (fp1>0)
      {
          lambda*=2.0;
      }    
      if (fp1<0)
      {
          lambda*=.5;
      }    

    }  

   // cout << " h = " << fp1 << "  "  << fp2 << "  " 
     //    << -fp1/fp2 <<  " lambda = " << lambda << endl;
    if (lambda<0) 
    {
      lambda=0.5*lambda_save;
    }
    yv(0)=get_yield_at_multiplier2(lambda,F,nyears,n1,pmature,
      alpha,beta,fsh,Bmsy,sBmsy);
    if (!sflag)
    {
      yv(0)=get_yield_at_multiplier2(lambda,F,nyears,n1,pmature,
        alpha,beta,fsh,Bmsy,sBmsy);
    }
    else
    {
      yv(0)=get_yield_at_multiplier2(lambda,F,nyears,n1,pmature,
        alpha,beta,fsh,Bmsy,sBmsy,FS,psmature);
    }
    oldiff = yv(0)-oy; 
    oy=value(yv(0));
//    if (fsh.pmsd)
//      cout << "A5  " <<  " species " << fsh.pmsd->current_species 
//         << "   " << yv(0) << endl;
    
    icount++;
  }
  while (fabs(value(h))>1.e-11 && icount<50);  
  
  fsh.old_log_lambda=value(log(lambda));
//  cout << yv(0) << endl;
  Cn=yv(0);
  
  dvariable Fmmsy=lambda;

  //--- collect some stuff for report
  fsh.p_Fmmsy=value(Fmmsy); 
  fsh.p_MSY=value(Cn);
  fsh.p_Bmsy=value(Bmsy);
  fsh.p_sBmsy=value(sBmsy);
  
  //---------------------------------------------------------------------
  // calc penalty on target Fmmsy, BBmsy, or sBsBmsy
  //---------------------------------------------------------------------
  dvariable pen = 0.0;

  if(fsh.age_flags(165)) 
  {
    //double wt=get_fmsy_pen_wt(fsh.age_flags(166));
    MY_DOUBLE_TYPE wt=get_flag(fsh.age_flags(166), 1000.);

    //fsh.p_Fmmsy=value(Fmmsy);
    if(fsh.age_flags(167) == 0) 
    {
      // compare Fmmsy (which is Fmsy/F) to target
      //cout<<"af165 "<<targetratio<< "pen wt = " << wt <<endl;
      //cout<<"Fmmsy/af165 "<<targetratio<<endl;
      pen += wt*square(log(Fmmsy/targetratio));
      fsh.p_Fmmsy_pen=value(pen);

      cout << "yieldbhp.cpp"<<endl;
      cout << "Target Fmsy/F ratio" 
           << setfixed() << setprecision(4) << setw(8) 
           << targetratio
           << " Actual Fmsy/F ratio " 
           << setfixed() << setprecision(4) << setw(8) 
           << Fmmsy <<endl;
      cout << " Target Penalty "   
           << setfixed() << setprecision(12) << setw(8) 
           << pen
           << endl;
      cout << "MSY: "<< Cn << endl;
    } 
    else if(fsh.age_flags(167) == 1)
    {
      // compare B/Bmsy to target
      //double tmsy=fsh.age_flags(165)/100.0;               
    //pen=biomass_target_calcs(fsh);
    pen=biomass_target_calcs(fsh,Bmsy,sBmsy,pmature);

      cout<<"====== yieldbhp ====="<<endl;
      cout << "Target B/B(msy) ratio" 
           << setfixed() << setprecision(4) << setw(8) 
           << targetratio
           << " Actual B/B(msy) ratio " 
           << setfixed() << setprecision(4) << setw(8) 
           << fsh.p_BBmsy <<endl;
      cout << " Target Penalty "       // PK 6-07
           << setfixed() << setprecision(12) << setw(8) 
           << pen
           << endl;
    }
    else if(fsh.age_flags(167) == 2)
    {
      // compare sB/sBmsy to target
      //double tmsy=fsh.age_flags(165)/100.0;               
      //pen=biomass_target_calcs(fsh);
      pen=biomass_target_calcs(fsh,Bmsy,sBmsy,pmature);

      cout<<"====== yieldbhp ====="<<endl;
      cout << "Target sB/sB(msy) ratio" 
           << setfixed() << setprecision(4) << setw(8) 
           << targetratio
           << " Actual sB/sB(msy) ratio " 
           << setfixed() << setprecision(4) << setw(8) 
           << fsh.p_sBsBmsy <<endl;
      cout << " Target Penalty "       // PK 6-07
           << setfixed() << setprecision(12) << setw(8) 
           << pen
           << endl;
    }
    else
    {
      cout << "target: " <<fsh.age_flags(167)<< " not implemented yet"<< endl;
      ad_exit(1);
    }
  }
  else
  {
    biomass_target_calcs(fsh,Bmsy,sBmsy,pmature); // for plot.rep file
  }
  return pen;
} 

dvariable get_yield_at_multiplier2(const dvariable& lambda, dvar_vector& F,
  int nyears,dvar_vector& n1,const dvar_vector& pmature,dvariable& alpha,
  dvariable& beta,dvar_len_fish_stock_history& fsh,
  dvariable& Bmsy,dvariable& sBmsy)
{
  int nage=fsh.nage;
  int numflag=fsh.age_flags(149);

  dvar_vector Z=lambda*F+exp(fsh.get_nat_mort_species());

  n1(1)=1.0;
  for (int j=2;j<=fsh.nage-1;j++)  
  {
    n1(j)=n1(j-1)*exp(-Z(j-1));
  }
  n1(nage)=n1(nage-1)*exp(-Z(nage-1))/(1.0-exp(-Z(nage)));
  
  //  Biomass
  dvariable ba, bt;
  int rr=fsh.get_region_bounds()(1);
  if (numflag)
  {
    bt=sum(n1);
    ba=pmature*n1;
  }
  else
  {
    bt=n1*fsh.mean_weight_yr(rr,1)/1000.;
    ba=pmature*(elem_prod(n1,fsh.mean_weight_yr(rr,1)))/1000.;
  }
  //-- bt is total biomass or nos. per recruit
  //-- ba is adult biomass or nos. per recruit
  dvariable recruitment;
  if (fsh.age_flags(161))
  {
    recruitment=alpha*mfexp(0.5*fsh.bh_variance)-beta/ba;
  }
  else
  {
    recruitment=alpha-beta/ba;
  }
// Adjustment for annualised BH-SRR   NMD_3July2015
  if ((fsh.age_flags(94)==3 || fsh.age_flags(182)) && fsh.age_flags(57)>1)
  {
    recruitment/=fsh.age_flags(57);
  }
//
  //cout<<" at QQQ: Bmsy "<<bt*recruitment<<"  sBmsy "<<ba*recruitment<<endl;
  Bmsy=bt*recruitment;  // equilib. total biomass (not necessarily Bmsy)
  sBmsy=ba*recruitment; // equilib. adult biomass (not necessarily sBmsy)

  dvar_vector C1=elem_prod(elem_div(lambda*F,Z),
     elem_prod(1.0-exp(-Z),n1));

  dvariable Cn;
  if (numflag)
    Cn=recruitment*sum(C1);
  else
    Cn=recruitment*(C1*fsh.mean_weight_yr(rr,1))/1000.;

  return Cn;
}

dvariable get_yield_at_multiplier2(const dvariable& lambda, 
  dvar_vector& _F,int nyears,
  dvar_vector& n1,const dvar_vector& _pmature,dvariable& alpha,
  dvariable& beta,
  dvar_len_fish_stock_history& fsh,dvariable& Bmsy,dvariable& sBmsy,
  dvar_vector& _FS,const dvar_vector& _psmature)
{
  ivector sp2=column(fsh.pmsd->species_flags,2);
  int mmin=sp2.indexmin();
  int mmax=sp2.indexmax();
  int isflag=0;
  for (int i=mmin;i<=mmax;i++)
  {
    if (sp2(i))
      isflag=i;
  }
  if (isflag==0)
  {
    cerr << "this can't happen" << endl;
    ad_exit(1);
  }
  dvar_vector F;
  dvar_vector FS;
  dvar_vector pmature;
  dvar_vector psmature;
  dvar_vector Z;
  dvar_vector ZS;

  int rr=0;
  int rs=0;
  switch (isflag)
  {
  case 1:
    F=_F;
    pmature=_pmature;
    FS=_FS;
    psmature=_psmature;
    Z=lambda*F+exp(fsh.get_nat_mort_species(1)(1));
    ZS=lambda*FS+exp(fsh.get_nat_mort_species(2)(1));
    rr=fsh.get_region_bounds(1)(1);
    rs=fsh.get_region_bounds(2)(1);
    break;
  case 2:
    F=_FS;
    pmature=_psmature;
    FS=_F;
    psmature=_pmature;
    Z=lambda*F+exp(fsh.get_nat_mort_species(2)(1));
    ZS=lambda*FS+exp(fsh.get_nat_mort_species(1)(1));
    rr=fsh.get_region_bounds(2)(1);
    rs=fsh.get_region_bounds(1)(1);
    break;
  default:
    cerr << "this can't happen" << endl;
    ad_exit(1);
  }
  int nage=fsh.nage;
  int numflag=fsh.age_flags(149);

  dvar_vector ns1(1,nage);


  n1(1)=1.0;
  for (int j=2;j<=fsh.nage-1;j++)  
  {
    n1(j)=n1(j-1)*exp(-Z(j-1));
  }
  n1(nage)=n1(nage-1)*exp(-Z(nage-1))/(1.0-exp(-Z(nage)));
  
  ns1(1)=1.0;
  for (int j=2;j<=fsh.nage-1;j++)  
  {
    ns1(j)=ns1(j-1)*exp(-ZS(j-1));
  }
  ns1(nage)=ns1(nage-1)*exp(-ZS(nage-1))/(1.0-exp(-ZS(nage)));
  
  //  Biomass
  dvariable ba, bt;
  dvariable bsa, bst;
  if (numflag)
  {
    bt=sum(n1);
    ba=pmature*n1;
    bst=sum(ns1);
    bsa=pmature*ns1;
  }
  else
  {
    bt=n1*fsh.mean_weight_yr(rr,1)/1000.;
    ba=pmature*(elem_prod(n1,fsh.mean_weight_yr(rr,1)))/1000.;
    bst=ns1*fsh.mean_weight_yr(rs,1)/1000.;
    bsa=psmature*(elem_prod(ns1,fsh.mean_weight_yr(rs,1)))/1000.;
  }
  //-- bt is total biomass or nos. per recruit
  //-- ba is adult biomass or nos. per recruit
  dvariable recruitment;
  if (fsh.age_flags(161))
  {
    recruitment=alpha*exp(0.5*fsh.bh_variance)-beta/ba;
  }
  else
  {
    recruitment=alpha-beta/ba;
  }
  //cout<<" at QQQ: Bmsy "<<bt*recruitment<<"  sBmsy "<<ba*recruitment<<endl;
  Bmsy=bt*recruitment;  // equilib. total biomass (not necessarily Bmsy)
  sBmsy=ba*recruitment; // equilib. adult biomass (not necessarily sBmsy)

  dvar_vector C1=elem_prod(elem_div(lambda*F,Z),
     elem_prod(1.0-exp(-Z),n1));

  dvar_vector CS1=elem_prod(elem_div(lambda*FS,ZS),
     elem_prod(1.0-exp(-ZS),ns1));

  dvariable Cn;
  if (numflag)
    Cn=recruitment*sum(C1);
  else
    Cn=recruitment*(C1*fsh.mean_weight_yr(rr,1))/1000.;

  dvariable CSn;
  if (numflag)
    CSn=recruitment*sum(CS1);
  else
    CSn=recruitment*(CS1*fsh.mean_weight_yr(rs,1))/1000.;

  return Cn+CSn;
}


dvariable biomass_target_calcs(dvar_len_fish_stock_history& fsh,
          dvariable & Bmsy, dvariable & sBmsy,const dvar_vector& pmature)
{
  int numflag=fsh.age_flags(149);
  int nyears=fsh.nyears;
  int navg=fsh.age_flags(148);
  int nomit=fsh.age_flags(155);
  int tot_or_adult=fsh.age_flags(167);
  int nstart=nyears-navg+1;
  int nend=nyears-nomit;
  //const MY_DOUBLE_TYPE targetratio=double(fsh.age_flags(165))/100.0;
  const MY_DOUBLE_TYPE targetratio=double(fsh.age_flags(165))/get_flag(fsh.age_flags(164),100.);
  dvar_vector abio(nstart,nend);  // Adult biomass
  dvar_vector tbio(nstart,nend);  // Total biomass
  tbio.initialize();
  abio.initialize();
  // calc avg. total nominal biomass over years given by af 148 & 155
  int iy,ir;
  int rmin=fsh.get_region_bounds()(1);
  int rmax=fsh.get_region_bounds()(2);
  for (ir=rmin;ir<=rmax;ir++)
  {
    for (iy=nstart;iy<=nend; iy++) 
    {
      // numflag determines if spawning stock is in numbers or weight
      if (numflag)  // --numbers--
      {
        abio(iy)+=sum(elem_prod(pmature,exp(fsh.N(ir,iy))));
        tbio(iy)+=sum(exp(fsh.N(ir,iy)));
      }
      else          // --weight--
      {
        abio(iy)+=pmature*(elem_prod(exp(fsh.N(ir,iy)),
          fsh.mean_weight_yr(ir,iy)))/1000.;
        tbio(iy)+=sum(elem_prod(exp(fsh.N(ir,iy)),
          fsh.mean_weight_yr(ir,iy)))/1000.;
      }
    }
  }
  dvariable Bt, Ba, BBmsy, sBsBmsy; 
  Bt = sum(tbio)/double(1.e-20+(nend-nstart+1));
  BBmsy = Bt/(1.e-20+Bmsy);
  Ba = sum(abio)/double(1.e-20+(nend-nstart+1));
  sBsBmsy = Ba/(1.e-20+sBmsy);

  //out<<" at SSS, nstart,nend: "
  //   << nstart<<" "<<nend
  //   <<"  Bt, Ba: "
  //   << setfixed() << setprecision(4) << setw(8) 
  //   <<Bt<<" "<<Ba
  //   <<"  Bmsy,sBmsy: "<<fsh.p_Bmsy<<" "<<fsh.p_sBmsy
  //   <<"  BBmsy,sBsBmsy "<<BBmsy<<" "<<sBsBmsy<<endl;
  //compare B/Bmsy to target

  MY_DOUBLE_TYPE wt=get_flag(fsh.age_flags(166), 1000.);
  dvariable pen=0;
  if(fsh.age_flags(167)==1) {
    pen = wt*square(log((targetratio+BBmsy)/(2.0*targetratio))); 
  }
  else if(fsh.age_flags(167)==2){
    pen = wt*square(log((targetratio+sBsBmsy)/(2.0*targetratio)));
  }
 
  // get values for report to par file and plot.rep file
  fsh.p_BBmsy=value(BBmsy); 
  fsh.p_sBsBmsy=value(sBsBmsy);
  //fsh.p_Fmmsy=value(Fmmsy); // this already done
  if(fsh.age_flags(167)==1) {
    fsh.p_BBmsy_pen=value(pen);
    fsh.p_sBsBmsy_pen=0.;
  }
  else if(fsh.age_flags(167)==2)  {
    fsh.p_BBmsy_pen=0.;
    fsh.p_sBsBmsy_pen=value(pen);
  }
  //cout<<" AT END OF biomass_target_calcs. BBmsy,sBsBmsy: "
  //    <<fsh.p_BBmsy<<" "<<fsh.p_sBsBmsy<<endl;

  return pen;
}
