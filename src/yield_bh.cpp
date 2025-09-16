/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

//% message to emacs: -*- mode: C++; fill-column: 80 -*-

#include "all.hpp"

//---This is a (mostly) constant version of get_yield_at_multiplier2 in 
//    yieldbhp.cpp.  Used here for report generaton----
dvariable get_yield_at_multiplier2(MY_DOUBLE_TYPE lambda, dvar_vector& _F,int nyears,
  dvar_vector& n1, const dvar_vector& _pmature,dvariable& alpha,dvariable& beta,
  dvar_len_fish_stock_history& fsh,dvariable& Bmsy,dvariable& sBmsy,
  dvar_vector& _FS,const dvar_vector& _psmature);

dvariable get_yield_at_multiplier2(MY_DOUBLE_TYPE lambda, dvar_vector& F,int nyears,
  dvar_vector& n1, const dvar_vector& pmature,dvariable& alpha,dvariable& beta,
  dvar_len_fish_stock_history& fsh,dvariable& Bmsy,dvariable& sBmsy)
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
    // this is equilibrium recruitment
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

dvariable yield_analysis_bh(dvar_len_fish_stock_history& fsh,
 ofstream * pof,dvariable alpha,dvariable beta,
 const dvar_vector& _pmature)
{
  ADUNCONST(dvar_vector,pmature)
  if(!(pof)) {
    cout<<"SHOULD NOT GET TO yield_analysis_bh W/O pof ASSIGNED"<<endl;
    ad_exit(1);
  }
  int nage=fsh.nage;
  int numflag=fsh.age_flags(149);
  int sflag=0;
  if (fsh.pmsd)
  {
    sflag=sum(column(fsh.pmsd->species_flags,2));
  }

  int nyears=fsh.nyears;
  dvar_vector n1(1,nage);

  (*pof) << "# Beverton-Holt yield analysis report" << endl;   //NMD 20 Feb 2014

  const int nsteps=5000;
  dvector tmp1(0,nsteps);
  dvar_vector tmp2(0,nsteps);
  dvar_vector tmp3(0,nsteps);  //JH 21/03/02 - equilib. adult biomass index
  dvar_vector tmp4(0,nsteps);  //JH 21/03/02 - equilib. total biomass index
  dvector lambda_values(0,nsteps); 
  //dmatrix tmp5;
  tmp1.initialize();
  tmp2.initialize();
  tmp3.initialize();
  tmp4.initialize();
  int i;
  int ii=1;
  int lastFmlt=nsteps;
  int navg=fsh.age_flags(148);
  int nomit=fsh.age_flags(155);
  int tmult=fsh.age_flags(57);
  if (!navg)navg=tmult;
  if (!navg)navg=1;


  if ((navg-nomit)<1)
  {
    cerr << "illegal values for af(148) and af(155)" << endl
         << "using different nomit value" << endl;
    nomit=navg-1;
  }
  if (!sflag)
  {
    fsh.get_fishing_mortality_by_age_by_year();  
  }
  else
  {
    fsh.get_fishing_mortality_by_age_by_year(1);  
    fsh.get_fishing_mortality_by_age_by_year(2);  
  }

  int cs=fsh.get_current_species();
  int ng=fsh.get_current_nage();
  dvar_vector F(1,nage);
  F.initialize();
  dvar_vector FS;
  dvar_vector psmature;
  /*   //NMD_19Dec2013   - do later when optimising yields
  if (sflag)
  {
    FS.allocate(1,nage);
    FS.initialize();
    for (i=nomit+1;i<=navg;i++)
      FS+=fsh.pmsd->F_by_age_by_year(2,nyears-i+1);
    FS/=double(navg-nomit);
    psmature=fsh.get_pmature_species(2);
  }
  */
  if (!fsh.pmsd || fsh.pmsd->current_species==1 )
  {
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

  // JH 21/03/02  This is to compute equil. biomass at F=0
  // --------------------------------------------------------
  dvar_vector M(1,ng);
  M=exp(fsh.get_nat_mort_species());
  n1(1)=1.0;
  for (int j=2;j<=ng-1;j++)  
  {
    n1(j)=n1(j-1)*exp(-M(j-1));
  }
  n1(ng)=n1(ng-1)*exp(-M(ng-1))/(1.0-exp(-M(ng)));
    
  dvariable b0;  // This is the adult biomass
  dvariable tb0; // This is for total biomass
  int rr=fsh.get_region_bounds()(1);
  if (numflag)
  {
    b0=pmature*n1;
    tb0=sum(n1); 
  }
  else
  {
    b0=pmature*(elem_prod(n1,fsh.mean_weight_yr(rr,1)))/1000.;
    tb0=sum(elem_prod(n1,fsh.mean_weight_yr(rr,1)))/1000.; 
  }
  dvariable n;              //NMD 21Dec2011
  if (fsh.age_flags(161))
  {
    n=alpha*mfexp(0.5*fsh.bh_variance)-beta/b0;
  }
  else
  {
    n=alpha-beta/b0;
  }                        //NMD 21Dec2011
// Adjustment for annualised BH-SRR   NMD_3July2015
  if ((fsh.age_flags(94)==3 || fsh.age_flags(182)) && fsh.age_flags(57)>1)
  {
    n/=fsh.age_flags(57);
  }
//
  tmp3(0)=b0*n;  // This scales up to the true equil. biomass
  tmp4(0)=tb0*n;
  // --------------------------------------------------------

  lambda_values(0) = 1.e-6;
  for (i=1;i<=nsteps;i++)
  {
    MY_DOUBLE_TYPE lambda=i/100.0;
    lambda_values(i) = lambda;
    //double lambda=i/2.0;
    
    MY_DOUBLE_TYPE mypen=0.0;
   
    dvar_vector Z;
    dvar_vector TF;
    if (fsh.age_flags(115)==0 && fsh.age_flags(92)!=3) // baranov
    {
      Z=lambda*F+exp(fsh.get_nat_mort_species());
    }
    else   // SS2
    {
      TF=lambda*F;
      for (int j=1;j<=fsh.nage;j++)  
      {
        if (TF(j)>0.9999) 
        {
          TF(j)=0.9999;
        }
      }

      Z=-log(1.0-TF)+exp(fsh.get_nat_mort_species());
    }
    int mmin=Z.indexmin();
    int mmax=Z.indexmax();
    for (int i1=mmin;i1<=mmax;i1++)
    {
#if !defined(NO_MY_DOUBLE_TYPE)
      Z(i1)=50.0-posfun(50.0-value(Z(i1)),1.0L,mypen);
#else
      Z(i1)=50.0-posfun(50.0-value(Z(i1)),1.0,mypen);
#endif
    }

    n1(1)=1.0;   // n1 is n-at-age per recruit (little phi in guide)
    for (int j=2;j<=ng-1;j++)  
    {
      n1(j)=n1(j-1)*exp(-Z(j-1));
    }
    n1(ng)=n1(ng-1)*exp(-Z(ng-1))/(1.0-exp(-Z(ng)));
    
    dvariable b1;   // Adult biomass
    dvariable tb1;  // This is for total biomass
    if (numflag)
    {
      b1=pmature*n1;
	    tb1=sum(n1);  // JH 21/03/02
    }
    else
    {
      b1=pmature*(elem_prod(n1,fsh.mean_weight_yr(rr,1)))/1000.;
      tb1=sum(elem_prod(n1,fsh.mean_weight_yr(rr,1)))/1000.; // JH 21/03/02
    }
    dvariable n; // n is recruitment
    if (fsh.age_flags(161))
    {
      n=alpha*mfexp(0.5*fsh.bh_variance)-beta/b1;
    }
    else
    {
      if (value(b1)>1.e-100)
        n=alpha-beta/b1;
      else
        n=alpha-beta*1.e+100;
    }
//  Adjustment for annualised BH-SRR   NMD_3July2015
    if ((fsh.age_flags(94)==3 || fsh.age_flags(182)) && fsh.age_flags(57)>1)
    {
      n/=fsh.age_flags(57);
    }
//
    dvar_vector C1;  // C1 is catch at age per recruit
    if (fsh.age_flags(115)==0 && fsh.age_flags(92)!=3) // baranov
    {
      C1=elem_prod(elem_div(lambda*F,Z),
        elem_prod(1.0-exp(-Z),n1));
    }
    else
    {
      C1=elem_prod(TF,elem_prod(exp(-0.5*exp(fsh.get_nat_mort_species())),n1));
    }

    dvariable Cn;
    dvariable TBn;
    dvariable ABn;
    int rrr=fsh.get_region_bounds()(1);
    if (numflag)
      Cn=n*sum(C1);
    else
      Cn=n*(C1*fsh.mean_weight_yr(rrr,1))/1000.;

    TBn=n*tb1;  // This is the scaled-up total biomass
    ABn=n*b1;   // and adult biomass

    if (Cn>0.0)  lastFmlt = ii; 
    if (fsh.age_flags(169)<0 && Cn<0.0)  
    {
      lastFmlt = ii-1; // for later output of just positive part
      break;
    }
                                       // yield curve PK 2003-06-10
    tmp1(ii)=lambda;         // deactivate truncation... PK 2003-06-10
    tmp3(ii)=ABn;
    tmp4(ii)=TBn;
    tmp2(ii++)=Cn;
  }
  get_biomass_ratio_by_time(fsh,numflag,pmature,tmp2,tmp3,tmp4);
  // JH 27/03/02 Calculates ratios of total and adult biomass (or numbers)
  //             to total and adult biomass (or numbers) at MSY
  // ---------------------------------------------------------------------
  dvar_vector abio(1,fsh.nyears);  // Adult biomass
  dvar_vector tbio(1,fsh.nyears);  // Total biomass
  dvar_vector tcat(1,fsh.nyears);  // Total catch
  abio.initialize();
  tbio.initialize();
  int rmin=fsh.get_region_bounds()(1);
  int rmax=fsh.get_region_bounds()(2);
  for (int ir=rmin;ir<=rmax;ir++)
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
    

  dvariable MSY=fsh.p_MSY+1.e-10;                // MSY
  fsh.MMSY=MSY;
  //double MSY=tmp2(imsy)+1.e-10;                // MSY

  dvar_vector F_agg(1,fsh.nyears);         // Aggregate F by year
  //double Fmsy=MSY/(fsh.p_Bmsy+1.e-10); // Aggregate F (exploitation rate) at MSY (F=catch/biomass)  
  fsh.p_Fmsy=value(MSY)/(fsh.p_Bmsy+1.e-10); // Aggregate F (exploitation rate) at MSY (F=catch/biomass)

  fsh.tb_ratio.deallocate();
  fsh.ab_ratio.deallocate();
  fsh.tb_ratio.allocate(1,fsh.nyears);
  fsh.ab_ratio.allocate(1,fsh.nyears);

  fsh.tb_ratio=tbio/(fsh.p_Bmsy+1.e-10); //  total bio to total biomass at MSY
  fsh.ab_ratio=abio/(fsh.p_sBmsy+1.e-10); //  adult bio to adult biomass at MSY
  fsh.F_ratio.deallocate();
  fsh.F_ratio.allocate(1,fsh.nyears);
  if (numflag)  // If yield calculations in number, calculate F_ratio in terms of
                // numbers
  {
    fsh.calculate_the_catch_numbers();
    F_agg=elem_div(fsh.catch_numbers,tbio+1.e-10);
    fsh.F_ratio=F_agg/fsh.p_Fmsy.to_double();
  }
  else  // If yield calculations in weight, calculate F_ratio in terms of weight
  {
    fsh.calculate_the_mean_weight();
    fsh.calculate_the_catch_biomass();
    F_agg=elem_div(fsh.catch_biomass,tbio+1.e-10);
    fsh.F_ratio=F_agg/fsh.p_Fmsy.to_double();
  }
  //cout<< "just before new stuff"<<endl; //xyxxx

  //-----------------------------------------------------------------------------------
  //NMD_19Dec2013 - set up in order of species before optimistion
  if (sflag)
  {
    FS.allocate(1,nage);
    FS.initialize();
    for (i=nomit+1;i<=navg;i++)
      FS+=fsh.pmsd->F_by_age_by_year(2,nyears-i+1);
    FS/=double(navg-nomit);
    psmature=fsh.get_pmature_species(2);
    F.initialize();
    for (i=nomit+1;i<=navg;i++)
      F+=fsh.F_by_age_by_year(nyears-i+1);
    F/=double(navg-nomit);
    pmature=fsh.pmature;
  }
  //NMD_19Dec2013

  //if (pof) 
  {
    dvar_vector eq_yield(0,nsteps);
    dvar_vector Beq(0,nsteps);
    dvar_vector sBeq(0,nsteps);
    Beq.initialize();
    sBeq.initialize();
    eq_yield.initialize();

    int lastFmult=nsteps;
    dvariable tt, stt;

    if (!sflag)
    {
      for (i=0;i<=nsteps;i++)
      {
        eq_yield(i)=get_yield_at_multiplier2(lambda_values(i),
          F,fsh.nyears,n1,pmature, alpha,beta,fsh,tt,stt);
        Beq(i)=tt;
        sBeq(i)=stt;
        if(eq_yield(i)<0.0) 
        {
          lastFmult=i-1;
          break;
        }     
      }
    }
    else
    {
      for (i=0;i<=nsteps;i++)
      {
        eq_yield(i)=get_yield_at_multiplier2(lambda_values(i),
          F,fsh.nyears,n1,pmature, alpha,beta,fsh,tt,stt,
          FS,psmature);
        Beq(i)=tt;
        sBeq(i)=stt;
        if(eq_yield(i)<0.0) 
        {
          lastFmult=i-1;
          break;
        }     
      }
    }

    fsh.predicted_yield_bh.deallocate();
    fsh.predicted_yield_bh.allocate(0,lastFmlt);  // lastFmlt replaces ii-1 -- PK 2003-06-10
    fsh.predicted_yield_bh=tmp2(0,lastFmlt);

    fsh.predicted_eqbio_bh.deallocate();        // JH 21/03/02 - This is in case we want
    fsh.predicted_eqbio_bh.allocate(0,lastFmlt);    // JH 21/03/02 - to include these in the
    fsh.predicted_eqbio_bh=tmp3(0,lastFmlt);        // JH 21/03/02 - SD report

    fsh.predicted_eqtotbio_bh.deallocate();     // JH 21/03/02
    fsh.predicted_eqtotbio_bh.allocate(0,lastFmlt); // JH 21/03/02
    fsh.predicted_eqtotbio_bh=tmp4(0,lastFmlt);     // JH 21/03/02

    fsh.predicted_yield_bh_x.deallocate();
    fsh.predicted_yield_bh_x.allocate(0,lastFmlt);
    fsh.predicted_yield_bh_x=tmp1(0,lastFmlt);

    {
      int cs=1;
      if (fsh.pmsd)
      {
        cs=fsh.pmsd->current_species;
      }
      dvariable alpha;
      dvariable beta;

      if (cs==1)
      {
        alpha=fsh.alpha;
        beta=fsh.beta;
      }
      else
      {
        alpha=fsh.pmsd->alpha(cs);
        beta=fsh.pmsd->beta(cs);
      }

      int sflag=0;    //NMD_18Dec2013
      if (fsh.pmsd)
      {
        sflag=sum(column(fsh.pmsd->species_flags,2));
      }
      if (sflag) fsh.pmsd->current_species=1;   //NMD_18Dec2013
   
      yield_analysis_bh_penalty(fsh,alpha,beta,pmature);
    }

    print_yield_stuff(pof,fsh,lastFmlt, F_agg, eq_yield,Beq,sBeq,lastFmult);
  }
  dvariable pen = 0.0;
  return pen;
} 

void print_yield_stuff(ofstream* pof,dvar_len_fish_stock_history& fsh,
                            int lastFmlt,  dvar_vector& F_agg, 
           dvar_vector& eq_yield, dvar_vector& Beq, dvar_vector& sBeq,int lastFmult)
{

  (*pof) << "# MSY" << endl;
  (*pof) << setw(10) << fsh.p_MSY << endl;
  //(*pof) << setw(10) << MSY << endl;

  (*pof) << "# F multiplier at MSY" << endl;
  (*pof) << setw(10) << fsh.p_Fmmsy << endl;
  //(*pof) << setw(10) << tmp1(imsy) << endl;

  (*pof) << "# F at MSY" << endl;

  if (fsh.p_Bmsy)
    (*pof) << setw(10) <<  fsh.p_MSY/fsh.p_Bmsy << endl;
  else
    (*pof) << setw(10) <<  "********"  << endl;

  //(*pof) << setw(10) <<  Fmsy << endl;

  (*pof) << "# Total biomass at MSY" << endl;
  (*pof) << setw(10) << fsh.p_Bmsy << endl;
  //(*pof) << setw(10) << tbio_at_msy << endl;

  (*pof) << "# Adult biomass at MSY" << endl;
  (*pof) << setw(10) << fsh.p_sBmsy  << endl;
  //(*pof) << setw(10) << abio_at_msy << endl;

  (*pof) << "# current Total Biomass to Total biomass at MSY" << endl;
  (*pof) << setw(10) << fsh.p_BBmsy << endl;

  (*pof) << "# current Adult Biomass to Adult Biomass at MSY" << endl;
  (*pof) << setw(10) << fsh.p_sBsBmsy << endl;

  (*pof) << "# Effort multiplier" << endl;
  (*pof) << setw(10) << fsh.predicted_yield_bh_x << endl;
  //(*pof) << setw(10) << tmp1(0,lastFmlt) << endl;

  (*pof) << "# Equilibrium yield" << endl;
  if (eq_yield.indexmin()==0)
    (*pof) << setw(10) << value(eq_yield(0,lastFmult)) << endl;
  else
    (*pof) << setw(10) << "  ********************  " << endl;
  //(*pof) << setw(10) << eq_yield << endl;
  //(*pof) << setw(10) << fsh.predicted_yield_bh << endl;
  //(*pof) << setw(10) << tmp2(0,lastFmlt) << endl;

  if (allocated(fsh.yield_eq_region))
  {
    (*pof) << "# Fmult values for Equilibrium yield by region" << endl;
    (*pof) << setw(10) << fsh.yield_eq_region_x(0,lastFmlt) << endl;
    (*pof) << "# Equilibrium yield by Fmult (across) and  region (down)" << endl;
    for (int i=1;i<=fsh.num_regions;i++)
    {
      (*pof) << setw(10) << column(fsh.yield_eq_region,i)(0,lastFmlt) << endl;
    }
  }

  (*pof) << "# Equilibrium adult biomass" << endl;
  if (sBeq.indexmin()==0)
    (*pof) << setw(10) << value(sBeq(0,lastFmult)) << endl;
  else
    (*pof) << setw(10) << "  ********************  " << endl;
  //(*pof) << setw(10) << fsh.predicted_eqbio_bh << endl;
  //(*pof) << setw(10) << tmp3(0,lastFmlt) << endl;

  (*pof) << "# Equilibrium total biomass" << endl;
  if (Beq.indexmin()==0)
    (*pof) << setw(10) << value(Beq(0,lastFmult)) << endl;
  else
    (*pof) << setw(10) << "  ********************  " << endl;

  //(*pof) << setw(10) << fsh.predicted_eqtotbio_bh << endl;
  //(*pof) << setw(10) << tmp4(0,lastFmlt) << endl;

  (*pof) << "# Adult biomass over adult biomass at MSY" << endl;
  (*pof) << setw(10) << fsh.ab_ratio << endl;

  (*pof) << "# Total biomass over total biomass at MSY" << endl;
  (*pof) << setw(10) << fsh.tb_ratio << endl;

  (*pof) << "# Aggregate F over F at MSY" << endl;
  (*pof) << setw(10) << fsh.F_ratio << endl;

  (*pof) << "# Aggregate F" << endl;
  (*pof) << setw(10) << F_agg << endl;
}

void get_biomass_ratio_by_time(dvar_len_fish_stock_history& fsh,
  int numflag,const dvar_vector& pmature,dvector& tmp2,dvector& tmp3,
  dvector& tmp4)
{
  // JH 27/03/02 Calculates ratios of total and adult biomass (or numbers)
  //             to total and adult biomass (or numbers) at MSY
  // ---------------------------------------------------------------------
  dvector abio(1,fsh.nyears);  // Adult biomass
  dvector tbio(1,fsh.nyears);  // Total biomass
  dvector tcat(1,fsh.nyears);  // Total catch
  abio.initialize();
  tbio.initialize();
  ivector rb=fsh.get_region_bounds();
  dvector cpmature=value(pmature);
  for (int ir=rb(1);ir<=rb(2);ir++)
  {
    for (int iy=1;iy<=fsh.nyears;iy++)  
    {
      // numflag determines if spawning stock is in numbers or weight
      if (numflag)
      {
        abio(iy)+=sum(elem_prod(cpmature,exp(value(fsh.N(ir,iy)))));
        tbio(iy)+=sum(exp(value(fsh.N(ir,iy))));
      }
      else
      {
        abio(iy)+=cpmature*(elem_prod(exp(value(fsh.N(ir,iy))),
          value(fsh.mean_weight_yr(ir,iy))))/1000.;
        tbio(iy)+=sum(elem_prod(exp(value(fsh.N(ir,iy))),
          value(fsh.mean_weight_yr(ir,iy))))/1000.;
      }
    }
  }

  int imsy=max_index(tmp2);                // gives the index of max(tmp2)
  if (imsy==0) imsy=1;   // can't use the 0 index
  MY_DOUBLE_TYPE tbio_at_msy=tmp4(imsy);        // equil. total biomass at MSY
  MY_DOUBLE_TYPE abio_at_msy=tmp3(imsy);        // equil. adult biomass at MSY
  MY_DOUBLE_TYPE MSY=tmp2(imsy)+1.e-10;                // MSY
  dvector F_agg(1,fsh.nyears);         // Aggregate F by year
  MY_DOUBLE_TYPE Fmsy=MSY/(tbio_at_msy+1.e-10); // Aggregate F at MSY (F=catch/biomass)
  fsh.tb_ratio.deallocate();
  fsh.ab_ratio.deallocate();
  fsh.tb_ratio.allocate(1,fsh.nyears);
  fsh.ab_ratio.allocate(1,fsh.nyears);
  fsh.tb_ratio=tbio/(tbio_at_msy+1.e-10); //  total bio to total biomass at MSY
  fsh.ab_ratio=abio/(abio_at_msy+1.e-10); //  adult bio to adult biomass at MSY
  fsh.F_ratio.deallocate();
  fsh.F_ratio.allocate(1,fsh.nyears);
  if (numflag)  // If yield calculations in number, calculate F_ratio in terms of
                // numbers
  {
    fsh.calculate_the_catch_numbers();
    F_agg=elem_div(value(fsh.catch_numbers),tbio+1.e-10);
    fsh.F_ratio=F_agg/Fmsy;
  }
  else  // If yield calculations in weight, calculate F_ratio in terms of weight
  {
    fsh.calculate_the_mean_weight();
    fsh.calculate_the_catch_biomass();
    F_agg=elem_div(value(fsh.catch_biomass),tbio+1.e-10);
    fsh.F_ratio=F_agg/Fmsy;
  }
  //cout<< "just before new stuff"<<endl; //xyxxx
}


void get_biomass_ratio_by_time(dvar_len_fish_stock_history& fsh,
  int numflag,const dvar_vector& pmature,dvar_vector& tmp2,dvar_vector& tmp3,
  const dvar_vector& tmp4)
{
  // JH 27/03/02 Calculates ratios of total and adult biomass (or numbers)
  //             to total and adult biomass (or numbers) at MSY
  // ---------------------------------------------------------------------
  dvar_vector abio(1,fsh.nyears);  // Adult biomass
  dvar_vector tbio(1,fsh.nyears);  // Total biomass
  dvar_vector tcat(1,fsh.nyears);  // Total catch
  abio.initialize();
  tbio.initialize();
  int rmin=fsh.get_region_bounds()(1);
  int rmax=fsh.get_region_bounds()(2);
  for (int ir=rmin;ir<=rmax;ir++)
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

  int imsy=max_index(tmp2);                // gives the index of max(tmp2)
  if (imsy==0) imsy=1;   // can't use the 0 index
  dvariable tbio_at_msy=tmp4(imsy);        // equil. total biomass at MSY
  dvariable abio_at_msy=tmp3(imsy);        // equil. adult biomass at MSY
  dvariable MSY=tmp2(imsy)+1.e-10;                // MSY
  dvar_vector F_agg(1,fsh.nyears);         // Aggregate F by year
  dvariable Fmsy=MSY/(tbio_at_msy+1.e-10); // Aggregate F at MSY (F=catch/biomass)
  fsh.tb_ratio.deallocate();
  fsh.ab_ratio.deallocate();
  fsh.tb_ratio.allocate(1,fsh.nyears);
  fsh.ab_ratio.allocate(1,fsh.nyears);
  fsh.tb_ratio=tbio/(tbio_at_msy+1.e-10); //  total bio to total biomass at MSY
  fsh.ab_ratio=abio/(abio_at_msy+1.e-10); //  adult bio to adult biomass at MSY
  fsh.F_ratio.deallocate();
  fsh.F_ratio.allocate(1,fsh.nyears);
  if (numflag)  // If yield calculations in number, calculate F_ratio in terms of
                // numbers
  {
    fsh.calculate_the_catch_numbers();
    F_agg=elem_div(fsh.catch_numbers,tbio+1.e-10);
    fsh.F_ratio=F_agg/Fmsy;
  }
  else  // If yield calculations in weight, calculate F_ratio in terms of weight
  {
    fsh.calculate_the_mean_weight();
    fsh.calculate_the_catch_biomass();
    F_agg=elem_div(fsh.catch_biomass,tbio+1.e-10);
    fsh.F_ratio=F_agg/Fmsy;
  }
  //cout<< "just before new stuff"<<endl; //xyxxx
}

//*************************************************
//*************************************************
//*************************************************
dvariable get_yield_at_multiplier2(MY_DOUBLE_TYPE lambda, dvar_vector& _F,int nyears,
  dvar_vector& n1, const dvar_vector& _pmature,dvariable& alpha,dvariable& beta,
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
    bsa=psmature*ns1;
  }
  else
  {
    bt=n1*fsh.mean_weight_yr(rr,1)/1000.;
    ba=pmature*(elem_prod(n1,fsh.mean_weight_yr(rr,1)))/1000.;
    bst=ns1*fsh.mean_weight_yr(rs,1)/1000.;
//    ba=psmature*(elem_prod(ns1,fsh.mean_weight_yr(rs,1)))/1000.;
    bsa=psmature*(elem_prod(ns1,fsh.mean_weight_yr(rs,1)))/1000.; //NMD_23jun-17
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

