/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
//% message to emacs: -*- mode: C++; fill-column: 80 -*-


#include "all.hpp"

dvariable yield_analysis_bh_penalty(dvar_len_fish_stock_history& fsh,
 dvariable alpha,dvariable beta,dvector& pmature);
 //dvariable alpha,dvariable beta,dvector& pmature,dvar_matrix& Fay);
//dvariable alpha,dvariable beta,MY_DOUBLE_TYPE lwc,dvector& pmature,dvar_matrix& Fay);


dvariable calculate_the_biomass(int i,dvar_fish_stock_history& fsh,
  dvar3_array & N,int numflag,const dvar_vector& _pmature)
{
  ADUNCONST(dvar_vector,pmature)
  dvariable bio=0.0;
  dvar_vector * pm=0;
  int mmin,mmax;
  if (fsh.pmsd==0) 
  {
    mmin=1;
    mmax=fsh.num_regions;
    pm=&pmature;
  } 
  else 
  {
    int nrr=fsh.pmsd->num_real_regions;
    int cs=fsh.pmsd->current_species;
    mmin=fsh.pmsd->region_bounds(cs,1);
    mmax=fsh.pmsd->region_bounds(cs,2);
    if (cs>1)
      pm=&(fsh.pmsd->pmature(cs));
    else
      pm=&(pmature);
  }
  
  for (int ir=mmin;ir<=mmax;ir++) {
    if (numflag) {
      bio+=sum(elem_prod(*pm,exp(N(ir,i))));
    } else {
      bio+=*pm*(elem_prod(exp(N(ir,i)),
        fsh.mean_weight_yr(ir,i)))/1000.;
    }
  }
  return bio;
}


dvariable stock_recruit_bh(dvar_len_fish_stock_history& fsh,
 ofstream * pof,dvar_matrix& Fay,ivector * pq_flag)
{
  fsh.bh_numbers.initialize();
  if (fsh.seasonal_recruitment_flag)
  {
    cerr << "seasonal recruitment in stock-recruit_bh not"
            " implemented yet" << endl;
    ad_exit(1);
  }
  int lag=fsh.age_flags(147);
  int numflag=fsh.age_flags(149);
  dvar_vector rec(1,fsh.last_real_year);
  dvar_vector bio(1,fsh.last_real_year);
  bio.initialize();
  rec.initialize();
  dvar_vector vvar=fsh.global_vars;
  MY_DOUBLE_TYPE lwc=fsh.len_wt_coff;
  if(lwc==0) lwc=1.0;
  dvar_vector& pmature = fsh.pmature;
  if (!value(sum(pmature)))
  {
    for (int j=1;j<=fsh.nage;j++)
    {
      pmature=1.0;
    }
  } 
  //dvariable N1=0.0;
  dvar3_array * pN=0;
  if (fsh.af170q0==0)
  {
    pN=&(fsh.N);
  }
  else
  {
    pN=&(fsh.N_q0);
  }
  for (int i=1;i<=fsh.last_real_year;i++)  
  {
    // numflag determines if spawning stock is in numbers or weight
    bio(i)+= calculate_the_biomass(i,fsh,*pN,numflag,pmature);

    for (int ir=1;ir<=fsh.num_regions;ir++)
    {
      rec(i)+=fsh.exp_N(ir,i,1);
      fsh.bh_numbers(i)+=sum(exp(value((*pN)(ir,i))));
    }
  }
  fsh.bh_reproductive_biomass(1,fsh.last_real_year)=value(bio);
  dvar_vector b=bio(1,fsh.last_real_year-lag);
  fsh.bh_bio.initialize();
  fsh.bh_bio(1,fsh.last_real_year-lag)=value(b);
  //cout <<"vrbioclc.cpp " << b(31) << endl;
  dvar_vector r=rec(1+lag,fsh.last_real_year).shift(1);
  dvar_vector lr=log(r);
  dvar_vector lb=log(b)-log((fsh.sv(21)+0.001)*b(1)+b);
  dvariable mlr=mean(lr);
  dvariable mlb=mean(lb);
  dvariable loga=mlr-mlb+fsh.sv(26);
//  MY_DOUBLE_TYPE pen_wt=fsh.age_flags(145);
//  pen_wt = pen_wt/1000.0;  //NMD 31Jan2013 testing
  MY_DOUBLE_TYPE pen_wt;
  if (fsh.age_flags(145)>0) pen_wt=fsh.age_flags(145);    //NMD_23Dec2013
  if (fsh.age_flags(145)<0) pen_wt=pow(10.0,fsh.age_flags(145));   //NMD_23Dec2013
  dvariable f;
  fsh.bh_predicted_recruits(1,fsh.last_real_year-lag)=exp(value(loga)+value(lb));
  dvar_vector  tmpdif=lr-loga-lb;
  fsh.bh_recr_devs(1,fsh.last_real_year-lag)=tmpdif;
  dvariable tmp=norm2(tmpdif);
  f=pen_wt*tmp;
  if (!fsh.af170q0 && !fsh.do_fishery_projections_flag //NMD_11Aug2022
      && !fsh.projection_sim_index)
  {
    fsh.bh_variance=tmp/(fsh.last_real_year-lag);
  }


   dvariable fpen=0.0;
   dvariable alpha=0.0;
   dvariable beta=0.0;
   dvariable steepness=0.0;  
   dvariable phi=0.0;  
    if (!fsh.zeroed_fisheries_flag)
    {
      alpha=exp(loga);
      beta=b(1)*(fsh.sv(21)+0.001);
      fsh.alpha=alpha;
      fsh.beta=beta;
      dvar_vector tmp=elem_div(alpha*b,beta+b);
      dvariable recpen=0.01*norm2(r-tmp);
      // calculate phi with alpha, beta, and mortality and
      phi=0.0;
      dvariable n0=1.0;
      dvar_vector M=exp(fsh.nat_mort(1));
      for (int j=1;j<=fsh.nage-1;j++)
      {
        if (numflag)          // spawning stock in numbers
          phi+=n0*pmature(j);
        else                  // spawning stock in weight
          phi+=n0*pmature(j)*fsh.mean_weight_yr(1,1,j)/1000.;
        n0*=exp(-M(j));
      }
      if (numflag)          // spawning stock in numbers
        phi+=n0/(1.000001-exp(-M(fsh.nage)))*pmature(fsh.nage);
      else                  // spawning stock in weight
        phi+=n0/(1.000001-exp(-M(fsh.nage)))*pmature(fsh.nage)*
             fsh.mean_weight_yr(1,1,fsh.nage)/1000.;
    
      // now use phi to finish the steepness calculation
      // steepeness is the ratio of equil. rec. at 20% unexploited spawning stock
      // to equil. rec. at the unexploited spawning stock level
      fsh.phi=phi;
      steepness=alpha*phi/(4.0*beta+alpha*phi);  
      steepness=0.2+posfun(steepness-0.2,.02,fpen);
      fsh.steepness=steepness;
      fsh.sv(29)=steepness;
    }
    else
    {
      alpha=fsh.alpha;
      beta=fsh.beta;
      phi=fsh.phi;
      steepness=fsh.steepness;
    }
   { 
    //cout << "&alpha = " << &alpha << endl;
    //cout << "alpha = " << alpha << endl;
    //cout << "beta  = " << beta  << endl;
    cout << "phi   = " << phi   << endl;
    // This is for a lognormal prior on the steepness
    /*
    if (fsh.age_flags(153)>0)
    {
      MY_DOUBLE_TYPE pen_wt_steepness=fsh.age_flags(153);
      MY_DOUBLE_TYPE steepness_prior=double(fsh.age_flags(154))/10.;
      dvariable steep_pen=pen_wt_steepness*square(log(steepness/steepness_prior));
      f+=steep_pen;
      cout << "steep = " << steepness << " Penalty = " << steep_pen << endl;
    }
    */
    // This is for a beta prior on the steepness
    // The PDF is f(x)=(x-lb)^(a-1)*(ub-x)^(b-1)/(B(a,b)*(ub-lb)^(a+b-1)
    // where a,b are parameters of the beta dist, lb is the lower bound,
    // ub is the upper bound, and B(a,b) is the integral of the beta function
    // evaluated over the range 0-1.
    if (fsh.age_flags(153)>0)
    {
      MY_DOUBLE_TYPE beta_a=double(fsh.age_flags(153))/10.; // a parameter for beta dist.
      MY_DOUBLE_TYPE beta_b=double(fsh.age_flags(154))/10.; // b parameter for beta dist.
      MY_DOUBLE_TYPE lb=0.2;
      MY_DOUBLE_TYPE ub=1.0;
      MY_DOUBLE_TYPE bnorm=0.0;
      MY_DOUBLE_TYPE x;
   // This evaluates the integral of Beta(a,b) from 0 to 1
      for (int dt=0;dt<=100;dt++)
      {
        x=dt/100.;
        bnorm+=pow(x,beta_a-1)*pow(1-x,beta_b-1);
      }
      MY_DOUBLE_TYPE beta_mode=(beta_a-1+lb)/(beta_a+beta_b-2); // prior mode
      MY_DOUBLE_TYPE beta_sd=sqrt(beta_a*beta_b/(beta_a+beta_b)/(beta_a+beta_b)/
                     (beta_a+beta_b+1))*(ub-lb);       // priod sd
    // Below is -log(beta(beta_a,beta_b,lb,ub))
#if !defined(NO_MY_DOUBLE_TYPE)
      dvariable steep_pen=-(beta_a-1.0L)*log(steepness-lb)-
#else
      dvariable steep_pen=-(beta_a-1.0)*log(steepness-lb)-
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
                          (beta_b-1.0L)*log(ub-steepness)+
#else
                          (beta_b-1.0)*log(ub-steepness)+
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
                          log(bnorm)+(beta_a+beta_b-1.0L)*log(ub-lb);
#else
                          log(bnorm)+(beta_a+beta_b-1.0)*log(ub-lb);
#endif
      f+=steep_pen;
      fsh.steep_penalty=value(steep_pen);
      cout << "steepness =  " << steepness << " Penalty =  " << steep_pen << endl;
      cout << "prior mode = " << beta_mode << " prior sd = " << beta_sd << endl;
    }
  
    int ii=0;
    const int nsteps=2000;
    dvector tmp1(1,nsteps);
    dvar_vector tmp2(1,nsteps);
  
    if (pof)
    {
      (*pof) << "# Beverton-Holt stock-recruitment relationship report" << endl;
      (*pof) << "# alpha = " << alpha << "  beta =  " << beta 
             << "  steepness = " << steepness << endl;
      (*pof) << "# Observed spawning Biomass" << endl;
      (*pof) << b << endl;
      (*pof) << "# Observed recruitment" << endl;
      (*pof) << r << endl;
      MY_DOUBLE_TYPE maxb=1.2*max(value(b));
      MY_DOUBLE_TYPE diff=maxb/100.0;
      MY_DOUBLE_TYPE fi=0; 
      do
      {
         if (ii==nsteps)
         {
           cerr << "error in beverton Holt analysis -- need to increase"
                   " nsteps " << endl;
           break;
         }
         tmp1(++ii)=fi;
         tmp2(ii)=alpha*fi/(beta+fi);
         fi+=diff;
      }
      while (fi<0.999*maxb);
      fsh.predicted_recruitment_bh.deallocate();
      fsh.predicted_recruitment_bh.allocate(1,ii);
      fsh.predicted_recruitment_bh=tmp2(1,ii);
      fsh.predicted_recruitment_bh_x.deallocate();
      fsh.predicted_recruitment_bh_x.allocate(1,ii);
      fsh.predicted_recruitment_bh_x=tmp1(1,ii);
      
      (*pof) << "# Spawning Biomass" << endl;
      (*pof) << tmp1(1,ii) << endl;
      (*pof) << "# Predicted recruitment" << endl;
      (*pof) << tmp2(1,ii) << endl;
      //yield_analysis_bh(fsh,pof,alpha,beta,pmature,Fay); //xyxxx
    }
  }
   //cout<< "just before call to yield_analysis_bh"<<endl; //xyxxx
  // DF  03 DEC 03
  if ( (pof) || fsh.age_flags(165)>0)
  {
    if (fsh.age_flags(140)==0) // af140 is for region-specific yield
    {
        if (fsh.age_flags(165)>0)  // Don't do this if only here for report
        {                          // generation.
          f += yield_analysis_bh_penalty(fsh,alpha,beta,pmature);
        }
        if( pof ) {
          if (fsh.age_flags(165)==0) // do once if not done above
               yield_analysis_bh_penalty(fsh,alpha,beta,pmature);
          // now print to report file
          yield_analysis_bh(fsh,pof,alpha,beta,pmature);
        }
    }
    else
    {
      if (pof)  { // --- this does region-specific yields in report (plot.rep)
        //cout<<"ABOUT TO ENTER folly IN vrbioclc:stock_recruit_bh"<<endl;
        //yield_analysis_bh(fsh,0,alpha,beta,pmature,Fay);
        yield_analysis_bh_daves_folly(fsh,pof,alpha,beta,lwc,
          pmature);
      }
    }
  }
  return f+1000.0*fpen;
} 

int max_index(const dvar_vector & v) 
// Returns the index of v that corresponds to max(v) (JH 27/03/02)
{ 
  int mmin=v.indexmin(); 
  int mmax=v.indexmax(); 
  int imax=mmin; 
  MY_DOUBLE_TYPE dmax=value(v(mmin)); 
  for (int i=mmin;i<=mmax;i++) 
  { 
    if (dmax<value(v(i))) 
    { 
      dmax=value(v(i)); 
      imax=i; 
    } 
  } 
  return imax; 
}

int max_index(const dvector & v) 
// Returns the index of v that corresponds to max(v) (JH 27/03/02)
{ 
  int mmin=v.indexmin(); 
  int mmax=v.indexmax(); 
  int imax=mmin; 
  MY_DOUBLE_TYPE dmax=v(mmin); 
  for (int i=mmin;i<=mmax;i++) 
  { 
    if (dmax<v(i)) 
    { 
      dmax=v(i); 
      imax=i; 
    } 
  } 
  return imax; 
}

int min_index(const dvar_vector & v) 
// Returns the index of v that corresponds to max(v) (JH 27/03/02)
{ 
  int mmin=v.indexmin(); 
  int mmax=v.indexmax(); 
  int imin=mmin; 
  MY_DOUBLE_TYPE dmin=value(v(mmin)); 
  for (int i=mmin;i<=mmax;i++) 
  { 
    if (dmin>value(v(i))) 
    { 
      dmin=value(v(i)); 
      imin=i; 
    } 
  } 
  return imin; 
}

int min_index(const dvector & v) 
// Returns the index of v that corresponds to max(v) (JH 27/03/02)
{ 
  int mmin=v.indexmin(); 
  int mmax=v.indexmax(); 
  int imin=mmin; 
  MY_DOUBLE_TYPE dmin=v(mmin); 
  for (int i=mmin;i<=mmax;i++) 
  { 
    if (dmin>v(i)) 
    { 
      dmin=v(i); 
      imin=i; 
    } 
  } 
  return imin; 
}

void yield_per_recruit_analysis(dvar_len_fish_stock_history& fsh,
 ofstream * pof,dvar_matrix& Fay)
{
  fsh.pmsd_error();
  int nage=fsh.nage;
  int numflag=fsh.age_flags(149);
  int nyears=fsh.nyears;
  dvar_vector n1(1,nage);
  (*pof) << "# Yield per recruit report" << endl;
  const int nsteps=2000;
  dvector tmp1(0,nsteps);
  dvar_vector tmp2(0,nsteps);
  int i;
  int ii=1;
  tmp1(0)=0.0;
  tmp2(0)=0.0;
  int navg=fsh.age_flags(148);
  int tmult=fsh.age_flags(57);
  if (!navg)navg=tmult;
  if (!navg)navg=1;

  dvar_vector F(1,nage);
  F.initialize();
  for (i=1;i<=navg;i++)
    F+=Fay(fsh.last_real_year-i+1);
  F/=double(navg);

  for (i=1;i<=nsteps;i++)
  {
    MY_DOUBLE_TYPE lambda=i/10.0;
    dvar_vector Z;
    dvar_vector TF;
    if (fsh.age_flags(115)==0) // baranov
    {
      Z=lambda*F+exp(fsh.nat_mort(nyears));
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

      Z=-log(1.0-TF)+exp(fsh.nat_mort(nyears));
    }
    n1(1)=1.0;
    for (int j=2;j<=fsh.nage-1;j++)  
    {
      n1(j)=n1(j-1)*exp(-Z(j-1));
    }
    n1(nage)=n1(nage-1)*exp(-Z(nage-1))/(1.0-exp(-Z(nage)));
    
    dvar_vector C1=elem_prod(elem_div(lambda*F,Z),
       elem_prod(1.0-exp(-Z),n1));

    dvariable Ypr;
    if (numflag)
      Ypr=sum(C1);
    else
      Ypr=(C1*fsh.mean_weight_yr(1,1))/1000.;

    if (value(Ypr)>0.0)
    {
      tmp1(ii)=lambda;
      tmp2(ii++)=Ypr;
    } 
    else 
    {
      tmp1(ii)=lambda;
      tmp2(ii++)=0.0;
      break;
    }
  }
  (*pof) << "# Effort multiplier" << endl;
  (*pof) << setw(10) << tmp1(0,ii-1) << endl;
  (*pof) << "# Yield per recruit" << endl;
  (*pof) << setw(10) << tmp2(0,ii-1) << endl;
} 


dvar_vector sbcalc(dvar_len_fish_stock_history& fsh)
{
  dvar_vector bio(1,fsh.last_real_year);
  bio.initialize();
  dvar_vector vvar=value(fsh.global_vars);
  MY_DOUBLE_TYPE lwc=fsh.len_wt_coff;
  if(lwc==0) lwc=1.0;
  dvar_vector & pmature = fsh.pmature;
  if (!value(sum(pmature)))
  {
    for (int j=1;j<=fsh.nage;j++)
    {
      pmature=1.0;
    }
  } 
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int i=1;i<=fsh.last_real_year;i++)  
    {
      bio(i)+=pmature*elem_prod(exp(fsh.N(ir,i)),
        fsh.mean_weight_yr(ir,i))/1000.;
    }
  }
  return bio;
} 

dvariable ratio_first_last_biocalc
  (dvar_len_fish_stock_history& fsh)
{
  dvar_vector bio(1,2);
  bio.initialize();
  int nav=1;
  if (fsh.age_flags(99)>0)nav=fsh.age_flags(99);
  int i;
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    for (i=1;i<=nav;i++)
    {
      bio(1)+=exp(fsh.N(ir,i))*fsh.mean_weight_yr(ir,i);
    }
    for (i=fsh.last_real_year-nav+1;i<=fsh.last_real_year;i++)
    {
      bio(2)+=exp(fsh.N(ir,i))*fsh.mean_weight_yr(ir,i);
    }
  }
  return bio(2)/bio(1);
} 

dvar_vector vbiocalc(dvar_len_fish_stock_history& fsh,int q0flag)
{
  dvar3_array * pN=0;
  if (q0flag==0)
  {
    pN=&fsh.N;
  }
  else
  {
    pN=&fsh.N_q0;
  }
  dvar_vector bio(1,fsh.last_real_year);
  bio.initialize();
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int i=1;i<=fsh.last_real_year;i++)  
    {
      bio(i)+=exp((*pN)(ir,i))*fsh.mean_weight_yr(ir,i);
    }
  }
  return bio;
} 

dvar_vector adult_vbiocalc(dvar_len_fish_stock_history& fsh,int q0flag)
{
  dvar_vector bio(1,fsh.last_real_year);
  bio.initialize();
  dvar3_array * pN=0;
  if (q0flag==0)
  {
    pN=&fsh.N;
  }
  else
  {
    pN=&fsh.N_q0;
  }
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int i=1;i<=fsh.last_real_year;i++)  
    {
      bio(i)+=elem_prod(fsh.pmature,exp((*pN)(ir,i)))
        *fsh.mean_weight_yr(ir,i);
    }
  }
  return bio;
} 

dvar_matrix reg_vbiocalc(dvar_len_fish_stock_history& fsh,int q0flag)
{
  dvar3_array * pN=0;
  if (q0flag==0)
  {
    pN=&fsh.N;
  }
  else
  {
    pN=&fsh.N_q0;
  }
  dvar_matrix bio(1,fsh.num_regions,1,fsh.last_real_year);
  bio.initialize();
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int i=1;i<=fsh.last_real_year;i++)  
    {
      bio(ir,i)=exp((*pN)(ir,i))*fsh.mean_weight_yr(ir,i);
    }
  }
  dvariable avgbio=sum(bio)/size_count(bio);
  bio=bio/avgbio;
  return bio;
} 

dvar_matrix unnormalized_reg_vbiocalc(dvar_len_fish_stock_history& fsh,
  int q0flag)
{
  dvar3_array * pN=0;
  if (q0flag==0)
  {
    pN=&fsh.N;
  }
  else
  {
    pN=&fsh.N_q0;
  }
  dvar_matrix bio(1,fsh.num_regions,1,fsh.last_real_year);
  bio.initialize();
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int i=1;i<=fsh.last_real_year;i++)  
    {
      bio(ir,i)=exp((*pN)(ir,i))*fsh.mean_weight_yr(ir,i);
    }
  }
  return bio;
} 

dvar_matrix unnormalized_adult_reg_vbiocalc(dvar_len_fish_stock_history& fsh,
  int q0flag)
{
  dvar3_array * pN=0;
  if (q0flag==0)
  {
    pN=&fsh.N;
  }
  else
  {
    pN=&fsh.N_q0;
  }
  dvar_matrix bio(1,fsh.num_regions,1,fsh.last_real_year);
  bio.initialize();
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int i=1;i<=fsh.last_real_year;i++)  
    {
      bio(ir,i)=elem_prod(fsh.pmature,exp((*pN)(ir,i)))
        *fsh.mean_weight_yr(ir,i);
    }
  }
  return bio;
} 

