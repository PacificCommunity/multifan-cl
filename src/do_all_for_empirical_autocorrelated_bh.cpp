/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
//% message to emacs: -*- mode: C++; fill-column: 80 -*-


#include "all.hpp"

//dvariable yield_analysis_bh(dvar_len_fish_stock_history& fsh,
// ofstream * pof,dvariable alpha,dvariable beta,
// dvector& pmature);

extern int debug_flag;


dvariable stock_recruit_bh_steep(dvar_len_fish_stock_history& fsh,
 ofstream * pof,ivector * pq_flag)
{
  int cs=1;
  if (fsh.pmsd)
    cs=fsh.pmsd->current_species;
  dvar_vector * pm=0;
  int mmin,mmax;
  if (fsh.pmsd==0)
  {
    mmin=1;
    mmax=fsh.num_regions;
    pm=&fsh.pmature;
  }
  else
  {
    int nrr=fsh.pmsd->num_real_regions;
    cs=fsh.pmsd->current_species;
    mmin=fsh.pmsd->region_bounds(cs,1);
    mmax=fsh.pmsd->region_bounds(cs,2);
    if (cs>1)
      pm=&(fsh.pmsd->pmature(cs));
    else
      pm=&(fsh.pmature);
  }
  dvar_vector & pmature = *pm;

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

  int ns=fsh.age_flags(57);
  int rem=fsh.last_real_year%ns;
  int na=fsh.last_real_year/ns;
  if (rem) na++;
  int na1=(fsh.last_real_year-1)/fsh.age_flags(57)+1;
  if (na != na1)
  {
    cout << "ERROR: stock-recruit_bh annual number of years" << endl;
    cout << "check setting of age_flags(57) - exiting " << endl;
    ad_exit(1);
  }
    
  bio.initialize();
  rec.initialize();
  dvar_vector vvar;
  MY_DOUBLE_TYPE lwc;
  if (cs==1)
  {
    vvar=fsh.global_vars;
    lwc=fsh.len_wt_coff;
  }
  else
  {
    vvar=fsh.pmsd->global_vars(cs);
    cout << "vvar(1) = " << vvar(1) << endl;
    lwc=fsh.pmsd->len_wt_coff(cs);
  }
  if(lwc==0) lwc=1.0;
  if (!value(sum(pmature)))
  {
    for (int j=1;j<=fsh.nage;j++)
    {
      pmature=1.0;
    }
  } 
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

    for (int ir=mmin;ir<=mmax;ir++)
    {
      rec(i)+=fsh.exp_N(ir,i,1);
      fsh.bh_numbers(i)+=sum(exp(value((*pN)(ir,i))));
    }
  }
  dvar_vector annual_rec(1,na);
  dvar_vector annual_bio(1,na);
  annual_rec.initialize();
  annual_bio.initialize();
  dvariable equilibrium_annual_rec=0;
 
  int ii=0;
  int i;
  for (i=1;i<=na;i++)
  {
    int ic=0;
    for (int j=1;j<=ns;j++)
    {
      ++ii;
      if (ii+lag>fsh.last_real_year) break;
      ic++;
      annual_rec(i)+=rec(ii+lag);
      if (i==1)
      {
        equilibrium_annual_rec+=rec(ii);
      }
      annual_bio(i)+=bio(ii);
    }
    if (ic>1)
      annual_bio(i)/=ic;
  }
  fsh.annual_phi=annual_bio(1)/equilibrium_annual_rec;

  cout << "phi = " << fsh.annual_phi << endl;
  if (fsh.age_flags(94)==3  || fsh.age_flags(182) )
  {
    fsh.bh_reproductive_biomass=value(annual_bio);
  }
  else
  {
    fsh.bh_reproductive_biomass(1,fsh.last_real_year)=value(bio);
  }
  dvar_vector b;
  if (fsh.age_flags(94)==3  || fsh.age_flags(182) )
  {
    b=annual_bio;
  }
  else
  {
    b=bio(1,fsh.last_real_year-lag);
  }
  fsh.bh_bio.initialize();
  if (fsh.age_flags(94)==3  || fsh.age_flags(182) )
  {
    fsh.bh_bio=value(b);
  }
  else
  {
    fsh.bh_bio(1,fsh.last_real_year-lag)=value(b);
  }
  dvar_vector r;
  if (fsh.age_flags(94)==3  || fsh.age_flags(182) )
  {
    r=annual_rec;
  }
  else
  {
    r=rec(1+lag,fsh.last_real_year).shift(1);
  }
  dvar_vector lr=log(r);
  dvariable sv21;
  if (cs==1)
  {
    sv21=fsh.sv(21);
  }
  else
  {
    sv21=fsh.pmsd->sv(cs,21);
  }
  dvariable phi=0.0;  
  // **********************************************************
  // **********************************************************
  if (!fsh.zeroed_fisheries_flag)
  {
    // calculate phi with alpha, beta, and mortality and
    phi=0.0;
    dvariable n0=1.0;
    dvar_vector M=exp(fsh.get_nat_mort_species(cs)(1));
    int rr=1;
    if (fsh.pmsd && fsh.pmsd->current_species>1)
      rr=fsh.pmsd->region_bounds(cs,1);
                 

    for (int j=1;j<=fsh.nage-1;j++)
    {
      if (numflag)          // spawning stock in numbers
        phi+=n0*pmature(j);
      else                  // spawning stock in weight
        phi+=n0*pmature(j)*fsh.mean_weight_yr(rr,1,j)/1000.;
      n0*=exp(-M(j));
    }
    if (numflag)          // spawning stock in numbers
      phi+=n0/(1.000001-exp(-M(fsh.nage)))*pmature(fsh.nage);
    else                  // spawning stock in weight
      phi+=n0/(1.000001-exp(-M(fsh.nage)))*pmature(fsh.nage)*
           fsh.mean_weight_yr(rr,1,fsh.nage)/1000.;
    // if using annual stock recruitment relationship
    // divide by the number of seasons because we want the
    // equilibrium biomass corresponding to a constant recruitment of
    // one fish per year and above is one fish per season.
    if (fsh.age_flags(182))
    {
      phi/=fsh.age_flags(57);
    }
  
    // now use phi to finish the steepness calculation
    // steepeness is the ratio of equil. rec. at 20% unexploited spawning stock
    // to equil. rec. at the unexploited spawning stock level
    if (cs==1)
    {
      fsh.phi=phi;
    }
    else
    {
      fsh.pmsd->phi(cs)=phi;
    }
  }
  else
  {
    // get phi calculated before effort was turned off
    // DF 25/04/08  !!!!!!!!!! check the logic
    if (cs==1)
    {
      phi=fsh.phi;
    }
    else
    {
      phi=fsh.pmsd->phi(cs);
    }

  }
  // **********************************************************
  // **********************************************************
  dvariable steepness=0.0;
  dvariable beta=0.0;
  dvariable alpha=0.0;
  dvariable loga=0.0;
  dvar_vector lb;

  int old_param=1;
  if (fsh.age_flags(94)==3)
    old_param=0;

  if (old_param)
  {
    if (cs==1)
    {
      steepness=fsh.sv(29);
      beta=b(1)*(fsh.sv(21)+0.001);
    }
    else
    {
      steepness=fsh.pmsd->sv(cs,29);
      beta=b(1)*(fsh.pmsd->sv(cs,21)+0.001);
    }
  
    dvariable logdelta=log(4.0*steepness/(phi*(1.0-steepness)));
    loga=logdelta+log(beta);
    alpha=exp(loga);
    lb=log(b)-log(beta+b);
  }
  else
  {
    if (cs==1)
    {
      steepness=fsh.sv(29);
      alpha=4.0*steepness*r(1)/(5.0*steepness-1);
#if !defined(NO_MY_DOUBLE_TYPE)
      beta=(r(1)*fsh.annual_phi*(1.0-steepness))/(5.0*steepness-1.0L);
#else
      beta=(r(1)*fsh.annual_phi*(1.0-steepness))/(5.0*steepness-1.0);
#endif
    }
    else
    {
      steepness=fsh.pmsd->sv(cs,29);
      alpha=4.0*steepness*r(1)/(5.0*steepness-1);
#if !defined(NO_MY_DOUBLE_TYPE)
      beta=(5.0*steepness-1.0L)/(r(1)*fsh.annual_phi*(1.0-steepness));
#else
      beta=(5.0*steepness-1.0)/(r(1)*fsh.annual_phi*(1.0-steepness));
#endif
    }
    loga=log(alpha);
    lb=log(b)-log(beta+b);
  }
  MY_DOUBLE_TYPE pen_wt=0.0;
  if (fsh.age_flags(145)>0) pen_wt=fsh.age_flags(145);    //NMD_23Dec2013
  if (fsh.age_flags(145)<0) pen_wt=pow(10.0,fsh.age_flags(145));   //NMD_23Dec2013
  dvariable f;
  dvar_vector  tmpdif=lr-loga-lb;
  // *********************************************************
  // *********************************************************
  // test stuff for using gamma regression for beverton holt
  dvariable tmp;
  int gamma_regression_bh_flag=0;
  // do gamma regression intead of log-normal
  if (fsh.age_flags(142))
  {
    gamma_regression_bh_flag=1;
  }
  
  int nsize;
  if (gamma_regression_bh_flag)
  {
    dvar_vector mean_r=exp(loga+lb);
    dvar_vector obsr=exp(lr);
    MY_DOUBLE_TYPE s2=0.5/pen_wt;
    MY_DOUBLE_TYPE es2=exp(s2);
#if !defined(NO_MY_DOUBLE_TYPE)
    MY_DOUBLE_TYPE v=1.0/((es2-1.0L)*es2);
#else
    MY_DOUBLE_TYPE v=1.0/((es2-1.0)*es2);
#endif
    cout << "enter v " << endl;
    cin >> v;
    // put into the gamma dist alpha,beta parameterization
    MY_DOUBLE_TYPE alpha=v;
    dvar_vector beta=v/mean_r;
    int mmin;
    int mmax;
    if (fsh.age_flags(199)<=0)  // Use tmpdif over full model period ND 2Mar2011
    {
      mmin=obsr.indexmin();
      mmax=obsr.indexmax();
    }
    else
    {
      mmin=fsh.last_real_year-fsh.age_flags(199)+1;    //ND 8Mar2011
      mmax=fsh.last_real_year-fsh.age_flags(200);    //ND 8Mar2011
    }
    nsize=mmax-mmin+1;
    tmp=norm2(tmpdif(mmin,mmax));

    dvariable ll=0.0;
    for (i=mmin;i<=mmax;i++)
    {
      ll+=log_gamma_density(obsr(i),alpha,beta(i));
    } 
    f=-ll;  // negative log-likelihood
  }
  else
  {
    if (fsh.age_flags(94)==3  || fsh.age_flags(182) )
    {
      if (!fsh.af170q0 && !fsh.do_fishery_projections_flag //NMD_11Aug2022
          && !fsh.projection_sim_index)
      {
        fsh.bh_recr_devs=tmpdif;
      }
    }
    else
    {
      if (!fsh.af170q0 && !fsh.do_fishery_projections_flag //NMD_11Aug2022
          && !fsh.projection_sim_index)
      {
        fsh.bh_recr_devs(1,fsh.last_real_year-lag)=tmpdif;
      }
    }
     // This is to restrict computation of SRR to a range of years
    int ibh,jbh;
    MY_DOUBLE_TYPE s2=0.5/pen_wt;
    dvar_vector w;
    dvariable rhoest=0;
    if (fsh.age_flags(199)<=0)  // Use tmpdif over full model period ND 2Mar2011
    {
      ibh=tmpdif.indexmin();
      jbh=tmpdif.indexmax();
      nsize=tmpdif.indexmax()-tmpdif.indexmin()+1;
      w=tmpdif;
    }
    else
    {
      int mmax=tmpdif.indexmax();
      ibh=fsh.last_real_year-fsh.age_flags(199)+1;    //ND 8Mar2011
      jbh=fsh.last_real_year-fsh.age_flags(200);    //ND 8Mar2011
      if (fsh.age_flags(94)==3  || fsh.age_flags(182) )
      {
        ibh=(ibh-1)/fsh.age_flags(57)+1;
        jbh=(jbh-1)/fsh.age_flags(57)+1;
      }
      if (jbh > mmax || ibh < 1)
      {
        cerr << "Interval for BH-SRR deviates incorrect"
             << " check age_flags(199,200) values " << endl;
        ad_exit(1);
      }
      nsize=jbh-ibh+1;
      w=tmpdif(ibh,jbh);
    }
    if (pof)  // report the moment estimate of autocorrelation coefficient
    {
      dvariable num;  //NMD
      dvariable den;  //NMD
      dvar_vector mw=w-mean(w);
      num=mw(ibh,jbh-1)*mw(ibh+1,jbh).shift(ibh)/(jbh-ibh);
      den=norm2(mw)/(jbh-ibh+1);
      dvariable t=num/den;
      dvariable fpen=0.0;
      dvariable rhoest =  0.99-posfun(0.99-t,.01,fpen);
      (*pof) << "# BH autocorrelation moment estimator  " << rhoest << endl;  //NMD_26jul2018
      dvariable rho=fsh.species_pars(2,cs);
      (*pof) << "# BH autocorrelation parameter estimate  " << rho << endl;  //NMD_26jul2018
    }
    if (!fsh.age_flags(135))
    {
      f=pen_wt*norm2(w);
    }
    else
    {
      if (nsize !=w.indexmax()-w.indexmin()+1)
      {
        cerr << "bad nsize" << endl;
        ad_exit(1);
      }
      dvariable rho=fsh.species_pars(2,cs);
      dvar_vector v=chinv_multiply(rho,w);
      MY_DOUBLE_TYPE var=0.5/pen_wt;

      f=0.5*nsize*log(var)-0.5*log(1.0-rho*rho)
        +0.5*((1-rho*rho)/var)*(v*v);

      cout << "# BH autocorrelation: " << rho << endl;
      cout << "# AR1 penalty: " << f << endl;
    }
    tmp=norm2(w);
  }
  fsh.bh_variance=tmp/nsize;

   dvariable fpen=0.0;
   if (!fsh.zeroed_fisheries_flag)
   {
     if (cs==1)
     {
       fsh.alpha=alpha;
       fsh.beta=beta;
       fsh.steepness=steepness;
       fsh.phi=phi;
     }
     else
     {
       fsh.pmsd->alpha(cs)=alpha;
       fsh.pmsd->beta(cs)=beta;
       fsh.pmsd->steepness(cs)=steepness;
       fsh.pmsd->phi(cs)=phi;
     }
   }
   else
   {
     if (cs==1)
     {
       alpha=fsh.alpha;
       beta=fsh.beta;
       phi=fsh.phi;
       steepness=fsh.steepness;
     }
     else
     {
       alpha=fsh.pmsd->alpha(cs);
       beta=fsh.pmsd->beta(cs);
       phi=fsh.pmsd->phi(cs);
       steepness=fsh.pmsd->steepness(cs);
     }
   }
   int ly=fsh.last_real_year-lag;
   if (fsh.age_flags(94)==3  || fsh.age_flags(182) )
   {
     if (cs==1) //NMD_12Apr2021
     {
       fsh.bh_predicted_recruits=exp(value(loga)+value(lb));
     }
     else
     {
       fsh.pmsd->bh_predicted_recruits(cs)=exp(value(loga)+value(lb));
     }	 
   }
   else
   {
     if (cs==1) //NMD_12Apr2021
     {
       fsh.bh_predicted_recruits(1,ly)=exp(value(loga)+value(lb));
     }
     else
     {
       fsh.pmsd->bh_predicted_recruits(cs)(1,ly)=exp(value(loga)+value(lb));
     }
   }
   /*
   dvector cbio=value(bio(1,ly));
   dvector op=elem_div(value(alpha)*cbio,value(beta)+cbio);

   cout << "newbh.cpp: "<<norm2(fsh.bh_predicted_recruits(1,fsh.last_real_year-lag)-op)
        << endl;
    ofstream ofs("bh");
    ofs << "newbh.cpp: "
        <<norm2(fsh.bh_predicted_recruits(1,fsh.last_real_year-lag)-op) << endl;
    ofs << "newbh.cpp: "
        <<fsh.bh_predicted_recruits(1,fsh.last_real_year-lag)-op << endl;
    */
   { 
    //cout << "&alpha = " << &alpha << endl;
    //cout << "alpha = " << alpha << endl;
    //cout << "beta  = " << beta  << endl;
    cout << "phi   = " << phi   << endl;
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
      (*pof) << "# Beverton-Holt stock-recruitment relationship report";
      if (fsh.pmsd)
      (*pof) << " for species " << fsh.pmsd->current_species;
      (*pof) << endl;
      (*pof) << "# alpha = " << alpha << "  beta =  " << beta 
             << "  steepness = " << steepness
             << "  variance = " << fsh.bh_variance << endl;
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
         if (fsh.age_flags(161))
         {
           tmp2(ii)*=mfexp(0.5*fsh.bh_variance);
         }
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
      //yield_analysis_bh(fsh,pof,alpha,beta,lwc,pmature,Fay); //xyxxx
    }
  }
  //cout<< "just before call to yield_analysis_bh"<<endl; //xyxxx
  // DF  03 DEC 03
  if ( (pof) || fsh.age_flags(165)>0)
  {
    if (fsh.age_flags(140)==0) // af140 is for region-specific yield
    {
      if (fsh.age_flags(165)>0)  // Don't do this if only here for report
      {                        // generation.
        f += yield_analysis_bh_penalty(fsh,alpha,beta,pmature);
      }
      if (pof) 
      {
        if (fsh.age_flags(165)==0) // do once if not done above
        {
          int sflag=0;    //NMD_18Dec2013
          if (fsh.pmsd)
          {
            sflag=sum(column(fsh.pmsd->species_flags,2));
          }
          if (sflag) fsh.pmsd->current_species=1;   //NMD_18Dec2013
          yield_analysis_bh_penalty(fsh,alpha,beta,pmature);
          if (fsh.pmsd)
          cout  << " species " << fsh.pmsd->current_species;
          cout  << "    fsh.p_MSY "  <<  fsh.p_MSY << endl;
        }
        // now print to report file
        if (fsh.pmsd) //NMD_18Dec2013
        {
          if (sum(column(fsh.pmsd->species_flags,2)))  //multi-sex
          {
            ivector sf2=column(fsh.pmsd->species_flags,2);
            for (int i=1; i<=fsh.pmsd->num_species; i++)
            {
              if(sf2(i)) fsh.pmsd->current_species=i; //get female NMD_18Dec2013
            }
          }
        }
        yield_analysis_bh(fsh,pof,alpha,beta,pmature);
      }
    }
    else
    {
      if (pof)  { // --- this does region-specific yields in report (plot.rep)
        //cout<<"ABOUT TO ENTER folly IN newbh"<<endl;
        //yield_analysis_bh(fsh,0,alpha,beta,pmature,Fay);
        yield_analysis_bh_daves_folly(fsh,pof,alpha,beta,lwc,
          pmature);
      }
    }
  }
  return f+1000.0*fpen;
} 
