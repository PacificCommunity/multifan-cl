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

#if !defined(ZFEV)
extern dvar_len_fish_stock_history * pcfsh;
#endif
  extern dvar_vector * psv;


/*
MY_DOUBLE_TYPE rmax_penalty_weight=1.e+5; 
MY_DOUBLE_TYPE rmax_multiplier= 0.96;

void check_sanity(ivector &t)
{
  int mmin=t.indexmin();
  int mmax=t.indexmax();
  for (int i=mmin;i<mmax;i++)
  {
    if (t(i) != t(i+1))
    {
      cerr << "sanity error in move_periods" << endl;
      ad_exit(1);
    }
  }
}
*/
void setup_diffusion(dvar3_array& Dad,int nage,int num_regions,
  imatrix& Dflags, dvar_vector& diff_coffs,dvar_vector& diff_coffs2,
  dvar_vector& diff_coffs3,ivector& age_flags,pmulti_species_data & pmsd);
/*
void check_sanity(ivector &t, dvar3_array& tg, ivector& rip,int it,
  dvar_fish_stock_history& fsh)
{
  ivector& itp=fsh.initial_tag_period(it);
  check_sanity(t);
# if !defined(OPT_LIB)
    
    int mmin=tg.indexmin();
    int mmax=tg.indexmax();
    int err=0;
    int i;
    for (i=mmin;i<mmax;i++)
    {
      if (tg(i,rip(i)).indexmin() != tg(i+1,rip(i+1)).indexmin())
      {
        err=1;
        break;
      }
    }
    if (err)
    {
      cerr << "sanity error in youngest age class for tag group " 
           << it << endl
           << "for regions " << i << " and " << i+1
           << " mimum ages are "
           <<  tg(i,rip(i)).indexmin() << " and "
           <<  tg(i+1,rip(i+1)).indexmin() 
           << " fishing periods are " << endl
           <<  rip << endl
           << " initial tag periods for this tag group are "
           << itp << endl;
   
      cerr << "years are" << endl;
      for (i=mmin;i<=mmax;i++)
      {
        cerr << fsh.year(i,rip(i)) << " ";
      }
      cerr << endl;
      cerr << "months are" << endl;
      for (i=mmin;i<=mmax;i++)
      {
        cerr << fsh.month(i,rip(i)) << " ";
      }
      cerr << endl;

      ad_exit(1);
    }

#  endif
}
  int sip;
 */ 
void dvar_len_fish_stock_history::catch_equations_calc1(dvar_vector& sv,
  ivector* pq_flag,const prevariable& fzpen)
{
  //greport("beginning catch_equations_calc");
  tmprecr.initialize();
  tmpinitpop.initialize();
  MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
  if (af170q0==0)
  {
    tot_mort.initialize();
  }
  else
  {
    tot_mort_q0.initialize();
  }
  // calculate the totalfishing mortality and survival rates for each
  // fishing period in each region
  //cout <<"rshort3.cpp " <<  num_fish_periods << endl;
  int ir;
  dvar_vector * ps=0;
  dvar_vector * ptm=0;
  dvar_matrix * pfm=0;
  for (ir=1;ir<=num_regions;ir++)
  {
    //cout <<"rshort3.cpp " << fraction(ir) << endl;
    int tmp_nfp;
    if (do_fishery_projections_flag==0)
      tmp_nfp=num_real_fish_periods(ir);
    else
      tmp_nfp=num_fish_periods(ir);

    for (int ip=1;ip<=tmp_nfp;ip++)  // Loop over fishing periods
    {
      if (af170q0==0)
      {
        ptm=&(tot_mort(ir,ip));
        pfm=&(fish_mort(ir,ip));
        ps=&(survival(ir,ip));
      }
      else
      {
        ptm=&(tot_mort_q0(ir,ip));
        pfm=&(fish_mort_q0(ir,ip));
        ps=&(survival_q0(ir,ip));
      }
      

      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        (*ptm)+=mfexp((*pfm)(fi));
      }
      if (age_flags(115))   // put bounds on Z (acutally total F)
      {                     // for SS2 parameterization
        bounds_for_ss2(*ptm,*pfm,fzpen);
      }
      
      int af115=get_age_flags_region(ir,115);
      if (!pmsd)
      {
        if (af115==0)   
        {                     
          (*ptm)+=mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
        }
        else
        {  // for SS2 parameterization
           // this is Z=-log(S)=-log(1-sum F)+M
           (*ptm)=-log(1-*ptm)+mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
        }
      }
      else
      {
        if (af115==0)   
        {                     
          (*ptm)+=
            mfexp(get_nat_mort_region(ir)(year(ir,ip))+fraction(ir,ip));
        }
        else
        {  // for SS2 parameterization
           // this is Z=-log(S)=-log(1-sum F)+M
           (*ptm)=-log(1-*ptm)
             +mfexp(get_nat_mort_region(ir)(year(ir,ip))+fraction(ir,ip));
        }
      }
    
      (*ps)=mfexp(-(*ptm));
    }
  }
  do_fish_mort_intermediate_calcs();
}

void dvar_len_fish_stock_history::catch_equations_calc2(dvar_vector& sv,
  ivector* pq_flag,const prevariable& fzpen)
{
  if (num_regions>1 && !age_flags(114))
  {
    setup_diffusion();
  }
     if (ff263flag==0)
     {
#if !defined(ZFEV)
       if (!parest_flags(143))
#endif
       {
         if (!parest_flags(175))
           mean_length_calc(0);
         else
           mean_length_calcx();
       }
     }
#if !defined(ZFEV)
     else
     {
       // do this once at the beginning of the run
       if (parest_flags(144))
       {
         pcfsh->mean_length_calc();
         pcfsh->variance_calc();
       }
     }
#endif

  if (age_flags(177))
  {
    totpop=totpop_coff;
    if (pmsd)
    {
      pmsd->totpop=pmsd->totpop_coff;
    }

    if (pq_flag)
    {
      age_flags(192)=1;
    }
    get_initial_population(sv,0,pq_flag);
    age_flags(192)=0;
  }
  else
  {
    if (af170q0==0)
    {
      get_population_multipliers(sv,pq_flag);
      get_initial_population(sv,0,pq_flag);
    }
    else
    {
      if (age_flags(171)==0)
      {
        get_population_multipliers(sv,pq_flag);
      }
      // get the initial age structure
      if (!age_flags(94))
      {
        get_initial_age_structure(totpop,sv);
      }
      else
      {  
        xget_initial_age_structure_equilibrium();
      }
    }
  }
  if (projection_sim_index>0 && do_fishery_projections_flag)
  {
    //ofstream ofs("projected_randomized_catches");
    int num_times=age_flags(20);
    random_number_generator rng(217+projection_sim_index);
    randrec.allocate(num_real_years,nyears);
    MY_DOUBLE_TYPE sd=age_flags(21)/100.;
    gradient_structure::set_NO_DERIVATIVES(); 
    {
      //sd=0.0;
      randrec.fill_randn(rng);
      randrec=exp(sd*randrec);
  
      get_numbers_at_age(sv,pq_flag);
     
      if(!pq_flag) *proj_output_files[0] << "#simulation " << projection_sim_index << endl;  //NMD 11Nov 2011
      if(!pq_flag) *proj_output_files[0] << "#catch " << endl;   //NMD 11Nov 2011
      for (int ir=1;ir<=num_regions;ir++)
      {
        if(!pq_flag) *proj_output_files[0]  << "#region " << ir << endl;      //NMD 11Nov 2011
        //for (int ip=num_real_fish_periods(ir)+1;ip<=num_fish_periods(ir);ip++)
        for (int ip=1;ip<=num_fish_periods(ir);ip++)
        {
          if (ip==num_real_fish_periods(ir)+1)
            if(!pq_flag) *(proj_output_files[0]) << "#projected catch" << endl;   //NMD 11Nov 2011
          for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)
          {
            //(*pc)(ir,ip,fi)=(*pfmc)(ir,ip,fi)+(*pnum_fish)(ir,ip);
            if(!pq_flag) *(proj_output_files[0]) << exp(catch(ir,ip,fi)) << endl;  //NMD 11Nov 2011
          }
        }
      }
      if(!pq_flag) *proj_output_files[0] << "#num fish in pop " << endl;    //NMD 11Nov 2011
      for (int ir=1;ir<=num_regions;ir++)
      {
        if(!pq_flag) *(proj_output_files[0]) << "#region " << ir << endl;    //NMD 11Nov 2011
        //for (int ip=num_real_fish_periods(ir)+1;ip<=num_fish_periods(ir);ip++)
        for (int ip=1;ip<=num_fish_periods(ir);ip++)
        {
          if (ip==num_real_fish_periods(ir)+1 && !pq_flag)        //NMD 11Nov 2011
             *(proj_output_files[0]) << "#projected numbers" << endl;
          if(!pq_flag) *(proj_output_files[0])  << exp(num_fish(ir,ip)) << endl;  //NMD 11Nov 2011
        }
      }
    }
  }

  get_numbers_at_age(sv,pq_flag);
      

  if (age_flags(180))
  {   
    (*clogf)  <<  "num_fish" << endl;
    (*clogf)  <<  num_fish << endl;
  }
  if (age_flags(180))
  {   
    (*clogf)  <<  "catch" << endl;
    (*clogf)  <<  catch << endl;
  }
}

