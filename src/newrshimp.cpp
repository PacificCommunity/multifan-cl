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

ofstream * pofs_num_fish=0;

void check_sanity(ivector &t);

void check_sanity(ivector &t, dvar3_array& tg, ivector& rip,int it,
  dvar_fish_stock_history& fsh);
extern int sip;

static void set_xxx(int & ipold,int & ip)
{ 
  ipold=ip;
}

static void newxxx(void)
{
}

void dvar_len_fish_stock_history::check_implict_catch_options
  (int ir,int ip,dvariable& ffpen,dvar_vector& fmlev,ivector& mci,
  dvar_matrix& fmlevdevs,dvariable& ffpen1)
{
  // just do this for now 
  /*
 if (num_fish_incidents(ir,ip)>0)
 {
   do_newton_raphson(ir,ip,ffpen,fm_level(ir,ip),
    missing_catch_counter,fm_level_devs,fpen1);
 }

  if (implicit_flag) // do some form of implcit stuff if required
  {
    // decide whether implicit catch stuff should be done
    if (num_missing_effort_by_region(ir,ip)>0)
    {
      // there are catches with no effort so fit these to observed total
      // catch via newton raphson
      if (num_missing_effort_by_region(ir,ip)==num_fish_incidents(ir,ip))
      {
        // all effort is missing for this fishing period so fit iall catchesd
        // via  newton-raphson 
        //
        num_present_effort_by_region(ir,ip)++;
      }
      else
      {
        // effort is missing for only some fisheries so check for
        // missing catches
      }
    }
  }
  */
}
  

  
void dvar_len_fish_stock_history::new_catch_equations_calc_implicit
  (dvar_vector& sv, ivector* pq_flag, dvariable & ffpen)
{
  if (pofs_num_fish==0)
  {
   // pofs_num_fish=new ofstream("num_fish");
  }
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
  int ir;
  if (num_regions>1 && !age_flags(114))
  {
    setup_diffusion();
  }
#if !defined(ZFEV)
     if (!parest_flags(143))
#endif
     {
       if (!parest_flags(175))
         mean_length_calc(0);
       else
         mean_length_calcx();
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
      get_initial_population(sv,0,0);
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
  // do we have the variance for this?
  calculate_the_mean_weight();

  dvar_matrix xvtmp;
  xvtmp=new_get_numbers_at_age_implicit(sv,pq_flag,ffpen);

  get_implicit_catchability(*this);

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

imatrix dvar_fish_stock_history::get_survey_samples_flags(void)
{
  imatrix m(1,num_regions,1,num_fish_periods);
  m.initialize();
  ivector ff92=column(fish_flags,92);
  for (int ir=1;ir<=num_regions;ir++)
  {
    ivector nfp;    //NMD_28jun2022
//    if (sum(data_fish_flags(2)))
    if (sum(data_fish_flags(2)) && !projection_sim_index) //NMD_21Aug2023
    {
      nfp=num_real_fish_periods;
    }
    else
    {
      nfp=num_fish_periods;
    }

//    for (int ip=1;ip<=num_fish_periods(ir);ip++) 
    for (int ip=1;ip<=nfp(ir);ip++) 
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)
      {
        int pi=parent(ir,ip,fi);
        if (ff92(pi))
        {
          m(ir,ip)=1;
        }
      }
    }
  }
  return m;   
}

dvar_matrix dvar_len_fish_stock_history::new_get_numbers_at_age_implicit
  (dvar_vector& sv,void * pq_flag ,dvariable& ffpen)
{
  ivector ff92=column(fish_flags,92);
  imatrix have_survey_samples;
  if (ff92sum)
  {
    have_survey_samples=get_survey_samples_flags();
  }
  ivector nft;    //NMD_28jun2022
//  if (sum(data_fish_flags(2)))
  if (sum(data_fish_flags(2)) && !projection_sim_index)  //NMD_21aug2023
  {
    nft=num_real_fish_times;
  }
  else
  {
    nft=num_fish_times;
  }
  if (ff92sum)
  {
//    if (!allocated(vul_at_age))
    if (!allocated(vul_at_age) || projection_sim_index)
    {
      if (projection_sim_index) vul_at_age.deallocate();
      vul_at_age.allocate(1,num_fisheries);
      for (int fi=1;fi<=num_fisheries;fi++)
      {
        if (ff92(fi))
        { 
          vul_at_age(fi).allocate(1,nft(fi),1,nage);
        }
      }
    }
    vul_at_age.initialize();
  }

  Zmax_flag=0;
  Zmax_fish=1.0;
  if (age_flags(116)!=0)
  {
    Zmax_fish=age_flags(116)/100.;
  }

  dvariable fpen1=0.0;
  missing_catch_counter.initialize();
  dvar_vector average_recruitment;
  dvar_matrix average_recruitment_by_season;
  dvar_matrix average_recruitment_by_season_for_regions;
  dvar_matrix recruitment_scalar_by_season_for_regions;
  dvar_vector average_recruitment_for_projections;

  if (age_flags(190)>0)
  {
    if (!af170q0)
    {
      average_recruitment_for_projections=
        get_average_recruitment_for_projections
       (N,num_regions,last_real_year,age_flags(190),age_flags(191));
    }
    else
    {
      average_recruitment_for_projections=
        get_average_recruitment_for_projections
       (N_q0,num_regions,last_real_year,age_flags(190),age_flags(191));
    }

    if (parest_flags(228)>0)
      average_recruitment_for_projections=1.e-8;
  }

  int year_change_flag=0;
  dvar3_array * pnum_fish=0;
  dvar3_array * pN=0;
  dvar3_array * ptm=0;
  dvar3_array * psurv=0;
  //if (!pq_flag)
  if (af170q0==0)
  {
    pnum_fish=&num_fish;
    ptm=&tot_mort;
    psurv=&survival;
    pN=&N;
  }
  else
  {
    ptm=&tot_mort_q0;
    pnum_fish=&num_fish_q0;
    psurv=&survival_q0;
    pN=&N_q0;
  }
  int ir;
  int current_year=1;
  ivector rip(1,num_regions);
  rip=1;
  int finished_flag=1;
  int break_flag=0;
  ivector tmp_mp(1,num_regions);
  ivector tmp_yr(1,num_regions);
  ivector tmp_mn(1,num_regions);
  ivector tmp_ip(1,num_regions);
  tmp_mp.initialize();
  tmp_yr.initialize();
  tmp_mn.initialize();
  tmp_ip.initialize();
  dvar_matrix vtmp;    //NMD_8jun2022
  if (!pmsd)
  {
    vtmp.allocate(1,1,1,last_real_year+1);   //NMD_8jun2022
  }
  else
  {
    vtmp.allocate(1,pmsd->num_species,1,last_real_year+1);  //NMD_8jun2022
  }
  vtmp.initialize();

  for (ir=1;ir<=num_regions;ir++)
  {
    (*pnum_fish)(ir,1)=(*pN)(ir,1);
  }
  dvariable alpha=0.0;
  dvariable beta=0.0;
  int have_ab_flag=0;

  int no_movement_flag=0;
  int year_save=0;
  int year_fill=0;  //NMD_8jun2022
  random_number_generator * rng2=0;
  if (age_flags(196)>0 && projection_sim_index>0)
  {
    rng2=new random_number_generator(age_flags(197)+projection_sim_index);
  }

  ivector fin_cal_year(1,num_regions);
  fin_cal_year.initialize();
  if (parest_flags(377)==0 && parest_flags(378)==0)
  {
    if (age_flags(57)==1)
    {
      fin_cal_year=last_real_year;
    }
    else
    {
      fin_cal_year=really_true_year(1,last_real_year);
    }
  }
  
  int irold=0;
  int ipold=0;
  do
  {
    finished_flag=1;
    year_save=year(1,rip(1));
    no_movement_flag=0;
    for (int ir=1;ir<=num_regions;ir++)
    {
      int ng=nage_by_region(ir);
      break_flag=0;
      int& ip=rip(ir);
      do
      {
        if (ip==ipold && ir==irold)
        {
          ip++;
          if (ip>num_fish_periods(ir)) break;
        }
        irold=ir;
        set_xxx(ipold,ip);
        if ( ip==num_real_fish_periods(ir))
        {
          //cout << "HERE" << endl;
        }

        dvariable srr;
        if ( ip==num_real_fish_periods(ir) && do_fishery_projections_flag==1 )
        {
          if (have_ab_flag==0)
          {
            if(!age_flags(190) || age_flags(195)>0)   //NMD_12jul2022
            {
              get_bh_alpha_and_beta(*this,alpha,beta);
              have_ab_flag=1;
            }
            if(af170q0 == 0)     //NMD_2Jul2013
            {
              average_recruitment=get_average_recruitment(N,num_regions,
                last_real_year);
            } else {
              average_recruitment=get_average_recruitment(N_q0,num_regions,
                last_real_year);
            }                    //NMD_2Jul2013
            average_recruitment/=sum(average_recruitment);
            average_recruitment=log(1.e-20+average_recruitment);
            {
              if(af170q0 == 0)     //NMD_2Jul2013
              {
                average_recruitment_by_season=
                  get_average_recruitment_by_season(N,num_regions,
                  last_real_year,age_flags(57));
              }
              else
              {
                average_recruitment_by_season=
                  get_average_recruitment_by_season(N_q0,num_regions,
                  last_real_year,age_flags(57));
              }
              MY_DOUBLE_TYPE ns=age_flags(57);
              if (ns==0) ns=1; 
              // normalize to sum to ns over all reagions
              // and seasons.
              average_recruitment_by_season/=
                sum(average_recruitment_by_season)/ns;
               
              average_recruitment_by_season=
                log(1.e-20+average_recruitment_by_season);

              average_recruitment_by_season_for_regions=
                get_average_recruitment_by_season_for_region
                (average_recruitment_by_season,num_regions,age_flags(57));

              recruitment_scalar_by_season_for_regions=
                get_scalar_recruitment_by_season_for_region
                (average_recruitment_by_season_for_regions,
                 num_regions,age_flags(57));
            }
          }
          // get recruitment for first projection period
          int lag=age_flags(147);
          dvariable r=get_bh_recruitment_for_projections
            (year(ir,ip+1)-lag,*this,alpha,beta);
          if (age_flags(182)>0 && age_flags(57)>1)
          {
            r/=age_flags(57);   //NMD25Aug2015
          }
          srr = r;
          if ( age_flags(190)==0) // Use SRR for projected recruitment
          {
            // says average_recruitment but means log_normalized_recruitment
            if (age_flags(183)==0)
            {
              (*pN)(ir,year(ir,ip+1),1)=log(r)+average_recruitment(ir);
            }
            else
            {
              int is=(year(ir,ip+1)-1)%age_flags(57)+1;
              (*pN)(ir,year(ir,ip+1),1)=log(r)+
                average_recruitment_by_season(is,ir);
            }
          }
          // Use average recruitment for periods defined by af(190), af(191)
          if (age_flags(190)>0 && age_flags(195)==0)
          {
            if (age_flags(183)==0)
            {
              (*pN)(ir,year(ir,ip+1),1)=log(1.e-20+average_recruitment_for_projections(ir));
            }
            else
            {
              int is=(year(ir,ip+1)-1)%age_flags(57)+1;
              (*pN)(ir,year(ir,ip+1),1)=log(1.e-20+average_recruitment_for_projections(ir))
		+recruitment_scalar_by_season_for_regions(is,ir);
            }
//            (*pN)(ir,year(ir,ip+1),1)=log(1.e-20+average_recruitment_for_projections(ir));
          }
          // Use SRR total rec but regional distribution according to
          // average distribution for periods defined by af(190), af(191)
          if (age_flags(190)>0 && age_flags(195)>0)
          {
            (*pN)(ir,year(ir,ip+1),1)=log(r)+
               log((1.e-20+average_recruitment_for_projections(ir))/
               (1.e-20+sum(average_recruitment_for_projections)));
          }
        }
        if (ip>num_real_fish_periods(ir) && do_fishery_projections_flag==1)
        {
          //check_implict_catch_options(ir,ip,ffpen,fm_level(ir,ip),
          //       missing_catch_counter,fm_level_devs,fpen1);

          if (num_fish_incidents(ir,ip)>0)
          {
            {
              do_newton_raphson(ir,ip,ffpen,fm_level(ir,ip),
                 missing_catch_counter,fm_level_devs,fpen1);
            }
          }
          else
          {
            if (!pmsd)
            {
               (*ptm)(ir,ip)=mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
            }
            else
            {
               (*ptm)(ir,ip)=mfexp(get_nat_mort_region(ir)(year(ir,ip))
                  +fraction(ir,ip));
            }
            (*psurv)(ir,ip)=mfexp(-(*ptm)(ir,ip));
          }
          if (allocated(have_survey_samples) && have_survey_samples(ir,ip))
          {
            dvar_vector nmid;
            nmid = (*pnum_fish)(ir,ip) -0.5*(*ptm)(ir,ip);
            for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)
            {
              int pi=parent(ir,ip,fi);
              if (ff92(pi))
              {
                // calculate the vulnerable numbers at age
                int ft=fish_times(ir,ip,fi);
                vul_at_age(pi,ft)=exp(incident_sel(ir,ip,fi)+nmid);
              }
            }   
          }
        }
        if ( ip>=num_fish_periods(ir) ||
          ( ip>num_real_fish_periods(ir)  && do_fishery_projections_flag==0 ) )
        {
          if (num_fish_incidents(ir,ip)>0)
          {
            do_newton_raphson(ir,ip,ffpen,fm_level(ir,ip),
               missing_catch_counter,fm_level_devs,fpen1);
          }
          else
          {
            if (!pmsd)
            {
              (*ptm)(ir,ip)=mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
            }
            else
            {
              (*ptm)(ir,ip)=mfexp(get_nat_mort_region(ir)(year(ir,ip))
                +fraction(ir,ip));
            }
            // *************************************************
            // *************************************************

            dvar_vector temp_tm(1,nage);
            temp_tm=(*ptm)(ir,ip);
            int jmin=(*ptm)(ir,ip).indexmin();
            int jmax=(*ptm)(ir,ip).indexmax();
            for (int j=jmin;j<=jmax;j++)
            {
              if ( (*ptm)(ir,ip)(j)>Zmax_fish)
              {
                dvariable u=(*ptm)(ir,ip)(j)-Zmax_fish;
                dvariable v=Zmax_fish + u/(1.0 + 3.0*u/Zmax_fish);
                temp_tm(j)=v;
                Zmax_flag=1;
              }
              MY_DOUBLE_TYPE pen=100.0;
              double Zcut_mult=0.8;
              if (age_flags(189)>0)
              {
                Zcut_mult=age_flags(189)/100.0;
              }
              double Zcut=Zcut_mult*Zmax_fish;
              if ( (*ptm)(ir,ip)(j)>Zcut)
              {
                MY_DOUBLE_TYPE pen=100.0;
                if (parest_flags(382))
                {
                  pen=parest_flags(382);
                }
                //ffpen1+=pen*square((*ptm)(ir,ip)(j)-Zcut);
                ffpen+=pen*log(1.0+square((*ptm)(ir,ip)(j)-Zcut));
              }
            }
            (*ptm)(ir,ip)=temp_tm;
            (*psurv)(ir,ip)=mfexp(-(*ptm)(ir,ip));
          }
          if (allocated(have_survey_samples) && have_survey_samples(ir,ip))
          {
            dvar_vector nmid;
            nmid = (*pnum_fish)(ir,ip) -0.5*(*ptm)(ir,ip);
            for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)
            {
              int pi=parent(ir,ip,fi);
              if (ff92(pi))
              {
                // calculate the vulnerable numbers at age
                int ft=fish_times(ir,ip,fi);
                vul_at_age(pi,ft)=exp(incident_sel(ir,ip,fi)+nmid);
              }
            }   
          }
          no_movement_flag=1; // we are at thenb end so ip does not
                              // get incremented so don't move the fish
                              // df july 13 05
          break;
        }
    
        finished_flag=0;
        if (year(ir,ip+1)==year(ir,ip))
        {
          if (num_fish_incidents(ir,ip)>0)
          {
            do_newton_raphson(ir,ip,ffpen,fm_level(ir,ip),
               missing_catch_counter,fm_level_devs,fpen1);
          }
          else
          {
            if (!pmsd)
            {
              (*ptm)(ir,ip)=mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
            }
            else
            {
              (*ptm)(ir,ip)=mfexp(get_nat_mort_region(ir)(year(ir,ip))
                 +fraction(ir,ip));
            }
            (*psurv)(ir,ip)=mfexp(-(*ptm)(ir,ip));
          }
          (*pnum_fish)(ir,ip+1)=(*pnum_fish)(ir,ip)-(*ptm)(ir,ip);

          if ((parest_flags(377)==0 && parest_flags(378)==0) &&
              do_fishery_projections_flag==1 && generate_report==1) //NMD_7may2024
          {
            if (age_flags(57)==1 && year(ir,ip)==fin_cal_year(ir))
            {
              put_terminal_catchability(ir,ip);
            }
            else if (age_flags(57)!=1 && really_true_year(ir,ip)==fin_cal_year(ir))
            {
              put_terminal_catchability(ir,ip);
            }
          }
	  
          if (pofs_num_fish)
          {
            if (ip==1)
            {
              *pofs_num_fish << ip << endl;
              *pofs_num_fish << exp((*pnum_fish)(ir,ip))<< endl;
            }
            *pofs_num_fish << ip+1 << endl;
            *pofs_num_fish << exp((*pnum_fish)(ir,ip+1)-(*ptm)(ir,ip))<< endl;
            *pofs_num_fish << " exp((*ptm)(ir,ip))  " 
              << exp((*ptm)(ir,ip))<< endl;
            *pofs_num_fish << " (*ptm)(ir,ip)  " << (*ptm)(ir,ip)<< endl;
          }

          // fit to survey samples based on ff92p
          if (allocated(have_survey_samples) && have_survey_samples(ir,ip))
          {
            dvar_vector nmid;
            nmid = (*pnum_fish)(ir,ip) -0.5*(*ptm)(ir,ip);
            for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)
            {
              int pi=parent(ir,ip,fi);
              if (ff92(pi))
              {
                // calculate the vulnerable numbers at age
                int ft=fish_times(ir,ip,fi);
                vul_at_age(pi,ft)=exp(incident_sel(ir,ip,fi)+nmid);
              }
            }   
          }

        }
        else
        {
          if (ip>num_real_fish_periods(ir))
          {
            // get recruitment for projection period
            int lag=age_flags(147);
            dvariable r=get_bh_recruitment_for_projections
              (year(ir,ip+1)-lag,*this,alpha,beta);
            if (age_flags(182)>0 && age_flags(57)>1) 
              r/=age_flags(57);   //NMD_12jul2022
            srr = r;
            if ( age_flags(190)==0) // Use SRR for projected recruitment
            {
              if (age_flags(183)==0)
              {
                (*pN)(ir,year(ir,ip+1),1)=log(r)+average_recruitment(ir);
              }
              else
              {
                int is=(year(ir,ip+1)-1)%age_flags(57)+1;
                (*pN)(ir,year(ir,ip+1),1)=log(r)+
                  average_recruitment_by_season(is,ir);
              }
              if (age_flags(196)>0)
              {
                MY_DOUBLE_TYPE sd=age_flags(196)/100.;
                MY_DOUBLE_TYPE bias=0.5*sd*sd;
                MY_DOUBLE_TYPE eps=randn(*rng2);
                (*pN)(ir,year(ir,ip+1),1)+=sd*eps-bias;

              }
            }
            // Use average recruitment for periods defined by af(190), af(191)
            if (age_flags(190)>0 && age_flags(195)==0)
            {
              if (age_flags(183)==0)
              {
                (*pN)(ir,year(ir,ip+1),1)=log(1.e-20+average_recruitment_for_projections(ir));
              }
              else
              {
                int is=(year(ir,ip+1)-1)%age_flags(57)+1;
                (*pN)(ir,year(ir,ip+1),1)=log(1.e-20+average_recruitment_for_projections(ir))
		  +recruitment_scalar_by_season_for_regions(is,ir);
              }
//              (*pN)(ir,year(ir,ip+1),1)=log(1.e-20+average_recruitment_for_projections(ir));
            }
            // Use SRR total rec but regional distribution according to
            // average distribution for periods defined by af(190), af(191)
            if (age_flags(190)>0 && age_flags(195)>0)
            {
              (*pN)(ir,year(ir,ip+1),1)=log(r)+
                 log((1.e-20+average_recruitment_for_projections(ir))/
                 (1.e-20+sum(average_recruitment_for_projections)));
            }

            if (age_flags(20)>0 && projection_sim_index>0)
            {
              // put in recruitment uncertainty for projections
              if(parest_flags(239) == 0)
              {  //NMD 7Mar2013
                if(age_flags(182) == 0)     //NMD 16Aug2022
                {
                  (*pN)(ir,year(ir,ip+1),1)=
                    log(simulated_recruitments(projection_sim_index,ir,year(ir,ip+1)));
                }
                else if (age_flags(182)>0 && age_flags(57)>1 && age_flags(183)>0)
                {
                  double tmprcr;
                  tmprcr=simulated_recruitments(projection_sim_index,ir,year(ir,ip+1));
                  int is=(year(ir,ip+1)-1)%age_flags(57)+1;
                  (*pN)(ir,year(ir,ip+1),1)=log(tmprcr)
                    +average_recruitment_by_season_for_regions(is,ir);
                }
                else
                {
                  cerr << "This can't happen in stochastic projections" << endl;
                  ad_exit(1);
                }     //NMD_16Aug2022
              }
              else if(parest_flags(239) == 1)
              {
                dvar_vector avrecr;    //NMD12Dec2012
                dvariable tot_avrecr;

                if(age_flags(199)>0)
                {
                  if (!af170q0)  //NMD_9Feb2015
      	          {
                    tot_avrecr=get_average_log_recruitment_for_projections
                      (N,num_regions,last_real_year,age_flags(199),age_flags(200),
                       age_flags(57),parest_flags(232),parest_flags(233),
                       age_flags(182));  //NMD_1Sep2015
                  }
                  else
       	          {
                    tot_avrecr=get_average_log_recruitment_for_projections
                      (N_q0,num_regions,last_real_year,age_flags(199),age_flags(200),
                       age_flags(57),parest_flags(232),parest_flags(233),
                       age_flags(182));  //NMD_1Sep2015
                  }

                  if (!allocated(simyears))
                  {
                    cerr << "Error -- trying to do simulations that requirer simyears"
                            " but it is not allocated " << endl;
                    ad_exit(1);
                  }
                  int psi=projection_sim_index;
                  int calyr=(year(ir,ip+1)-1)/age_flags(57)+1;
                  int is=(year(ir,ip+1)-1)%age_flags(57)+1;
                  dvariable deviate1;
                  if (age_flags(182))     //NMD_7May2018
                  {
                    deviate1=bh_recr_devs(simyears(psi,calyr));
                  }
                  else
                  {
                    deviate1=bh_recr_devs(simyears(psi,year(ir,ip+1)));		    
                  }         //NMD_7May2018

                  if(age_flags(161))
                  {
                    (*pN)(ir,year(ir,ip+1),1)= 
                         exp(average_recruitment_by_season(is,ir)
                         +log(srr) + deviate1 - bh_variance/2);
                  }
                  else
                  {
                    (*pN)(ir,year(ir,ip+1),1)= 
                         exp(average_recruitment_by_season(is,ir)
                         +log(srr) + deviate1);
                  }
                }
                else
                {
                  tot_avrecr=totpop;
                  int psi=projection_sim_index;
                  int calyr=(year(ir,ip+1)-1)/age_flags(57)+1;
                  int is=(year(ir,ip+1)-1)%age_flags(57)+1;
                  dvariable deviate1;
                  if (age_flags(182))     //NMD_7May2018
                  {
                    deviate1=bh_recr_devs(simyears(psi,calyr));
                  }
                  else
                  {
                    deviate1=bh_recr_devs(simyears(psi,year(ir,ip+1)));		    
                  }         //NMD_7May2018
                  if(age_flags(161))
                  {
                    (*pN)(ir,year(ir,ip+1),1)= 
                         exp(average_recruitment_by_season(is,ir)
                         +log(srr) + deviate1 - bh_variance/2);
                  }
                  else
                  {
                    (*pN)(ir,year(ir,ip+1),1)= 
                         exp(average_recruitment_by_season(is,ir)
                         +log(srr) + deviate1);
                  }				
                  if (age_flags(182)==1 && age_flags(57)>1)  //NMD_27Aug2015
                  {
                    cerr << "Error simulation recruitments don't support this option"
                         << endl;
                    ad_exit(1);
                  }
                }
              }
            }

          }
          if (ip<=num_real_fish_periods(ir))
          {
            if (num_fish_incidents(ir,ip)>0)
            {
              do_newton_raphson(ir,ip,ffpen,fm_level(ir,ip),
                missing_catch_counter,fm_level_devs,fpen1);
            }
            else
            {
              if (!pmsd)
              {
                (*ptm)(ir,ip)=mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
              }
              else
              {
                (*ptm)(ir,ip)=mfexp(get_nat_mort_region(ir)(year(ir,ip))+fraction(ir,ip));
              }
              (*psurv)(ir,ip)=mfexp(-(*ptm)(ir,ip));
            }

            if ((parest_flags(377)==0 && parest_flags(378)==0) &&
                 do_fishery_projections_flag==1 && generate_report==1) //NMD_7may2024
            {
              if (age_flags(57)==1 && year(ir,ip)==fin_cal_year(ir))
              {
                put_terminal_catchability(ir,ip);
              }
              else if (age_flags(57)!=1 && really_true_year(ir,ip)==fin_cal_year(ir))
              {
                put_terminal_catchability(ir,ip);
              }
            }


            if (allocated(have_survey_samples) && have_survey_samples(ir,ip))
            {
              dvar_vector nmid;
              nmid = (*pnum_fish)(ir,ip) -0.5*(*ptm)(ir,ip);
              for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)
              {
                int pi=parent(ir,ip,fi);
                if (ff92(pi))
                {
                  // calculate the vulnerable numbers at age
                  int ft=fish_times(ir,ip,fi);
                  vul_at_age(pi,ft)=exp(incident_sel(ir,ip,fi)+nmid);
                }
              }   
            }
          }

          // age the fish
          if (ip==num_real_fish_periods(ir) 
            && do_fishery_projections_flag==1
            && projection_sim_index>0)
          {
            if(parest_flags(239) == 0)
            {     //NMD 15Mar2012
              if(age_flags(182) == 0)     //NMD 16Aug2022
              {
                (*pnum_fish)(ir,ip+1,1)=
                  log(simulated_recruitments(projection_sim_index,ir,year(ir,ip+1)));
              }
              else if (age_flags(182)>0 && age_flags(57)>1 && age_flags(183)>0)
              {
                double tmprcr;
                tmprcr=simulated_recruitments(projection_sim_index,ir,year(ir,ip+1));
                int is=(year(ir,ip+1)-1)%age_flags(57)+1;
                (*pnum_fish)(ir,ip+1,1)=log(tmprcr)
                  +average_recruitment_by_season_for_regions(is,ir);
              }
              else
              {
                cerr << "This can't happen in stochastic projections" << endl;
                ad_exit(1);
              }     //NMD_16Aug2022
//              (*pnum_fish)(ir,ip+1,1)= log(
//                simulated_recruitments(projection_sim_index,ir,year(ir,ip+1)));
            }
            else if(parest_flags(239) == 1)     //NMD 15Mar2012
            {
              dvar_vector avrecr;    //NMD12Dec2012
              dvariable tot_avrecr;				
              if(age_flags(199)>0)
              {
                if (!af170q0)  //NMD_9Feb2015
                {
                  tot_avrecr=get_average_log_recruitment_for_projections
                    (N,num_regions,last_real_year,age_flags(199),age_flags(200),
                     age_flags(57),parest_flags(232),parest_flags(233),
                     age_flags(182));  //NMD_1Sep2015
                }
                else
                {
                  tot_avrecr=get_average_log_recruitment_for_projections
                    (N_q0,num_regions,last_real_year,age_flags(199),
                     age_flags(200),
                     age_flags(57),parest_flags(232),parest_flags(233),
                     age_flags(182));  //NMD_1Sep2015
                }
                int psi=projection_sim_index;
                int calyr=(year(ir,ip+1)-1)/age_flags(57)+1;
                int is=(year(ir,ip+1)-1)%age_flags(57)+1;
                dvariable deviate1;
                if (!allocated(simyears))
                {
                  cerr << "Error -- trying to do simulations that require"
                          " simyears but it is not allocated " << endl;
                  ad_exit(1);
                }
                if (age_flags(182))     //NMD_7May2018
                {
                  deviate1=bh_recr_devs(simyears(psi,calyr));
                }
                else
                {
                  deviate1=bh_recr_devs(simyears(psi,year(ir,ip+1)));
                }         //NMD_7May2018

                if(age_flags(161))
                {
                  (*pnum_fish)(ir,ip+1,1)= 
                       average_recruitment_by_season(is,ir)+
                       log(srr) + deviate1 - bh_variance/2;
                }
                else
                {
                  (*pnum_fish)(ir,ip+1,1)= 
                       average_recruitment_by_season(is,ir)+
                       log(srr) + deviate1;
                }
              }
              else
              {
                tot_avrecr=totpop;
                int psi=projection_sim_index;
                int calyr=(year(ir,ip+1)-1)/age_flags(57)+1;
                int is=(year(ir,ip+1)-1)%age_flags(57)+1;
                dvariable deviate1;
                if (age_flags(182))     //NMD_7May2018
                {
                  deviate1=bh_recr_devs(simyears(psi,calyr));
                }
                else
                {
                  deviate1=bh_recr_devs(simyears(psi,year(ir,ip+1)));
                }         //NMD_7May2018

                if(age_flags(161))
                {
                  (*pnum_fish)(ir,ip+1,1)= 
                       average_recruitment_by_season(is,ir)+
                       log(srr) + deviate1 - bh_variance/2;
                }
                else
                {
                  (*pnum_fish)(ir,ip+1,1)= 
                       average_recruitment_by_season(is,ir)+
                       log(srr) + deviate1;
                }
                if (age_flags(182)==1 && age_flags(57)>1)  //NMD_27Aug2015
                {
                  cerr << "Error simulation recruitments don't support this option"
                       << endl;
                  ad_exit(1);
                }
              }    //NMD12Dec2012
            }

            if (ng>2)   //NMD 22Feb2012
            {
              if(af170q0==1)
              {
                (*pnum_fish)(ir,ip+1)(2,ng)=
                  log(simulated_numbers_at_age_noeff(projection_sim_index,ir));
              }
              else
              {
                (*pnum_fish)(ir,ip+1)(2,ng)=
                  log(simulated_numbers_at_age(projection_sim_index,ir));
              }
            }    //NMD 22Feb2012
             
            (*pN)(ir,year(ir,ip+1))=(*pnum_fish)(ir,ip+1);


          }
          else if (ip>num_real_fish_periods(ir) 
            && do_fishery_projections_flag==1
            && projection_sim_index>0)
          {
            if(parest_flags(239) == 0)
            {     //NMD 15Mar2012
              if(age_flags(182) == 0)     //NMD 16Aug2022
              {
                (*pnum_fish)(ir,ip+1,1)=
                  log(simulated_recruitments(projection_sim_index,ir,year(ir,ip+1)));
              }
              else if (age_flags(182)>0 && age_flags(57)>1 && age_flags(183)>0)
              {
                double tmprcr;
                tmprcr=simulated_recruitments(projection_sim_index,ir,year(ir,ip+1));
                int is=(year(ir,ip+1)-1)%age_flags(57)+1;
                (*pnum_fish)(ir,ip+1,1)=log(tmprcr)
                  +average_recruitment_by_season_for_regions(is,ir);
              }
              else
              {
                cerr << "This can't happen in stochastic projections" << endl;
                ad_exit(1);
              }     //NMD_16Aug2022
//              (*pnum_fish)(ir,ip+1,1)=
//                log(simulated_recruitments(projection_sim_index,ir,year(ir,ip+1)));
            }
            else if(parest_flags(239) == 1)     //NMD 15Mar2012
            {
              dvar_vector avrecr;    //NMD12Dec2012
              dvariable tot_avrecr;
              if(age_flags(199)>0)
              {
                if (!af170q0)  //NMD_9Feb2015
                {
                  tot_avrecr=get_average_log_recruitment_for_projections
                    (N,num_regions,last_real_year,age_flags(199),age_flags(200),
                     age_flags(57),parest_flags(232),parest_flags(233),
                     age_flags(182));  //NMD_1Sep2015
                }
                else
                {
                  tot_avrecr=get_average_log_recruitment_for_projections
                    (N_q0,num_regions,last_real_year,age_flags(199),
                       age_flags(200),
                       age_flags(57),parest_flags(232),parest_flags(233),
                       age_flags(182));  //NMD_1Sep2015
                }
                int psi=projection_sim_index;
                int calyr=(year(ir,ip+1)-1)/age_flags(57)+1;
                int is=(year(ir,ip+1)-1)%age_flags(57)+1;
                dvariable deviate1;
                if (age_flags(182))     //NMD_7May2018
                {
                  deviate1=bh_recr_devs(simyears(psi,calyr));
                }
                else
                {
                  deviate1=bh_recr_devs(simyears(psi,year(ir,ip+1)));
                }         //NMD_7May2018

                if(age_flags(161))
                {
                  (*pnum_fish)(ir,ip+1,1)= 
                       average_recruitment_by_season(is,ir)+
                       log(srr) + deviate1 - bh_variance/2;
                }
                else
                {
                  (*pnum_fish)(ir,ip+1,1)= 
                       average_recruitment_by_season(is,ir)+
                       log(srr) + deviate1;
                }
              }
              else
              {
                tot_avrecr=totpop;
                int psi=projection_sim_index;
                int calyr=(year(ir,ip+1)-1)/age_flags(57)+1;
                int is=(year(ir,ip+1)-1)%age_flags(57)+1;
                dvariable deviate1;
                if (age_flags(182))     //NMD_7May2018
                {
                  deviate1=bh_recr_devs(simyears(psi,calyr));
                }
                else
                {
                  deviate1=bh_recr_devs(simyears(psi,year(ir,ip+1)));
                }         //NMD_7May2018
                if(age_flags(161))
                {
                  (*pnum_fish)(ir,ip+1,1)= 
                       average_recruitment_by_season(is,ir)+
                       log(srr) + deviate1 - bh_variance/2;
                }
                else
                {
                  (*pnum_fish)(ir,ip+1,1)= 
                       average_recruitment_by_season(is,ir)+
                       log(srr) + deviate1;
                }
                if (age_flags(182)==1 && age_flags(57)>1)  //NMD_27Aug2015
                {
                  cerr << "Error simulation recruitments don't support this option"
                       << endl;
                  ad_exit(1);
                }
              }    //NMD12Dec2012 
            }
            if (ng>2)
            {
              --(*pnum_fish)(ir,ip+1)(2,ng-1)=
                (*pnum_fish)(ir,ip)(1,ng-2)-(*ptm)(ir,ip)(1,ng-2);
            }

            (*pnum_fish)(ir,ip+1,ng)=log(1.e-10 
                + mfexp((*pnum_fish)(ir,ip,ng-1)-(*ptm)(ir,ip,ng-1))
                + mfexp((*pnum_fish)(ir,ip,ng)-(*ptm)(ir,ip,ng)) );

            (*pN)(ir,year(ir,ip+1))=(*pnum_fish)(ir,ip+1); //NMD19Aug2015
    
          }
          else
          {

            (*pnum_fish)(ir,ip+1,1)=(*pN)(ir,year(ir,ip+1),1);

            if (nage>2)
            --(*pnum_fish)(ir,ip+1)(2,nage-1)=
              (*pnum_fish)(ir,ip)(1,nage-2)-(*ptm)(ir,ip)(1,nage-2);

            (*pnum_fish)(ir,ip+1,nage)=
              log(1.e-10 + mfexp((*pnum_fish)(ir,ip,nage-1)-(*ptm)(ir,ip,nage-1))
                + mfexp((*pnum_fish)(ir,ip,nage)-(*ptm)(ir,ip,nage)) );

            // fit to survey samples based on ff92p
            if (allocated(have_survey_samples) && have_survey_samples(ir,ip)
                && ip<=num_real_fish_periods(ir))
            {
              dvar_vector nmid;
              nmid = (*pnum_fish)(ir,ip) -0.5*(*ptm)(ir,ip);
              for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)
              {
                int pi=parent(ir,ip,fi);
                if (ff92(pi))
                {
                  // calculate the vulnerable numbers at age
                  int ft=fish_times(ir,ip,fi);
                  vul_at_age(pi,ft)=exp(incident_sel(ir,ip,fi)+nmid);
                }
              }   
            }

            (*pN)(ir,year(ir,ip+1))=(*pnum_fish)(ir,ip+1);

          }
        }
        if (move_flags(ir,ip+1))
        {
          tmp_ip(ir)=ip+1;
          tmp_yr(ir)=year(ir,ip+1);
          tmp_mn(ir)=month(ir,ip+1);
          tmp_mp(ir)=move_index(ir,ip+1);
          
          if (ir==num_regions)
          {
            int diff_flag=0;
            for (int i=2;i<=num_regions;i++)
            {
              if (tmp_mp(i) !=tmp_mp(1))
                diff_flag=1;
            }

            if (diff_flag)
            {
              cout << "sanity error" << endl;
              cout << "map = " << tmp_mp(1,ir) << endl;
              cout << "year = " << tmp_yr(1,ir) << endl;
              cout << "month = " << tmp_mn(1,ir) << endl;
              ad_exit(1);
            }
          }
          break_flag=1;
        }
        else if (year(ir,ip+1)> year(ir,ip))
        {
          if ( af170q0==1 && age_flags(171)==1)
          {
            break_flag=1;
          }
        }
        ip++;
      }
      while (!break_flag); 
    }
    // move the fish
    if ( af170q0==1 && age_flags(171)==1
         && rip(1)<num_fish_periods(1) &&  year(1,rip(1)) <= last_real_year)
    {
      if (year(1,rip(1))-year_save>=2)
      {
        cerr << "This can't happen" << endl;
        ad_exit(1);
      }
      if (year_fill>=vtmp(1).indexmax())
      {
        cerr << "Warning year_fill would be out of bounds in lesmatrix"
             << endl;
        cout << "year_fill: " << year_fill << "max: " << vtmp(1).indexmax() << endl;
        year_fill=year_save;
      }
      else
      {
        year_fill=year_save+1;
      }
      if (year(1,rip(1))-year_save==1)
      {
//      Deal with annual recruitments and multi-species cases
        if ((age_flags(94)==3  || age_flags(182))
             && get_annual_recruitment_flag)  //NMD 7jun2022
        {
          if (!pmsd)  //- single species case   NMD_13Apr2021
          {
            dvariable tmp=get_bh_annual_recruitment_multiplier(N_q0,
                            year_save+1);
            vtmp(1,year_fill)=tmp;
          }
          else if (!sum(column(pmsd->species_flags,2))) // - multi-species
          {
            for (int is=1;is<=pmsd->num_species;is++)
            {
              pmsd->current_species=is;
              dvariable tmp=get_bh_annual_recruitment_multiplier(N_q0,
                              year_save+1);
              vtmp(is,year_fill)=tmp;
            }
          }
          else        // multi-sex
          {
            ivector sf2=column(pmsd->species_flags,2);
            int fem=0;
            for (int is=1;is<=pmsd->num_species;is++)
            {
              if (sf2(is))   // get female
              {
                fem=is;
                pmsd->current_species=is;
                dvariable tmp=get_bh_annual_recruitment_multiplier(N_q0,
                                year_save+1);
                vtmp(is,year_fill)=tmp;
              }
            }
            for (int is=1;is<=pmsd->num_species;is++)
            {
              if (!sf2(is))   // get male
              {
                pmsd->current_species=is;
                vtmp(is,year_fill)=vtmp(fem,year_fill);  // share recr.mult.
              }
            }
          }    //NMD_13Apr2021
        }
        else if (get_annual_recruitment_flag)    //NMD_7jun2022
        {
          if (!pmsd)  //NMD_13Apr2021
          {
            dvariable tmp=get_bh_recruitment_multiplier(N_q0,
                            year_save+1);
            vtmp(1,year_fill)=tmp;
          }
          else if (!sum(column(pmsd->species_flags,2))) // - multi-species
          {
            for (int is=1;is<=pmsd->num_species;is++)
            {
              pmsd->current_species=is;
              dvariable tmp=get_bh_recruitment_multiplier(N_q0,
                              year_save+1);
              vtmp(is,year_fill)=tmp;
            }
          }
          else        // multi-sex
          {
            ivector sf2=column(pmsd->species_flags,2);
            int fem=0;
            for (int is=1;is<=pmsd->num_species;is++)
            {
              if (sf2(is))   // get female
              {
                fem=is;
                pmsd->current_species=is;
                dvariable tmp=get_bh_recruitment_multiplier(N_q0,
                                year_save+1);
                vtmp(is,year_fill)=tmp;
              }
            }
            for (int is=1;is<=pmsd->num_species;is++)
            {
              if (!sf2(is))   // get male
              {
                pmsd->current_species=is;
                vtmp(is,year_fill)=vtmp(fem,year_fill);  // share recr.mult.
              }
            }
          }    //NMD_13Apr2021
        }

        for (ir=1;ir<=num_regions;ir++)
        {
          if (!pmsd)  //NMD_13Apr2021
          {
            (*pnum_fish)(ir,rip(ir),1)=vtmp(1,year_fill)
              +(*pN)(ir,year(1,rip(1)),1);
          }
          else
          {
            int is=pmsd->region_species_pointer(ir);
            (*pnum_fish)(ir,rip(ir),1)=vtmp(is,year_fill)
              +(*pN)(ir,year(1,rip(1)),1);
          }
        }

        for (ir=1;ir<=num_regions;ir++)
        {
          (*pN)(ir,year(1,rip(1)))=(*pnum_fish)(ir,rip(ir));
        }
      }
    }

    {
      if (num_regions>1)
      {
        if (no_movement_flag==0)
        {
          int tp=tmp_mp(1);
          int ir;
          check_sanity(tmp_mp);
          dvar_matrix tmp;
          if (!pmsd) 
            tmp=fast_diffusion_calcs(nage,num_regions,
             (*pnum_fish),Dad(tp),rip);
          else
            tmp=fast_diffusion_calcs(nage,num_regions,
              (*pnum_fish),Dad(tmp_mp(1)),rip,0,pmsd);
          for (ir=1;ir<=num_regions;ir++)
          {
            (*pnum_fish)(ir,rip(ir))=log(5.e-10+tmp(ir));
          }
        } 
      }
    }
  } // need to decide when to quit
  while (!finished_flag);
  
  do_fish_mort_intermediate_calcs();

  if (parest_flags(360)>0)
  {
    do_tag_loss_fish_mort_intermediate_calcs();
  }

  if (do_fishery_projections_flag==1 )
  {
    do_fish_mort_intermediate_projection_calcs();
  }

  dvar4_array * pc =0;
  dvar4_array * pfmc =0;
  if (af170q0==0)
  {
    pc = &catch;
    pfmc = &fish_mort_calcs;
  }
  else
  {
    pc = &catch_q0;
    pfmc = &fish_mort_calcs_q0;
  }
  
  for (ir=1;ir<=num_regions;ir++)
  {
    int tmp_nfp;
    if (do_fishery_projections_flag==0)
    {
      tmp_nfp=num_real_fish_periods(ir);  
    }
    else
    {
      tmp_nfp=num_fish_periods(ir);  
    }
   
    for (int ip=1;ip<=tmp_nfp;ip++)  
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {        
        (*pc)(ir,ip,fi)=(*pfmc)(ir,ip,fi)+(*pnum_fish)(ir,ip);
      }
    }
  }
 
 /*
  ofstream ofs("www");
  ofs << setprecision(15) << exp(*pnum_fish) << endl;
  ad_exit(1);

  int ip;
  for (ip=1;ip<=10;ip++)  
  {
    for (ir=1;ir<=num_regions;ir++)
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {
        if (parent(ir,ip,fi)==1)
        {
          cout << endl << "ip =" << ip << " ir = " << ir << endl;
          cout << sum(exp((*pc)(ir,ip,fi))) << " " 
              << (*pc)(ir,ip,fi) << " " << (*pfmc)(ir,ip,fi) << endl;
        }
      }
    }
  }
  ad_exit(1);
 */
 cout << "PPPPP   Tot mort penalty " << fpen1 << "  Zmax_flag = "
      << Zmax_flag << endl;
 ffpen+=fpen1;
 return vtmp;
}

void dvar_len_fish_stock_history::zero_out_catches(void)
{
  // zero out the catches from fisheries that are to be turned off;
  ivector ff55=column(fish_flags,55);
  
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++) 
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {        
        if (ff55(parent(ir,ip,fi)) && obs_tot_catch(ir,ip,fi)>-0.5L)
        {
          obs_tot_catch(ir,ip,fi)=1.e-3;
        }
      }
    }
  }
}

void dvar_len_fish_stock_history::zero_out_log_effort(void)
{
  // zero out the effort from fisheries that are to be turned off;
  ivector ff55=column(fish_flags,55);
  
  for (int fi=1;fi<=num_fisheries;fi++)
  {
    for (int nt=1;nt<=num_fish_times(fi);nt++) 
    {
      if (ff55(fi))
      {
        log_effort_by_fishery(fi,nt)=-20.0;
      }
    }
  }
}
