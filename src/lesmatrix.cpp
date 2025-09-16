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
#endif

extern dvar_len_fish_stock_history * pcfsh;
ofstream pofs_recrdbg("recruitment_debug");

void sanity_check(ivector& v,imatrix& year)
{  
  int mmin=v.indexmin();
  int mmax=v.indexmax();

  int yr=year(1,v(1));
  for (int i=2;i<=mmax;i++)
  {
    if (year(i,v(i)) != yr)
    {
      cerr << "year sanity error  yr = " << yr << "y(" << i << ")=" 
           << year(i,v(i)) << endl;
      cerr << " ip " << v << endl;
      ad_exit(1);
    }
  }
}      

void print_movement_stuff(ofstream& xpofs, ivector& tmp_mp, ivector& tmp_yr,
  ivector& tmp_mn,ivector& tmp_ip,int num_regions)
{
  xpofs << "map = " << tmp_mp(1,num_regions) << endl;
  xpofs << "year = " << tmp_yr(1,num_regions) << endl;
  xpofs << "month = " << tmp_mn(1,num_regions) << endl;
  xpofs << "periods = " << tmp_ip(1,num_regions) << endl;
  xpofs << endl;
}
int cobb_douglas_flag=0;

dvar_matrix dvar_len_fish_stock_history::get_numbers_at_age(dvar_vector& sv,
  void * pq_flag)
{
  ofstream ofss1("biomasses");
  if (age_flags(125)>0)
    cobb_douglas_flag=1;

  Zmax_flag=0;    //NMD_18jan2022
  Zmax_fish=1.0;
  if (age_flags(116)!=0)
  {
    Zmax_fish=age_flags(116)/100.;
  }
  
  dvariable ffpen1=0.0;
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
  int year_change_flag=0; dvar3_array * pnum_fish=0; dvar3_array * pN=0;
  dvar3_array * ptm=0;
  
  if (af170q0==0)
  {
    pnum_fish=&num_fish; ptm=&tot_mort; pN=&N;
  }
  else
  {
    ptm=&tot_mort_q0; pnum_fish=&num_fish_q0; pN=&N_q0;
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
  dvar_matrix vtmp;
  if (!pmsd)
  {
    vtmp.allocate(1,1,1,last_real_year+1);   //NMD_13Apr2021
  }
  else
  {
    vtmp.allocate(1,pmsd->num_species,1,last_real_year+1);   //NMD_13Apr2021
  }
  vtmp.initialize();
  int ns=age_flags(57);
  int af94=age_flags(94);
  if (af94!=3)
  {
    for (ir=1;ir<=num_regions;ir++)
    {
      (*pnum_fish)(ir,1)=(*pN)(ir,1);
    }
  }
  else
  {
    for (ir=1;ir<=num_regions;ir++)
    {
      for (int is=1;is<=ns;is++)
      {
        (*pnum_fish)(ir,is)=(*pN)(ir,is);
      }
    }
  }
 
  dvariable alpha=0.0;
  dvariable beta=0.0;
  int have_ab_flag=0;

  int no_movement_flag=0;
  int year_save=0;
  int year_fill=0;
  random_number_generator * rng2=0;
  if (age_flags(196)>0 && projection_sim_index>0)
  {
    rng2=new random_number_generator(age_flags(197)+projection_sim_index);
  }
  dvar_vector lbio1;
  dvar_vector lcurrent_rbio;
  

  if (cobb_douglas_flag)
  {
    current_biomass_by_year.initialize();
    lbio1=log(1.e-20+pcfsh->calculate_the_current_biomass(rip));
    lcurrent_rbio.allocate(1,num_regions);
    lcurrent_rbio.initialize();
    ofss1 << exp(lbio1)  << endl;
  }

  do
  {
    finished_flag=1;
    year_save=year(1,rip(1));
    no_movement_flag=0;
    if ( af170q0==1 && age_flags(171)==1
         && year(1,rip(1)) <= last_real_year )
    {

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
      if (age_flags(94)==3  || age_flags(182) )
      {
        if (get_annual_recruitment_flag)
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
      }
      else
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
    }

    for (int ir=1;ir<=num_regions;ir++)
    {
      int ng=nage_by_region(ir);
      break_flag=0;
      int& ip=rip(ir);
      do
      {
        if (cobb_douglas_flag)
        {
          if (af170q0==0)
          {
            get_fishing_and_total_mortality_for_this_period
              (ir,ip,lcurrent_rbio);
            do_fish_mort_intermediate_calcs(ir,ip);
          }
          else
          {
            get_fishing_and_total_mortality_for_this_period_q0
              (ir,ip,lcurrent_rbio);
            do_fish_mort_intermediate_calcs_q0(ir,ip);
          }
        }

        dvariable srr;   //NMD 15Mar2012
        if ( ip==num_real_fish_periods(ir) && do_fishery_projections_flag==1 )
        {
          if (have_ab_flag==0)
          {
            if(!age_flags(190) || age_flags(195)>0)   //NMD_25Sep2018
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
          
          if (age_flags(182)>0 && age_flags(57)>1 ) 
          {
            r/=age_flags(57);   //NMD25Aug2015
          } 
          pofs_recrdbg << ir << " " << ip+1 << " " << year(ir,ip+1)  
            << " " << setscientific() << r << endl;
          srr = r;     //NMD 15Mar2012
          if ( age_flags(190)<=0) // Use SRR for projected recruitment
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
        dvariable ffpen=0.0;
        if (ip>num_real_fish_periods(ir) && do_fishery_projections_flag==1)
        {
          if (num_fish_incidents(ir,ip)>0)
          {
            if (missing_catch_for_period_flag(ir,ip)==0 )
            {
              do_newton_raphson_with_totcatch(ir,ip,ffpen,
                fm_level(ir,ip),ffpen1);
            }
            else
            {
              do_newton_raphson_missing_totcatch(ir,ip,ffpen,
               fm_level(ir,ip));
            }
          }
        }
        if ( ip>=num_fish_periods(ir) ||
          ( ip>num_real_fish_periods(ir)  && do_fishery_projections_flag==0 ) )
        {
          no_movement_flag=1; // we are at thenb end so ip does not
                              // get incremented so don't move the fish
                              // df july 13 05
          break;
        }
    
        finished_flag=0;
        if (year(ir,ip+1)==year(ir,ip))
        {
          if (af94!=3 || year(ir,ip)>=ns)
            (*pnum_fish)(ir,ip+1)=(*pnum_fish)(ir,ip)-(*ptm)(ir,ip);
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
              r/=age_flags(57);   //NMD25Aug2015
            srr = r;     //NMD 15Mar2012
            pofs_recrdbg << ir << " " << ip+1 << " " << year(ir,ip+1)  
              << " " << setscientific() << r << endl;
            if ( age_flags(190)<=0) // Use SRR for projected recruitment
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
                if(age_flags(182) == 0)     //NMD 05Oct2022
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
//                (*pN)(ir,year(ir,ip+1),1)=
//                  simulated_recruitments(projection_sim_index,ir,year(ir,ip+1));
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
          if (af170q0 && year(ir,ip)<=last_real_year &&
            get_annual_recruitment_flag) //NMD_9Feb2015
          {
            if (year_fill && get_annual_recruitment_flag)
            {
              if (year_fill != year(ir,ip+1))
              {
                 cerr << "rethink" << year_fill << " " << year(ir,ip+1)
                      << endl;
                 (*pnum_fish)(ir,ip+1,1)=Rsave(ir,year(ir,ip+1));
              }
              else
              {
                if (!pmsd)  //NMD_13Apr2021
                {
                  (*pnum_fish)(ir,ip+1,1)=vtmp(1,year_fill)
                    +Rsave(ir,year(ir,ip+1));
                }
                else
                {
                  int is=pmsd->region_species_pointer(ir);
                  (*pnum_fish)(ir,ip+1,1)=vtmp(is,year_fill)
                    +Rsave(ir,year(ir,ip+1));		  
                }  //NMD_13Apr2021
              }
            } 
          }
          else
          {
            (*pnum_fish)(ir,ip+1,1)=(*pN)(ir,year(ir,ip+1),1);
            Rsave(ir,year(ir,ip+1))=(*pN)(ir,year(ir,ip+1),1);
          }          

          // age the fish
          if (ip==num_real_fish_periods(ir) 
            && do_fishery_projections_flag==1
            && projection_sim_index>0)
          {
            if(parest_flags(239) == 0)
            {     //NMD 15Mar2012
              if(age_flags(182) == 0)     //NMD_05Oct2022
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
              }     //NMD_05Oct2022
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

            if(parest_flags(239) == 0){     //NMD 15Mar2012
              if(age_flags(182) == 0)     //NMD_05Oct2022
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
              }     //NMD_05Oct2022
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
            if (af94!=3 || year(ir,ip)>=ns)
            {
              if (af170q0==0 || !age_flags(171))
                (*pnum_fish)(ir,ip+1,1)=(*pN)(ir,year(ir,ip+1),1);

              if (ng>2)
                --(*pnum_fish)(ir,ip+1)(2,ng-1)=
                  (*pnum_fish)(ir,ip)(1,ng-2)-(*ptm)(ir,ip)(1,ng-2);
      
              (*pnum_fish)(ir,ip+1,ng)=
                log(1.e-10 + mfexp((*pnum_fish)(ir,ip,ng-1)-(*ptm)(ir,ip,ng-1))
                  + mfexp((*pnum_fish)(ir,ip,ng)-(*ptm)(ir,ip,ng)) );

              (*pN)(ir,year(ir,ip+1))=(*pnum_fish)(ir,ip+1);
            }
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
         && year(1,rip(1)) <= last_real_year )
    {
      if (year(1,rip(1))-year_save>=2)
      {
        cerr << "This can't happen" << endl;
        ad_exit(1);
      }
    }

    {
      if (num_regions>1)
      {
        if (af94!=3 || year(1,rip(1))>ns)
        {
          sanity_check(rip,year);
          if (no_movement_flag==0)
          {
            check_sanity(tmp_mp);
            dvar_matrix tmp=fast_diffusion_calcs(nage,num_regions,(*pnum_fish),
              Dad(tmp_mp(1)),rip,0,pmsd);
            for (int ir=1;ir<=num_regions;ir++)
            {
              (*pnum_fish)(ir,rip(ir))=log(5.e-10+tmp(ir));
            }
          }
        } 
      }
    }
    if (cobb_douglas_flag)
    {
      dvar_vector current_biomass=pcfsh->calculate_the_current_biomass(rip);
      lcurrent_rbio=log(current_biomass)-lbio1;
      if (year(1,rip(1))>year(1,rip(1)-1))
      {
        current_biomass_by_year(year(1,rip(1)))=lcurrent_rbio;
      }
      ofss1 << current_biomass << endl;
    }
  } // need to decide when to quit
  while (!finished_flag);

  if (parest_flags(360)>0)  //NMD_18Feb2022
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
  if (do_fishery_projections_flag)
  {
    ofstream ofs("projpop");
    for (ir=1;ir<=num_regions;ir++)
    {
      int nrfp=num_real_fish_periods(ir);
      int nfp=num_fish_periods(ir);  
      if (nfp>nrfp)
      {
        int ip=nrfp+1;
        int fi;
        for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {        
          ofs << setw(11) << (*pfmc)(ir,ip,fi) << endl;;
        }
        ofs << "catch" << endl;
        for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {        
          ofs << setw(11) << exp((*pc)(ir,ip,fi)) << endl;;
        }
        ofs << "num_fish" << endl;
        ofs << setw(11) << exp((*pnum_fish)(ir,ip)) << endl;;
        ofs << "fish mort" << endl;
        for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {        
          ofs << setw(11) << exp(fish_mort(ir,ip,fi)) << endl;;
        }
      }
    }
  }
 
  if (rng2)
  {
    delete rng2;
    rng2=0;
  } 
  return vtmp;
}


dvar_vector get_average_recruitment(dvar3_array& N,
  int num_regions,int last_real_year)
{
  dvar_vector avg_recruitment(1,num_regions);
  avg_recruitment.initialize();
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int iy=1;iy<=last_real_year;iy++)
    {
      avg_recruitment(ir)+=mfexp(N(ir,iy,1));
    }
  }
  avg_recruitment/=last_real_year;
  return avg_recruitment;
}   

dvar_matrix get_average_recruitment_by_season(dvar3_array& N,
  int num_regions,int last_real_year,int af57)
{
  dvar_matrix avg_recruitment_by_season(1,af57,1,num_regions);
  avg_recruitment_by_season.initialize();
  dmatrix nis(1,af57,1,num_regions);
  nis.initialize();
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int iy=1;iy<=last_real_year;iy++)
    {
      int is=(iy-1)%af57+1;
      avg_recruitment_by_season(is,ir)+=mfexp(N(ir,iy,1));
      nis(is,ir)+=1.0;
    }
  }
  avg_recruitment_by_season=elem_div(avg_recruitment_by_season,nis);
  return avg_recruitment_by_season;
}   

dvar_matrix get_average_recruitment_by_season_for_region(dvar_matrix& average_recruitment_by_season,
  int num_regions,int af57)
{
  dvar_matrix avg_recruitment(1,af57,1,num_regions);
  avg_recruitment.initialize();
  for (int ir=1;ir<=num_regions;ir++)
  {
    dvariable tmp=0;
    for (int is=1;is<=af57;is++)
    {
      tmp+=mfexp(average_recruitment_by_season(is,ir));
    }
    for (int is=1;is<=af57;is++)
    {
      avg_recruitment(is,ir)=mfexp(average_recruitment_by_season(is,ir))/tmp;
    }
  }
  avg_recruitment=log(1.e-20+avg_recruitment);
  return avg_recruitment;
}   

dvar_matrix get_scalar_recruitment_by_season_for_region
  (dvar_matrix& average_recruitment_by_season_for_regions,
  int num_regions,int af57)
{
  dvar_matrix avg_recruitment(1,af57,1,num_regions);
  avg_recruitment.initialize();
  for (int ir=1;ir<=num_regions;ir++)
  {
    dvariable tmp=0;
    for (int is=1;is<=af57;is++)
    {
      tmp+=mfexp(average_recruitment_by_season_for_regions(is,ir));
    }
    tmp/=af57;
    for (int is=1;is<=af57;is++)
    {
      avg_recruitment(is,ir)=
	mfexp(average_recruitment_by_season_for_regions(is,ir))/tmp;
    }
  }
  avg_recruitment=log(1.e-20+avg_recruitment);
  return avg_recruitment;
}   


dvar_vector get_average_recruitment_for_projections(dvar3_array& N,
  int num_regions,int last_real_year, int imin, int imax)
{
  dvar_vector avg_recruitment_for_projections(1,num_regions);
  avg_recruitment_for_projections.initialize();
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int iy=last_real_year-imin+1;iy<=last_real_year-imax;iy++)
    {
      avg_recruitment_for_projections(ir)+=mfexp(N(ir,iy,1));
    }
  }
  avg_recruitment_for_projections/=(imin-imax);
  return avg_recruitment_for_projections;
}   

//NMD13Dec2012
//dvariable get_average_log_recruitment_for_projections(dvar3_array& N,
//  int num_regions, int last_real_year, int imin, int imax)
//NMD_1Sep2015
dvariable get_average_log_recruitment_for_projections(dvar3_array& N,
  int num_regions, int last_real_year, int imin, int imax, int af57,
  int pf232, int pf233, int af182)
{
  dvariable avg_log_recruitment_for_projections;
  avg_log_recruitment_for_projections=0.0;

  if (af182==1 && af57>1)  //NMD_1Sep2015
  {
    int ns=af57;
    int rem=last_real_year%ns;
    int na=last_real_year/ns;
    if (rem) na++;
    dvar_matrix Nann(1,num_regions,1,na);
    for (int ir=1;ir<=num_regions;ir++)
    {
      int ii=0;
      for (int i=1;i<=na;i++)
      {
        Nann(ir,i)=0.0;
        for (int j=1;j<=ns;j++)
        {
          ++ii;
          Nann(ir,i)+=mfexp(N(ir,ii,1));
        }
      }
    }
    int yrfirst=pf232/af57+1;    //first calendar year of BH-SRR period
    int yrlast=(pf233-1)/af57+1;       //last calendar year of BH-SRR period
    for (int i=yrfirst;i<=yrlast;i++)
    {
      dvariable sum=0.0;
      for (int ir=1;ir<=num_regions;ir++)
      {
        sum+=Nann(ir,i);
      }
     // !!!!!!!!XXXX
      avg_log_recruitment_for_projections+=log(sum/num_regions);
    }
    avg_log_recruitment_for_projections/=(yrlast-yrfirst+1);
  }
  else
  {
    for (int iy=last_real_year-imin+1;iy<=last_real_year-imax;iy++)
    {
      dvariable tmp;
      tmp=0.0;
      for (int ir=1;ir<=num_regions;ir++)
      {
        tmp+=mfexp(N(ir,iy,1));
      }
      avg_log_recruitment_for_projections+=log(tmp);
    }
    avg_log_recruitment_for_projections=avg_log_recruitment_for_projections/(imin-imax);
  }
  return avg_log_recruitment_for_projections;
}   


