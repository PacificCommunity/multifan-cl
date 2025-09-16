/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

void  dvar_len_fish_stock_history::do_simulation_stuff(void)
{
  remove("test_lw_sim_alt");
  remove("test_lw_sim");
  remove("effort_sim");
  remove("effort_sim_true");  //NMD_29jun2020
  remove("catch_sim");
  remove("cpue_sim");  //NMD_31aug2023
  remove("cpue_sim_true");  //NMD_30jul2024
  set_option_flag("-length_seed",simulation_seeds(1),ad_argc,ad_argv);
  set_option_flag("-weight_seed",simulation_seeds(2),ad_argc,ad_argv);
  set_option_flag("-effort_seed",simulation_seeds(3),ad_argc,ad_argv);
  set_option_flag("-catch_seed",simulation_seeds(4),ad_argc,ad_argv);
  set_option_flag("-tag_seed",simulation_seeds(5),ad_argc,ad_argv);
  set_option_flag("-cpue_seed",simulation_seeds(6),ad_argc,ad_argv);  //NMD_31aug2023

  read_tag_simulation_info();
}


void  dvar_len_fish_stock_history::read_tag_simulation_info(void)
{
  adstring tmpstring=ad_root+adstring(".tag_sim");
  cifstream cifs((char *)(tmpstring));

  if (!cifs)
  {
    *clogf << "Error trying to open file " << tmpstring << endl;
    cerr << "Error trying to open file " << tmpstring << endl;
    exit(1);
  }
  //*clogf << "BD" << endl;
  cifs >> sim_num_tag_releases;

  sim_num_tags_at_length.allocate(1,sim_num_tag_releases);
  sim_tag_region.allocate(1,sim_num_tag_releases);
  sim_true_tag_year.allocate(1,sim_num_tag_releases);
  sim_tag_year.allocate(1,sim_num_tag_releases);
  sim_tag_month.allocate(1,sim_num_tag_releases);
  sim_true_tag_month.allocate(1,sim_num_tag_releases);
  sim_tag_fishery.allocate(1,sim_num_tag_releases);
  sim_tag_numbers_released.allocate(1,sim_num_tag_releases);
  sim_initial_tag_period.allocate(1,sim_num_tag_releases,1,num_regions);
  sim_tag_incident.allocate(1,sim_num_tag_releases);
  sim_initial_tag_release_by_age.allocate(1,sim_num_tag_releases,1,nage);
  sim_initial_tag_recruitment_period.allocate(1,sim_num_tag_releases,
    1,num_regions);
  sim_tag_fish_rep.allocate(1,num_fisheries,1,sim_num_tag_releases);
  sim_tag_region.initialize();
  sim_true_tag_year.initialize();
  sim_tag_year.initialize();
  sim_tag_month.initialize();
  sim_tag_fishery.initialize();
  sim_tag_numbers_released.initialize();
  sim_initial_tag_period.initialize();
  sim_tag_incident.initialize();
  sim_initial_tag_release_by_age.initialize();

  if (pmsd)
  {
    pmsd->sim_tag_species_flag.allocate(1,sim_num_tag_releases,
      1,pmsd->num_species);
    pmsd->sim_num_tag_release_by_species.allocate(1,pmsd->num_species);
  }
  
  sim_initial_tag_year.allocate(1,sim_num_tag_releases);
  
  for (int it=1;it<=sim_num_tag_releases;it++)
  {
    cifs >> sim_tag_region(it) >> sim_tag_year(it) >> sim_tag_month(it);
    if (pmsd)
    {
      cifs >> pmsd->sim_tag_species_flag(it);
      pmsd->sim_num_tag_release_by_species+=pmsd->sim_tag_species_flag(it);
      if (sum(pmsd->sim_tag_species_flag(it))>0)
        pmsd_error();
    }
    
    sim_true_tag_year(it)=sim_tag_year(it);

    sim_true_tag_month(it)=sim_tag_month(it);



    {
      normalize_tag_dates(sim_tag_year,sim_tag_month,sim_true_tag_year,
        month_1,year1,direction_flag,month_factor,first_time,it);
    }


    cifs >> sim_tag_fishery(it);
    cifs >> sim_tag_numbers_released(it);

    int p=sim_tag_fishery(it);
    for (int nt=1;nt<=num_fish_times(p);nt++)
    {
      int rr=realization_region(p,nt);
      int rp=realization_period(p,nt);
      int ty=really_true_year(rr,rp)+year1-1;
      int tm=really_true_month(rr,rp);
      int inc=realization_incident(p,nt);
      if (ty==sim_true_tag_year(it) && tm==sim_true_tag_month(it)) 
      {
        sim_tag_year(it)=year(rr,rp);
        sim_tag_incident(it)=inc;
        break;
      }
    }
  }
  for (int it=1;it<=sim_num_tag_releases;it++)
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      int ip;
      int no_match_flag=1;
      for (ip=1;ip<=num_fish_periods(ir);ip++) 
      {
        if (sim_tag_year(it) < year(ir,ip)) 
        {
          no_match_flag=0;
          break;
        }
        if ( sim_tag_year(it) == year(ir,ip) 
           && sim_tag_month(it) <= month(ir,ip))
        {
          no_match_flag=0;
          break;
        }
      }
      if (no_match_flag)
      {
        *clogf << endl << "Sim Tag group " << it 
               << " released after last fishery "
             << "remove it from tag file " << endl << endl;
        cerr << endl << "Tag group " << it << " released after last fishery "
             << "remove it from tag file " << endl << endl;
        ad_exit(1);
      }
      sim_initial_tag_period(it,ir)=ip;
    }
  }
  if (age_flags(96))
  {
    sim_terminal_tag_period.allocate(1,sim_num_tag_releases,1,num_regions);

    for (int it=1;it<=sim_num_tag_releases;it++)
    {
      for (int ir=1;ir<=num_regions;ir++)
      {
        int ip;
        for (ip=1;ip<=num_fish_periods(ir);ip++) 
        {
          if (sim_tag_year(it)+age_flags(96) < year(ir,ip)) break;
          if ( sim_tag_year(it)+age_flags(96) == year(ir,ip) 
           && sim_tag_month(it) <= month(ir,ip)) break; 
        }
        sim_terminal_tag_period(it,ir)=min(ip,num_fish_periods(ir));
      }
    }
  }

  sim_minimum_initial_tag_period.allocate(1,num_regions);
  sim_initial_tag_year.allocate(1,sim_num_tag_releases);
  sim_initial_tag_year=sim_tag_year;
  int bad_tag_flag=0;
    
  for (int ir=1;ir<=num_regions;ir++)
  {
    sim_minimum_initial_tag_period(ir)=sim_initial_tag_period(1,ir);
    for (int it=2;it<=sim_num_tag_releases;it++)
    {
      if (sim_minimum_initial_tag_period(ir)>sim_initial_tag_period(it,ir)) 
        sim_minimum_initial_tag_period(ir)=sim_initial_tag_period(it,ir);
    }
  }

  sim_min_tag_age.allocate(1,sim_num_tag_releases,1,num_regions,
      sim_initial_tag_period,
      imatrix(1,sim_num_tag_releases,num_fish_periods));

    
  for (int it=1;it<=sim_num_tag_releases;it++)
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      int min_age=1;
      sim_min_tag_age(it,ir,sim_initial_tag_period(it,ir))=min_age;
      for (int ip=sim_initial_tag_period(it,ir);ip<=num_fish_periods(ir);ip++)
      {
        sim_min_tag_age(it,ir,ip)=
          min(nage,min_age+year(ir,ip)-sim_initial_tag_year(it));
      }
    }
  }

  sim_min_tag_age1.allocate(1,sim_num_tag_releases,1,num_regions,
      index_type(sim_initial_tag_year),nyears);

  for (int it=1;it<=sim_num_tag_releases;it++)
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      int min_age=1;
      sim_min_tag_age1(it,ir,sim_initial_tag_year(it))=min_age;
      for (int iy=sim_initial_tag_year(it)+1;iy<=nyears;iy++)
      {
        if (min_age<nage) min_age++;
        sim_min_tag_age1(it,ir,iy)=min_age;
      }
    }
  }
 
    
  sanity_check_2(sim_initial_tag_year,year,sim_initial_tag_period);
    
  if (!age_flags(96)) 
  {
    cout << " sim_initial_tag_year(1) "  << sim_initial_tag_year(1)  << endl;
    sim_tagN.allocate(1,sim_num_tag_releases,1,num_regions,
      index_type(sim_initial_tag_year),nyears,sim_min_tag_age1,nage);
    sim_tagnum_fish.allocate(1,sim_num_tag_releases,1,num_regions,
      sim_initial_tag_period,
      imatrix(1,sim_num_tag_releases,num_fish_periods),sim_min_tag_age,nage);
  }
  else
  {
    int it;
    sim_terminal_tag_year.allocate(1,sim_num_tag_releases);
    sim_terminal_tag_year=sim_initial_tag_year+age_flags(96);
    for (int it=1;it<=sim_num_tag_releases;it++)
    {
      if (sim_terminal_tag_year(it)>nyears)
        sim_terminal_tag_year(it)=nyears;
      if (sim_terminal_tag_year(it)<sim_initial_tag_year(it))
      {
        cout << "Error sim_terminal tag year for sim tag group " << it
             << " is less than the sim_initial tag year " << endl;
        ad_exit(1);
      }
    }


    sim_min_tag_age5.allocate(1,sim_num_tag_releases,1,num_regions,
      sim_initial_tag_period,sim_terminal_tag_period);

    for (int it=1;it<=sim_num_tag_releases;it++)
    {
      for (int ir=1;ir<=num_regions;ir++)
      {
        int min_age=1;
        sim_min_tag_age5(it,ir,sim_initial_tag_period(it,ir))=min_age;
        for (int ip=sim_initial_tag_period(it,ir)+1;
          ip<=sim_terminal_tag_period(it,ir);ip++)
        {
          if (year(ir,ip)>year(ir,ip-1)) 
          {
            if (min_age<nage) min_age++;
          }
          sim_min_tag_age5(it,ir,ip)=min_age;
        }
      }
    }
    sim_min_tag_age6.allocate(1,sim_num_tag_releases,1,num_regions,
      sim_initial_tag_year,sim_terminal_tag_year);

    for (int it=1;it<=sim_num_tag_releases;it++)
    {
      for (int ir=1;ir<=num_regions;ir++)
      {
        int min_age=1;
        sim_min_tag_age6(it,ir,sim_initial_tag_year(it))=min_age;
        for (int ip=sim_initial_tag_year(it)+1;
          ip<=sim_terminal_tag_year(it);ip++)
        {
          if (year(ir,ip)>year(ir,ip-1)) 
          {
            if (min_age<nage) min_age++;
          }
          sim_min_tag_age6(it,ir,ip)=min_age;
        }
      }
    }

  *clogf << "BM" << endl;

    sim_tagN.allocate(1,sim_num_tag_releases,1,num_regions,
      index_type(sim_initial_tag_year),index_type(sim_terminal_tag_year),
      sim_min_tag_age6,nage);
    sim_tagnum_fish.allocate(1,sim_num_tag_releases,1,num_regions,
      sim_initial_tag_period,sim_terminal_tag_period,sim_min_tag_age5,nage);
    sim_minttp.allocate(1,num_regions);
    int ir;
    for (ir=1;ir<=num_regions;ir++)
    {
      sim_minttp(ir)=sim_terminal_tag_period(1,ir);
      for (int it=2;it<=sim_num_tag_releases;it++)
      {
        if (sim_minttp(ir)>sim_terminal_tag_period(it,ir))
          sim_minttp(ir)=sim_terminal_tag_period(it,ir);
      }
    }
    sim_min_tag_year=min(sim_tag_year);
    sim_pooledtagN.allocate(1,num_regions,sim_min_tag_year,nyears,1,nage);

    sim_pooled_tagnum_fish.allocate(1,num_regions,sim_minttp+1,num_fish_periods,
       1,nage);
    sim_epooled_tagnum_fish_recr.allocate(1,num_regions,sim_minttp+1,
       num_fish_periods,1,nage);
  }

  sim_min_init_tag_period.allocate(1,num_regions);
  for (int ir=1;ir<=num_regions;ir++)
  {
    sim_min_init_tag_period(ir)=sim_initial_tag_period(1,ir);
    for (int it=2;it<=sim_num_tag_releases;it++)
    {
      if (sim_min_init_tag_period(ir)>sim_initial_tag_period(it,ir))
        sim_min_init_tag_period(ir)=sim_initial_tag_period(it,ir);
    }
  }


  if (!age_flags(96)) 
  {
    sim_num_alltagfish_incidents.allocate(1,sim_num_tag_releases,1,num_regions,
      sim_initial_tag_period,imatrix(1,sim_num_tag_releases,num_fish_periods));
    sim_num_tagfish_incidents.allocate(1,sim_num_tag_releases,1,num_regions,
      sim_initial_tag_period,imatrix(1,sim_num_tag_releases,num_fish_periods));

  }
  else
  {
    sim_num_alltagfish_incidents.allocate(1,sim_num_tag_releases,1,num_regions,
      sim_initial_tag_period,sim_terminal_tag_period);
    sim_num_tagfish_incidents.allocate(1,sim_num_tag_releases,1,num_regions,
      sim_initial_tag_period,sim_terminal_tag_period);
    sim_num_pooledtagfish_incidents.allocate(1,num_regions,
      sim_minttp+1,num_fish_periods);
  }

  *clogf << "BO" << endl;
  for (int it=1;it<=sim_num_tag_releases;it++)
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      int ip;
         
      if (sim_initial_tag_period(it,ir))
      {
        if (!age_flags(96)) 
          for (ip=sim_initial_tag_period(it,ir);ip<=num_fish_periods(ir);ip++)
          {
            sim_num_alltagfish_incidents(it,ir,ip)=num_fish_incidents(ir,ip);
            sim_num_tagfish_incidents(it,ir,ip)=num_fish_incidents(ir,ip);
          }
          else
            for (ip=sim_initial_tag_period(it,ir);ip<=sim_terminal_tag_period(it,ir);ip++)
            {
              sim_num_alltagfish_incidents(it,ir,ip)=num_fish_incidents(ir,ip);
              sim_num_tagfish_incidents(it,ir,ip)=num_fish_incidents(ir,ip);
            }
        }
      }
    }
    if (age_flags(96)) 
      for (int ir=1;ir<=num_regions;ir++)
        for (int ip=sim_minttp(ir)+1;ip<=num_fish_periods(ir);ip++)
          sim_num_pooledtagfish_incidents(ir,ip)=num_fish_incidents(ir,ip);

    get_sim_initial_tag_recruitment_period();

  *clogf << "BP" << endl;

  if (!age_flags(96)) 
  {
    sim_tot_tag_catch.allocate(1,sim_num_tag_releases,1,num_regions,
      sim_initial_tag_period,imatrix(1,sim_num_tag_releases,num_fish_periods),
      1,sim_num_alltagfish_incidents);

    {
      ofstream ofs("sim_obstag");
      ofs << "sim_num_tag_releases" << endl;
      ofs << sim_num_tag_releases << endl;
      ofs << "num_regions" << endl;
      ofs << num_regions << endl;
      ofs << "num_fish_periods" << endl;
      ofs << num_fish_periods << endl;
      ofs << "sim_num_alltagfish_incidents" << endl;
      ofs << sim_num_alltagfish_incidents << endl;
    }
  }
  else
  {
    sim_tot_tag_catch.allocate(1,sim_num_tag_releases,1,num_regions,
      sim_initial_tag_period,sim_terminal_tag_period,1,
      sim_num_alltagfish_incidents);
  }


  sim_initial_tag_release_by_age.allocate(1,sim_num_tag_releases,1,nage);


  if (age_flags(96)) 
  {
    sim_min_tag_age4.allocate(1,sim_num_tag_releases,1,num_regions,
       sim_initial_tag_period,sim_terminal_tag_period,1,
       sim_num_tagfish_incidents);
  
    for (int it=1;it<=sim_num_tag_releases;it++)
    {
      for (int ir=1;ir<=num_regions;ir++)
      {
        int min_age=1;
        sim_min_tag_age4(it,ir,sim_initial_tag_period(it,ir))=min_age;
        for (int ip=sim_initial_tag_period(it,ir)+1;
          ip<=sim_terminal_tag_period(it,ir);ip++)
        {
          if (year(ir,ip)>year(ir,ip-1)) 
          {
            if (min_age<nage) min_age++;
          }
          for (int fi=1;fi<=sim_num_tagfish_incidents(it,ir,ip);fi++)
          {
            sim_min_tag_age4(it,ir,ip,fi)=min_age;
          }
        }
      }
    }
    
  *clogf << "BQ" << endl;

    sim_tagcatch.allocate(1,sim_num_tag_releases,1,num_regions,
      sim_initial_tag_period,sim_terminal_tag_period,
      1,sim_num_tagfish_incidents,sim_min_tag_age4,nage);

    sim_obs_tagcatch.allocate(1,sim_num_tag_releases,1,num_regions,
      sim_initial_tag_period,sim_terminal_tag_period,
      1,sim_num_tagfish_incidents,sim_min_tag_age4,nage);

    sim_pooled_tagcatch.allocate(1,num_regions,sim_minttp+1,
      num_fish_periods,
      1,sim_num_pooledtagfish_incidents,1,nage);

    sim_pooledtot_tag_catch.allocate(1,num_regions,sim_minttp+1,
      num_fish_periods,
      1,sim_num_pooledtagfish_incidents);

    sim_tagcatch.initialize();

    sim_pooled_tagcatch.initialize();

    sim_pooledtot_tag_catch.initialize();
  }
  cout << "Reading sim_tag_fish_rep" << endl;
  cifs >> sim_tag_fish_rep;
  if (!cifs) 
  {
    cerr << "Error reading sim_tag_fish_rep" << endl;
    ad_exit(1);
  }
  *clogf << "BR" << endl;
}

void  dvar_fish_stock_history::sim_get_initial_tag_population(void)
{
  for (int it=1;it<=sim_num_tag_releases;it++)
  {
    int ir=sim_tag_region(it);
    int ip=sim_initial_tag_period(it,ir);
    int fi=sim_tag_incident(it);
    int yr1=sim_tag_year(it);

    //    sim_initial_tag_release_by_age(it)=
    //      prop(ir,ip,fi)*sim_tag_numbers_released(it);
    dvar_vector ecatch=exp(catch(ir,ip,fi));
    sim_initial_tag_release_by_age(it)=
      ecatch/(sum(ecatch)+1.0e-15)
        *sim_tag_numbers_released(it);

    
    sim_tagN(it).initialize();
    sim_tagN(it,ir,yr1)=log(sim_initial_tag_release_by_age(it)+1.e-12);

    // here is where we "diffuse " the tagged fish to move them between
    // the regions
    if (num_regions>1) do_the_diffusion(yr1,sv,sim_tagN(it));
  }
}

  
  
void dvar_fish_stock_history::get_sim_initial_tag_recruitment_period(void)
{
  int irc=initial_recruitment_count;
  for (int it=1;it<=sim_num_tag_releases;it++)
  {  
    for (int ir=1;ir<=num_regions;ir++)
    {
      int ip=sim_initial_tag_period(it,ir);
      int i=0;
      for (i=1;i<=age_flags(93);i++)
      {
        if (month(ir,ip)< 1000)
           break;
        if (rec_times(i)<1000)
           break;
        if (month(ir,ip)< rec_times(i))
           break;
      }
      sim_initial_tag_recruitment_period(it,ir)=irc+age_flags(93)*(year(ir,ip)-2.)
           +i-1;
    }
  }
}
