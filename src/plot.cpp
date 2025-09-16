/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#define HOME_VERSION
#include "all.hpp"

#include "variable.hpp"
extern adstring full_input_parfile_path;
extern adstring full_output_parfile_path;
extern adstring full_datafile_path;
static void  xxx(dvar_vector& u){}
int debug_flag=0;

void show_version(void);
void print_version(ofstream& ofs);

void print_identifier_stuff(ofstream& ofs)
{
  ofs << "# Frq file = " << full_datafile_path << "\n";
  ofs << "# Input par file = " << full_input_parfile_path << "\n";
  ofs << "# Output par file = " << full_output_parfile_path << "\n";
}
void print_length_ssmult_stuff
  (dvar_len_fish_stock_history& fsh,const ofstream & ofs);

void get_length_sel_at_age
  (dvar_len_fish_stock_history& fsh)
{
  for (int i=1;i<=fsh.num_fisheries;i++)
  {
    int checkis=0;
    int checkib=0;
    for (int j=1;j<=fsh.num_fish_times(i);j++)
    {
      int rr=fsh.realization_region(i,j);
      int rp=fsh.realization_period(i,j);
      int ri=fsh.realization_incident(i,j);
      int is=fsh.sseason(rr,rp,ri);
      int ib=fsh.bblock(rr,rp,ri);
      if (is > checkis || ib > checkib)
      {
        fsh.bstempsel_afl(i,is,ib)=fsh.incident_sel(rr,rp,ri);
	if (is > checkis) checkis++;
	if (ib > checkib) checkib++;
      }
    }
  }
//  cout << "Finished bstempsel_afl" << endl;
}


void ests_write1(ofstream& of1,dvar_len_fish_stock_history& fsh,void * pq_flag)
{
  if (!pq_flag && (fsh.parest_flags(139)>8 || fsh.parest_flags(141)>8))
  {
    fsh.print_ssmult_stuff();
  }
  if(fsh.parest_flags(186)){
    ofstream ofs("fishmort");
    ofs << fsh.fish_mort << endl;

    ofstream ofs2("fishmort2");
    ofs2 << "Region" << "   " << "Period" << "   " << "Incident" << "   "
         << "fish_mort(reg,pd,fi): ages" << endl;
    for (int ir=1;ir<=fsh.num_regions;ir++)
    {
      for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)
      {
        for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
        {
          ofs2 << ir << "  "  << ip << "  " << fi << "  "
               << fsh.fish_mort(ir,ip,fi) << endl;
        }
      }
    }
  }
  int i;
  int j;
  int ir;
  int iy;

  if(fsh.parest_flags(187)) fsh.test_tag_report(pq_flag);

  of1 << "# MULTIFAN-CL Viewer" << endl;
  //  if(!fsh.pmsd)   //NMD 4Apr2012
  if(fsh.parest_flags(200) < 1043)   //NMD 20Feb2014
  {
    of1 << "# 2.0" << endl;
  } 
  else if(fsh.parest_flags(200) < 1052)
  {
    of1 << "# 3.0" << endl;
  }   //NMD 4Apr2012
  else if(fsh.parest_flags(200) >= 1052)
  {
    of1 << "# 4.0" << endl;
  }   //NMD 20Dec2016

  print_identifier_stuff(of1);
  print_version(of1); //NMD_14Dec2016
  of1 << "# Number of time periods" << endl << "  " << fsh.nyears << endl;
  of1 << "# Year 1" << endl << "  " << fsh.year1 << endl;
  of1 << "# Number of regions" << endl << "  " << fsh.num_regions << endl;
  if(fsh.parest_flags(200) > 1042)   //NMD 20Feb2014
  {
    if(!fsh.pmsd)   //NMD 4Apr2012
    {
       of1 << "# Number of species " << endl << "  " << " 1 " << endl;
    } else {
       of1 << "# Number of species" << endl << "  " << fsh.pmsd->num_species << endl;
    }   //NMD 4Apr2012
  }     //NMD 20Feb2014

  if(fsh.parest_flags(200) > 1051)   //NMD 11Jan2017
  {

    if(!fsh.pmsd)
    {
      int species_flag_pointer;
      species_flag_pointer=1;
      of1 << "# Multi species pointer " << endl << "  " 
          << species_flag_pointer  << endl;
      species_flag_pointer=-1;
      of1 << "# Species sex pointer " << endl << "  " 
          << species_flag_pointer  << endl;
    } else {
      ivector species_flag_pointer(1,fsh.pmsd->num_species);
      species_flag_pointer=column(fsh.pmsd->species_flags,1);
      of1 << "# Multi species pointer " << endl << "  " 
          << species_flag_pointer  << endl;
      species_flag_pointer=column(fsh.pmsd->species_flags,2);
      of1 << "# Species sex pointer " << endl << "  " 
          << species_flag_pointer  << endl;
    }
  }     //NMD 11Jan2017



  if(!fsh.pmsd)   //NMD 4Apr2012
    {
       of1 << "# Number of age classes" << endl << "  " << fsh.nage << endl;
    } else {
       of1 << "# Number of age classes" << endl << "  " << fsh.nage << 
	 "  " << fsh.pmsd->nage << endl;
    }   //NMD 4Apr2012
  //  if(!fsh.pmsd)   //NMD 4Apr2012
  if(fsh.parest_flags(200) > 1042)   //NMD 20Feb2014
  {
    if(!fsh.pmsd)   //NMD20Feb2014
    {
      ivector rsp(1,fsh.num_regions);
      rsp = 1;
      of1 << "# Regions species pointer " << endl << "  " << rsp  << endl;
    } else {
      of1 << "# Regions species pointer " << endl << "  " << 
             fsh.pmsd->region_species_pointer << endl;
    }   //NMD 4Apr2012
  }     //NMD 20Feb2014
  if(fsh.parest_flags(200) > 1051)   //NMD 20Dec2016
  {
    ivector fishery_species_pointer(1,fsh.num_fisheries);
    if(!fsh.pmsd)
    {
      fishery_species_pointer=1;
      of1 << "# Fishery species pointer " << endl << "  " 
          << fishery_species_pointer  << endl;
    } else {
      int mmin=fsh.fishery_regions.indexmin();
      int mmax=fsh.fishery_regions.indexmax();
      for(int ir=mmin;ir<=mmax;ir++)
      {
        int ireg=fsh.fishery_regions(ir);
        fishery_species_pointer(ir)=fsh.pmsd->region_species_pointer(ireg);
      }
      of1 << "# Fishery species pointer " << endl << "  " 
          << fishery_species_pointer << endl;
    }
  }     //NMD 20Dec2016


  of1 << "# Number of length bins" << endl << "  " << fsh.fmid.size() << endl;
  of1 << "# Number of recruitments per year" << endl << "  " << 
          fsh.age_flags(57)<< endl;
  of1 << "# Number of fisheries" << endl << "  " << fsh.num_fisheries << endl;
  of1 << "# Fishery selectivity seasons" << endl;   //NMD_20May2015
  of1 << column(fsh.fish_flags,74) << endl;
  of1 << "# Fishery selectivity time-blocks" << endl;   //NMD_20May2015
  of1 << column(fsh.fish_flags,71) + 1 << endl;
  of1 << "# Number of realizations per fishery" << endl
     << setw(4) << fsh.num_fish_times << endl;
  of1 << "# Region for each fishery" << endl
     << fsh.fishery_regions << endl;
  of1 << "# Time of each realization by fishery (down)" << endl;
  for (int i=1;i<=fsh.num_fisheries;i++)
  {
    for (int j=1;j<=fsh.num_fish_times(i);j++)
    {
      of1 << setfixed() << setprecision(3) << setw(10) 
          << fsh.fishing_incident_times(i,j);
    }
    of1 << endl;
  }

  //dvector ml(1,fsh.nage);
  int ns=1;
  if (fsh.pmsd)
   ns=fsh.pmsd->num_species;
  for (int is=1;is<=ns;is++)
  {
    int ng=fsh.get_nage_species(is);   //NMD_10Oct2018
    dvector vbc(1,4);
    dvector vc(1,2);
    if (is==1)
    {
      vbc=value(fsh.vb_coff);
      vc=value(fsh.var_coff);
    }
    else
    {
      vbc=value(fsh.pmsd->vb_coff(is));
      vc=value(fsh.pmsd->var_coff(is));
    }
    MY_DOUBLE_TYPE rho=exp(-vbc(3));
    MY_DOUBLE_TYPE dd=1-pow(rho,(ng-1));
    dvector ml(1,ng);
    dvector mw(1,ng);
    dvector sdl(1,ng);
    
    if (is==1)
    {
      ml=value(fsh.mean_length_yr(1,1));
      mw=value(fsh.mean_weight_yr(1,1));
      sdl=value(sqrt(fsh.global_vars));
    }
    else
    {
      int rb=fsh.pmsd->region_bounds(is,1);
      ml=value(fsh.mean_length_yr(rb,1));
      mw=value(fsh.mean_weight_yr(rb,1));
      sdl=value(sqrt(fsh.pmsd->global_vars(is)));
    }
    /*
    if (fsh.parest_flags(226)==0)        //NMD_10Oct2018
    {
      for (j=1;j<=ng;j++)
      {
        if (fsh.parest_flags(157))
          ml(j)=vbc(1)+(vbc(2)-vbc(1))*
                ((1-pow(rho,(j-1)))/dd);
        sdl(j)=vc(1)*exp(vc(2)*(-1.+2*(1-pow(rho,(j-1)))/dd));
      }
    }
    else
    {        //NMD_10Oct2018
      dvector scaled_ml=setm11(ml);   //scale the ml between -1 and 1
      sdl=vc(1)*exp(vc(2)*scaled_ml);

    }        //NMD_10Oct2018
    */
    of1 << "# Mean lengths at age" << endl;
    of1 << setfixed()<< setprecision(4) << setw(7) << ml << endl;
    of1 << "# SD of length at age" << endl;
    of1 << setfixed()<< setprecision(4) << setw(7) << sdl << endl;
    of1 << "# Mean weights at age" << endl;
    of1 << setfixed()<< setprecision(4) << setw(7) << mw << endl;
  }

  of1 << "# Natural mortality at age" << endl;
  for (int is=1;is<=ns;is++)
  {
//    dvector nm(1,fsh.nage);
    if (is==1)
    {
      dvector nm(1,fsh.nage);
      nm=exp(value(fsh.nat_mort(1)));
      of1 << setfixed()<< setprecision(4) << setw(7) 
        << nm << endl;
    }
    else
    {
      dvector nm(1,fsh.pmsd->nage(is));
      nm=exp(value(fsh.pmsd->nat_mort(is,1)));
      of1 << setfixed()<< setprecision(4) << setw(7) 
        << nm << endl;
    }    //NMD_2April2020    
//    if (is==1)
//    {
//      nm=exp(value(fsh.nat_mort(1)));
//    }
//    else
//    {
//      nm=exp(value(fsh.pmsd->nat_mort(is,1)));
//    }
//    of1 << setfixed()<< setprecision(4) << setw(7) 
//       << nm << endl;
  }
  of1 << "# Selectivity by age class (across) and fishery (down)" << endl;
  ivector ff16=column(fsh.fish_flags,16);
  ivector ff57=column(fsh.fish_flags,57);
  ivector ff74=column(fsh.fish_flags,74);
  ivector ff71=column(fsh.fish_flags,71);
  ivector ff75=column(fsh.fish_flags,75);
  ivector nd=column(fsh.fish_flags,3);
  d4_array sel(1,fsh.num_fisheries,1,ff74,1,ff71+1,1,fsh.nage);
  d3_array ms;

  if (fsh.ff263flag) get_length_sel_at_age(fsh);
  
  for (i=1;i<=fsh.num_fisheries;i++)
  {
    for (int is=1;is<=ff74(i);is++)
    {
      for (int ib=1;ib<=ff71(i)+1;ib++)
      {
        if (!fsh.ff263flag)
        {
          if (!ff75(i))
          {
            sel(i,is,ib)(1,nd(i))=exp(value(fsh.bstempsel(i,is,ib)));
          }
          else
          {
            sel(i,is,ib)(1,ff75(i))=0.0;
            sel(i,is,ib)(ff75(i)+1,nd(i))=exp(value(fsh.bstempsel(i,is,ib)));
          }
          if (nd(i)<fsh.nage)
          {
            if (ff57(i)>0)
            {
              sel(i,is,ib)(nd(i)+1,fsh.nage)=sel(i,is,ib,nd(i));
            }
            else
            {
              if (ff16(i)==1)
              {
                sel(i,is,ib)(nd(i)+1,fsh.nage)=sel(i,is,ib,nd(i));
              }
              else if (ff16(i)==2)
              {
                sel(i,is,ib)(nd(i),fsh.nage)=0.0;
              }
              else
              {
                cerr << "ERROR: incorrect fish_flags(16) setting, fishery: "
                    << i << endl;
                ad_exit(1);
              }
            }
          }
        }
        else
        {
          sel(i,is,ib)=exp(value(fsh.bstempsel_afl(i,is,ib)));
        }
      }
    }
  }
  if (fsh.pmsd) 
  {
    ivector nd=column(fsh.fish_flags,3);
    for (int i=1;i<=fsh.num_fisheries;i++)
    {
      if (nd(i)==0) nd(i)=fsh.nage;
    }

    int nrf=fsh.pmsd->num_real_fisheries;
    ms.allocate(1,fsh.num_fisheries,1,ff74,1,ff71+1);
    ms.initialize();
    int nage=fsh.nage;
    for (int i=1;i<=nrf;i++)
    {
      for (int is=1;is<=ff74(i);is++)
      {
        for (int ib=1;ib<=ff71(i)+1;ib++)
        {
          //ms(i,is,ib)=value(exp(max(fsh.bstempsel(i,is,ib,nd(i)),
          MY_DOUBLE_TYPE d1=max(sel(i,is,ib));
          MY_DOUBLE_TYPE d2=max(sel(i+nrf,is,ib));
          ms(i,is,ib)=max(d1,d2);
          ms(i+nrf,is,ib)=ms(i,is,ib);
//          cout << endl << ms(i,is,ib);      //NMD_24Sep2018
//          cout << endl << sel(i,is,ib) << endl;
//          cout << sel(i+nrf,is,ib) << endl;
//          cout << endl << sel(i,is,ib)/ms(i,is,ib) << endl;
//          cout << sel(i+nrf,is,ib)/ms(i,is,ib) << endl;      //NMD_24Sep2018
        }
      }
    }
  }
  ofstream of2("selectivity-multi-sex");
  int nrf=0;
  if (fsh.pmsd)
  {
    nrf=fsh.pmsd->num_real_fisheries;
  }
  for (i=1;i<=fsh.num_fisheries;i++)
  {
    for (int is=1;is<=ff74(i);is++)
    {
      for (int ib=1;ib<=ff71(i)+1;ib++)
      {
        if (fsh.pmsd  && fsh.pmsd->num_species==2 && (fsh.fish_flags(i,57)==3 || fsh.fish_flags(i,57)==1 )
          && fsh.age_flags(193) )
        {
          sel(i,is,ib)/=ms(i,is,ib);
        }
        else
        {
	  MY_DOUBLE_TYPE maxs=max(sel(i,is,ib));
          sel(i,is,ib)/=maxs;
        }
        of1 << setfixed()<< setprecision(4) << setw(7) 
            << sel(i,is,ib) << endl;
      }
    }
  }
  for (i=1;i<=nrf;i++)
  {
    for (int is=1;is<=ff74(i);is++)
    {
      for (int ib=1;ib<=ff71(i)+1;ib++)
      {
        if (fsh.pmsd  && fsh.pmsd->num_species==2 && (fsh.fish_flags(i,57)==3 || fsh.fish_flags(i,57)==1 )
          && fsh.age_flags(193) )
        {
          if ( i<=nrf )
          {
            of2 << "# fishery " << i << endl;
            of2 << setfixed()<< setprecision(4) << setw(7) 
                << sel(i,is,ib) << endl;
            of2 << setfixed()<< setprecision(4) << setw(7) 
              << sel(i+nrf,is,ib) << endl;
          }
        }
      }
    }
  }
  //if (fsh.ff263flag)
  {
    of1 << "# length bin mid-points" << endl;
	of1 << fsh.fmid << endl;
    of1 << "# length-based selectivity by fishery" << endl;
    dvector minusone(1,fsh.nlint);
    minusone=-1;
    for (int fi=1;fi<=fsh.num_fisheries;fi++)
    {
      if (fsh.fish_flags(fi,26)==3)
      {
//        MY_DOUBLE_TYPE fmax=max(value(fsh.splinesel(fi)));
//        of1 << value(fsh.splinesel(fi))/fmax << endl;
        for (int is=1;is<=ff74(fi);is++)
        {
          for (int ib=1;ib<=ff71(fi)+1;ib++)
          {
	    MY_DOUBLE_TYPE maxs=max(value(exp(fsh.bstempsel(fi,is,ib))));
	    of1 << value(exp(fsh.bstempsel(fi,is,ib)))/maxs << endl;
          }
        }
      }
      else
      {
        of1 << minusone << endl;
      }
    }
  }

  if (!fsh.age_flags(92))
  {
    of1 << "# Catchability by realization (across) by fishery (down)" << endl;
  }
  if (fsh.age_flags(92) && allocated(fsh.grouped_catchability_coffs))
  {
    of1 << "# Kalman filter smoothed catchability by realization (across) "
        << "by fishery (down)" << endl;
  }
  MY_DOUBLE_TYPE fzero=0.0;
  ivector grouping=column(fsh.fish_flags,29);
  for (i=1;i<=fsh.num_fisheries;i++)
  {
    for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
    {
      int rr=fsh.realization_region(i,nt);
      int rp=fsh.realization_period(i,nt);
      int ri=fsh.realization_incident(i,nt);
      if (fsh.age_flags(92)==0)
      {
        of1 << setscientific() << setprecision(4) << " " 
            << exp(value(fsh.catchability(rr,rp,ri)));
      }
      else
      {
        //if (fsh.age_flags(92))
        if (nt<=fsh.num_real_fish_times(i))
        {
          if (allocated(fsh.grouped_catchability_coffs))
          {
            of1 << setscientific() << setprecision(4) << " " 
                << exp(fsh.grouped_catchability_coffs
                   (grouping(i),fsh.gfish_ptr(i,nt)));
          }
        }
        else
        {
          of1 << setscientific() << setprecision(4) << " " 
              << fzero; 
        }
      }
    }
    of1 << endl;
  }
  if (!fsh.age_flags(92))
    of1 << "# Catchability+effort dev. by realization (across) "
        << "by fishery (down)" << endl;
  if (fsh.age_flags(92))
    of1 << "# Implicit catchability by realization (across) "
        << "by fishery (down)" << endl;
  for (i=1;i<=fsh.num_fisheries;i++)
  {
    for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
    {
      int rr=fsh.realization_region(i,nt);
      int rp=fsh.realization_period(i,nt);
      int ri=fsh.realization_incident(i,nt);
      if (!fsh.age_flags(92))
        of1 << setscientific() << setprecision(4) << " " 
            << exp(value(fsh.catchability(rr,rp,ri)))*
               exp(value(fsh.effort_devs(rr,rp,ri)));
      if (fsh.age_flags(92))
        of1 << setscientific() << setprecision(4) << " " 
            << exp(value(fsh.implicit_catchability(i,nt)));
    }
    of1 << endl;
  }
  // ***************************************************************
  // ***************************************************************
  //  new code for projections
  if (fsh.projection_sim_index>0 &&!pq_flag && fsh.parest_flags(191))  //NMD_8May2018 
  {
    ofstream& ofxx=*(fsh.proj_output_files[0]);
    ofxx << "#simulation number " << fsh.projection_sim_index << endl;
    if (!fsh.age_flags(92))
      ofxx << "# Catchability by realization (across) by fishery (down)" << endl;
    if (fsh.age_flags(92))
      ofxx << "# Kalman filter smoothed catchability by realization (across) "
          << "by fishery (down)" << endl;
    MY_DOUBLE_TYPE fzero=0.0;
    ivector grouping=column(fsh.fish_flags,29);
    for (i=1;i<=fsh.num_fisheries;i++)
    {
      for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
      {
        int rr=fsh.realization_region(i,nt);
        int rp=fsh.realization_period(i,nt);
        int ri=fsh.realization_incident(i,nt);
        if (fsh.age_flags(92)==0)
        {
          ofxx << setscientific() << setprecision(4) << " " 
              << exp(value(fsh.catchability(rr,rp,ri)));
        }
        else
        {
          if (allocated(fsh.grouped_catchability_coffs))
          {
          //if (fsh.age_flags(92))
          if (nt<=fsh.num_real_fish_times(i))
          {
            ofxx << setscientific() << setprecision(4) << " " 
                << exp(fsh.grouped_catchability_coffs
                      (grouping(i),fsh.gfish_ptr(i,nt)));
          }
          else
          {
            ofxx << setscientific() << setprecision(4) << " " 
                << fzero; 
          }
          }
        }
      }
      ofxx << endl;
    }
    if (!fsh.age_flags(92))
      ofxx << "# Catchability+effort dev. by realization (across) "
          << "by fishery (down)" << endl;
    if (fsh.age_flags(92))
      ofxx << "# Implicit catchability by realization (across) "
          << "by fishery (down)" << endl;
    for (i=1;i<=fsh.num_fisheries;i++)
    {
      for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
      {
        int rr=fsh.realization_region(i,nt);
        int rp=fsh.realization_period(i,nt);
        int ri=fsh.realization_incident(i,nt);
        if (!fsh.age_flags(92))
          ofxx << setscientific() << setprecision(4) << " " 
              << exp(value(fsh.catchability(rr,rp,ri)))*
                 exp(value(fsh.effort_devs(rr,rp,ri)));
        if (fsh.age_flags(92))
          ofxx << setscientific() << setprecision(4) << " " 
              << exp(value(fsh.implicit_catchability(i,nt)));
      }
      ofxx << endl;
    }
  }
  // *********************************************************
  // *********************************************************
  of1 << "# Fishing mortality by age class (across) and year (down)" << endl;

  fsh.get_fishing_mortality_by_age_by_year_by_region();
  if (fsh.pmsd==0)
  {
    if(fsh.parest_flags(200) > 1042) of1 << "# Species 1" << endl; //NMD_20Feb2014
    fsh.get_fishing_mortality_by_age_by_year();
    of1 << setprecision(6) << fsh.F_by_age_by_year << endl;
  }
  else
  {
    for (int is=1;is<=fsh.pmsd->num_species;is++)
    {
      of1 << "# Species " << is << endl;
      int rmin=fsh.pmsd->region_bounds(is,1);
      int rmax=fsh.pmsd->region_bounds(is,2);
      dmatrix tmp(1,fsh.nyears,1,fsh.nage);
      tmp.initialize();
      for (int ir=rmin;ir<=rmax;ir++)
      {
        tmp+=value(fsh.F_by_age_by_year_by_region(ir));
      }
      of1 << setprecision(6) << tmp << endl;
    }
  }

  of1 << "# Fishing mortality by age class (across), year (down) and region (block)" 
      << endl;
  for (ir=1;ir<=fsh.num_regions;ir++)
  {
    of1 << "# Region " << ir << endl;
    of1 << setprecision(6) << fsh.F_by_age_by_year_by_region(ir) << endl;
  }

 //********************************************************** 
 //********************************************************** 
 //********************************************************** 
  of1 << "# Population Number by age (across), year (down) and region" << endl;
  ofstream& of3=*(fsh.proj_output_files[3]);
  ofstream& of4=*(fsh.proj_output_files[4]);   //NMD 8Nov2011
  ofstream& of5=*(fsh.proj_output_files[5]);   //NMD 10Nov2011
  if (pq_flag==0)
  {
    if (fsh.projection_sim_index>0) 
    {
      of3 << "# Simulation number " << fsh.projection_sim_index << endl;
      of3 << "# Population Number by age (across), year (down) and region" 
          << endl;
      of4 << "# Simulation number " << fsh.projection_sim_index << endl;   //NMD 8Nov2011
      of4 << "# Fishing mortality by age class (across), year (down)"         //NMD 8Nov2011
          << endl;
      of4 << setprecision(6) << fsh.F_by_age_by_year << endl;                //NMD 8Nov2011
    }
//    int nage=fsh.nage;   //NMD_2Apr2020
    for (ir=1;ir<=fsh.num_regions;ir++)
    {
//      dvar_vector avgN(1,nage);
//      dvar_vector avgz(1,nage-1);
//      dvar_vector avgN_proj(1,nage);
//      dvar_vector avgz_proj(1,fsh.nage-1);
      dvar_vector avgN(1,fsh.nage_by_region(ir));
      dvar_vector avgz(1,fsh.nage_by_region(ir)-1);
      dvar_vector avgN_proj(1,fsh.nage_by_region(ir));
      dvar_vector avgz_proj(1,fsh.nage_by_region(ir)-1);
      avgN.initialize();
      avgN_proj.initialize();
      avgz.initialize();
      avgz_proj.initialize();
      of1 << "# Region " << ir << endl;
      if (fsh.projection_sim_index>0) 
        of3 << "# Region " << ir << endl;
      int iy;
      for (iy=1;iy<=fsh.last_real_year;iy++)
      {
        of1 <<  setscientific() << setprecision(5) 
            << exp(value(fsh.N(ir,iy))) << endl;
        if (fsh.projection_sim_index>0) 
          of3 <<  setscientific() << setprecision(5) 
              << exp(value(fsh.N(ir,iy))) << endl;
      }
      int icount=0;
      for (iy=1;iy<fsh.last_real_year;iy++)
      {
        avgN+=exp(fsh.N(ir,iy));
//        avgz+=fsh.N(ir,iy)(1,nage-1)-fsh.N(ir,iy+1)(2,nage).shift(1);
        avgz+=fsh.N(ir,iy)(1,fsh.nage_by_region(ir)-1)-
          fsh.N(ir,iy+1)(2,fsh.nage_by_region(ir)).shift(1);    //NMD_2Apr2020
        icount++;
      }
      avgN/=icount;
      avgz/=icount;
      if (fsh.last_real_year<fsh.nyears)
      {
        of1 << "#   Projected years " << endl;
        if (fsh.projection_sim_index>0) 
          of3 << "#   Projected years " << endl;
        for (iy=fsh.last_real_year+1;iy<=fsh.nyears;iy++)
        {
          of1 <<  setscientific() << setprecision(5) 
              << exp(value(fsh.N(ir,iy))) << endl;
          if (fsh.projection_sim_index>0) 
            of3 <<  setscientific() << setprecision(5) 
                << exp(value(fsh.N(ir,iy))) << endl;
        }
        icount=0;
        for (iy=fsh.last_real_year+1;iy<fsh.nyears;iy++)
        {
          avgN_proj+=exp(fsh.N(ir,iy));
//          avgz_proj+=fsh.N(ir,iy)(1,nage-1)-fsh.N(ir,iy+1)(2,nage).shift(1);
          avgz_proj+=fsh.N(ir,iy)(1,fsh.nage_by_region(ir)-1)-
            fsh.N(ir,iy+1)(2,fsh.nage_by_region(ir)).shift(1);    //NMD_2Apr2020
          icount++;
        }
        avgN_proj/=icount;
        avgz_proj/=icount;
      }
      ofstream ofs("avgz");
      ofs << setw(9) << setprecision(4) << avgN  << endl;
      ofs << setw(9) << setprecision(4) << avgN_proj << endl;
      ofs << setw(9) << setprecision(4) << avgz << endl;
      ofs << setw(9) << setprecision(4) << avgz_proj << endl;
    }
  }
  else
  {
    if (fsh.projection_sim_index>0)         //NMD 1Dec 2011
	{
      of5 << "# Simulation number " << fsh.projection_sim_index << endl;       //NMD 10Nov2011
      of5 << "# Population Number by age (across), year (down) and region"     //NMD 10Nov2011
        << endl;                                                             //NMD 10Nov2011
	}
    for (ir=1;ir<=fsh.num_regions;ir++)
    {
      of1 << "# Region " << ir << endl;
      if (fsh.projection_sim_index>0)         //NMD 10Nov 2011
        of5 << "# Region " << ir << endl;
      for (int iy=1;iy<=fsh.last_real_year;iy++)
      {
        of1 <<  setscientific() << setprecision(3) 
            << exp(value(fsh.N_q0(ir,iy))) << endl;
        if (fsh.projection_sim_index>0)      //NMD 10Nov 2011
          of5 <<  setscientific() << setprecision(3) 
              << exp(value(fsh.N_q0(ir,iy))) << endl;
      }
      if (fsh.last_real_year<fsh.nyears)
      {
        of1 << "#   Projected years " << endl;
        if (fsh.projection_sim_index>0)      //NMD 10Nov 2011
          of5 << "#   Projected years " << endl;
        for (int iy=fsh.last_real_year+1;iy<=fsh.nyears;iy++)
        {
          of1 <<  setscientific() << setprecision(3) 
              << exp(value(fsh.N_q0(ir,iy))) << endl;
          if (fsh.projection_sim_index>0)   //NMD 10Nov 2011
            of5 <<  setscientific() << setprecision(3) 
              << exp(value(fsh.N_q0(ir,iy))) << endl;
        }
      }
    }
  }

  //  new code for projections
  if (fsh.projection_sim_index>0 && fsh.parest_flags(191))  //NMD_8May2018  
  {
    if (pq_flag==0)
    {
      ofstream& of1=*(fsh.proj_output_files[0]);    //NMD 11Nov 2011
      of1 << "# Population Number by age (across), year (down) and region" << endl;       //NMD 11Nov 2011
      int nage=fsh.nage;
      for (ir=1;ir<=fsh.num_regions;ir++)
      {
        dvar_vector avgN(1,nage);
        dvar_vector avgz(1,nage-1);
        dvar_vector avgN_proj(1,nage);
        dvar_vector avgz_proj(1,fsh.nage-1);
        avgN.initialize();
        avgN_proj.initialize();
        avgz.initialize();
        avgz_proj.initialize();
        of1 << "# Region " << ir << endl;
        int iy;
        for (iy=1;iy<=fsh.last_real_year;iy++)
        {
          of1 <<  setscientific() << setprecision(3) 
              << exp(value(fsh.N(ir,iy))) << endl;
        }
        int icount=0;
        for (iy=1;iy<fsh.last_real_year;iy++)
        {
          avgN+=exp(fsh.N(ir,iy));
          avgz+=fsh.N(ir,iy)(1,nage-1)-fsh.N(ir,iy+1)(2,nage).shift(1);
          icount++;
        }
        avgN/=icount;
        avgz/=icount;
        if (fsh.last_real_year<fsh.nyears)
        {
          of1 << "#   Projected years " << endl;
          for (iy=fsh.last_real_year+1;iy<=fsh.nyears;iy++)
          {
            of1 <<  setscientific() << setprecision(3) 
                << exp(value(fsh.N(ir,iy))) << endl;
          }
          icount=0;
          for (iy=fsh.last_real_year+1;iy<fsh.nyears;iy++)
          {
            avgN_proj+=exp(fsh.N(ir,iy));
            avgz_proj+=fsh.N(ir,iy)(1,nage-1)-fsh.N(ir,iy+1)(2,nage).shift(1);
            icount++;
          }
          avgN_proj/=icount;
          avgz_proj/=icount;
        }
        ofstream ofs("avgz");
        ofs << setw(9) << setprecision(4) << avgN  << endl;
        ofs << setw(9) << setprecision(4) << avgN_proj << endl;
        ofs << setw(9) << setprecision(4) << avgz << endl;
        ofs << setw(9) << setprecision(4) << avgz_proj << endl;
      }
    }
    else
    {
     /*
      for (ir=1;ir<=fsh.num_regions;ir++)
      {
        of1 << "# Region " << ir << endl;
        for (int iy=1;iy<=fsh.last_real_year;iy++)
        {
          of1 <<  setscientific() << setprecision(3) 
              << exp(value(fsh.N_q0(ir,iy))) << endl;
        }
        if (fsh.last_real_year<fsh.nyears)
        {
          of1 << "#   Projected years " << endl;
          for (int iy=fsh.last_real_year+1;iy<=fsh.nyears;iy++)
          {
            of1 <<  setscientific() << setprecision(3) 
                << exp(value(fsh.N_q0(ir,iy))) << endl;
          }
        }
      }
    */
    }
  }
  imatrix m74(1,fsh.nyears,ff74);

  of1 << "# Exploitable population biomass by fishery (down) and by year-season  (across)"
      << endl;
  d3_array exp_pop_bio(1,fsh.nyears,1,fsh.num_fisheries,1,m74);
  exp_pop_bio.initialize();
  if (pq_flag==0)
  {
    for (iy=1;iy<=fsh.nyears;iy++)
    {
      for (i=1;i<=fsh.num_fisheries;i++)
      {
        for (int is=1;is<=ff74(i);is++)
        {
          int rr=fsh.realization_region(i,1);
          for (j=1;j<=fsh.nage;j++)
          {
            int ib=fsh.yearblock(i,iy);
            if (ib>0)
            {
              exp_pop_bio(iy,i,is)+=sel(i,is,ib,j)*
                exp(value(fsh.N(rr,iy,j)))*
                value(fsh.mean_weight_yr(rr,iy,j))/1000.;
            }
          }
        }
      }
    }
  }
  else
  {
    for (iy=1;iy<=fsh.nyears;iy++)
    {
      for (i=1;i<=fsh.num_fisheries;i++)
      {
        for (int is=1;is<=ff74(i);is++)
        {
          int rr=fsh.realization_region(i,1);
          int ib=fsh.yearblock(i,iy);
          for (j=1;j<=fsh.nage;j++)
          {
            if (ib>0)
            {
              exp_pop_bio(iy,i,is)+=sel(i,is,ib,j)*
                exp(value(fsh.N_q0(rr,iy,j)))*
                value(fsh.mean_weight_yr(rr,iy,j))/1000.;
            }
          }
        }
      }
    }
 }
 for (i=1;i<=fsh.num_fisheries;i++)
 {
   for (iy=1;iy<=fsh.nyears;iy++)
   {
     for (int is=1;is<=ff74(i);is++)
     {
       of1 <<  setscientific() << setprecision(3) << exp_pop_bio(iy,i,is) 
           << " ";
     }
   }
   of1 << endl; 
 }

  of1 << "# Exploitable population in same units as catch by fishery (down) and year-season (across)" << endl;
  d3_array exp_pop(1,fsh.nyears,1,fsh.num_fisheries,1,m74);
  exp_pop.initialize();
  for (iy=1;iy<=fsh.nyears;iy++)
  {
    for (i=1;i<=fsh.num_fisheries;i++)
    {
      for (int is=1;is<=ff74(i);is++)
      {
        int rr=fsh.realization_region(i,1);
        int ib=fsh.yearblock(i,iy);
        for (j=1;j<=fsh.nage;j++)
        {
          if (ib>0)
          {
            if (fsh.data_fish_flags(1,i)==0)
            {
              exp_pop(iy,i,is)+=sel(i,is,ib,j)*exp(value(fsh.N(rr,iy,j)));
            }  
            if (fsh.data_fish_flags(1,i)==1)
            {
              exp_pop(iy,i,is)+=sel(i,is,ib,j)*exp(value(fsh.N(rr,iy,j)))*
                 value(fsh.mean_weight_yr(rr,iy,j))/1000.;
            }
          }
        }
      }
    }
  }
  for (i=1;i<=fsh.num_fisheries;i++)
  {
    for (iy=1;iy<=fsh.nyears;iy++)
    {
      for (int is=1;is<=ff74(i);is++)
      {
        of1 <<  setscientific() << setprecision(3) << exp_pop(iy,i,is) 
            << " ";
      }
    }
    of1 << endl; 
  }
 /*
  for (iy=1;iy<=fsh.nyears;iy++)
  {
    for (int ifish=1;ifish<=fsh.num_fisheries;ifish++)
    {
      int rp=fsh.realization_period(i,nt);
      int ri=fsh.realization_incident(i,nt);
      int rr=fsh.realization_region(ifish,1);
      int ib=fsh.bblock(ifish,iy);
      if (fsh.data_fish_flags(1,ifish)==0)
      {
        for (j=1;j<=fsh.nage;j++)
        {
          exp_pop(iy,ifish)+=sel(ifish,1,ib,j)*exp(value(fsh.N(rr,iy,j)));
        }
      }
      if (fsh.data_fish_flags(1,ifish)==1)
      {
        for (j=1;j<=fsh.nage;j++)
        {
          exp_pop(iy,ifish)+=sel(ifish,j)*exp(value(fsh.N(rr,iy,j)))*
                           value(fsh.mean_weight_yr(rr,iy,j))/1000.;
        }
      }
    }
    of1 <<  setscientific() << setprecision(3) << exp_pop(iy) << endl;
  }
  */

  of1 << "# Absolute biomass by region (across) and year (down)" << endl;
  int mmin=fsh.N(1).rowmin();
  int mmax=fsh.N(1).rowmax();
  dmatrix recx(1,fsh.num_regions,mmin,mmax);
  dmatrix rbio(1,fsh.num_regions,mmin,mmax);
  dmatrix sbio(1,fsh.num_regions,mmin,mmax);
  dvector pmature = value(fsh.pmature);
  if (!sum(pmature))
  {
    for (int j=1;j<=fsh.nage;j++)
    {
      pmature=1.0;
    }
  } 
  if (fsh.af170q0ex==0)
  {
    //dvar_vector * pm=0;  //NMD 9May 2012
    for (ir=1;ir<=fsh.num_regions;ir++)
    {
      //pm=&(fsh.get_pmature_region(ir));   //NMD 9May 2012
      dvar_vector pmature = fsh.get_pmature_region(ir);   //NMD 9May 2012
	  if (!value(sum(pmature)))
      {
        for (int j=1;j<=fsh.nage;j++)
        {
          pmature=1.0;
        }
      }                                 //NMD 9May 2012
      for (iy=mmin;iy<=mmax;iy++)
      {
        recx(ir,iy)=exp(value(fsh.N(ir,iy,1)));
        rbio(ir,iy)=sum(elem_prod(exp(value(fsh.N(ir,iy))),
                    value(fsh.mean_weight_yr(ir,iy))))/1000.;
        sbio(ir,iy)=sum(elem_prod(value(pmature),elem_prod(exp(value(fsh.N(ir,iy))),
                    value(fsh.mean_weight_yr(ir,iy)))))/1000.;
      }
    }
  }
  else
  {
    //dvar_vector * pm=0;  //NMD 9May 2012
    for (ir=1;ir<=fsh.num_regions;ir++)
    {
      //pm=&(fsh.get_pmature_region(ir));   //NMD 9May 2012
      dvar_vector pmature=fsh.get_pmature_region(ir);   //NMD 9May 2012
      if (!value(sum(pmature)))
      {
        for (int j=1;j<=fsh.nage;j++)
        {
          pmature=1.0;
        }
      }                                 //NMD 9May 2012
      for (iy=mmin;iy<=mmax;iy++)
      {
        recx(ir,iy)=exp(value(fsh.N_q0(ir,iy,1)));
        rbio(ir,iy)=sum(elem_prod(exp(value(fsh.N_q0(ir,iy))),
                    value(fsh.mean_weight_yr(ir,iy))))/1000.;
        sbio(ir,iy)=sum(elem_prod(value(pmature),elem_prod(exp(value(fsh.N_q0(ir,iy))),
		value(fsh.mean_weight_yr(ir,iy)))))/1000.;
      }
    }
  }
  of1 << "# Recruitment" << endl;
  of1 << trans(recx) << endl;
  of1 << "# Total biomass" << endl;
  of1 << setprecision(4)  << trans(rbio) << endl; 
  of1 << "# Adult biomass" << endl;
  of1 << setprecision(4)  << trans(sbio) << endl; 

//  if (fsh.projection_sim_index>0 && !pq_flag)
  if (fsh.projection_sim_index>0)  //NMD_13Jun2018
  {
    if (fsh.af170q0ex==0)
    {
      if (fsh.projection_sim_index==1) // start simulations with new file
      {
        remove("projected_spawning_biomass");
      }
      ofstream ofs("projected_spawning_biomass",ios::app);
      ofs << "# Simulation number " << fsh.projection_sim_index << endl;
      ofs << "# Adult biomass by region (across), and year (down) "
            << endl;
      ofs << setprecision(4)  << trans(sbio) << endl;
    }
    else
    {
      if (fsh.projection_sim_index==1) // start simulations with new file
      {
        remove("projected_spawning_biomass_noeff");
      }
      ofstream ofs("projected_spawning_biomass_noeff",ios::app);
      ofs << "# Simulation number " << fsh.projection_sim_index << endl;
      ofs << "# Adult biomass by region (across), and year (down) "
            << endl;
      ofs << setprecision(4)  << trans(sbio) << endl;
    }
  }
  
  of1 << "# Relative biomass by region (across) and year (down)" << endl;
  rbio/=max(rbio)+1.e-10;
  of1 << setprecision(4)  << trans(rbio) << endl; 

  of1 << "# Observed catch by fishery (down) and time (across)" << endl;
  for (i=1;i<=fsh.num_fisheries;i++)
  {
    for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
    {
      int rr=fsh.realization_region(i,nt);
      int rp=fsh.realization_period(i,nt);
      int ri=fsh.realization_incident(i,nt);
      MY_DOUBLE_TYPE tmp=fsh.obs_tot_catch(rr,rp,ri);
#if !defined(NO_MY_DOUBLE_TYPE)
      if (tmp==-1.0L) tmp=0.0;
#else
      if (tmp==-1.0) tmp=0.0;
#endif
      of1 << setprecision(3) << setw(12) << tmp;
    }
    of1 << endl;
  }
  of1 << "# Predicted catch by fishery (down) and time (across)" << endl;
  ofstream * pofs2=0;
  if (!pq_flag && fsh.age_flags(92)==2)
  {
    pofs2= new ofstream("Missing_catch_report");
    (*pofs2) << "Missing_catch_report" << endl << endl;;
    (*pofs2) << " Fishery   Fish time  Predicted Catch  Region Period  Incident" << endl;
  }
  ofstream ofs("totalcatches");  //NMD_1jul2024
  for (i=1;i<=fsh.num_fisheries;i++)
  {
    int nt;
//    if (i==1)
//    {
//      ofstream ofs("totalcatches");
    ofs << "Fishery:  " << i << endl;
    ofs << " " << endl;
    for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
    {
      int rr=fsh.realization_region(i,nt);
      int rp=fsh.realization_period(i,nt);
      int ri=fsh.realization_incident(i,nt);
      ofs << fsh.tot_catch(rr,rp,ri) << " " << sum(fsh.exp_catch(rr,rp,ri)) 
          << exp(fsh.fish_mort_calcs(rr,rp,ri)) << " " 
          << exp(fsh.num_fish(rr,rp)) << endl;
    }
   //    }
    for (nt=1;nt<=fsh.num_fish_times(i);nt++)
    {
      if (i==24 && nt==1)
        cout << "HERE" << endl;
      //if ( nt==fsh.num_fish_times(i))
      //{
      //  //cout << "HERE" << endl;
      //}
      int rr=fsh.realization_region(i,nt);
      int rp=fsh.realization_period(i,nt);
      int ri=fsh.realization_incident(i,nt);
      if (!fsh.data_fish_flags(1,fsh.parent(rr,rp,ri))) 
      {
        if (pq_flag)
        {
          of1 << setprecision(3) << setw(12) 
              << sum(exp(fsh.catch_q0(rr,rp,ri)));
        }
        else
        {
          of1 << setprecision(3) << setw(12) << value(fsh.tot_catch(rr,rp,ri));
          if (!pq_flag && fsh.age_flags(92)==2 &&
            fsh.missing_catch_for_incident_flag(rr,rp,ri)) //NMD_9jun2021
          {
            (*pofs2) << setw(5) << i << setw(10)  << nt  << setw(16)
                     << value(fsh.tot_catch(rr,rp,ri)) << setw(11) << rr
                     << setw(8) 
                     << rp << setw(8) 
                     << ri << endl;
          }
        }
      }
      else
      {
        dvariable totcatch=0.0;
        dvariable sv27=fsh.get_sv_region(rr,27);
        dvariable sv28=fsh.get_sv_region(rr,28);
        dvector tmp;
        if (pq_flag)
        {
          tmp=exp(value(fsh.catch_q0(rr,rp,ri)));
        }
        else
        {
          tmp=value(fsh.exp_catch(rr,rp,ri));
        }
        if (value(sv28)==3.0)
        {
          //totcatch=fsh.len_wt_coff/1000.* (fsh.exp_catch(rr,rp,ri)*
          totcatch=fsh.len_wt_coff/1000.* (tmp *
            (pow(fsh.mean_length(rr,rp,ri),3)+
            3.0*elem_prod(fsh.mean_length(rr,rp,ri),fsh.vars(rr,rp,ri))));
        }
        else
        {
          for (int j=1;j<=fsh.nage;j++)
          {  
            //totcatch+=fsh.exp_catch(rr,rp,ri,j)*
            totcatch+=tmp(j)*
              normal_length_to_weight(0.5,-3.5,
              fsh.mean_length(rr,rp,ri,j),sqrt(fsh.vars(rr,rp,ri,j)),
              value(sv27),value(sv28));   
          }
          totcatch/=1000.;
          if (i==18 && nt==fsh.num_real_fish_times(i)+1)
          {
            ofstream ofs("compare");
            ofs << fsh.mean_length(rr,rp,ri) << endl;
            ofs << sqrt(fsh.vars(rr,rp,ri)) << endl;
            ofs << fsh.exp_catch(rr,rp,ri) << endl;
            ofs << totcatch << endl;
          }
        }
        of1 << setprecision(3) << setw(12) << value(totcatch);
        if (!pq_flag && fsh.age_flags(92)==2 &&
          fsh.missing_catch_for_incident_flag(rr,rp,ri)) //NMD_1jun2020
        {
          (*pofs2) << setw(5) << i << setw(10)  << nt  << setw(16)
                   << value(totcatch) << setw(11) << rr << setw(8) 
                   << rp << setw(8) 
                   << ri << endl;
        }
      }
    }
    of1 << endl;
  }
  if (pofs2)
  {
    delete pofs2;
    pofs2=0;
  }

  // ***************************************************************
  // ***************************************************************
  //  new code for projections
  if (fsh.projection_sim_index>0 && fsh.parest_flags(191))  //NMD_8May2018 
  {
    ofstream& of1=*(fsh.proj_output_files[1]);
    ofstream& of2=*(fsh.proj_output_files[2]);
    of1 << "# Simulation " << fsh.projection_sim_index << endl;
    of1 << "# Observed catch by fishery (down) and time (across)" << endl;
    of2 << "# Simulation " << fsh.projection_sim_index << endl;
    for (i=1;i<=fsh.num_fisheries;i++)
    {
      for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
      {
        int rr=fsh.realization_region(i,nt);
        int rp=fsh.realization_period(i,nt);
        int ri=fsh.realization_incident(i,nt);
        MY_DOUBLE_TYPE tmp=fsh.obs_tot_catch(rr,rp,ri);
        //if (fsh.realization_period(i,nt)>fsh.num_real_fish_periods(ir))
        {
#if !defined(NO_MY_DOUBLE_TYPE)
          if (tmp==-1.0L) tmp=0.0;
#else
          if (tmp==-1.0) tmp=0.0;
#endif
          of1 << setprecision(3) << setw(12) << tmp;
        }  
      }
      of1 << endl;
    }
    of1 << "# Predicted catch by fishery (down) and time (across)" << endl;
    of2 << "# Predicted catch at age by fishery " << endl;
    for (i=1;i<=fsh.num_fisheries;i++)
    {
      of2 << "# Fishery " << i << endl;
      int nt;
      if (i==1)
      {
        ofstream ofs("totalcatches");
        for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
        {
        
          int rr=fsh.realization_region(i,nt);
          int rp=fsh.realization_period(i,nt);
          int ri=fsh.realization_incident(i,nt);
          ofs << fsh.tot_catch(rr,rp,ri) << " " << sum(fsh.exp_catch(rr,rp,ri)) 
              << exp(fsh.fish_mort_calcs(rr,rp,ri)) << " " 
              << exp(fsh.num_fish(rr,rp)) << endl;
        }
      }
      for (nt=1;nt<=fsh.num_fish_times(i);nt++)
      {
        int rr=fsh.realization_region(i,nt);
        int rp=fsh.realization_period(i,nt);
        int ri=fsh.realization_incident(i,nt);
        of2 << setprecision(3) << setw(12) << value(fsh.exp_catch(rr,rp,ri));
        if (!fsh.data_fish_flags(1,fsh.parent(rr,rp,ri))) 
        {
          of1 << setprecision(3) << setw(12) << value(fsh.tot_catch(rr,rp,ri));
        }
        else
        {
          dvariable totcatch=0.0;
          dvariable sv27=fsh.get_sv_region(rr,27);
          dvariable sv28=fsh.get_sv_region(rr,28);
          if (fsh.sv(28)==3.0)
          {
            totcatch=fsh.len_wt_coff/1000.* (fsh.exp_catch(rr,rp,ri)*
              (pow(fsh.mean_length(rr,rp,ri),3)+
              3.0*elem_prod(fsh.mean_length(rr,rp,ri),fsh.vars(rr,rp,ri))));
          }
          else
          {
            for (int j=1;j<=fsh.nage;j++)
            {  
              totcatch+=fsh.exp_catch(rr,rp,ri,j)*
                normal_length_to_weight(0.5,-3.5,
                fsh.mean_length(rr,rp,ri,j),sqrt(fsh.vars(rr,rp,ri,j)),
//                value(fsh.sv(27)),value(fsh.sv(28)));
                value(sv27),value(sv28));
            }
            totcatch/=1000.;
            if (i==18 && nt==fsh.num_real_fish_times(i)+1)
            {
              ofstream ofs("compare");
              ofs << fsh.mean_length(rr,rp,ri) << endl;
              ofs << sqrt(fsh.vars(rr,rp,ri)) << endl;
              ofs << fsh.exp_catch(rr,rp,ri) << endl;
              ofs << totcatch << endl;
            }
          }
          if (totcatch<=0.0) 
          {
            //cout << "negative totcatch for ir " << rr << " ip " << rp 
              //   << "  fi " << ri << endl;
          }
          of1 << setprecision(3) << setw(12) << value(totcatch);
        }
        of2 << endl;
      }
      of1 << endl;
    }
  }

  // ***************************************************************
  // ***************************************************************
  // ***************************************************************
  {
    ofstream of1("testproj");
    of1 << "# Projected catch by fishery (down) and time (across)" << endl;
    for (i=1;i<=fsh.num_fisheries;i++)
    {
      for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
      {
        int rp=fsh.realization_period(i,nt);
        int rr=fsh.realization_region(i,nt);
        int ri=fsh.realization_incident(i,nt);
        if (rp > fsh.num_real_fish_times(rr))
        {
          if (!fsh.data_fish_flags(1,fsh.parent(rr,rp,ri))) 
          {
            of1 << setprecision(3) << setw(12) << value(fsh.tot_catch(rr,rp,ri));
            dvar_vector& u=fsh.fish_mort(rr,rp,ri);
            xxx(u);
          }
          else
          {
            dvariable totcatch=0.0;
            dvariable sv27=fsh.get_sv_region(rr,27);
            dvariable sv28=fsh.get_sv_region(rr,28);
            if (fsh.sv(28)==3.0)
            {
              totcatch=fsh.len_wt_coff/1000.* (fsh.exp_catch(rr,rp,ri)*
                (pow(fsh.mean_length(rr,rp,ri),3)+
                3.0*elem_prod(fsh.mean_length(rr,rp,ri),fsh.vars(rr,rp,ri))));
            }
            else
            {
              for (int j=1;j<=fsh.nage;j++)
              {  
                totcatch+=fsh.exp_catch(rr,rp,ri,j)*
                  normal_length_to_weight(0.5,-3.5,
                  fsh.mean_length(rr,rp,ri,j),sqrt(fsh.vars(rr,rp,ri,j)),
//                  value(fsh.sv(27)),value(fsh.sv(28)));
                  value(sv27),value(sv28));  //NMD
              }
              totcatch/=1000.;
            }
            if (totcatch<=0.0) 
            {
              //cout << "negative totcatch for ir " << rr << " ip " << rp 
                //   << "  fi " << ri << endl;
            }
            of1 << setprecision(3) << setw(12) << value(totcatch);
          }
        }
      }
      of1 << endl;
    }
  }

  of1 << "# Observed CPUE by fishery (down) and time (across)" << endl;
  for (i=1;i<=fsh.num_fisheries;i++)
  {
    int nft=fsh.num_fish_times(i);
    int nrft=nft;
    if (fsh.do_fishery_projections_flag==1)
    {
      nrft=fsh.num_real_fish_times(i);
    }
//    for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
    for (int nt=1;nt<=nft;nt++)
    {
      int rr=fsh.realization_region(i,nt);
      int rp=fsh.realization_period(i,nt);
      int ri=fsh.realization_incident(i,nt);
      MY_DOUBLE_TYPE tmp=fsh.obs_tot_catch(rr,rp,ri);
      MY_DOUBLE_TYPE tmp1=exp(fsh.effort(rr,rp,ri));
      MY_DOUBLE_TYPE tmp2=tmp/tmp1;
#if !defined(NO_MY_DOUBLE_TYPE)
      if (tmp==-1.0L) tmp=0.0;
#else
      if (tmp==-1.0) tmp=0.0;
#endif
//      if (fsh.fish_flags(i,92) != 0) //NMD_8mar_2022
      if (fsh.fish_flags(i,92) != 0 && nt <= nrft) //NMD_8mar_2022
      {
        tmp2=fsh.survey_cpue_obs(i,nt);
      }
      of1 << setprecision(3) << setw(12) << tmp2;
    }
    of1 << endl;
  }
  of1 << "# Predicted CPUE by fishery (down) and time (across)" << endl;
  for (i=1;i<=fsh.num_fisheries;i++)
  {
    int nft=fsh.num_fish_times(i);
    int nrft=nft;
    if (fsh.do_fishery_projections_flag==1)
    {
      nrft=fsh.num_real_fish_times(i);
    }
//    for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
    for (int nt=1;nt<=nft;nt++)
    {
     // **************************************************************
      int rr=fsh.realization_region(i,nt);
      int rp=fsh.realization_period(i,nt);
      int ri=fsh.realization_incident(i,nt);
      MY_DOUBLE_TYPE tmp=0.0;
      if (!fsh.data_fish_flags(1,fsh.parent(rr,rp,ri))) 
      {
        tmp=value(fsh.tot_catch(rr,rp,ri));
      }
      else
      {
        MY_DOUBLE_TYPE totcatch=0.0;
        dvariable sv27=fsh.get_sv_region(rr,27);
        dvariable sv28=fsh.get_sv_region(rr,28);
        if (fsh.sv(28)==3.0)
        {
          totcatch=fsh.len_wt_coff/1000.* 
            (value(fsh.exp_catch(rr,rp,ri))*
            (pow(value(fsh.mean_length(rr,rp,ri)),3)+
            3.0*elem_prod(value(fsh.mean_length(rr,rp,ri)),
            value(fsh.vars(rr,rp,ri)))));
        }
        else
        {
          for (int j=1;j<=fsh.nage;j++)
          {  
            totcatch+=value(fsh.exp_catch(rr,rp,ri,j)*
              normal_length_to_weight(0.5,-3.5,
              fsh.mean_length(rr,rp,ri,j),sqrt(fsh.vars(rr,rp,ri,j)),
//              value(fsh.sv(27)),value(fsh.sv(28))));
              value(sv27),value(sv28)));  //NMD
          }
          totcatch/=1000.;
        }
        tmp=totcatch;
      }

     // **************************************************************

      MY_DOUBLE_TYPE tmp1=exp(fsh.effort(rr,rp,ri));
      MY_DOUBLE_TYPE tmp2=tmp/tmp1;

      if(fsh.age_flags(92)==2)    // NMD_20may2021
      {
        if (fsh.fish_flags(i,92) == 0) //NMD_8mar_2022
        {
          if (fsh.missing_catch_for_incident_flag(rr,rp,ri)==0
            &&   fsh.missing_effort_by_region_flag(rr,rp,ri)==0)
          {
            tmp1=exp(fsh.log_pred_effort_by_fishery(i,nt));
          }
          else
          {
            tmp1=-1.0;
          }
          tmp2=tmp/tmp1;
        }
        else
        {
          if (nt <= nrft)
          {
            tmp2=fsh.survey_cpue_pred(i,nt);
          }
          else
          {
            if (fsh.missing_catch_for_incident_flag(rr,rp,ri)==0
              &&   fsh.missing_effort_by_region_flag(rr,rp,ri)==0)
            {
              tmp1=exp(fsh.log_pred_effort_by_fishery(i,nt));
            }
            else
            {
              tmp1=-1.0;
            }
              tmp2=tmp/tmp1;
	  }
        }
        /*	
        if (fsh.missing_catch_for_incident_flag(rr,rp,ri)==0
          &&   fsh.missing_effort_by_region_flag(rr,rp,ri)==0)
        {
          tmp1=exp(fsh.log_pred_effort_by_fishery(i,nt));
        }
        else
        {
          tmp1=-1.0;
        }
        tmp2=tmp/tmp1;
        */
      }
        of1 << setprecision(3) << setw(12) << tmp2;
    }
    of1 << endl;
  }
  of1 << "# Yield analysis option: 0=none, 1=Bev&Holt, 2=Pella Tomlinson"
      << endl;
  int yflag=0;
  if (fsh.age_flags(145))
    yflag=1;
  if (fsh.age_flags(150))
    yflag=2;
  of1 << yflag << endl;

  // add this kludge so that alpha and beta don't get changed
  // no that fisheries have been turned off DF feb10 05
  ivector kludge(1,1);
  kludge(1)=1;
  if (fsh.age_flags(145))
  {  
    if (fsh.age_flags(163)==0)
    {
      if (fsh.pmsd==0)
      {
        stock_recruit_bh_steep(fsh,&of1,&kludge);
      }
      else
      {
        // check if this is multi-species or multi-sex
        if (!sum(column(fsh.pmsd->species_flags,2)))
        {
          for (int i=1;i<=fsh.pmsd->num_species;i++)
          {
            fsh.pmsd->current_species=i;
            stock_recruit_bh_steep(fsh,&of1,&kludge);
            debug_flag=0;
          }
        }
        else
        {      // multi-sex
          ivector sf2=column(fsh.pmsd->species_flags,2); //NMD_18Dec2013
          for (int i=1;i<=fsh.pmsd->num_species;i++)
          {
            if (sf2(i))   // get female
            {
              fsh.pmsd->current_species=i;
              stock_recruit_bh_steep(fsh,&of1,&kludge);
            }
          }  //NMD_18Dec2013
        }
      }
    }
    else
    {
      stock_recruit_bh(fsh,&of1,fsh.F_by_age_by_year,&kludge);
    }
  }

  if (fsh.age_flags(150))
  {
    if (!fsh.pmsd)
    {
      fsh.yield_analysis_pt(&of1);
      if (fsh.pmsd)
      {
        for (int is=2;is<=fsh.pmsd->num_species;is++)
        {
          fsh.pmsd->current_species=is;
          fsh.yield_analysis_pt(&of1);
        }
      }
    }
  } 

  fsh.get_equilibrium_structure_for_yield();

  ////////////////////////////////////
  // YT 2018-06-07 Additional code to obtain Fmult from each projection simulation
  if (fsh.projection_sim_index>0 && !pq_flag)
  {
    ofstream& of6=*(fsh.proj_output_files[6]);
    of6<< "#simulation " << fsh.projection_sim_index ;  // YT 2018-04-25
    of6<< "  F_mult@MSY   currentSB/SB@MSY  currentSB" ; //YT 2018-04-25
    of6<< "  SB@MSY"<<endl;
    of6<< setw(10) << fsh.p_Fmmsy << " " <<  fsh.p_sBsBmsy << " ";
    of6<< fsh.p_sBsBmsy.to_double() * fsh.p_sBmsy <<" "<<fsh.p_sBmsy  << endl; // YT 2018-04-25
  }
  
  if (!fsh.pmsd)
    yield_per_recruit_analysis(fsh,&of1,fsh.F_by_age_by_year);

  if (fsh.num_tag_releases) 
  {
//    if(!fsh.age_flags(198))
//    {
      fsh.print_tagging_fit_info(of1);
      fsh.print_tag_return_by_time_at_liberty(of1);
//    }
//    else
//    {
//      fsh.new_print_tagging_fit_info(of1);
//      fsh.new_print_tag_return_by_time_at_liberty(of1);
//    }
  }
  if (!fsh.pmsd && fsh.num_regions > 1) 
  {
    fsh.print_movement_report(of1);
  }
  else if (fsh.pmsd)  //in case of multi-species NMD_jan28-19
  {
    if (fsh.pmsd->num_real_regions > 1)  //in case of single region NMD_jan28-19
    {
      fsh.print_movement_report(of1);
    }
  }  
  //output components of likelihood
  if (!fsh.parest_flags(145)) fsh.print_likelihood_components(of1); //NMD_11dec2023

  {
    int tmult=1;
    if (fsh.age_flags(57)) tmult= fsh.age_flags(57);
    fsh.calculate_the_mean_weight();
    fsh.calculate_the_catch_biomass();
    ofstream ofs3("catch.rep");
    int mmin=fsh.catch_biomass.indexmin();
    int mmax=fsh.catch_biomass.indexmax();
    dvector tmp(1,mmax);
    tmp.initialize();
    dmatrix tmp2(1,mmax,1,fsh.num_fisheries);    //NMD21Nov2013
    tmp2.initialize();
    if (!fsh.pmsd)    //NMD21Nov2013
    {
      for (int i=mmin;i<=mmax;i++)
      {
        tmp(i)=value(fsh.catch_biomass(i));
      }
    }
    else
    {
      for (int j=1; j<=fsh.pmsd->num_species; j++)
      {
        fsh.pmsd->current_species=j;
        fsh.calculate_the_catch_biomass();
        if (j==1)
	{
          for (int i=mmin;i<=mmax;i++)
          {
            tmp(i)=value(fsh.catch_biomass(i));
          }
        }
        else
	{
          for (int i=mmin;i<=mmax;i++)
          {
            tmp(i)+=value(fsh.catch_biomass(i));
          }
        }
        int fmin=(j-1)*fsh.pmsd->num_real_fisheries+1;
        int fmax=j*fsh.pmsd->num_real_fisheries;
        for (int fi=fmin; fi<=fmax; fi++)
	{
          for (int i=mmin;i<=mmax;i++)
	  {
            tmp2(i,fi)=value(fsh.catch_biomass_by_fishery(i,fi));
	  }
        }
      }
    }  //NMD21Nov2013
  
    ofs3 << "# Total catch by year"<<endl;
    ofs3 << tmp << endl;
    ofs3 << "# Catch by year (across) by fishery (down)"<<endl;
    for (int fi=1;fi<=fsh.num_fisheries;fi++)
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (!fsh.pmsd)  //NMD21Nov2013
	{
          ofs3 << value(fsh.catch_biomass_by_fishery(i,fi)) << " ";
        }
        else
	{
          ofs3 << tmp2(i,fi) << " ";
        }     //NMD21Nov2013
      }
      ofs3 << endl;
    }
  }
}

void  dvar_len_fish_stock_history::allocate_plotstuff(void)
{
  int nspp=1;  //NMD_jun27-17
  if (pmsd) nspp=pmsd->num_species;
  ppstf=new plotstuff(num_fisheries,num_regions,nspp,num_fish_times,
      num_fish_periods,num_fish_incidents,gfish_index,fish_flags);
}
void dvar_len_fish_stock_history::print_ssmult_stuff(void)
{
  ofstream ofs("ssmultstuff");
  if (tcs)     //NMD_jun10_19
  {

    ivector& len_rho_group_ptr=tcs->len_rho_group_ptr;
    ivector inv_len_rho_group_ptr=get_inv_group_ptr(len_rho_group_ptr);
    int lmaxg=max(len_rho_group_ptr);
    ofs << " length N  ";
    for (int ig=1;ig<=lmaxg;ig++)
    {
      int gp=inv_len_rho_group_ptr(ig);
      if (fish_flags(gp,67))
      ofs << " "  << exp(fish_pars(14,gp));
    }
    ofs << endl << " length rho  ";
    for (int ig=1;ig<=lmaxg;ig++)
    {
      int gp=inv_len_rho_group_ptr(ig);
      if (fish_flags(gp,67))
      ofs << " " << fish_pars(16,gp);
    }
    ofs << endl << " length var ";
    for (int ig=1;ig<=lmaxg;ig++)
    {
      int gp=inv_len_rho_group_ptr(ig);
      if (fish_flags(gp,82))
      ofs << " " << exp(fish_pars(18,gp));
    }
    ofs << endl << "length sample size covariate coefficient ";
    for (int ig=1;ig<=lmaxg;ig++)
    {
      int gp=inv_len_rho_group_ptr(ig);
      if (fish_flags(gp,85))
      ofs << "  " << exp(fish_pars(20,gp));
    }
    ofs << endl << "length RE correlation" << endl;
    for (int ig=1;ig<=lmaxg;ig++)
    {
      int gp=inv_len_rho_group_ptr(ig);
      if (fish_flags(gp,67) && wpcsa )
      {
        wpcsa[ig]->set_z(fish_pars(16,gp));
        banded_lower_triangular_dvar_matrix vbltd=
              pcsa[ig]->get_ltcholeski_inv(10);
        banded_symmetric_dvar_matrix vbsd2=mult_trans_mult(vbltd);
        print_correlation_matrices(vbsd2,ofs,gp,fish_pars(16,gp));
      }
    }
  }          //NMD_jun10_19
  if (wtcs)
  {
    ivector& wght_rho_group_ptr=wtcs->wght_rho_group_ptr;
    ivector inv_wght_rho_group_ptr=get_inv_group_ptr(wght_rho_group_ptr);
    int wmaxg=max(wght_rho_group_ptr);
    ofs << endl << " weight N  ";
    for (int ig=1;ig<=wmaxg;ig++)
    {
      int gp=inv_wght_rho_group_ptr(ig);
      if (fish_flags(gp,76))
      ofs << " "  << exp(fish_pars(15,gp));
    }
    ofs << endl << " weight rho  ";
    for (int ig=1;ig<=wmaxg;ig++)
    {
      int gp=inv_wght_rho_group_ptr(ig);
      if (fish_flags(gp,76))
      ofs << " " << fish_pars(17,gp);
    }
    ofs << endl << " weight var ";
    for (int ig=1;ig<=wmaxg;ig++)
    {
      int gp=inv_wght_rho_group_ptr(ig);
      if (fish_flags(gp,84))
      ofs << " " << exp(fish_pars(19,gp));
    }
    ofs << endl << "weight sample size covariate coefficient "; 
    for (int ig=1;ig<=wmaxg;ig++)
    {
      int gp=inv_wght_rho_group_ptr(ig);
      if (fish_flags(gp,86))
      ofs << "  " << exp(fish_pars(21,gp));
    }
    /*
    ofs << endl << "length RE correlation" << endl;
    for (int ig=1;ig<=lmaxg;ig++)
    {
      int gp=inv_len_rho_group_ptr(ig);
      if (fish_flags(gp,67))
      {
        wpcsa[ig]->set_z(fish_pars(16,gp));
        banded_lower_triangular_dvar_matrix vbltd=
              pcsa[ig]->get_ltcholeski_inv(10);
        banded_symmetric_dvar_matrix vbsd2=mult_trans_mult(vbltd);
        print_correlation_matrices(vbsd2,ofs,gp,fish_pars(16,gp));
      }
    }
    */
    ofs << "weight RE correlation" << endl;
    for (int ig=1;ig<=wmaxg;ig++)
    {
      int gp=inv_wght_rho_group_ptr(ig);
      if (fish_flags(gp,76) && wpcsa)
      {
        wpcsa[ig]->set_z(fish_pars(17,gp));
        banded_lower_triangular_dvar_matrix vbltd=
              wpcsa[ig]->get_ltcholeski_inv(10);
        banded_symmetric_dvar_matrix vbsd2=mult_trans_mult(vbltd);
        print_correlation_matrices(vbsd2,ofs,gp,fish_pars(17,gp));
      }
    }
  }
  print_length_ssmult_stuff(*this,ofs);
}




void print_length_ssmult_stuff
  (dvar_len_fish_stock_history& fsh,const ofstream & ofs)
{}
 
dmatrix dvar_len_fish_stock_history::get_time(void)
{
  dmatrix time(1,num_fisheries,1,num_fish_times);
  for (int i=1;i<=num_fisheries;i++)
  {
    for (int j=1;j<=num_fish_times(i);j++)
    {
      int ir=realization_region(i,j);
      int ip=realization_period(i,j);
      int yr;
      int mnth;
      if (month_factor!=0 && month_factor!=1)
      {         
        date_struc newdate=get_old_time(year(ir,ip),
          month(ir,ip),week(ir,ip),month_factor,first_time);
        yr=newdate.year+year1-1;
        mnth=newdate.month;
      }
      else
      {
        yr=year(ir,ip)+year1-1;
        mnth=month(ir,ip);
      }
      if (month_1>1)
      {
#if !defined(NO_MY_DOUBLE_TYPE)
        time(i,j)=yr+(mnth-(13-month_1)-0.5L)/12.;
#else
        time(i,j)=yr+(mnth-(13-month_1)-0.5)/12.;
#endif
      }
      else
      {
#if !defined(NO_MY_DOUBLE_TYPE)
        time(i,j)=yr+(mnth-0.5L)/12.;
#else
        time(i,j)=yr+(mnth-0.5)/12.;
#endif
      }
    }
  }
  return time;
}

dmatrix dvar_len_fish_stock_history::get_real_time(void)
{
  dmatrix time(1,num_fisheries,1,num_real_fish_times);
  for (int i=1;i<=num_fisheries;i++)
  {
    for (int j=1;j<=num_real_fish_times(i);j++)
    {
      int ir=realization_region(i,j);
      int ip=realization_period(i,j);
      int yr;
      int mnth;
      if (month_factor!=0 && month_factor!=1)
      {         
        date_struc newdate=get_old_time(year(ir,ip),
          month(ir,ip),week(ir,ip),month_factor,first_time);
        yr=newdate.year+year1-1;
        mnth=newdate.month;
      }
      else
      {
        yr=year(ir,ip)+year1-1;
        mnth=month(ir,ip);
      }
      if (month_1>1)
      {
#if !defined(NO_MY_DOUBLE_TYPE)
        time(i,j)=yr+(mnth-(13-month_1)-0.5L)/12.;
#else
        time(i,j)=yr+(mnth-(13-month_1)-0.5)/12.;
#endif
      }
      else
      {
#if !defined(NO_MY_DOUBLE_TYPE)
        time(i,j)=yr+(mnth-0.5L)/12.;
#else
        time(i,j)=yr+(mnth-0.5)/12.;
#endif
      }
    }
  }
  return time;
}
