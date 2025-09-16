/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
void write_bs_selcoff(par_ofstream& pof,dvar4_array& bsc,int nage);
void read_bs_selcoff(par_ofstream& pof,dvar4_array& bsc,int nage);

//void write_check_number(const dvar_fish_stock_history& fsh,par_ofstream & pof,char * s)
void write_check_number(const dvar_fish_stock_history& fsh,par_ofstream & pof,const char * s)
{
  if (fsh.parest_flags(195))
  {
    pof << s << endl;
  }
}

void write_bs_selcoff(par_ofstream& pof,dvar4_array& bsc,int nage)
{
  int imin=bsc.indexmin();
  int imax=bsc.indexmax();
  dvector tmp(1,nage);
  for (int i=imin;i<=imax;i++)
  {
    int jmin=bsc(i).indexmin();
    int jmax=bsc(i).indexmax();
    for (int j=jmin;j<=jmax;j++)
    {
      int kmin=bsc(i,j).indexmin();
      int kmax=bsc(i,j).indexmax();
      for (int k=kmin;k<=kmax;k++)
      {
         tmp.initialize();
         int lmin=bsc(i,j,k).indexmin();
         int lmax=bsc(i,j,k).indexmax();
         tmp(lmin,lmax)=value(bsc(i,j,k));
         pof << tmp;
      }
    }
  }
}
void read_bs_selcoff(cifstream& cif,dvar4_array& bsc,int nage)
{
  int imin=bsc.indexmin();
  int imax=bsc.indexmax();
  dvector tmp(1,nage);
  for (int i=imin;i<=imax;i++)
  {
    int jmin=bsc(i).indexmin();
    int jmax=bsc(i).indexmax();
    for (int j=jmin;j<=jmax;j++)
    {
      int kmin=bsc(i,j).indexmin();
      int kmax=bsc(i,j).indexmax();
      for (int k=kmin;k<=kmax;k++)
      {
         tmp.initialize();
         int lmin=bsc(i,j,k).indexmin();
         int lmax=bsc(i,j,k).indexmax();
         cif >> tmp;
         bsc(i,j,k)= tmp(lmin,lmax);
      }
    }
  }
}
void read_bs_selcoff(cifstream& cif,d4_array& bsc,int nage)
{
  int imin=bsc.indexmin();
  int imax=bsc.indexmax();
  dvector tmp(1,nage);
  for (int i=imin;i<=imax;i++)
  {
    int jmin=bsc(i).indexmin();
    int jmax=bsc(i).indexmax();
    for (int j=jmin;j<=jmax;j++)
    {
      int kmin=bsc(i,j).indexmin();
      int kmax=bsc(i,j).indexmax();
      for (int k=kmin;k<=kmax;k++)
      {
         tmp.initialize();
         int lmin=bsc(i,j,k).indexmin();
         int lmax=bsc(i,j,k).indexmax();
         cif >> tmp;
         bsc(i,j,k)= tmp(lmin,lmax);
      }
    }
  }
}
         

void write_p1(ostream& ofs, dvar_fish_stock_history& fsh);
void write_p2(ostream& ofs, dvar_fish_stock_history& fsh);
void write_p3(ostream& ofs, dvar_fish_stock_history& fsh);

void write_report(ostream& ofs,dvar_fish_stock_history& fsh)
{
  ofs << fsh;
  ofs << "obs_catch - catch" << endl <<endl;
  for (int i=fsh.obs_catch.slicemin();i<=fsh.obs_catch.slicemax();i++)
  {
    ofs << fsh.obs_catch(i) - exp(fsh.catch(i)) << endl;
  }
  ofs << endl;
  ofs << "obs_tot_catch - tot_catch" << endl <<endl;
  ofs << fsh.obs_tot_catch - fsh.tot_catch << endl << endl;
}

ostream& operator << (ostream& ofs, dvar_fish_stock_history& fsh)
{
  write_p1(ofs,fsh);
  write_p2(ofs,fsh);
  write_p3(ofs,fsh);
  return ofs;
}


    par_ofstream& operator << (par_ofstream& pof, dvar_fish_stock_history& fsh)
    {
      fsh.parest_flags(200)=max(fsh.parest_flags(200),1000);

      pof <<  "# age flags"<<endl;
      //pof << std::noshowpoint << "# age flags"<<endl;
      if (fsh.age_flags(52)==0)
      {
        pof << fsh.age_flags << endl;
        if (fsh.pmsd)
        {
          if (fsh.pmsd->num_species>1)
          {
            pof <<  "# multi-species age flags"<<endl;
            {
              pof <<  fsh.pmsd->age_flags.sub(2,fsh.pmsd->num_species) <<endl;
            }
          }
        }
      }
      else
      {
        ivector tmp(fsh.age_flags.indexmin(),fsh.age_flags.indexmax());
        tmp=fsh.age_flags;
        if (fsh.age_flags(41)==0)
        {
          tmp(41)=1;
          tmp(32)=1;
        }
        else
        {
          tmp(41)=0;
          tmp(32)=0;
        }
        tmp(52)=0;
        pof << tmp << endl;
      }
      write_check_number(fsh,pof, " 871" );
      pof << "# fish flags"<<endl;
      pof << setw(4) << fsh.fish_flags << endl;

      //pof  << "# kludged_equilib_coffs" << endl;
      //pof << setw(4) << fsh.kludged_equilib_coffs << endl;

      write_check_number(fsh,pof, " 872" );
      if (fsh.num_tag_releases)
      {
        pof << "# tag flags"<<endl;
        if (fsh.pmsd)
        {
          pof << fsh.true_tag_flags << endl;
        }
        else
        {
          pof << fsh.tag_flags << endl;
        }
        write_check_number(fsh,pof, " 1872");
        pof << "# tagmort"<<endl;
        pof << setscientific()<< setprecision(14) 
            << fsh.tagmort << endl;
        

        pof << setprecision(14) << endl;
        pof << "# tag fish rep" << endl;
        pof << setprecision(14) << endl;
        pof << fsh.tag_fish_rep << endl;
        pof << "# tag fish rep group flags" << endl;
        pof <<  fsh.tag_fish_rep_group_flags << endl;
        pof << "# tag_fish_rep active flags" << endl;
        pof << fsh.tag_fish_rep_active_flags << endl;
        pof << "# tag_fish_rep target" << endl;
        pof << fsh.tag_fish_rep_target << endl;
        pof << "# tag_fish_rep penalty" << endl;
        pof << fsh.tag_fish_rep_penalty << endl;

      }

      write_check_number(fsh,pof, " 2872");
      pof << "# region control flags "<<endl;
      pof << setscientific()<< setprecision(14) << fsh.region_flags << endl;

      write_check_number(fsh,pof, " 873");
      if (fsh.pmsd)
      {
        pof << "# species flags "<<endl;
        pof << setscientific()<< setprecision(14) << fsh.pmsd->species_flags << endl;
      }
      write_check_number(fsh,pof, " 1874");

      pof << "# percent maturity "<<endl;
      pof << setscientific()<< fsh.pmature <<endl;
      write_check_number(fsh,pof, " 874" );
      if (fsh.pmsd)
      {
        pof << "#multi-species percent maturity "<<endl;
        pof << setscientific()<< fsh.pmsd->pmature <<endl;
      }
      write_check_number(fsh,pof, " 875" );
      pof << "# total populations scaling parameter "<<endl;
      pof << setscientific()<< fsh.totpop_coff <<endl;
      if (fsh.pmsd)
      {
        pof << "# multi-species total populations scaling parameter "<<endl;
        pof << setscientific()<< fsh.pmsd->totpop_coff <<endl;
      }
      write_check_number(fsh,pof, " 876");
      pof << "# implicit total populations scaling parameter "<<endl;
      pof << setscientific()<< fsh.implicit_totpop_coff << endl;
      pof << "# rec init pop level difference "<<endl;
      pof << setscientific()<< fsh.rec_init_diff <<endl;
      if (fsh.pmsd)
      {
        pof << "# multi-species rec init pop level difference "<<endl;
        pof << setscientific() << fsh.pmsd->rec_init_diff << endl;
      }
      write_check_number(fsh,pof, " 877");
      pof << "# recruitment times "<<endl;
      pof << setscientific()<< fsh.rec_times <<endl;
      pof << "# relative recruitment "<<endl;
      if (fsh.age_flags(52)==0)
      {
        pof << "#" << setscientific() << setprecision(14) << endl;
        pof <<  exp(fsh.recr) << endl;
        if (fsh.pmsd && fsh.pmsd->num_species>1)
        {
          pof << "# multi-species relative recruitment "<<endl;
          for (int is=2;is<=fsh.pmsd->num_species;is++)
          {
            pof <<  exp(fsh.pmsd->recr(is)) << endl;
          }
        }
      }
      else
      {
        pof << "#" << setscientific() << setprecision(14) << endl;
        pof << exp(fsh.tmprecr) << endl;
        fsh.pmsd_error();
      }
      write_check_number(fsh,pof, " 878" );
      pof << "# Lambdas for augmented Lagrangian" << endl;
      pof << "#" << setscientific() << setprecision(14) << endl;
      if (fsh.parest_flags(399))
      {
        fsh.update_lagrange_lambda();
      }
      pof << fsh.lagrange_lambda << endl;
      pof << "# Other lambdas for augmented Lagrangian" << endl;
      pof << "#" << setscientific() << setprecision(14) << endl;
      if (fsh.parest_flags(349))
      {
        fsh.update_lagrange_avg_biomass();
      }
      pof << fsh.other_lagrange_lambda << endl;
      pof << "# Reporting rate dev coffs"<<endl;
      pof << "# " << setscientific() << setprecision(14) << endl;

      if (allocated(fsh.rep_dev_coffs))
        pof << fsh.rep_dev_coffs << endl;

      pof << "# availability coffs "<<endl;
      pof << "# " << setscientific() << setprecision(14) << endl;
      pof << fsh.avail_coff << endl;
      write_check_number(fsh,pof, " 879");
      if (fsh.parest_flags(155)>0)
      {
	pof << "# annual coffs for relative recruitment "<<endl;
        pof << "# " << setscientific() << setprecision(17) 
            << endl << fsh.yearly_recr_all << endl;
        write_check_number(fsh,pof, " 6875");
	pof << "# orthogonal poly coffs for relative recruitment "<<endl;
        pof << "# " << setscientific() << setprecision(17) 
            << endl << fsh.orth_recr_all << endl;
        write_check_number(fsh,pof, " 5875" );
        if (fsh.pmsd && fsh.pmsd->num_species>1)
        {
          for (int is=2;is<=fsh.pmsd->num_species;is++)
          {
            if (fsh.pmsd->parest_flags(is,155)>0)
            {
	      pof << "# annual coffs for relative recruitment for species "
                  << is << endl;
              pof << "# " << setscientific() << setprecision(17) 
                  << endl << fsh.pmsd->yearly_recr_all(is) << endl;
	      pof << "# orthogonal poly coffs for relative recruitment "<<endl;
              pof << "# " << setscientific() << setprecision(17) 
                  << endl << fsh.pmsd->orth_recr_all(is) << endl;
            }
          }
        }
        else
        {
          //ad_exit(1);
        }
      }
      write_check_number(fsh,pof, " 4875" );
      pof << "# relative initial population "<<endl;
      if (fsh.age_flags(52)==0)
      {
        pof << setscientific() << setprecision(14) 
            << endl <<  exp(fsh.initpop) << endl;
      }
      else
      {
        pof << setscientific() << setprecision(14) 
            << endl <<  exp(fsh.tmpinitpop) << endl;
      }
      write_check_number(fsh,pof, " 3875" );
      pof << "# fishery selectivity "<< endl;
      pof << setprecision(14);
      
      // ***********************************************************
      //pof << fsh.bs_selcoff << endl;
      write_bs_selcoff(pof,fsh.bs_selcoff,fsh.nage);
      // ***********************************************************
      pof << "# age-dependent component of fishery selectivity "<< endl;
      pof << setprecision(14);
      pof << fsh.ageselcoff << endl;
      pof << "# natural mortality coefficient "<<endl;
      pof << "# " << setscientific() << setprecision(14) << endl;
      pof << setscientific()<< setprecision(14);
      pof <<exp(fsh.nat_mort_coff) << endl;
//NMD 02Nov2011
      if (fsh.pmsd)
      {
        pof << "#multi-species natural mortality coefficient "<< endl;
        pof << exp(fsh.pmsd->nat_mort_coff) << endl;
      }
//NMD 02Nov2011

      pof << "# average catchability coefficients "<<endl;
      pof << "# " << endl;
      pof << setscientific()<< setprecision(14)<< endl;
      pof <<exp(fsh.q0) << endl;


      pof << "# initial trend in catchability coefficients "<<endl;
      pof << "# " << endl;
      pof << setscientific()<< setprecision(14) <<fsh.q1 << endl;

      pof << "# q0_miss " << endl;
      pof << "# " << endl;
      pof << setscientific()<< setprecision(14) << fsh.q0_miss << endl;

      write_check_number(fsh,pof, " 3876");
      if (fsh.missing_catch_flag)
      {
        pof << "# fm_level_devs "<< endl;
        pof << setscientific()<< setprecision(14) << fsh.fm_level_devs << endl;
      }

      fsh.Dad2=fsh.Dad;    //NMD_20Nov_2018  - store current Dad for output

      int af184=fsh.age_flags(184);
      fsh.rationalize_movement_coffs(af184);
      pof << "# movement map "<<endl;
      pof << setscientific()<< setprecision(14) <<fsh.move_map << endl;
      pof << "# diff_coffs movement coefficients "<<endl;
      pof <<  setprecision(14) <<fsh.diff_coffs << endl;
      pof << "# xdiff_coffs movement coefficients "<<endl;
      pof <<  setprecision(14) <<fsh.xdiff_coffs << endl;
      pof << "# y1diff_coffs movement coefficients "<<endl;
      pof <<  setprecision(14) <<fsh.y1diff_coffs << endl;
      pof << "# y2diff_coffs movement coefficients "<<endl;
      pof <<  setprecision(14) <<fsh.y2diff_coffs << endl;
      pof << "# zdiff_coffs movement coefficients "<<endl;
      pof <<  setprecision(14) <<fsh.zdiff_coffs << endl;
      pof << "# movement matrices "<<endl;

      int np=fsh.mo.num_periods();
      // ****************************************
      // ****************************************
      /*
      {

        ofstream ofs("dadpar");
        for (int i=1;i<=4;i++)
        {
          int j=i%4+1;
          ofs << fsh.Dad(j) << endl;
        }
        ofs << "Dad" << endl;
        ofs << fsh.Dad << endl; 
        ofs << "diff_coffs" << endl;
        ofs << fsh.diff_coffs << endl;
        ofs << "diff_coffs2" << endl;
        ofs << fsh.diff_coffs2 << endl;
        ofs << "diff_coffs3" << endl;
        ofs << fsh.diff_coffs3 << endl;
        //ad_exit(1);
      }
      */
      // ****************************************
      // ****************************************

      for (int ii=1;ii<=np;ii++)
      {
        for (int j=1;j<=fsh.nage;j++)
        {
          pof << "# Movement period " << ii << "  age class " << j << endl;
//          pof <<  setprecision(14) <<fsh.Dad(ii,j) << endl;
          pof <<  setprecision(14) <<fsh.Dad2(ii,j) << endl;  //NMD_20Nov2018
        }
      }
      pof << "# age dependent movement coefficients "<<endl;
      pof <<  setprecision(14) <<fsh.diff_coffs2 << endl;
      pof << "# nonlinear movement coefficients "<<endl;
      pof <<  setprecision(14) <<fsh.diff_coffs3 << endl;
      pof << "# Movement coefficients priors "<<endl;
      pof <<  setprecision(14) <<fsh.diff_coffs_prior << endl;
      pof << "# age dependent movement coefficients priors"<<endl;
      pof <<  setprecision(14) <<fsh.diff_coffs2_prior << endl;
      pof << "# nonlinear movement coefficients priors"<<endl;
      pof <<  setprecision(14) <<fsh.diff_coffs3_prior << endl;
      pof << "# regional recruitment variation "<<endl;
      pof << setprecision(14);
      pof << fsh.region_rec_diff_coffs << endl;
      pof << "# effort deviation coefficients "<<endl;
      pof << setscientific()<< setprecision(14);
      pof << fsh.effort_dev_coffs << endl;
      pof << "# correlation in selectivity deviations "<<endl;
      pof << fsh.corr_wy << endl;
      pof << "# extra fishery parameters "<<endl;
      pof << "# " << setprecision(14)<< endl;
      pof << setscientific()<< setprecision(14) << endl;
      pof << fsh.fish_pars << endl;

      pof << "# fsh.implicit_fm_level_regression_pars "<<endl;
      pof << "# " << setprecision(14)<< endl;
      pof << setscientific()<< setprecision(14) << endl;
      pof << fsh.implicit_fm_level_regression_pars << endl;

      pof << "# species parameters "<<endl;
      pof << setscientific()<< setprecision(14) << endl;
      pof << fsh.species_pars << endl;
      pof << "# seasonal_catchability_pars "<<endl;
      pof << fsh.seasonal_catchability_pars << endl;

      pof << "# age-class related parameters (age_pars)"<<endl;
      pof << "# " << endl;
      pof << setscientific()<< setprecision(14) 
          << fsh.age_pars << endl;
//NMD 04Nov2011
      if (fsh.pmsd)
      {
        for (int is=2;is<=fsh.pmsd->num_species;is++)
        {
          pof << "#multi-species age-class related parameters"
                 " (age_pars), species: " << is << endl;
          pof << "# " << endl;
          pof << fsh.pmsd->age_pars(is) <<endl;
		}
      }
//NMD 04Nov2011
      pof << "# region parameters "<<endl;
      pof << setprecision(14);
      pof << fsh.region_pars << endl;

      pof << "# catchability deviation coefficients "<<endl;
      pof << "# " << endl;
      pof << setprecision(14);
      pof << fsh.catch_dev_coffs << endl;
      {
	d3_array tarr(1,fsh.num_fisheries,2,fsh.num_fish_times,1,fsh.nage);
	for (int i=1;i<=fsh.num_fisheries;i++)
	{
	  for (int nt=2;nt<=fsh.num_fish_times(i);nt++)
	  {
	    for (int j=1;j<=fsh.nage;j++)
	    {
	      tarr(i,nt,j)=value(fsh.delta2(i,nt,min(j,
                max(1,fsh.fish_flags(i,3)))));
	    }
	  }
	}
	pof << "# selectivity deviation coefficients "<<endl;
        pof << "# " << setprecision(14) << endl;
	pof << setprecision(14) << endl;
	pof <<tarr << endl;
      }
      pof << "# sel_dev_coffs "<<endl;
      for (int i=1;i<=fsh.num_fisheries;i++)
      {
        for (int j=1;j<=fsh.num_fish_times(i);j++)
        {
          pof << setprecision(14);
          pof << fsh.sel_dev_coffs(i,j);
        }
        pof << endl;
      }
        
      pof << "# year_flags "<<endl;
      pof << fsh.year_flags << endl;
      pof << "# season_flags "<<endl;
      pof << fsh.season_flags << endl;
      return pof;
    }

void write_p1(ostream& ofs, dvar_fish_stock_history& fsh)
{
  ofs <<     "parest_flags" << endl <<endl;
  ofs <<     fsh.parest_flags << endl <<endl;
  ofs <<     "age_flags" << endl <<endl;
  ofs <<     fsh.age_flags << endl <<endl;
  ofs <<     "fish_flags" << endl <<endl;
  ofs <<     fsh.fish_flags << endl <<endl;
  ofs <<     "nage" << endl <<endl;
  ofs <<     fsh.nage << endl <<endl;
  ofs <<     "nyears" << endl <<endl;
  ofs <<     fsh.nyears << endl <<endl;
  ofs <<     "num_fisheries" << endl <<endl;
  ofs <<     fsh.num_fisheries << endl <<endl;
  ofs <<     "num_fish_periods" << endl <<endl;
  ofs <<     fsh.num_fish_periods << endl <<endl;
  ofs <<     "num_fish_incidents" << endl <<endl;  
  ofs <<     fsh.num_fish_incidents << endl <<endl;  
  ofs <<     "num_fish_times" << endl <<endl;      
  ofs <<     fsh.num_fish_times << endl <<endl;      
  ofs <<     "parent" << endl <<endl;
  ofs <<     fsh.parent << endl <<endl;              
  ofs <<     "year" << endl <<endl;                
  ofs <<     fsh.year << endl <<endl;                
  ofs <<     "fraction" << endl <<endl;            
  ofs <<     fsh.fraction << endl <<endl;            
  ofs <<     "totpop_coff" << endl <<endl;              
}

void write_p2(ostream& ofs, dvar_fish_stock_history& fsh)
{
  ofs <<     fsh.totpop_coff << endl <<endl;              
  ofs <<     "recr" << endl <<endl;                
  ofs <<     fsh.recr << endl <<endl;
  ofs <<     "initpop" << endl <<endl;
  ofs <<     fsh.initpop << endl <<endl;
  ofs <<     "num_fish" << endl <<endl;            
  ofs <<     fsh.num_fish << endl <<endl;
  ofs <<     "catch" << endl <<endl;
  ofs <<     exp(fsh.catch) << endl <<endl;
  ofs <<     "tot_catch" << endl <<endl;
  ofs <<     fsh.tot_catch << endl <<endl;
  ofs <<     "prop" << endl <<endl;
  ofs <<     fsh.prop << endl <<endl;
  ofs <<     "incident_sel" << endl <<endl;
  ofs <<     fsh.incident_sel << endl <<endl;
  ofs <<     "fishery_sel" << endl <<endl;
  ofs <<     fsh.selcoff << endl <<endl;
  ofs <<     "selcoff" << endl <<endl;
  ofs <<     fsh.selcoff << endl <<endl;
  ofs <<     "survival" << endl <<endl;
  ofs <<     fsh.survival << endl <<endl;
  ofs <<     "fish_mort" << endl <<endl;
  ofs <<     fsh.fish_mort << endl <<endl;
  ofs <<     "tot_mort" << endl <<endl;
  ofs <<     fsh.tot_mort << endl <<endl;
  ofs <<     "nat_mort" << endl <<endl;            
  ofs <<     fsh.nat_mort << endl <<endl;
  ofs <<     "catchability" << endl <<endl;
  ofs <<     fsh.catchability << endl <<endl;
  ofs <<     "effort" << endl <<endl;              
  ofs <<     fsh.effort << endl <<endl;              
  ofs <<     "month" << endl <<endl;
  ofs <<     fsh.month << endl <<endl;               
  ofs <<     "week" << endl <<endl;                
  ofs <<     fsh.week << endl <<endl;                
  ofs <<     "nat_mort_coff" << endl <<endl;       
  ofs <<     fsh.nat_mort_coff << endl <<endl;       
  ofs <<     "catch_init" << endl <<endl;          
  ofs <<     fsh.catch_init << endl <<endl;          
}

void write_p3(ostream& ofs, dvar_fish_stock_history& fsh)
{
  ofs <<     "q0" << endl <<endl;                  
  ofs <<     fsh.q0 << endl <<endl;                  
  ofs <<     "q1" << endl <<endl;                  
  ofs <<     fsh.q1 << endl <<endl;                  
  ofs <<     "effort_dev_coffs" << endl <<endl;
  ofs <<     fsh.effort_dev_coffs << endl <<endl;
  ofs <<     "effort_devs" << endl <<endl;         
  ofs <<     fsh.effort_devs << endl <<endl;
  ofs <<     "sel_dev_coffs" << endl <<endl;
  ofs <<     fsh.sel_dev_coffs << endl <<endl;       
  ofs <<     "sel_devs" << endl <<endl;       
  ofs <<     fsh.sel_devs << endl <<endl;
  ofs <<     "corr_wy" << endl <<endl;             
  ofs <<     fsh.corr_wy << endl <<endl;             
  ofs <<     "corr_by" << endl <<endl;             
  ofs <<     fsh.corr_by << endl <<endl;             
  ofs <<     "corr_wc" << endl <<endl;             
  ofs <<     fsh.corr_wc << endl <<endl;             
  ofs <<     "corr_eff" << endl <<endl;            
  ofs <<     fsh.corr_eff << endl <<endl;            
  ofs <<     "obs_tot_catch" << endl <<endl;
  ofs <<     fsh.obs_tot_catch << endl <<endl;
  ofs <<     "obs_prop" << endl <<endl;            
  ofs <<     fsh.obs_prop << endl <<endl;            
  ofs <<     "obs_catch" << endl <<endl;
  ofs <<     fsh.obs_catch << endl <<endl;
}

void dvar_fish_stock_history::update_lagrange_lambda(void)
{
  if (no_lagrangian_update==0)
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {
        lagrange_lambda(ir,ip,fi)-=lagrange_mu*value(lagrange_c(ir,ip,fi));
      }
    }
  }
}

void dvar_fish_stock_history::update_lagrange_avg_biomass(void)
{
  if (no_lagrangian_update==0)
  other_lagrange_lambda(1)-=
    other_lagrange_mu(1)*value(other_lagrange_c(1));
}

#undef HOME_VERSION
