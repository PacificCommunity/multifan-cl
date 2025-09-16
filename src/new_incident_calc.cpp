/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"

  
void block_error();
extern dvar_vector * const_vb;
extern dvar_vector * const_var;

void dvar_len_fish_stock_history::new_incident_selectivity_calc
  (d3_array& len_sample_size,dvar_vector& vb_coff,dvar_vector& var_coff)
{
  ivector nd=column(fish_flags,3);
  ivector ff16=column(fish_flags,16);
  ivector ff26=column(fish_flags,26);
  ivector ff57=column(fish_flags,57);
  ivector ff74=column(fish_flags,74);
  ivector ff71=column(fish_flags,71);
  ivector ff75=column(fish_flags,75);
  int i,j;
  for (i=1;i<=num_fisheries;i++)
  {
    if (nd(i)==0) nd(i)=nage_by_fishery(i);
  }
  ivector tmplength(1,num_fisheries);
  for (i=1;i<=num_fisheries;i++)
  {
    if (ff26(i)==3)
      tmplength(i)=nlint;
    else
      tmplength(i)=nd(i);
  }
  
  dvar_matrix tempsel(1,num_fisheries,1,tmplength);
  dvar4_array bs_tempsel(1,num_fisheries,1,ff74,1,ff71+1,1,tmplength);
  tempsel.initialize();
  dvariable rho=0.0;
  //dvar_matrix tlength(1,num_fisheries,1,nd);
  get_scaled_single_or_two_species_lengths(); 
  dvar_matrix lbsel(1,num_fisheries,1,nage_by_fishery);
  
  if (age_flags(185)==0)
  {
    lbsel_calc(lbsel,*this,vb_coff,var_coff);
  }
  else
  {
    lbsel_calc(lbsel,*this,*const_vb,*const_var);
  }
 
  logistic_sel_calc(*this);
  MY_DOUBLE_TYPE min_select=.000001;
  //double rnage=sqrt(double(nage));
  MY_DOUBLE_TYPE one_plus=1.+min_select;
  MY_DOUBLE_TYPE alpha=min_select/exp(1.);
  dvar_vector& sel_corr = corr_wy;
  dvar_vector root_sel_corr=sqr(1.0-square(sel_corr));
  int num_age=0;

  no_splines();

  splines();


  if (ff263flag)
  {
    all_length_dist_calcs_long();
    if (nwint>0)
    {
      all_weight_dist_calcs_long();
    }
//    {
//      ofstream ofs("aws");
//      ofs << age_weight_fishery_size_dist(2,1,1) << endl;
//      cout << "here" << endl;
//    }
  }
      
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)  // Loop over fishing
      {
        int pi=parent(ir,ip,fi);
        int is=sseason(ir,ip,fi);
        int ib=bblock(ir,ip,fi);
        if (!ff263flag)     //NMD_20dec2022
        {	  
          if (!ff75(pi))
          {
            int mmin=bstempsel(pi,is,ib).indexmin();
            int mmax=bstempsel(pi,is,ib).indexmax();
            //incident_sel(ir,ip,fi)(1,nd(pi))=bstempsel(pi,is,ib);
            // YYYY
            if (mmax > incident_sel(ir,ip,fi).indexmax())
            {
              incident_sel(ir,ip,fi).deallocate();
              incident_sel(ir,ip,fi).allocate(mmin,mmax);
            }
            incident_sel(ir,ip,fi)(mmin,mmax)=bstempsel(pi,is,ib);
          }
          else
          {
            incident_sel(ir,ip,fi)(ff75(pi)+1,nd(pi))=bstempsel(pi,is,ib);
            incident_sel(ir,ip,fi)(1,ff75(pi))=-20.;
          }
        }  // NMD_20dec2022
        if (nd(pi)<nage_by_region(ir))
        {
          if (!ff263flag)
          {
            if (ff57(pi)>0)
            {
              incident_sel(ir,ip,fi)(nd(pi),nage_by_region(ir))=
              bstempsel(pi,is,ib,nd(pi));
            }
            else
            {
              if (ff16(pi)==1)
              {
                incident_sel(ir,ip,fi)(nd(pi),nage_by_region(ir))=
                bstempsel(pi,is,ib,nd(pi));
              }
              else if (ff16(pi)==2)
              {
                incident_sel(ir,ip,fi)(nd(pi),nage_by_region(ir))= -20.0;
              }
              else
              {
                cerr << "ERROR: incorrect fish_flags(16) setting, fishery: "
                     << pi << endl;
                ad_exit(1);
              }
            }
          }
          else
          {
            if (ff16(pi)==0)
            {
              incident_sel(ir,ip,fi)(nd(pi),nage_by_region(ir))=
                incident_sel(ir,ip,fi,nd(pi));
            }
            else if (ff16(pi)==1)
            {
//            do nothing - managed by spline penalty
            }
            else if (ff16(pi)==2)
            {
              incident_sel(ir,ip,fi)(nd(pi),nage_by_region(ir))= -20.0;
            }
            else
            {
              cerr << "ERROR: incorrect fish_flags(16) setting, fishery: "
                   << pi << endl;
              ad_exit(1);
            }
          }
        }
////////////////
        if (age_flags(92)>0 && fish_flags(pi,19)>0)
        {
          int num_ages=min(fish_flags(pi,19),nage_by_region(ir));
	  int nt;
	  nt=fish_times(ir,ip,fi);
          if (len_sample_size(ir,ip,fi)>0 || wght_sample_size(ir,ip,fi)>0)
          {
            for (int jj=1;jj<=num_ages;jj++)
            {
              incident_sel(ir,ip,fi)+=sel_dev_coffs(pi,nt,jj);
            }
          }
        }
////////////////
      }    //ND_26feb2015
    }
  }
}

