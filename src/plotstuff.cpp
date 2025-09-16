/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
plotstuff::plotstuff(int num_fisheries,int num_regions,int nspp,
  const ivector& num_fish_times,
  const ivector& num_fish_periods,
  const imatrix& num_fish_incidents,
  const imatrix & gfish_index,
  const imatrix & fish_flags)
{
   lencontrib.allocate(1,num_fisheries);
   wghtcontrib.allocate(1,num_fisheries);
   lencontrib_by_realization.allocate(1,num_regions,
     1,num_fish_periods,1,num_fish_incidents);
   wghtcontrib_by_realization.allocate(1,num_regions,
     1,num_fish_periods,1,num_fish_incidents);
   tot_catch_like_by_realization.allocate(1,num_regions,
     1,num_fish_periods,1,num_fish_incidents);
   effort_dev_penalty_by_fishery.allocate(1,num_fisheries);
   survey_index_like_by_fishery.allocate(1,num_fisheries);
   catchability_dev_penalty_by_fishery.allocate(1,num_fisheries);
   bh_steep_contribution.allocate(1,nspp);
   orth_poly_penalty_by_level.allocate(1,nspp,1,4);
   spline_sel_bound_pen.allocate(1,num_fisheries);
   nonspline_sel_curve_pen.allocate(1,num_fisheries);

   ivector grouping=column(fish_flags,29);
   int tmp_ngroups=max(grouping);
   catchability_dev_penalty_by_group.allocate(1,tmp_ngroups);
   grouping=column(fish_flags,99);
   tmp_ngroups=max(grouping);
   survey_index_like_by_group.allocate(1,tmp_ngroups);

   lencontrib = 0.;
   wghtcontrib = 0.;
//   bh_steep_contribution=0.0;
   lencontrib_by_realization.initialize();
   wghtcontrib_by_realization.initialize();
   tot_catch_like_by_realization.initialize();
   effort_dev_penalty_by_fishery.initialize();
   survey_index_like_by_fishery.initialize();
   survey_index_like_by_group.initialize();
   catchability_dev_penalty_by_fishery.initialize();
   catchability_dev_penalty_by_group.initialize();
   bh_steep_contribution.initialize();  //NMD_jun27-17
   orth_poly_penalty_by_level.initialize();  //NMD_05dec2022
   spline_sel_bound_pen.initialize();   //NMD_21nov2023
   nonspline_sel_curve_pen.initialize();
   reset_pen.initialize();
   do_every_pen.initialize();
   reg_recr_pen.initialize();
   diff_coff_pen.initialize();
   seasonal_q_pen.initialize();
   tag_rep_rate_pen.initialize();
   region_recr_dev_pen.initialize();
   mean_region_recr_pen.initialize();
   init_recr_pen.initialize();
   init_agecomp_curve_pen.initialize();
   temporal_recr_dev_pen.initialize();
   norm_temporal_recr_pen.initialize();
   recr_trend_pen.initialize();
   recr_curve_pen.initialize();
   incid_sel_curve_pen.initialize();
   mean_sel_pen.initialize();
   time_block_sel_mean_pen.initialize();
   sel_devs_pen.initialize();
   sel_form_pen.initialize();
   kludge_surv_pen.initialize();
   impl_fm_level_regr_pars_pen.initialize();
   impl_fm_level_regr_pen.initialize();
   fmort_max_pen.initialize();
   fmort_agelag_max_pen.initialize();
   mean_len_constr_pen.initialize();
}

void dvar_len_fish_stock_history::print_likelihood_components
  (ofstream& of2)
{

  ofstream of1("test_plot_output");

  of1 << "# Kludged equilibrium survival penalty contribution" << endl  
      << " " << ppstf->kludge_surv_pen << endl;

  of1 << "# reset penalty contribution" << endl  
      << " " << ppstf->reset_pen << endl;

  of1 << "# do_everything_calc penalty contribution" << endl  
      << " " << setprecision(12) << ppstf->do_every_pen << endl;

  of1 << "# BH_steep contribution" << endl  
      << ppstf->bh_steep_contribution << endl;

  of1 << "# Regional recruitment distribution penalty contribution" << endl  
      << " " << ppstf->reg_recr_pen << endl;

  of1 << "# Movement diffusion coefficients penalty contribution" << endl  
      << " " << ppstf->diff_coff_pen << endl;

  of1 << "# Implicit fm_level regression coefficients penalty contribution" << endl  
      << " " << ppstf->impl_fm_level_regr_pars_pen << endl;

  of1 << "# Implicit fm_level regression penalty contribution" << endl  
      << " " << ppstf->impl_fm_level_regr_pen << endl;

  of1 << "# von Bertalanffy deviates penalty contribution" << endl  
      << " " << ppstf->vonb_devs_pen << endl;

  of1 << "# Seasonal catchability penalty contribution" << endl  
      << " " << ppstf->seasonal_q_pen << endl;

  of1 << "# Tagged fish reporting rate penalty contribution" << endl  
      << " " << ppstf->tag_rep_rate_pen << endl;

  of1 << "# Regional recruitment deviates penalty contribution" << endl  
      << " " << ppstf->region_recr_dev_pen << endl;

  of1 << "# Average regional recruitment distribution penalty contribution" << endl  
      << " " << ppstf->mean_region_recr_pen << endl;

  of1 << "# Initial population recruitment penalty contribution" << endl  
      << " " << ppstf->init_recr_pen << endl;

  of1 << "# Initial population age composition penalty contribution" << endl  
      << " " << ppstf->init_agecomp_curve_pen << endl;

  of1 << "# Temporal recruitment deviates penalty contribution" << endl  
      << " " << ppstf->temporal_recr_dev_pen << endl;

  of1 << "# Normal distribution temporal recruitment deviates penalty contribution" << endl  
      << " " << ppstf->norm_temporal_recr_pen << endl;

  of1 << "# Temporal recruitment deviates trend penalty contribution" << endl  
      << " " << ppstf->recr_trend_pen << endl;

  of1 << "# Temporal recruitment deviates curve penalty contribution" << endl  
      << " " << ppstf->recr_curve_pen << endl;

  of1 << "# Effort_dev_penalty_by_fishery" << endl 
      << ppstf->effort_dev_penalty_by_fishery << endl;

  of1 << "# Mean_length_constraint_penalty" << endl 
      << ppstf->mean_len_constr_pen << endl;

  of1 << "# Incident selectivity curvature penalty" << endl 
      << ppstf->incid_sel_curve_pen << endl;

  of1 << "# Overall mean selectivity penalty" << endl 
      << ppstf->mean_sel_pen << endl;

  of1 << "# Time-block mean selectivity penalty" << endl 
      << ppstf->time_block_sel_mean_pen << endl;

  of1 << "# Selectivity deviations penalty" << endl 
      << ppstf->sel_devs_pen << endl;

  of1 << "# Non-spline selectivity curvature penalty (all fisheries)" << endl 
      << sum(ppstf->nonspline_sel_curve_pen) << endl;

  of1 << "# Spline selectivity bound penalty (all fisheries)" << endl 
      << sum(ppstf->spline_sel_bound_pen) << endl;

  of1 << "# Selectivity form penalty (all fisheries)" << endl 
      << ppstf->sel_form_pen << endl;

  of1 << "# Maximum fishing mortality penalty (all fisheries)" << endl 
      << ppstf->fmort_max_pen << endl;

  of1 << "# Maximum age-lag fishing mortality penalty (all fisheries)" << endl 
      << ppstf->fmort_agelag_max_pen << endl;

  of1 << "# Orthogonal polynomial recruitment penalty by level: " << endl 
      << "#     year :: season :: region :: season-region" << endl 
      << "#     species 1:" << endl 
      << ppstf->orth_poly_penalty_by_level(1) << endl;
  if (pmsd)
  {
    int ns=pmsd->num_species;
    for (int is=2;is<=ns;is++)
    {
      of1 << "#      species " << is << ":" << endl 
      << ppstf->orth_poly_penalty_by_level(is) << endl;
    }
  }

  ivector ff92=column(fish_flags,92);
  ivector ff99=column(fish_flags,99);
  if (sum(ff92))
  {
    if(sum(ff99))
    { 
    of1 << "# Survey_index_like_by_group" << endl 
        << setprecision(12) << ppstf->survey_index_like_by_group << endl;  //NMD_nov15-21
    }
    else
    {
    of1 << "# Survey_index_like_by_fishery" << endl 
        << ppstf->survey_index_like_by_fishery << endl;  //NMD_aug17-21
    }
  }
  of1 << "# catchability_dev_penalty_by_fishery" << endl 
      << ppstf->catchability_dev_penalty_by_fishery << endl;

  of1 << "# catchability_dev_penalty_by_group" << endl 
      << ppstf->catchability_dev_penalty_by_group << endl;

  of1 << "# total length component of likelihood for each fishery" 
      << endl;

  of1 << setscientific() << setprecision(12) 
            << ppstf->lencontrib << endl;

  for (int i=1;i<=num_fisheries;i++)
  {
    of1 << "# length-sample components of likelihood for fishery " 
        << i << endl;
    for (int j=1;j<=num_fish_times(i);j++)
    {
      int ir=realization_region(i,j);
      int ip=realization_period(i,j);
      int fi=realization_incident(i,j);
      if (len_sample_size(ir,ip,fi)>0)
      {
        of1 << " " << setscientific() << setprecision(12) 
            << ppstf->lencontrib_by_realization(ir,ip,fi) << " ";
      }
    }
    of1 << endl;
  }

  of1 << "# total weight component of likelihood for each fishery" 
      << endl;
  of1 << setscientific() << setprecision(12) 
            << ppstf->wghtcontrib << endl;

  for (int i=1;i<=num_fisheries;i++)
  {
    of1 << "# weight-sample components of likelihood for fishery " 
        << i << endl;
    for (int j=1;j<=num_fish_times(i);j++)
    {
      int ir=realization_region(i,j);
      int ip=realization_period(i,j);
      int fi=realization_incident(i,j);
      if (wght_sample_size(ir,ip,fi)>0)
      {
        of1 << " " << setscientific() << setprecision(12) 
            << ppstf->wghtcontrib_by_realization(ir,ip,fi) << " ";
      }
    }
    of1 << endl;
  }

  for (int i=1;i<=num_fisheries;i++)
  {
    of1 << "# total catch components of likelihood for fishery " 
        << i << endl;
    for (int j=1;j<=num_fish_times(i);j++)
    {
      int ir=realization_region(i,j);
      int ip=realization_period(i,j);
      int fi=realization_incident(i,j);
      of1 << " " << setscientific() << setprecision(12) 
            << ppstf->tot_catch_like_by_realization(ir,ip,fi) << " ";
    }
    of1 << endl;
  }

  if (allocated(ppstf->grouped_tag_like))
  {
    of1 << "# Tag likelihood by tag release by fishery groups" << endl;
  
    int mmin=ppstf->grouped_tag_like.indexmin();
    int mmax=ppstf->grouped_tag_like.indexmax();
    ivector group_flags32=column(fish_flags,32);
    int gmax32=Max(group_flags32);
  
    for (int it=mmin;it<=mmax;it++)
    {
      of1 << "# tag release " << it << endl;
      for (int ig=1;ig<=gmax32;ig++)
      {
        if (ppstf->grouped_tag_like(it,ig).indexmin()>0)
        {
          of1 << "# fishery group " << ig << endl;
          of1 << ppstf->grouped_tag_like(it,ig) << endl;
        }
        else
        {
          of1 << "# XXXXX no tags stuff for tag release " << it 
              << " and tag group " << ig << " XXXXXXX" << endl;
        }
      }
    }
  }
  if (allocated(ppstf->grouped_pooled_tag_like))  //NMD_1Ooct2019
  {
    of1 << "# Tag likelihood for pooled tag release by fishery groups" << endl;
  
    ivector group_flags32=column(fish_flags,32);
    int gmax32=Max(group_flags32);
  
    of1 << "# pooled tag release " << endl;
    for (int ig=1;ig<=gmax32;ig++)
    {
      if (ppstf->grouped_pooled_tag_like(ig).indexmin()>0)
      {
        of1 << "# fishery group " << ig << endl;
        of1 << ppstf->grouped_pooled_tag_like(ig) << endl;
      }
      else
      {
        of1 << "# XXXXX no tags stuff for pooled tag release " 
            << " tag group " << ig << " XXXXXXX" << endl;
      }
    }
  }

  if (allocated(ppstf->ungrouped_tag_like))
  {
    of1 << "# Tag likelihood by tag release by fishery" << endl;
  
    int mmin=ppstf->ungrouped_tag_like.indexmin();
    int mmax=ppstf->ungrouped_tag_like.indexmax();
  
    for (int it=mmin;it<=mmax;it++)
    {
      of1 << "# tag release " << it << endl;
      for (int fi=1;fi<=num_fisheries;fi++)
      {
        if (ppstf->ungrouped_tag_like(it,fi).indexmin()>0)
        {
          of1 << "# fishery " << fi << endl;
          of1 << ppstf->ungrouped_tag_like(it,fi) << endl;
        }
        else
        {
          of1 << "# XXXXX no tags stuff for tag release " << it 
              << " and fishery " << fi << " XXXXXXX" << endl;
        }
      }
    }
  }
  if (allocated(ppstf->ungrouped_pooled_tag_like))  //NMD_1Ooct2019
  {
    of1 << "# Tag likelihood for pooled tag release by fishery" << endl;
  
    of1 << "# pooled tag release " << endl;
    for (int fi=1;fi<=num_fisheries;fi++)
    {
      if (ppstf->ungrouped_pooled_tag_like(fi).indexmin()>0)
      {
        of1 << "# fishery " << fi << endl;
        of1 << ppstf->ungrouped_pooled_tag_like(fi) << endl;
      }
      else
      {
        of1 << "# XXXXX no tags stuff for pooled tag release "
            << " and fishery " << fi << " XXXXXXX" << endl;
      }
    }
  }

  if (allocated(ppstf->age_length_like))
  {
    of1 << "# age length likelihood " << endl;
    of1 << ppstf->age_length_like << endl;
  }
}
