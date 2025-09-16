/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#ifndef _VARIABLE_HPP
  #define _VARIABLE_HPP
//#include <cifstrem.h>
//#include <adstring.hpp>
// The class fish_stock_history contains all the relevant data and
// parameters which determine the exploitation of the stock
// the data are ordered by fishing period. A fishing period is
// a time (month and week of month for length frequency purposes)
// during which at least one fishing incident occurs. A fishing
// incident is an instance of a fishery.  A historical fishery
// consists of a time series of fishing incidents. Fishing incdients
// are grouped into fisheries because fisheries are supposed to reflect
// regularities in the fishing process which means that different
// fishing incidents within a fishery share parameters -- the idea is
// to reduce the number of free parameters to be estimated to a minimum.

// class par_uostream;
// class fishery_catch_at_age_record_array;

  class print_double
  {
    MY_DOUBLE_TYPE x;
    //print_double( const MY_DOUBLE_TYPE _x) { x=_x;}
  public:
    print_double(void) { x=0.0;}
    void operator = (MY_DOUBLE_TYPE _x) {x=_x;}
    print_double pprvalue (const prevariable & x);
    MY_DOUBLE_TYPE to_double(void){ return x;}
  };


  print_double prvalue (const prevariable & x);
  ostream & operator << (const ostream& _os, const print_double& _v);

  class mf_pvm_manager
  {
  public:
    pvmhostinfo *ad_hostp;
    mf_timer  mft;
    int ad_nhost;
    ivector ad_stid;
    int pvm_switch;
    int pvm_save_switch;
    int minspawn;
    int maxspawn;
    int narch;
    ivector hostmatch;
    ivector multiplier;
    mf_pvm_manager(void);
    void setup_pvm(void);
  };

 /*
  class mfcl_thread_mamanger
  {
  public:
    ad_master_slave_synchro * ppsynch;
    mfcl_thread_mamanger(void)
    {
      ppsynch=new ad_master_slave_synchro();
    }
  };
 */

  class plotstuff
  {
  public:
    dmatrix ungrouped_pooled_tag_like;
    dmatrix grouped_pooled_tag_like;
    dvector age_length_like;
    d3_array grouped_tag_like;
    d3_array ungrouped_tag_like;
    d3_array tot_catch_like_by_realization;
    d3_array tag_like;
    dmatrix orth_poly_penalty_by_level;  //NMD_05dec2022
    dvector lencontrib;
    dvector wghtcontrib;
    dvector effort_dev_penalty_by_fishery; 
    dvector survey_index_like_by_fishery;   //NMD_aug17-21
    dvector survey_index_like_by_group;   //NMD_nov15-21
    d3_array lencontrib_by_realization;
    d3_array wghtcontrib_by_realization;
    dvector catchability_dev_penalty_by_fishery;
    dvector catchability_dev_penalty_by_group;
    dvector bh_steep_contribution;  //NMD_jun-27-17
    dvector spline_sel_bound_pen;
    dvector nonspline_sel_curve_pen;
    dvariable reset_pen;  //NMD_21nov2023
    dvariable do_every_pen;
    dvariable reg_recr_pen;
    dvariable diff_coff_pen;
    dvariable vonb_devs_pen;
    dvariable seasonal_q_pen;
    dvariable tag_rep_rate_pen;
    dvariable region_recr_dev_pen;
    dvariable mean_region_recr_pen;
    dvariable init_recr_pen;
    dvariable init_agecomp_curve_pen;
    dvariable temporal_recr_dev_pen;
    dvariable norm_temporal_recr_pen;
    dvariable recr_trend_pen;
    dvariable recr_curve_pen;
    dvariable incid_sel_curve_pen;
    dvariable mean_sel_pen;
    dvariable time_block_sel_mean_pen;
    dvariable sel_devs_pen;
    dvariable sel_form_pen;
    dvariable kludge_surv_pen;
    dvariable impl_fm_level_regr_pen;
    dvariable impl_fm_level_regr_pars_pen;
    dvariable fmort_max_pen;
    dvariable fmort_agelag_max_pen;
    dvariable mean_len_constr_pen;

    //plotstuff(int num_fisheries,int num_regions,const ivector& num_fish_times,
    //  const ivector &  num_fish_periods,
    //  const imatrix& num_fish_incidents,
    //  const imatrix& gfish_index,
    //  const imatrix& fish_flags);
    plotstuff(int num_fisheries,int num_regions,int nspp,
      const ivector& num_fish_times,
      const ivector &  num_fish_periods,
      const imatrix& num_fish_incidents,
      const imatrix& gfish_index,
      const imatrix& fish_flags);

  };
  class xxxtags;

  class smooth_minimizer
  {
    MY_DOUBLE_TYPE eps;
    dvector scoffs;
  public:
    MY_DOUBLE_TYPE poly(MY_DOUBLE_TYPE& u)
    {
      return scoffs(1)+u*(scoffs(2)+u*(scoffs(3)+u*scoffs(4)));
    }
    smooth_minimizer(MY_DOUBLE_TYPE _eps=0.05L) : eps(_eps), scoffs(1,4)
    {
      MY_DOUBLE_TYPE eps2=eps*eps;
      MY_DOUBLE_TYPE eps3=eps*eps2;
      dmatrix A(1,4,1,4);
      A(1,1)=1.0;
      A(2,1)=1.0;
      A(3,1)=0.0;
      A(4,1)=0.0;
    
      A(1,2)=-eps;
      A(2,2)=eps;
      A(3,2)=1.0;
      A(4,2)=1.0;
    
      A(1,3)=eps2;
      A(2,3)=eps2;
      A(3,3)=-2.0*eps;
      A(4,3)=2.0*eps;
    
      A(1,4)=-eps3;
      A(2,4)=eps3;
      A(3,4)=3.0*eps2;
      A(4,4)=3.0*eps2;
    
      dvector v(1,4);
      v(1)=0;
      v(2)=1;
      v(3)=0;
      v(4)=0;
    
      scoffs=solve(A,v);
    
    }
    MY_DOUBLE_TYPE operator () (MY_DOUBLE_TYPE x,MY_DOUBLE_TYPE y)
    {
      MY_DOUBLE_TYPE diff=x-y;
      if (diff>eps)
        return y;
      else if (diff<-eps)
        return x;
  
      //cout << "A " << diff << " " << poly(diff) << endl;
      return x + poly(diff)*(y-x);
    }
  };

  class vsmooth_minimizer
  {
    MY_DOUBLE_TYPE eps;
    dvector scoffs;
  public:
    dvariable poly(const prevariable& u)
    {
      return scoffs(1)+u*(scoffs(2)+u*(scoffs(3)+u*scoffs(4)));
    }
    vsmooth_minimizer(MY_DOUBLE_TYPE _eps=0.05L) : eps(_eps), scoffs(1,4)
    {
      MY_DOUBLE_TYPE eps2=eps*eps;
      MY_DOUBLE_TYPE eps3=eps*eps2;
      dmatrix A(1,4,1,4);
      A(1,1)=1.0;
      A(2,1)=1.0;
      A(3,1)=0.0;
      A(4,1)=0.0;
    
      A(1,2)=-eps;
      A(2,2)=eps;
      A(3,2)=1.0;
      A(4,2)=1.0;
    
      A(1,3)=eps2;
      A(2,3)=eps2;
      A(3,3)=-2.0*eps;
      A(4,3)=2.0*eps;
    
      A(1,4)=-eps3;
      A(2,4)=eps3;
      A(3,4)=3.0*eps2;
      A(4,4)=3.0*eps2;
    
      dvector v(1,4);
      v(1)=0;
      v(2)=1;
      v(3)=0;
      v(4)=0;
    
      scoffs=solve(A,v);
    
    }
    dvariable operator () (const prevariable& x,const prevariable& y)
    {
      dvariable diff=x-y;
      if (diff>eps)
        return y;
      else if (diff<-eps)
        return x;
  
      //cout << "A " << diff << " " << poly(diff) << endl;
      return x + poly(diff)*(y-x);
    }
  };

  class smooth_maximizer
  {
    MY_DOUBLE_TYPE eps;
    dvector scoffs;
  public:
    MY_DOUBLE_TYPE poly(MY_DOUBLE_TYPE& u)
    {
      return scoffs(1)+u*(scoffs(2)+u*(scoffs(3)+u*scoffs(4)));
    }
    smooth_maximizer(MY_DOUBLE_TYPE _eps=0.05L) : eps(_eps), scoffs(1,4)
    {
      MY_DOUBLE_TYPE eps2=eps*eps;
      MY_DOUBLE_TYPE eps3=eps*eps2;
      dmatrix A(1,4,1,4);
      A(1,1)=1.0;
      A(2,1)=1.0;
      A(3,1)=0.0;
      A(4,1)=0.0;
    
      A(1,2)=-eps;
      A(2,2)=eps;
      A(3,2)=1.0;
      A(4,2)=1.0;
    
      A(1,3)=eps2;
      A(2,3)=eps2;
      A(3,3)=-2.0*eps;
      A(4,3)=2.0*eps;
    
      A(1,4)=-eps3;
      A(2,4)=eps3;
      A(3,4)=3.0*eps2;
      A(4,4)=3.0*eps2;
    
      dvector v(1,4);
      v(1)=0;
      v(2)=1;
      v(3)=0;
      v(4)=0;
    
      scoffs=solve(A,v);
    
    }
    MY_DOUBLE_TYPE operator () (MY_DOUBLE_TYPE x,MY_DOUBLE_TYPE y)
    {
      MY_DOUBLE_TYPE diff=x-y;
      if (diff>eps)
        return x;
      else if (diff<-eps)
        return y;
      return y + poly(diff)*(x-y);
    }
  };

  class vsmooth_maximizer
  {
    MY_DOUBLE_TYPE eps;
    dvector scoffs;
  public:
    dvariable poly(const prevariable& u)
    {
      return scoffs(1)+u*(scoffs(2)+u*(scoffs(3)+u*scoffs(4)));
    }
    vsmooth_maximizer(MY_DOUBLE_TYPE _eps=0.05L) : eps(_eps), scoffs(1,4)
    {
      MY_DOUBLE_TYPE eps2=eps*eps;
      MY_DOUBLE_TYPE eps3=eps*eps2;
      dmatrix A(1,4,1,4);
      A(1,1)=1.0;
      A(2,1)=1.0;
      A(3,1)=0.0;
      A(4,1)=0.0;
    
      A(1,2)=-eps;
      A(2,2)=eps;
      A(3,2)=1.0;
      A(4,2)=1.0;
    
      A(1,3)=eps2;
      A(2,3)=eps2;
      A(3,3)=-2.0*eps;
      A(4,3)=2.0*eps;
    
      A(1,4)=-eps3;
      A(2,4)=eps3;
      A(3,4)=3.0*eps2;
      A(4,4)=3.0*eps2;
    
      dvector v(1,4);
      v(1)=0;
      v(2)=1;
      v(3)=0;
      v(4)=0;
    
      scoffs=solve(A,v);
    
    }
    dvariable operator () (const prevariable& x,const prevariable& y)
    {
      dvariable diff=x-y;
      if (diff>eps)
        return x;
      else if (diff<-eps)
        return y;
      return y + poly(diff)*(x-y);
    }
  };

 
  class group_manager_1;
  class  blocked_orthogonal_design;

  class dvar_fish_stock_history
  {
  public:
    xgroup_manager tag_mortality_group;
    dvar_vector tagmort;
    ivector first_survey_time;
    int ff92sum;
    dmatrix survey_cpue_pred; // predicted survey fishery indices NMD_8mar2022
    dmatrix survey_cpue_obs;  // observed survey fishery indices NMD_8mar2022
    dmatrix survey_index;    // for fitting survey fisheries by number
    dvar3_array vul_at_age;    // for fitting survey fisheries by number
    dvariable newfpen;
    blocked_orthogonal_design * pBOD;
    int initial_orthp_estimate_flag;
    dvar_matrix Ninit_standard;
    dvar_matrix Ninit_orthogonal;
    dvar_matrix Npred;
    i3_array movement_ptr;
    dmatrix psqrt_hess;
    dmatrix psqrt_hess_inv;
    dvector pos_hess_scale;
    dvariable kludged_selmean_square;
    dmatrix Xchk;
    d4_array csurv;
    dvar4_array surv;
    d5_array csurv_chk;
    int iloop;
    int loop_flag;
    dvar4_array mean_weight;
    dvariable debug_par;
    MY_DOUBLE_TYPE pamin1;
    MY_DOUBLE_TYPE pamin2;
    int Zmax_flag;
    MY_DOUBLE_TYPE Zmax_fish;
    imatrix implicit_fml_bounds;
    group_flag_manager * pfml_group;

    d4_array fml_designvec;

    d3_array fml_R;
    d3_array fml_Q;
    d3_array fml_M;
    imatrix fml_columns;

    dmatrix ptr_fml_Mbs_mnth_finyr;  //NMD_4jul2022
    d3_array fml_Mbs_finyr;  //NMD_4jul2022
    ivector fshtms_finyr;  //NMD_4jul2022
    ivector eff_proj_fshry;  //NMD_11jul2022
    ivector catch_proj_fshry;
    dmatrix q_level_finyr; //NMD_7may2024
    
    int print_implicit_effort_flag;
    int implicit_flag;
    int censored_gamma_report_flag;
    dvector censor_report;
    int nb_tag_report_flag;
    dmatrix orthogonal_diffusion_matrix;
    dmatrix new_orthogonal_diffusion_matrix;
    int first_unfixed_year;
    ivector first_unfixed_fish_time;
    d3_array allprob;
    d3_array alltagsurv;
    d3_array prob;
    imatrix sel_pointer;
    imatrix months_used;
    dvector eta_hat;
    dvariable testpar;
    ccubic_spline_array * pccsa;
    cubic_spline_array ** pcsa;
    cubic_spline_array ** wpcsa;
    //au_au_cubic_spline_stuff * pau_au_css;
    int no_lagrangian_update;
    ~dvar_fish_stock_history();
    vsmooth_maximizer smax;
    vsmooth_minimizer smin;
    d3_array len_sample_size;
    d3_array wght_sample_size;
    MY_DOUBLE_TYPE bl_sel_scaling; // scaling factor for selectivity pars
    int get_annual_recruitment_flag;
    ivector projected_simulated_data_flags;
    ivector simulation_seeds;
    dvector average_effort;
    
    const int max_group_ptr;
    imatrix group_ptr;   // group each grouped obect is in
    imatrix inv_group_ptr;   // first grouped object in each group
    dmatrix lSinv;
    d4_array len_Sinv;
    dmatrix wtSinv;
    d4_array wt_Sinv;
    imatrix better_sbb;
    imatrix yearblock;
    imatrix selblks_ff16;
    imatrix selseas_ff16;
    dvar4_array bstempsel;   // selectivity of an age class j fish or length
    dvar4_array bswtempsel;  // selectivity of an age class j fish or weight
    dvar4_array bstempsel_afl; // selectivity of an age class j fish from length
    i3_array sseason;
    i3_array bblock;
    dvar4_array kludged_surv; 
    dvar4_array bs_selcoff;   // probability that an age class j fish
    dvar3_array bs_selmean;   // probability that an age class j fish
    ivector num_blocked_seasons;
    ivector num_blocks;
    ivector num_breaks;
    i3_array sel_seasons;
    int na;
    dvariable annual_phi;
    dvariable log_length_variance;
    dvariable log_length_dof;
    dvariable length_rho;
    dvariable length_exp;
    dvariable length_tot_exp;
    dvariable length_psi;
    dvariable log_weight_variance;
    dvariable log_weight_dof;
    dvariable weight_rho;
    dvariable weight_exp;
    dvariable weight_tot_exp;
    dvariable weight_psi;
    dvar3_array blbsel;
    int sum_ff71_flag;
    int sum_ff48_flag;
    int sum_ff79_flag;
    ivector block_ptr;
    ivector blocked_fishery_index;
    ivector sel_block_index;
    ivector sel_block_fisheries;
    imatrix sel_block_breaks; 
    imatrix break_block; 
    int threaded_tag_flag;
    int set_global_vars_flag;
    pmulti_species_data  pmsd;
    ofstream * tagofs;
    ofstream * proj_output_files[10];
    plotstuff * ppstf;
    ivector q_flag;
    int tag_fish_switch;
    int ff263flag;
    int  projection_sim_index;
    int  ny_begin_yr;
    int  ny_begin_reg;
    int  ny_begin_ses;
    int  ny_begin_ses_reg;
    int  nd_begin_yr;
    int  nd_begin_reg;
    int  nd_begin_ses;
    int  nd_begin_ses_reg;
    int  ny_end_yr;
    int  ny_end_reg;
    int  ny_end_ses;
    int  ny_end_ses_reg;
    int  nd_end_yr;
    int  nd_end_reg;
    int  nd_end_ses;
    int  nd_end_ses_reg;
    dmatrix OR;
    int recr_degree_yr;
    int recr_degree_reg;
    int recr_degree_ses;
    int num_new_weights;
    dvar_vector new_orth_recr;
    dmatrix recr_polys_yr;
    dmatrix recr_polys_reg;
    dmatrix recr_polys_ses;
    dmatrix recr_polys_ses_reg;
    d3_array simulated_numbers_at_age;
    d3_array simulated_numbers_at_age_noeff; //NMD 22Feb2012
    d3_array simulated_recruitments;
    d3_array simulated_effort_devs;
    dvar_vector rec_delta;

    //mfcl_thread_mamanger * pmfcltm;
    fmmt1 * pfmin1;
    dvariable imppen;
    int ss2_flag;
    d3_array orth_recr_basis;
    dmatrix grouped_catchability_coffs;
    ivector numcomp;
    int degree_yr;
    int degree_reg;
    int degree_ses;
    int degree_ses_reg;
    ivector ses_reg_rsum;
    ivector ses_reg_csum;
    imatrix ses_reg_ractives;
    imatrix ses_reg_cactives;
    imatrix ses_reg_active;
    int num_ses_reg_active;
    int have_projection_periods_flag;
    int first_data_month;
    dvar_vector sv;
    imatrix year_flags;
    int num_seasons;
    int seasonal_recruitment_flag;
    imatrix season_flags;
    dvariable bh_variance;
    dvariable save_bh_variance;			//NMD20Jan2012
    dvar3_array mean_weight_yr;
    dvar3_array mean_weight_yr_alternative;
    dvar_matrix mean_weight_yr_proj;
    int af170q0;
    int af170q0ex;
    int zeroed_fisheries_flag;
    dvariable alpha;
    dvariable beta;
    dvariable steepness;
    dvariable phi;
    MY_DOUBLE_TYPE old_log_lambda;
    int total_num_obs;
    MY_DOUBLE_TYPE likecomponent[10];
    movement_info mo; 
    int direction_flag;
    int generate_report;
    int age_age1;
    int age_nage;
    int     month_1;       // the month in which recruitment occurs
    int frq_file_version;
    ivector      parest_flags;        // the old control flags from the standard multifan model so they can be accessed from the class members// The standard flags from the old multifan model
    ivector      historical_parest_flags;        // the old control flags from the standard multifan model so they can be accessed from the class members// The standard flags from the old multifan model
    ivector      age_flags;           // The new flags for the age structured part of the model
    ivector      historical_age_flags;           // The new flags for the age structured part of the model
    imatrix      fish_flags;          // Fishery specific control flags
    imatrix      historical_fish_flags;   // Fishery specific control flags
    imatrix      data_fish_flags;          // Fishery specific control flags
    imatrix      old_fish_flags;          // Fishery specific control flags
    imatrix      old_zero_effdev_flag;
    imatrix      tag_flags;           // 
    imatrix      true_tag_flags;           // 
    dvar_vector  tag_return_probability;
    imatrix      move_flags;           // 
    ivector      move_map;           // 
    // fish_flags(i,1) :  the number of age classes with selecitvity in fisehry i
    int          nage;                // The number of age classes
    dvar4_array  nrsurv;
    dvar_vector  biomass;
    int min_tag_group;
    int max_tag_group;
    dvar_vector  adult_biomass;
    dvar_vector  bh_recr_devs;
    dvector  bh_predicted_recruits;
    dvector  randrec;
    dvector  bh_numbers;
    dvector bh_reproductive_biomass;
    dvector  bh_bio;

    i3_array grouped_fishery_projection_flag;
    i3_array fishery_projection_flag;
    int  do_fishery_projections_flag;

    imatrix ses_reg_recr_flags;
    dvar_matrix  biomass_by_region;
    dvar_matrix  rel_biomass_by_region;
    dvar_vector  catch_biomass;
    dvar_matrix  catch_biomass_by_fishery;
    dvar_vector  selmean;
    dvar_vector  pop_delta;
    dvar_vector  epop_delta;
    dvar_vector  ageselmean;
    dmatrix  effort_weight_by_fishery;
    dvar_vector  catch_numbers;       // JH 27/03/02
    dvar3_array  F_by_age_by_year_by_region;
    dvar_matrix  F_by_age_by_year;
    dvar5_array  nrfm;
    dvar4_array  nrtm;
 // fisheries grouping stuff used in catchability calculations
    imatrix gfish_ptr;
    imatrix real_gfish_ptr;
    ivector num_grouped_fish_times;
    ivector num_real_grouped_fish_times;
    ivector gp_year;
    ivector gp_month;
    ivector projection_year;
    ivector projection_month;
    int ngroups;
    imatrix global_fishing_periods;  // need for grouping effort for cobbs-douglas production function
    ivector effective_len_size;
    ivector effective_weight_size;
    dmatrix grouped_effort;
    dmatrix grouped_fish_time;
    imatrix grouped_year;
    imatrix grouped_month;
    imatrix grouped_week;
    imatrix grouped_true_year;
    imatrix grouped_true_month;
    imatrix grouped_true_week;
    dmatrix grouped_between_times;
    dvar_matrix grouped_catchability;
    imatrix gfish_index;
    dvar_matrix grouped_catch_dev_coffs;

    int          nyears;              // The number of years over which fishing occurred
    int          last_real_year;              // The number of years over which fishing occurred
    int          num_real_years;              // The number of years over which fishing occurred before month "doubling"
    int          num_tag_releases; // the number of sets of releases of tagged fish
    int          sim_num_tag_releases; // the number of sets of releases of tagged fish
    imatrix      sim_num_tags_at_length; 
    ivector          sim_tag_fishery; // the number of sets of releases of tagged fish
    int          old_num_tag_releases; // the number of sets of releases of tagged fish
    int          true_num_tag_releases; // the number of sets of releases of tagged fish
    int          min_tag_year;
    int          sim_min_tag_year;
    ivector      sim_tag_region;        // the region in which each tag set was released
    ivector      sim_tag_numbers_released;        // the region in which each tag set was released
    ivector      tag_region;        // the region in which each tag set was released
    ivector      tag_year;          // the year in which each tag set was released
    ivector      sim_tag_year;          // the year in which each tag set was released
    ivector      sim_true_tag_year;          // the year in which each tag set was released
    ivector      true_tag_year;     // the year in which each tag set was released
    ivector      tag_month;         // the month in which each tag set was released
    ivector      sim_tag_month;         // the month in which each tag set was released
    imatrix      sim_initial_tag_period;         // the month in which each tag set was released
    ivector      sim_tag_incident;         // the month in which each tag set was released
    ivector      true_tag_month;    // the month in which each tag set was released
    ivector      sim_true_tag_month;    // the month in which each tag set was released
    ivector      itr;               // the number of differenct periods in which tags were returned fomr that tag group
    ivector      itind;            // the first realization with tags present
    ivector      initial_tag_year;   // the  each tag set was released
    ivector      sim_initial_tag_year;   // the  each tag set was released
    ivector      terminal_tag_year;   // the  each tag set was released
    ivector      sim_terminal_tag_year;   // the  each tag set was released
    imatrix      initial_tag_period;   // the  each tag set was released
    imatrix      terminal_tag_period;   // the  each tag set was released
    imatrix      sim_terminal_tag_period;   // the  each tag set was released
    imatrix      initial_tag_recruitment_period;   // the  each tag set was released
    imatrix      sim_initial_tag_recruitment_period;   // the  each tag set was released
    dmatrix      initial_tag_release_by_length; // the tag release data by length intervals
    dvar_matrix      initial_tag_release_by_age; // the tag release data by length intervals
    dvar_matrix      sim_initial_tag_release_by_age; // the tag release data by age
    dmatrix      csim_initial_tag_release_by_age; // the tag release data by age
    i3_array     tag_recaptures_by_length;
    //dvar_vector region_rec_diff_colmeans;
    i4_array min_tag_age2;
    i4_array min_tag_age4;
    i4_array sim_min_tag_age4;
    i3_array min_tag_age1;
    i3_array sim_min_tag_age1;
    i3_array min_tag_age3;
    i3_array min_tag_age5;
    i3_array sim_min_tag_age5;
    i3_array min_tag_age6;
    i3_array sim_min_tag_age6;
    i3_array min_tag_age;
    i3_array sim_min_tag_age;
    int          tag_shlen;       // The number of fisheries
    MY_DOUBLE_TYPE       tag_filen;       // The number of fisheries
    int          tag_nlint;       // The number of fisheries

    int          first_time;       // used as date offset for renumbering data dates
    int          month_factor;       // used as date offset for renumbering data dates
    int          num_fisheries;       // The number of fisheries
    ivector      minttp;
    ivector      sim_minttp;
    ivector      offset;
    ivector      min_init_tag_period;
    ivector      sim_min_init_tag_period;
    ivector      maxttp;
    ivector      minimum_initial_tag_period;
    ivector      sim_minimum_initial_tag_period;
    dvar_vector  pmature;
    dvector      calc_pmature_at_age;
    int          num_recruitment_periods;
    int          initial_recruitment_count;
    ivector      num_fish_periods;    // The number of fishing periods
    imatrix      tag_num_fish_periods;    // The number of fishing periods
    ivector      num_real_fish_periods;    // The number of fishing periods
    int          num_fish_data_recs;  // The total number of fishery data records
    int          num_regions;  // The total number of fishery data records
    int          num_real_regions;  // The total number of fishery data records
    ivector      fishery_regions;  // The region in which each fishery occurs
    imatrix      num_fish_incidents;  // The number of fishing incidents which occurred during each fishing period in each region
    imatrix      recruitment_period;  // 
    i3_array     num_tagfish_incidents;  // The number of fishing incidents which occurred during each fishing period in each region for each tag group
    i3_array     sim_num_tagfish_incidents;  // The number of fishing incidents which occurred during each fishing period in each region for each tag group
    i3_array     num_alltagfish_incidents;  // The number of fishing incidents which occurred during each fishing period in each region for each tag group
    i3_array     sim_num_alltagfish_incidents;  // The number of fishing incidents which occurred during each fishing period in each region for each tag group
    imatrix      num_pooledtagfish_incidents;  // The number of fishing incidents which occurred during each fishing period in each region
    imatrix      sim_num_pooledtagfish_incidents;  // The number of fishing incidents which occurred during each fishing period in each region
    ivector      num_fish_times;      // the number of times a particular fishery occurred
    ivector      num_real_fish_times;      // the number of real times a particular fishery occurred ie not including projections
    ivector      fishing_period;      // the fishing period during which an instance of a fishery occurs
    dvar_vector totalcatch_by_numbers;
    dvector obstotalcatch_by_numbers;
    ivector numtotalcatch_by_numbers;
    dvar_vector totalcatch_by_weight;
    dvector obstotalcatch_by_weight;
    ivector numtotalcatch_by_weight;
    dmatrix      effort_by_fishery;
    dmatrix      true_effort_by_fishery;
    d3_array     really_true_effort;
    dmatrix      log_true_effort_by_fishery;
    int missing_catch_flag;
    ivector missing_catch_counter;
    ivector missing_catch_by_fishery_flag;
    i3_array missing_effort_by_region_flag;
    imatrix num_missing_effort_by_region;
    imatrix num_present_effort_by_region;
    ivector 	missing_effort_flag;
    imatrix 	missing_effort_by_realization_flag;
    dmatrix      log_effort_by_fishery;
    dmatrix      normalized_log_effort_by_fishery;
    dmatrix      log_pred_effort_by_fishery;  //NMD_20may2021
    ivector      fishing_region;      // the fishing period during which an instance of a fishery occurs
    ivector      fishing_incident;    // the fishing_incident corresponding to an instance of a fishery
    i3_array     parent;              // Point to the fishery containing a fishing incident
    i3_array     fish_times;              // the order of this incident by fishery
    imatrix      year;                // The year during which the fishing incident occurred
    imatrix      movement_period;     // The year during which the fishing incident occurred
    int          year1;                // The first year a fishing incident occurred
    imatrix      realization_period;  // The fishing period in which each realization of a fishery occurred
    i3_array      fishery_realization_index;  // The fish_time for the occurence of this fishery 
    imatrix      realization_incident;// The fishing incident to which each realization of a fishery corresponded
    imatrix      realization_region;// The fishing incident to which each realization of a fishery corresponded
    i3_array      header_record_index; // The index of the (sorted) header record which corresponds to a particular fishing period and fishing incident
    dvar_matrix region_pars;
    dvar_matrix current_biomass_by_year;
    imatrix region_flags;
    dvar_vector  gml;
    dvar_vector predicted_yield_bh;
	dvar_vector predicted_eqbio_bh;    // JH 21/03/02 - equil. adult biomass
	dvar_vector predicted_eqtotbio_bh; // JH 21/03/02 - equil. total biomass
    dvar_vector tb_ratio;              // JH 27/03/02 - total biomass over TB at MSY
    dvar_vector ab_ratio;              // JH 27/03/02 - adult biomass over AB at MSY
    dvar_vector F_ratio;               // JH 27/03/02 - aggregate F over F at MSY
    dvector predicted_yield_bh_x;
    dvar_vector predicted_yield_pt;
    dvector predicted_yield_pt_x;
    dvar_vector predicted_recruitment_bh;
    dvector predicted_recruitment_bh_x;
    dvector  region_area;
    dvector  thread_xsave;
    MY_DOUBLE_TYPE thread_f;
    dvar_matrix region_rec_diffs;
    dvar_matrix region_rec_diff_coffs;
    dvar_vector region_rec_diff_sums;
    dvar_matrix      fraction;            // The fraction of the total natural mortality occurring during this fishing period
    ivector      regmin;
    ivector      regmax;
    dvariable        totpop;              // population size scaling parameter
    dvariable        totpop_coff;              // population size scaling parameter
    dvariable        implicit_totpop_coff;              // population size scaling parameter
    dvariable        recmean; 
    dvariable        initmean; 
    dvariable    rec_init_diff;           // diff level between rec and init pop
    dvar_matrix      D;                   // moves the fish around
    dvar4_array      Dad;                 // moves the fish around
    dvar4_array      Dad2;                 // moves the fish around
    imatrix      Dflags;                  // moves the fish around
    //dvar_matrix      diff_coff_species_region_means;   // parameters in the diffusion matrix
    //dvar_vector      diff_coff_species_region_means_sum;   // parameters in the diffusion matrix
    dvar_matrix      diff_coffs;          // parameters in the diffusion matrix
    dvar_matrix      xdiff_coffs;          // parameters in the diffusion matrix
    dvar_matrix      y1diff_coffs;    // sum parameters in the diffusion matrix
    dvar_matrix      y2diff_coffs;    // sum proportion parameters in the diffusion matrix
    dvar_matrix      zdiff_coffs;          // parameters in the diffusion matrix
    dvar_matrix      diff_coffs2;         // parameters in the diffusion matrix
    dvar_matrix      diff_coffs3;         // parameters in the diffusion matrix
    dvar_matrix      diff_coffs_prior;          // parameters in the diffusion matrix
    dvar_matrix      diff_coffs2_prior;         // parameters in the diffusion matrix
    dvar_matrix      diff_coffs3_prior;         // parameters in the diffusion matrix
    //dvar_vector      xdiff_coffs;       // parameters in the diffusion matrix
    //dvar_vector      xdiff_coffs2;      // parameters in the diffusion matrix
    //dvar_vector      xdiff_coffs3;      // parameters in the diffusion matrix
    dvar_vector      recr;                // the relative recruitment levels
    ivector        rec_times;                // the relative recruitment levels
    dvector        rec_covars;                // environmenital factors affecting recruitment 
    dvar_vector      tmprecr;                // the relative recruitment levels
    dvar_vector      avail_coff;         // orthogonal polys which  determine
    dvar_matrix      orth_recr;         // orthogonal polys which  determine
					// the relative recruitment levels
    dvar_matrix      orth_recr_all;    // orthogonal poly coffs which  determine
    dvar_matrix      yearly_recr_all;    // orthogonal poly coffs which  determine
    dvar_vector      orth_recr_tot;       // orthogonal poly coffs which  determine
					// overall recruitment levels
    dvar_matrix      orth_recr_season;    // orthogonal poly coffs which  determine
					// seasonal component of recruitment levels
    dvar_matrix      orth_recr_region;    // orthogonal poly coffs which  determine
					// regional component of recruitment levels
    dvar3_array      orth_recr_ses_reg;    // orthogonal poly coffs which  determine
					// regional component of recruitment levels
    dmatrix         recr_polys;
    dvar_matrix      initpop;             // the relative inital population at age levels
    dvar_vector      tmpinitpop;             // the relative inital population at age levels
    dvar_vector      actual_recruit;           // The total number of fish recruiting to the population each year
    dvar_vector      actual_init;           // The initial number of fish in the population in year 1
    dvar3_array      num_fish;            // The number of fish in each age class in the population at the beginning of each fishing period
    dvar3_array      num_fish_q0;            // The number of fish in each age class in the population at the beginning of each fishing period
    dvar3_array      num_fish0;            // The number of fish in each age class in the population at the beginning of each fishing period in the absence of fishing
    dvar3_array      total_num_fish;            // The number of fish in in the population at the beginning of each fishing period
    dvar4_array      tagnum_fish;            // The number of fish in each tag group in age class in the population at the beginning of each fishing period
    dvar4_array      sim_tagnum_fish;            // The number of fish in each tag group in age class in the population at the beginning of each fishing period
    dvar3_array      pooled_tagnum_fish;            // The number of fish in each tag group in age class in the population at the beginning of each fishing period
    dvar3_array      sim_pooled_tagnum_fish;            // The number of fish in each tag group in age class in the population at the beginning of each fishing period
    dvar3_array      epooled_tagnum_fish_recr;            // The number of fish in each tag group in age class in the population at the beginning of each fishing period
    dvar3_array      sim_epooled_tagnum_fish_recr;            // The number of fish in each tag group in age class in the population at the beginning of each fishing period
    dmatrix       tag_release_by_length;
    dvar_matrix   lbsel;
    dvar3_array      pooledtagN;            // The number of fish in each age class in the population at the beginning of each Year
    dvar3_array      sim_pooledtagN;            // The number of fish in each age class in the population at the beginning of each Year
    dvar4_array      tagN;            // The number of fish in each age class in the population at the beginning of each Year
    dvar4_array      sim_tagN;            // The number of fish in each age class in the population at the beginning of each Year
    dvar3_array      tagrelease;            // The number of fish in each age class in the population at the beginning of each Year
    dvar3_array      N;            // The number of fish in each age class in the population at the beginning of each Year
    dvar3_array      Nsave;            // The number of fish in each age class in the population at the beginning of each Year
    dvar_matrix      Rsave;            // The number of fish in each age class in the population at the beginning of each Year
    dvar3_array      N_q0;            // The number of fish in each age class in the population at the beginning of each Year
    //dvar3_array      N0;            // The number of fish in each age class in the population at the beginning of each Year in the absence of fishing
    dvar3_array      exp_N;            // The number of fish in each age class in the population at the beginning of each Year
    dvar4_array     catch;               // The number of fish in the catch from each age class, fishing period, and each fishery
    dvar4_array     catch_q0;               // The number of fish in the catch from each age class, fishing period, and each fishery
    dvar3_array     pred_totcatch;               // The number of fish in the catch from each age class, fishing period, and each fishery
    dvar4_array     exp_catch;               // The number of fish in the catch from each age class, fishing period, and each fishery
    dvar5_array     tagcatch;               // The number of tagged fish in the catch from each age class, fishing period, and each fishery
    d5_array     probtagcatch;               // The number of tagged fish in the catch from each age class, fishing period, and each fishery
    dvar5_array     sim_tagcatch;               // The number of tagged fish in the catch from each age class, fishing period, and each fishery
    dvar5_array     sim_obs_tagcatch;               // The number of tagged fish in the catch from each age class, fishing period, and each fishery
    dvar4_array     pooled_tagcatch;               // The number of tagged fish in the catch from each age class, fishing period, and each fishery
    dvar4_array     sim_pooled_tagcatch;               // The number of tagged fish in the catch from each age class, fishing period, and each fishery
    dvar5_array     obstagcatch;               // The number of tagged fish in the catch from each age class, fishing period, and each fishery
    dvar4_array     pooledobstagcatch;               // The number of tagged fish in the catch from each age class, fishing period, and each fishery
    dvar4_array     sim_pooledobstagcatch;               // The number of tagged fish in the catch from each age class, fishing period, and each fishery
    dvar5_array obstagcatch1;
    d4_array     tot_tag_catch;               // The number of tagged fish in the catch from each age class, fishing period, and each fishery
    d4_array     sim_tot_tag_catch;               // The number of tagged fish in the catch from each age class, fishing period, and each fishery
    d3_array     region_tot_tag_catch;               // The number of tagged fish in the catch from each age class, fishing period, and each fishery
    d3_array     pooledtot_tag_catch;               // The number of tagged fish in the catch from each age class, fishing period, and each fishery
    d3_array     sim_pooledtot_tag_catch;               // The number of tagged fish in the catch from each age class, fishing period, and each fishery
    d5_array     obstagcatch_by_length;     // The number of tagged fish in the catch from each age class, fishing period, and each fishery
    d4_array     pooledobstagcatch_by_length;     // The number of tagged fish in the catch from each age class, fishing period, and each fishery
    dvar3_array      tot_catch;           // The total number of fish in the catch fishing period, and each fishery
    dvar4_array      prop;                // The number of proportion in the catch from each age class, fishing period, and each fishery
    dvar4_array      incident_sel;        // The selectivity parameters for each fishing incident
    dvar5_array      age_len_incident_sel;        // The selectivity parameters for each fishing incident
    dvar5_array      age_weight_incident_sel;        // The selectivity parameters for each fishing incident
  dvar4_array      age_len_incident_sel_yr;    
    dvar3_array      mean_incident_sel;        // The selectivity parameters for each fishing incident
    dvar3_array     delta2;             // The random variations in selectivity 
    dmatrix      between_times;        // The  times between realizations of a fishery
    imatrix      imp_back;        // How far back to go in implicit robust mean calculations for catchability
    imatrix      imp_forward;      // How far forward to go in implicit robust mean calculations for catchability
    dvar_matrix      fishery_sel;         // The selectivity parameters for each fishery
    dvar_matrix      selcoff;             // The selectivity parameters for each fishery
    dvar_matrix      ageselcoff;             // The age dependent part of selectivity parameters for each fishery
    dvar3_array      survival;            // The survival rate for each fishing period
    dvar3_array      survival_q0;            // The survival rate for each fishing period
    dvar4_array      fish_mort;           // The fishing mortality rate for each fishing period, fishing incident and age class
    dvar4_array      fish_mort_q0;           // The fishing mortality rate for each fishing period, fishing incident and age class
    dvar4_array      fish_mort_calcs;           // The fishing mortality rate for each fishing period, fishing incident and age class
    dvar5_array      tag_fish_mort_calcs;           // The fishing mortality rate for each fishing period, fishing incident and age class
    dvar4_array      fish_mort_calcs_q0;           // The fishing mortality rate for each fishing period, fishing incident and age class
    dvar4_array      tag_tot_mort;            // The total mortality rate for each fishing period
    dvar3_array      tot_mort;            // The total mortality rate for each fishing period
    dvar3_array      e_nat_mort_by_period;            // The natural mortality rate for each fishing period
    dvar3_array      tot_mort_q0;            // The total mortality rate for each fishing period
    dvar3_array      tot_mort0;            // The total mortality rate for each fishing period
    dvar_matrix      nat_mort;            // The annual mortality rate for each year and age class
    dvar3_array      catchability;        // The catchability by fishing period by fishing incident
    dvar3_array      catchability_q0;        // The catchability by fishing period by fishing incident
    dvar3_array      FF;        // The catchability by fishing period by fishing incident
    d3_array      effort;              // the fishing effort by fishing period by fishing incident
    d3_array      effort_normalization_factor;              // the fishing effort by fishing period by fishing incident
    dvar3_array      fm_level;              // the fishing effort by fishing period by fishing incident
    d3_array      effort_weight;              // the fishing effort by fishing period by fishing incident
    dvar_matrix		fm_level_devs;  // when there are time-series changes in implcit catchability
    imatrix          month;               // the month of the year during which each fishery occurs
    imatrix          true_month;               // the month of the year during which each fishery occurs
    imatrix          really_true_month;               // the month of the year during which each fishery occurs
    imatrix          true_week;               // the month of the year during which each fishery occurs
    imatrix          true_year;               // the month of the year during which each fishery occurs
    imatrix          really_true_year;               // the month of the year during which each fishery occurs
    imatrix          week;                // The week of the month during which each fishery occurs
    imatrix          move_index;          // index of the parameters used to parameterize movement
    dvariable        nat_mort_coff;       // determines the natural mortality
    MY_DOUBLE_TYPE           catch_init;          // Initial value for the catchability
    dvar_vector      q0;                  // overall catchability coefficient
    dvar_vector      q0_miss;             // overall catchability coefficient
                                          // for catch-conditioned calculations
                                          // with missing catches

    dvar_vector      q1;                  // time dependent trend in catchability coefficient
    dvar_matrix      effort_dev_coffs;
    imatrix          zero_effdev_flag;
    imatrix          zero_catch_for_period_flag; 
    i3_array         zero_catch_for_incident_flag;
    imatrix          missing_catch_for_period_flag; 
    i3_array         missing_catch_for_incident_flag;
    imatrix         missing_catch_for_realization;
    dvar_matrix      catch_dev_coffs;    // determine random walk deviations in catchability time trend
    dvar_matrix      rep_dev_coffs;    // determine random walk deviations in reporting rate time trend
    dvar3_array      rep_rate;          // Tag reporting rate by region, fishing period, fishing incident
    dvar_matrix      implicit_catchability;    //
    dvar_matrix      implicit_catchability_ccond;    //
    dvar_matrix      new_fm_level;    //
    dvar_matrix      fish_pars;    // extra fisheries related parameters
    dvar_matrix      implicit_fm_level_regression_pars;  // extra fisheries related parameters

    
    dvar_matrix      species_pars;    // extra fisheries related parameters
  // neweswtuff for tag * fishery reporting rates
  // ***************************************************************
    dvar_matrix      tag_fish_rep; // tag * fish reporting rates 
    dmatrix      sim_tag_fish_rep; // tag * fish reporting rates 
    imatrix      tag_fish_rep_group_flags; // tag * fish reporting rates 
    imatrix      tag_fish_rep_active_flags; // tag * fish reporting rates 
    dmatrix      tag_fish_rep_target; // tag fish reporting rate target value
    dmatrix      tag_fish_rep_penalty; // tag fish reporting rate target value
                                      // penalty


    dvariable tag_rep_rate(int it,int ir,int rp,int fi); // function to get
  // ***************************************************************
    


    dvar_matrix      seasonal_catchability_pars;    // explicit seasonal catchability as opposed to trig function
    imatrix seasonal_catchability_pars_index;
    dmatrix seasonal_catchability_pars_mix;
    dvar_matrix      age_pars;    // extra age-class related parameters
    dmatrix      biomass_index;    // extra fisheries related parameters
    dvar3_array      effort_devs;         // deviations in effort fishing mortality relationship
    dvar3_array      sel_dev_coffs;       //determine deviations in selectivity
    dvar4_array      sel_devs;       //determine deviations in selectivity
    dvar_vector      corr_wy;             // determines within year correlation in sel deviations
    dvar_vector      corr_by;             // determines between year correlation in sel deviations
    dvar_vector      corr_wc;             // determines within cohort correlation in sel deviations
    dvar_vector      corr_eff;            // determines  // The data
    d3_array         obs_tot_catch;
    dmatrix         obs_region_tot_catch;
    d4_array     obs_catch;               // The number of fish in the catch from each age class, fishing period, and each fishery
    d4_array         obs_prop;            // The estimated proportion in the catch from each age class, fishing period, and each fishery
    dvar_matrix qmu;
    dvar3_array lagrange_c;
    d3_array lagrange_lambda;
    MY_DOUBLE_TYPE lagrange_mu;
    dvar_vector other_lagrange_c;
    dvector other_lagrange_lambda;
    dvector other_lagrange_mu;
    dvar_matrix qvar;
    dmatrix qstudent;
    dvar4_array kludged_equilib_coffs;
    dvar_vector kludged_equilib_level_coffs;
    ivector get_min_time_indices(ivector& rip);
    dvar_vector get_bh_recruitment(dvar3_array& NP,int iy);
    dvariable get_bh_recruitment_multiplier(dvar3_array& NP,int iy);
    void make_fishing_period_report(void);
    int have_fish_periods(ivector& rip);
    void read_movement_coffs(adstring & movement_coffs_file);
    void set_fishery_projection_flags(void);
    void allocate_fishery_projection_flags(void);
    void set_grouped_fishery_projection_flags(void);
    void allocate_grouped_fishery_projection_flags(void);
    void un_normalize_year(void);
    void do_grouped_between_time_calculations(void);
    void set_catchabilities_to_zero(ivector * pq_flag);
    void do_the_diffusion_yf(int year,dvar_matrix& DINV);
    void fishing_mortality_calc(void);
    void fishing_mortality_calc(d3_array&,d3_array&);
    void fishing_mortality_calc(d3_array&,int);
    dvar_fish_stock_history(int ng,int nfp, ivector& nfi, int nfsh, int nyrs);
    dvar_fish_stock_history(int ng,int nfp, ivector& nfi, int nfsh,
      int nyrs,ivector& flags);
    void seasonal_catchability_calc(void);
    void fishery_selectivity_calc(void);
    void incident_selectivity_calc(void);
    dvariable get_initial_population(dvar_vector& sv,int,void *);
    void get_initial_population_test(dvar_vector& sv);
    void get_pop_delta(void);
    void get_initial_tag_population(dvar_vector& sv, int it);
    void do_the_diffusion(int year,dvar_vector& sv,dvar3_array& N,ofstream * pofs=NULL);
    void add_cs_selectivity_coffs(void);
    void add_equilibrium_selectivity_coffs(void);
    void do_the_diffusion(ivector& rip,dvar_vector& sv,dvar3_array& N,ofstream * pofs=NULL);
    void do_the_diffusion(int it,int year,dvar_vector& sv,dvar3_array& N,ofstream * pofs=NULL);
    void get_exp_N(void);
    void proportion_at_age_calc(void);
    void total_fish_mort_and_survival(int it,int ip,int yr,
      const dvar_vector & tmp);
    int fishing_selectivity_interface_nvar(void);
    void fishing_selectivity_interface
      (dvar_vector& x, int& ii,const prevariable& pen);
    void fishing_selectivity_interface_inv
      (dvector& x, int& ii);
    void fishery_selcoff_init();
    void get_fishing_and_total_mortality_for_this_period_q0(int ir,int ip,
      dvar_vector& lcurrent_rbio);
    void get_fishing_and_total_mortality_for_this_period(int ir,int ip,
      dvar_vector& lcurrent_rbio);
    void rep_rate_devs_calc(void);
    void do_everything_calc();
    void do_everything_calc(d3_array&);
    void get_initial_parameter_values();
    dvariable fit_implicit_q(void);
    void transform_in();
    void transform_out();
    //have_data_this_year(void);
    dvariable have_data_this_year(int it,int ir,int& ip,int cy);
    void print_tag_accounting_info(void);
    void calculate_tag_catches(int it);
    void have_no_data_this_year(int it,int ir,int cy);
    void initial_population_profile_calc();
    void natural_mortality_calc(void);
    void natural_mortality_splines(void);
    void natural_mortality_calc2(void);
    //void catchability_init();
    void natural_mortality_init();
    void catchability_calc(void);
    void catchability_devs_calc(void);
    void set_zero_effdev_flag(void);
    void set_missing_totcatch_flags(void);
    void set_zero_catchdev_flag(void);
    void get_observed_total_catch();
    //void catch_equations_calc(dvar_vector& sv);
    void xcatch_equations_calc(dvar_vector& sv);
    void catch_equations_calc_movement(dvar_vector& sv);
    void tag_catch_equations_calc_mc(dvar_vector& sv);
    void catch_equations_calc_avail();
    dvariable reset(dvar_vector& x,int& ii);
    dvariable reset(dvar_vector& x,int& ii,d3_array& len_sample_size);
    void sel_dev_all_comp();
    void fmcomp(void);
    void effort_devs_calc();
    void do_fish_mort_intermediate_calcs(void);
    void do_fish_mort_intermediate_projection_calcs(void);
    void do_fish_mort_intermediate_calcs(int ir,int ip);
    void do_fish_mort_intermediate_calcs_q0(int ir,int ip);
      //(d3_array& len_sample_size,dvar_vector& vb_coff,dvar_vector& var_coff);
    void albacore_control_switches(void);
    dvariable calculate_the_totalbiomass_depletion(void);
    dvariable calculate_the_relative_biomass_depletion(void);
    dvariable calculate_the_average_biomass(void);
    dvector totalbiomass_depletion_for_plot_file(void);  // PK6-28-05
    void set_control_switches(void);
    void xinit(dvector& x);
    void xinit(dvector& x,int& ii);
    void xinit(dvector& x,int& ii,d3_array& len_sample_size);
    void xinit(dvector& x,int& ii,d3_array& len_sample_size,ofstream& xof);
    void do_recruitment_period_calculations(void);
    void get_initial_tag_recruitment_period(void);
    void get_sim_initial_tag_recruitment_period(void);
    int nvcal(d3_array& len_sample_size);
    void get_initial_age_structure(void);
    void get_initial_age_structure(
      const dvariable& totpop,dvar_vector& sv);
    dvar_fish_stock_history(int ntg,int nregions, int ng,ivector& nfp,
      imatrix& nfi,int nfsh,int nyrs,ivector& fl ,ivector& par_fl, ivector& nft,
      ivector& _regmin,ivector& _regmax,ivector& _dataswitch,imatrix& _Dflags,
      movement_info& mo,int df,int _mfactor,imatrix& ses_reg_flags,
      pmulti_species_data & pmsd,ivector & nage_by_region,
      ivector & nage_by_fishery);
    void get_initial_age_structure_equilibrium(void);
    void xget_initial_age_structure_equilibrium(void);
    void put_X_in_P(dvar_vector& X,dvar_matrix& P);
    void put_P_in_N(dvar_matrix& P);
    void put_P_in_N(dvar3_array& EN);
    dvar_matrix get_equilibrium_survival_rate(void);
    void print_tagging_fit_info(ofstream& of);
    void print_movement_report(ofstream& of);
    dvariable tag_catch_equations_calc_pooled(dvar_vector& sv);
    dvariable xtag_catch_equations_calc_pooled(dvar_vector& sv);
    void get_initial_tag_fishery_realization_index(void);
    void pooled_tag_catch_equations_calc(dvar_vector& sv);
    void xpooled_tag_catch_equations_calc(dvar_vector& sv);
    //void sim_xpooled_tag_catch_equations_calc(dvar_vector& sv);
    void grouped_catchability_calc(void);
    void do_grouped_tags_between_time_calculations(void);
    dvariable robust_kalman_filter_for_catchability(void);
  
    void process_effort_weight_data(void);
    void calculate_nat_mort_by_period(void);
    void effort_multiplier_for_cobb_douglas_production_function(void);
    void read_recruitment_env_data(adstring& root);
    void print_tag_return_by_time_at_liberty(ofstream& of);
    MY_DOUBLE_TYPE get_rep_rate_correction(int it,int ir,int ip,int fi);
    void explicit_seasonal_catchability_calc(void);
    dvar_matrix do_the_diffusion(dvar_matrix& N);
    dvariable normalize_seasonal_catchability(void);
    dvar_matrix catch_equations_calc_for_equilibrium
      (dvar_vector& sv, int navg,MY_DOUBLE_TYPE lambda,dvar_matrix& tmpcatch);
    dvariable calculate_reproductive_biomass(dvar3_array& NP,int iy);
    dvar_matrix catch_equations_calc_for_equilibrium
      (dvar_vector& sv, int navg,const prevariable & lambda,dvar_matrix& tmpcatch);
    dvar_matrix get_equilibrium_age_structure
      (int navg,int nomit,const MY_DOUBLE_TYPE lambda,dvar_matrix& tmpcatch);
    dvar_matrix get_equilibrium_age_structure
      (int navg,int nomit,const prevariable& lambda,dvar_matrix& tmpcatch);
    dvariable fit_implicit_catchability(void);
    void do_newton_raphson_for_tags(int it,int ir,int ip,dvariable& ffpen);
    void do_newton_raphson_for_tags_fix_zeroes(int it,int ir,int ip,dvariable& ffpen);
    void kludge_the_date(int & year,int & month);
    int is_projection_period(int ir,int ip);
    void set_year_flags_for_recruitment(void);
    int get_season(int i,int j);
    int get_season2(int i);
    void test_initial_equilibrium_code(ivector *);
    ivector  get_equilibrium_movements_per_season(void);
    ivector  kludged_get_equilibrium_movements_per_season(void);
    dvar_matrix get_initial_equilibrium_population1(dvar3_array& EN);
    dvar_vector get_initial_equilibrium_A(ivector& mps,
      dvar_matrix& ipop,imatrix& emi);
    dvar4_array  get_initial_equilibrium_survival(ivector & mps,int nya,
      ivector* pq_flag);
    dvar4_array  get_initial_equilibrium_survival_for_catch_conditioned
      (ivector & mps,int nya,ivector * pq_flag);
    dvar_matrix get_initial_equilibrium_B(ivector& mps,imatrix& emi);
      imatrix  get_initial_equilibrium_movement(ivector &mps);
      dvar_matrix get_initial_equilibrium_recruitment(int nya);
    dvar3_array get_equilibrium_cohorts(ivector& eqmi,dvar_matrix&,
      imatrix& emi);
    void get_equilibrium_cohorts_1(ivector& eqmi,dvar_matrix&,
      dvar4_array& tmort,imatrix& emi,dvar3_array &);
    void get_orthogonal_polynomial_weights(void);
      void get_orth_recr_info(void);
      void get_orthogonal_recruitment_poly_info(void);
      //void  get_new orthogonal_season_region_recruitment_poly_info(void);
      void new_get_orthogonal_recruitment_poly_info(void);
    void get_f_normal(MY_DOUBLE_TYPE & f,dmatrix& Y,dvector& x,
      dvector& g,banded_symmetric_dmatrix& H,int Hessflag,int robflag,
      dmatrix& resids,dmatrix & wts,int wtflag,ivector& ,imatrix&,MY_DOUBLE_TYPE c=0.05L);
    void get_f(MY_DOUBLE_TYPE & f,dmatrix& Y,dvector& x,
      dvector& g,banded_symmetric_dmatrix& H,int Hessflag,int robflag,
      dmatrix& resids,dmatrix & wts,int wtflag,ivector& ,imatrix&,MY_DOUBLE_TYPE c=0.05L);
    dvariable get_vf(dvar_matrix& Y,dvector& x,ivector&,imatrix&);
    dvariable get_vf_normal(dvar_matrix& Y,dvector& x,ivector&,imatrix&,MY_DOUBLE_TYPE);
    dvariable grouped_implicit_catchability_deviations_calc(void);
  
    dvariable grouped_implicit_catchability_deviations_calc_part2_pvm(void);
  
    MY_DOUBLE_TYPE grouped_implicit_catchability_deviations_calc_part1_pvm(int n,dvector& x,
      ivector&,imatrix&,MY_DOUBLE_TYPE);
    MY_DOUBLE_TYPE grouped_implicit_catchability_deviations_calc_part1(int n,dvector& x,
      ivector&,imatrix&,MY_DOUBLE_TYPE);
    void grouped_implicit_catchability_deviations_calc_part1_thread(void);
    dvariable grouped_implicit_catchability_deviations_calc_part2(dvector&,
      ivector&,imatrix&,MY_DOUBLE_TYPE);
    void grouped_implicit_catchability_deviations_calc_part1_pvm(void);
    dvariable grouped_implicit_catchability_deviations_calc_part2_thread(void);
    dvariable grouped_implicit_catchability_deviations_calc_part2_thread_with_signal(void);
    void bounds_for_ss2(dvar_vector & TF,const dvar_matrix & LF,
      const prevariable & fpen);
    dvar_vector ddnrf_interface(int ir,int ip,int nfi1,ivector& wtflag);
    dvar_vector ddnrf_interface_mc(int ir,int ip,ivector& wtflag);
    void cimplicit_seasonal_catchability_calc(dmatrix& Y);
    void vimplicit_seasonal_catchability_calc(dvar_matrix& Y);
    void set_grouped_true_months(void);
    dvar_vector get_seasonal_catchability_effect(int i);
    void do_newton_raphson_for_tags2(int it,int ir,int ip,dvariable& _ffpen);
    void new_do_newton_raphson_for_tags(int it,int ir,int ip,dvariable& _ffpen);
    void get_fishery_realization_index(void);
    void allocate_some_tag_stuff(void);
    void new_print_tagging_fit_info(ofstream& of);
    void new_print_tag_return_by_time_at_liberty(ofstream& of);
    dvariable get_sv_region(int is,int svind);
    dvariable get_sv_species(int is,int svind);
    dvar_vector get_nat_mort_species(void);
    dvar_matrix get_nat_mort_species(int);
    dvar_matrix get_nat_mort_region(int ir);
    dvariable get_nat_mort_coff_region(int ir);
    int get_age_flags_region(int ir,int i);
    dvar4_array  get_initial_equilibrium_survival_ms
      (ivector & mps,int nya,ivector * pq_flag);
    int  get_nage_region(int ir);
    int  get_nage(int is);
    dvar_vector get_age_pars_species(int is,int i);
    void pmsd_error();
    void indepvars_error();
    dvar_vector get_pmature_region(int ir);    //NMD 9May2012
    int get_current_species(void);
    int get_current_nage(void);
    ivector get_region_bounds(void);
    ivector get_region_bounds(int sp);
    dvar_vector get_pmature_species(int is);
    void test_tag_report(void *);  //NMD 14Feb2013
    void rescale_initial_tag_population(void);
    void set_selectivity_breaks(void);
    void add_blocked_sel_devs(void);
    dvar_matrix equilibrium_calcs(const dvar3_array & B,dvar_matrix & R);
    dvar3_array get_equilibrium_cohorts_0
      (ivector& mps,dvar_matrix& eqrec,dvar4_array& surv,imatrix& emi);
    dvar_matrix  get_initial_R(void);
    imatrix simyears;
    void init_recruit_equilibrium_code(void);
    void check_flag_stuff(void);
    void update_lagrange_lambda(void);
    void update_lagrange_avg_biomass(void);
    void read_fish_flags(cifstream *pinfile);
    ivector calculate_blocked_sel_sizes(void);
    void set_sel_seasons(void);
    dvariable get_bh_annual_recruitment_multiplier
      (dvar3_array& NP,int iy);
    void fish_flags_sanity_check(void);
    void age_flags_sanity_check(void);
    void tag_flags_sanity_check(void);
    void put_values_in_fish_mort(void);
    int new_fishing_selectivity_interface_nvar(void);
    void new_fishing_selectivity_interface_reset
      (int& ii,dvar_vector& x,dvariable& fpen,MY_DOUBLE_TYPE lb,MY_DOUBLE_TYPE ub,MY_DOUBLE_TYPE scale);
    void new_fishing_selectivity_interface_xinit
      (int& ii,dvector& x,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s);
    dmatrix  get_learner_code_coffs_low(void);
    int check_total_catch_pooling_flags(void);
    void  sim_get_initial_tag_population(void);
    void sim_xpooled_pd(dmatrix & prob,int & id,const char * s);
    dmatrix test_orthogonal(void);
    void expand_orthogonal_coefficients(void);
    void get_corresponding_orthogonal_coefficients(void);
    void setup_sel_ptr(void);
    void get_real_simtag_probs(dvar_vector& sv,int it);
    MY_DOUBLE_TYPE get_simulation_reporting_rate(int it,int ir,int ip,int fi,
      int real_flag);
    void modify_Dad(void);
    void modify_diff_coffs(dvar_matrix&);
    void modify_move_map(void);
    void modify_all_diff_coffs(void);
    int get_nage_species(int is);
    ivector nage_by_region;
    ivector nage_by_tag_release;
    ivector nage_by_fishery;
    imatrix tag_region_bounds;
    int get_nage_tag(int it);
    void get_first_unfixed_year(void);
    void setup_diffusion(void);
    void setup_diffusion2(void);
    void get_xdiff_parameterization(void);
    void get_diff_parameterization(void);
    dmatrix get_n_season_diffusion_matrix(void);
    dmatrix get_diffusion_matrix(void);
    int check_flat_tag_time_periods(int it);

//    dvar_vector censored_gamma(const prevariable & tau,dvar_vector& o,
//      const dvar_vector& mu,MY_DOUBLE_TYPE e);
    //group_manager_1* pgroup_manager_1;
    std::shared_ptr<group_manager_1> pgroup_manager_1;
    dvar4_array get_kludged_initial_equilibrium_survival
     (ivector & mps,ivector * pq_flag);
    int get_numyears_average(void);
    dvariable kludged_initial_survival_equilibrium_penalty(ivector *);
    void check_kludge_bounds(void);
    void get_the_initial_age_structure_options(ivector * pq_flag);
    int get_num_for_nr(int ir,int ip);
    int get_non_zero_effort(int ir,int ip);  //NMD_12mar2025
    dvar_matrix get_initial_equilibrium_plus_group
      (const dvar_matrix & M, ivector * pq_flag);
    void get_new_xdiff_parameterization(void);
    void get_old_diff_parameterization(void);
    void get_pos_hess_scale(int nvar);
    void read_psqrt_hess(int nvar);
    void test_movement_pointer_2(void);
    void get_movement_pointer(void);
    void test_movement_pointer_0(void);
    void test_movement_pointer_2_inv(void);
    void do_fragment_stuff(void);
    void assign_movement_coffs_option_0(void);
    void old_assign_movement_coffs_option_0(void);
    void assign_movement_coffs_option_2(void);
    void movement_coffs_option_0_to_2(void);
    void movement_coffs_option_2_to_0(void);
    void movement_coffs_option_1_to_0(void);
    void movement_coffs_option_0_to_1(void);
    void assign_movement_coffs_option_1(void);
    void rationalize_all_coffs(void);
    void rationalize_movement_coffs(int i);
    void get_corresponding_orthogonal_initial_recruitment_pars(void);
    void recinpop_standard(ivector * pq_flag,dvar3_array * pN,int override=0);
    dvariable recinpop_orth(dvar3_array *pN,int override);
    //void recinpop_orth(dvar3_array * pN,dvariable& fpen);
    dvector get_new_orthpolys(dvector& oldsel,dvector& tmpsel,dvector& xx,dvector& x1);
    MY_DOUBLE_TYPE orth_poly_fit(dvar_vector& vx);
    void recinpop_orth_xinit(ofstream& xof,dvector& x,int& ii);
    void recinpop_orth_reset(dvar_vector& x,int& ii);
    int recinpop_orth_size_count(void);
    void orth_pf155_less_than_zero_why(dvar3_array * pN);
    void set_totpop(void);
    void  copy_recruitment_to_population(dvar_matrix& Ninit,dvar3_array *pN);
    void fish_flags_sanity_check(int index,int group_index);
    void fish_flags_sanity_check(int index,int group_index,int override_index); //NMD_11apr2023
    void data_fish_flags_sanity_check(int index,int group_index);
    void orthpoly_recr_flags_sanity_check(void); //NMD_12dec2023
    void check_for_no_catch_or_effort(void);
    dvariable movement_coffs_option_2_penalty(void);
    dvariable recinpop_orth_penalties(void);
    imatrix get_survey_samples_flags(void);
    dvariable fit_survey_fishery_indices(ivector& ff92,int ps);
    void calculate_survey_indices(void);
    void get_tot_mort_with_tag_shedding(void);
    void do_tag_loss_fish_mort_intermediate_calcs(int ir,int ip);
    void do_tag_loss_fish_mort_intermediate_calcs(void);
    i3_array fishery_group_ptr;
    void get_fishery_group_pointers(void);
    void make_fishery_group_pointer(int i);
    dvar4_array get_initial_equilibrium_survival_natmort
      (ivector & mps,ivector * pq_flag);

    void print_survey_fishery_fits(ivector& ff92);
    void set_fish_mort_q0(double u);
    void do_all_sanity_checks(void);

    dvector indep_var;
    dvector indep_var_lo;
    dvector indep_var_hi;
    void construct_indep_vars(int nvar);
    void get_indep_vars_for_report(dvar_vector& x,int& ii, int& runsetv);
    void recinpop_orth_reset_indepvar(dvar_vector& x,int& ii,int& runsetv);

    int historical_version;
    int current_version;

    MY_DOUBLE_TYPE get_disagg_catch(int ir, int ip, int fi, int wtf,
      dvar_vector enf);
    
  }; // end class dvar_fish_stock_history

  // new stuff to manage changin grouping flags during fitting
  class group_manager_1
  {
  public:
    imatrix old_fish_flags;
    imatrix old_zero_effdev_flags;
    ivector key;
    const ivector active_flags;
    const imatrix mflags;
    const ivector group;
    ivector active_group;
    int maxg;
    int groupsum;
    group_manager_1(void);
    group_manager_1(ivector & ff4,imatrix& mflags,
      const ivector& group);
    ~group_manager_1(); 
  };

    //fish_stock_history value(dvar_fish_stock_history& dfsh);
    ostream& operator << (ostream& ofs, dvar_fish_stock_history&);

class tail_compression_stuff
{
public:
  ivector sizes;
  ivector num_samples;
  ivector len_rho_group_ptr;
  ivector inv_len_rho_group_ptr;
  dvector average_sample_size;
  int tot_nfi;
  int minsize;
  ivector reg;
  ivector per;
  ivector finc;
  ivector left_bound;
  ivector right_bound;
  imatrix group_size_flag;
  ivector is;
 /*
  tail_compression_stuff(int num_fisheries): sizes(0), num_samples(;
   len_rho_group_ptr;
   inv_len_rho_group_ptr;
   average_sample_size;
   tot_nfi;
   minsize;
   reg;
   per;
   finc;
   left_bound;
   right_bound;
  */
};

class wght_tail_compression_stuff
{
public:
  ivector sizes;
  ivector num_samples;
  ivector wght_rho_group_ptr;
  ivector inv_wght_rho_group_ptr;
  dvector average_sample_size;
  int tot_nfi;
  int minsize;
  ivector reg;
  ivector per;
  ivector finc;
  ivector left_bound;
  ivector right_bound;
  imatrix group_size_flag;
  ivector is;
 /*
  tail_compression_stuff(int num_fisheries): sizes(0), num_samples(;
   len_rho_group_ptr;
   inv_len_rho_group_ptr;
   average_sample_size;
   tot_nfi;
   minsize;
   reg;
   per;
   finc;
   left_bound;
   right_bound;
  */
};


class age_length_record
{
public:
  int year;
  int month;
  int fishery;
  int species;
  int region;
  int period;
  int incident;
  
  dmatrix age_length;
  void allocate(int nlint,int nage)
  {
    age_length.allocate(1,nlint,1,nage);
  }
};

class age_length_record_array
{
public:
  dmatrix obs_sample_size_by_length;
  MY_DOUBLE_TYPE sumwght;
  dvector effective_sample_size;
  int * ncopies;
  age_length_record * ptr;
  int mmin;
  int mmax;
  //age_length_record_array(int mmin,int mmax,int nlint,int nage);
  void allocate(int mmin,int mmax,int nlint,int nage);
  age_length_record_array(const age_length_record_array & alarray);
  age_length_record_array(void);
  ~age_length_record_array();
  int indexmin() { return mmin; }
  int indexmax() { return mmax; }
  age_length_record & operator () (int i)
  {
    if (i<mmin || i>mmax)
    {
      cerr << "bounds error in age length record array" << endl;
      ad_exit(1);
    }
    return ptr[i];
  }
};

  class dvar_len_fish_stock_history: public dvar_fish_stock_history
  {
  public:
    dmatrix fishing_incident_times;
    dmatrix fishing_incident_real_times;  //NMD_27jun2022
    dvector cpmature_at_length;
    tail_compression_stuff * tcs;
    wght_tail_compression_stuff * wtcs;
    dmatrix m_estimator_coffs;
    dvar_matrix tlength;
    d4_array lestp;
    d4_array wtestp;
    d4_array lestp1;
    d4_array wtestp1;
    dvar4_array bswsplinesel;
    dvar4_array bsplinesel;
    age_length_record_array alarray;
    dvar5_array qij;   // probability that an age class j fish
                     // has length in length interval i
    print_double p_Fmmsy, p_Fmmsy_pen; // for target F to Fmsy
    print_double p_BBmsy, p_BBmsy_pen; // for target B to Bmsy
    print_double p_sBsBmsy, p_sBsBmsy_pen; // for target spawning B to Bmsy
    //print_double p_Bmsy,p_sBmsy; 
    MY_DOUBLE_TYPE p_Bmsy,p_sBmsy; //  
    print_double p_Fmsy; // exploitation rate at MSY
    MY_DOUBLE_TYPE p_MSY;  
    //print_double p_MSY;
    dvariable MMSY;  //    
    dmatrix yield_eq_region; //for output to plot.rep, calc. in yield_analysis_bh_daves_folly
    dvector yield_eq_region_x; //corresponding Fmult values for output to plot.rep, 
    MY_DOUBLE_TYPE steep_penalty;      // steepness penalty
    ivector nlintv;
    //dmatrix freq;
    dvar_matrix splinesel;
    dvar_matrix wsplinesel;
    dvector fmid;
    dvector wmid;
    dvector wm2;
    dvector realwmid;
    dvar_vector wljac;
    dvar_matrix sdevs_yr;
    //dvar5_array age_weight_sel;
    dvar5_array age_weight_fishery_sel;
    dvar5_array age_length_fishery_sel;
    dvar5_array age_weight_fishery_size_dist;
    dvar5_array age_length_fishery_size_dist;
    dvar4_array age_fishery_sel_from_length;
    dvar4_array age_fishery_sel_from_weight;
    dvar4_array mean_length;
    dvar_matrix mean_length_yr_proj;
    dvar3_array mean_length_yr;
    MY_DOUBLE_TYPE len_wt_coff;
    dvar_matrix length_month_yr;
    dvar_matrix length_sel;
    dvar3_array relative_bio;
    dvar3_array      len_dist;      
    dvar3_array      RB;      
    d3_array      ORB;      
    dvar4_array vars;
    dvar_vector global_vars;
    dvar4_array sdevs;
    MY_DOUBLE_TYPE shlen;
    MY_DOUBLE_TYPE wshlen;
    MY_DOUBLE_TYPE wfilen;
    MY_DOUBLE_TYPE filen;
    dvar_vector vb_coff;
    dvar_vector vb_bias;
    ivector common_vb_bias;
    dvar_vector common_vb_bias_coffs;
    dvar_vector growth_dev;
    dvar_vector var_coff;
    dvar4_array tprob;
    dvar4_array tc_tprob;
    dvar4_array wtprob;
    dvar4_array tc_wtprob;
    dvar_matrix lengthbsel;   
    d4_array len_freq;
    d4_array len_eps;
    d4_array tc_len_freq;
    d4_array tc_wght_freq;

    d4_array wght_eps;
    d4_array wght_freq;
    d4_array age_freq;
    d3_array max_len_obs;
    d3_array max_wght_obs;
    d3_array age_len_sample_size;
    d3_array age_sample_size;
    int nlint;
    int nwint;
    MY_DOUBLE_TYPE fmin1;
    MY_DOUBLE_TYPE fmax1;
    MY_DOUBLE_TYPE fminl;
    MY_DOUBLE_TYPE fmaxl;
    MY_DOUBLE_TYPE vmin1;
    MY_DOUBLE_TYPE vmax1;
    MY_DOUBLE_TYPE vminl;
    MY_DOUBLE_TYPE vmaxl;
    MY_DOUBLE_TYPE rhomin;
    MY_DOUBLE_TYPE rhomax;
    int nfmbound;
    ivector ifper;
    ivector ifinc;
    ivector iffish;
    ivector iageclass;
    dvector fmmin;
    dvector fmmax;
    dvector pdown;
    dvector pup;

    friend dvariable objective_function(dvar_len_fish_stock_history& fsh);

  public:
   void  big_fish_catch(void);
    void do_everything_calc(dvariable&,ivector* pq=NULL);
    void mean_length_calc(void);
    void mean_length_calcx(void);
    void variance_calc(void);
    void predicted_frequency_calc(void);
    void fast_pred_frequency_calc(void);
    void fast_weight_pred_frequency_calc(void);
    void fast_pred_frequency_calc_len_based(void);
    dvariable reset(dvar_vector& x);
    int nvcal(dvector&x);
    friend par_uostream& operator << (par_uostream& pof,
			   dvar_len_fish_stock_history& fsh);
    friend par_ofstream& operator << (par_ofstream& pof,
			   dvar_len_fish_stock_history& fsh);
    dvariable robust_fit(void);
    dvariable fmeanpen(void);
    void obs_length_moments_calc(
      dvar_matrix& obs_catch_moment1,dvar_matrix& obs_catch_moment2);
    void set_global_variance(void);
    void pred_length_moments_calc(dvar_matrix& pred_catch_moment1,
      dvar_matrix& pred_catch_moment2);
    void main_length_calcs(dvar_vector& tprobb,dvar_vector& propp,
      dvar_vector& mean_len,dvar_vector& sigg);
    void main_length_calcs2(dvar_vector& tprobb,
      dvar_vector& propp,dvar_vector& mean_len,dvar_vector& sigg,
      dvar_vector& relative_num_at_length,dvar_vector& length_sel,
      dvar_vector& nnum_fish,dvar_vector& llen_dist,int fi);
    dvar_vector calculate_the_current_biomass(ivector& rip);
    void main_length_calcs_len_based
    (dvar_vector& tprobb,dvar_vector& propp,dvar_vector& mean_len,
       dvar_vector& sigg,dvar_vector& nnum_fish,dvar_vector& llendist,
       int fi);
    d3_array report_catch_by_year(void);
    d3_array report_estimated_catch_by_year(void);

    void calculate_the_biomass_by_region(int ir,int i);

    dvariable fit_tag_returns_parallel(void);
    void main_length_calcs_print(uostream& ofs,
      dvar_vector& propp,dvar_vector& mean_len,
      dvar_vector& sigg);

    dvector get_mean_length(int k)
    {
      return  value(mean_length( fishing_region(k),fishing_period(k), 
        fishing_incident(k)) );
    }

    dvector get_props(int k)
    {
      return  value(prop( fishing_region(k),fishing_period(k), 
        fishing_incident(k) ));
    }

    dvector get_vars(int k)
    {
      return  value(vars(fishing_region(k),fishing_period(k), 
        fishing_incident(k) ));
    }
    d4_array vb_length_calc(void);
    dvector vb_length_calc(int,int,int);

    void print_pred_frequencies(ofstream& ofs);
    void print_pred_frequencies(uostream& ofs);

    dvector main_length_calcs_print(ofstream& ofs,
      dvar_vector& propp,dvar_vector& mean_len, dvar_vector& sigg,
      int ir,int ip,int fi,int pi);
    int nvcal(void);

    //void xinit(dvector& x,int& ii,d3_array& len_sample_size,ofstream& xof);

    void set_Takeuchi_san_control_switches(void);
    dvariable fit_tag_returns(void);
    dvariable fit_tag_returns_like_ss3(void);
    dvariable fit_pooled_tag_returns_like_ss3(void);
    dvariable fit_tag_returns2(void);
    dvariable fit_tag_returns_survival_analysis(void);
    dvariable fit_tag_returns_sqrt(void);
    void read_tagging_data(adstring& root,pmulti_species_data & pmsd);
    void set_control_switches(void);
    void set_shark_control_switches(void);
    void set_shark_no_vonb_control_switches(void);
    void set_shark_no_vonb_var_control_switches(void);
    void set_Yukio_control_switches(void);
    void set_gridsearch_control_switches(void);
    void albacore_control_switches(void);
    void xinit(dvector& x);
    dvar_len_fish_stock_history(dvar_fish_stock_history& fsh,
       int _nlint, MY_DOUBLE_TYPE& _shlen,MY_DOUBLE_TYPE& _filen );
    dvar_len_fish_stock_history(const dvar_fish_stock_history& fsh,
       int _nlint, MY_DOUBLE_TYPE& _shlen,MY_DOUBLE_TYPE& _filen,i3_array& x,
       i3_array& u, i3_array& v, int _nwint);
    void cons_convert_tag_lengths_to_age(void);
    void var_convert_tag_lengths_to_age(void);
    void cfast_pred_frequency_calc(dvar_len_fish_stock_history& cfsh);
    void set_some_flags(ivector& dataswitch);
    void set_some_flags(ivector& dataswitch,int min_year);
    void cons_observed_tag_catch_calc(void);
    void observed_tag_catch_calc(int direction_flag);
    void observed_tag_catch_by_length_calc(void);
    void observed_tags_by_age_from_length(void);
    void observed_tags_by_age_from_length_pooled(void);
    void tot_tags_catch(void);
    void print_tag_data(int ir,int ip,ofstream &ofs);
    dvariable fit_pooled_tag_returns(void);
    dvariable fit_pooled_tag_returns_ss3(void);
    dvariable grouped_tag_reporting_rate_penalty(void);
    dvariable biomass_dynamics_pt(void);
    void yield_analysis_pt(ofstream * pof);
    void calculate_the_biomass(void);
    void calculate_the_biomass_by_region(void);
    void calculate_the_catch_biomass(void);
    void calculate_the_catch_numbers(void);    // JH 27/03/02
    void get_fishing_mortality_by_age_by_year_by_region(void);
    void get_fishing_mortality_by_age_by_year(void);
    void get_fishing_mortality_by_age_by_year(int i);
    void calculate_the_mean_weight(void);
    void mean_lengths_by_year_calc(void);
    void mean_weights_by_year_calc(void);
    void get_equilibrium_structure_for_yield(void);
    dvariable fit_tag_returns_mix(void);
    dvariable fit_pooled_tag_returns_mix(void);
    void tag_returns_report(const ofstream &);
    void main_weight_calcs(dvar_vector& tprobb,dvar_vector& propp,
      dvar_vector& mean_len,dvar_vector& sigg);
    void setup_some_stuff(MY_DOUBLE_TYPE tmp_len_wt_coff,
      MY_DOUBLE_TYPE _wshlen,MY_DOUBLE_TYPE _wfilen,int nwint,int _num_fisheries,imatrix _dff,
      par_cifstream * _pinfile,int _month_1,int _first_time,int _mfactor,
      ivector& _v,int nn,pmulti_species_data & pmsd);
    void print_pred_wght_frequencies(ofstream& ofs);
    dvector main_wght_calcs_print(ofstream& ofs,
      dvar_vector& propp,dvar_vector& mean_len,dvar_vector& sigg,
      int ir,int ip,int fi,int pi);   //NMD_24jan2023
    void set_effective_length_and_weight_sizes(void);
    void do_everything_calc(d3_array& len_sample_size,
      dvar_vector& vb_coff,dvar_vector& var_coff,dvar_vector& sv,
      dvariable&,ivector* pq=NULL);
//    dvariable fit_totals(int acf);
    dvariable fit_totals(int print_switch,int acf); //NMD_11dec2023
    void get_population_multipliers(dvar_vector& sv,void *);
    //void catch_equations_calc(dvar_vector& sv,ivector * pq_flag);
    void catch_equations_calc2(dvar_vector& sv,ivector * pq_flag,
      const prevariable& );
    void catch_equations_calc1(dvar_vector& sv,ivector * pq_flag,
      const prevariable& );
    void catch_equations_calc(dvar_vector& sv,ivector * pq_flag,
      const prevariable& );
    void catch_equations_calc_implicit(dvar_vector& sv,dvariable& ffpen);
//    void new_get_numbers_at_age_implicit(dvar_vector& sv,void * pq_flag,
//      dvariable& ffpen);
    dvar_matrix new_get_numbers_at_age_implicit(dvar_vector& sv,void * pq_flag,
      dvariable& ffpen);
    void new_catch_equations_calc_implicit_experiment_loop(dvar_vector& sv,
      ivector* pq_flag,dvariable& ffpen);
    void new_catch_equations_calc_implicit_experiment(dvar_vector& sv,
      ivector* pq_flag,dvariable& ffpen);
    void new_catch_equations_calc_implicit(dvar_vector& sv,ivector* pq_flag,
      dvariable& ffpen);
    void dvar_len_fish_stock_history::q0debug_report(ivector* pq_flag,
      ofstream& ofs);
      // NMD_18jun2020
    //void new_catch_equations_calc_implicit(dvar_vector& sv,dvariable& ffpen);
    void do_newton_raphson(int ir,int ip,dvariable& ffpen,dvar_vector& q,ivector&,dvar_matrix&,dvariable& fpen1);
    void do_newton_raphson_with_totcatch(int ir,int ip,dvariable& ffpen,
      dvar_vector&);
    void do_newton_raphson_with_totcatch(int ir,int ip,dvariable& ffpen,
      dvar_vector&,dvariable& ffpen1);

    void do_newton_raphson_missing_totcatch2(int ir,int ip,dvariable& ffpen,
      dvar_vector& q);
    void do_newton_raphson_missing_totcatch(int ir,int ip,dvariable& ffpen);
    void do_newton_raphson_missing_totcatch(int ir,int ip,dvariable& ffpen,
      dvar_vector& q,ivector&,dvar_matrix&);
    void do_newton_raphson_missing_totcatch_fml(int ir,int ip,dvariable& ffpen,
      dvar_vector& fmlev,ivector& mci,dvar_matrix& fmlevdevs);

    void catch_equations_calc_implicit_mc(dvar_vector& sv,dvariable& ffpen);
    dvar_vector get_kludged_catch_with_totcatch(int ir,int ip,
      dvariable& ffpen,dvar_vector& enf);
    dvariable get_mean_weight(int ir,int ip,int fi);
//    dvar_vector get_numbers_at_age(dvar_vector& sv,void * pq_flag);
    dvar_matrix get_numbers_at_age(dvar_vector& sv,void * pq_flag); //NMD_13Apr2021
    void mean_lengths_by_year_for_projections(void);
    void mean_weights_by_year_for_projections(void);
    void get_initial_projection_period(void);
    void mean_length_calc(int mode_flag);


    void mean_weights_calc(int mode_flag);
    void  resize_if_necessary(void);
    void create_catchability_deviation_thread(void);
    dvariable tag_catch_equations_calc(dvar_vector& sv);
    void do_newton_raphson_for_tags_with_totcatch(int it,int ir,int ip,
      dvariable& ffpen, dvar_vector& qimp);
    void zero_out_catches(void);
    void zero_out_log_effort(void);   //NMD_13jul2022
    void do_newton_raphson(int ir,int ip,dvariable& ffpen,dvar_vector& fmlev,
      ivector& mci,dvar_matrix& fmlevdevs,ivector * pq_flag,dvariable& ffpen1);
    void do_newton_raphson_missing_totcatch(int ir,int ip,dvariable& ffpen,
      dvar_vector& fmlev,ivector& mci,dvar_matrix& fmlevdevs,ivector * pq_flag);
    void incident_selectivity_calc(d3_array&);
    void incident_selectivity_calc(d3_array&,dvar_vector&,dvar_vector&);
    void new_incident_selectivity_calc(d3_array&,dvar_vector&,dvar_vector&);
  void length_dist_calcs_optimized
  (dvar_vector& tprobb,dvar_vector& propp,dvar_vector& mean_len,
    dvar_vector& sigg,dvar_matrix& r);
  dvar_matrix length_dist_calcs
   (dvar_vector& mean_len,dvar_vector& sigg);
  dvar_matrix length_dist_calcs(void);
  void all_length_dist_calcs(void);
  void all_length_dist_calcs_long(void);
  void all_weight_dist_calcs_long(void);
  void allocate_optional_stuff(void);
  void main_length_calcs_len_based(dvar_vector& tprobb,dvar_vector& propp,
    dvar_matrix& alis);
  dvar_matrix weight_dist_calcs(dvar_vector& mean_len,dvar_vector& sigg);
  void main_weight_calcs_len_based(dvar_vector& tprobb,
    dvar_vector& propp,dvar_matrix& alis);
  dvar_vector get_global_vars(int cs);
  dvar_vector get_global_vars_region(int cs);
  dvar_vector get_vb_coff_region(int cs);
  dvar_vector get_vb_coff(int is);
  MY_DOUBLE_TYPE  get_len_wt_coff(int cs);
  dmatrix get_multi_wm2(void);
  dmatrix get_multi_wmid(void);
  dvar_vector get_var_coff_species(int is);
  dvar_matrix get_sdevs_yr_species(int is);
  MY_DOUBLE_TYPE get_len_wt_coff_region(int ir);
  xxxtags make_tagstuff(int sno);
  int send_variables_to_tag_slave1(int slave_number,int flag);
  int send_variables_to_tag_slave3(int slave_number,int flag);
  int send_variables_to_tag_slave2(int slave_number,int flag);
  void send_constant_data(int sno);
  void send_constant_data_3(int sno);
  void all_age_length_calcs(void);
  void age_length_calcs(dvar_matrix& q,dvar_vector& mean_len,
    dvar_vector& sigg);
  void get_age_length_periods(void);
  dvariable fit_age_length_data(void);
  void print_likelihood_components(ofstream& of1);
  void allocate_plotstuff(void);
  void simulate_length_and_weight_frequencies();
  void generate_simulated_cpue();  //NMD_24aug2023
  void make_tail_compressed_samples(void);
  void make_tail_compressed_weight_samples(void);
  void tail_compress_predicted_length_frequencies(void);
  void tail_compress_predicted_weight_frequencies(void);
  int check_freq_pooling_flags(int data_type);

  dvariable square_fita(d3_array& tf,d3_array & sample_size,d4_array & of,
    dvar4_array & pf,dvector& contrub,d3_array & contrib_by,int data_type,
    int print_switch);
  dvariable square_fita_modified(d3_array& tf,
    d3_array & sample_size,d4_array & of,dvar4_array & pf,
    dvector& contrib, d3_array& contrib_by_realization);
  void splines(void);   // i'th fishery
  void splines_common_length(void);   // i'th fishery
  void no_splines(void);
  void new_set_selectivity_breaks(void);
  dvariable len_self_scaling_multinomial_re(dvar_len_fish_stock_history& fsh);
  //dvariable wt_self_scaling_multinomial_re(dvar_len_fish_stock_history& fsh);
  void print_ssmult_stuff(void);
  void natural_mortality_lorenzen(void);
  void natural_mortality_double_normal(void);
  void do_simulation_stuff(void);
  void read_tag_simulation_info(void);
  dvar_vector getmin(int fi); 
  void get_scaled_two_species_lengths(void); 
  void get_scaled_lengths(void); 
  void length_logistic(void); 
  void get_scaled_single_or_two_species_lengths(void);
  void generate_simulated_tags(void);
  void generate_simulated_real_tags(void);
  void generate_simulated_tags_all(void);
  dvariable get_simtag_pd(int it, int iage, int irr,
    dmatrix & prob,int & id,const char * s);
  dvariable get_simtag_pd_all(int it, dvector& v, int irr,
    dmatrix & prob,int & id,const char * s);
//  void  make_tag_data_report(dmatrix& data,int ii,char s[]);
  void  make_tag_data_report(dmatrix& data,int ii,const char * s);
//  void  make_real_tag_data_report(dmatrix& data,int ii,char s[]);
  void  make_real_tag_data_report(dmatrix& data,int ii,const char * s);
//  ivector get_sim_num_tags_at_length(int it,ivector sim_num_tags_at_age);
  ivector get_sim_num_tags_at_length(int it,ivector sim_num_tags_at_age,
    random_number_generator& rng);   //NMD_20dec2019
  dmatrix get_initial_tag_release_age_dist(void);
  dvariable get_real_simtag_pd2
    (int it, int iage, int irr,dmatrix & prob,int & id,const char * s);
  void generate_real_simulated_tags2(void);
  void generate_real_simulated_tagsx(void);
  MY_DOUBLE_TYPE get_simulated_real_tag_fish_age(void);
  int check_if_tag_exists(MY_DOUBLE_TYPE d,random_number_generator &rng);
  void check_if_tag_survived_and_was_caught(int it,int iage,
    random_number_generator& rng,int il,const dmatrix & _data1,int & ii);
  MY_DOUBLE_TYPE get_simulated_real_tag_fish_age(int il,
    random_number_generator& rng);
  dvar_vector maturity_length_to_age_simple_spline(int isp);
  dvar_vector maturity_length_to_age_weighted_spline(int isp);
  void  make_fake_recapture(int it, dmatrix & _data1,int & ii);
  void check_for_ss3_structure(void);
  dvar_vector censored_gamma_or_lognormal
    (const prevariable & tau,dvar_vector& o,const dvar_vector& mu,MY_DOUBLE_TYPE e,
     int pf111);
  dvar_vector dvar_len_fish_stock_history::censored_gamma_mix(
    const prevariable & tau,dvar_vector& o,const dvar_vector& mu,MY_DOUBLE_TYPE e);
  dvar_vector dvar_len_fish_stock_history::censored_lognormal_mix(
    const prevariable & tau,dvar_vector& o,const dvar_vector& mu,MY_DOUBLE_TYPE e);
  dvar_vector censored_gamma(const prevariable & tau,dvar_vector& o,
      const dvar_vector& mu,MY_DOUBLE_TYPE e);

  void do_newton_raphson_missing_effort(int ir,int ip,dvariable& ffpen);

  void check_implict_catch_options(int ir,int ip,dvariable& ffpen,
   dvar_vector& fmlev,ivector& mci,dvar_matrix& fmlevdevs,
   dvariable& ffpen1);
   dmatrix get_time(void);
   dmatrix get_real_time(void);  //NMD_27jun2022
   void new_build_implicit_catch_effort_design_matrix(void);
   dvariable old_build_implicit_catch_effort_design_matrix(void);
   dvariable implicit_catch_effort_relationship(void);
   void do_old_build_part(int fi,dvar_vector& O,int& irow);
   void new_do_build_part(int fi,dmatrix& M,int& irow,int& icol);
   void new_do_build_part_for_projections1(void);   //NMD_30jun2022
   void new_do_build_part_for_projections2(void);   //NMD_4jul2022
   void do_col_bounds_count(int fi,int& icol);
   void do_row_bounds_count(int fi,int& irow);
   imatrix get_bounds_for_implicit_catch_effort_design_matrix(void);
   void catch_conditioned_popes_approx(int ir,int ip,dvariable& ffpen);
   void do_popes_approximation(int ir,int ip,dvariable& ffpen,
     dvariable& ffpen1);
   dvar_vector get_mean_weight_at_age(int ir,int ip,int fi);
   dvar4_array mean_weight_calcs(void);
   void do_newton_raphson_missing_totcatch(int ir,int ip,
     dvariable& ffpen,dvar_vector& fmlev);
   void new_trunc_newton_code(void);
   void remove_extra_columns_from_fml_designvec(int fi);

   void put_terminal_catchability(int ir,int ip);  //NMD_7may2024

   void get_indep_vars_for_report(dvar_vector& x);
 };  // end class dvar_len_fish_stock_history

dvariable prob_calc(MY_DOUBLE_TYPE& fmidd,prevariable& mean_len,prevariable& sigg);
dvariable age_at_length_calc(const prevariable& v,const dvar_vector& vb_coff,
  int nage);
dvariable daves_kludge(prevariable& x);
MY_DOUBLE_TYPE daves_kludge(MY_DOUBLE_TYPE x);

  MY_DOUBLE_TYPE fcomp(const dvar_len_fish_stock_history& fsh,const dvar_vector& x,
    int nvar,int print_switch,const dvector& gbest,int grad_switch,
    ivector * pq_flag=NULL);

void do_robust_mean_spread_calcs(dvar_fish_stock_history& fsh);

void  calc_FF_mean_partial(dvar_fish_stock_history& fsh,int print_switch=0);

class date_struc
{
  public:
  int year;
  int month;
  int week;
};

date_struc get_new_time(int oy,int om,int ow,int mult,int first_time);
date_struc get_old_time(int oy,int om,int ow,int mult,int first_time);

  dvariable age_at_length_calc(MY_DOUBLE_TYPE v,const dvar_vector& vb_coff,int nage);
  MY_DOUBLE_TYPE age_at_length_calc(MY_DOUBLE_TYPE v,const dvector& vb_coff,int nage);

  int check_pos(const dvar_vector& qq,int n);

  void slave_fitting_procedure(dvar_len_fish_stock_history& fsh);
  void  mfget_dv3_sum_from_slaves(const dvar3_array& w);

void set_value_partial(const dvar3_array& _w,const dvar_vector& x,
  const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,const ivector& flags,const ivector& group,
  MY_DOUBLE_TYPE scale,ivector& block_ptr,ivector& ff57,int nage,ivector& ff61,
  ivector& sel_block_index);
  
void set_value_inv_partial(const dvar3_array& _w,const dvector& x,
  const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const ivector& flags,const ivector& group,
  MY_DOUBLE_TYPE scale,ivector& block_ptr,ivector& ff57,int nage,ivector& ff61,
  ivector& sel_block_index);
  
int size_count_partial(const dvar3_array& _w,
  const ivector& flags,const ivector& group,
  ivector& block_ptr,ivector& ff57,int nage,ivector& ff61,
  ivector& sel_block_index);


#endif
