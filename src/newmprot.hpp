/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

    dvariable square_fita(dvar_len_fish_stock_history& fsh,d3_array& total_freq);

    dvariable square_fit(dvar_len_fish_stock_history& fsh,d3_array& total_freq);

    dvariable square_fit0(dvar_len_fish_stock_history& fsh,d3_array& total_freq);
    dvariable square_fit0a(dvar_len_fish_stock_history& fsh,d3_array& total_freq);

    dvariable mean_fit(dvar_len_fish_stock_history& fsh,
      d3_array& total_freq);
    dvariable objective_function(dvar_fish_stock_history& fsh);
    dvariable objective_function(dvar_len_fish_stock_history& fsh,
      dvariable& mean_length_constraints,int print_switch);
    dvariable robust_freq_fit(d3_array & obs, dvar3_array& pred,
      dvar3_array& weights);
dvariable robust_freq_fit(d3_array & obs, dvar3_array& pred,
  dvar3_array& weights, MY_DOUBLE_TYPE& min_sig);
#ifdef __MSVC32__
  MY_DOUBLE_TYPE mmiinn(const MY_DOUBLE_TYPE& x,const MY_DOUBLE_TYPE y);
#else
  MY_DOUBLE_TYPE mmiinn(const MY_DOUBLE_TYPE& x,const MY_DOUBLE_TYPE y);
#endif
dvariable shmodel_fit(dvar_len_fish_stock_history& fsh,int print_switch);
dvariable robust_freq_partial_fit(d4_array & obs, dvar4_array& pred,
  dvar4_array& weights,dvar_len_fish_stock_history& fsh);

dvariable eff_dev_penalty(dvar_fish_stock_history& fsh,int print_switch);
dvariable implicit_eff_dev_penalty(dvar_fish_stock_history& fsh,int print_switch);
dvariable selectivity_form_penalty(dvar_len_fish_stock_history& fsh,
  int print_switch);
dvariable selectivity_deviations_penalty(dvar_fish_stock_history& fsh,
  int print_switch);
dvariable incident_sel_curvature_penalty(dvar_fish_stock_history& fsh,
  int print_switch);

dvariable catch_at_age_fit(dvar_len_fish_stock_history& fsh,int print_switch);
dvariable catch_at_age_fit_dirichlet_multinomial(dvar_len_fish_stock_history& fsh,int print_switch);
dvariable catch_at_length_fit(dvar_len_fish_stock_history& fsh,
  int print_switch);
d3_array total_freq_calc(dvar_len_fish_stock_history& fsh,int print_switch);
dvar4_array weights_calc(dvar_len_fish_stock_history& fsh,d3_array& total_freq,
  int print_switch);
dvariable recr_initpop_sum_penalty(dvar_fish_stock_history& fsh,
  int print_switch);
dvariable first_length_bias_penalty(dvar_len_fish_stock_history& fsh,
  int print_switch);
dvariable no_penalties(dvar_len_fish_stock_history& fsh,
  dvariable& mean_length_constraints,int print_switch,ofstream* of_pen);
dvariable good_penalties(dvar_len_fish_stock_history& fsh,
  dvariable& mean_length_constraints,int print_switch,ofstream* of_pen);
dvariable call_penalties(dvar_len_fish_stock_history& fsh,
  dvariable& mean_length_constraints,int print_switch,ofstream* of_pen);
dvariable total_catch_fit(dvar_fish_stock_history& fsh,
  int print_switch,int);
dvariable total_weight_fit(dvar_fish_stock_history& fsh,int print_switch,int);

dvariable ratio_first_last_biocalc(dvar_len_fish_stock_history& fsh);
void length_calc(dvariable& tt,prevariable& vb1,
    prevariable& vb2,prevariable& t1,prevariable& t2);
void write_the_parfile(const adstring& fopp);
dvariable multinomial_fit(dvar_len_fish_stock_history& fsh, 
  d3_array& total_freq);
dvariable dirichlet_multinomial_fit(dvar_len_fish_stock_history& fsh, 
  d3_array& total_freq);
dvariable dirichlet_multinomial_mixture_fit(dvar_len_fish_stock_history& fsh, 
  d3_array& total_freq);
dvariable logistic_normal_fit(dvar_len_fish_stock_history& fsh, 
  d3_array& total_freq);
dvariable tail_compressed_logistic_normal_fit_heteroscedastic(dvar_len_fish_stock_history& fsh, 
  d3_array& total_freq);
dvariable tail_compressed_weight_logistic_normal_fit_heteroscedastic(dvar_len_fish_stock_history& fsh, 
  d3_array& total_freq);
dvariable logistic_normal_fit_heteroscedastic(dvar_len_fish_stock_history& fsh, 
  d3_array& total_freq);
dvariable logistic_normal_weight_fit(dvar_len_fish_stock_history& fsh, 
  d3_array& total_freq);
dvariable logistic_normal_weight_fit_heteroscedastic
  (dvar_len_fish_stock_history& fsh,d3_array& total_freq);

void  get_implicit_catchability(dvar_fish_stock_history& fsh);
void  get_implicit_catchability_catch_conditioned(dvar_fish_stock_history& fsh);

extern struct pvmhostinfo *ad_hostp;
extern int ad_nhost;
dvar_vector choleski_factor_solve(_CONST dvar_matrix& MM,
    const dvar_vector& vv,const prevariable& det,const int& sgn);
void set_option_flag(const char * s,int& ss,int argc,char * argv[]);
