/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

class xxxmspd 
{
public:
  // data
  int num_species;
  int num_real_regions;
  ivector region_species_pointer;
  ivector tag_species_pointer;
  imatrix region_bounds;
  // variables
  dvar3_array nat_mort; 

  dvar_matrix vb_coff;
  dvar_matrix var_coff;
};

typedef xxxmspd * pxxxmspd;

dvar_matrix fast_diffusion_calcs(int nage,int num_regions,dvar3_array& N,
  dvar3_array& Dad,ivector& rip,int maxage,pxxxmspd & pmsd);

dvar_matrix fast_diffusion_calcs(int nage,int num_regions,dvar3_array& N,
  dvar3_array& Dad,int year,pxxxmspd & pmsd);

//double mmaxZ[5];

dvar_matrix get_jac2(int nfi,int jmin,int nage,dvar_vector q,
  dvar_matrix sel,const dvar_vector& N,dvar_vector& M,MY_DOUBLE_TYPE rmax);

class xxxtags
{
  int slave_number;
  //DATA
  ivector parest_flags;  
  int num_regions;
  int min_tag_group;
  int max_tag_group;
  int num_tag_releases;
  ivector tag_year;
  imatrix initial_tag_period;
  imatrix move_flags;
  imatrix move_index;
  ivector num_fish_periods;
  imatrix num_fish_incidents;
  imatrix tag_flags;
  imatrix year;
  ivector tag_region; 
  int nage;
  i3_array     parent; 

  imatrix      realization_period;  // The fishing period in which each realiz
  imatrix      realization_region;
  imatrix      realization_incident;

  ivector      age_flags;
  d4_array tot_tag_catch; 
  imatrix      fish_flags;  
  int		num_fisheries;
  dvar5_array   obstagcatch;
  imatrix      terminal_tag_period; 
  ivector      num_real_fish_periods;
  i3_array fishery_projection_flag;
  i3_array grouped_fishery_projection_flag;

  ivector min_init_tag_period;
  ivector      minttp;

  int          tag_shlen;    
  MY_DOUBLE_TYPE       tag_filen; 
  int          tag_nlint; 
  dmatrix      initial_tag_release_by_length;

  int pq_flag;

  imatrix      num_pooledtagfish_incidents;

  int     min_tag_year;
  int     nyears;
  d4_array     pooledobstagcatch_by_length;    
  i4_array min_tag_age4;
  i3_array     num_tagfish_incidents;
  // VARIABLES
  dvar3_array rep_rate;
  dvar_matrix nat_mort; 
  dvar_matrix fraction;
  dvar4_array incident_sel;
  dvar3_array tot_mort; 
  dvar4_array tagnum_fish;
  dvar4_array tagN;
  dvar4_array  nrtm;
  dvar4_array Dad; 
  dvar5_array  nrfm;
  dvar4_array  nrsurv;
  dvar5_array  tagcatch; 
  dvar4_array  fish_mort_calcs; 
  dvar_matrix  initial_tag_release_by_age;
  dvar_matrix      fish_pars;
  dvar_matrix      tag_fish_rep;
  dvar4_array     pooled_tagcatch; 
  dvar4_array     pooledobstagcatch;

  dvar_vector vb_coff;
  dvar_vector  gml;

  dvar3_array      epooled_tagnum_fish_recr;
  dvar3_array      pooledtagN; 

  dvar3_array  pooled_tagnum_fish;
  dvar5_array obstagcatch1;
  dvar_vector var_coff;
  d5_array     obstagcatch_by_length;  
  d3_array     pooledtot_tag_catch;   

  ivector initial_tag_year;
  i3_array min_tag_age1;
  ivector  terminal_tag_year; 
  i3_array min_tag_age5;
  i3_array min_tag_age6;
  i3_array min_tag_age;



public:
  pxxxmspd pmsd;
  void allocate(void);
  xxxtags(void);
 /*
  xxxtags(int _sno, int _num_regions, int _min_tag_group,int __max_tag_group, int _num_tag_releases, ivector & _tag_year, imatrix & _initial_tag_period, imatrix & _move_flags, imatrix & _move_index, ivector _num_fish_periods, imatrix & _num_fish_incidents, imatrix & _tag_flags, imatrix & _year, ivector _tag_region, int _nage, i3_array &  _parent, imatrix &_realization_period,imatrix &_realization_region, imatrix &_realization_incident, ivector & _age_flags, d4_array& _tot_tag_catch,imatrix & fish_flags,
   ivector & _num_real_fish_periods,
   imatrix & _terminal_tag_period,
   ivector & _parest_flags,
   ivector min_init_tag_period,
   ivector minttp,
   imatrix      num_pooledtagfish_incidents,
   int          tag_shlen,
   int          min_tag_year,
   int          nyears,
   MY_DOUBLE_TYPE       tag_filen,
   int          tag_nlint,
   imatrix      initial_tag_release_by_length,
   int          num_fisheries,
   i3_array&    num_tagfish_incidents,
   d5_array&    obstagcatch_by_length,
   d4_array&    pooledobstagcatch_by_length,    
   d3_array&    pooledtot_tag_catch,
   ivector& initial_tag_year,
   i3_array& min_tag_age1,
   ivector&  terminal_tag_year,
   i4_array& min_tag_age4,
   i3_array& min_tag_age5,
   i3_array& min_tag_age6,
   i3_array& min_tag_age,
   i3_array grouped_fishery_projection_flag);
  */

  dvariable tag_catch_equations_calc(dvar_vector& sv);
  dvariable tag_catch_equations_calc_loop(void);
  dvariable tag_catch_equations_calc_loop_1(void);
  dvariable tag_catch_equations_calc_loop_2(void);
  void get_initial_tag_population(dvar_vector& sv, int it);
  void do_newton_raphson_for_tags2(int it,int ir,int ip,dvariable& _ffpen);
  void do_newton_raphson_for_tags(int it,int ir,int ip,dvariable& ffpen);
  void do_the_diffusion(int year,dvar_vector& sv,dvar3_array& N,
    ofstream * pofs=0);
  dvar_matrix get_nat_mort_region(int ir);
  int get_variables_from_master1(void);
  int get_variables_from_master_1A(int nt);
  int get_variables_from_master3(void);
  int get_variables_from_master2(int tn);
  dvariable fit_tag_returns(void);
  dvariable tag_rep_rate(int it,int ir,int rp,int fi);
  dvariable grouped_tag_reporting_rate_penalty(void);
  dvariable fit_pooled_tag_returns(void);
  void var_convert_tag_lengths_to_age(void);
  dvariable xtag_catch_equations_calc_pooled(dvar_vector& sv);
  void xpooled_tag_catch_equations_calc(dvar_vector& sv);
  void thread_mortality_calcs(dvar_vector& sv,ivector* pq_flag,
    const prevariable& fzpen);
  int send_variables_to_tag_slave1(int flag,int slave_number);
  int send_variables_to_tag_slave_1A(int flag,int slave_number);
  int send_variables_to_tag_slave2(int flag,int slave_number);
  void observed_tags_by_age_from_length_pooled(void);
  void observed_tags_by_age_from_length(void);
  dvar_vector get_vb_coff_region(int cs);
  dvar_vector get_var_coff_species(int is);
  void tot_tags_catch(void);
  void get_constant_data(int nt);
};
class x2mspd
{
public:
    //i4_array fisc; 
    i3_array fisn; 
    //i4_array sp_in_catch;
    i4_array reg_in_catch;
    //dvar_matrix rec_delta;
    //dvar_vector recmean;
    //dvar_vector totpop_coff;
};


class  xxtotcatwt
{
public:
  // constant
  imatrix      fish_flags; 
  ivector      num_real_fish_periods;
  int          num_regions; 
  ivector      parest_flags; 
  d3_array     obs_tot_catch;
  imatrix      num_fish_incidents; 
  ivector      age_flags; 
  dvector obstotalcatch_by_weight;
  ivector numtotalcatch_by_weight;
  int nage;
  dvector obstotalcatch_by_numbers;
  ivector numtotalcatch_by_numbers;
  i3_array     parent;  
  imatrix      data_fish_flags;     
  MY_DOUBLE_TYPE len_wt_coff;
  int avg_calc_flag;
  ivector      num_fish_periods; 

  // variables
  dvar3_array      tot_catch;
  dvar_vector totalcatch_by_weight;
  dvar3_array     pred_totcatch; 
  dvar4_array     catch;   
  dvar4_array mean_length;
  dvar4_array vars;
  dvar_vector totalcatch_by_numbers;
  dvar4_array     exp_catch;    
  dvar_vector      q0;  
  x2mspd * pmsd;
  void pmsd_error();
  dvariable get_sv_region(int is,int svind);
  dvar_vector sv;
  void allocate(void);
  void get_constant_data_3(int tn);
  void catch_wt_loop_3(void);
  //int get_variables_from_master_catch_wt2(int);
  int xxtotcatwt::get_variables_from_master3(int);
};

void admb_catch_wt(void *ptr);
