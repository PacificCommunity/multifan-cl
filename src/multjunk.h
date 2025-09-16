/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
extern ofstream DERCH_ofs;
extern ofstream FMM_ofs;
extern dvector mean_weights_kludge;
extern adstring * _current_path;
//extern int _sd_flag;
extern int _version_no;
extern int _file_version_no;
adstring delayed_infile;
adstring delayed_outfile;

#include "variable.hpp"
#ifdef __i860
  #include "timer.h"
  event_timer fmin_timer("fmc.fmin(f,x,g)");
  event_timer fcomp_timer("fcomp(fsh,x,nvar,0)");
  event_timer gradcalc_timer("gradcalc(nvar,g)");
  ofstream event_timer::report_file;

#endif
dvar_len_fish_stock_history * pcfsh;


int check_mean_length_bounds(dvar_len_fish_stock_history& fsh);


    void fitting_procedure(dvar_len_fish_stock_history& fsh,int& quit_flag);
    void master_fitting_procedure(dvar_len_fish_stock_history& fsh,int& quit_flag);
    dvariable objective_function(dvar_fish_stock_history& fsh);
    MY_DOUBLE_TYPE fcomp(dvar_len_fish_stock_history& fsh,dvar_vector& x,
      int print_switch);
    void write_report(ostream& ofs,dvar_fish_stock_history& fsh);
    void normalize_year(fishery_freq_record_array& fra);
#include <float.h>

#ifdef HERE
  #undef HERE
#endif

