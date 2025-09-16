/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#ifndef __MF_GLOBALS__
#define __MF_GLOBALS__

 const int MAX_DEPVAR_LABELS = 1000;
  class _mf_cl_globals
  {
  public:
    dvector dep_vars_values;
    adstring_array dep_var_labels;
    int max_depvar_labels;
    _mf_cl_globals();
  };
  
  extern _mf_cl_globals mfglobals;
#endif
