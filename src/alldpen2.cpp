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
 
dvariable avail_coff_penalty(dvar_fish_stock_history& fsh) 
{
  dvariable pen=0.0;
  for (int i=1;i<=fsh.age_flags(54)-1;i++)
  {
    pen+=.1*square(fsh.avail_coff(i));
  }
  return pen;
}

