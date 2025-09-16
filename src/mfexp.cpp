/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


#include <fvar.hpp>

  dvariable mfexp(_CONST prevariable& x)
  {
    MY_DOUBLE_TYPE b=60;
    if (x<b) 
    {
      return exp(x);
    }
    else
    {
      dvariable y=x-b;
      return exp(b)*(1.+2.*y)/(1.+y);
    }
  }
  dvariable mfexp(prevariable& x)
  {
    MY_DOUBLE_TYPE b=60;
    if (x<b) 
    {
      return exp(x);
    }
    else
    {
      dvariable y=x-b;
      return exp(b)*(1.+2.*y)/(1.+y);
    }
  }

