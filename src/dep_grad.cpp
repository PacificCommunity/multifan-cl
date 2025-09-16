/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
#include "here.h"
  dvar_vector sbcalc(dvar_len_fish_stock_history& fsh);

  dvariable dep_gradients_calc(dvar_len_fish_stock_history& fsh,
    int gradient_switch)
  {
    cout << "In dep_gradients_calc " << endl;
    dvariable f=0.;
    switch(gradient_switch)
    {
      case(1):
      {
    cout << "In dep_gradients_calc " << endl;
	dvar_vector rbio=vbiocalc(fsh);
	f=log(rbio(fsh.nyears)/max(rbio));
	break;
      }
      case(2):
      {
    cout << "In dep_gradients_calc " << endl;
	dvar_vector rbio=vbiocalc(fsh);
	f=log(-log(0.99*rbio(fsh.nyears)/max(rbio)));
	break;
      }
      case(3):
      {
    cout << "In dep_gradients_calc " << endl;
	dvar_vector rbio=vbiocalc(fsh);
	f=log(rbio(fsh.nyears)/mean(rbio));
	break;
      }
      case(4):
      {
    cout << "In dep_gradients_calc " << endl;
	dvar_vector rbio=sbcalc(fsh);
	if (fsh.age_flags(95)>0)
	  f=log(rbio(fsh.nyears)/mean(rbio(1,fsh.age_flags(95))));
        else
	  f=log(rbio(fsh.nyears)/rbio(1));
	break;
      }
      default:
      {
	cerr << "Error illegal switch value in dep_gradients_calc"
	     << endl;
      }
    }
    return f;
  }

  dvariable biomass_gradients_calc(dvar_len_fish_stock_history& fsh,
    int gradient_switch)
  {
    dvariable f=0.;
    dvar_vector rbio=vbiocalc(fsh);
    f=log(rbio(gradient_switch-3));
    return f;
  }

  dvariable recruitment_gradients_calc(dvar_len_fish_stock_history& fsh,
    int gradient_switch)
  {
    dvariable f=0.;
    int this_year=gradient_switch-fsh.nyears-3;
    MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
    // get the initial numbers of fish at age
    int yr=fsh.nyears; 
    if (fsh.age_flags(41)==0)
    { 
      cerr << "error  -- can't do standard devs of recruitment for"
              "age_flags(41)=0"<< endl;
      f=0.;
      return f;
    }
    else
    {
      f=fsh.recr(this_year)-fsh.recmean+fsh.rec_init_diff+fsh.totpop;
      return f;
    }
  }
