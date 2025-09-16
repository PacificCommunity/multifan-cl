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

//extern adstring * _current_path;

  void getpaths(adstring& current_path,adstring& root);

void dep_grad_routine(independent_variables& x,dvar_len_fish_stock_history& fsh,
  int nvar,int grad_switch)
{
  MY_DOUBLE_TYPE f;
  adstring current_path;
  adstring root;
  root=ad_root;
  adstring fname=root + adstring(".dep");
  adstring fname2=root + adstring(".dp2");
  cout <<"gradrout.cpp " << fname << endl;
  uostream ofs(fname);
  if (!ofs)
  {
    cerr << "Error trying to open file " << fname << endl;
    exit(1);
  }
  HERE
  dvector g(1,nvar);
  gradient_structure::set_YES_DERIVATIVES();
  cout << " Calling fcomp for the dependent variables " << endl;
  f=fcomp(fsh,dvar_vector(x),nvar,0,g,grad_switch);
  cout << " Finished fcomp for the dependent variables " << endl;
  HERE
  //ofs << nvar << "  "  << int(f) << endl;
  ofs << nvar << int(f) ;
  if (fsh.parest_flags(229)>0)
  {
    ofs << fsh.parest_flags(229);
    ofs << fsh.parest_flags(230);
  }
  HERE
  if(fsh.age_flags(20))  //NMD_3May2018
  {
    if(fsh.parest_flags(237))
    {
      cout << " Calculating derivatives for the dependent variables " << endl;
      jacobcalc(nvar,ofs);
    }
  }
  else
  {
    cout << " Calculating derivatives for the dependent variables " << endl;
    jacobcalc(nvar,ofs);
  }     //NMD_3May2018
//  jacobcalc(nvar,ofs);
  HERE

  ivector q_flag=column(fsh.fish_flags,55);
  if (sum(q_flag))
  {
    uostream ofs2(fname2);
    cout << " Calling fcomp no effort for the dependent variables " << endl;
//
// apply same no-fishing settings as used in minimising_routine NMD_7Sep2015
    fsh.af170q0ex=1;
    fsh.af170q0=1;
    fsh.age_flags(170)=1;
    fsh.allocate_optional_stuff();
// NMD_7Sep2015
    fsh.zero_out_catches();  //NMD_7Apr2022
// zero out the effort from fisheries that are to be turned off;
    if (fsh.age_flags(92)==2 && sum(fsh.data_fish_flags(2)))
    {
      fsh.zero_out_log_effort();
    }
//    f=fcomp(fsh,dvar_vector(x),nvar,0,g,1,&q_flag);
    f=fcomp(fsh,dvar_vector(x),nvar,0,g,grad_switch,&q_flag);   //NMD 22Feb2012
    cout << " Finished fcomp no effort for the dependent variables " << endl;
    ofs2 << nvar << int(f) ;

    if(fsh.age_flags(20))  //NMD_3May2018
    {
      if(fsh.parest_flags(237))
      {
        cout << " Calculating derivatives for the dependent variables " << endl;
        jacobcalc(nvar,ofs2);
      }
    }
    else
    {
      cout << " Calculating derivatives for the dependent variables " << endl;
      jacobcalc(nvar,ofs2);
    }     //NMD_3May2018
//    jacobcalc(nvar,ofs2);
  }
  gradient_structure::set_NO_DERIVATIVES();
}

#undef HOME_VERSION

