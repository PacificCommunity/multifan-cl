/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#if !defined(__OLD_SET__)
#define __OLD_SET__
void set_value(const dvar_vector& x,_CONST dvar_vector& v, const int& _ii,
   double fmin,double fmax,const double s);
#endif 
void set_value(const dvar3_array& _w,const dvar_vector& x,const int& ii,
  double fmin,double fmax,const dvariable& pen,double s,ivector&  flags,
  const ivector& group);
#include "allfiles"
