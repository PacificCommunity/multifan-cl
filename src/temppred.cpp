/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#define HOME_VERSION
#include <fvar.hpp>

//#ifdef __ZTC__
//void DF_prob_calc(fpos_t filepos);
//#else
void DF_prob_calc(void);
//#endif

dvariable prob_calc(MY_DOUBLE_TYPE& fmidd,prevariable& mean_len,prevariable& sigg)
{
  dvariable out;
  MY_DOUBLE_TYPE tmp=(fmidd-value(mean_len))/value(sigg);
  MY_DOUBLE_TYPE prob=exp(-tmp*tmp/2.)/value(sigg);
  value(out)=prob;
//    save_identifier_string("ACE2");
  const char * str1;
  str1="ACE2";
  char* strx1=const_cast <char*> (str1);
  save_identifier_string(strx1);
  save_double_value(tmp);
  sigg.save_prevariable_value();
  mean_len.save_prevariable_position();
  sigg.save_prevariable_position();
  out.save_prevariable_position();
  out.save_prevariable_value();
//    save_identifier_string("ACE1");
  const char * str2;
  str2="ACE1";
  char* strx2=const_cast <char*> (str2);
  save_identifier_string(strx2);
  fpos_t filepos;
  //#ifdef __ZTC__
  //  gradient_structure::GRAD_STACK1->set_gradient_stack(DF_prob_calc,filepos);
 // #else
    gradient_structure::GRAD_STACK1->set_gradient_stack(DF_prob_calc);
  //#endif
  return out;
}

//#ifdef __ZTC__
//void DF_prob_calc(fpos_t filepos)
//#else
void DF_prob_calc(void)
//#endif
{
  verify_identifier_string("ACE1");
  MY_DOUBLE_TYPE prob=restore_prevariable_value();
  prevariable_position prob_pos=restore_prevariable_position();
  prevariable_position sigg_pos=restore_prevariable_position();
  prevariable_position mean_len_pos=restore_prevariable_position();
  MY_DOUBLE_TYPE sigg=restore_prevariable_value();
  MY_DOUBLE_TYPE tmp=restore_double_value();
  MY_DOUBLE_TYPE dfprob=restore_prevariable_derivative(prob_pos);
  verify_identifier_string("ACE2");
  // prob=exp(-tmp*tmp/2.)/sigg;
  MY_DOUBLE_TYPE dftmp=-dfprob*prob*tmp;
  MY_DOUBLE_TYPE dfsigg=-dfprob*prob/sigg;
  // MY_DOUBLE_TYPE tmp=(value(fmidd)-value(mean_len))/value(sigg);
  MY_DOUBLE_TYPE dfmean_len=dftmp/sigg;
  dfsigg-=dftmp*tmp/sigg;
  save_double_derivative(dfsigg,sigg_pos);
  save_double_derivative(dfmean_len,mean_len_pos);
}

