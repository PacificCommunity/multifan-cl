/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"
#ifdef __MSC__
  MY_DOUBLE_TYPE mfexp(const MY_DOUBLE_TYPE& );
  dvariable mfexp(const prevariable& );
  dvar_vector mfexp(const dvar_vector& );
  dvector mfexp(const dvector& );
#else
  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
  dvariable mfexp(prevariable& );
  //dvar_vector mfexp(dvar_vector& );
  //dvector mfexp(dvector& );
#endif
void ageflags_error(int n,const char * s)
{
  cerr << "Error yo need to set pf(" << n << ") to a number > 0" << endl;
  cerr << " see function " << s << endl;
  ad_exit(1);
}

void dvar_fish_stock_history::new_get_orthogonal_recruitment_poly_info(void)
{
  if (pmsd==0)
  {
    if (age_flags(137)==0)
    {
      ageflags_error(137,"new_get_orthogonal_recruitment_poly_info");
    }
    if (age_flags(138)==0)
    {
      ageflags_error(138,"new_get_orthogonal_recruitment_poly_info");
    }
    if (age_flags(139)==0)
    {
      ageflags_error(139,"new_get_orthogonal_recruitment_poly_info");
    }
    recr_degree_yr=age_flags(137)-1;
    recr_degree_reg=age_flags(138)-1;
    recr_degree_ses=age_flags(139)-1;
    OR=test_orthogonal();
  }
  else
  {
    if (!sum(column(pmsd->species_flags,2)))
    {
      for (int is=1;is<=pmsd->num_species;is++)
      {
        pmsd->current_species=is;
        if (pmsd->current_species==1)
        {
          recr_degree_yr=age_flags(137)-1;
          recr_degree_reg=age_flags(138)-1;
          recr_degree_ses=age_flags(139)-1;
          OR=test_orthogonal();
        }
        else
        {
          pmsd->recr_degree_yr(is)=pmsd->age_flags(is,137)-1;
          pmsd->recr_degree_reg(is)=pmsd->age_flags(is,138)-1;
          pmsd->recr_degree_ses(is)=pmsd->age_flags(is,139)-1;
          pmsd->OR(is)=test_orthogonal();
        }
      } 
    }
    else
    {
      ivector sf2=column(pmsd->species_flags,2);
      if (pmsd->num_species !=2 || sum(sf2)!=1 )
      {
        cerr << "only works at present for 2 species model" << endl;
        ad_exit(1);
      }
      for (int is=1;is<=pmsd->num_species;is++)
      {
        pmsd->current_species=is;
        if (is==1)
        {
          recr_degree_yr=age_flags(137)-1;
          recr_degree_reg=age_flags(138)-1;
          recr_degree_ses=age_flags(139)-1;
          OR=test_orthogonal();
        }
        else
        {
          pmsd->recr_degree_yr(is)=pmsd->age_flags(is,137)-1;
          pmsd->recr_degree_reg(is)=pmsd->age_flags(is,138)-1;
          pmsd->recr_degree_ses(is)=pmsd->age_flags(is,139)-1;
          pmsd->OR(is)=OR;
        }
      }
    }
  }
}

  
#undef HOME_VERSION


