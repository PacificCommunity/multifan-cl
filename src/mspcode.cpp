/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

//% message to emacs: -*- mode: C++; fill-column: 80 -*-


#include "all.hpp"

dvar_vector  dvar_len_fish_stock_history::get_global_vars_region(int ir)
{
  if (!pmsd)
    return global_vars;
  else
  {
    int is=pmsd->region_species_pointer(ir);
    if (is==1)
      return global_vars;
    else
      return pmsd->global_vars(is);
  }
}
dvar_vector  dvar_len_fish_stock_history::get_vb_coff_region(int ir)
{
  if (!pmsd)
    return vb_coff;
  else
  {
    int is=pmsd->region_species_pointer(ir);
    if (is==1)
      return vb_coff;
    else
      return pmsd->vb_coff(is);
  }
}

dvar_vector  dvar_len_fish_stock_history::get_vb_coff(int is)
{
  if (!pmsd)
    return vb_coff;
  else
  {
    if (is==1)
      return vb_coff;
    else
      return pmsd->vb_coff(is);
  }
}

dvar_vector  dvar_len_fish_stock_history::get_global_vars(int cs)
{
  if (cs==1)
    return global_vars;
  else
    return pmsd->global_vars(cs);
}

int  dvar_fish_stock_history::get_nage_region(int ir)
{
  if (!pmsd)
    return nage;
  else
  {
    int is=pmsd->region_species_pointer(ir);
    if (is==1)
      return nage;
    else
      return pmsd->nage(is);
  }
}
int  dvar_fish_stock_history::get_nage(int is)
{
  if (!pmsd)
    return nage;
  else
  {
    if (is==1)
      return nage;
    else
      return pmsd->nage(is);
  }
}

MY_DOUBLE_TYPE  dvar_len_fish_stock_history::get_len_wt_coff_region(int ir)
{
  if (!pmsd)
    return len_wt_coff;
  else
  {
    int is=pmsd->region_species_pointer(ir);
    if (is==1)
    {
      return len_wt_coff;
    }
    else
    {
      return pmsd->len_wt_coff(is);
    }
  }
}

int dvar_fish_stock_history::get_current_species(void)
{
  if (!pmsd)
    return 1;
  else
    return pmsd->current_species;
}

int dvar_fish_stock_history::get_current_nage(void)
{
  if (!pmsd || pmsd->current_species==1)
    return nage;
  else
    return pmsd->nage(pmsd->current_species);
}

dvariable dvar_fish_stock_history::get_sv_region(int ir,int svind)
{
  if (!pmsd)
    return sv(svind);
  else
  {
    int is=pmsd->region_species_pointer(ir);
    if (is==1)
      return sv(svind);
    else
      return pmsd->sv(is,svind);
  }
}
    
dvariable dvar_fish_stock_history::get_sv_species(int is,int svind)
{
  if (!pmsd)
    return sv(svind);
  else
  {
    if (is==1)
      return sv(svind);
    else
      return pmsd->sv(is,svind);
  }
}
    
dmatrix dvar_len_fish_stock_history::get_multi_wmid(void)
{
  dmatrix tmp(2,pmsd->num_species);

  for (int is=2;is<=pmsd->num_species;is++)
  {
    tmp(is)=pow(realwmid/value(get_sv_species(is,27)),1.0/value(get_sv_species(is,28)));
  }
  return tmp;
}

    
dmatrix dvar_len_fish_stock_history::get_multi_wm2(void)
{
  dmatrix tmp(2,pmsd->num_species);

  for (int is=2;is<=pmsd->num_species;is++)
  {
#if !defined(NO_MY_DOUBLE_TYPE)
    tmp(is)=pow(pmsd->wmid(is),value(get_sv_species(is,28)-1.0L));
#else
    tmp(is)=pow(pmsd->wmid(is),value(get_sv_species(is,28)-1.0));
#endif
  }
  return tmp;
}

extern int _NUMSV;

multi_species_data::multi_species_data(int _num_species,int _num_regions): 
  num_species(_num_species),num_real_regions(_num_regions),
  totpop_coff(2,_num_species),
  //totpop(2,_num_species),
  recr(2,_num_species),
  recr1(2,_num_species),
  recmean(2,_num_species),
  nat_mort_coff(2,_num_species),   //NMD 02Nov2011
  parest_flags(2,_num_species,1,400),
  age_flags(2,_num_species,1,200),
  historical_age_flags(2,_num_species,1,200), //NMD_23Apr2015
  historical_parest_flags(2,_num_species,1,400), //NMD_23Apr2015
  nage(2,_num_species),
  num_tag_groups(1,_num_species),
  sv(2,_num_species,1,_NUMSV),
  rowoffset(1,_num_species),
  region_bounds(1,_num_species,1,2),
  species_by_region(1,_num_species,1,_num_regions),
  vb_coff(2,_num_species,1,4),
  fmin1(2,_num_species),
  fmax1(2,_num_species),
  fminl(2,_num_species),
  fmaxl(2,_num_species),
  rhomin(2,_num_species),
  rhomax(2,_num_species),
  var_coff(2,_num_species,1,2),
  vmin1(2,_num_species),
  vmax1(2,_num_species),
  vminl(2,_num_species),
  vmaxl(2,_num_species),
  len_wt_coff(2,_num_species),
  global_vars(2,_num_species),
  alpha(2,_num_species),
  beta(2,_num_species),
  bh_predicted_recruits(2,num_species),  //NMD_12Apr2021
  steepness(2,_num_species),
  species_flags(1,_num_species,1,10),
  combined_tags_flag(0),
  bh_recr_devs(2,_num_species),
  phi(2,_num_species),
  recr_degree_yr(2,_num_species),
  recr_degree_reg(2,_num_species),
  recr_degree_ses(2,_num_species), 
  OR(2,_num_species),
  num_new_weights(2,_num_species),
  new_orth_recr(2,_num_species),
  numcomp(2,_num_species,1,4), 
  orth_recr_basis(2,_num_species),
  degree_yr(2,_num_species),
  degree_reg(2,_num_species),
  degree_ses(2,_num_species),
  degree_ses_reg(2,_num_species),
  recr_polys_yr(2,_num_species),
  recr_polys_reg(2,_num_species),
  recr_polys_ses(2,_num_species),
  recr_polys_ses_reg(2,_num_species),
  ny_begin_yr(2,_num_species),
  ny_begin_reg(2,_num_species),
  ny_begin_ses(2,_num_species),
  ny_begin_ses_reg(2,_num_species),
  nd_begin_yr(2,_num_species),
  nd_begin_reg(2,_num_species),
  nd_begin_ses(2,_num_species),
  nd_begin_ses_reg(2,_num_species),
  ny_end_yr(2,_num_species),
  ny_end_reg(2,_num_species),
  ny_end_ses(2,_num_species),
  ny_end_ses_reg(2,_num_species),
  nd_end_yr(2,_num_species),
  nd_end_reg(2,_num_species),
  nd_end_ses(2,_num_species),
  nd_end_ses_reg(2,_num_species),
  yearly_recr_all(2,_num_species),
  orth_recr_all(2,_num_species)
{ 
  orth_recr_all.allocate(2,num_species);
  numcomp.initialize();
  num_new_weights.initialize();
  parest_flags.initialize();
}



multi_species_data::multi_species_data(void) { }

void multi_species_data::allocate(int _num_species,int _num_regions)
{
  num_species=_num_species;
  num_real_regions=_num_regions;
  num_tag_groups.allocate(1,num_species);
  species_by_region.allocate(1,num_species,1,num_real_regions);
}

dvar_vector dvar_fish_stock_history::get_nat_mort_species(void)
{
  if (!pmsd)
    return nat_mort(nyears);
  else
  {
    int is=pmsd->current_species;
    if (is==1)
      return nat_mort(nyears);
    else
      return pmsd->nat_mort(is,nyears);
  }
}

ivector dvar_fish_stock_history::get_region_bounds(void)
{
  if (!pmsd)
  {
    ivector tmp(1,2);
    tmp(1)=1;
    tmp(2)=num_regions;
    return tmp;
  }
  else
  {
    return pmsd->region_bounds(pmsd->current_species);
  }
}

ivector dvar_fish_stock_history::get_region_bounds(int sp)
{
  if (!pmsd)
  {
    ivector tmp(1,2);
    tmp(1)=1;
    tmp(2)=num_regions;
    return tmp;
  }
  else
  {
    return pmsd->region_bounds(sp);
  }
}

dvar_vector dvar_len_fish_stock_history::get_var_coff_species(int is)
{
  if (!pmsd)
    return var_coff;
  else
  {
    if (is==1)
      return var_coff;
    else
      return pmsd->var_coff(is);
  }
}
int dvar_fish_stock_history::get_nage_tag(int it)
{
  if (!pmsd)
    return nage;
  else
  {
    int is=pmsd->tag_species_pointer(it);
    if (is==1)
      return nage;
    else
      return pmsd->nage(is);
  }
}
int dvar_fish_stock_history::get_nage_species(int is)
{
  if (!pmsd)
    return nage;
  else
  {
    if (is==1)
      return nage;
    else
      return pmsd->nage(is);
  }
}
dvar_matrix dvar_fish_stock_history::get_nat_mort_species(int is)
{
  if (!pmsd)
    return nat_mort;
  else
  {
    if (is==1)
      return nat_mort;
    else
      return pmsd->nat_mort(is);
  }
}

dvar_vector dvar_fish_stock_history::get_age_pars_species(int is,int i)
{
  if (!pmsd)
    return age_pars(i);
  else
  {
    if (is==1)
      return age_pars(i);
    else
      return pmsd->age_pars(is,i);
  }
}

dvariable dvar_fish_stock_history::get_nat_mort_coff_region(int ir)
{
  if (!pmsd)
    return nat_mort_coff;
  else
  {
    int is=pmsd->region_species_pointer(ir);
    if (is==1)
      return nat_mort_coff;
    else
      return pmsd->nat_mort_coff(is);
  }
}

int dvar_fish_stock_history::get_age_flags_region(int ir,int i)
{
  if (!pmsd)
    return age_flags(i);
  else
  {
    int is=pmsd->region_species_pointer(ir);
    if (is==1)
      return age_flags(i);
    else
      return pmsd->age_flags(is,i);
  }
}

dvar_matrix dvar_fish_stock_history::get_nat_mort_region(int ir)
{
  if (!pmsd)
    return nat_mort;
  else
  {
    int is=pmsd->region_species_pointer(ir);
    if (is==1)
      return nat_mort;
    else
      return pmsd->nat_mort(is);
  }
}

dvar_vector dvar_fish_stock_history::get_pmature_region(int ir)
{
  if (!pmsd)
    return pmature;
  else
  {
    int is=pmsd->region_species_pointer(ir);
    if (is==1)
      return pmature;
    else
      return pmsd->pmature(is);
  }
}

dvar_vector dvar_fish_stock_history::get_pmature_species(int is)
{
  if (!pmsd)
    return pmature;
  else
  {
    if (is==1)
      return pmature;
    else
      return pmsd->pmature(is);
  }
}



