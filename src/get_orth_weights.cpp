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

  
void dvar_fish_stock_history::get_orthogonal_polynomial_weights(void)
{
  ivector pf;
  int * pny_begin_yr=0;
  int * pny_begin_reg=0;
  int * pny_begin_ses=0;
  int * pny_begin_ses_reg=0;
  int * pnd_begin_yr=0;
  int * pnd_begin_reg=0;
  int * pnd_begin_ses=0;
  int * pnd_begin_ses_reg=0;
  int * pny_end_yr=0;
  int * pny_end_reg=0;
  int * pny_end_ses=0;
  int * pny_end_ses_reg=0;
  int * pnd_end_yr=0;
  int * pnd_end_reg=0;
  int * pnd_end_ses=0;
  int * pnd_end_ses_reg=0;

  int ns=1;
  if (pmsd)
  {
    ns=pmsd->num_species;
  }
  int is=1;
  if (pmsd) is=pmsd->current_species;
  if (!pmsd || is==1)
  {
    pf=parest_flags;
    pny_begin_yr=&ny_begin_yr;
    pny_begin_reg=&ny_begin_reg;
    pny_begin_ses=&ny_begin_ses;
    pny_begin_ses_reg=&ny_begin_ses_reg;
    pnd_begin_yr=&nd_begin_yr;
    pnd_begin_reg=&nd_begin_reg;
    pnd_begin_ses=&nd_begin_ses;
    pnd_begin_ses_reg=&nd_begin_ses_reg;
    pny_end_yr=&ny_end_yr;
    pny_end_reg=&ny_end_reg;
    pny_end_ses=&ny_end_ses;
    pny_end_ses_reg=&ny_end_ses_reg;
    pnd_end_yr=&nd_end_yr;
    pnd_end_reg=&nd_end_reg;
    pnd_end_ses=&nd_end_ses;
    pnd_end_ses_reg=&nd_end_ses_reg;
  }
  else
  {
    pf=pmsd->parest_flags(is);
    pny_begin_yr=&pmsd->ny_begin_yr(is);
    pny_begin_reg=&pmsd->ny_begin_reg(is);
    pny_begin_ses=&pmsd->ny_begin_ses(is);
    pny_begin_ses_reg=&pmsd->ny_begin_ses_reg(is);
    pnd_begin_yr=&pmsd->nd_begin_yr(is);
    pnd_begin_reg=&pmsd->nd_begin_reg(is);
    pnd_begin_ses=&pmsd->nd_begin_ses(is);
    pnd_begin_ses_reg=&pmsd->nd_begin_ses_reg(is);
    pny_end_yr=&pmsd->ny_end_yr(is);
    pny_end_reg=&pmsd->ny_end_reg(is);
    pny_end_ses=&pmsd->ny_end_ses(is);
    pny_end_ses_reg=&pmsd->ny_end_ses_reg(is);
    pnd_end_yr=&pmsd->nd_end_yr(is);
    pnd_end_reg=&pmsd->nd_end_reg(is);
    pnd_end_ses=&pmsd->nd_end_ses(is);
    pnd_end_ses_reg=&pmsd->nd_end_ses_reg(is);
  }

  int dy,dr,ds,dsr;
  ivector nc;
  if (is==1)
  {
    nc=numcomp;
    dy=degree_yr;
    dr=degree_reg;
    ds=degree_ses;
    dsr=degree_ses_reg;
  }
  else
  {
    nc=pmsd->numcomp(is);
    dy=pmsd->degree_yr(is);
    dr=pmsd->degree_reg(is);
    ds=pmsd->degree_ses(is);
    dsr=pmsd->degree_ses_reg(is);
  }

  int nry=num_real_years;
//  int nry=(nyears-1)/num_seasons+1;
 
  *pny_begin_yr=pf(183);
  // pf(183) is the default for the number of begining years with
  // orthogonal polynomials
  *pny_begin_reg=pf(183);
  // change the default begi*pning year for regions
  if (pf(204)>0)
  {
    *pny_begin_reg=pf(204);
  }
  if (pf(204)==-1)
  {
    // turn off begi*pning year for regions
    *pny_begin_reg=0;
  }
  *pny_begin_ses=pf(183);
  if (pf(206)>0)
  {
    *pny_begin_ses=pf(206);
  }
  if (pf(206)==-1)
  {
    // turn off begi*pning year for regions
    *pny_begin_ses=0;
  }
  *pny_begin_ses_reg=pf(183);
  if (pf(208)>0)
  {
    *pny_begin_ses_reg=pf(208);
  }
  if (pf(208)==-1)
  {
    // turn off begi*pning year for regions
    *pny_begin_ses_reg=0;
  }
 
  *pnd_begin_yr=pf(201);
  // pf(201) is the default for the  begining lower degree
  // of the orthogonal polynomials
  *pnd_begin_reg=pf(201);
  // change the default begi*pning degree for regions
  if (pf(205)>0)
  {
    *pnd_begin_reg=pf(205);
  }
  if (pf(205)==-1)
  {
    // turn off begi*pning year for regions
    *pnd_begin_reg=0;
  }
  *pnd_begin_ses=pf(201);
  // change the default begi*pning degree for seasons
  if (pf(207)>0)
  {
    *pnd_begin_ses=pf(207);
  }
  if (pf(207)==-1)
  {
    // turn off begi*pning year for regions
    *pnd_begin_ses=0;
  }
  *pnd_begin_ses_reg=pf(201);
  // change the default begi*pning degree for seasons
  if (pf(209)>0)
  {
    *pnd_begin_ses_reg=pf(209);
  }
  if (pf(209)==-1)
  {
    // turn off begi*pning year for regions
    *pnd_begin_ses_reg=0;
  }
 
  *pny_end_yr=pf(202);
  // pf(202) is the default for the number of ending years with
  // lower degree orthogonal polynomials
 
  *pny_end_reg=pf(202);
  // change the default for the number of ending years with
  // lower degree orthogonal polynomials for regions
  if (pf(210)>0)
  {
    *pny_end_reg=pf(210);
  }
  if (pf(210)==-1)
  {
    // turn off ending year for regions
    *pny_end_reg=0;
  }
  *pny_end_ses=pf(202);
  // change the default for the number of ending years with
  // lower degree orthogonal polynomials for seasons
  if (pf(212)>0)
  {
    *pny_end_ses=pf(212);
  }
  if (pf(212)==-1)
  {
    // turn off ending year for regions
    *pny_end_ses=0;
  }
  *pny_end_ses_reg=pf(202);
  // change the default for the number of ending years with
  // lower degree orthogonal polynomials for season-regions
  if (pf(214)>0)
  {
    *pny_end_ses_reg=pf(214);
  }
  if (pf(214)==-1)
  {
    // turn off ending year for regions
    *pny_end_ses_reg=0;
  }
 
  *pnd_end_yr=pf(203);
  // pf(203) is the default for the  ending lower degree
  // of the orthogonal polynomials
  *pnd_end_reg=pf(203);
  // change the default ending degree for regions
  if (pf(211)>0)
  {
    *pnd_end_reg=pf(211);
  }
  if (pf(211)==-1)
  {
    // turn off begi*pning year for regions
    *pnd_end_reg=0;
  }
  *pnd_end_ses=pf(203);
  // change the default ending degree for seasons
  if (pf(213)>0)
  {
    *pnd_end_ses=pf(213);
  }
  if (pf(213)==-1)
  {
    // turn off begi*pning year for regions
    *pnd_end_ses=0;
  }
  *pnd_end_ses_reg=pf(203);
  // change the default ending degree for region-seasons
  if (pf(215)>0)
  {
    *pnd_end_ses_reg=pf(215);
  }
  if (pf(215)==-1)
  {
    // turn off begi*pning year for regions
    *pnd_end_ses_reg=0;
  }
  if (ny_begin_yr==0)
    *pny_begin_yr=1;
  if (ny_end_yr==0)
    *pny_end_yr=1;
    //ny_end_yr=num_real_years;
 
  if (ny_begin_reg==0)
    *pny_begin_reg=1;
  if (ny_end_reg==0)
    *pny_end_reg=1;
    //ny_end_reg=num_real_years;
 
  if (ny_begin_ses==0)
    *pny_begin_ses=1;
  if (ny_end_ses==0)
    *pny_end_ses=1;
    //ny_end_ses=num_real_years;
 
  if (ny_begin_ses_reg==0)
    *pny_begin_ses_reg=1;
  if (ny_end_ses_reg==0)
    *pny_end_ses_reg=1;
    //ny_end_ses_reg=num_real_years;
 
 // note tht the number of parameters is one more than the
 // degree of the orthogonal polynomial so we subtract one
  if (is==1)
  {
    if (dy>0)
    { 
      recr_polys_yr=
        orthpoly_constant_begin_end(nry,dy-1,*pnd_begin_yr,
        *pny_begin_yr,*pnd_end_yr,*pny_end_yr);
    }
         
    if (dr>0)
    { 
      recr_polys_reg=
        orthpoly_constant_begin_end(nry,dr-1,*pnd_begin_reg,
        *pny_begin_reg,*pnd_end_reg,*pny_end_reg);
    }
         
    if (ds>0)
    { 
      recr_polys_ses=
        orthpoly_constant_begin_end(nry,ds-1,*pnd_begin_ses,
        *pny_begin_ses,*pnd_end_ses,*pny_end_ses);
    }
   
    if (dsr>0)
    { 
      recr_polys_ses_reg=
        orthpoly_constant_begin_end(nry,dsr-1,*pnd_begin_ses_reg,
        *pny_begin_ses_reg,*pnd_end_ses_reg,*pny_end_ses_reg);
    }
  }
  else
  {
    if (dy>0)
    { 
      pmsd->recr_polys_yr(is)=
        orthpoly_constant_begin_end(nry,dy-1,*pnd_begin_yr,
        *pny_begin_yr,*pnd_end_yr,*pny_end_yr);
    }
         
    if (dr>0)
    { 
      pmsd->recr_polys_reg(is)=
        orthpoly_constant_begin_end(nry,dr-1,*pnd_begin_reg,
        *pny_begin_reg,*pnd_end_reg,*pny_end_reg);
    }     

    if (ds>0)
    { 
      pmsd->recr_polys_ses(is)=
        orthpoly_constant_begin_end(nry,ds-1,*pnd_begin_ses,
        *pny_begin_ses,*pnd_end_ses,*pny_end_ses);
    }
   
    if (dsr>0)
    { 
      pmsd->recr_polys_ses_reg(is)=
        orthpoly_constant_begin_end(nry,dsr-1,*pnd_begin_ses_reg,
        *pny_begin_ses_reg,*pnd_end_ses_reg,*pny_end_ses_reg);
    }
  }
}
