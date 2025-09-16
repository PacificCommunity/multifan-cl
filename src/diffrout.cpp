/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
//#include <f:\linad99\fvar.hpp>
#include "all.hpp"


dvar_matrix fast_diffusion_calcs(int nage,int num_regions,dvar3_array& N,
  dvar3_array& Dad,ivector& rip,int maxage)
{ 
  const MY_DOUBLE_TYPE delta=1.e-6;
  const MY_DOUBLE_TYPE nd=num_regions*delta;
  int j;
  if (maxage) nage=min(nage,maxage);
  int jmin=N(1,rip(1)).indexmin();
  dvar_matrix EN(jmin,nage,1,num_regions);
  dvar_matrix tmp(jmin,nage,1,num_regions);
  tmp.initialize();

  for (j=jmin;j<=nage;j++)
  {
    for (int is=1;is<=num_regions;is++)
    {
      EN(j,is)=mfexp(N(is,rip(is),j));
    }
  }
  for (j=jmin;j<=nage;j++)
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int is=1;is<=num_regions;is++)
      {
        tmp(j,ir)+=Dad(j)(ir,is)*EN(j,is);
      }
    }
    dvariable ss=sum(tmp(j));
    tmp(j)+=delta;
    tmp(j)*=ss/(ss+nd);
  }
  return trans(tmp);
}


dvar_matrix fast_diffusion_calcs(int nage,int num_regions,dvar_matrix& N,
  dvar3_array& Dad,int maxage)
{ 
  const MY_DOUBLE_TYPE delta=1.e-6;
  const MY_DOUBLE_TYPE nd=num_regions*delta;
  int j;
  if (maxage) nage=min(nage,maxage);
  int jmin=N(1).indexmin();
  dvar_matrix EN(jmin,nage,1,num_regions);
  dvar_matrix tmp(jmin,nage,1,num_regions);
  tmp.initialize();

  for (j=jmin;j<=nage;j++)
  {
    for (int is=1;is<=num_regions;is++)
    {
      EN(j,is)=mfexp(N(is,j));
    }
  }
  for (j=jmin;j<=nage;j++)
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int is=1;is<=num_regions;is++)
      {
        tmp(j,ir)+=Dad(j)(ir,is)*EN(j,is);
      }
    }
    dvariable ss=sum(tmp(j));
    tmp(j)+=delta;
    tmp(j)*=ss/(ss+nd);
  }
  return trans(tmp);
}

dvar_vector fast_diffusion_calcs_for_plusgroup_biomass(int nage,
  int num_regions,dvar_matrix& B, dvar3_array& Dad,ivector& rip)
{ 
  const MY_DOUBLE_TYPE delta=1.e-6;
  const MY_DOUBLE_TYPE nd=num_regions*delta;
  int j;
  dvar_vector EN(1,num_regions);
  dvar_vector tmp(1,num_regions);
  tmp.initialize();

  for (int is=1;is<=num_regions;is++)
  {
    EN(is)=B(is,rip(is));
  }
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int is=1;is<=num_regions;is++)
    {
      tmp(ir)+=Dad(nage)(ir,is)*EN(is);
    }
  }
  return tmp;
}

dvar_matrix fast_diffusion_calcs(int nage,int num_regions,dvar3_array& N,
  dvar3_array& Dad,int year)
{ 
  const MY_DOUBLE_TYPE delta=1.e-6;
  const MY_DOUBLE_TYPE nd=num_regions*delta;
  dvar_matrix EN(1,nage,1,num_regions);
  dvar_matrix tmp(1,nage,1,num_regions);
  tmp.initialize();
  int j;
  for (j=1;j<=nage;j++)
  {
    for (int is=1;is<=num_regions;is++)
    {
      EN(j,is)=mfexp(N(is,year,j));
    }
  }
  for (j=1;j<=nage;j++)
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int is=1;is<=num_regions;is++)
      {
        tmp(j,ir)+=Dad(j)(ir,is)*EN(j,is);
      }
    }
    dvariable ss=sum(tmp(j));
    tmp(j)+=delta;
    tmp(j)*=ss/(ss+nd);
  }
  return trans(tmp);
}

dvar_matrix fast_diffusion_calcs(int nage,int num_regions,dvar3_array& N,
  dvar3_array& Dad,int year,pmulti_species_data & pmsd)
{ 
  const MY_DOUBLE_TYPE delta=1.e-6;
  //const MY_DOUBLE_TYPE nd=num_regions*delta;   //NMD 2Oct2012
  if (!pmsd)
  { 
    const MY_DOUBLE_TYPE nd=num_regions*delta;  //NMD 2Oct2012
    dvar_matrix EN(1,nage,1,num_regions);
    dvar_matrix tmp(1,nage,1,num_regions);
    tmp.initialize();
    int j;
    for (j=1;j<=nage;j++)
    {
      for (int is=1;is<=num_regions;is++)
      {
        EN(j,is)=mfexp(N(is,year,j));
      }
    }
    for (j=1;j<=nage;j++)
    {
      for (int ir=1;ir<=num_regions;ir++)
      {
        for (int is=1;is<=num_regions;is++)
        {
          tmp(j,ir)+=Dad(j)(ir,is)*EN(j,is);
        }
      }
      dvariable ss=sum(tmp(j));
      tmp(j)+=delta;
      tmp(j)*=ss/(ss+nd);
    }
    return trans(tmp);
  }
  else
  { 
    const MY_DOUBLE_TYPE nd=pmsd->num_real_regions*delta;  //NMD 2Oct2012
    int it=pmsd->tag_index;
    int nrr=pmsd->num_real_regions;
    int ns=pmsd->num_species;
    dvar_matrix tmp1(1,num_regions,1,pmsd->nage_by_region);
    int isp;
    int ng;
    for (isp=1;isp<=ns;isp++)
    {
      if (it==0 || isp==pmsd->tag_species_pointer(it))
      {
        if (isp==1)
          ng=nage;
        else
          ng=pmsd->nage(isp);
        int offset=(isp-1)*nrr;
        dvar_matrix EN(1,ng,1,nrr);
        dvar_matrix tmp(1,ng,1,nrr);
        tmp.initialize();
        int j;
        for (j=1;j<=ng;j++)
        {
          for (int is=1;is<=nrr;is++)
          {
            EN(j,is)=mfexp(N(is+offset,year,j));
          }
        }
        for (j=1;j<=ng;j++)
        {
          dvar_matrix td=Dad(j).sub(1+offset,nrr+offset).shift(1);
          for (int ir=1;ir<=nrr;ir++)
          {
            for (int is=1;is<=nrr;is++)
            {
              tmp(j,ir)+=td(ir,is)*EN(j,is);
            }
          }
          dvariable ss=sum(tmp(j));        
          tmp(j)+=delta;
          tmp(j)*=ss/(ss+nd);
        }
        tmp1.sub(1+offset,nrr+offset).shift(1)=trans(tmp);
      }
    }
    return tmp1;
  }
}


dvar_matrix fast_diffusion_calcs(int nage,int num_regions,dvar3_array& N,
  dvar3_array& Dad,ivector& rip,int maxage,pmulti_species_data & pmsd)
{ 
  const MY_DOUBLE_TYPE delta=1.e-6;
  //const MY_DOUBLE_TYPE nd=num_regions*delta;    //NMD 4Oct2012

  int j;
  if (!pmsd)
  { 
    const MY_DOUBLE_TYPE nd=num_regions*delta;  //NMD 2Oct2012
    if (maxage) nage=min(nage,maxage);
    int jmin=N(1,rip(1)).indexmin();
    dvar_matrix EN(jmin,nage,1,num_regions);
    dvar_matrix tmp(jmin,nage,1,num_regions);
    tmp.initialize();
  
    for (j=jmin;j<=nage;j++)
    {
      for (int is=1;is<=num_regions;is++)
      {
        EN(j,is)=mfexp(N(is,rip(is),j));
      }
    }
    for (j=jmin;j<=nage;j++)
    {
      for (int ir=1;ir<=num_regions;ir++)
      {
        for (int is=1;is<=num_regions;is++)
        {
          tmp(j,ir)+=Dad(j)(ir,is)*EN(j,is);
        }
      }
      dvariable ss=sum(tmp(j));
      tmp(j)+=delta;
      tmp(j)*=ss/(ss+nd);
    }
    return trans(tmp);
  }
  else
  { 
    const MY_DOUBLE_TYPE nd=pmsd->num_real_regions*delta;  //NMD 2Oct2012
    int it=pmsd->tag_index;
    if (maxage) nage=min(nage,maxage);
    int jmin=N(N.indexmin(),rip(rip.indexmin())).indexmin();
    int nrr=pmsd->num_real_regions;
    int ns=pmsd->num_species;
    dvar_matrix tmp1(1,num_regions,jmin,pmsd->nage_by_region);
    tmp1.initialize();
    int isp;
    int ng;
    for (isp=1;isp<=ns;isp++)
    {
      if (it==0 || isp==pmsd->tag_species_pointer(it))
      {
        if (isp==1)
          ng=nage;
        else
          ng=pmsd->nage(isp);
        int offset=(isp-1)*nrr;
        dvar_matrix EN(jmin,ng,1,nrr);
        dvar_matrix tmp(jmin,ng,1,nrr);
        tmp.initialize();
      
        for (j=jmin;j<=ng;j++)
        {
          for (int is=1;is<=nrr;is++)
          {
            EN(j,is)=mfexp(N(is+offset,rip(is+offset),j));
          }
        }
        for (j=jmin;j<=ng;j++)
        {
          dvar_matrix td=Dad(j).sub(1+offset,nrr+offset).shift(1);
          for (int ir=1;ir<=nrr;ir++)
          {
            for (int is=1;is<=nrr;is++)
            {
              tmp(j,ir)+=td(ir,is)*EN(j,is);
            }
          }
          dvariable ss=sum(tmp(j));
          tmp(j)+=delta;
          tmp(j)*=ss/(ss+nd);
        }
        tmp1.sub(1+offset,nrr+offset).shift(1)=trans(tmp);
      }
    }
    return tmp1;
  }
}

/*
dvar_matrix dvar_fish_stock_history::do_the_diffusion(dvar_matrix& N)
{ 
  const MY_DOUBLE_TYPE delta=1.e-4;
  const MY_DOUBLE_TYPE nd=num_regions*delta;
  dvar_matrix EN(1,nage,1,num_regions);
  dvar_matrix tmp(1,nage,1,num_regions);
  tmp.initialize();
  for (int j=1;j<=nage;j++)
  {
    for (int is=1;is<=num_regions;is++)
    {
      EN(j,is)=mfexp(N(is,j));
    }
  }
  for (j=1;j<=nage;j++)
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int is=1;is<=num_regions;is++)
      {
        tmp(j,ir)+=Dad(j)(ir,is)*EN(j,is);
      }
    }
    dvariable ss=sum(tmp(j));
    tmp(j)+=delta;
    tmp(j)*=ss/(ss+nd);
  }
  return trans(tmp);
}

*/

