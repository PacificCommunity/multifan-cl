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
  dvariable mfexp(BOR_CONST prevariable& );
  //dvar_vector mfexp(dvar_vector& );
  //dvector mfexp(dvector& );
#endif

#ifdef __MSVC32__
  dvariable age_at_length_calcxx(MY_DOUBLE_TYPE& v,dvar_vector& gml,int nslots);
#else
  dvariable age_at_length_calcxx(const dvariable& v,dvar_vector& gml,int nslots);
#endif

extern mf_pvm_manager * mf_pvm;

void dvar_len_fish_stock_history::observed_tags_by_age_from_length(void)
{
  obstagcatch1.initialize();
  dvar_vector vbc(1,vb_coff.indexmax());
  dvar_vector vc(1,var_coff.indexmax());
  vbc=vb_coff;
  vc=var_coff;
  dvar_vector sigma(1,nage);
  dvariable rho=exp(-vb_coff(3));
  dvariable vbdiff=vb_coff(2)-vb_coff(1);
  dvariable vbtmp=1-pow(rho,nage-1);
  int ng=nage;   //NMD_27Sep2018

  for (int j=1;j<=nage;j++)
  {
    dvariable tmp=(1.-pow(rho,j-1))/(1.-pow(rho,nage-1));
    sigma(j)=var_coff(1)*exp(var_coff(2)*(-1+2*tmp));
  }
  int mmin,mmax;
  if (mf_pvm->pvm_switch)
  {
    mmin=min_tag_group;
    mmax=max_tag_group;
  }
  else
  {
    mmin=1;
    mmax=num_tag_releases;
  }
  int isold=1;
  for (int it=mmin;it<=mmax;it++)
  {
    int tr=tag_region(it);
    int lb=1;
    int ub=num_regions;
    int is=1;
    if (pmsd)
    {
      is=pmsd->region_species_pointer(tr);
      lb=pmsd->region_bounds(is,1);
      ub=pmsd->region_bounds(is,2);
    }
    if (is !=isold)
    {
      isold=is;
      vbc=get_vb_coff_region(tr);
      rho=exp(-vbc(3));
      vbdiff=vbc(2)-vbc(1);
      vc=get_var_coff_species(is);
      vbtmp=1-pow(rho,nage-1);
      ng=get_nage_species(is);   //NMD_27Sep2018
      for (int j=1;j<=nage;j++)
      {
        dvariable tmp=(1.-pow(rho,j-1))/(1.-pow(rho,nage-1));
        sigma(j)=vc(1)*exp(vc(2)*(-1+2*tmp));
      }
    }
    for (int ir=lb;ir<=ub;ir++)
    {
      for (int ip=initial_tag_period(it,ir);ip<=num_fish_periods(ir);ip++)
      {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {
          if(tot_tag_catch(it,ir,ip,fi))
          {
            int tyear=year(ir,ip);
            int diff=tyear-tag_year(it);
            const ivector& ocbl=ivector(obstagcatch_by_length(it,ir,ip,fi));
            for (int il=1;il<=tag_nlint;il++)
            {
              if (ocbl(il)>0)
              {
                int num_recaps=ocbl(il);
                MY_DOUBLE_TYPE len=tag_shlen+(il-0.5L)*tag_filen;

//                dvariable age=age_at_length_calc(len,rho,vbdiff,
//                  vbtmp,vbc,nage,parest_flags);
                dvariable age=age_at_length_calc(len,rho,vbdiff,
                  vbtmp,vbc,ng,parest_flags);   //NMD_27Sep2018

                MY_DOUBLE_TYPE cage=value(age);
                if (cage<=1)
                {
                  int ind=1+diff;
                  if (ind>nage) ind=nage;
                  obstagcatch1(it,ir,ip,fi,ind)+=num_recaps;
                }
                else if (cage+diff>=nage_by_tag_release(it))
                {
                  obstagcatch1(it,ir,ip,fi,nage_by_tag_release(it))
                    +=num_recaps;
                }
                else
                {
                  int jj=int(cage)+diff;
                  dvariable sf;
                  sf=daves_kludge1(age);
                  obstagcatch1(it,ir,ip,fi,jj)+=(1.-sf)*num_recaps;
                  obstagcatch1(it,ir,ip,fi,jj+1)+=sf*num_recaps;
                }
              }
            }
          }
          obstagcatch(it,ir,ip,fi)=obstagcatch1(it,ir,ip,fi);
        }
      }
    }
  }
}

void dvar_len_fish_stock_history::observed_tags_by_age_from_length_pooled(void)
{
  obstagcatch1.initialize();
  dvar_vector vbc=vb_coff;
  dvar_vector vc=var_coff;
  dvar_vector sigma(1,nage);
  dvariable rho=exp(-vb_coff(3));
  int ng=nage;    //NMD_27Sep2018
  for (int j=1;j<=nage;j++)
  {
    dvariable tmp=(1.-pow(rho,j-1))/(1.-pow(rho,nage-1));
    sigma(j)=vc(1)*exp(vc(2)*(-1+2*tmp));
  }
  int mmin,mmax;
  if (mf_pvm->pvm_switch)
  {
    mmin=min_tag_group;
    mmax=max_tag_group;
  }
  else
  {
    mmin=1;
    mmax=num_tag_releases;
  }
  int isold=1;
  for (int it=mmin;it<=mmax;it++)
  {
    int tr=tag_region(it);
    int lb=1;
    int ub=num_regions;
    int is=1;
    if (pmsd)
    {
      int is=pmsd->region_species_pointer(tr);
      lb=pmsd->region_bounds(is,1);
      ub=pmsd->region_bounds(is,2);
    }
    if (is !=isold)
    {
      isold=is;
      vbc=get_vb_coff_region(tr);
      rho=exp(-vbc(3));
      //vbdiff=vbc(2)-vbc(1);
      vc=get_var_coff_species(is);
      //vbtmp=1-pow(rho,nage-1);
      ng=get_nage_species(is);   //NMD_27Sep2018
      for (int j=1;j<=nage;j++)
      {
        dvariable tmp=(1.-pow(rho,j-1))/(1.-pow(rho,nage-1));
        sigma(j)=vc(1)*exp(vc(2)*(-1+2*tmp));
      }
    }
    for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
    {
      for (int ip=initial_tag_period(it,ir);ip<=terminal_tag_period(it,ir);ip++)
      {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {
          if(tot_tag_catch(it,ir,ip,fi))
          {
            int tyear=year(ir,ip);
            int diff=tyear-tag_year(it);
            const ivector& ocbl=ivector(obstagcatch_by_length(it,ir,ip,fi));
            for (int il=1;il<=tag_nlint;il++)
            {
              if (ocbl(il)>0)
              {
                int num_recaps=ocbl(il);
                MY_DOUBLE_TYPE len=tag_shlen+(il-0.5L)*tag_filen;
//                dvariable age=age_at_length_calc(len,vbc,nage);
                dvariable age=age_at_length_calc(len,vbc,ng,parest_flags);  //NMD_27Sep2018
                MY_DOUBLE_TYPE cage=value(age);
                if (cage<=1)
                {
                  int ind=1+diff;
                  if (ind>nage) ind=nage;
                  obstagcatch1(it,ir,ip,fi,ind)+=num_recaps;
                }
                else if (cage+diff>=nage)
                {
                  obstagcatch1(it,ir,ip,fi,nage)+=num_recaps;
                }
                else
                {
                  int jj=int(cage)+diff;
                  dvariable sf;
                  sf=daves_kludge(age);
                  obstagcatch1(it,ir,ip,fi,jj)+=(1.-sf)*num_recaps;
                  obstagcatch1(it,ir,ip,fi,jj+1)+=sf*num_recaps;
                }
              }
            }
          }
          obstagcatch(it,ir,ip,fi)=obstagcatch1(it,ir,ip,fi);
        }
      }
    }
  }
  pooledobstagcatch.initialize();

  for (int ir=1;ir<=num_regions;ir++)
  {
    int ipmax=num_fish_periods(ir);  
    for (int ip=minttp(ir)+1;ip<=ipmax;ip++)
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {
        if (pooledtot_tag_catch(ir,ip,fi))
        {
          int tyear=year(ir,ip);
          const ivector& ocbl=ivector(pooledobstagcatch_by_length(ir,ip,fi));
          for (int il=1;il<=tag_nlint;il++)
          {
            if (ocbl(il)>0)
            {
              int num_recaps=ocbl(il);
              pooledobstagcatch(ir,ip,fi,nage)+=num_recaps;
            }
          }
        }
      }
    }
  }
}
