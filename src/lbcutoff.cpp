/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"

MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
dvariable mfexp(dvariable& v1);
dvar_vector mfexp(dvar_vector& );
dvector mfexp(const dvector& );

dvar_vector cutoff_sel(dvar_fish_stock_history& fsh, dvar_vector& _vb_coff,
  dvar_vector& var_coff)
{
  dvar_vector sigma(1,fsh.nage);
  dvar_vector csel(1,fsh.num_fisheries);
  int ip;
  MY_DOUBLE_TYPE month;
  int i;
  for (i=1;i<=fsh.num_fisheries;i++)  // Loop over fisheries
  {

    int ir=fsh.realization_region(i,1);
    int isp=fsh.pmsd->region_species_pointer(ir);
    dvar_vector vb_coff;
    if (isp==1)
    {
      vb_coff=_vb_coff;
    }
    else
    {
      vb_coff=fsh.pmsd->vb_coff(isp);
    }
    dvariable rho=exp(-vb_coff(3));

    ip=fsh.realization_period(i,1);
    month=(fsh.month(1,ip)-1.)/12.;
    for (int j=1;j<=fsh.nage;j++)
    {
      dvariable tmp=(1.-pow(rho,j-1+month))/(1.-pow(rho,fsh.nage-1));
      sigma(j)=var_coff(1)*exp(var_coff(2)*(-1+2*tmp));
    }
  }  
  dvar_matrix esel=exp(fsh.selcoff);
  for (i=1;i<=fsh.num_fisheries;i++)  // Loop over fisheries
  {
    int ir=fsh.realization_region(i,1);
    int isp=fsh.pmsd->region_species_pointer(ir);
    dvar_vector vb_coff;
    if (isp==1)
    {
      vb_coff=_vb_coff;
    }
    else
    {
      vb_coff=fsh.pmsd->vb_coff(isp);
    }
    if (fsh.fish_flags(i,31))
    {
      dvariable alm1;
      dvariable alml;
      int ff3=fsh.fish_flags(i,3);
      MY_DOUBLE_TYPE clength=fsh.fish_flags(i,31)/10.;
      switch (fsh.fish_flags(i,26))
      {
      case 1:
        alm1=age_at_length_calc(clength,vb_coff,fsh.nage);
        alml=age_at_length_calc(vb_coff(2),vb_coff,fsh.nage);
        break;
      case 2:
        alm1=age_at_length_calc(clength,vb_coff,fsh.nage,sigma);
        alml=age_at_length_calc(vb_coff(2),vb_coff,fsh.nage,sigma);
        break;
      default:
        cerr << "Illegal value for fish_flag(" << i << ",26)"
           << endl;
      }
      int im1=static_cast<int>(value(alm1));   
      int iml=static_cast<int>(value(alml));   
      dvariable b=daves_kludge(alm1-im1);
      dvariable bl=daves_kludge(alml-iml);
      csel(i)=(1-b)*esel(i,min(ff3,im1))+b*esel(i,min(ff3,im1+1));
      dvariable tmp=(1-bl)*esel(i,min(ff3,iml))+bl*esel(i,min(ff3,iml+1));
      csel(i)/=tmp;
    }
    else
    {
      csel(i)=0;
    }
  }
  return csel;
}

dvar_matrix cutoff_sel_report(dvar_fish_stock_history& fsh, dvar_vector& _vb_coff,
  dvar_vector& var_coff)
{
  int nint=20;
  dvar_vector sigma(1,fsh.nage);
  dvar_matrix csel(0,fsh.num_fisheries,0,nint);
  int ip;
  MY_DOUBLE_TYPE month;
  int i;
  dvar_vector vb_coff;
  for (i=1;i<=fsh.num_fisheries;i++)  // Loop over fisheries
  {
    int ir=fsh.realization_region(i,1);
    int isp=fsh.pmsd->region_species_pointer(ir);
    if (isp==1)
    {
      vb_coff=_vb_coff;
    }
    else
    {
      vb_coff=fsh.pmsd->vb_coff(isp);
    }
    dvariable rho=exp(-vb_coff(3));
    ip=fsh.realization_period(i,1);
    month=(fsh.month(1,ip)-1.)/12.;
    for (int j=1;j<=fsh.nage;j++)
    {
      dvariable tmp=(1.-pow(rho,j-1+month))/(1.-pow(rho,fsh.nage-1));
      sigma(j)=var_coff(1)*exp(var_coff(2)*(-1+2*tmp));
    }
  }  
  dvar_matrix esel=exp(fsh.selcoff);
  dvariable delta=(vb_coff(2)-vb_coff(1))/nint;
  for (int j=0;j<=nint;j++)
  {
    dvariable clength=vb_coff(1)+j*delta;
    csel(0,j)=clength;
  }
    
  for (i=1;i<=fsh.num_fisheries;i++)  // Loop over fisheries
  {
    int ir=fsh.realization_region(i,1);
    int isp=fsh.pmsd->region_species_pointer(ir);
    dvar_vector vb_coff;
    if (isp==1)
    {
      vb_coff=_vb_coff;
    }
    else
    {
      vb_coff=fsh.pmsd->vb_coff(isp);
    }
    dvariable alm1;
    int ff3=fsh.fish_flags(i,3);
    for (int j=0;j<=nint;j++)
    {
      dvariable clength=vb_coff(1)+j*delta;
      switch (fsh.fish_flags(i,26))
      {
      case 1:
        alm1=age_at_length_calc(clength,vb_coff,fsh.nage);
        break;
      case 2:
        alm1=age_at_length_calc(clength,vb_coff,fsh.nage,sigma);
        break;
      default:
        cerr << "Illegal value for fish_flag(" << i << ",26)"
           << endl;
      }
      int im1=static_cast<int>(value(alm1));   
      dvariable b=daves_kludge(alm1-im1);
      csel(i,j)=(1-b)*esel(i,min(ff3,im1))+b*esel(i,min(ff3,im1+1));
    }
  }
  for (i=1;i<=fsh.num_fisheries;i++)  // Loop over fisheries
  {
    csel(i)/=csel(i,nint);
  }
  return csel;
}

