/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"

void dvar_len_fish_stock_history::set_global_variance(void)
{
  int j;
  //if (!set_global_vars_flag)
  dvariable rho=exp(-vb_coff(3));
  dvariable temp2=1.-pow(rho,nage-1);
  dvariable xn1inv=1.0/temp2;
  if (parest_flags(226)==0)
  {
    dvar_vector u(1,nage);
    dvar_vector local_sdevs(1,nage);
    //dvariable sd=var_coff(1);
    MY_DOUBLE_TYPE v;
    for (j=1;j<=nage;j++)         // Loop over age classes
    {
      if (parest_flags(34) == 0)
      {
        dvariable tmp=(1.-pow(rho,j-1))*xn1inv;
        u(j)=(-1.0+2.0*tmp);
      }
      else if (parest_flags(34) == 1)
      {
        u(j)= ( 2.*(j-1.e0)-(nage-1.0))/(nage-1.e0);
      }
      else if (parest_flags(34) == 2)
      {
        cout << "WARNING: incorrect scaled mean-length calculated!!!" << endl;
        v= j-1.e0;
        u(j)= -1.e0+2.e0*(1.0-pow(vb_coff(3),v))/
           (1.0-pow(vb_coff(3),nage-1.0));
      }
      else
      {
        cerr << "ERROR: incorrect option for parest_flags(34)" << endl;
        ad_exit(1);
      }
    }
    dvariable tt;
    dvariable v1=var_coff(1);
    dvariable v2=var_coff(2);
    for (j=1;j<=nage;j++)         // Loop over age classes
    {
      if (parest_flags(35)==0)  //log-linear relationship
      {
        local_sdevs(j)=v1*exp(v2*u(j));
        global_vars(j)=square(local_sdevs(j));
      }
      else if (parest_flags(35)==1)   //normal linear relationship
      {
        local_sdevs(j)=v1+(v2*u(j));
        global_vars(j)=square(local_sdevs(j));
      }
      else
      {
        cerr << "ERROR: incorrect option for parest_flags(35)" << endl;
        ad_exit(1);
      }
    }
  }
  else     //NMD_5Oct2018
  {
    dvariable T;
    if (parest_flags(226)==1)
      T=exp(vb_coff(4));   // change sign to make simpler
    else
      T=-exp(vb_coff(4));
    dvar_vector u(1,nage);
    dvar_vector local_sdevs(1,nage);
    dvariable c1=pow(vb_coff(1),1.0/T);
    dvariable cN=pow(vb_coff(2),1.0/T);
    dvariable diff=cN-c1;
    dvariable rho=exp(-vb_coff(3));
    dvariable temp2=1.-pow(rho,nage-1);
    dvariable xn1inv=1.0/temp2;
    for (int j=1;j<=nage;j++)
    {
      dvariable tmp=(1.-pow(rho,j-1))*xn1inv;
      dvariable tt=c1+diff*tmp;
      u(j)=pow(tt,T);
    }
    dvariable v1=var_coff(1);
    dvariable v2=var_coff(2);
    dvar_vector scaled_ml=setm11(u);   //scale the ml between -1 and 1
    local_sdevs=v1*exp(v2*scaled_ml);
    global_vars=square(local_sdevs);
  }     //NMD_5Oct2018


  if (pmsd)
  {

    for (int is=2;is<=pmsd->num_species;is++)
    {
      if (parest_flags(226)==0)
      {
        dvar_vector vbc=pmsd->vb_coff(is);
        dvar_vector vc=pmsd->var_coff(is);
        int ng=pmsd->nage(is);
        if(!allocated(pmsd->global_vars(is)))
          pmsd->global_vars(is).allocate(1,ng);
        dvar_vector u(1,ng);
        dvar_vector local_sdevs(1,ng);
        //dvariable sd=var_coff(1);
        MY_DOUBLE_TYPE v;
        for (j=1;j<=ng;j++)         // Loop over age classes
        {
//NMD_26oct2023
          if (parest_flags(34) == 0)
          {
            dvariable tmp=(1.-pow(rho,j-1))*xn1inv;
            u(j)=(-1.0+2.0*tmp);
          }
          else if (parest_flags(34) == 1)
          {
            u(j)= ( 2.*(j-1.e0)-(ng-1.0))/(ng-1.e0);
          }
          else if (parest_flags(34) == 2)
          {
            cout << "WARNING: incorrect scaled mean-length calculated!!!" << endl;
            v= j-1.e0;
            u(j)= -1.e0+2.e0*(1.0-pow(vb_coff(3),v))/
               (1.0-pow(vb_coff(3),ng-1.0));
          }
          else
          {
            cerr << "ERROR: incorrect option for parest_flags(34)" << endl;
            ad_exit(1);
          }
//NMD_26oct2023
        }
        dvariable tt;
        dvariable v1=vc(1);
        dvariable v2=vc(2);
        for (j=1;j<=ng;j++)         // Loop over age classes
        {
          local_sdevs(j)=v1*exp(v2*u(j));
          pmsd->global_vars(is,j)=square(local_sdevs(j));
        }
//        cout << "pmsd->global_vars(2,1) = " << pmsd->global_vars(2,1) << endl;
      }
      else
      {
        dvar_vector vbc=pmsd->vb_coff(is);
        dvar_vector vc=pmsd->var_coff(is);
        int ng=pmsd->nage(is);
        if(!allocated(pmsd->global_vars(is)))
          pmsd->global_vars(is).allocate(1,ng);
        dvariable T;
        if (parest_flags(226)==1)
          T=exp(vbc(4));   // change sign to make simpler
        else
          T=-exp(vbc(4));
        dvar_vector u(1,ng);
        dvar_vector local_sdevs(1,ng);
        dvariable c1=pow(vbc(1),1.0/T);
        dvariable cN=pow(vbc(2),1.0/T);
        dvariable diff=cN-c1;
        dvariable rho=exp(-vbc(3));
        dvariable temp2=1.-pow(rho,ng-1);
        dvariable xn1inv=1.0/temp2;
        for (int j=1;j<=ng;j++)
        {
          dvariable tmp=(1.-pow(rho,j-1))*xn1inv;
          dvariable tt=c1+diff*tmp;
          u(j)=pow(tt,T);
        }
        dvariable v1=vc(1);
        dvariable v2=vc(2);
        dvar_vector scaled_ml=setm11(u);   //scale the ml between -1 and 1
        local_sdevs=v1*exp(v2*scaled_ml);
        pmsd->global_vars(is)=square(local_sdevs);
      }
    }
  }
  set_global_vars_flag=1;
}
