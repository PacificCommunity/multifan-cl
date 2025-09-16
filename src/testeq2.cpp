/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#define HOME_VERSION
#include "all.hpp"
      void check_sanity(imatrix& mpy);
      void check_sanity_trans(imatrix& mpy);
  dvar_matrix mycolumn(const dvar_matrix& _M,int i,int j);

 dvar_matrix matrix_calculations(const dvar3_array  & B)
 {
   int ns=B.indexmax();
   int nr=B(1).indexmax();
   int n=ns*nr;
   dvar_matrix T(1,n,1,n);
   T.initialize();
 
   for (int i=1;i<=n;i++)
   {
     T(i,i)=-1.0;
   }
 
   int off=0;
   int off1=-nr;
   for (int i=1;i<=ns;i++)
   {
     off+=nr;
     off1+=nr;
     if (off+nr>n) off=0;
     dvar_matrix ST=mycolumn(T.sub(off+1,off+nr),off1+1,off1+nr);
     ST.rowshift(1);
     ST.colshift(1);
     ST=B(i);
   }
   return T;
 }    

 void dvar_fish_stock_history::get_equilibrium_cohorts_1
   (ivector& mps,dvar_matrix& initial_R,dvar4_array& surv,imatrix& emi,
   dvar3_array& EN)
 {
   int ns=age_flags(57);
   dvar3_array TRANSPM(1,ns,1,num_regions,1,num_regions);
   int is,im;
 
   dvar_matrix PEN(1,ns,1,num_regions);
   for (int ir=1;ir<=num_regions;ir++)
   {
     PEN.initialize();
     if (!pmsd)
     {
       for (is=1;is<=ns;is++)
       {
         PEN(is,ir)=1;
         //EN(1,is)=initial_R(is);
         dvar_vector tmp(1,num_regions);
         //for (int j=1;j<=nage;j++)
         int j=nage;
         {
           int si=get_season(is,j);
           tmp=PEN(is);
           for (im=1;im<=mps(si);im++)
           {
             tmp = elem_prod(tmp,surv(si,im,j));
             // now do themovement 
             tmp=Dad(emi(si,im),j)*tmp;
           } 
           //if (j>=nage-3)
           //  cout <<"testeq1.cpp " << j << endl;
           PEN(is)=tmp;
         }
         TRANSPM(is,ir)=PEN(is);
       }
     }
     else   // mult species code
     {
       cerr << "Not doen for multi species" << endl;
       ad_exit(1);
      /*
       int nrr=pmsd->num_real_regions;
       int numsp=pmsd->num_species;
   
       dvar3_array TEN(1,nage,1,ns,1,nrr);
       TEN.initialize();
       dvar_vector tmp(1,nrr);
       for (int isp=1;isp<=numsp;isp++)  // loop over species
       {
         int offset=(isp-1)*nrr;
         for (is=1;is<=ns;is++)
         {
           EN(1,is)(1+offset,nrr+offset).shift(1)
             =initial_R(is)(1+offset,nrr+offset).shift(1);
           TEN(1,is)=initial_R(is)(1+offset,nrr+offset).shift(1);
           dvar_vector tmp(1,nrr);
           for (int j=1;j<nage-1;j++)
           {
             int si=get_season(is,j);
             tmp=TEN(j,is);
             for (im=1;im<=mps(si);im++)
             {
               tmp = 
                 elem_prod(tmp,surv(si,im,j).sub(1+offset,nrr+offset).shift(1));
               // now do themovement 
               tmp=
                 Dad(emi(si,im),j).sub(1+offset,nrr+offset).shift(1)*tmp;
             } 
             //if (j>=nage-3)
             //  cout <<"testeq1.cpp " << j << endl;
             TEN(j+1,is)=tmp;
             EN(j+1,is).sub(1+offset,nrr+offset).shift(1)=tmp;
               //cout << tmp << endl;
           }
         }
       }
     */
     }
   }
   dvar3_array PMM(1,ns);
   for (is=1;is<=ns;is++)
   {
     PMM(is)=trans(TRANSPM(is));
   }
 
   dvar_matrix EPLUS=equilibrium_calcs(PMM,EN(nage));
   
   EN(nage)=EPLUS;
 }
 // 
 // 
 dvar_matrix dvar_fish_stock_history::equilibrium_calcs(
   const dvar3_array & B,dvar_matrix & w)
 {
   int ns=age_flags(57);
   
   dvar_matrix T=matrix_calculations(B);
   int nr=num_regions;
 
   dvar_vector ww(1,nr*ns);
   int off=-nr;
   for (int i=1;i<=ns;i++)
   {
     off+=nr;
     ww(off+1,off+nr).shift(1)=w(i);
   }
   dvar_vector solution=-solve(T,ww);
   dvar_matrix EPLUS(1,ns,1,nr);
   off=0;
   for (int i=1;i<=ns;i++)
   {
     EPLUS(i)=solution(off+1,off+nr).shift(1);
     off+=nr;
   }
   return EPLUS;
 }
void dvar_fish_stock_history::init_recruit_equilibrium_code(void)
{
  ivector *pq_flag=0;
  int nya=month_factor; // num_years_for_average
  if (nya==0)
  {
    cerr << "Error month factor must be >0 " << endl;
    ad_exit(1);
  }

  // get the number of movmement periods for each season
  ivector  mps = get_equilibrium_movements_per_season();

  imatrix equilib_move_index = get_initial_equilibrium_movement(mps);

  dvar4_array  surv;
  if (!pmsd)
  {
    surv=get_initial_equilibrium_survival(mps,nya,pq_flag);
  }
  else
  {
    surv=get_initial_equilibrium_survival_ms(mps,nya,pq_flag);
  }

  dvar_matrix  eqrec=get_initial_R();

  dvar3_array EN=get_equilibrium_cohorts(mps,eqrec,equilib_move_index);

  get_equilibrium_cohorts_1(mps,eqrec,surv,equilib_move_index,EN);

  put_P_in_N(EN);

}
 
void dvar_fish_stock_history::put_P_in_N(dvar3_array& EN)
{
  //dvar3_array EN(1,nage,1,ns,1,num_regions);
  int ns=age_flags(57);
  int ir;
  if (af170q0==0)
  {
    for (ir=1;ir<=num_regions;ir++)
    {
      for (int is=1;is<=ns;is++)
      {
        for (int j=1;j<=nage;j++)
        {
          N(ir,is,j)=log(1.e-10+EN(j,is,ir));
        }
      }
    }
  }
  else
  {
    for (ir=1;ir<=num_regions;ir++)
    {
      for (int is=1;is<=ns;is++)
      {
        for (int j=1;j<=nage;j++)
        {
          N_q0(ir,is,j)=log(1.e-10+EN(j,is,ir));
        }
      }
    }
  }
}

