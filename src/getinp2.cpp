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

void dvar_len_fish_stock_history::get_population_multipliers(dvar_vector& sv,
  void * pq_flag)
{
  totalcatch_by_numbers.initialize();
  obstotalcatch_by_numbers.initialize();
  numtotalcatch_by_numbers.initialize();
  totalcatch_by_weight.initialize();
  obstotalcatch_by_weight.initialize();
  numtotalcatch_by_weight.initialize();
  dvariable A;
  dvariable B;
  if (!pq_flag)
  {
    totpop=0.0;
  }
  pmsd_error(); //NMD 10Apr2012
  get_initial_population(sv,1,pq_flag);
 
  get_numbers_at_age(sv,pq_flag);
  proportion_at_age_calc();
//  fit_totals(0,1);
  fit_totals(0,1);  //NMD_11dec2023
  A=sum(totalcatch_by_numbers)+1.e-3;
  B=sum(totalcatch_by_weight)+1.e-3;

  if (!pq_flag)
  {
    dvariable tmp1=sum(obstotalcatch_by_numbers)/A;
    dvariable tmp2=sum(obstotalcatch_by_weight)/B;
    MY_DOUBLE_TYPE s1=sum(numtotalcatch_by_numbers);
    MY_DOUBLE_TYPE s2=sum(numtotalcatch_by_weight);
    dvariable mult=(s1*tmp1+s2*tmp2)/(s1+s2);
    if (mult<=0.0)
    {
      cerr << "Mult<=0.0 " << mult << endl;
      ad_exit(1);
    }
    totpop=log(mult)+totpop_coff;  //+totpop_coff;
  }
}

dvariable dvar_fish_stock_history::get_initial_population(dvar_vector& sv,
  int avg_calc_flag,void * _pq_flag)
{  
  ivector * pq_flag = (ivector *) (_pq_flag);
  dvariable fpen=0.0;
  dvar3_array * pN=0;
  if (af170q0==0)
  {
    pN=&N;
  }
  else
  {
    pN=&N_q0;
  }
  if (parest_flags(155)==0)
  {
    recinpop_standard(pq_flag,pN,4);
  }
  else if (parest_flags(155)<0)
  {
    orth_pf155_less_than_zero_why(pN);
  }
  else if (initial_orthp_estimate_flag==0)
  {
    fpen+=recinpop_orth(pN,4);
  }
    get_the_initial_age_structure_options(pq_flag);
  return fpen;
}

void dvar_fish_stock_history::get_the_initial_age_structure_options
  (ivector * pq_flag)
{
  //dvar_vector& sv, int avg_calc_flag,void * _pq_flag)//
  // get the initial age structure
  if (!age_flags(94))
  {
    get_initial_age_structure(totpop,sv);
  }
  else
  {  
    if (num_regions==1)
    {
      /*
      if (!age_flags(92))   //NMD_2jul2021
      {
        xget_initial_age_structure_equilibrium();
      }
      else
      {
      */
        // **************************************
        // **************************************
         if (age_flags(94)==1 || age_flags(94)==2)   //NMD_2jul2021
           test_initial_equilibrium_code(pq_flag);
         else if (age_flags(94)==3)
           init_recruit_equilibrium_code();
         else
         {
           cerr << "Invalid option for age_flags(94)" << endl;
           ad_exit(1);
         }
        // **************************************
        // **************************************
       //}
    }
    else
    {
      if (age_flags(94)==1 || age_flags(94)==2)   //NMD_2jul2021
        test_initial_equilibrium_code(pq_flag);
      else if (age_flags(94)==3)
        init_recruit_equilibrium_code();
      else
      {
        cerr << "Invalid option for age_flags(94)" << endl;
        ad_exit(1);
      }
    }
  }
}

static MY_DOUBLE_TYPE dot(const dmatrix&A,const dmatrix&B)
{
  int imin=A.indexmin();
  int imax=A.indexmax();
  
  int jmin=A(imin).indexmin();
  int jmax=A(imin).indexmax();

  MY_DOUBLE_TYPE sum=0.0;
  for (int i=imin;i<=imax;i++)
  { 
    for (int j=jmin;j<=jmax;j++)
    { 
      sum+=A(i,j)*B(i,j);
    }
  }
  return sum;
}

void dvar_fish_stock_history::get_pop_delta(void)
{
  dvar_vector tmp(1,num_regions);
  tmp=region_pars(1)+1.e-12;
  if (!pmsd || pmsd->num_species==1)
  {
    tmp/=sum(tmp);
  }
  else
  {
    int ns=pmsd->num_species;
    for (int is=1;is<=ns;is++)
    {
      int lb=pmsd->region_bounds(is,1);
      int ub=pmsd->region_bounds(is,2);
      tmp(lb,ub)/=sum(tmp(lb,ub));
    }
  }
  epop_delta=tmp;
  pop_delta=log(tmp);
}
  
void dvar_fish_stock_history::orth_pf155_less_than_zero_why(dvar3_array * pN)
{
  int ns=age_flags(57);
  int rem=last_real_year%ns;
  int na=last_real_year/ns;
  if (rem) na++;
  int ny=na;
  int nr=get_region_bounds(1)(2)-get_region_bounds(1)(1)+1;
  {
    int mmin=OR.rowmin();
    int mmax=OR.rowmax();
    int cmin=OR(mmin).indexmin();
    int cmax=OR(mmin).indexmax();
    dvar_vector RR(cmin,cmax);
    RR.initialize();
    //ofstream ofs("orthstuff");
    for (int i=mmin;i<=mmax;i++)
    {
      RR+=new_orth_recr(i)*OR(i);
      //ofs << i << " " << new_orth_recr(i) << " " << OR(i) << endl;
    }
    //ofs.close();
    int mm=0;
    dvar3_array Rec(1,ny,1,ns,1,nr);
    for (int ii=1;ii<=ny;ii++)
    {
      for (int jj=1;jj<=ns;jj++)
      {
        for (int kk=1;kk<=nr;kk++)
        {
          Rec(ii,jj,kk)=RR(++mm);
        }
      }
    }
    for (int ir=1;ir<=nr;ir++)
    {
      for (int iy=1;iy<=nyears;iy++)
      {
        int ty=(iy-1)/num_seasons+1;
        int ts=(iy-1)%num_seasons+1;
        {
          (*pN)(ir,iy,1)=15.0+Rec(ty,ts,ir);
        }
      }
    }
  }        
  if (pmsd)
  {
    for (int is=2;is<=pmsd->num_species;is++)
    {
      int rlow=get_region_bounds(is)(1);
      int nr=get_region_bounds(is)(2)-get_region_bounds(is)(1)+1;
      dmatrix& POR=pmsd->OR(is);
      int mmin=POR.rowmin();
      int mmax=POR.rowmax();
      int cmin=POR(mmin).indexmin();
      int cmax=POR(mmin).indexmax();
      dvar_vector RR(cmin,cmax);
      RR.initialize();
      for (int i=mmin;i<=mmax;i++)
      {
        RR+=pmsd->new_orth_recr(is,i)*POR(i);
      }
      int mm=0;
      dvar3_array Rec(1,ny,1,ns,1,nr);
      for (int ii=1;ii<=ny;ii++)
      {
        for (int jj=1;jj<=ns;jj++)
        {
          for (int kk=1;kk<=nr;kk++)
          {
            Rec(ii,jj,kk)=RR(++mm);
          }
        }
      }
      for (int ir=1;ir<=nr;ir++)
      {
        for (int iy=1;iy<=nyears;iy++)
        {
          int ty=(iy-1)/num_seasons+1;
          int ts=(iy-1)%num_seasons+1;
          {  
            (*pN)(ir+rlow-1,iy,1)=15.0+Rec(ty,ts,ir);
          }
        }
      }
    }
  }        
}
#undef HOME_VERSION


