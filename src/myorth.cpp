/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#define USE_DD_NOT
#include "all.hpp"

#if defined(USE_DD)
void ddvector::fill_seqadd(MY_DOUBLE_TYPE x,MY_DOUBLE_TYPE y)
{
  int mmin=indexmin();
  int mmax=indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    (*this)(i)=x+(i-mmin)*y;
  }
}
#endif

dmatrix dvar_fish_stock_history::test_orthogonal(void)
{
  int ny;
  int nr;
  int yd;
  int rd;
  int sd;
  int yd1;
  int sd1;
  int rd1;
  int ns;

  if (pmsd==0  || pmsd->current_species==1 )
  {
    ns=age_flags(57);
    int rem=last_real_year%ns;
    int na=last_real_year/ns;
    if (rem) na++;
    ny=na;
    nr=get_region_bounds()(2)-get_region_bounds()(1)+1;

    yd=recr_degree_yr;
    rd=recr_degree_reg;
    sd=recr_degree_ses;
    yd1=yd+1;
    sd1=sd+1;
    rd1=rd+1;
    num_new_weights=yd1*sd1*rd1;
    if (!allocated(new_orth_recr))
    {
      new_orth_recr.allocate(1,num_new_weights);
    }
  }
  else
  {
    int cs=pmsd->current_species;
    ns=age_flags(57);
    int rem=last_real_year%ns;
    int na=last_real_year/ns;
    ny=na;
    nr=get_region_bounds()(2)-get_region_bounds()(1)+1;
    yd=pmsd->recr_degree_yr(cs);
    rd=pmsd->recr_degree_reg(cs);
    sd=pmsd->recr_degree_ses(cs);
    yd1=yd+1;
    sd1=sd+1;
    rd1=rd+1;
    pmsd->num_new_weights(cs)=yd1*sd1*rd1;
    if (!allocated(pmsd->new_orth_recr(cs)))
    {
      pmsd->new_orth_recr(cs).allocate(1,num_new_weights);
    }
  }

  dmatrix Y(0,yd,1,ny);
  dmatrix S(0,sd,1,ns);
  dmatrix R(0,rd,1,nr);

  Y(0)=1.0;
  S(0)=1.0;
  R(0)=1.0;
 
#if !defined(NO_MY_DOUBLE_TYPE)
  if (yd>0) Y(1).fill_seqadd(-1.0L,2.0/(ny-1.0L));
#else
  if (yd>0) Y(1).fill_seqadd(-1.0,2.0/(ny-1.0));
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
  if (sd>0) S(1).fill_seqadd(-1.0L,2.0/(ns-1.0L));
#else
  if (sd>0) S(1).fill_seqadd(-1.0,2.0/(ns-1.0));
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
  if (rd>0) R(1).fill_seqadd(-1.0L,2.0/(nr-1.0L));
#else
  if (rd>0) R(1).fill_seqadd(-1.0,2.0/(nr-1.0));
#endif

  for (int i=2;i<=yd;i++)
  {
    Y(i)=elem_prod(Y(i-1),Y(1));
  }
  for (int i=2;i<=sd;i++)
  {
    S(i)=elem_prod(S(i-1),S(1));
  }
  for (int i=2;i<=rd;i++)
  {
    R(i)=elem_prod(R(i-1),R(1));
  }

  int ll=0;
  dmatrix O(1,yd1*sd1*rd1,1,ny*nr*ns);
  int N=yd1*sd1*rd1;
  //ddmatrix O(1,1000,1,1000);
  O.initialize();
  cout << Y(0) << endl;
  cout << S(0) << endl;
  cout << R(0) << endl;
  for (int i=0;i<=yd;i++)
  {
    for (int j=0;j<=sd;j++)
    {
      for (int k=0;k<=rd;k++)
      {
        ll++;
        int mm=0;
        for (int ii=1;ii<=ny;ii++)
        {
          for (int jj=1;jj<=ns;jj++)
          {
            for (int kk=1;kk<=nr;kk++)
            {
              O(ll,++mm)=Y(i,ii)*S(j,jj)*R(k,kk);
            }
          }
        }
      }
    }
  }
  for (int i=1;i<=N;i++)
  {
    O(i)/=norm(O(i));
    for (int j=i+1;j<=N;j++)
    {
      O(j)-=(O(j)*O(i))*O(i);
    }
  }
#if defined(USE_DD)
  return make_dmatrix(O);
#else
  return O;
#endif
}
   

void dvar_fish_stock_history::expand_orthogonal_coefficients(void)
{
  if (historical_age_flags(137)> age_flags(137)
    || historical_age_flags(138)> age_flags(138)
    || historical_age_flags(139)> age_flags(139))
  {
     cerr << "Can not reduce degree of orthogonal recruitment" << endl;
     ad_exit(1);
  }
  if (historical_age_flags(137)< age_flags(137)
    || historical_age_flags(138)< age_flags(138)
    || historical_age_flags(139)< age_flags(139))
  {
    int ydold=historical_age_flags(137)-1;
    int rdold=historical_age_flags(138)-1;
    int sdold=historical_age_flags(139)-1;
    int yd=age_flags(137)-1;
    int rd=age_flags(138)-1;
    int sd=age_flags(139)-1;
    int yd1=yd+1;
    int rd1=rd+1;
    int sd1=sd+1;
    int nl=0;
    int ol=0;
    dvar_vector tmp(1,yd1*sd1*rd1);
    tmp.initialize();
    for (int i=0;i<=yd;i++)
    {
      for (int j=0;j<=sd;j++)
      {
        for (int k=0;k<=rd;k++)
        {
          nl++;
          if (i<=ydold && j<= sdold && k<=rdold)
          {
            tmp(nl)=new_orth_recr(++ol);
          }
        }
      }
    }
    cout << (new_orth_recr.shape->get_ncopies()) << endl;
    if ( (new_orth_recr.shape->get_ncopies()) >1)
    {
      cerr << "Expanding orthogonal-polynomial recruitment parameters" << endl;
      cerr << "exceeds bounds - exiting" << endl;
      ad_exit(1);
    }
    new_orth_recr.deallocate();
    new_orth_recr=tmp;
  }
  if (pmsd && pmsd->num_species>1)
  {
    for (int is=2;is<=pmsd->num_species;is++)
    {
      if (pmsd->historical_age_flags(is,137)> pmsd->age_flags(is,137)
        || pmsd->historical_age_flags(is,138)> pmsd->age_flags(is,138)
        || pmsd->historical_age_flags(is,139)> pmsd->age_flags(is,139))
      {
         cerr << "Can not reduce degree of orthogonal recruitment" << endl;
         ad_exit(1);
      }
      if (pmsd->historical_age_flags(is,137)< pmsd->age_flags(is,137)
        || pmsd->historical_age_flags(is,138)< pmsd->age_flags(is,138)
        || pmsd->historical_age_flags(is,139)< pmsd->age_flags(is,139))
      {
        int ydold=pmsd->historical_age_flags(is,137)-1;
        int rdold=pmsd->historical_age_flags(is,138)-1;
        int sdold=pmsd->historical_age_flags(is,139)-1;
        int yd=pmsd->age_flags(is,137)-1;
        int rd=pmsd->age_flags(is,138)-1;
        int sd=pmsd->age_flags(is,139)-1;
        int yd1=yd+1;
        int rd1=rd+1;
        int sd1=sd+1;
        int nl=0;
        int ol=0;
        dvar_vector tmp(1,yd1*sd1*rd1);
        tmp.initialize();
        for (int i=0;i<=yd;i++)
        {
          for (int j=0;j<=sd;j++)
          {
            for (int k=0;k<=rd;k++)
            {
              nl++;
              if (i<=ydold && j<= sdold && k<=rdold)
              {
                tmp(nl)=pmsd->new_orth_recr(is,++ol);
              }
            }
          }
        }
        cout << pmsd->new_orth_recr(is).shape->get_ncopies() << endl;
        if ( (pmsd->new_orth_recr(is).shape->get_ncopies()) >1)
        {
          cerr << "Error: orthogonal-polynomial recruitment parameters out"
	       << "of bounds" << endl;
          ad_exit(1);
        }
        pmsd->new_orth_recr(is).deallocate();
        pmsd->new_orth_recr(is)=tmp;
      }
    }
  }
}
