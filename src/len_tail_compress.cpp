/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
#include "phi_stuff.h"
#include "tridiagonal_dmatrix.h"
void check_grouping_flags(ivector v);
ivector get_inv_group_ptr(ivector v);

extern int stupid_print_switch;
extern "C" void adfloat_except(int k);


void len_tail_compress(dvar_len_fish_stock_history& fsh)
{
  ivector & sizes =fsh.tcs->sizes;
  ivector & num_samples =fsh.tcs->num_samples;
  ivector& len_rho_group_ptr=fsh.tcs->len_rho_group_ptr;
  ivector& inv_len_rho_group_ptr=fsh.tcs->inv_len_rho_group_ptr;
  dvector& average_sample_size=fsh.tcs->average_sample_size;
  int& tot_nfi=fsh.tcs->tot_nfi;
  int& minsize=fsh.tcs->minsize;
  ivector& reg=fsh.tcs->reg;
  ivector& per=fsh.tcs->per;
  ivector& finc=fsh.tcs->finc;
  ivector& left_bound=fsh.tcs->left_bound;
  ivector& right_bound=fsh.tcs->right_bound;
  imatrix& group_size_flag=fsh.tcs->group_size_flag;
  ivector& is=fsh.tcs->is;
  minsize=1;
  if (fsh.parest_flags(320)>0)
    minsize=fsh.parest_flags(320);
  int nlint=fsh.nlint;
  ivector ff68=column(fsh.fish_flags,68);
  if (sum(ff68)==0)
    ff68.fill_seqadd(1,1);
  len_rho_group_ptr=ff68;
  inv_len_rho_group_ptr=get_inv_group_ptr(len_rho_group_ptr);
  // get number of fish incidents
  tot_nfi=0;
  if(!allocated(group_size_flag))
    group_size_flag.allocate(1,max(len_rho_group_ptr),1,nlint);
  group_size_flag.initialize();
  int maxg=max(len_rho_group_ptr);
  if (!allocated(average_sample_size))
    average_sample_size.allocate(1,maxg);
  if (!allocated(num_samples))
    num_samples.allocate(1,maxg);
   
  MY_DOUBLE_TYPE leftmin=0.0;
  MY_DOUBLE_TYPE rightmin=0.0;
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    int ntimes;
    if (fsh.parest_flags(142)==0)
    {
      ntimes=fsh.num_fish_periods(ir);
    }
    else 
    {
      ntimes=min(fsh.parest_flags(142),fsh.num_fish_periods(ir));
    }

    for (int ip=1;ip<=ntimes;ip++)
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
      {
        int lb=0;
        int rb=0;
        if (fsh.len_sample_size(ir,ip,fi)>0)
        {
          dvector OP=fsh.len_freq(ir,ip,fi);
          MY_DOUBLE_TYPE ssum=0.0;
          for (int k=OP.indexmin();k<=OP.indexmax();k++)
          {
            ssum+=OP(k);
            if (ssum>leftmin) 
            {
              lb=k;
              break;
            }
          }
          ssum=0.0;
          for (int k=OP.indexmax();k>=OP.indexmin();k--)
          {
            ssum+=OP(k);
            if (ssum>rightmin) 
            {
              rb=k;
              break;
            }
          }
          if (rb-lb+1>=minsize)
          tot_nfi++;
          int p=fsh.parent(ir,ip,fi);
          int gp=len_rho_group_ptr(p);
          group_size_flag(gp,rb-lb+1)=1;
        }
      }
    }
  }
  if (!allocated(num_samples))
    num_samples.allocate(1,maxg);
  if (!allocated(reg))
    reg.allocate(1,tot_nfi);
  if (!allocated(per))
    per.allocate(1,tot_nfi);
  if (!allocated(finc))
    finc.allocate(1,tot_nfi);
  if (!allocated(left_bound))
    left_bound.allocate(1,tot_nfi);
  if (!allocated(right_bound))
    right_bound.allocate(1,tot_nfi);

  average_sample_size.initialize();
  num_samples.initialize();
  reg.initialize();
  per.initialize();
  finc.initialize();
  left_bound.initialize();
  right_bound.initialize();
  int ii=0;
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    int ntimes;
    if (fsh.parest_flags(142)==0)
    {
      ntimes=fsh.num_fish_periods(ir);
    }
    else 
    {
      ntimes=min(fsh.parest_flags(142),fsh.num_fish_periods(ir));
    }
    for (int ip=1;ip<=ntimes;ip++)
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
      {
        int lb=0.0;
        int rb=0.0;
        int pp=fsh.parent(ir,ip,fi);
        int gp=len_rho_group_ptr(pp);  //group this fishery is in   
        if (fsh.len_sample_size(ir,ip,fi)>0)
        {
          average_sample_size(gp)+=fsh.len_sample_size(ir,ip,fi);
          num_samples(gp)+=1;
          dvector OP=fsh.len_freq(ir,ip,fi);
          MY_DOUBLE_TYPE ssum=0.0;
          for (int k=OP.indexmin();k<=OP.indexmax();k++)
          {
            ssum+=OP(k);
            if (ssum>leftmin) 
            {
              lb=k;
              break;
            }
          }
          ssum=0.0;
          for (int k=OP.indexmax();k>=OP.indexmin();k--)
          {
            ssum+=OP(k);
            if (ssum>rightmin) 
            {
              rb=k;
              break;
            }
          }
          if (rb-lb+1>=minsize)
          {
            ii++;
            left_bound(ii)=lb;
            right_bound(ii)=rb;
            reg(ii)=ir;
            per(ii)=ip;
            finc(ii)=fi;
          }
        }
      }
    }
  }
  for (int i=1;i<=maxg;i++)
  {
    if (num_samples(i)>0)
      average_sample_size(i)/=num_samples(i);
  }
  sizes=right_bound-left_bound+1;
  int maxsize=max(sizes);

  if (min(sizes)<fsh.m_estimator_coffs.indexmin())
  {
     cerr << "Minsize too small for tail compression" << endl
          << "Need to increase pf320 to at least  " 
          << fsh.m_estimator_coffs.indexmin() << endl;
     ad_exit(1);
  }

  ivector ss=sort(sizes);

  int mmin=sizes.indexmin();
  int mmax=sizes.indexmax();
  ivector tmp(mmin,mmax);
  tmp(mmin)=ss(mmin);
  int ic=mmin;
  for (int i=sizes.indexmin()+1;i<=sizes.indexmax();i++)
  {
    if (ss(i)>ss(i-1))
    {
      tmp(++ic)=ss(i);
    }
  }
  is=tmp(mmin,ic);

}
