/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


#define HOME_VERSION
#include "all.hpp"

#include "variable.hpp"
#include "newmprot.hpp"
int dvar_len_fish_stock_history::check_freq_pooling_flags(int data_type)
{
  int no_pool_flag=0;
  if (data_type==1)
  {
    if (parest_flags(343))
      no_pool_flag=1;
  }
  else if (data_type==2)
  {
    if (parest_flags(344))
      no_pool_flag=1;
  }
  else 
  {
    cerr << "Illegal value for data_type" << endl;
    ad_exit(1);
  }
  return no_pool_flag;
}

int dvar_fish_stock_history::check_total_catch_pooling_flags(void)
{
  // need to figure out options for this function
  int no_pool_flag=0;
  if (parest_flags(345)) no_pool_flag=1;
  return no_pool_flag;
}

dvariable dvar_len_fish_stock_history::square_fita(d3_array& tf,
  d3_array & sample_size,d4_array & of,dvar4_array & pf,
  dvector& contrib, d3_array& contrib_by_realization,int data_type,
  int print_switch)
{
  ofstream * ppofs = 0;
  if (print_switch)
  {
    if (data_type==1)
    {
      ppofs = new ofstream("lsizemult");
    }
    else if (data_type==2)
    {
      ppofs = new ofstream("wsizemult");
    }
  }

  i3_array fisn_lw;
  i4_array reg_in_lw;
  if (pmsd && data_type==1)
  {
    fisn_lw=pmsd->fisn_lf;
    reg_in_lw=pmsd->reg_in_lf;
  }
  else if (pmsd && data_type==2)
  {
    fisn_lw=pmsd->fisn_wf;
    reg_in_lw=pmsd->reg_in_wf;
  }
  int no_pool_flag=check_freq_pooling_flags(data_type);
  dvariable tot_tmp=0.;
  dvariable tmp_sigma=0.;
  contrib=0.;
  for (int ir=1;ir<=num_regions;ir++)
  {
    int ntimes;
    if (parest_flags(142)==0)
    {
      ntimes=num_fish_periods(ir);
    }
    else 
    {
      ntimes=min(parest_flags(142),num_fish_periods(ir));
    }

    for (int ip=1;ip<=ntimes;ip++)
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)
      {
        if (sample_size(ir,ip,fi)>0)
        {
          dvector varQ;
          int mmax=(of(ir,ip,fi).indexmax());
          int mmin=(of(ir,ip,fi).indexmin());
          int nint=mmax-mmin+1;
          MY_DOUBLE_TYPE eps=1.0/nint;
          if (parest_flags(193))
            eps=parest_flags(193)/(100.*nint);
          dvariable len_tmp=0.;

          if (ppofs)
          {
            // effective sample size calcs
            // assumes that of and pf are vectors of proportions, i.e.
            // that they both sum to 1.

            dvector cpf=value(pf(ir,ip,fi));
            dvector denom=elem_prod(cpf,1.0-cpf);
            dvector numer=square(of(ir,ip,fi)-cpf);
            int pp=parent(ir,ip,fi);
            dvector tmp(1,nint);
            tmp.initialize();
            int ioff=0;
//            for (int i=1;i<=nint;i++)
            for (int i=mmin;i<=mmax;i++)     //NMD_7jun2019
            {
              if (denom(i) !=0)
              {
                ioff++;
                tmp(ioff)=numer(i)/denom(i);
              }
            }
            MY_DOUBLE_TYPE ess;
	    if (sum(tmp) < 1e-10)   //NMD_7jun2019
            {
              ess=0.0;
            }
	    else
            {
              ess=ioff/sum(tmp);
            }   //NMD_7jun2019
            (*ppofs) << ir << " " << ip << " " << fi << " "
                  << pp << " "
//                  << ioff/sum(tmp) << " " << nint << " "   //NMD_7jun2019
                  << ess << " " << nint << " " 
                  <<  sample_size(ir,ip,fi) << " "
                  << sum(tmp) << endl;
          }


          if (no_pool_flag ||!pmsd || fisn_lw(ir,ip,fi)==1)
          {
            MY_DOUBLE_TYPE tau2=1.0/tf(ir,ip,fi);
            dvector esquigle=elem_prod(1.0-of(ir,ip,fi),
                           of(ir,ip,fi));
            varQ=(esquigle+eps)*tau2;;
            dvar_vector r2=elem_div(
              square(of(ir,ip,fi)-pf(ir,ip,fi)),
              2*varQ);
            len_tmp -= sum(log(exp(-r2) + 0.001));
          }
          else
          {
            if (ir==reg_in_lw(ir,ip,fi,1))
            {
              int mmin=pf(ir,ip,fi).indexmin();
              int mmax=pf(ir,ip,fi).indexmax();
              dvar_vector tmp(mmin,mmax);
              MY_DOUBLE_TYPE tot=tf(ir,ip,fi);

// NMD_23jun-17
              if (!parest_flags(350))  // original case
              {
                tmp=pf(ir,ip,fi);
                int n=fisn_lw(ir,ip,fi);
                for (int i=2;i<=n;i++)
                {
                  int rr=reg_in_lw(ir,ip,fi,i);
                  tmp+=pf(rr,ip,fi);
                  // need to deal with this
                  //tot+=tf(rr,ip,fi);
                }
                tmp/=double(n);
              }
              else  // new case applying sex ratio of predicted catch
              {
                dvar_vector sp_ratio(1,fisn_lw(ir,ip,fi)); // YT 2017-06-24 
                dvariable tmpc=tot_catch(ir,ip,fi);
                int n=fisn_lw(ir,ip,fi);
                for (int i=2;i<=n;i++)
                {
                  int rr=reg_in_lw(ir,ip,fi,i);
                  tmpc+=tot_catch(rr,ip,fi);
                }
                for (int i=1;i<=n;i++)
                {
                  int rr=reg_in_lw(ir,ip,fi,i);
                  sp_ratio(i)=tot_catch(rr,ip,fi)/tmpc;
                }
// apply sex ratio to predicted proportions
                tmp=pf(ir,ip,fi)*sp_ratio(1);
                for (int i=2;i<=n;i++)
                {
                  int rr=reg_in_lw(ir,ip,fi,i);
                  tmp+=pf(rr,ip,fi)*sp_ratio(i);
                }
              }
// NMD_23jun-17
              MY_DOUBLE_TYPE tau2=1.0/tot;
              dvector esquigle=elem_prod(1.0-of(ir,ip,fi),
                           of(ir,ip,fi));
              varQ=(esquigle+eps)*tau2;;
              dvar_vector r2=elem_div(
                square(of(ir,ip,fi)-tmp),
                2*varQ);
              len_tmp -= sum(log(exp(-r2) + 0.001));
            }
          }
          if (parest_flags(161)==0)
          {
            if (allocated(varQ))
              len_tmp+=0.5*sum(log(6.283186*varQ));
          }

          //contrib by fishery ....PK jun26-07
          contrib(parent(ir,ip,fi)) += value(len_tmp);
          contrib_by_realization(ir,ip,fi) = value(len_tmp);
          tot_tmp += len_tmp;
        }
      }
    }
  }
  return tot_tmp;
}


