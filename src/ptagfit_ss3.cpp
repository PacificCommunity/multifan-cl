/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
#include <admodel.h>
#if !defined(linux)
#include <windows.h>
#endif
#define  __declspec(dllexport) 
#ifdef __MSC__
  MY_DOUBLE_TYPE mfexp(const MY_DOUBLE_TYPE& );
  dvariable mfexp(const prevariable& );
  dvar_vector mfexp(const dvar_vector& );
  dvector mfexp(const dvector& );
#else
  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
  dvariable mfexp(BOR_CONST prevariable& );
#endif
dvariable sum(dvar3_array& M);
extern mf_pvm_manager * mf_pvm;
dvariable grouped_dirichlet_multinomial_ll(dmatrix & obs,dvar_matrix & pred,
  dvariable & s);

int check_zero(dvar_vector & v);

dvariable dvar_len_fish_stock_history::fit_pooled_tag_returns_like_ss3(void)
{
  pmsd_error();
  // check if we can exchange loop order
  check_for_ss3_structure();
  dvariable f=0.0;
  dvariable gp_pen=0.0;
  int fi;
  dvar_matrix obsgroupedcatch;
  dvar_matrix groupedcatch;

  ivector group_flags32=column(fish_flags,32);
  int gsum32=sum(group_flags32);
  int gmax32=Max(group_flags32);
  ivector gp_fish32(1,gmax32);
  int rmin=1;
  int rmax=num_regions;
  if (gmax32)
  {
    gp_fish32.initialize();
    for (fi=1;fi<=num_fisheries;fi++)
      gp_fish32(group_flags32(fi))=fi;
    if (allocated(obsgroupedcatch))
      obsgroupedcatch.deallocate();
    obsgroupedcatch.allocate(1,gmax32,1,nage);
    if (allocated(groupedcatch))
      groupedcatch.deallocate();
    groupedcatch.allocate(1,gmax32,1,nage);
  }
  else
  {
    if (allocated(obsgroupedcatch))
      obsgroupedcatch.deallocate();
    obsgroupedcatch.allocate(1,num_fisheries,1,nage);
    if (allocated(groupedcatch))
      groupedcatch.deallocate();
    groupedcatch.allocate(1,num_fisheries,1,nage);
  }
  dvariable totobstags=sum(obstagcatch);
  dvariable totpredtags=0.0;

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
  if (!gsum32)   // no fishery grouping
  {
    imatrix jflag(rmin,rmax,1,gmax32);
    {
      for (int ip=minttp(rmin)+1;ip<=num_real_fish_periods(rmin);ip++)
      {
        jflag=nage;
        ivector nfi=column(num_fish_incidents,ip); 
        dvar3_array pred(rmin,rmax,1,nfi,1,nage);
        dvar3_array obs(rmin,rmax,1,nfi,1,nage);
        dvariable tobs=0.0;
        dvariable tpred=0.0;
        obs.initialize();
        pred.initialize();
        for (int ir=rmin;ir<=rmax;ir++)
        {
          for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
          {
            dvar_vector rtc;
              rtc=rep_rate(ir,ip,fi)*pooled_tagcatch(ir,ip,fi);
            int jmin=rtc.indexmin();
            jflag(ir,fi)=min(jflag(ir,fi),rtc.indexmin());
            pred(ir,fi)(rtc.indexmin(),rtc.indexmax())=rtc;
            obs(ir,fi)(rtc.indexmin(),rtc.indexmax())=pooledobstagcatch(ir,ip,fi); 
            tobs+=sum(pooledobstagcatch(ir,ip,fi)); 
            tpred+=sum(rtc);
          }
        }
        dvariable a=10.0;  // need to get a somewhere
        dvariable multiplier=1.0/tpred;
        for (int ir=rmin;ir<=rmax;ir++) {
          for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) {
            pred(ir,fi)*=multiplier;
          }
        }

        int onetime=0;
        for (int ir=rmin;ir<=rmax;ir++)
        {
          for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
          {
            //if(iflag(fi))
            {
              if (!onetime)
              {
#if !defined(NO_MY_DOUBLE_TYPE)
                f-= gammln(tobs+1.0L);
#else
                f-= gammln(tobs+1.0);
#endif
                f-= gammln(a);
                f+= gammln(tobs+a);
                onetime=1;
              }
              int mmin=jflag(ir,fi);
              int mmax=pred(ir,fi).indexmax();
              f-= sum(gammln(obs(ir,fi)(mmin,mmax)+a*pred(ir,fi)(mmin,mmax)));
              f+= sum(gammln(a*pred(ir,fi))(mmin,mmax));
#if !defined(NO_MY_DOUBLE_TYPE)
              f+= sum(gammln(obs(ir,fi)(mmin,mmax)+1.0L));
#else
              f+= sum(gammln(obs(ir,fi)(mmin,mmax)+1.0));
#endif
            }
          }
        }
      }
    }
  }
  else   // have grouping
  {
    int rmin=1;
    int rmax=num_regions;
    imatrix iflag(rmin,rmax,1,gmax32);
    imatrix jflag(rmin,rmax,1,gmax32);
    for (int ip=minttp(rmin)+1;ip<=num_real_fish_periods(rmin);ip++)
    {
      iflag.initialize();
      jflag=nage;
      dvar3_array pred(rmin,rmax,1,gmax32,1,nage);
      dvar3_array obs(rmin,rmax,1,gmax32,1,nage);
      dvariable tobs=0.0;
      dvariable tpred=0.0;
      obs.initialize();
      pred.initialize();
      for (int ir=rmin;ir<=rmax;ir++)
      {
        ivector& pi=parent(ir,ip);
        for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {
          int pp1=pi(fi);
          int pp=group_flags32(pp1);
          iflag(ir,pp)=1;
          dvar_vector rtc;
          rtc=rep_rate(ir,ip,fi)*pooled_tagcatch(ir,ip,fi);
          if (check_zero(rtc))
          {
            cout << rep_rate(ir,ip,fi) << endl;
          }

          int jmin=rtc.indexmin();
          jflag(ir,pp)=min(jmin,jflag(ir,pp));
          pred(ir,pp)(rtc.indexmin(),rtc.indexmax())+=rtc;
          obs(ir,pp)(rtc.indexmin(),rtc.indexmax())+=pooledobstagcatch(ir,ip,fi); 
          tobs+=sum(pooledobstagcatch(ir,ip,fi)); 
          tpred+=sum(rtc);
        }
      }
      dvariable a=10.0;  // need to get a somewhere
      dvariable multiplier=1.0/tpred;
      for (int ir=rmin;ir<=rmax;ir++) {
        for (int ig=1;ig<=gmax32;ig++) {
          pred(ir,ig)*=multiplier;
        }
      }
      int onetime=0;
      for (int ir=rmin;ir<=rmax;ir++)
      {
        for (int ig=1;ig<=gmax32;ig++)
        {
          if(iflag(ir,ig))
          {
            if (grouped_fishery_projection_flag(ir,ip,ig)==0)
            {
              if (!onetime)
              {
#if !defined(NO_MY_DOUBLE_TYPE)
                f-= gammln(tobs+1.0L);
#else
                f-= gammln(tobs+1.0);
#endif
                f-= gammln(a);
                f+= gammln(tobs+a);
                onetime=1;
              }
              int mmin=jflag(ir,ig);
              int mmax=pred(ir,ig).indexmax();
              f-= sum(gammln(obs(ir,ig)+a*pred(ir,ig)(mmin,mmax)));
              f+= sum(gammln(a*pred(ir,ig)(mmin,mmax)));
              f+= sum(gammln(obs(ir,ig)+1));
            } 
          }
        }
      }
    }
  }
  return f;
}
