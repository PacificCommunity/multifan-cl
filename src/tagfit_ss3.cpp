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

int check_zero(dvar_vector & v)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    if (fabs(value(v(i))<1.e-50)) return 1;
  }
  return 0;
}

dvariable sum(dvar3_array& M)
{
  dvariable tmp=0;
  int mmin=M.indexmin();
  int mmax=M.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    tmp+=sum(M(i));
  }
  return tmp;
}

extern mf_pvm_manager * mf_pvm;

dvariable grouped_dirichlet_multinomial_ll(dmatrix & obs,dvar_matrix & pred,
  dvariable & s)
{
  int cmin=obs.indexmin();
  int cmax=obs.indexmax();
  if (cmin !=pred.indexmin() ||
    cmax !=pred.indexmax())
  {
    cerr << "shape error in grouped_dirichlet_multinomial_ll" << endl;
    ad_exit(1);
  }
  MY_DOUBLE_TYPE tobs=0.0;
  dvariable tpred=0.0;
  dvariable f=0.0;
  for (int i=cmin;i<=cmax;i++)
  {
    tobs+=sum(obs(i));
    tpred+=sum(pred(i));
  }
  dvariable a=s/tpred;
#if !defined(NO_MY_DOUBLE_TYPE)
  f+=gammln(tobs+1.0L);
#else
  f+=gammln(tobs+1.0);
#endif
  f+=gammln(s);
  f-=gammln(tobs+s);
  for (int i=cmin;i<=cmax;i++)
  {
    dvar_vector alpha=a*pred(i);
    f+=sum(gammln(obs(i)+alpha));
#if !defined(NO_MY_DOUBLE_TYPE)
    f-=sum(gammln(obs(i)+1.0L));
#else
    f-=sum(gammln(obs(i)+1.0));
#endif
    f-=sum(gammln(alpha));
  }
  return f;
}  

dvariable dvar_len_fish_stock_history::fit_tag_returns_like_ss3(void)
{
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
  for (int it=mmin;it<=mmax;it++)
  {
    int flatflag=check_flat_tag_time_periods(it);
    if (flatflag)
    {
      cout << "non flatfor tag release " << it << endl;
    }
  }
  //if (!gsum32)   // no fishery grouping
//  if (1)   // no fishery grouping
  if (!gsum32)   // no fishery grouping
  {
    for (int it=mmin;it<=mmax;it++)
    {
      int ng=nage_by_tag_release(it);
      int rmin=tag_region_bounds(1,it);
      int rmax=tag_region_bounds(2,it);
      int ub;
      if (!age_flags(96)) ub=num_fish_periods(rmin);
      else
        ub=terminal_tag_period(it,rmin);
      // DF  july 14 05
      if (ub>num_real_fish_periods(rmin)) {
        ub=num_real_fish_periods(rmin);
      }
      imatrix iflag(rmin,rmax,1,gmax32);
      imatrix jflag(rmin,rmax,1,gmax32);
      for (int ip=initial_tag_period(it,rmin);ip<=ub;ip++)
      {
        iflag.initialize();
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
            if (age_flags(198)) {
              rtc=tag_rep_rate(it,ir,ip,fi)*tagcatch(it,ir,ip,fi);
            }  
            else {
              rtc=rep_rate(ir,ip,fi)*tagcatch(it,ir,ip,fi);
            }

            int jmin=rtc.indexmin();
            pred(ir,fi)(rtc.indexmin(),rtc.indexmax())=rtc;
            obs(ir,fi)(rtc.indexmin(),rtc.indexmax())=obstagcatch(it,ir,ip,fi); 
            tobs+=sum(obstagcatch(it,ir,ip,fi)); 
            tpred+=sum(rtc);
            iflag(ir,fi)=1;   
            jflag(ir,fi)=jmin;   
          }
        }
        dvariable a=10.0;  // need to get a somewhere
        dvariable multiplier=1.0/tpred;
        for (int ir=rmin;ir<=rmax;ir++) {
          for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) {
            pred(ir,fi)*=multiplier;
          }
        }
        //parest_flags(111)=5;

        int onetime=0;
        for (int ir=rmin;ir<=rmax;ir++)
        {
          for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
          {
            if(iflag(ir,fi))
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
              f+= sum(gammln(a*pred(ir,fi)(mmin,mmax)));
#if !defined(NO_MY_DOUBLE_TYPE)
              f+= sum(gammln(obs(ir,fi)+1.0L));
#else
              f+= sum(gammln(obs(ir,fi)+1.0));
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
    for (int it=mmin;it<=mmax;it++)
    {
      if (pmsd)
      {
        int cs=pmsd->tag_species_pointer(it);
        rmin=pmsd->region_bounds(cs,1);
        rmax=pmsd->region_bounds(cs,2);
      }
      int ub;
      if (!age_flags(96))
        ub=num_real_fish_periods(rmin);
      else
        ub=terminal_tag_period(it,rmin);
      // DF  july 14 05
      if (ub>num_real_fish_periods(rmin))
      {
        ub=num_real_fish_periods(rmin);
      }
      imatrix iflag(rmin,rmax,1,gmax32);
      imatrix jflag(rmin,rmax,1,gmax32);
      for (int ip=initial_tag_period(it,rmin);ip<=ub;ip++)
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
            dvar_vector rtc;
            if (age_flags(198)) {
              rtc=tag_rep_rate(it,ir,ip,fi)*tagcatch(it,ir,ip,fi);
            }  
            else {
              rtc=rep_rate(ir,ip,fi)*tagcatch(it,ir,ip,fi);
            }
            int jmin=rtc.indexmin();
            pred(ir,pp)(rtc.indexmin(),rtc.indexmax())+=rtc;
            obs(ir,pp)(rtc.indexmin(),rtc.indexmax())+=obstagcatch(it,ir,ip,fi); 
            tobs+=sum(obstagcatch(it,ir,ip,fi)); 
            tpred+=sum(rtc);
            iflag(ir,pp)=1;   
            jflag(ir,pp)=min(jmin,jflag(ir,pp));
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
                f-= sum(gammln(obs(ir,ig)(mmin,mmax)+a*pred(ir,ig)(mmin,mmax)));
                f+= sum(gammln(a*pred(ir,ig)(mmin,mmax)));
#if !defined(NO_MY_DOUBLE_TYPE)
                f+= sum(gammln(obs(ir,ig)+1.0L));
#else
                f+= sum(gammln(obs(ir,ig)+1.0));
#endif
              } 
            }
          }
        }
      }
    }

    cout << " tagfit before pooling " << f << endl;
  }
  const dvar_vector& a=fish_pars(4)+50.0001;
  //cout << "a = " << a << endl;
  dvar_vector& q=fish_pars(5);
  //cout << "q = " << q << endl;
  if (sum(column(fish_flags,34)))
    gp_pen=grouped_tag_reporting_rate_penalty();
  if (age_flags(105)>0)
  {
    dvariable tot_pen=age_flags(105)/100.*square(totobstags-totpredtags);
    cout << "Total observed tags = " << totobstags << endl;
    cout << "Total predicted tags = " << totpredtags << endl;
    cout << "Total tags penalty =   " << tot_pen << endl;
    f+=tot_pen;
  }

  if (age_flags(96))
  {
    if (mf_pvm->pvm_switch == 0 || mf_pvm->pvm_switch == 1)
    {
      dvariable tmp=fit_pooled_tag_returns_like_ss3();
      cout << " tagfit pooling " << tmp << endl;
      f+=tmp;
    }
  }
  f+=gp_pen;
  return f;
}
void check_equal_bounds(ivector bds,const char s[])
{
  int mmin=bds.indexmin();
  int mmax=bds.indexmax();
  int bs=bds(mmin);
  for (int i=mmin+1;i<=mmax;i++)
  {
    if (bs != bds(i))
    {
      cerr << s << endl;
      ad_exit(1);
    }
  }
}

void dvar_len_fish_stock_history::check_for_ss3_structure(void)
{
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
  for (int it=mmin;it<=mmax;it++)
  {
    int ng=nage_by_tag_release(it);
    int rmin=tag_region_bounds(1,it);
    int rmax=tag_region_bounds(2,it);
    ivector region_lower_bound(rmin,rmax);
    ivector region_upper_bound(rmin,rmax);
    for (int ir=rmin;ir<=rmax;ir++)
    {
      int ub;
      if (!age_flags(96))
        ub=num_fish_periods(ir);
      else
        ub=terminal_tag_period(it,ir);
      // DF  july 14 05
      if (ub>num_real_fish_periods(ir))
      {
        ub=num_real_fish_periods(ir);
      }
      region_lower_bound=initial_tag_period(it,ir);
      region_upper_bound=ub;

      check_equal_bounds(region_lower_bound,"region_lower_bound");
      check_equal_bounds(region_upper_bound,"region_upper_bound");
    }
  }
}

