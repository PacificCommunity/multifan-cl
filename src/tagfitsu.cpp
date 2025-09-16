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
extern mf_pvm_manager * mf_pvm;

dvariable dvar_len_fish_stock_history::fit_tag_returns_survival_analysis(void)
{
  dvariable f=0.0;
  dvariable gp_pen=0.0;
  int fi;
  dvar_vector& rep_rate1=fish_pars(3);

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
  // here terminal tag period should be when the tag is caught
  // and reported
  //if (!gsum32)   // no fishery grouping
  //ofstream ofs("tagpred");
  {
    for (int it=mmin;it<=mmax;it++)
    {
      if (itr(it)>0)
      {
        //ofs << "tag group " << it << endl;
        dvariable totpredtags=0.0;
        for (int ir=1;ir<=num_regions;ir++)
        {
          int ub=terminal_tag_period(it,ir);
          for (int ip=initial_tag_period(it,ir);ip<=ub;ip++)
          {
            for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
            {
              int mmin=tagcatch(it,ir,ip,fi).indexmin();
              for (int j=mmin;j<=nage;j++)
              {
                if (tagcatch(it,ir,ip,fi,j)<0)
                {
                  cerr << tagcatch(it,ir,ip,fi,j) << endl;
                }
              }
              const dvar_vector& rtc=rep_rate(ir,ip,fi)*tagcatch(it,ir,ip,fi);
              //ofs << sum(rtc) << endl;
              // totpredtags should be the prob that the tag is caught
              // and reported between the time when it was released and
              // when it was caught and reported
              // this sum is summing over age classes
              totpredtags+=sum(rtc);
            }
          }
        }
        //ofs << " tot_pred_tags = " << totpredtags << endl;
  
        // at this point we need ir,ip,ir for when the tag
        // was caught and we add the likelihood component
        // this is the region for this fishery
        // 
        int fishery=tag_recaptures_by_length(it,1,2);
        int irr=fishery_regions(fishery);
        int ipp=terminal_tag_period(it,irr);
        int found_fishery=0;
        {
          //ofstream ofs("tagsur");
          //ofs << obstagcatch(it) << endl;
        }
        for (int fi=1;fi<=num_fish_incidents(irr,ipp);fi++)
        {
          if (parent(irr,ipp,fi)==fishery)
          {
            dvariable tsum=sum(obstagcatch(it,irr,ipp,fi));
            if (tsum<=0.0)
            {
              cerr << "Error in survival analysis tag times"
                   << " tsum=0 in routine "
                   << " fit_tag_returns_survival_analysis(void)" << endl;
              ad_exit(1);
            }
            //tag_return_probability(it)=
            //  sum(obstagcatch(it,irr,ipp,fi))/totpredtags;
            const dvar_vector& rtc=rep_rate(irr,ipp,fi)*tagcatch(it,irr,ipp,fi);
            tag_return_probability(it)=
              sum(rtc)/totpredtags+1.e-5;
            //ofs << " return prob = " << sum(rtc) << endl;
            f-=log(tag_return_probability(it));
  
            found_fishery=1;
            break;
          }
        }
        if (found_fishery==0)
        {
          cerr << "Error -- did not find a fishery corresponding to this"
               " returned tag " << endl;
        }
      }  
    }
  }
  //ofs.close();
  //exit(1);
  return f;
}

