/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"
#include <admodel.h>

#ifdef __MSC__
  MY_DOUBLE_TYPE mfexp(const MY_DOUBLE_TYPE& );
  dvariable mfexp(const prevariable& );
  dvar_vector mfexp(const dvar_vector& );
  dvector mfexp(const dvector& );
#else
  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
  dvariable mfexp(BOR_CONST prevariable& );
  //dvar_vector mfexp(dvar_vector& );
  //dvector mfexp(dvector& );
#endif

void dvar_len_fish_stock_history::tag_returns_report(const ofstream & ofs)
{
  ivector group_flags=column(fish_flags,44);
  int gsum44=sum(group_flags);
  int gmax44=Max(group_flags);
  dvariable f=0.0;
  dvariable gp_pen=0.0;
  int fi;
  dvar_vector& rep_rate1=fish_pars(3);
  ivector gp_fish44(1,gmax44);
  dvar_matrix obsgroupedcatch;
  dvar_matrix groupedcatch;
  ivector iflag(1,gmax44);
  ivector vsize;
  if (!gsum44)
    vsize.allocate(1,num_fisheries);
  else
    vsize.allocate(1,gmax44);
  vsize.initialize();
  if (gmax44)
  {
    gp_fish44.initialize();
    for (fi=num_fisheries;fi>=1;fi--)
      gp_fish44(group_flags(fi))=fi;
    obsgroupedcatch.allocate(1,gmax44,1,nage);
    groupedcatch.allocate(1,gmax44,1,nage);
  }
  else
  {
    obsgroupedcatch.allocate(1,num_fisheries,1,nage);
    groupedcatch.allocate(1,num_fisheries,1,nage);
  }
  // count up size of vectors
  int it;
  for (it=1;it<=num_tag_releases;it++)
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      int ub;
      if (!age_flags(96))
        ub=num_fish_periods(ir);
      else
        ub=terminal_tag_period(it,ir);
      for (int ip=initial_tag_period(it,ir);ip<=ub;ip++)
      {
        for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {
          if (!gsum44)
          {
            // don't group fisheries
            int pp1=parent(ir,ip,fi);
            const dvar_vector& rtc=rep_rate(ir,ip,fi)*tagcatch(it,ir,ip,fi);
            vsize(pp1)+=size_count(rtc);
            //obstagcatch(it,ir,ip,fi)
          }
          else
          {
            iflag.initialize();
            ivector& pi=parent(ir,ip);
            for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
            {
              int pp1=pi(fi);
              int pp=group_flags(pp1);
              iflag(pp)=1;
            }

            for (fi=1;fi<=gmax44;fi++)
            {
           //   cout << "Group " << fi << endl;
              if(iflag(fi))
              {
                vsize(fi)+=size_count(groupedcatch(fi));
                //obsgroupedcatch(fi)
                // -groupedcatch(fi)),
              }
            } 
          }
        }
      }
    }
  }
  d3_array all_predicted_and_observed_tags(1,vsize.indexmax(),1,vsize,1,2);
  ivector ioffset(1,vsize.indexmax());
  ioffset=1;

  
  for (it=1;it<=num_tag_releases;it++)
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      int ub;
      if (!age_flags(96))
        ub=num_fish_periods(ir);
      else
        ub=terminal_tag_period(it,ir);
      for (int ip=initial_tag_period(it,ir);ip<=ub;ip++)
      {
        for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {
          if (!gsum44)
          {
            // don't group fisheries
            int pp1=parent(ir,ip,fi);
            const dvar_vector& rtc=rep_rate(ir,ip,fi)*tagcatch(it,ir,ip,fi);
            //obstagcatch(it,ir,ip,fi)
            for (int j=1;j<=nage;j++)
            {
              all_predicted_and_observed_tags(pp1,ioffset(pp1),1)=value(rtc(j));
              all_predicted_and_observed_tags(pp1,ioffset(pp1)++,2)
                = value(obstagcatch(it,ir,ip,fi,j));
            }
          }
          else
          {
            obsgroupedcatch.initialize();
            groupedcatch.initialize();
            iflag.initialize();
            ivector& pi=parent(ir,ip);
            for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
            {
              int pp1=pi(fi);
              int pp=group_flags(pp1);
              const dvar_vector& rtc=rep_rate(ir,ip,fi)*tagcatch(it,ir,ip,fi);
              iflag(pp)=1;
              obsgroupedcatch(pp)+=obstagcatch(it,ir,ip,fi); 
              groupedcatch(pp)+=rtc;
            }

            for (fi=1;fi<=gmax44;fi++)
            {
              cout << "Group " << fi << endl;
              if(iflag(fi))
              {
                    //obsgroupedcatch(fi)
                    // -groupedcatch(fi)),
                for (int j=1;j<=nage;j++)
                {
                  all_predicted_and_observed_tags(fi,ioffset(fi),1)
                     =value(groupedcatch(fi,j));
                  all_predicted_and_observed_tags(fi,ioffset(fi)++,2)
                    = value(obsgroupedcatch(fi,j));
                }
              }
            } 
          }
        }
      }
    }
  }
  ofstream ofs1("tagtmp");
  ofs1 << " tags " << endl;
  for (fi=1;fi<=vsize.indexmax();fi++)
  {
    ofs1 << "#  fishery (group) "  << fi << endl;
    ofs1 << sort(all_predicted_and_observed_tags(fi),1) << endl;
  }
}

