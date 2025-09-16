/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

void dvar_fish_stock_history::test_tag_report(void * pq_flag)
{
  // shape of tagcatch
  //tagcatch.allocate(1,num_tag_releases,1,num_regions,initial_tag_period,
  //  terminal_tag_period,
  //   1,num_tagfish_incidents,min_tag_age4,nage);
  ofstream * pofs=0;
  if (pq_flag==0)
  {
    pofs = new ofstream("temporary_tag_report");
  }
  else
  {
    pofs=new ofstream("temporary_tag_report_q0");
  }
  if (!pmsd)
  {
    for (int it=1;it<=num_tag_releases;it++)
    {
      (*pofs) << "release group " << it << endl 
              << initial_tag_release_by_age(it) << endl;
    }
  }
  else
  {
    for (int it=1;it<=num_tag_releases;it++)
    {
      (*pofs) << "release group " <<  " " 
              << pmsd->tag_species_index(it,1) << "  species "
              << pmsd->tag_species_pointer(it) << " " << endl
              << initial_tag_release_by_age(it) << endl;
    }
  }

  if (!pmsd)
  {
    for (int it=1;it<=num_tag_releases;it++)
    {
      (*pofs) << "release group " << it << endl;
      for (int ir=1;ir<=num_regions;ir++)
      {
        (*pofs) << "predicted recapture in region  " << ir << endl;
        int mmin=tagcatch(it,ir).indexmin();
        int mmax=tagcatch(it,ir).indexmax();
        for (int ip=mmin;ip<=mmax;ip++)
        {
          int imax=tagcatch(it,ir,ip).indexmax();
          for (int fi=1;fi<=imax;fi++)
          {
            dvector rtc;
            if (age_flags(198))
                  rtc=value(tag_rep_rate(it,ir,ip,fi))*
                    value(tagcatch(it,ir,ip,fi));
                else
                  rtc=value(rep_rate(ir,ip,fi))*
                    value(tagcatch(it,ir,ip,fi));
  
            (*pofs) << "fishery " << parent(ir,ip,fi) 
                << "  year " << "  "  << really_true_year(ir,ip)  
                << "  month " << "  "  << really_true_month(ir,ip)  << endl;
            //(*pofs) << sum(rtc) <<  "  "  
            //    <<  sum(obstagcatch(it,ir,ip,fi))   << endl;
            (*pofs) << rtc <<  endl  
                <<  obstagcatch(it,ir,ip,fi)   << endl;
          }
        }
      }
    }
  }
  else
  {
    for (int it=1;it<=num_tag_releases;it++)
    {
      (*pofs) << "release group " <<  " " 
              << pmsd->tag_species_index(it,1) << "  species "
              << pmsd->tag_species_pointer(it) << " " << endl;
      for (int ir=1;ir<=num_regions;ir++)
      {
        (*pofs) << "predicted recapture in region  " << ir << endl;
        int mmin=tagcatch(it,ir).indexmin();
        int mmax=tagcatch(it,ir).indexmax();
        for (int ip=mmin;ip<=mmax;ip++)
        {
          int imax=tagcatch(it,ir,ip).indexmax();
          for (int fi=1;fi<=imax;fi++)
          {
            dvector rtc;
            if (age_flags(198))
                  rtc=value(tag_rep_rate(it,ir,ip,fi))*
                    value(tagcatch(it,ir,ip,fi));
                else
                  rtc=value(rep_rate(ir,ip,fi))*
                    value(tagcatch(it,ir,ip,fi));
  
            (*pofs) << "fishery " << parent(ir,ip,fi) 
                << "  year " << "  "  << really_true_year(ir,ip)  
                << "  month " << "  "  << really_true_month(ir,ip)  << endl;
            //(*pofs) << sum(rtc) <<  "  "  
            //    <<  sum(obstagcatch(it,ir,ip,fi))   << endl;
            (*pofs) << rtc <<  endl  
                <<  obstagcatch(it,ir,ip,fi)   << endl;
          }
        }
      }
    }
  }
  delete pofs;
  pofs=0;
}


    


