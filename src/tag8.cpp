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
  dvariable mfexp(BOR_CONST prevariable& );
  //dvar_vector mfexp(dvar_vector& );
  //dvector mfexp(dvector& );
#endif

#ifdef __MSVC32__
  dvariable age_at_length_calcxx(MY_DOUBLE_TYPE& v,dvar_vector& gml,int nslots);
#else
  dvariable age_at_length_calcxx(const dvariable& v,dvar_vector& gml,int nslots);
#endif

void dvar_len_fish_stock_history::observed_tag_catch_by_length_calc(void)
{
  obstagcatch_by_length.initialize();
  int ignore_flag=0;
  for (int it=1;it<=num_tag_releases;it++)
  {
    imatrix& trl=tag_recaptures_by_length(it);
    int nrecaps=itr(it);
    for (int itr=1;itr<=nrecaps;itr++)
    {
      // get the fishing period and fishing incident for this recapture
      int fishery=trl(itr,2);   // assumes that fishery is the second field in
      // the *.tag file.
      int recap_year= trl(itr,3);
      if (recap_year >nyears)
      {
        cerr << "Ignoring tag caught after last fishery" << endl;
        cerr << " from recapture " << itr << " in  tag release "
             << it << endl;
      }
      else
      {
        int recap_month= trl(itr,4);
        int num_recaps=trl(itr,5);
        int il=int((trl(itr,1)-tag_shlen)/tag_filen)+1;
        if (il<1 || il > tag_nlint)
        {
          cerr << "error converting tag recapture data from length to age"
               << endl;
          exit(1);
        }
        int nt=0;
        int tyear=0;
        int tmonth=0;
        int ir;
        int ip;
        int fi;
        for (nt=1;nt<=num_fish_times(fishery);nt++)
        {
          ir=realization_region(fishery,nt);
          ip=realization_period(fishery,nt);
          fi=realization_incident(fishery,nt);
          tyear=year(ir,ip);
          tmonth=month(ir,ip);
            if (tyear==recap_year && tmonth==recap_month) break;
            // add this for the moment to deal with incorrectly
            // grouped data
            if (nt <num_fish_times(fishery))
            {
            int ir1=realization_region(fishery,nt+1);
            int ip1=realization_period(fishery,nt+1);
            int tyear1=year(ir1,ip1);
            int tmonth1=month(ir1,ip1);
              int tm=12*recap_year+recap_month;
              if (12*tyear+tmonth<=12*recap_year+recap_month &&
                12*tyear1+tmonth1>12*recap_year+recap_month)
              {
              if (initial_tag_period(it,ir)>ip)
              {
                cerr << "Initial_tag_period for recapture " << itr
                     << " for tag group "
                     << it << " in region " << ir << endl
                     << " is greater than the tag period being "
                     << " assigned " << endl;
                cerr << " I think this tag was caught before it was released "
                     << endl;
                ignore_flag=1;
              }
              break;
              }
          }
        }
        if (nt>num_fish_times(fishery))
        {
          cout << "tag recaptures for fishery " << fishery
             << " in month " << recap_month
             << " in year " << recap_year+year1-1 << endl
             << " do not match any fishery" << endl;
          ignore_flag=1;
        }
        if (!ignore_flag)
        {
          if (!age_flags(96))
            obstagcatch_by_length(it,ir,ip,fi,il)+=num_recaps;
          else
          {
            if (ip<=terminal_tag_period(it,ir))
              obstagcatch_by_length(it,ir,ip,fi,il)+=num_recaps;
            else
              pooledobstagcatch_by_length(ir,ip,fi,il)+=num_recaps;
          }
        }
        ignore_flag=0;
      }
    }
  }
}
