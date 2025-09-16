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


void dvar_len_fish_stock_history::observed_tag_catch_calc(int direction_flag)
{
  ofstream ofs("tags.err");
  obstagcatch.initialize();
  int ignore_flag=0;
  int minj=nage;
  for (int it=1;it<=num_tag_releases;it++)
  {
    int isp=1;
    int offset=0.0;
    if (pmsd)
    {
      isp=pmsd->tag_species_pointer(it);
      offset=4*(isp-2);
    }


    imatrix& trl=tag_recaptures_by_length(it);
    int nrecaps=itr(it);
    for (int irr=1;irr<=nrecaps;irr++)
    {
      // get the fishing period and fishing incident for this recapture
      int fishery=trl(irr,2);   // assumes that fishery is the second field in
      // the *.tag file.
      int recap_year= trl(irr,3);
      int recap_month= trl(irr,4);

      if (recap_year >nyears)
      {
        ofs << "Ignoring tag from fishery " << fishery 
             << " in month " << recap_month
             << " in year " << recap_year << endl
         << " caught after last fishery" << endl;
        ofs << " from recapture " << irr << " in  tag release "
             << it << endl;
        cerr << "Ignoring tag caught after last fishery" << endl;
        cerr << " from recapture " << irr << " in  tag release "
             << it << endl;
      }
      else
      {
  
        int num_recaps=trl(irr,5);
        int il=int((trl(irr,1)-tag_shlen)/tag_filen)+1;
        if (il<1 || il > tag_nlint)
        {
          ofs << "error converting tag recapture data from length to age"
               << endl;
          cerr << "error converting tag recapture data from length to age"
               << endl;
          exit(1);
        }
        int nt=0;
        int tyear=0;
        int tmonth=0;
        int ir;
        int ip;
        for (nt=1;nt<=num_fish_times(fishery);nt++)
        {
          ir=realization_region(fishery,nt);
          ip=realization_period(fishery,nt);
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
              if (12*tyear+tmonth<12*recap_year+recap_month &&
                12*tyear1+tmonth1>12*recap_year+recap_month)
              {
                cerr << " Tag recapture " << irr << " from tag group "
                  << it << " was caught in year " << recap_year << " month "
                  << recap_month << endl
                  << " in fishery " << fishery << " but there was no"
                     " realization of this fishery at that time " << endl;

              }
              if (12*tyear+tmonth<=12*recap_year+recap_month &&
                12*tyear1+tmonth1>12*recap_year+recap_month)
              {
              if (initial_tag_period(it,ir)>ip)
              {
                cerr << " tyear = " << tyear << " tmonth = " << tmonth << endl;
                cerr << " tyear1 = " << tyear1 << " tmonth1 = " << tmonth1 << endl;
                ofs  << "Initial_tag_period for recapture " << irr
                     << " for tag group "
                     << it << " in region " << ir << endl
                     << " is greater than the tag period being "
                     << " assigned " << endl;
                ofs  << " I think this tag was caught before it was released "
                     << endl;
                cerr << "Initial_tag_period for recapture " << irr
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
          ofs << "tag recaptures for fishery " << fishery
             << " in month " << recap_month
             << " in year " << recap_year+year1-1 << endl
             << " do not match any fishery" << endl;
          cout << "tag recaptures for fishery " << fishery
             << " in month " << recap_month
             << " in year " << recap_year+year1-1 << endl
             << " do not match any fishery" << endl;
          ignore_flag=1;
        }
        if (!ignore_flag)
        {
          MY_DOUBLE_TYPE len=tag_shlen+(il-0.5L)*tag_filen;
          MY_DOUBLE_TYPE age;
          if (isp==1)
          {
//            age=value(age_at_length_calc(len,vb_coff,nage));
            age=value(age_at_length_calc(len,vb_coff,nage,parest_flags));	     //NMD_27Sep2018    
          }
          else
          {
            dvar_vector tvb=pmsd->vb_coff(isp);
//            age=value(age_at_length_calc(len,tvb,pmsd->nage(isp)));
            age=value(age_at_length_calc(len,tvb,pmsd->nage(isp),parest_flags));   //NMD_27Sep2018
          }
          int diff=tyear-tag_year(it);
          int fi=realization_incident(fishery,nt);
          if (age<=1)
          {
            int ind=1+diff;
            if (ind>nage) ind=nage;
            if (!age_flags(96))
              obstagcatch(it,ir,ip,fi,ind)+=num_recaps;
            else
            {
              if (ip<=terminal_tag_period(it,ir))
                obstagcatch(it,ir,ip,fi,ind)+=num_recaps;
              else
                pooledobstagcatch(ir,ip,fi,ind)+=num_recaps;
            }
          }

          else if (age+diff>=nage_by_tag_release(it) )
          {
            if (!age_flags(96))
              obstagcatch(it,ir,ip,fi,nage_by_tag_release(it))+=num_recaps;
            else
            {
              if (ip<=terminal_tag_period(it,ir))
                obstagcatch(it,ir,ip,fi,nage_by_tag_release(it))+=num_recaps;
              else
                pooledobstagcatch(ir,ip,fi,nage_by_tag_release(it))+=num_recaps;
            }
          }
          else
          {
            MY_DOUBLE_TYPE sf;
            int jj=int(age)+diff;
            /*
            if (f<=0.5L)
            {
              sf=2.*f*f;
            }
            else
            {
              sf=1.0-2.*square(1.-f);
            }
            */
            sf=daves_kludge(age);
            dvariable tp= sf*num_recaps;
            if (!age_flags(96))
            {
              obstagcatch(it,ir,ip,fi,jj)+=num_recaps;
              obstagcatch(it,ir,ip,fi,jj)-=tp;
              obstagcatch(it,ir,ip,fi,jj+1)+=tp;
            }
            else
            {
              if (ip<=terminal_tag_period(it,ir))
              {
                obstagcatch(it,ir,ip,fi,jj)+=num_recaps;
                obstagcatch(it,ir,ip,fi,jj)-=tp;
                obstagcatch(it,ir,ip,fi,jj+1)+=tp;
              }
              else
              {
                minj=min(jj,minj);
                pooledobstagcatch(ir,ip,fi,jj)+=num_recaps;
                pooledobstagcatch(ir,ip,fi,jj)-=tp;
                pooledobstagcatch(ir,ip,fi,jj+1)+=tp;
              }
            }
          }
        }
        ignore_flag=0;
      }
    }
  }
  // do we need this anyway?  dave f 9 feb 03
  /*
  if (minj<nage)
  {
    cerr << "Need to increase age_flags(96) to "
         << age_flags(96)+nage-minj << endl;
    exit(1);
  }
  */
}
