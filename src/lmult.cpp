/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#define HOME_VERSION
#include "all.hpp"

#ifndef __GNUC__
//#include <mf_menu.h>
#endif

  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
  extern int _NUMSV;
  const MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;

void klxx(dvector& u){;}

imatrix grouping_flags(const ivector& v)
{
  int i;
  int minv=min(v);
  int maxv=max(v);
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  ivector icount(minv,maxv);
  icount.initialize();
  for (i=mmin;i<=mmax;i++) icount(v(i))+=1;
  imatrix tmp(minv,maxv,1,icount);
  icount.initialize();
  for (i=mmin;i<=mmax;i++) 
  {
    icount(v(i))+=1;
    tmp(v(i),icount(v(i)))=i;
  }
  tmp.rowshift(1);
  return tmp;
}
  
    
    
    
    

  dvar_len_fish_stock_history::dvar_len_fish_stock_history
    (const dvar_fish_stock_history& fsh,
    int _nlint, MY_DOUBLE_TYPE& _shlen,MY_DOUBLE_TYPE& _filen,i3_array& len_shape,
    i3_array& wght_shape,
    i3_array& age_shape,int _nwint) :
    dvar_fish_stock_history((dvar_fish_stock_history&)fsh) ,
    nlintv(1,sum(num_fish_incidents)) ,
    //freq(1,num_fish_periods,1,num_fish_incidents) ,
    splinesel(1,num_fisheries) ,
    wsplinesel(1,num_fisheries) ,
    fmid(1,_nlint) ,
    //mean_length_yr(1,fsh.num_regions,1,fsh.nyears,1,nage),
    mean_length_yr_proj(1,12,1,nage),
    //mean_weight_yr(1,fsh.num_regions,1,fsh.nyears,1,nage),
    mean_length(1,fsh.num_regions,1,num_fish_periods,1,num_fish_incidents,
      1,nage),
    length_sel(1,num_fisheries,1,_nlint),
    relative_bio(1,fsh.num_regions,1,num_fish_periods,1,num_fish_incidents),
    len_dist(1,fsh.num_regions,1,num_fish_periods,1,_nlint) ,
    RB(1,fsh.num_regions,1,num_fish_periods,1,num_fish_incidents) ,
    ORB(1,fsh.num_regions,1,num_fish_periods,1,num_fish_incidents) ,
    vars(1,fsh.num_regions,1,num_fish_periods,1,num_fish_incidents,1,nage) ,
    global_vars(1,nage) ,
    sdevs(1,fsh.num_regions,1,num_fish_periods,1,num_fish_incidents,1,nage) ,
    vb_coff(1,4) ,
    var_coff(1,2) ,
    vb_bias(1,num_fisheries) ,
    common_vb_bias(0,num_fisheries) ,
    common_vb_bias_coffs(1,num_fisheries) ,
    tprob(1,fsh.num_regions,1,num_fish_periods,1,num_fish_incidents,1,_nlint) ,
    cpmature_at_length(1,_nlint),
    lengthbsel(1,num_fisheries,1,_nlint),
    wtprob(1,fsh.num_regions,1,num_fish_periods,1,num_fish_incidents,1,_nwint) ,
    age_sample_size(1,fsh.num_regions,1,num_fish_periods,1,num_fish_incidents) ,
    ifper(1,30) ,
    iffish(1,30) ,
    ifinc(1,30) ,
    iageclass(1,30) ,
    fmmin(1,30) ,
    fmmax(1,30) ,
    pdown(1,30) ,
    pup(1,30)
  {
    if(fsh.pmsd)
    {
      mean_length_yr.allocate(1,fsh.num_regions);
      int ng=nage;
      for (int ir=1;ir<=num_regions;ir++)
      {
        int  is=fsh.pmsd->region_species_pointer(ir);
        if (is>1) ng=fsh.pmsd->nage(is);
        mean_length_yr(ir).allocate(1,fsh.nyears,1,ng);
      }
    }
    else
    {
      mean_length_yr.allocate(1,fsh.num_regions,1,fsh.nyears,1,nage);
    }
    if (allocated(len_shape))   //NMD_jun10_19
    {    
      tcs=new tail_compression_stuff();
    }
    else
    {
      tcs=0;
    }                          //NMD_jun10_19
    if (allocated(wght_shape)) 
    {
      wtcs=new wght_tail_compression_stuff();
    }
    else
    {
      wtcs=0;
    }
    wght_eps.allocate(1,fsh.num_regions,
        1,num_fish_periods,1,num_fish_incidents);
    //len_sample_size.allocate(1,fsh.num_regions,1,num_fish_periods,1,
    //  num_fish_incidents); 
    //wght_sample_size.allocate(1,fsh.num_regions,1,num_fish_periods,1,
    //  num_fish_incidents);

    if (allocated(wght_shape)) 
    {
      wght_freq.allocate(1,fsh.num_regions,
        1,num_fish_periods,1,num_fish_incidents,1,wght_shape);

      /*
//      if (parest_flags(301>0))
      if (parest_flags(301)>0)     //  NMD_7jun-19
      {
        tc_wtprob.allocate(1,fsh.num_regions,
          1,num_fish_periods,1,num_fish_incidents);
        tc_wght_freq.allocate(1,fsh.num_regions,
          1,num_fish_periods,1,num_fish_incidents);
        max_wght_obs.allocate(1,fsh.num_regions,
          1,num_fish_periods,1,num_fish_incidents);
      }
      */
    }
    len_eps.allocate(1,fsh.num_regions,
      1,num_fish_periods,1,num_fish_incidents);

    if (allocated(len_shape))
    {
      len_freq.allocate(1,fsh.num_regions,
        1,num_fish_periods,1,num_fish_incidents,1,len_shape);

      /*
//      if (parest_flags(311>0))
      if (parest_flags(311)>0)     //  NMD_7jun-19
      {
        tc_tprob.allocate(1,fsh.num_regions,
          1,num_fish_periods,1,num_fish_incidents);
        tc_len_freq.allocate(1,fsh.num_regions,
          1,num_fish_periods,1,num_fish_incidents);
        max_len_obs.allocate(1,fsh.num_regions,
          1,num_fish_periods,1,num_fish_incidents);
      }
      */
    }

    if (allocated(age_shape)) age_freq.allocate(1,fsh.num_regions,
      1,num_fish_periods,1,num_fish_incidents,1,age_shape);

    growth_dev.allocate(1,fsh.nyears);
    shlen=_shlen;
    filen=_filen;
    nlint = _nlint;
    if (allocated(fmid)) fmid.fill_seqadd(shlen+.5*filen,filen);
  }

  dvar_len_fish_stock_history fishery_freq_record_array::
    get_history_data(int ntg,int num_regions,int nage, ivector& parest_flags,
    ivector& regmin,ivector& regmax,ivector& dataswitch,imatrix& Dflags,
    dvector& _region_area,int month_1, int _mfactor,int& _first_time,
    ivector& move_weeks,int direction_flag,imatrix & _season_region_flags,
    pmulti_species_data & pmsd)
  {
    fishery_header_record_array fhra(indexmin(),indexmax());
    int i;
    for (i=indexmin();i<=indexmax();i++)
    {
      fhra(i)=elem(i);
    }

   const dvar_fish_stock_history rfsh= 
	   fhra.get_history_data(ntg,num_regions,nage,
     parest_flags,regmin,regmax,dataswitch,Dflags,month_1,_mfactor,
     _first_time,move_weeks,direction_flag,_season_region_flags,pmsd);

     dvar_fish_stock_history fsh= 
  	   (dvar_fish_stock_history&) rfsh;
 
    i=1;
    i3_array len_shape;
        
    if (nlint)
    {
      len_shape.allocate(1,fsh.num_regions,1,fsh.num_fish_periods,
        1,fsh.num_fish_incidents);
      int ir;
      for (ir=1;ir<=fsh.num_regions;ir++)
      {
        for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)  
        {
          for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) 
          {     
            if (allocated(elem(i).freq))
            {
              if (sum(elem(i).freq)>0)
              {
                fsh.total_num_obs+=nlint;
                len_shape(ir,ip,fi)=nlint;
              }
              else
              {
                len_shape(ir,ip,fi)=1;
              }
            }
            else
            {
              len_shape(ir,ip,fi)=1;
            }
            i++;
          }
        }
      }
    }
    i=1;
    i3_array wght_shape;
        
    if (nwint)
    {
      int ir;
      wght_shape.allocate(1,fsh.num_regions,1,fsh.num_fish_periods,
        1,fsh.num_fish_incidents);
      for (ir=1;ir<=fsh.num_regions;ir++)
      {
        for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)  
        {
          for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) 
          {     
            if (allocated(elem(i).wfreq))
            {
              dvector& tmp= elem(i).wfreq;
              if (sum(elem(i).wfreq)>0)
              {
                fsh.total_num_obs+=nwint;
                wght_shape(ir,ip,fi)=nwint;
              }
              else
              {
                wght_shape(ir,ip,fi)=1;
              }
            }
            else
            {
              wght_shape(ir,ip,fi)=1;
            }
            i++;
          }
        }
      }
    }

    i=1;
    i3_array age_shape;
    if (fsh.age_nage)
    {
      int ir;
      age_shape.allocate(1,fsh.num_regions,1,fsh.num_fish_periods,
        1,fsh.num_fish_incidents);
      for (ir=1;ir<=fsh.num_regions;ir++)
      {
        for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)  
        {
          for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) 
          {     
            if (allocated(elem(i).age_freq))
            {
              dvector& tmp= elem(i).age_freq;
              if (sum(elem(i).age_freq)>0)
              {
                fsh.total_num_obs+=fsh.age_nage;
                age_shape(ir,ip,fi)=fsh.age_nage;
              }
              else
              {
                age_shape(ir,ip,fi)=-1;
              }
            }
            else
            {
              age_shape(ir,ip,fi)=-1;
            }
            i++;
          }
        }
      }
    }

    dvar_len_fish_stock_history lfsh(fsh,nlint,shlen,filen,len_shape,
      wght_shape,age_shape,nwint);

    lfsh.len_freq.initialize();
    i=1;
    int ir;
    for (ir=1;ir<=fsh.num_regions;ir++)
    {
      for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)  
      {
        for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) 
        {     
          if (allocated(elem(i).freq))
          {
            MY_DOUBLE_TYPE lss=sum(elem(i).freq);
#if !defined(NO_MY_DOUBLE_TYPE)
            if (lss<1.0L)
#else
            if (lss<1.0)
#endif
            {
              cerr << "This happened " << lss << endl;
              cerr << elem(i).freq << endl;
            }
            lfsh.len_sample_size(ir,ip,fi)=lss;
          }
          else
          {
            lfsh.len_sample_size(ir,ip,fi)=0;
          }

          if (lfsh.len_sample_size(ir,ip,fi))
          {
            lfsh.len_freq(ir,ip,fi)=elem(i).freq;
            {
              lfsh.len_freq(ir,ip,fi)=lfsh.len_freq(ir,ip,fi)/
                lfsh.len_sample_size(ir,ip,fi);
            }
          }
          i++;
        }
      }
    }
    if (nwint)
    {
      lfsh.wght_freq.initialize();
      i=1;
      for (ir=1;ir<=fsh.num_regions;ir++)
      {
        for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)  
        {
          for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) 
          {     
            if (allocated(elem(i).wfreq))
              lfsh.wght_sample_size(ir,ip,fi)=sum(elem(i).wfreq);
            else
              lfsh.wght_sample_size(ir,ip,fi)=0;

            if (lfsh.wght_sample_size(ir,ip,fi))
            {
              dvector& tmp1=lfsh.wght_freq(ir,ip,fi);
              dvector& tmp2=elem(i).wfreq;
              tmp1=tmp2;
              lfsh.wght_freq(ir,ip,fi)=elem(i).wfreq;
              {
                lfsh.wght_freq(ir,ip,fi)=lfsh.wght_freq(ir,ip,fi)/
                  lfsh.wght_sample_size(ir,ip,fi);
              }
            }
            i++;
          }
        }
      }
    }

    if (lfsh.age_nage)
    {
      lfsh.age_freq.initialize();
      i=1;
      for (ir=1;ir<=fsh.num_regions;ir++)
      {
        for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)  
        {
          for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) 
          {     
            if (allocated(elem(i).age_freq))
              lfsh.age_sample_size(ir,ip,fi)=sum(elem(i).age_freq);
            else
              lfsh.age_sample_size(ir,ip,fi)=0;


            if (lfsh.age_sample_size(ir,ip,fi))
            {
              lfsh.age_freq(ir,ip,fi)=elem(i).age_freq;
              {
                lfsh.age_freq(ir,ip,fi)=lfsh.age_freq(ir,ip,fi)/
                  lfsh.age_sample_size(ir,ip,fi);
              }
            }
            i++;
          }
        }
      }
    }
    lfsh.region_area.allocate(1,num_regions);
    lfsh.region_area=_region_area;
    lfsh.parest_flags=parest_flags;
    return lfsh;
  }


void get_global_fish_periods(dvar_fish_stock_history& fsh)
{
  int ir,ip;    
  dvector dates(1,fsh.num_regions);
  if (allocated(fsh.global_fishing_periods))
    fsh.global_fishing_periods.deallocate(); 
  fsh.global_fishing_periods.allocate(1,fsh.num_regions,
    1,fsh.num_fish_periods);
  ivector periods(1,fsh.num_regions);
  int current_period=1;
  periods=1;
  for (ir=1;ir<=fsh.num_regions;ir++)
  {
    dates(ir)=get_node_date_for_region_and_period(ir,periods(ir),fsh);
  }
  do
  {
    int rmin=min_index(dates);
    MY_DOUBLE_TYPE datemin=dates(rmin);
    if (datemin>1.e+60) break;

    for (ir=1;ir<=fsh.num_regions;ir++)
    {
      if (dates(ir)==datemin)
      {
        fsh.global_fishing_periods(ir,periods(ir))=current_period;
        periods(ir)++;
        if (periods(ir)>fsh.num_fish_periods(ir))
          dates(ir)=1.e+61;
        else
          dates(ir)=get_node_date_for_region_and_period(ir,periods(ir),fsh);
      }
    }
    current_period++;
  }
  while(1);
}

  
  void normalize_effort_data(dvar_fish_stock_history& fsh)
  {
    fsh.missing_effort_flag.initialize();
    MY_DOUBLE_TYPE effset;
    int i;
    int lastset;
    ivector ff29=column(fsh.fish_flags,29);
    ivector oldff29=column(fsh.old_fish_flags,29);
    fsh.average_effort.initialize();
    for (i=1;i<=fsh.num_fisheries;i++)
    {
      fsh.missing_effort_by_realization_flag(i)=1;
      MY_DOUBLE_TYPE ssum=0.0;
      int icount=0;
      int nt;
      for (nt=1;nt<=fsh.num_fish_times(i);nt++)
      {
        int rp=fsh.realization_period(i,nt);
        int ri=fsh.realization_incident(i,nt);
        int rr=fsh.realization_region(i,nt);
      
        if (fsh.effort(rr,rp,ri) == 0.0) 
        {
          if (fsh.num_fish_times(i)>fsh.num_real_fish_times(i))
          {
            // add a small number for 0 projections so that
            // log does not blow up
            fsh.effort(rr,rp,ri)=1.e-10;
          }
          else
          {
            cerr << "Zero effort for fishery " << i
               << " realization " << nt 
               << " region " << rr 
               << " fishing period " << rp
               << " fishing incident " << ri
               << endl;
            ad_exit(1);
          }
        }
     
        fsh.really_true_effort(rr,rp,ri)=fsh.effort(rr,rp,ri); 
        fsh.true_effort_by_fishery(i,nt)=fsh.effort(rr,rp,ri); 
#if !defined(NO_MY_DOUBLE_TYPE)
        if (fsh.effort(rr,rp,ri) > -0.5L) 
#else
        if (fsh.effort(rr,rp,ri) > -0.5) 
#endif
        {
          ssum+=fsh.true_effort_by_fishery(i,nt);
          icount++;
        }
      }
      if (icount>0)
      {
        fsh.average_effort(i)=ssum/icount;
      }
      
      for (nt=1;nt<=fsh.num_fish_times(i);nt++)
      {
        int rp=fsh.realization_period(i,nt);
        int ri=fsh.realization_incident(i,nt);
        int rr=fsh.realization_region(i,nt);
#if !defined(NO_MY_DOUBLE_TYPE)
        if (fsh.effort(rr,rp,ri) > -0.5L) 
#else
        if (fsh.effort(rr,rp,ri) > -0.5) 
#endif
        {
          fsh.true_effort_by_fishery(i,nt)/=fsh.average_effort(i); 
          fsh.log_true_effort_by_fishery(i,nt) 
            = log(fsh.true_effort_by_fishery(i,nt));
        }
        else
        {
          fsh.log_true_effort_by_fishery(i,nt)=0.0;
        }
      }
      //cout <<"lmult.cpp " << mean(fsh.log_true_effort_by_fishery(i))<<endl;
    }
    fsh.num_missing_effort_by_region.initialize();
    fsh.num_present_effort_by_region.initialize();
    for (i=1;i<=fsh.num_fisheries;i++)
    {
      effset=-1;
      lastset=0;
      for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
      {
        int rp=fsh.realization_period(i,nt);
        int ri=fsh.realization_incident(i,nt);
        int rr=fsh.realization_region(i,nt);
#if !defined(NO_MY_DOUBLE_TYPE)
        if (fsh.effort(rr,rp,ri) > -0.5L) 
#else
        if (fsh.effort(rr,rp,ri) > -0.5) 
#endif
        {
          // Here zero is "off" i.e. means no missing effort 
          fsh.missing_effort_by_region_flag(rr,rp,ri)=0;
          fsh.num_present_effort_by_region(rr,rp)++;
        }
        else
        {
          // Here 1 is "on" i.e. means missing effort 
          fsh.missing_effort_by_region_flag(rr,rp,ri)=1;
          fsh.num_missing_effort_by_region(rr,rp)++;
        }
 
#if !defined(NO_MY_DOUBLE_TYPE)
        if (fsh.effort(rr,rp,ri) > -0.5L) 
#else
        if (fsh.effort(rr,rp,ri) > -0.5) 
#endif
        {
          effset=fsh.effort(rr,rp,ri);
          for (int ii=lastset+1;ii<=nt-1;ii++)
          {
            int rp=fsh.realization_period(i,ii);
            int ri=fsh.realization_incident(i,ii);
            int rr=fsh.realization_region(i,ii);
            fsh.effort(rr,rp,ri)=effset;
          }
          lastset=nt;
        }
        else
        {
          fsh.missing_effort_flag(i)=1;
          // Here zero is "on" brcause it turns off parameters
          // in mewmau5a.cpp
          fsh.missing_effort_by_realization_flag(i,nt)=0;
          if (effset>0)
          {
            fsh.effort(rr,rp,ri)=effset;
            lastset=nt;
          }         
        }    
      }
      if (lastset==0)
      {
        for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
        {
          int rp=fsh.realization_period(i,nt);
          int ri=fsh.realization_incident(i,nt);
          int rr=fsh.realization_region(i,nt);
          fsh.effort(rr,rp,ri)=1;
        }
      }
    }

    int projflag=sum(fsh.data_fish_flags(4));
    if (!sum(ff29))
    {
      for (i=1;i<=fsh.num_fisheries;i++)
      {
        MY_DOUBLE_TYPE sum=0.;
        int icount=0;
        int nt;
        for (nt=1;nt<=fsh.num_fish_times(i);nt++)
        {
          int rp=fsh.realization_period(i,nt);
          int ri=fsh.realization_incident(i,nt);
          int rr=fsh.realization_region(i,nt);
          if (projflag)
          {
            int yr=fsh.year(rr,rp);
            int mn=fsh.month(rr,rp);
            int pyr=fsh.data_fish_flags(4,i);
            int pmn=fsh.data_fish_flags(5,i);
             
            if (yr>=pyr && mn >=pmn)
            {
              cout << yr << " " << mn << " " << pyr << " "
                  << pmn << "  " << endl;
              break;
            }
          }
          icount++;
          sum+=fsh.effort(rr,rp,ri);
        }
        if (icount) sum/=icount;
        if (sum>0)
        {
          for (nt=1;nt<=fsh.num_fish_times(i);nt++)
          {
            int rp=fsh.realization_period(i,nt);
            int ri=fsh.realization_incident(i,nt);
            int rr=fsh.realization_region(i,nt);
            fsh.effort_normalization_factor(rr,rp,ri)=sum;
            fsh.effort(rr,rp,ri)/=fsh.effort_normalization_factor(rr,rp,ri);
          }
        }
      }
    }
    else 
    {  
      dvector rsum(min(ff29),max(ff29));
      dvector isum(min(ff29),max(ff29));
      dvector average_grouped_region_area(min(ff29),max(ff29));
      rsum.initialize();
      isum.initialize();
      for (i=1;i<=fsh.num_fisheries;i++)
      {
        rsum(ff29(i))+=fsh.region_area(fsh.realization_region(i,1));
        isum(ff29(i))+=1;
      }
      average_grouped_region_area=elem_div(rsum,isum);
      if (!sum(oldff29) || !norm2(oldff29-ff29))
      {   //either first run or grouped the same as last time
        dvector ssum(min(ff29),max(ff29));
        dvector ntsum(min(ff29),max(ff29));
        dvector average_grouped_effort(min(ff29),max(ff29));
        ssum.initialize();
        ntsum.initialize();
        
        for (i=1;i<=fsh.num_fisheries;i++)
        {
          for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
          {
            int rp=fsh.realization_period(i,nt);
            int ri=fsh.realization_incident(i,nt);
            int rr=fsh.realization_region(i,nt);
            
            if (projflag)
            {
              int yr=fsh.year(rr,rp);
              int mn=fsh.month(rr,rp);
              int pyr=fsh.data_fish_flags(4,i);
              int pmn=fsh.data_fish_flags(5,i);
              
              if (yr>=pyr && mn >=pmn)
              {
                cout << yr << " " << mn << " " << pyr << " "
                    << pmn << "  " << endl;
                break;
              }
            }
            
            ntsum(ff29(i))++;
            ssum(ff29(i))+=fsh.effort(rr,rp,ri);
          }
        } 
        average_grouped_effort=elem_div(ssum,ntsum);
        //ofstream xcout("tt");
        dmatrix tmp(1,fsh.num_fisheries,1,fsh.num_fish_times);


        for (i=1;i<=fsh.num_fisheries;i++)
        {
          for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
          {
            int rp=fsh.realization_period(i,nt);
            int ri=fsh.realization_incident(i,nt);
            int rr=fsh.realization_region(i,nt);

            MY_DOUBLE_TYPE tmp1=fsh.effort(rr,rp,ri)/
              (fsh.region_area(rr)/rsum(ff29(i))*isum(ff29(i))
              *ssum(ff29(i))/ntsum(ff29(i)));


            fsh.effort_normalization_factor(rr,rp,ri)=
              fsh.region_area(rr)/rsum(ff29(i))*isum(ff29(i))
              *ssum(ff29(i))/ntsum(ff29(i));

            fsh.effort(rr,rp,ri)/=fsh.effort_normalization_factor(rr,rp,ri);

            if (fabs(fsh.effort(rr,rp,ri)-tmp1)>1.e-20)
            {
              cout << "ERROR: normalising effort data - out of bounds" << endl;
              cout << "Exiting..." << endl;
              ad_exit(1);
            }

            tmp(i,nt)=fsh.effort(rr,rp,ri);
            //xcout << fsh.effort(rr,rp,ri) << " ";
          }
          //xcout <<"lmult.cpp " << endl;
        }
        int mmin=min(ff29);
        int mmax=max(ff29);
        
       /*
        for (int ig=mmin;ig<=mmax;ig++)
        { 
          xcout << "Group " << ig << endl;
          for (i=1;i<=fsh.num_fisheries;i++)
          {
            
            if (ff29(i)==ig)
            {
              xcout << "fishery " << i << endl;
              xcout << max(tmp(i))<< " "  << mean(tmp(i)) 
                   << " " << sort(tmp(i)) << endl;
            }
          }
        }
       */
      }
      else
      {
        //grouping of catchability has been changed
        dvector ssum(min(ff29),max(ff29));
        dvector ntsum(min(ff29),max(ff29));
        dvar_vector oldssum(min(oldff29),max(oldff29));
        dvector oldntsum(min(oldff29),max(oldff29));
        oldssum.initialize();
        oldntsum.initialize();
        ssum.initialize();
        ntsum.initialize();
        for (i=1;i<=fsh.num_fisheries;i++)
        {
          for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
          {
            int rp=fsh.realization_period(i,nt);
            int ri=fsh.realization_incident(i,nt);
            int rr=fsh.realization_region(i,nt);
            ssum(ff29(i))+=fsh.effort(rr,rp,ri);
            oldssum(oldff29(i))+=fsh.effort(rr,rp,ri);
          }
          ntsum(ff29(i))+=fsh.num_fish_times(i);
          oldntsum(oldff29(i))+=fsh.num_fish_times(i);
        } 

        for (i=1;i<=fsh.num_fisheries;i++)
        {
          fsh.q0(i)*=ssum(ff29(i))/oldssum(oldff29(i))*oldntsum(oldff29(i))/
            ntsum(ff29(i));
        }

        //ofstream xcout("tt");
        for (i=1;i<=fsh.num_fisheries;i++)
        {
          for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
          {
            int rp=fsh.realization_period(i,nt);
            int ri=fsh.realization_incident(i,nt);
            int rr=fsh.realization_region(i,nt);
            fsh.effort(rr,rp,ri)/=
              fsh.region_area(rr)*ssum(ff29(i))/ntsum(ff29(i));
            //xcout << fsh.effort(rr,rp,ri) << " ";
          }
          //xcout <<"lmult.cpp " << endl;
        }
      }
    }
  }

  void get_effort_data_by_fishery(dvar_fish_stock_history& fsh)
  {
    for (int i=1;i<=fsh.num_fisheries;i++)
    {
      MY_DOUBLE_TYPE meanf=0.0;
      int icount=0.0;
      for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
      {
        int rp=fsh.realization_period(i,nt);
        int ri=fsh.realization_incident(i,nt);
        int rr=fsh.realization_region(i,nt);
        if (fsh.effort(rr,rp,ri) > 0) 
        {
          fsh.effort_by_fishery(i,nt)=fsh.effort(rr,rp,ri);
          fsh.log_effort_by_fishery(i,nt)=log(fsh.effort(rr,rp,ri));
          meanf+=fsh.log_effort_by_fishery(i,nt);
          icount ++;
        }
        else
        {
          if (fsh.age_flags(92)>0)
          {
            cerr << "Error -- can not deal with missing effort in"
              " catch-conditioned option yet" << endl;
            ad_exit(1);
          }
        }
      }
      meanf/=icount;
      fsh.normalized_log_effort_by_fishery(i)=
        fsh.log_effort_by_fishery(i)-meanf;
    }
  }

  void dvar_fish_stock_history::process_effort_weight_data(void)
  {
    int i,nt;
    ivector ff66=column(fish_flags,66);  // fish flags for effort weights
    ff66.initialize(); // set all flags to 0
    for (i=1;i<=num_fisheries;i++)
    {
      for (nt=1;nt<=num_fish_times(i);nt++)
      {
        int rp=realization_period(i,nt);
        int ri=realization_incident(i,nt);
        int rr=realization_region(i,nt);
        if (parest_flags(225)>0)
        {
#if !defined(NO_MY_DOUBLE_TYPE)
          if (effort_weight(rr,rp,ri) < -0.5L) 
#else
          if (effort_weight(rr,rp,ri) < -0.5) 
#endif
            effort_weight(rr,rp,ri)=0.1; 
        }
  
        effort_weight_by_fishery(i,nt)=effort_weight(rr,rp,ri); 
#if !defined(NO_MY_DOUBLE_TYPE)
        if (effort_weight(rr,rp,ri) > -0.5L) 
#else
        if (effort_weight(rr,rp,ri) > -0.5) 
#endif
        {
          ff66(i)=1;
        }
      }
    }
      
    for (i=1;i<=num_fisheries;i++)
    {
      // now check effort weights for consistency
      if (ff66(i)==1)  // can't have effort_weight = -1
      {
        for (nt=1;nt<=num_fish_times(i);nt++)
        {
          int rp=realization_period(i,nt);
          int ri=realization_incident(i,nt);
          int rr=realization_region(i,nt);
#if !defined(NO_MY_DOUBLE_TYPE)
          if (effort(rr,rp,ri) > -0.5L) 
#else
          if (effort(rr,rp,ri) > -0.5) 
#endif
          {
#if !defined(NO_MY_DOUBLE_TYPE)
            if (effort_weight(rr,rp,ri) < -0.5L) 
#else
            if (effort_weight(rr,rp,ri) < -0.5) 
#endif
            {
              cerr << "Consistency error in effort weights for fishery "
                   << i << " some but not all effort weights are != -1"
                   << endl  << " value is " << effort_weight(rr,rp,ri)  
                   << " year = " << really_true_year(rr,rp)+year1-1 
                   << " month = " << really_true_month(rr,rp) << endl;
              ad_exit(1);
            }
          }
        }
      }
    }
  }
  
  void dvar_fish_stock_history::set_sel_seasons(void)
  {
    ivector ff74=column(fish_flags,74);

    for (int fi=1;fi<=num_fisheries;fi++)
    {
      if (ff74(fi)>1)
      {
        for (int i=1;i<=num_fish_times(fi);i++)
        {
          int ir=realization_region(fi,i);
          int ip=realization_period(fi,i);
          int ri=realization_incident(fi,i);

          sel_seasons(ir,ip,ri)=(really_true_month(ir,ip)-1)%ff74(fi)+1;
//          cout << sel_seasons(ir,ip,ri) << " ";
        }
//        cout << endl;
      }
    }
  }

   
