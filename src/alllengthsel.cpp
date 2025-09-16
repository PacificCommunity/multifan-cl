/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"
static void uuuu(void) {;}

dvar_matrix dvar_len_fish_stock_history::length_dist_calcs
(dvar_vector& mean_len,dvar_vector& sigg)
{
  const MY_DOUBLE_TYPE a2 = 35.91908881;
  const MY_DOUBLE_TYPE a3 = 785.858644;
  const MY_DOUBLE_TYPE a22 = 2*35.91908881;
  const MY_DOUBLE_TYPE a33 = 3*785.858644;
  dvariable prob=0.0;
  int js=1;
  int break_flag=0;
  int i;
  dvar_matrix r(1,nage,1,nlint);
  r.initialize();
  for (i=1; i<=nlint; i++) 
  {
    MY_DOUBLE_TYPE& fmidd=fmid(i);
   
    for (int j=1; j<=nage; j++) 
    {
      dvariable t=(fmidd-mean_len(j))/sigg(j);
      if ( fabs(t) <= 3.03)
      {
        if ( fabs(t) <= 3.00)
          prob=exp(-t*t/2.e0)/sigg(j);
        else if (t > 0.e0)
        {
          dvariable u=t-3.03;
          prob=(u*u*(a2+a3*u))/sigg(j);
        }
        else if (t < 0.e0)
        {
          dvariable u=t+3.03;
          prob=(u*u*(a2-a3*u))/sigg(j);
        }
        r(j,i) =  prob;
      }
    }
  }
  {
//   ofstream ofs("rvals");
//   ofs << r << endl;
//   cout << "here" << endl;
  }
  return r;
}


dvar_matrix dvar_len_fish_stock_history::weight_dist_calcs
(dvar_vector& mean_len,dvar_vector& sigg)
{
  const MY_DOUBLE_TYPE a2 = 35.91908881;
  const MY_DOUBLE_TYPE a3 = 785.858644;
  const MY_DOUBLE_TYPE a22 = 2*35.91908881;
  const MY_DOUBLE_TYPE a33 = 3*785.858644;
  dvariable prob=0.0;
  int js=1;
  int break_flag=0;
  int i;
  dvar_matrix r(1,nage,1,nwint);
  r.initialize();
  cout << r.indexmin() << endl;
  cout << r.indexmax() << endl;
  //  wmid=pow(realwmid/value(sv(27)),1.0/value(sv(28)));
  for (i=1; i<=nwint; i++) 
  {
    //dvariable fmidd=twmid(i);
    dvariable fmidd=wmid(i);
   
    for (int j=1; j<=nage; j++) 
    {
      dvariable t=(fmidd-mean_len(j))/sigg(j);
      if ( fabs(t) <= 3.03)
      {
        if ( fabs(t) <= 3.00)
          prob=exp(-t*t/2.e0)/sigg(j);
        else if (t > 0.e0)
        {
          dvariable u=t-3.03;
          prob=(u*u*(a2+a3*u))/sigg(j);
        }
        else if (t < 0.e0)
        {
          dvariable u=t+3.03;
          prob=(u*u*(a2-a3*u))/sigg(j);
        }
        r(j,i) =  prob;
      }
      if (parest_flags(93)<2)
      {
        r(j,i) += 1.e-6;  // NMD_14feb2023
      }
    }
  }

  for (int j=1; j<=nage; j++) 
  {
    if (sum(r(j))<1.e-8)
    {
// NMD_14feb2023
// - check for zeroes in all size intervals
      if (parest_flags(93)<2)
      {
        r(j)(1,nwint)=1.e-6;
      }
//       cout << "Error: sum(r_ji) matrix over age is zero: "
//            << r(j) << "  age: " << j << endl;
    }
    r(j)=elem_div(r(j),wm2);
    r(j)/=sum(r(j));
  }
  {
//    ofstream ofs("rvals");
//    ofs << r << endl;
//    cout << "here" << endl;
  }
  return r;
}

void dvar_len_fish_stock_history::all_weight_dist_calcs_long
  (void)
{
  //dvar_matrix etsel=mfexp(tsel);
  ivector ff75=column(fish_flags,75);
  dvar4_array r(1,12,1,4,1,nage,1,nwint);
  r.initialize();
  imatrix month_week_weight_flag(1,12,1,4);
  i3_array month_week_fishery_weight_flag(1,12,1,4,1,num_fisheries);
  month_week_weight_flag.initialize();
  month_week_fishery_weight_flag.initialize();
  if (!allocated(age_weight_fishery_sel))
  {
    //age_weight_fishery_sel.allocate(1,12,1,4,1,num_fisheries,1,nage,1,nwint);
    age_weight_fishery_sel.allocate(1,12,1,4,1,num_fisheries,1,nage,1,0);
  }
  age_weight_fishery_sel.initialize();

  if (!allocated(age_weight_fishery_size_dist))
  {
    //age_weight_fishery_size_dist.allocate(1,12,1,4,1,num_fisheries,1,nage,1,nwint);
    age_weight_fishery_size_dist.allocate(1,12,1,4,1,num_fisheries,1,nage,
      1,0);
  }
  age_weight_fishery_size_dist.initialize();

  if (!allocated(age_fishery_sel_from_weight))
  {
    age_fishery_sel_from_weight.allocate(1,12,1,4,1,num_fisheries,1,nage);
  }

  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
    {
      int mn=month(ir,ip);
      int wk=week(ir,ip);
      if (month_week_weight_flag(mn,wk)==0)
      {
        month_week_weight_flag(mn,wk)=1;
        dvar_vector& mean_len=mean_length(ir,ip,1);
        dvar_vector& sigg=sdevs(ir,ip,1);
        //const dvar_vector& _mean_len=dvar_vector(value(mean_length(ir,ip,1)));
        //const dvar_vector& _sigg=dvar_vector(value(sdevs(ir,ip,1)));
        //dvar_vector& mean_len=(dvar_vector&)(_mean_len);
        //dvar_vector& sigg=(dvar_vector&)(_sigg);
        r(mn,wk)=weight_dist_calcs(mean_len,sigg);
      }
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)  // Loop over fishing
      {
        int pi=parent(ir,ip,fi);
        int is=sseason(ir,ip,fi);   //NMD_20dec2022
        int ib=bblock(ir,ip,fi);
//        if (allocated(wsplinesel(pi)))
        if (allocated(bswtempsel(pi,is,ib)))   //NMD_20dec2022
        {
          if (month_week_fishery_weight_flag(mn,wk,pi)==0)
          {
            dvar_vector ebswtempsel(1,bswtempsel(pi,is,ib).indexmax());
            ebswtempsel.initialize();
            ebswtempsel=exp(bswtempsel(pi,is,ib));
            month_week_fishery_weight_flag(mn,wk,pi)=1;
            for (int j=1;j<=nage;j++)
            {
              // this is the selectivity at length for each age class
//              age_weight_fishery_sel(mn,wk,pi,j)=
//                elem_prod(r(mn,wk,j),wsplinesel(pi));
              age_weight_fishery_sel(mn,wk,pi,j)=
                elem_prod(r(mn,wk,j),ebswtempsel);
              // mean selectivity at age after 
              age_fishery_sel_from_weight(mn,wk,pi,j)=
                  1.e-20+sum(age_weight_fishery_sel(mn,wk,pi,j));
              // length distribution of age j fish in the catch of
              //  fishery pi at this month and week
              age_weight_fishery_size_dist(mn,wk,pi,j)=
                age_weight_fishery_sel(mn,wk,pi,j)/
                  age_fishery_sel_from_weight(mn,wk,pi,j);
            }
          }
        }
      }
    }
  }
// NMD_21Feb2023  - assign selwt to all incidents for fisheries with WF
  ivector fshry_wtflag(1,num_fisheries);
  fshry_wtflag.initialize();
  for (int pi=1;pi<=num_fisheries;pi++)
  {
    for (int ii=1;ii<=num_fish_times(pi);ii++)
    {
      int ir=realization_region(pi,ii);
      int ip=realization_period(pi,ii);
      int fi=realization_incident(pi,ii);
      if (wght_sample_size(ir,ip,fi)>0.0)    //NMD_20jan2023
      {
        fshry_wtflag(pi)=1;
        break;
      }
    }
  }
  for (int pi=1;pi<=num_fisheries;pi++)
  {
    for (int ii=1;ii<=num_fish_times(pi);ii++)
    {
      int ir=realization_region(pi,ii);
      int ip=realization_period(pi,ii);
      int fi=realization_incident(pi,ii);
      int mn=month(ir,ip);
      int wk=week(ir,ip);
      if (parest_flags(94)==1)
      {
        if (fshry_wtflag(pi)>0)    //NMD_20jan2023
        {
//          incident_sel(ir,ip,fi)=age_fishery_sel_from_weight(mn,wk,pi);
          if (!ff75(pi))
          {
            incident_sel(ir,ip,fi)=age_fishery_sel_from_weight(mn,wk,pi);
          }
          else
          {
            incident_sel(ir,ip,fi)(ff75(pi)+1,nage)=
              age_fishery_sel_from_weight(mn,wk,pi)(ff75(pi)+1,nage);
            incident_sel(ir,ip,fi)(1,ff75(pi))=-20.;
          }
        }
      }
    }
  }
}

void dvar_len_fish_stock_history::all_length_dist_calcs_long
  (void)
{
  //dvar_matrix etsel=mfexp(tsel);
  ivector ff75=column(fish_flags,75);
  dvar4_array r(1,12,1,4,1,nage,1,nlint);
  r.initialize();
  imatrix month_week_length_flag(1,12,1,4);
  i3_array month_week_fishery_length_flag(1,12,1,4,1,num_fisheries);
  month_week_length_flag.initialize();
  month_week_fishery_length_flag.initialize();
  if (!allocated(age_length_fishery_sel))
  {
    age_length_fishery_sel.allocate(1,12,1,4,1,num_fisheries,
      1,nage,1,0);
      //1,nage,1,nlint);
  }
  age_length_fishery_sel.initialize();

  if (!allocated(age_length_fishery_size_dist))
  {
    age_length_fishery_size_dist.allocate(1,12,1,4,1,num_fisheries,
      1,nage,1,0);
      //1,nage,1,nlint);
  }
  age_length_fishery_size_dist.initialize();

  if (!allocated(age_fishery_sel_from_length))
  {
    age_fishery_sel_from_length.allocate(1,12,1,4,1,num_fisheries,1,nage);
  }

  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
    {
      int mn=month(ir,ip);
      int wk=week(ir,ip);
      if (month_week_length_flag(mn,wk)==0)
      {
        month_week_length_flag(mn,wk)=1;
        dvar_vector& mean_len=mean_length(ir,ip,1);
        dvar_vector& sigg=sdevs(ir,ip,1);
        r(mn,wk)=length_dist_calcs(mean_len,sigg);
      }
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)  // Loop over fishing
      {
        int pi=parent(ir,ip,fi);
        int is=sseason(ir,ip,fi);   //NMD_20dec2022
        int ib=bblock(ir,ip,fi);
//        if (allocated(splinesel(pi)))
        if (allocated(bstempsel(pi,is,ib)))   //NMD_20dec2022
        {
          if (month_week_fishery_length_flag(mn,wk,pi)==0)
          {
            dvar_vector ebstempsel(1,bstempsel(pi,is,ib).indexmax());
            ebstempsel.initialize();
            ebstempsel=exp(bstempsel(pi,is,ib));
            month_week_fishery_length_flag(mn,wk,pi)=1;
            for (int j=1;j<=nage;j++)
            {
              // this is the selectivity at length for each age class
              age_length_fishery_sel(mn,wk,pi,j)=
                elem_prod(r(mn,wk,j),ebstempsel);
//                elem_prod(r(mn,wk,j),splinesel(pi));
              // mean selectivity at age after 
              age_fishery_sel_from_length(mn,wk,pi,j)=
                  sum(age_length_fishery_sel(mn,wk,pi,j));
              // length distribution of age j fish in the catch of
              //  fishery pi at this month and week
              age_length_fishery_size_dist(mn,wk,pi,j)=
                age_length_fishery_sel(mn,wk,pi,j)/
                  age_fishery_sel_from_length(mn,wk,pi,j);
            }
          }
          if (!ff75(pi))
          {
            incident_sel(ir,ip,fi)=age_fishery_sel_from_length(mn,wk,pi);
          }
          else
          {
            incident_sel(ir,ip,fi)(ff75(pi)+1,nage)=
              age_fishery_sel_from_length(mn,wk,pi)(ff75(pi)+1,nage);
            incident_sel(ir,ip,fi)(1,ff75(pi))=-20.;
          }
        }
      }
    }
  }
}
