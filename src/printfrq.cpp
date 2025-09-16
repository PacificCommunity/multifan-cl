/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"

  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
 long int diff(long int&,const long int&);

void dvar_len_fish_stock_history::print_pred_frequencies(ofstream& ofs)
{
// Print version number and pointers     //NMD_13Jan2017
  if(parest_flags(200) >= 1052)
  {
    ofs << "# FIT 3" << endl;
  }  //NMD 13Jan2017

  //d3_array tprob(1,num_fish_periods, 1,num_fish_incidents, 1, nlint);
  tprob.initialize();
  //tprob+=1.e-5;
  int i;
  int j;
  dvariable temp;
  dvariable fdiff2;
  dvariable prob;
  dvariable u;
  MY_DOUBLE_TYPE a2 = 35.91908881;
  MY_DOUBLE_TYPE a3 = 785.858644;

  dvar_vector rnuml(1,nlint);
  ivector num_freqs(1,num_fisheries);
  for (i=1;i<=num_fisheries;i++)
  {
    int ii=0;
    for (j=1;j<=num_fish_times(i);j++)
    {
      int ir=realization_region(i,j);
      int ip=realization_period(i,j);
      int fi=realization_incident(i,j);
      if (len_sample_size(ir,ip,fi)>0)
      ii++;
    }
    num_freqs(i)=ii;
  }
  dvector dt(1,10);
  dt.initialize();
  dt(1)=0;
  dt(3)=20;
  ofs << dt << endl;

  ofs << nlint <<  " " << shlen << " " << filen << endl;
  ofs << num_fisheries+1 << endl;
  ofs << num_freqs << " " << num_fisheries << endl;
  ofs << nage << endl;

// Print pointer for species in each fishery
  if(parest_flags(200) > 1051)   //NMD 13Jan2017
  {
    ivector fishery_species_pointer(1,num_fisheries);
    if(!pmsd)
    {
      fishery_species_pointer=1;
      ofs << "# Fishery species pointer " << endl << "  " 
          << fishery_species_pointer  << endl;
    } else {
      int mmin=fishery_regions.indexmin();
      int mmax=fishery_regions.indexmax();
      for(int ir=mmin;ir<=mmax;ir++)
      {
        int ireg=fishery_regions(ir);
        fishery_species_pointer(ir)=pmsd->region_species_pointer(ireg);
      }
      ofs << "# Fishery species pointer " << endl << "  " 
          << fishery_species_pointer << endl;
    }
  }     //NMD 13Jan2017


  dmatrix tolf(1,num_fisheries,1,nlint);
  tolf.initialize();
  dmatrix tplf(1,num_fisheries,1,nlint);
  tplf.initialize();
  for (i=1;i<=num_fisheries;i++)
  {
    int ii=0;
    ofs << "# fishery " << i << endl;
    for (j=1;j<=num_fish_times(i);j++)
    {
      int ir=realization_region(i,j);
      int ip=realization_period(i,j);
      int fi=realization_incident(i,j);
      int yr;
      int mnth;
      int wk;
      if (len_sample_size(ir,ip,fi)>0)
      {
        if (month_factor!=0 && month_factor!=1)
        {         
          date_struc newdate=get_old_time(year(ir,ip),
            month(ir,ip),week(ir,ip),month_factor,first_time);
          yr=newdate.year+year1-1;
          mnth=newdate.month;
          wk=newdate.week;
        }
        else
        {
          yr=year(ir,ip);
          mnth=month(ir,ip);
          wk=week(ir,ip);
        }
        ofs << yr << " " << mnth << " " << wk << endl;
        if (pmsd)  //NMD_17Mar2017
        {
          ofs << pmsd->fisc_lf(ir,ip,fi) << endl; 
        }
        else
        {
          int tmp = 1;
          ofs << tmp << endl;
        }  //NMD_17Mar2017
        ofs << len_sample_size(ir,ip,fi) << endl;    //NMD_21Mar2017
        int pi=parent(ir,ip,fi);
        if (fish_flags(pi,57) ==3 && fish_flags(pi,26) ==3)
        {
          int mn=month(ir,ip);
          int wk=week(ir,ip);
          //prob(j,i)= value(propp(j)) 
          //  * value(age_length_fishery_size_dist(mn,wk,pi,j,i));
//   old code //NMD_24jan2023
          dvector mnl=value(age_length_fishery_size_dist(mn,wk,pi))*fmid;
          ofs << (mnl-shlen)/filen + 1.0 << endl;
        }
        else
        {
          ofs << (mean_length(ir,ip,fi)-shlen)/filen + 1.0 << endl;
        }
        //ofs << mean_length(ir,ip,fi) << endl;
        MY_DOUBLE_TYPE ssum=sum(len_freq(ir,ip,fi));
        ofs << ssum << endl;
        ofs << len_freq(ir,ip,fi)/ssum << endl;
        dvar_vector& propp=prop(ir,ip,fi);
        dvar_vector& mean_len=mean_length(ir,ip,fi);
        dvar_vector& sigg=sdevs(ir,ip,fi);
        tolf(i)+=len_sample_size(ir,ip,fi)*len_freq(ir,ip,fi); 
        tplf(i)+=len_sample_size(ir,ip,fi)*
           main_length_calcs_print(ofs,propp,mean_len,sigg,ir,ip,fi,i);
        ofs << endl;
      }
    }
  }
  ofs << "# fishery " << "totals" << endl;
  for (i=1;i<=num_fisheries;i++)
  {
    ofs << i << " " << " 1" << " " << " 1" << endl;
    int ia;
    for (ia=1;ia<=nage;ia++)
      ofs << " -1.000";
    ofs << endl;
    ofs << " 1.0000" << endl;
    ofs << tolf(i) << endl;
    ofs << tplf(i) << endl;
    ofs << endl;
    for (ia=1;ia<=nage;ia++)
    {
      for (int ilen=1;ilen<=nlint;ilen++)
        ofs << " 0.0000";
      ofs << endl;
    }
  }
}

void dvar_len_fish_stock_history::print_pred_frequencies(uostream& ofs)
{
  //d3_array tprob(1,num_fish_periods, 1,num_fish_incidents, 1, nlint);
  tprob.initialize();
  //tprob+=1.e-5;
  int i;
  int j;
  dvariable temp;
  dvariable fdiff2;
  dvariable prob;
  dvariable u;
  MY_DOUBLE_TYPE a2 = 35.91908881;
  MY_DOUBLE_TYPE a3 = 785.858644;

  dvar_vector rnuml(1,nlint);
  ivector num_freqs(1,num_fisheries);
  for (i=1;i<=num_fisheries;i++)
  {
    int ii=0;
    for (j=1;j<=num_fish_times(i);j++)
    {
      int ir=realization_region(i,j);
      int ip=realization_period(i,j);
      int fi=realization_incident(i,j);
      if (len_sample_size(ir,ip,fi)>0)
      ii++;
    }
    num_freqs(i)=ii;
  }
  dvector dt(1,10);
  dt.initialize();
  dt(1)=year1;
  ofs << dt;

  ofs << nlint <<  " " << shlen << " " << filen;
  ofs << num_fisheries;
  ofs << num_freqs;
  ofs << nage;

  for (i=1;i<=num_fisheries;i++)
  {
    int ii=0;
    for (j=1;j<=num_fish_times(i);j++)
    {
      int ir=realization_region(i,j);
      int ip=realization_period(i,j);
      int fi=realization_incident(i,j);
      if (len_sample_size(ir,ip,fi)>0)
      {
        ofs << year(ir,ip) << " " << month(ir,ip) << " " << week(ir,ip);
        ofs << (mean_length(ir,ip,fi)-shlen)/filen + 1.0;
        MY_DOUBLE_TYPE ssum=sum(len_freq(ir,ip,fi));
        ofs << ssum;
        ofs << len_freq(ir,ip,fi)/ssum;
        dvar_vector& propp=prop(ir,ip,fi);
        dvar_vector& mean_len=mean_length(ir,ip,fi);
        dvar_vector& sigg=sdevs(ir,ip,fi);
        main_length_calcs_print(ofs,propp,mean_len,sigg);
      }
    }
  }
}

dvector dvar_len_fish_stock_history::main_length_calcs_print(ofstream& ofs,
  dvar_vector& propp,dvar_vector& mean_len, dvar_vector& sigg,
  int ir,int ip,int fi,int pi)
{
  const MY_DOUBLE_TYPE a2 = 35.91908881;
  const MY_DOUBLE_TYPE a3 = 785.858644;
  const MY_DOUBLE_TYPE a22 = 2*35.91908881;
  const MY_DOUBLE_TYPE a33 = 3*785.858644;
  int js=1;
  int break_flag=0;
  dvector tprob(1,nlint);
  dmatrix prob(1,nage,1,nlint);
  tprob.initialize();
  prob.initialize();
  int i=1; 
  int j;
  for (i=1; i<=nlint; i++) 
  {
    MY_DOUBLE_TYPE& fmidd=fmid(i);
   
    for (j=js; j<=nage; j++) 
    {
      MY_DOUBLE_TYPE t=(fmidd-value(mean_len(j)))/value(sigg(j));
      if ( fabs(t) <= 3.03)
      {
        if ( fabs(t) <= 3.00)
          prob(j,i)=exp(-t*t/2.e0)/value(sigg(j));
        else if (t > 0.e0)
        {
          MY_DOUBLE_TYPE u=t-3.03;
          prob(j,i)=(u*u*(a2+a3*u))/value(sigg(j));
        }
        else if (t < 0.e0)
        {
          MY_DOUBLE_TYPE u=t+3.03;
          prob(j,i)=(u*u*(a2-a3*u))/value(sigg(j));
        }
        if (fish_flags(pi,57) ==3 && fish_flags(pi,26) ==3)
        {
          int mn=month(ir,ip);
          int wk=week(ir,ip);
          prob(j,i)= value(propp(j)) 
            * value(age_length_fishery_size_dist(mn,wk,pi,j,i));
        }
        else
        {
          prob(j,i)= value(propp(j)) * prob(j,i);
        }
        tprob(i) += prob(j,i);
      }
      else
      {
        if (mean_len(j)>fmidd)
          break_flag=1;
        else
          js=j;
        if (break_flag ==1) break;
      }
      if (break_flag ==1) break;
    }
    break_flag=0;
  }
  
  MY_DOUBLE_TYPE ssum=0.0;
  ssum=sum(tprob);
  if (ssum<=0)
  {
    ssum=1.e-20;
  }

  tprob/=ssum;
  ofs << setprecision(4) << setfixed()<< tprob << endl << endl;
  for (j=1; j<=nage; j++) 
  {
    prob(j)/=ssum;
    ofs << setprecision(4) << setfixed()<< prob(j) << endl;
  }
  return tprob;
}

void dvar_len_fish_stock_history::main_length_calcs_print(uostream& ofs,
  dvar_vector& propp,dvar_vector& mean_len,
  dvar_vector& sigg)
{
  const MY_DOUBLE_TYPE a2 = 35.91908881;
  const MY_DOUBLE_TYPE a3 = 785.858644;
  const MY_DOUBLE_TYPE a22 = 2*35.91908881;
  const MY_DOUBLE_TYPE a33 = 3*785.858644;
  int js=1;
  int break_flag=0;
  dvector tprob(1,nlint);
  dmatrix prob(1,nage,1,nlint);
  tprob.initialize();
  prob.initialize();
  int i=1; 
  int j;
  for (i=1; i<=nlint; i++) 
  {
    MY_DOUBLE_TYPE& fmidd=fmid(i);
   
    for (j=js; j<=nage; j++) 
    {
      MY_DOUBLE_TYPE t=(fmidd-value(mean_len(j)))/value(sigg(j));
      if ( fabs(t) <= 3.03)
      {
        if ( fabs(t) <= 3.00)
          prob(j,i)=exp(-t*t/2.e0)/value(sigg(j));
        else if (t > 0.e0)
        {
          MY_DOUBLE_TYPE u=t-3.03;
          prob(j,i)=(u*u*(a2+a3*u))/value(sigg(j));
        }
        else if (t < 0.e0)
        {
          MY_DOUBLE_TYPE u=t+3.03;
          prob(j,i)=(u*u*(a2-a3*u))/value(sigg(j));
        }
        prob(j,i)= value(propp(j)) * prob(j,i);
        tprob(i) += prob(j,i);
      }
      else
      {
        if (mean_len(j)>fmidd)
          break_flag=1;
        else
          js=j;
        if (break_flag ==1) break;
      }
      if (break_flag ==1) break;
    }
    break_flag=0;
  }
  
  MY_DOUBLE_TYPE ssum=0.0;
  ssum=sum(tprob);
  if (ssum<=1.e-100)
  {
    cout << "this can't happen" << endl;
    ssum+=1.e-20;
  }
     
  tprob/=ssum;
  ofs << tprob;
  for (j=1; j<=nage; j++) 
  {
    prob(j)/=ssum;
    ofs << prob(j);
  }
}



void dvar_len_fish_stock_history::print_pred_wght_frequencies(ofstream& ofs)
{
//
// Print version number and pointers     //NMD_13Jan2017
  if(parest_flags(200) >= 1052)
  {
    ofs << "# FIT 3" << endl;
  }  //NMD 13Jan2017

  calculate_the_mean_weight();
  wtprob.initialize();
  int i;
  int j;
  dvariable temp;
  dvariable fdiff2;
  dvariable prob;
  dvariable u;
  MY_DOUBLE_TYPE a2 = 35.91908881;
  MY_DOUBLE_TYPE a3 = 785.858644;

  dvar_vector rnuml(1,nwint);
  ivector num_freqs(1,num_fisheries);
  for (i=1;i<=num_fisheries;i++)
  {
    int ii=0;
    for (j=1;j<=num_fish_times(i);j++)
    {
      int ir=realization_region(i,j);
      int ip=realization_period(i,j);
      int fi=realization_incident(i,j);
      if (wght_sample_size(ir,ip,fi)>0)
      ii++;
    }
    num_freqs(i)=ii;
  }
  dvector dt(1,10);
  dt.initialize();
  dt(1)=0;
  dt(3)=20;
  ofs << dt << endl;

  ofs << nwint <<  " " << wshlen << " " << wfilen << endl;
  ofs << num_fisheries+1 << endl;
  ofs << num_freqs << " " << num_fisheries << endl;
  ofs << nage << endl;

// Print pointer for species in each fishery
  if(parest_flags(200) > 1051)   //NMD 13Jan2017
  {
    ivector fishery_species_pointer(1,num_fisheries);
    if(!pmsd)
    {
      fishery_species_pointer=1;
      ofs << "# Fishery species pointer " << endl << "  " 
          << fishery_species_pointer  << endl;
    } else {
      int mmin=fishery_regions.indexmin();
      int mmax=fishery_regions.indexmax();
      for(int ir=mmin;ir<=mmax;ir++)
      {
        int ireg=fishery_regions(ir);
        fishery_species_pointer(ir)=pmsd->region_species_pointer(ireg);
      }
      ofs << "# Fishery species pointer " << endl << "  " 
          << fishery_species_pointer << endl;
    }
  }     //NMD 13Jan2017

  dmatrix tolf(1,num_fisheries,1,nwint);
  tolf.initialize();
  dmatrix tplf(1,num_fisheries,1,nwint);
  tplf.initialize();
  for (i=1;i<=num_fisheries;i++)
  {
    int ii=0;
    ofs << "# fishery " << i << endl;
    for (j=1;j<=num_fish_times(i);j++)
    {
      int ir=realization_region(i,j);
      int ip=realization_period(i,j);
      int fi=realization_incident(i,j);
      int yr;
      int mnth;
      int wk;
      if (wght_sample_size(ir,ip,fi)>0)
      {
        if (month_factor!=0 && month_factor!=1)
        {         
          date_struc newdate=get_old_time(year(ir,ip),
            month(ir,ip),week(ir,ip),month_factor,first_time);
          yr=newdate.year+year1-1;
          mnth=newdate.month;
          wk=newdate.week;
        }
        else
        {
          yr=year(ir,ip);
          mnth=month(ir,ip);
          wk=week(ir,ip);
        }
        ofs << yr << " " << mnth << " " << wk << endl;
        if (pmsd)  //NMD_17Mar2017
        {
          ofs << pmsd->fisc_wf(ir,ip,fi) << endl; 
        }
        else
        {
          int tmp = 1;
          ofs << tmp << endl;
        }  //NMD_17Mar2017
        ofs << wght_sample_size(ir,ip,fi) << endl;    //NMD_21Mar2017
        int pi=parent(ir,ip,fi);
        if (fish_flags(pi,57) ==3 && fish_flags(pi,26) ==3)
        {
          int mn=month(ir,ip);
          int wk=week(ir,ip);
          dvector mnw=value(age_weight_fishery_size_dist(mn,wk,pi))*realwmid;
          ofs << (mnw-wshlen)/wfilen + 1.0 << endl;
        }
        else
        {
          ofs << (mean_weight(ir,ip,fi)-wshlen)/wfilen + 1.0 << endl;
        }
        //ofs << mean_length(ir,ip,fi) << endl;
        MY_DOUBLE_TYPE ssum=sum(wght_freq(ir,ip,fi));
        if (ssum==0.0)
        {
          cout <<"printfrq.cpp " << wght_sample_size(ir,ip,fi) << endl;
          cout <<"printfrq.cpp " << ssum << endl;
        }
        ofs << "#wghtsum" << endl;
        ofs << "1" << endl;
        ofs << "#wghtfrq" << endl;
        ofs << wght_freq(ir,ip,fi)/ssum << endl;
        dvar_vector& propp=prop(ir,ip,fi);
        dvar_vector& mean_len=mean_length(ir,ip,fi);
        dvar_vector& sigg=sdevs(ir,ip,fi);
        tolf(i)+=wght_sample_size(ir,ip,fi)*wght_freq(ir,ip,fi); 
        tplf(i)+=wght_sample_size(ir,ip,fi)*
	  main_wght_calcs_print(ofs,propp,mean_len,sigg,ir,ip,fi,pi);
        ofs << endl;
      }
    }
  }
  ofs << "# fishery " << "totals" << endl;
  int ia;
  for (i=1;i<=num_fisheries;i++)
  {
    ofs << i << " " << " 1" << " " << " 1" << endl;
    for (ia=1;ia<=nage;ia++)
      ofs << " -1.000";
    ofs << endl;
    ofs << " 1.0000" << endl;
    ofs << tolf(i) << endl;
    ofs << tplf(i) << endl;
    ofs << endl;
    for (ia=1;ia<=nage;ia++)
    {
      for (int ilen=1;ilen<=nwint;ilen++)
        ofs << " 0.0000";
      ofs << endl;
    }
  }
}

dvector dvar_len_fish_stock_history::main_wght_calcs_print(ofstream& ofs,
  dvar_vector& propp,dvar_vector& mean_len,
  dvar_vector& sigg, int ir, int ip, int fi, int pi)
{
  const MY_DOUBLE_TYPE a2 = 35.91908881;
  const MY_DOUBLE_TYPE a3 = 785.858644;
  const MY_DOUBLE_TYPE a22 = 2*35.91908881;
  const MY_DOUBLE_TYPE a33 = 3*785.858644;
  int js=1;
  int break_flag=0;
  dvector tprob(1,nwint);
  dmatrix prob(1,nage,1,nwint);
  tprob.initialize();
  prob.initialize();
  int i=1; 
  int j;
  for (i=1; i<=nwint; i++) 
  {
    MY_DOUBLE_TYPE& fmidd=wmid(i);
   
    for (j=js; j<=nage; j++) 
    {
      MY_DOUBLE_TYPE t=(fmidd-value(mean_len(j)))/value(sigg(j));
      if ( fabs(t) <= 3.03)
      {
        if ( fabs(t) <= 3.00)
          prob(j,i)=exp(-t*t/2.e0)/value(sigg(j));
        else if (t > 0.e0)
        {
          MY_DOUBLE_TYPE u=t-3.03;
          prob(j,i)=(u*u*(a2+a3*u))/value(sigg(j));
        }
        else if (t < 0.e0)
        {
          MY_DOUBLE_TYPE u=t+3.03;
          prob(j,i)=(u*u*(a2-a3*u))/value(sigg(j));
        }
        prob(j,i)= value(propp(j)) * prob(j,i);
        tprob(i) += prob(j,i);
      }
      else
      {
        if (mean_len(j)>fmidd)
          break_flag=1;
        else
          js=j;
        if (break_flag ==1) break;
      }
      if (break_flag ==1) break;
    }
    break_flag=0;
  }
//   new code //NMD_24jan2023
  dvar_vector tprobb(1,nwint);
  tprobb.initialize();
  if (fish_flags(pi,57) ==3 && fish_flags(pi,26) ==3)
  {
     int mn=month(ir,ip);
     int wk=week(ir,ip);
     dvar_vector& propp=prop(ir,ip,fi);
     dvar_matrix& aws=age_weight_fishery_size_dist(mn,wk,pi);
     main_weight_calcs_len_based(tprobb,propp,aws);
     if (pmsd==0 || pmsd->current_species==1)
       tprobb=elem_div(tprobb,wm2);
     else
       tprobb=elem_div(tprobb,pmsd->wm2(pmsd->current_species));     
     tprobb/=(1.e-20+sum(tprobb));
     //NMD_17mar2023
    prob.initialize();
    for (j=1;j<=nage;j++)
    {
      prob(j)=value(propp(j)*aws(j));
    }
  }

  MY_DOUBLE_TYPE ssum=0.0;
  tprob=elem_div(tprob,wm2);
  ssum=sum(tprob);
  tprob/=ssum;
  ofs << "#predicted weight distribution" << endl;
  if (fish_flags(pi,57) ==3 && fish_flags(pi,26) ==3)
  {
    ofs << setprecision(4) << setfixed()<< tprobb << endl << endl;
  }
  else
  {
    ofs << setprecision(4) << setfixed()<< tprob << endl << endl;
  }
  ofs << "#predicted weight distribution by age" << endl;
  for (j=1; j<=nage; j++) 
  {
    prob(j)=elem_div(prob(j),wm2);
    prob(j)/=ssum;
    ofs << setprecision(4) << setfixed()<< prob(j) << endl;
  }
  //NMD_17Mar2023
  if (fish_flags(pi,57) ==3 && fish_flags(pi,26) ==3)
  {
    tprob=value(tprobb);
  }
  return tprob;
}

