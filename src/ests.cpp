/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


#define HOME_VERSION
#include "all.hpp"

#include "variable.hpp"

dvector cbiocalc(dvar_len_fish_stock_history& fsh);
dvector cbiocalc2(dvar_len_fish_stock_history& fsh);
dvector cbiocalc_wt(dvar_len_fish_stock_history& fsh);
extern dvector mean_weights_kludge;

dvariable calculate_overall_exploitation();
void get_availability_props(dvar_vector& avail_prop,int nfa,
  dvar_vector& avail_coff, dvar_vector& unavail_prop);

MY_DOUBLE_TYPE max(_CONST dmatrix&);

dvector get_generic_mean_lengths(dvar_len_fish_stock_history& fsh)
{	
  int xip;
  int xir=1;
  int lb=1;
  int ub=fsh.num_regions;
  int cs=1;
  if (fsh.pmsd)
  {
    cs=fsh.pmsd->current_species;
    lb=fsh.pmsd->region_bounds(cs,1);
    ub=fsh.pmsd->region_bounds(cs,2);
  }
  
  int break_flag=0;
  for (xir=lb;xir<=ub;xir++)
  {	     
    for (xip=1;xip<=fsh.num_fish_periods(xir);xip++)
    {
      if (fsh.num_fish_incidents(xir,xip)>0) 
      {
        break_flag=1;
        break;
      }
    }	    
    if (break_flag==1)break;
  }  
  return value(fsh.mean_length(xir,xip,1));
}  

dvar_vector vget_generic_mean_lengths(dvar_len_fish_stock_history& fsh)
{	
  int xip;
  int xir=1;
  int lb=1;
  int ub=fsh.num_regions;
  int cs=1;
  if (fsh.pmsd)
  {
    cs=fsh.pmsd->current_species;
    lb=fsh.pmsd->region_bounds(cs,1);
    ub=fsh.pmsd->region_bounds(cs,2);
  }
  
  int break_flag=0;
  for (xir=lb;xir<=ub;xir++)
  {	     
    for (xip=1;xip<=fsh.num_fish_periods(xir);xip++)
    {
      if (fsh.num_fish_incidents(xir,xip)>0) 
      {
        break_flag=1;
        break;
      }
    }	    
    if (break_flag==1)break;
  }  
  return fsh.mean_length(xir,xip,1);
}  

void ests_write(ofstream& of,dvar_len_fish_stock_history& fsh)
{
  dvector ml;
  //cout << "entering ests_write" << endl;
  int ccond=0;
  if (fsh.age_flags(92)==2) ccond=1;
  int i;
  {
    int xip,xir;
    int break_flag=0;
    for (xir=1;xir<=fsh.num_regions;xir++)
    {	     
      for (xip=1;xip<=fsh.num_fish_periods(xir);xip++)
      {
        if (fsh.num_fish_incidents(xir,xip)>0) 
        {
          break_flag=1;
          break;
        }
      }	    
      if (break_flag==1)break;
    }  

    ml=get_generic_mean_lengths(fsh);
    dvector sig=sqrt(value(fsh.global_vars));
    dvector dif(1,fsh.nage-1);

    if (sum(column(fsh.fish_flags,31)))
    {
      fsh.pmsd_error();
      dvar_matrix csr=cutoff_sel_report(fsh,fsh.vb_coff,fsh.var_coff);
      of << "Selectivity at length" << endl;
      of << setw(8) << setprecision(3) << csr << endl;
    }

  of << "exp(fsh.N(1,fsh.nyears,1) = " << exp(fsh.N(1,fsh.nyears,1)) << endl;

    for (int j=1;j<fsh.nage;j++)
    {
      dif(j)=(ml(j+1)-ml(j))/(sig(j)+sig(j+1));
    }
    of << "the mean lengths" << endl;
    of << ml << endl;
    of << --ml(2,fsh.nage)-ml(1,fsh.nage-1) << endl;
    of << "Scaled rao distance of the mean lengths" << endl;
    of << dif << endl;
  }
   of << "Natural mortality in first year" << endl
	   << exp(value(fsh.nat_mort(1))) << endl;

    for (i=1;i<=fsh.num_fisheries;i++)
    {
      of << endl << "Fishery " << i << endl;
      for (int nt =1;nt<=fsh.num_fish_times(i);nt++)
      {
        of << "Fishery Realization" << nt << endl;
        int rr=fsh.realization_region(i,nt);
        int rp=fsh.realization_period(i,nt);
        int ri=fsh.realization_incident(i,nt);
        of << setfixed()<< setprecision(4) << setw(8) << "  " 
           << exp(value(fsh.fish_mort(rr,rp,ri)))/
              max(exp(value(fsh.fish_mort(rr,rp,ri)))) << endl;
      }
    }
    int nt;

    if (!fsh.age_flags(92))
    {
      for (i=1;i<=fsh.num_fisheries;i++)
      {
        of << endl << "Fishery " << i << endl;
        for ( nt =1;nt<=fsh.num_fish_times(i);nt++)
        {
          int rr=fsh.realization_region(i,nt);
          int rp=fsh.realization_period(i,nt);
          int ri=fsh.realization_incident(i,nt);
          {
            of << setscientific() << setprecision(4) << " " 
               << exp(value(fsh.catchability(rr,rp,ri)));
          }
          if (!(nt%10)) of << endl;
        }
        of << endl << endl << endl;
        for (nt =1;nt<=fsh.num_fish_times(i);nt++)
        {
          int rr=fsh.realization_region(i,nt);
          int rp=fsh.realization_period(i,nt);
          int ri=fsh.realization_incident(i,nt);
          {
            of << setscientific() << setprecision(4) << " " 
               << exp(value(fsh.catchability(rr,rp,ri)))*
                  exp(value(fsh.effort_devs(rr,rp,ri)));
          }
          if (!(nt%10)) of << endl;
        }
        of << endl << endl;
      }
    }
    else
    {
      of << endl << "Raw implicit catchabilities" << i << endl;
      for (i=1;i<=fsh.num_fisheries;i++)
      {
        of << endl << "Fishery " << i << endl;
        for (int nt =1;nt<=fsh.num_fish_times(i);nt++)
        {
          {
            of << setscientific() << setprecision(4) << " " 
               << exp(value(fsh.implicit_catchability(i,nt)));
          }
          if (!(nt%10)) of << endl;
        }
        of << endl << endl << endl;
      }
    }
    of << endl;

    if (!ccond)
    {
      dvector ac(1,fsh.num_fisheries);
      dvector v(1,fsh.num_fisheries);
      ac.initialize();
      v.initialize();
      for (i=1;i<=fsh.num_fisheries;i++)
      {
        of << endl << "Fishery " << i << endl;
        for (int nt =1;nt<fsh.num_fish_times(i);nt++)
        {
          int rr=fsh.realization_region(i,nt);
          int rp=fsh.realization_period(i,nt);
          int ri=fsh.realization_incident(i,nt);
          int rr1=fsh.realization_region(i,nt+1);
          int rp1=fsh.realization_period(i,nt+1);
          int ri1=fsh.realization_incident(i,nt+1);
       
          ac(i)+= exp(value(fsh.effort_devs(rr,rp,ri)))
            *exp(value(fsh.effort_devs(rr1,rp1,ri1)));
          v(i)+= exp(value(fsh.effort_devs(rr,rp,ri)))
            *exp(value(fsh.effort_devs(rr,rp,ri)));
        }
      }
      of << " The effort deviation lag 1 auto correlations " << endl;
      of << setfixed() << setprecision(3) << setw(6) 
         << elem_div(ac,v) << endl;
    } 


    for (i=1;i<=fsh.num_fisheries;i++)
    {
      of << endl << "Fishery " << i << endl;
      for (int nt =1;nt<=fsh.num_fish_times(i);nt++)
      {
        of << "Fishery Realization" << nt << endl;
        int rr=fsh.realization_region(i,nt);
        int rp=fsh.realization_period(i,nt);
        int ri=fsh.realization_incident(i,nt);
        of << setfixed()<< setprecision(4) << setw(8) << "  " 
           << exp(value(fsh.fish_mort(rr,rp,ri))) << endl;
      }
    }

    // #ifndef __NDPX__
    //lhmodel_fit(fsh,of);
    int ny=1;
    int ir=1;

    of << endl << " Predicted numbers of fish" << endl;
    for ( ir=1;ir<=fsh.num_regions;ir++)
    {
      of << "Group " << ir << endl;
      of <<  setscientific() << setprecision(3) << exp(value(fsh.N(ir))) << endl;
      
      of << "Last Year " << ir << endl;
      of <<  exp(value(fsh.N(ir,fsh.nyears))) << endl;
    }

      //fsh.pmsd_error();  //NMD11Apr2012
      int mmin=fsh.N(1).rowmin();
      int mmax=fsh.N(1).rowmax();
      dmatrix rbio(1,fsh.num_regions,mmin,mmax);
      dvector vvar=value(fsh.global_vars);
      MY_DOUBLE_TYPE lwc=fsh.len_wt_coff;
      if (lwc==0) lwc=1.0;
      for (ir=1;ir<=fsh.num_regions;ir++)
      {
        for (int iy=mmin;iy<=mmax;iy++)
        {
          rbio(ir,iy)=lwc*sum(elem_prod(exp(value(fsh.N(ir,iy))),
            pow(ml,3)+
            3.*elem_prod(ml,vvar)));
        }
      }
      of << "Absolute biomass by area" << endl;
      of << setprecision(4)  << trans(rbio) << endl; 

      rbio=rbio/max(rbio);
      
      of << "Relative biomass by area" << endl;
      of << setprecision(4)  << trans(rbio) << endl; 


    dvector csum(1,fsh.num_fisheries);
    csum.initialize();
    dvector nsum(1,fsh.num_fisheries);
    nsum.initialize();
    if (fsh.age_flags(121))
    {
      of << "Survival analysis results" << endl;

      for (int it=1;it<=fsh.num_tag_releases;it++)
      {
        int fishery=fsh.tag_recaptures_by_length(it,1,2);
        int irr=fsh.fishery_regions(fishery);
        int ipp=fsh.terminal_tag_period(it,irr);
         of <<fsh. tag_return_probability(it) << setw(12)
              << fishery << setw(12)
              << ipp << setw(12) << endl;
      }
    }
   
    for (i=1;i<=fsh.num_fisheries;i++)
    {
      of << endl << "Fishery " << i <<  endl << 
          "  fishery  realiz. realiz. ratio      predicted     observed   observed    effort      number of" <<  endl <<
          "  realiz.  region  period             total cat     tot cat    effort      devs       fish in pop" << endl;

      for (int nt =1;nt<=fsh.num_fish_times(i);nt++)
      {
        of << "  " << setw(4) << nt;
        int rr=fsh.realization_region(i,nt);
        int rp=fsh.realization_period(i,nt);
        int ri=fsh.realization_incident(i,nt);
       
        MY_DOUBLE_TYPE ratio=value(fsh.pred_totcatch(rr,rp,ri))/ 
                (.1+fsh.obs_tot_catch(rr,rp,ri)); 
        if (fsh.parest_flags(41)==0)
        {
          if ( (ratio > 2.0 || ratio < 0.5L) && 
	       (value(fsh.pred_totcatch(rr,rp,ri))>1.1) )
            of << " ***";
          else
            of << "    ";
          of << setw(4) << rr << " " << setw(4) << rp << "  " 
             << setfixed()<< setprecision(2) << setw(12)  
             << ratio << "  " 
             << setfixed()<< setprecision(1)  << setw(12)    
             << value(fsh.pred_totcatch(rr,rp,ri)) <<  "  "
             << setfixed()<< setprecision(1)   << setw(12)  
             << fsh.obs_tot_catch(rr,rp,ri) << setw(12)
             << fsh.effort(rr,rp,ri) << setw(12) 
             << fsh.effort_devs(rr,rp,ri) << setw(12)
             << sum(exp(fsh.num_fish(rr,rp)))/1000.0;
          //if (ratio > 2.0 || ratio < 0.5L) 
          //   fsh.print_tag_data(rr,rp,of);
          of << endl; 
        }
        else
        {
          of << setfixed()<< setprecision(4) << setw(8) << "  " 
             << exp(value(fsh.catch(rp,ri)))*mean_weights_kludge <<  "  "
             << fsh.obs_tot_catch(rp,ri) << endl;
        }
      }
    }
  
 {
   of << "The recruitment" << endl;
   int af57=fsh.age_flags(57);
   dvector pr;
   dvector tmp(1,fsh.nyears);
   tmp.initialize();
   ivector numpr;
   if (af57)
   {
     pr.allocate(1,af57);
     numpr.allocate(1,af57);
     pr.initialize();
     numpr.initialize();
   }
   for (i=1;i<=fsh.nyears;i++)
   {
     for (int ii=1;ii<=fsh.num_regions;ii++)
     {
       tmp(i)+=exp(value(fsh.N(ii,i,1)));
     }
     if (af57)
     {
       int isub=((i-1)%af57)+1;
       pr(isub)+=tmp(i);
       numpr(isub)+=1;
     }

     if (!af57)
     {
       of << setw(5) << i << " " << tmp(i) << " " << endl;
     }
     else
     {
       of << setw(5) << setprecision(1) << setfixed()
          << i/double(af57) << " " << tmp(i) << " " << endl;
     }
   }
   of << endl;
   if (af57)
   {
     of << "Average recruitment by period " << endl;
     of << elem_div(pr,numpr) << endl;
     for (int ii=1;ii<=fsh.num_regions;ii++)
     of << "Recruitment by period" << endl;
     for (int iper=1;iper<=af57;iper++)
     {
       for (int i=iper;i<=fsh.nyears;i+=af57)
       {
         of << setw(5) << setprecision(1) << setfixed()
            << i/double(af57) << " " << tmp(i) << " " << endl;
       }
       of << endl;
     }
   }
 }
 {
   dvector ssum(1,fsh.nage);
   dvector szsum(1,fsh.nage);
   ssum.initialize();
   szsum.initialize();
  
   for (ir=1;ir<=fsh.num_regions;ir++)
   {
     for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++) 
     {
       for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) 
       {                                 
         szsum+=exp(value(fsh.num_fish(ir,ip)));
         ssum+=elem_prod(exp(value(fsh.fish_mort(ir,ip,fi))),
           exp(value(fsh.num_fish(ir,ip))));
       }
     }
   }
   ssum=elem_div(ssum,szsum);
   MY_DOUBLE_TYPE mult=double(fsh.num_fish_data_recs)/double(fsh.nyears);
   of << "# Average fish mort per year by age class is " << endl 
      << "#" << setprecision(6) << ssum*mult << endl;
 }
 {
   dmatrix ssum(1,fsh.nyears,1,fsh.nage);
   dmatrix szsum(1,fsh.nyears,1,fsh.nage);
   ssum.initialize();
   szsum.initialize();
  
   for (ir=1;ir<=fsh.num_regions;ir++)
   {
     for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++) 
     {
       for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) 
       {                                 
         szsum(fsh.year(ir,ip))+=exp(value(fsh.num_fish(ir,ip)));
         ssum(fsh.year(ir,ip))+=elem_prod(exp(value(fsh.fish_mort(ir,ip,fi))),
           exp(value(fsh.num_fish(ir,ip))));
       }
     }
   }
   ssum=elem_div(ssum,1.e-10+szsum);
   of << "# Average Fish mort by age class is " << endl 
      << "#" << endl << setprecision(6) << colsum(ssum)/fsh.nyears << endl;
   of << "# Fish mort by year by age class is " << endl 
      << "#" << endl << setprecision(6) << ssum << endl;
  }
  if (fsh.nage>5)
  {
    dmatrix exploitation_by_year(1,fsh.nyears,1,fsh.nage);
    dvector age_summed_exploitation_by_year(1,fsh.nyears);
    dmatrix totnum(1,fsh.nyears,1,fsh.nage);
    d3_array regional_totnum(1,fsh.num_regions,1,fsh.nyears,1,fsh.nage);
    d3_array exploitation_by_region_and_year(1,fsh.num_regions,
      1,fsh.nyears,1,fsh.nage);
    dmatrix ssum(1,fsh.nyears,1,fsh.nage);
    d3_array regional_ssum(1,fsh.num_regions,1,fsh.nyears,1,fsh.nage);
    dmatrix totcatch(1,fsh.nyears,1,fsh.nage);
    d3_array regional_totcatch(1,fsh.num_regions,1,fsh.nyears,1,fsh.nage);
    regional_totnum.initialize();
    regional_totcatch.initialize();
    totnum.initialize();
    totcatch.initialize();
  
    for (ir=1;ir<=fsh.num_regions;ir++)
    {
      for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++) 
      {
        for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) 
        {                                 
          // totcatch is the total catch by year and age
          totcatch(fsh.year(ir,ip))+=exp(value(fsh.catch(ir,ip,fi)));
          // regional totcatch is the total catch by region year and age
          regional_totcatch(ir,fsh.year(ir,ip))+=
            exp(value(fsh.catch(ir,ip,fi)));
        }
      }
    }
    for (i=1;i<=fsh.nyears;i++)
    {
      for (ir=1;ir<=fsh.num_regions;ir++)
      {
        // totnum is the total number of fish at the beginning of year i
        totnum(i)+=exp(value(fsh.N(ir,i)));
        // totnum is the total number of fish in each region 
        // at the beginning of year i
        regional_totnum(ir,i)=exp(value(fsh.N(ir,i)));
      }
    }
    of << "Initial population by age" << endl << totnum(1) << endl;
    if (fsh.nyears>4)
    {
      for (i=1;i<=fsh.nyears;i++)
      {
        exploitation_by_year(i)=elem_div(totcatch(i),1.e-10+totnum(i));
        age_summed_exploitation_by_year(i)=
           sum(totcatch(i)(5,fsh.nage))/(1.e-10+sum(totnum(i)(5,fsh.nage)));
      }
    }
  
   //if (fsh.direction_flag==1)
    {
      d3_array catch_by_fishery_by_year= fsh.report_catch_by_year();
      of << "Predicted catch by fishery by year and age"  << endl;
      int ny=catch_by_fishery_by_year(1).indexmax();
      for (i=1;i<=fsh.num_fisheries;i++)
      {
        of << "Fishery " << i << endl;
        for (int j=1;j<=ny;j++)
        {
         //of << setw(8) << ivector(catch_by_fishery_by_year(i,j)+0.5L) << endl;
          of << setw(8) << setfixed() << setprecision(0)
             << catch_by_fishery_by_year(i,j)+0.5 << endl;
        }
      }
      if (allocated(fsh.age_freq))
      {
        d3_array estimated_catch_by_fishery_by_year=
          fsh.report_estimated_catch_by_year();
        of << "Predicted estimated_catch by fishery by year and age"  << endl;
        ny=estimated_catch_by_fishery_by_year(1).indexmax();
        for (i=1;i<=fsh.num_fisheries;i++)
        {
          of << "Fishery " << i << endl;
          for (int j=1;j<=ny;j++)
          {
            //of << setw(8) << ivector(100.*
            //   estimated_catch_by_fishery_by_year(i,j)+0.5L) << endl;
            of << setw(8) << setfixed() << setprecision(0)
               <<  (100.*estimated_catch_by_fishery_by_year(i,j)+0.5L) << endl;
          }
        }
      }
    }
  
    MY_DOUBLE_TYPE tc=0;
    MY_DOUBLE_TYPE tn=0;
    if (fsh.nyears>4)
    {
      for (i=fsh.nyears-4;i<=fsh.nyears;i++)
      {
        tc+=sum(value(totcatch)(i)(5,fsh.nage));
        tn+=sum(value(totnum)(i)(5,fsh.nage));
      }
    }
    of << "Other exploitation rate = " << tc/tn << endl; 

    for (i=1;i<=fsh.nyears;i++)
    {
      for (ir=1;ir<=fsh.num_regions;ir++)
      {
        exploitation_by_region_and_year(ir,i)=
          elem_div(regional_totcatch(ir,i),1.e-10+regional_totnum(ir,i));
      }
    }

    for (ir=1;ir<=fsh.num_regions;ir++)
    {
      dvector avg_age= colsum(regional_totnum(ir));
      avg_age/=sum(avg_age); 
      of << "# Average age structure in region " << ir << endl 
         << "#" << endl << setfixed()<< setprecision(4) 
         << avg_age << endl;
    }
    of << endl;


    for (ir=1;ir<=fsh.num_regions;ir++)
    {
      of << "# Annual exploitation rate by age class in"  
            " region " << ir << endl 
         << "#" << endl << setprecision(4)  << setw(6)
         << exploitation_by_region_and_year(ir) << endl;
    }
    of << endl;

    of << "# Overall Annual exploitation rate by age class"  
       << endl 
       << "#" << endl << setprecision(4)  << setw(6)
       << exploitation_by_year << endl;
    of << endl;

    of << "# Age-summed (5+) overall Annual exploitation rate"  
       << endl 
       << "#" << endl << setprecision(4)  << setw(6)
       << age_summed_exploitation_by_year << endl;
    of << endl;
  }

  if (!ccond)
  {
    for (i=1;i<=fsh.num_fisheries;i++)
    {
      of << "# Catchability for fishery " << i << endl;
      for (int j=1;j<=fsh.num_fish_times(i);j++)
      {
        int ir=fsh.realization_region(i,j);
        int ip=fsh.realization_period(i,j);
        int fi=fsh.realization_incident(i,j);
        int year=fsh.year(ir,ip);
        int month=fsh.month(ir,ip);
        int week=fsh.week(ir,ip);
        MY_DOUBLE_TYPE q=exp(value(fsh.catchability(ir,ip,fi)));
        MY_DOUBLE_TYPE ed=exp(value(fsh.effort_devs(ir,ip,fi)));
        of << setw(4) << year << " " << setw(4) << month << " " << setw(4) 
           << week << " " 
           << setscientific() << setprecision(4) << setw(10) << q << " " 
           << setscientific() << setprecision(4) << setw(10) << ed << endl;
      }
    }
  }
  
  if (fsh.age_flags(61))
  {
    of << endl; 
    for (i=1;i<=fsh.num_fisheries;i++)
    {
      of << "Relative Big Fish biomass Fishery " << i << endl; 
      for (int j=1;j<=fsh.num_fish_times(i);j++)
      {
        int rr=fsh.realization_region(i,j);
        int rp=fsh.realization_period(i,j);
        int ri=fsh.realization_incident(i,j);
                
        if (fsh.biomass_index(i,j)>0)
        {
          MY_DOUBLE_TYPE wt=fsh.age_flags(61);
          const prevariable& x1=fsh.RB(rr,rp,ri);
          MY_DOUBLE_TYPE& x2=fsh.biomass_index(i,j);
          of << j << " " << x1 << " " << x2 << endl; 
        }
      }
      of << endl; 
    }
  }
  of << "# The relative distribution of recruitment by region is"
     << endl << setfixed()<< setprecision(3) << setw(5) << fsh.epop_delta 
     << endl;
  if (fsh.num_regions>1)
  {
    if (fsh.age_flags(89))
    {
      of << "Diffusion by age :" << endl;
      for (int i=1;i<=fsh.nage;i++)
      {
        of << "Age class " << i << endl << setfixed()<< setprecision(3)
             << setw(7) << inv(fsh.Dad(1,i)) << endl;
      }
    }
  
    
    if (fsh.age_flags(62))
    {
      if (fsh.num_regions>1) fsh.do_the_diffusion(1,fsh.sv,fsh.N,&of);
    }
  } 
  //int im=fsh.num_tag_releases;
  /*    //NMD_5jun2025
  if (fsh.parest_flags(189))
  {
    ofstream ofs("length.fit");
    fsh.print_pred_frequencies(ofs);

    if (fsh.nwint && !fsh.parest_flags(181))
    {
      ofstream ofs1("weight.fit");
      fsh.print_pred_wght_frequencies(ofs1);
    }
  }
  */    //NMD_5jun2025

  ofstream oftag("tag.rep");
  if (fsh.num_tag_releases) fsh.print_tagging_fit_info(oftag);

  
  if (fsh.age_flags(96))
  {
    oftag << "pooled_tagnum_fish" << endl;
    oftag << setw(10) << setprecision(1) << fsh.pooled_tagnum_fish << endl;
  }
  if (fsh.num_tag_releases) fsh.print_tag_return_by_time_at_liberty(oftag);

//  if (fsh.num_regions > 1) fsh.print_movement_report(oftag);
  if (!fsh.pmsd && fsh.num_regions > 1) //NMD_jan28-19 
  {
    fsh.print_movement_report(oftag);
  }
  else if (fsh.pmsd)  //in case of multi-species NMD_jan28-19
  {
    if (fsh.pmsd->num_real_regions > 1)  //in case of single region NMD_jan28-19
    {
      fsh.print_movement_report(oftag);
    }
  }  

}
