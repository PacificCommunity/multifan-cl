/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

//int __nage;

dvariable dvar_fish_stock_history::tag_catch_equations_calc_pooled(dvar_vector& sv)
{
  //__nage=nage;
  dvariable ffpen=0.0;

  int Bigbreak=0;
  epooled_tagnum_fish_recr.initialize();
  pooledtagN=-20;
  //for (int it=1;it<=15;it++)
  for (int it=1;it<=num_tag_releases;it++)
  {
    get_initial_tag_population(sv,it);

    int current_year=tag_year(it);
    //yreport(tagN,it,current_year,num_regions);
    ivector rip(1,num_regions);
    rip=initial_tag_period(it);
    do
    {
      for (int ir=1;ir<=num_regions;ir++)
      {
        int& ip=rip(ir);
        if (year(ir,ip)==current_year &&
            ip<=num_fish_periods(ir) &&
            ip<=terminal_tag_period(it,ir))
        {
          if (ip<num_fish_periods(ir))
            ffpen+=have_data_this_year(it,ir,ip,current_year);
          else  // check to see if we need newton raphson
          {
            if (tag_flags(it,1) && 
              ip < initial_tag_period(it,ir)+tag_flags(it,1))
            {
              do_newton_raphson_for_tags(it,ir,ip,ffpen);
            }
          }
        }
        else   // there were no fisheries in this region for this year
        {
          have_no_data_this_year(it,ir,current_year);
        }
      } // loop over regions
      if (Bigbreak) break;
      if(current_year<terminal_tag_year(it) )
      {
      // Changed af(57) to af(53) J.H. 27 Nov 01
        if (age_flags(53))
        {
          if ( !((current_year+1)%age_flags(53)) )
          {
            if (num_regions>1) do_the_diffusion(current_year+1,sv,tagN(it));
          }
        }
        else
        {
          if (num_regions>1) do_the_diffusion(current_year+1,sv,tagN(it));
        }
      }
  
      current_year++;
      if (Bigbreak)
      {
        Bigbreak=0;
        break;
      }
    }
    while (current_year<=terminal_tag_year(it));
    
    calculate_tag_catches(it);
  
  }

# if defined(COUNT_FISH)
  print_tag_accounting_info();
#endif

  return ffpen;
}

dvariable dvar_fish_stock_history::have_data_this_year(int it,int ir,int& ip,
  int current_year)
{
  tagnum_fish(it,ir,ip);
  tagN(it,ir,current_year);
  tagnum_fish(it,ir,ip)=tagN(it,ir,current_year);
  dvariable ffpen=0.0;
          
  do
  {
    if (ip>=num_fish_periods(ir)) break;
    if (ip>terminal_tag_period(it,ir)) break;
    if (year(ir,ip+1)==current_year)
    { 
      if (ip==terminal_tag_period(it,ir))
      {
        cerr << "Grouped tags within a year " << endl;
      }
      if (!tag_flags(it,1) || 
        ip >= initial_tag_period(it,ir)+tag_flags(it,1))
      {
        tagnum_fish(it,ir,ip+1)=tagnum_fish(it,ir,ip)-tot_mort(ir,ip);
      }
      else
      {
        do_newton_raphson_for_tags(it,ir,ip,ffpen);
        tagnum_fish(it,ir,ip+1)=tagnum_fish(it,ir,ip)-nrtm(it,ir,ip);
      }
      tagnum_fish(it,ir,ip+1)=tagnum_fish(it,ir,ip)-tot_mort(ir,ip);
      ip++;
    }
    else
    {
      if (ip== terminal_tag_period(it,ir))
      {
        if (!tag_flags(it,1) || 
        ip >= initial_tag_period(it,ir)+tag_flags(it,1))
        {
          if (nage>2)
            --epooled_tagnum_fish_recr(ir,ip+1)(2,nage-1)
              +=mfexp(tagnum_fish(it,ir,ip)(1,nage-2)
               -tot_mort(ir,ip)(1,nage-2));
       
          epooled_tagnum_fish_recr(ir,ip+1,nage) 
            += mfexp(tagnum_fish(it,ir,ip,nage)-tot_mort(ir,ip,nage))
            + mfexp(tagnum_fish(it,ir,ip,nage-1)-tot_mort(ir,ip,nage-1));
        }
        else
        {
          do_newton_raphson_for_tags(it,ir,ip,ffpen);

          if (nage>2)
            --epooled_tagnum_fish_recr(ir,ip+1)(2,nage-1)
             +=mfexp(tagnum_fish(it,ir,ip)(1,nage-2)
              -nrtm(it,ir,ip)(1,nage-2));
       
          epooled_tagnum_fish_recr(ir,ip+1,nage) 
           += mfexp(tagnum_fish(it,ir,ip,nage)-nrtm(it,ir,ip,nage))
           + mfexp(tagnum_fish(it,ir,ip,nage-1)-nrtm(it,ir,ip,nage-1));
        }
               
        if (current_year<nyears)
        {
          pooledtagN(ir,current_year+1)=log(1.e-20+epooled_tagnum_fish_recr(ir,ip+1));
        }

        break;
      }

      tagnum_fish(it,ir,ip+1,1)=-15.0;

      if (!tag_flags(it,1) || 
      ip >= initial_tag_period(it,ir)+tag_flags(it,1))
      {
        if (nage>2)
          --tagnum_fish(it,ir,ip+1)(2,nage-1)=
            tagnum_fish(it,ir,ip)(1,nage-2)-tot_mort(ir,ip)(1,nage-2);

        tagnum_fish(it,ir,ip+1,nage)=
          log(1.e-12 + mfexp(tagnum_fish(it,ir,ip,nage-1)
            -tot_mort(ir,ip,nage-1))
            + mfexp(tagnum_fish(it,ir,ip,nage)-tot_mort(ir,ip,nage)) );
      }
      else
      {
        do_newton_raphson_for_tags(it,ir,ip,ffpen);

        if (nage>2)
          --tagnum_fish(it,ir,ip+1)(2,nage-1)=
          tagnum_fish(it,ir,ip)(1,nage-2)-nrtm(it,ir,ip)(1,nage-2);

        tagnum_fish(it,ir,ip+1,nage)=
          log(1.e-10 + mfexp(tagnum_fish(it,ir,ip,nage-1)
            -nrtm(it,ir,ip,nage-1))
            + mfexp(tagnum_fish(it,ir,ip,nage)-nrtm(it,ir,ip,nage)) );
      }
              
      if (current_year<nyears)
      {
        tagN(it,ir,current_year+1)=tagnum_fish(it,ir,ip+1);
      }
      ip++;
      break;
    }
  }
  while (1);
  return ffpen;
}

void dvar_fish_stock_history::have_no_data_this_year(int it,int ir,
  int current_year)
{
  int ng=nage;
  if (pmsd) 
   ng=get_nage_region(ir);
  if (current_year<terminal_tag_year(it))
  {
    if (!pmsd)
    {
      for (int j=1;j<ng-1;j++)      // Loop over age classes
      {
        tagN(it,ir,current_year+1,j+1)=tagN(it,ir,current_year,j)-
          exp(nat_mort(current_year,j));
      }
    }
    else
    {
      const dvar_matrix& nm = get_nat_mort_region(ir);
      for (int j=1;j<ng-1;j++)      // Loop over age classes
      {
        tagN(it,ir,current_year+1,j+1)=tagN(it,ir,current_year,j)-
          exp(nm(current_year,j));
      }
    }
    tagN(it,ir,current_year+1,ng)=
      log (
        mfexp(tagN(it,ir,current_year,ng-1)
          - exp(nat_mort(current_year,ng-1)))
       +
        mfexp(tagN(it,ir,current_year,ng)
          - exp(nat_mort(current_year,ng)))
          );
  }
}

void dvar_fish_stock_history::calculate_tag_catches(int it)
{  
  const MY_DOUBLE_TYPE one_plus=1.e0+1.e-12;
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=initial_tag_period(it,ir);ip<=terminal_tag_period(it,ir);ip++)
    {
      if (!tag_flags(it,1) || 
        ip >= initial_tag_period(it,ir)+tag_flags(it,1))
      {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {        
          tagcatch(it,ir,ip,fi)=
            exp(fish_mort(ir,ip,fi)-log(1.e-10+tot_mort(ir,ip))+
            log(one_plus-survival(ir,ip))+tagnum_fish(it,ir,ip));
        }
      }
      else
      {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {        
          tagcatch(it,ir,ip,fi)=
            exp(nrfm(it,ir,ip,fi)-log(1.e-10+nrtm(it,ir,ip))+
            log(one_plus-nrsurv(it,ir,ip))+tagnum_fish(it,ir,ip));
        }
      }
    }
  }
}

void dvar_fish_stock_history::print_tag_accounting_info(void)
{
  dvar_matrix tt(1,num_tag_releases,
    index_type(initial_tag_year),index_type(terminal_tag_year));
  tt.initialize();
     
  int it;
  for (it=1;it<=num_tag_releases;it++) 
  {
   for (int iy=initial_tag_year(it);iy<=terminal_tag_year(it);iy++)
    {
      for (int ir=1;ir<=num_regions;ir++)
        for (int j=1;j<=nage;j++)
         tt(it,iy)+=exp(tagN(it,ir,iy,j));
    }
  }

  ofstream ofs("tagrep");
  //for (it=1;it<=10;it++) 
  for (it=1;it<=num_tag_releases;it++) 
  {
    ofs << "tag group " << it << endl;
    for (int iy=initial_tag_year(it);iy<=terminal_tag_year(it);iy++)
    {
      ofs << "year " << iy << "  ";
      ofs << tt(it,iy) << endl;
    }
  }
  //pooledtagN.allocate(1,num_regions,min_tag_year,nyears,1,nage);
  dvar_vector ttt(min_tag_year,nyears);
  ttt.initialize();
  int iy;
  for (iy=min_tag_year;iy<=nyears;iy++)
  {
    for (int ir=1;ir<=num_regions;ir++)
      for (int j=1;j<=nage;j++)
        ttt(iy)+=mfexp(pooledtagN(ir,iy,j));
  }
  //ofstream ofs("ptagrep");
  ofs << "pooled tags" << endl;
  for (iy=min_tag_year;iy<=nyears;iy++)
  {
    ofs << "year " << iy << "  ";
    ofs << ttt(iy) << endl;
  }
  {
    int it,ir;
    /*
    ofstream ofslog("log");
    for (it=1;it<=num_tag_releases;it++)
    {
      ofslog << "ag group " << it << endl;
      for (ir=1;ir<=num_regions;ir++)
      {
        ofslog << "Region " << ir << endl;
        ofslog << setfixed()<< setprecision(0)
           << setw(4) << exp(tagN(it,ir)) << endl << endl;
      }
      ofslog << endl;
    }

    for (ir=1;ir<=num_regions;ir++)
      ofslog << setfixed()<< setprecision(0) << setw(4)
             << mfexp(pooledtagN(ir)) << endl << endl;
      ofslog << endl;
    */
  }
}

