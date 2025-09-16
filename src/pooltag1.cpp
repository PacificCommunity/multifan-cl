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
  dvariable mfexp(prevariable& );
  //dvar_vector mfexp(dvar_vector& );
  //dvector mfexp(dvector& );
#endif

  extern dvar_vector * psv;


void dvar_fish_stock_history::pooled_tag_catch_equations_calc
  (dvar_vector& sv)
{
  MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
  int current_year=min_tag_year;
  ivector rip(1,num_regions);
  rip=minttp+1;
  int ir;
  pooledtagN=-30;

  int icounter=0;
  int jcounter=0;
  ivector add_flag(1,num_regions);
  add_flag.initialize();

  pooled_tagnum_fish=-30;
  do
  {
    ++icounter;
    jcounter=0;
    for (ir=1;ir<=num_regions;ir++)
    {
      int& ip=rip(ir);
      int tip=ip;
      int yr=year(ir,ip);
      //cout << icounter << " " << ++jcounter << endl;
      //if (icounter==37 && jcounter==3)
        //cout << "trap1" << endl;
      if (yr==current_year &&
            ip<num_fish_periods(ir)) 
      {
        if (!add_flag(ir))
        {
          pooled_tagnum_fish(ir,ip)=log(1.e-20+epooled_tagnum_fish_recr(ir,ip)
            +mfexp(pooled_tagnum_fish(ir,ip)));
          add_flag(ir)=1;
        }
        do
        {
          if (ip>=num_fish_periods(ir)) break;
          if (year(ir,ip+1)==current_year)
          {
            pooled_tagnum_fish(ir,ip+1)=
	      log(mfexp(pooled_tagnum_fish(ir,ip)-tot_mort(ir,ip))+
                epooled_tagnum_fish_recr(ir,ip+1));
            //cout <<"pooltag1.cpp " << sum(exp(pooled_tagnum_fish(ir,ip+1))) << endl;

            ip++;
          }
          else
          {
            pooled_tagnum_fish(ir,ip+1,1)=
              log(1.e-20+epooled_tagnum_fish_recr(ir,ip+1,1));

            pooled_tagnum_fish(ir,ip+1)(2,nage-1)=
              log(1.e-10+
              ++(mfexp(pooled_tagnum_fish(ir,ip)(1,nage-2)-
	          tot_mort(ir,ip)(1,nage-2)))
              + epooled_tagnum_fish_recr(ir,ip+1)(2,nage-1));

            pooled_tagnum_fish(ir,ip+1,nage)=
              log(1.e-10 
               + mfexp(pooled_tagnum_fish(ir,ip,nage-1)-tot_mort(ir,ip,nage-1))
               + mfexp(pooled_tagnum_fish(ir,ip,nage)-tot_mort(ir,ip,nage)) 
               + epooled_tagnum_fish_recr(ir,ip+1,nage));
  
            //cout << ir << " " << sum(exp(pooled_tagnum_fish(ir,ip+1))) << endl;

            if (current_year<nyears)
            {
              pooledtagN(ir,current_year+1)(2,nage)=pooled_tagnum_fish(ir,ip+1)(2,nage);
            }

            ip++;
            break;
          }
        }
        while (1);
      }
      else   // there were no fisheries in this region for this year
      {
        if (current_year<nyears)
        {
          for (int j=1;j<nage-1;j++)      // Loop over age classes
          {
            if (!pmsd)
            {
              pooledtagN(ir,current_year+1,j+1)=
                pooledtagN(ir,current_year,j)-
                exp(nat_mort(current_year,j));
            }
            else
            {
              pooledtagN(ir,current_year+1,j+1)=
                pooledtagN(ir,current_year,j)-
                exp(get_nat_mort_region(ir)(current_year,j));
            }
          }
          if (!pmsd)
          {
            pooledtagN(ir,current_year+1,nage)=
              log(1.e-20
                +mfexp(pooledtagN(ir,current_year,nage-1)
                  -exp(nat_mort(current_year,nage-1)))
                +mfexp(pooledtagN(ir,current_year,nage)
                  -exp(nat_mort(current_year,nage))));
          }
          else
          {
            pooledtagN(ir,current_year+1,nage)=
              log(1.e-20
                +mfexp(pooledtagN(ir,current_year,nage-1)
                  -exp(get_nat_mort_region(ir)(current_year,nage-1)))
                +mfexp(pooledtagN(ir,current_year,nage)
                  -exp(get_nat_mort_region(ir)(current_year,nage))));
          }
        }
      }
    }
    if(current_year<nyears)
    {
      // Changed af(57) to af(53) J.H. 27 Nov 01
      if (age_flags(53))
      {
        if ( !((current_year+1)%age_flags(53)) )
        {
          if (num_regions>1) do_the_diffusion(current_year+1,sv,pooledtagN);
        }
      }
      else
      {
        if (num_regions>1) do_the_diffusion(current_year+1,sv,pooledtagN);
      }
    }
    current_year++;
  }
  while (current_year<=nyears);
  //greport("B catch_equations_calc");

  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=minttp(ir)+1;ip<=num_fish_periods(ir);ip++)
    {
      for (int fi=1;fi<=num_pooledtagfish_incidents(ir,ip);fi++)
      {        
        pooled_tagcatch(ir,ip,fi)=mfexp(fish_mort(ir,ip,fi)-log(1.e-10+tot_mort(ir,ip))+
          log(one_plus-survival(ir,ip))+pooled_tagnum_fish(ir,ip));
      }
    }
  }
  //greport("leaving pooled_tag_catch_equations_calc");
#if defined(COUNT_FISH)
  {
    dvar_vector ttt(min_tag_year,nyears);
    int iy;
    for (iy=min_tag_year;iy<=nyears;iy++)
    {
      for (int ir=1;ir<=num_regions;ir++)
        for (int j=1;j<=nage;j++)
          ttt(iy)+=exp(pooledtagN(ir,iy,j));
    }
    ofstream ofs("ptagrep");
    ofs << "pooled tags" << endl;
    for (iy=min_tag_year;iy<=nyears;iy++)
    {
      ofs << "year " << iy << "  ";
      ofs << ttt(iy) << endl;
    }
    dvariable ttc=0.0;
    for (int i1=tagcatch.indexmin();i1<=tagcatch.indexmax();i1++)
      for (int j1=1;j1<=num_regions;j1++)
        for (int k1=tagcatch(i1,j1).indexmin();k1<=tagcatch(i1,j1).indexmax();k1++)
         ttc+=sum(tagcatch(i1,j1,k1));
    
    ofs << "total number of tags caught " << endl;
    ofs << ttc << endl;
    ad_exit(1);
  }
#endif
}

