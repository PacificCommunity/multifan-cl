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

void dvar_fish_stock_history::xpooled_tag_catch_equations_calc
  (dvar_vector& sv)
{
  MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
  ivector rip(1,num_regions);
  rip=minttp+1;
  int ir;
  ivector bug_tmp_mp(1,num_regions);
  ivector tmp_mp(1,num_regions);
  ivector tmp_yr(1,num_regions);
  ivector tmp_mn(1,num_regions);
  ivector tmp_ip(1,num_regions);
  bug_tmp_mp.initialize();
  tmp_yr.initialize();
  tmp_mn.initialize();
  tmp_ip.initialize();
  
  //ofstream xpofs("xpool_movestuff");
  //pooledtagN=-30;
  {
    int mmin=pooledtagN.indexmin();
    int mmax=pooledtagN.indexmax();
    for (int i1=mmin;i1<=mmax;i1++)
    {
      if (allocated(pooledtagN(i1)))
      {
        int mmin=pooledtagN(i1).indexmin();
        int mmax=pooledtagN(i1).indexmax();
        for (int i2=mmin;i2<=mmax;i2++)
        {
          if (allocated(pooledtagN(i1,i2)))
          {
            pooledtagN(i1,i2)=-30.;
          }
        }
      }
    }            
  }
  int finished_flag=0;
  int break_flag=0;
  int nomoveflag=0;

  ivector add_flag(1,num_regions);
  add_flag.initialize();
  {
    int mmin=pooled_tagnum_fish.indexmin();
    int mmax=pooled_tagnum_fish.indexmax();
    for (int i1=mmin;i1<=mmax;i1++)
    {
      if (allocated(pooled_tagnum_fish(i1)))
      {
        int mmin=pooled_tagnum_fish(i1).indexmin();
        int mmax=pooled_tagnum_fish(i1).indexmax();
        for (int i2=mmin;i2<=mmax;i2++)
        {
          if (allocated(pooled_tagnum_fish(i1,i2)))
          {
            pooled_tagnum_fish(i1,i2)=-30.;
          }
        }
      }
    }            
  }
  do
  {
    finished_flag=1;
    for (ir=1;ir<=num_regions;ir++)
    {
      break_flag=0;
      int& ip=rip(ir);
      do
      {

        if (ip>=num_fish_periods(ir)) 
        {
          nomoveflag=1;
          break;
        }
        finished_flag=0;
        if (!add_flag(ir))
        {
          pooled_tagnum_fish(ir,ip)=log(1.e-20+epooled_tagnum_fish_recr(ir,ip)
            +mfexp(pooled_tagnum_fish(ir,ip)));
          add_flag(ir)=1;
        }
        if (year(ir,ip+1)==year(ir,ip))
        {
          pooled_tagnum_fish(ir,ip+1)=
	    log(1.e-20+mfexp(pooled_tagnum_fish(ir,ip)-tot_mort(ir,ip))+
              epooled_tagnum_fish_recr(ir,ip+1));
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
         } 
         if (move_flags(ir,ip+1)) 
         {
           tmp_ip(ir)=ip+1;
           tmp_yr(ir)=year(ir,ip+1);
           tmp_mn(ir)=month(ir,ip+1);
           tmp_mp(ir)=move_index(ir,ip+1);
           bug_tmp_mp(ir)=move_index(ir,ip);
          
           /*
           if (ir==num_regions)
           {
             {
               xpofs << "map = " << tmp_mp(1,ir) << endl;
               xpofs << "year = " << tmp_yr(1,ir) << endl;
               xpofs << "month = " << tmp_mn(1,ir) << endl;
               xpofs << "periods = " << tmp_ip(1,ir) << endl;
               xpofs << endl;
             }
           }
           */
           break_flag=1;
         } 
         ip++;
       }
       while (!break_flag); 
     }
    // move the tags
    {
      if (num_regions>1&& !nomoveflag)
      {
        check_sanity(tmp_mp);
//        if (parest_flags(356))
//          tp=bug_tmp_mp(1);

        dvar_matrix tmp=fast_diffusion_calcs(nage,num_regions,
          pooled_tagnum_fish,
              Dad(tmp_mp(1)),rip,0,pmsd);
        //dvar_matrix tmp=fast_diffusion_calcs(nage,num_regions,
        //  pooled_tagnum_fish,Dad(tp),rip);
        for (int ir=1;ir<=num_regions;ir++)
        {
          pooled_tagnum_fish(ir,rip(ir))=log(5.e-10+tmp(ir));
        }
      }
      else
      {
        cout << "Avoided move" << endl;
      }
    }
  } // need to decide when to quit
  while (!finished_flag);

  /*
  {
    ofstream ofs("xpnumfish");
    dvariable ssum=0.0;
    for (int ir=1;ir<=num_regions;ir++)
    {
      ofs << "region " << ir << endl;
      for (int ip=minttp(ir)+1;ip<=num_fish_periods(ir);ip++)
      {
        ofs << setw(3) << ip;
        ofs << setw(3) << move_flags(ir,ip);
        ofs << setw(5) << year(ir,ip) << setw(4) << month(ir,ip);
        ofs << setw(5) << num_fish_incidents(ir,ip);
        ofs << setw(7) << setfixed() << setprecision(1) 
            << sum(exp(pooled_tagnum_fish(ir,ip)));
        ofs << " " << setw(6) << setfixed() << setprecision(1) 
            << exp(pooled_tagnum_fish(ir,ip)) << endl;
      }
    }
    exit(0);
  }
  */

  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=minttp(ir)+1;ip<=num_fish_periods(ir);ip++)
    {
      for (int fi=1;fi<=num_pooledtagfish_incidents(ir,ip);fi++)
      {        
        pooled_tagcatch(ir,ip,fi)=
          mfexp(fish_mort_calcs(ir,ip,fi)+pooled_tagnum_fish(ir,ip));
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

