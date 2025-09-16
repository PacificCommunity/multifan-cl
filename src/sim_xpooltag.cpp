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

//#define COUNT_FISH
extern MY_DOUBLE_TYPE multiplier;


void dvar_fish_stock_history::sim_xpooled_pd(dmatrix & prob,int & id,
  const char * s)
{
  MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
  ivector rip(1,num_regions);
  rip=sim_minttp+1;
  int ir;
  
   ofstream xpofs("getxsp_movestuff");
  //sim_pooledtagN=-30;
  {
    int xmmin=sim_pooledtagN.indexmin();
    int xmmax=sim_pooledtagN.indexmax();
    for (int i1=xmmin;i1<=xmmax;i1++)
    {
      if (allocated(sim_pooledtagN(i1)))
      {
        int mmin=sim_pooledtagN(i1).indexmin();
        int mmax=sim_pooledtagN(i1).indexmax();
        for (int i2=mmin;i2<=mmax;i2++)
        {
          if (allocated(sim_pooledtagN(i1,i2)))
          {
            sim_pooledtagN(i1,i2)=-30.;
          }
        }
      }
    }            
  }
  int finished_flag=0;
  int break_flag=0;
  int nomoveflag=0;

  ivector tmp_mp(1,num_regions);
  ivector tmp_yr(1,num_regions);
  ivector tmp_mn(1,num_regions);
  ivector tmp_ip(1,num_regions);
  tmp_mp.initialize();
  tmp_yr.initialize();
  tmp_mn.initialize();
  tmp_ip.initialize();

  ivector add_flag(1,num_regions);
  sim_pooled_tagnum_fish.initialize();
  add_flag.initialize();
  {
    int xmmin=sim_pooled_tagnum_fish.indexmin();
    int xmmax=sim_pooled_tagnum_fish.indexmax();
    for (int i1=xmmin;i1<=xmmax;i1++)
    {
      if (allocated(sim_pooled_tagnum_fish(i1)))
      {
        int mmin=sim_pooled_tagnum_fish(i1).indexmin();
        int mmax=sim_pooled_tagnum_fish(i1).indexmax();
        for (int i2=mmin;i2<=mmax;i2++)
        {
          if (allocated(sim_pooled_tagnum_fish(i1,i2)))
          {
            sim_pooled_tagnum_fish(i1,i2)=-30.;
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
          sim_pooled_tagnum_fish(ir,ip)=log(1.e-20+sim_epooled_tagnum_fish_recr(ir,ip)
            +mfexp(sim_pooled_tagnum_fish(ir,ip)));
          add_flag(ir)=1;
        }
        if (year(ir,ip+1)==year(ir,ip))
        {
          sim_pooled_tagnum_fish(ir,ip+1)=
	    log(mfexp(sim_pooled_tagnum_fish(ir,ip)-tot_mort(ir,ip))+
              sim_epooled_tagnum_fish_recr(ir,ip+1));
        }
        else
        {
          sim_pooled_tagnum_fish(ir,ip+1,1)=
            log(1.e-20+sim_epooled_tagnum_fish_recr(ir,ip+1,1));

          sim_pooled_tagnum_fish(ir,ip+1)(2,nage-1)=
            log(1.e-10+
            ++(mfexp(sim_pooled_tagnum_fish(ir,ip)(1,nage-2)-
	        tot_mort(ir,ip)(1,nage-2)))
            + sim_epooled_tagnum_fish_recr(ir,ip+1)(2,nage-1));

          sim_pooled_tagnum_fish(ir,ip+1,nage)=
            log(1.e-10 
             + mfexp(sim_pooled_tagnum_fish(ir,ip,nage-1)-tot_mort(ir,ip,nage-1))
             + mfexp(sim_pooled_tagnum_fish(ir,ip,nage)-tot_mort(ir,ip,nage)) 
             + sim_epooled_tagnum_fish_recr(ir,ip+1,nage));
         } 
         if (move_flags(ir,ip+1))
         {
           tmp_ip(ir)=ip+1;
           tmp_yr(ir)=year(ir,ip+1);
           tmp_mn(ir)=month(ir,ip+1);
           tmp_mp(ir)=move_index(ir,ip+1);
          
           if (ir==num_regions)
           {
             print_movement_stuff(xpofs,tmp_mp,tmp_yr,tmp_mn,tmp_ip,
               num_regions);
             xpofs << endl;

             int diff_flag=0;
             for (int i=2;i<=num_regions;i++)
             {
               if (tmp_mp(i) !=tmp_mp(1))
                 diff_flag=1;
             }

             if (diff_flag)
             {
               cout << "sanity error" << endl;
               cout << "map = " << tmp_mp(1,ir) << endl;
               cout << "year = " << tmp_yr(1,ir) << endl;
               cout << "month = " << tmp_mn(1,ir) << endl;
               ad_exit(1);
             }
           }
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
        dvar_matrix tmp=fast_diffusion_calcs(nage,num_regions,
          sim_pooled_tagnum_fish,Dad(tmp_mp(1)),rip);
        for (int ir=1;ir<=num_regions;ir++)
        {
          sim_pooled_tagnum_fish(ir,rip(ir))=log(5.e-10+tmp(ir));
        }
      }
      else
      {
        cout << "Avoided move" << endl;
      }
    }
  } // need to decide when to quit
  while (!finished_flag);

  cout << sum(mfexp(sim_pooled_tagnum_fish)) << endl;
  /*
  {
    ofstream ofs("xpnumfish");
    dvariable ssum=0.0;
    for (int ir=1;ir<=num_regions;ir++)
    {
      ofs << "region " << ir << endl;
      for (int ip=sim_minttp(ir)+1;ip<=num_fish_periods(ir);ip++)
      {
        ofs << setw(3) << ip;
        ofs << setw(3) << move_flags(ir,ip);
        ofs << setw(5) << year(ir,ip) << setw(4) << month(ir,ip);
        ofs << setw(5) << num_fish_incidents(ir,ip);
        ofs << setw(7) << setfixed() << setprecision(1) 
            << sum(exp(sim_pooled_tagnum_fish(ir,ip)));
        ofs << " " << setw(6) << setfixed() << setprecision(1) 
            << exp(sim_pooled_tagnum_fish(ir,ip)) << endl;
      }
    }
    exit(0);
  }
  */

  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=sim_minttp(ir)+1;ip<=num_fish_periods(ir);ip++)
    {
      for (int fi=1;fi<=sim_num_pooledtagfish_incidents(ir,ip);fi++)
      {        
        sim_pooled_tagcatch(ir,ip,fi)=
          mfexp(fish_mort_calcs(ir,ip,fi)+sim_pooled_tagnum_fish(ir,ip));
        dvector tc=value(sim_pooled_tagcatch(ir,ip,fi));
        int jmin=tc.indexmin();
        MY_DOUBLE_TYPE pp=sum(tc)/multiplier;
        if (pp>1.e-6)
        {
          int nn=0;
          for (int i=jmin;i<=nage;i++)
          {
            if (tc(i)>1.e-6) nn=i;
          }
          prob(++id,1)=pp;
          prob(id,2)=ir;
          prob(id,3)=ip;
          prob(id,4)=fi;
          prob(id,5)=nn;
        }
      }
    }
  }
  ofstream ofscc(s);
  for (int i=1;i<=id;i++)
  {
    ofscc  << " " << setw(10) << setprecision(5) << setscientific() << prob(i,1) 
           << " " << setw(4) << setprecision(5) << setfixed() << ivector(prob(i))(2,5) 
         << endl;
  }
  //cout <<  sum(sim_pooled_tagcatch) << endl;
  //cout <<  sum(sim_pooled_tagcatch) << endl;
  //ad_exit(1);
  //greport("leaving pooled_tag_catch_equations_calc");
#if defined(COUNT_FISH)
  {
    dvar_vector ttt(sim_min_tag_year,nyears);
    ttt.initialize();
    int iy;
    for (iy=sim_min_tag_year;iy<=nyears;iy++)
    {
      for (int ir=1;ir<=num_regions;ir++)
        for (int j=1;j<=nage;j++)
          ttt(iy)+=exp(sim_pooledtagN(ir,iy,j));
    }
    ofstream ofs("sim_ptagrep");
    ofs << "sim_pooled tags" << endl;
    for (iy=sim_min_tag_year;iy<=nyears;iy++)
    {
      ofs << "year " << iy << "  ";
      ofs << ttt(iy) << endl;
    }
    dvariable ttc=0.0;
    for (int i1=sim_tagcatch.indexmin();i1<=sim_tagcatch.indexmax();i1++)
      for (int j1=1;j1<=num_regions;j1++)
        for (int k1=sim_tagcatch(i1,j1).indexmin();k1<=sim_tagcatch(i1,j1).indexmax();k1++)
         ttc+=sum(sim_tagcatch(i1,j1,k1));
    
    ofs << "total number of tags caught " << endl;
    ofs << ttc << endl;
    ad_exit(1);
  }
#endif
}

