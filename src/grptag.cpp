/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"


imatrix group_pointers(ivector& v)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  int vmin=min(v);
  int vmax=max(v);
  int i,j;
  if (vmin !=1)
  {
    cerr << "Error in grouping vector" << endl;
    exit(1);
  }
  ivector tmp(vmin,vmax);
  tmp.initialize();
  for (i=vmin;i<=vmax;i++)
  {
    for (j=mmin;j<=mmax;j++)
    {
      if (v(j)==i)
        tmp(i)++;
    }
  }
  for (i=vmin;i<=vmax;i++)
  {
    if (!tmp(i))
    {
      cerr << "Error in grouping vector" << endl;
      exit(1);
    }
  }

  imatrix M(1,vmax,1,tmp);
  tmp.initialize();
  for (i=vmin;i<=vmax;i++)
  {
    for (j=mmin;j<=mmax;j++)
    {
      if (v(j)==i)
      {
        tmp(i)++;
        M(i,tmp(i))=j;
      }
    }
  }
  return M;
}

/*
int  get_minimum_time(ivector& yr,ivector&  mo,ivector&  wk,ivector& nt)
{
  int nfish=yr.indexmax();
  dvector times(1,nfish);
  int i=1;
  for (i=1;i<=nfish;i++)
  {
    times(i)=336*yr(i)+28*mo(i)+7*wk(i);
  }
  int imin=1;
  for (i=2;i<=nfish;i++)
  {
    if (times(i)>times(i-1))
      imin=i;
  }
  return imin;
}



void dvar_fish_stock_history::do_grouped_tags_between_time_calculations(void)
{
  ivector grouping=column(fish_flags,32);
  imatrix gtr=group_pointers(grouping);
  int num_gtr_fisheries=gtr.indexmax();
  MY_DOUBLE_TYPE weekdays=7.0;
  MY_DOUBLE_TYPE monthdays=7.0*4.0;
  MY_DOUBLE_TYPE yeardays=12.0*7.0*4.0;
  int tmult=1;
  if (age_flags(57)) tmult= age_flags(57);
  // count the number of realizations of grouped fisheries
  for (int i=1;i<=num_gtr_fisheries;i++)
  {
    int j;
    int ic=0;
    // nf is the number of fisheries in  the current group
    int nf=gtr(i).indexmax();
    ivector rr(1,nf);
    ivector nt(1,nf);
    ivector rp(1,nf);
    ivector yr(1,nf);
    ivector mo(1,nf);
    ivector wk(1,nf);
    nt=1;
    do
      for (j=1;j<=nf;j++)
      {
        // get date for each fishing relaization
        if (nt(j)<=num_fish_times(j))
        {
          rr(j)=realization_region(i,nt(j));
          rp(j)=realization_period(i,nt(j));
          yr(j)=year(rr(j),rp(j));
          mo(j)=month(rr(j),rp(j));
          wk(j)=week(rr(j),rp(j));
        }
        else
        {
  	yr(j)=1000000;
  	mo(j)=0;
  	wk(j)=0;
        }
      }
      // this fishery occurred first
      int jmin=get_minimum_time(yr,mo,wk,nt);
      for (j=1;j<=nf;j++)
      {
        if (same_calendar_date(j,jmin)
        {
          nt(j)++
        }
      }
      ic++;
      // check to see if we are finished
      int cont_flag=0;
      for (j=1;j<=nf;j++)
      {
        if (nt(j)<=num_fish_time(j))
        cont_flag=1;
        break;
      }  	
      if (!cont_flag) break;
    }
    while(1);  
    num_grouped_tag_fish_times(i)=ic;
  }  


}


    ivector nt(1,nf);
    nt=2;  
    if (fish_flags(i,23)==0)
    {
      	
      do
      {
        int rr2=fsh.realization_region(i,nt);
        int rr1=fsh.realization_region(i,nt-1);
        int rp2=fsh.realization_period(i,nt);
        int rp1=fsh.realization_period(i,nt-1);
        fsh.between_times(i,nt)=
          yeardays*(fsh.year(rr2,rp2)-fsh.year(rr1,rp1));
        fsh.between_times(i,nt)+=
          monthdays*(fsh.month(rr2,rp2)-fsh.month(rr1,rp1));
        fsh.between_times(i,nt)+=
          weekdays*(fsh.week(rr2,rp2)-fsh.week(rr1,rp1));
        // inverse of the time period for use in penalty function
        // results are in 1/yr
        fsh.between_times(i,nt)=(tmult*yeardays)/fsh.between_times(i,nt);
      }
      while(1);
    }
    else
    {
      MY_DOUBLE_TYPE tmp=0.0;
      MY_DOUBLE_TYPE t_interval = fsh.fish_flags(i,23)*tmult*monthdays;
      for (int nt=2;nt<=fsh.num_fish_times(i);nt++)
      {
        int rp2=fsh.realization_period(i,nt);
        int rp1=fsh.realization_period(i,nt-1);
        int rr2=fsh.realization_region(i,nt);
        int rr1=fsh.realization_region(i,nt-1);
        tmp+=yeardays*(fsh.year(rr2,rp2)-fsh.year(rr1,rp1));
        tmp+=monthdays*(fsh.month(rr2,rp2)-fsh.month(rr1,rp1));
        tmp+=weekdays*(fsh.week(rr2,rp2)-fsh.week(rr1,rp1));
        // inverse of the time period for use in penalty function
        // results are in 1/yr
        if (tmp>t_interval)
        {
          fsh.between_times(i,nt)=(tmult*yeardays)/tmp;
          tmp=0.0;
        }
        else
        {
          fsh.between_times(i,nt)=0.0;
        }
      }
    }
  }
  ofstream ofs("tt.1");
  ofs << fsh.between_times << endl;
}
*/
