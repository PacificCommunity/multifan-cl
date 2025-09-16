/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


#include "all.hpp"
date_struc normalize_projection_time(int year,int month,
  int week,int month_1,int direction_flag);

MY_DOUBLE_TYPE something_date(int year,int month,int week)
{
  const MY_DOUBLE_TYPE weekdays=7.0;
  const MY_DOUBLE_TYPE monthdays=7.0*4.0;
  const MY_DOUBLE_TYPE yeardays=12.0*7.0*4.0;
  MY_DOUBLE_TYPE it=yeardays*year
    + monthdays*month
    + weekdays*week;
  return it;
}

MY_DOUBLE_TYPE get_node_date(int i,int nt,dvar_fish_stock_history& fsh)
{
  const MY_DOUBLE_TYPE weekdays=7.0;
  const MY_DOUBLE_TYPE monthdays=7.0*4.0;
  const MY_DOUBLE_TYPE yeardays=12.0*7.0*4.0;
  int rr=fsh.realization_region(i,nt);
  int rp=fsh.realization_period(i,nt);
  return something_date(fsh.year(rr,rp),fsh.month(rr,rp),
    fsh.week(rr,rp));
}

MY_DOUBLE_TYPE get_node_date_for_region_and_period(int rr,int rp,
  const dvar_fish_stock_history& fsh)
{
  const MY_DOUBLE_TYPE weekdays=7.0;
  const MY_DOUBLE_TYPE monthdays=7.0*4.0;
  const MY_DOUBLE_TYPE yeardays=12.0*7.0*4.0;
  MY_DOUBLE_TYPE it=yeardays*fsh.year(rr,rp)
    + monthdays*fsh.month(rr,rp)
    + weekdays*fsh.week(rr,rp);
  return it;
}

void check_between_times(const ivector& iv,const imatrix& gf)
{
  int mmin=gf.indexmin();
  int mmax=gf.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    int jmin=gf(i).indexmin();
    int jmax=gf(i).indexmax();
    for (int j=jmin+1;j<=jmax;j++)
    {
      if (iv(gf(i,jmin))!=iv(gf(i,j)))
      cerr << "grouping fish flags 29 inconsistent with between times flags "
              "23" << endl;
    }
  }
}

void dvar_fish_stock_history::do_grouped_between_time_calculations(void)
{
  const MY_DOUBLE_TYPE weekdays=7.0;
  const MY_DOUBLE_TYPE monthdays=7.0*4.0;
  const MY_DOUBLE_TYPE yeardays=12.0*7.0*4.0;
  int tmult=1;
  if (age_flags(57)) tmult= age_flags(57);
  // calculate the nodes for grouped fisheries
  // num_in_group(i) is the number of fisheries in each group
  // gfish_index(i) is a vector of the fisheries in the i'th grouping
  //  ivector nt(1,num_in_group(i)) points to the index of the current 
  // realization of each fishery in the group
  ivector grouping=column(fish_flags,29);
  //cout <<"grpcatch.cpp " << grouping << endl;
  int min_group=min(grouping);
  if (min_group !=1)
  {
    cerr << "min group for fisheries must be 1" << endl;
    ad_exit(1);
  }
  ngroups=max(grouping);
  int i;
  ivector num_in_group(1,ngroups);
  num_in_group.initialize();

  for (i=1;i<=num_fisheries;i++)
    num_in_group(grouping(i))++;

  for (i=1;i<=ngroups;i++)
  {
    if (num_in_group(i)==0)
    {
      cout << "empty group from grouping ff(*,29) for group " 
        << i << endl;
      cout << "the grouping flag values are" << endl  << grouping << endl;
      ad_exit(1);
    }
  }
  //cout <<"grpcatch.cpp " << num_in_group << endl;
  ivector ii_tmp(1,ngroups);
  ivector real_ii_tmp(1,ngroups);
  ii_tmp=1;
  gfish_index.allocate(1,ngroups,1,num_in_group);
  for (i=1;i<=num_fisheries;i++)
    gfish_index(grouping(i),ii_tmp(grouping(i))++ )=i;
    
  //cout <<"grpcatch.cpp " << gfish_index << endl;

  dvector tmp_grouped_fish_time(1,ngroups);
  ii_tmp=0;
  real_ii_tmp=0;
  gfish_ptr.allocate(1,num_fisheries,1,num_fish_times);
  real_gfish_ptr.allocate(1,num_fisheries,1,num_real_fish_times);
  ivector ff23=column(fish_flags,23);
  if (sum(ff23)) check_between_times(ff23,gfish_index);
  for (i=1;i<=ngroups;i++)
  {
    // get the first time period
    dvector init_time_vector(1,num_in_group(i));
    init_time_vector.initialize();
    ivector nt(1,num_in_group(i));
    int j;
    for (j=1;j<=num_in_group(i);j++) nt(j)=1;

    int break_flag;
    // now we do all this twice to get the sizes for the dmatrix
    // grouped_fish_time which holds the times for the nodes for
    // grouped fisheries
    do 
    {
      break_flag=0;
      init_time_vector=1.e+60;
      for (j=1;j<=num_in_group(i);j++)
      {
        if (nt(j)<=num_fish_times(gfish_index(i,j)))
        {
          break_flag=1;
          init_time_vector(j)=get_node_date(gfish_index(i,j),nt(j),*this);
        }
      }
      
      if (break_flag)
      {
        tmp_grouped_fish_time(i)=min(init_time_vector);
        ++ii_tmp(i);
        for (j=1;j<=num_in_group(i);j++)
        {
          if (nt(j)<=num_fish_times(gfish_index(i,j)))
          {
            break_flag=1;
            if (get_node_date(gfish_index(i,j),nt(j),*this)==
              tmp_grouped_fish_time(i))
            {
              gfish_ptr(gfish_index(i,j),nt(j)++)=ii_tmp(i);
            }    
          }
        }
      }
    }
    while(break_flag);
  }
  for (i=1;i<=ngroups;i++)
  {
    // get the first time period
    dvector init_time_vector(1,num_in_group(i));
    init_time_vector.initialize();
    ivector nt(1,num_in_group(i));
    int j;
    for (j=1;j<=num_in_group(i);j++) nt(j)=1;

    int break_flag;
    // now we do all this twice to get the sizes for the dmatrix
    // grouped_fish_time which holds the times for the nodes for
    // grouped fisheries
    do 
    {
      break_flag=0;
      init_time_vector=1.e+60;
      for (j=1;j<=num_in_group(i);j++)
      {
        if (nt(j)<=num_real_fish_times(gfish_index(i,j)))
        {
          break_flag=1;
          init_time_vector(j)=get_node_date(gfish_index(i,j),nt(j),*this);
        }
      }
      
      if (break_flag)
      {
        tmp_grouped_fish_time(i)=min(init_time_vector);
        ++real_ii_tmp(i);
        for (j=1;j<=num_in_group(i);j++)
        {
          if (nt(j)<=num_real_fish_times(gfish_index(i,j)))
          {
            break_flag=1;
            if (get_node_date(gfish_index(i,j),nt(j),*this)==
              tmp_grouped_fish_time(i))
            {
              real_gfish_ptr(gfish_index(i,j),nt(j)++)=real_ii_tmp(i);
            }    
          }
        }
      }
    }
    while(break_flag);
  }
  //cout <<"grpcatch.cpp " << ii_tmp << endl;
  num_real_grouped_fish_times.allocate(1,ngroups);
  num_grouped_fish_times.allocate(1,ngroups);
  num_grouped_fish_times=ii_tmp;
  num_real_grouped_fish_times=real_ii_tmp;

  grouped_catchability_coffs.allocate(1,ngroups,1,num_real_grouped_fish_times);
  grouped_catchability_coffs.initialize();
  grouped_catchability.allocate(1,ngroups,1,num_grouped_fish_times);
  grouped_between_times.allocate(1,ngroups,2,num_grouped_fish_times);
  grouped_catch_dev_coffs.allocate(1,ngroups,2,num_grouped_fish_times);
  grouped_fish_time.allocate(1,ngroups,1,num_grouped_fish_times);
  grouped_year.allocate(1,ngroups,1,num_grouped_fish_times);
  grouped_month.allocate(1,ngroups,1,num_grouped_fish_times);
  grouped_week.allocate(1,ngroups,1,num_grouped_fish_times);
  grouped_true_year.allocate(1,ngroups,1,num_grouped_fish_times);
  grouped_true_month.allocate(1,ngroups,1,num_grouped_fish_times);
  grouped_true_week.allocate(1,ngroups,1,num_grouped_fish_times);

  gp_year.allocate(1,ngroups);
  gp_month.allocate(1,ngroups);
  dvector grouped_projection_date(1,ngroups);
  /*
  for (i=1;i<=num_fisheries;i++)
  {
    int tyr=data_fish_flags(2,i);
    int tym=data_fish_flags(3,i);
    if (tyr==0)
    {
      tyr=year(realization_region(i,1),realization_period(i,1));
      tym=month(realization_region(i,1),realization_period(i,1));
      data_fish_flags(2,i)=tyr+year1-1;
      data_fish_flags(3,i)=tym;
    }
  }    
  */

  for (i=1;i<=ngroups;i++)
  {
    int mmin=gfish_index(i).indexmin();
    int mmax=gfish_index(i).indexmax();
    gp_year(i)=data_fish_flags(2,gfish_index(i,mmin))-year1+1;
    gp_month(i)=data_fish_flags(3,gfish_index(i,mmin));

    date_struc nd1=normalize_projection_time(gp_year(i),
      gp_month(i),1,month_1,direction_flag);

    gp_year(i)=nd1.year;
    gp_month(i)=nd1.month;

    for (int j=mmin+1;j<=mmax;j++)
    {
      if (data_fish_flags(2,gfish_index(i,j))>gp_year(i))
        gp_year(i)=data_fish_flags(2,gfish_index(i,j))-year1+1;

      if (data_fish_flags(3,gfish_index(i,j))>gp_month(i))
        gp_month(i)=data_fish_flags(3,gfish_index(i,j));
    }
    int tmult=1;
    if (age_flags(57)) tmult= age_flags(57);
    date_struc nd=
      get_new_time(gp_year(i),gp_month(i),1,tmult,first_time);
    gp_year(i)=nd.year;
    gp_month(i)=nd.month;
  }
  ii_tmp=0;
  for (i=1;i<=ngroups;i++)
  {
    // get the first time period
    dvector init_time_vector(1,num_in_group(i));
    init_time_vector.initialize();
    ivector nt(1,num_in_group(i));
    int j;
    for (j=1;j<=num_in_group(i);j++) nt(j)=1;

    int break_flag;
    do 
    {
      break_flag=0;
      init_time_vector=1.e+60;
      for (j=1;j<=num_in_group(i);j++)
      {
        if (nt(j)<=num_fish_times(gfish_index(i,j)))
        {
          break_flag=1;
          init_time_vector(j)=get_node_date(gfish_index(i,j),nt(j),*this);
        }
      }
      if (break_flag)
      {
        grouped_fish_time(i,++ii_tmp(i))=min(init_time_vector);
        for (j=1;j<=num_in_group(i);j++)
        {
          if (nt(j)<=num_fish_times(gfish_index(i,j)))
          {
            break_flag=1;
            if (get_node_date(gfish_index(i,j),nt(j),*this)==
              grouped_fish_time(i,ii_tmp(i)))
            {
              gfish_ptr(gfish_index(i,j),nt(j)++)=ii_tmp(i);
            }    
          }
        }
      }
    }
    while(break_flag);
  }

  for (i=1;i<=num_fisheries;i++)
  {
    int rr=realization_region(i,1);
    for (int nt=1;nt<=num_fish_times(i);nt++)
    {
      int rr=realization_region(i,nt);
      int rp=realization_period(i,nt);
      grouped_year(grouping(i),gfish_ptr(i,nt))=year(rr,rp);
      grouped_month(grouping(i),gfish_ptr(i,nt))=month(rr,rp);
      grouped_week(grouping(i),gfish_ptr(i,nt))=week(rr,rp);
    }
  }

  for (i=1;i<=ngroups;i++)
  {
    //cout <<"grpcatch.cpp " << fish_flags(i,23) << endl; 
    if (fish_flags(gfish_index(i,1),23)==0)
    {
      for (int nt=2;nt<=ii_tmp(i);nt++)
      {
        grouped_between_times(i,nt)=
          grouped_fish_time(i,nt)
           -grouped_fish_time(i,nt-1);
          
        // inverse of the time period for use in penalty function
        // results are in 1/yr
        grouped_between_times(i,nt)
          =(tmult*yeardays)/grouped_between_times(i,nt);
      }
    }
    else
    {
      MY_DOUBLE_TYPE tmp=0.0;
      MY_DOUBLE_TYPE t_interval = fish_flags(gfish_index(i,1),23)*tmult*monthdays;
      for (int nt=2;nt<=ii_tmp(i);nt++)
      {
        tmp+=grouped_fish_time(i,nt)-grouped_fish_time(i,nt-1);
        // inverse of the time period for use in penalty function
        // results are in 1/yr
        int yr=grouped_year(i,nt);
        int mn=grouped_month(i,nt);
        int pyr=gp_year(i);
        int pmn=gp_month(i);
        if (pyr>0)
        {
          if ( (tmp>t_interval) && ( (yr<pyr) || ( (yr==pyr) && (mn<pmn))) ) 
          {
            grouped_between_times(i,nt)=(tmult*yeardays)/tmp;
            tmp=0.0;
          }
          else
          {
            grouped_between_times(i,nt)=0.0;
          }
        }
        else
        {
          if ( tmp>t_interval) 
          {
            grouped_between_times(i,nt)=(tmult*yeardays)/tmp;
            tmp=0.0;
          }
          else
          {
            grouped_between_times(i,nt)=0.0;
          }
        }
      }
    }
  }
  //ofstream ofs("tt.1");
  //ofs << grouped_between_times << endl;
}

  void dvar_fish_stock_history::grouped_catchability_calc(void)
  {
    int i;
    ivector grouping=column(fish_flags,29);
 
    for (i=1;i<=ngroups;i++)
    {
      grouped_catchability(i,1)=0.0; 
      for (int nt=2;nt<=num_grouped_fish_times(i);nt++)
      {
        grouped_catchability(i,nt)=grouped_catchability(i,nt-1)
          + grouped_catch_dev_coffs(i,nt);
      }
    }
    const MY_DOUBLE_TYPE tpi=2.e0*3.14159e0;
    dvar3_array * pc=0;
    if (af170q0)
    {
      pc=&catchability_q0;
    }
    else
    {
      pc=&catchability;
    }
    for (i=1;i<=num_fisheries;i++)
    {
      int rr=realization_region(i,1);
      (*pc)(rr,realization_period(i,1),realization_incident(i,1))
        =q0(i);
      for (int nt=2;nt<=num_fish_times(i);nt++)
      {
        int rr=realization_region(i,nt);
        int rp=realization_period(i,nt);
        (*pc)(rr,rp,realization_incident(i,nt))
          =q0(i)+grouped_catchability(grouping(i),gfish_ptr(i,nt));
      }
    }
  }

  
date_struc normalize_projection_time(int year,int month,
  int week,int month_1,int direction_flag)
{
  date_struc newdate;
  if (month_1 !=1)
  {
    //if (min_month>=month_1)
    if (direction_flag==-1)
    {
      if (month >= month_1)
      {
        month-=(month_1-1);
      }
      else
      {
          year-=1;
          month+=12;
          month-=(month_1-1);
      }
    }      
    else
    {
      //direction_flag=1;
      month+=(13-month_1);
      if (month>12)
      {
        year+=1;
        month-=12;
      }
    }
  }
  newdate.year=year;
  newdate.month=month;
  newdate.week=week;
  return newdate;
}



void dvar_fish_stock_history::set_grouped_true_months(void)
{
  int i,nt;
  ivector grouping=column(fish_flags,29);
  for (i=1;i<=num_fisheries;i++)
  {
    int rr=realization_region(i,1);
    for (int nt=1;nt<=num_fish_times(i);nt++)
    {
      int rr=realization_region(i,nt);
      int rp=realization_period(i,nt);
      grouped_true_year(grouping(i),gfish_ptr(i,nt))=true_year(rr,rp);
      grouped_true_month(grouping(i),gfish_ptr(i,nt))=true_month(rr,rp);
      grouped_true_week(grouping(i),gfish_ptr(i,nt))=true_week(rr,rp);
    }
  }
}
