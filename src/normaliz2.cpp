/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include <admodel.h>
#if !defined(__MSVC32__)
#  pragma implementation "newmult.hpp"
#  pragma implementation "variable.hpp"
#endif
#include "all.hpp"
      static void xxx(int mn,int yr){}

  imatrix normalize_year(fishery_freq_record_array& fra,int month_1,
    int& direction_flag,int& min_year,int& first_data_month)
  {
    int i;
    int maxfish=0;
    for (i=fra.indexmin();i<=fra.indexmax();i++)
    {
      maxfish=max(fra[i].fishery,maxfish);
    }
    imatrix month_detector(1,maxfish,1,12);
    month_detector.initialize();
    for (i=fra.indexmin();i<=fra.indexmax();i++)
    {
      month_detector(fra[i].fishery,fra[i].month)+=1;
    }

    min_year=fra[fra.indexmin()].year;
    int min_month=12;
    direction_flag=0;
    for (i=fra.indexmin()+1;i<=fra.indexmax();i++)
    {
      min_year=min(fra[i].year,min_year);
    }
    // find the minimum month in the first year
    // !!! dave aug 28 01
    // change next line 
    //   old line for (i=fra.indexmin()+1;i<=fra.indexmax();i++)
    //  to this
    for (i=fra.indexmin();i<=fra.indexmax();i++)
    {
      if(min_year==fra[i].year)
      {
        min_month=min(min_month,fra[i].month);
      }
    }
    first_data_month=min_month;

    if (min_year>1)
    {
      int itemp=min_year-1;
      for (i=fra.indexmin();i<=fra.indexmax();i++)
      {
        fra[i].year-=itemp;
      }
    }
    if (month_1 !=1)
    {
      // !!! dave june 5 01
      //if (min_month>month_1)
      if (min_month>=month_1)
      {
        direction_flag=-1;
        for (i=fra.indexmin();i<=fra.indexmax();i++)
        {
          int yr=fra[i].year;
          int mn=fra[i].month;
          if (fra[i].month >= month_1)
          {
            fra[i].month-=(month_1-1);
          }
          else
          {
            fra[i].year-=1;
            fra[i].month+=12;
            fra[i].month-=(month_1-1);
          }
          int yr1=fra[i].year;
          int mn1=fra[i].month;
          xxx(mn1,yr1);
          xxx(mn,yr);
        }
      }        
      else
      {
        direction_flag=1;
        for (i=fra.indexmin();i<=fra.indexmax();i++)
        {
          int yr=fra[i].year;
          int mn=fra[i].month;
          fra[i].month+=(13-month_1);
          if (fra[i].month>12)
          {
            fra[i].year+=1;
            fra[i].month-=12;
          }
          int yr1=fra[i].year;
          int mn1=fra[i].month;
          xxx(mn1,yr1);
          xxx(mn,yr);
        }
      }
    }
    
    //cout <<"normaliz2, min_year min_month: " << min_year << "  " << min_month << endl;
    return month_detector;
  }

  void month_doubling_kludge(fishery_header_record_array& fra,
    int mult,ivector& regmin,int& first_time,pmulti_species_data & pmsd)
  {
    int minr=regmin.indexmin();
    int maxr=regmin.indexmax();
    //if (pmsd)
    //{
    //  maxr=pmsd->num_real_regions;
    //}
    int i1=fra.indexmin();
    int oy=fra[i1].year;
    int om=fra[i1].month;
    int ow=fra[i1].week;
    int minoffset1=48*(oy-1)+4*(om-1)+(ow-1);
    int i;
    for (i=minr;i<=maxr;i++)
    {
      i1=regmin(i);
      oy=fra[i1].year;
      om=fra[i1].month;
      ow=fra[i1].week;
      //cout << oy << " " << om << " " << ow << endl;
      int offset1=48*(oy-1)+4*(om-1)+(ow-1);
      if (offset1<minoffset1) minoffset1=offset1;
    }
    first_time=minoffset1;
  
    int movemult=2;
    imatrix tester(fra.indexmin(),fra.indexmax(),1,2);
    tester.initialize();
    date_struc newdate_movement=get_new_time(oy,om,ow,movemult*mult,first_time);
    date_struc newdate=get_new_time(oy,om,ow,mult,first_time);
    for (i=fra.indexmin();i<=fra.indexmax();i++)
    {
      int oy=fra[i].year;
      int om=fra[i].month;
      int ow=fra[i].week;
      date_struc newdate=get_new_time(oy,om,ow,mult,first_time);
      date_struc newdate_movement=get_new_time(oy,om,ow,movemult*mult,
        first_time);
      fra[i].year=newdate.year;
      fra[i].movement_period=newdate_movement.year;
      fra[i].month=newdate.month;
      fra[i].week=newdate.week;
      tester(i,1)=newdate.year;
      tester(i,2)=newdate_movement.year;
    }
    //ofstream ofs("tester");
    //ofs << tester << endl;
  }

  date_struc get_new_time(int oy,int om,int ow,int mult,int first_time)
  {
    if (!mult) mult=1;
    date_struc newtime;
    int offset=48*(oy-1)+4*(om-1)+(ow-1);
    int newoff=mult*(offset-first_time)+first_time;;
    int year=newoff/48;
    int month=fmod(newoff,48.)/4;
    int week=fmod(fmod(double(newoff),48),4.);
    newtime.year=year+1;
    newtime.month=month+1;
    newtime.week=week+1;
    return newtime;
  }

  date_struc get_old_time(int oy,int om,int ow,int mult,int first_time)
  {
    date_struc newtime;
    int offset=48*(oy-1)+4*(om-1)+(ow-1);
    int newoff=(offset-first_time)/mult+first_time;;
    int year=newoff/48;
    int month=fmod(newoff,48.)/4;
    int week=fmod(fmod(double(newoff),48),4.);
    newtime.year=year+1;
    newtime.month=month+1;
    newtime.week=week+1;
    return newtime;
  }

  void normalize_projection_year(int month_1,int direction_flag,int min_year,
    imatrix& fdat,dvar_fish_stock_history& fsh)
  {
    int i;

    int num_fisheries=fdat(2).indexmax();
    ivector yr=fdat(2);
    ivector mn=fdat(3);
    ivector newyr=fdat(4);
    ivector newmn=fdat(5);
    if (min_year>1)
    {
      int itemp=min_year-1;
      for (i=1;i<=num_fisheries;i++)
      {
        newyr(i)=yr(i)-itemp;
      }
    }
    if (month_1 !=1)
    {
      if (direction_flag==-1)
      {
        for (i=1;i<=num_fisheries;i++)
        {
          if (mn(i) >= month_1)
          {
            newmn(i)=mn(i)-(month_1-1);
          }
          else
          {
            newyr(i)-=1;
            newmn(i)=mn(i)+12;
            newmn(i)-=(month_1-1);
          }
        }
      }        
      else
      {
        for (i=1;i<=num_fisheries;i++)
        {
          newmn(i)=mn(i)+13-month_1;
          if (newmn(i)>12)
          {
            newyr(i)+=1;
            newmn(i)-=12;
          }
        }
      }
    }
    else
    {
      for (i=1;i<=num_fisheries;i++)
      {
        newmn(i)=mn(i);
      }
    }
    if (fsh.month_factor!=0 && fsh.month_factor!=1)
    {         
      date_struc newdate=get_new_time(newyr(1),newmn(1),1,fsh.month_factor,
          fsh.first_time);
        newyr(1)=newdate.year;
        newmn(1)=newdate.month;
      for (i=2;i<=num_fisheries;i++)
      {
        date_struc newdate=get_new_time(newyr(i),newmn(i),1,fsh.month_factor,
          fsh.first_time);
        newyr(i)=newdate.year;
        newmn(i)=newdate.month;
        if (newyr(1) !=newyr(1) || newmn(1) !=newmn(i))
        {
          cerr << "for projections projection years and months must now"
                   " be the same for all fisheries" << endl;
          ad_exit(1);
        } 
      }
    }
    int ip,ir;
    int adjustflag=1;
    for (ir=1;ir<=fsh.num_regions;ir++)
    {
      for (ip=1;ip<=fsh.num_fish_periods(ir);ip++)
      {
        if (fsh.year(ir,ip)==newyr(1) && fsh.month(ir,ip)<newmn(1) ) 
        {
          adjustflag=0;
        }
        if (    fsh.year(ir,ip)> newyr(1) 
          || ( fsh.year(ir,ip)==newyr(1) && fsh.month(ir,ip)>= newmn(1) ) )
        {
          fsh.num_real_fish_periods(ir)= ip-1;
          break;
        }
      }
    }
    fsh.projection_year=newyr;
    fsh.projection_month=newmn;
    fsh.last_real_year=newyr(1)-adjustflag;

    fsh.num_real_fish_times.initialize();


    for (int fi=1;fi<=fsh.num_fisheries;fi++) 
    {
      for (int nt=1;nt<=fsh.num_fish_times(fi);nt++) 
      {
        int ir=fsh.realization_region(fi,nt);
        int ip=fsh.realization_period(fi,nt);
        if ( ip > fsh.num_real_fish_periods(ir))
          break;
        fsh.num_real_fish_times(fi)++;
      }
    } 
  }

  imatrix months_analyzer(imatrix& month_detector)
  {
    int numfish=month_detector.indexmax();
    imatrix months_used(1,numfish);
    for (int ii=1;ii<=numfish;ii++)
    {
      int used=0;
      for (int i=1;i<=12;i++)
      {
        if (month_detector(ii,i)) used++; 
      }
      months_used(ii).allocate(0,used);
      used=0;
      for (int i=1;i<=12;i++)
      {
        if (month_detector(ii,i))  months_used(ii,++used)=i;
      }
      // check intervals
      if (used>1)
      {
        ivector diff(1,used-1);
        diff=months_used(ii)(2,used).shift(1) - months_used(ii)(1,used-1);
        int space=diff(1);
        int badspace_flag=0;
        for (int i=2;i<=used-1;i++)
        {
          if (diff(i) != space)
          {
            badspace_flag=1;
            break;
            cerr << "Error ! months not equally space. This option is not"
                    " implemented yet"  << endl;
            ad_exit(1);
          }
        }
        if (!badspace_flag)
          months_used(ii,0)=1;
        else
          months_used(ii,0)=-1;
      }
      else
      {
        months_used(ii,0)=0;
      }
    }
    return months_used;
  }
