/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"
void check_region_number(int tag_region,int nu,int it,int yr,int mn,int nr);
void check_for_early_returns(const i3_array& tr,
  const ivector& true_tag_year,const ivector& true_tag_month,int ntr,
  ivector& itr);

int within_interval(int ir,int l,int u)
{
  if (ir>=l && ir<=u)
    return 1;
  else
    return 0;
}

void sanity_check_2(ivector& initial_tag_year,imatrix& year,
  imatrix& initial_tag_period)
{
  int mmin=initial_tag_year.indexmin();
  int mmax=initial_tag_year.indexmax();
  //int rmin=year.indexmin();
  //int rmax=year.indexmax();
  for (int it=mmin;it<=mmax;it++)
  {
    int rmin=initial_tag_period(it).indexmin();
    int rmax=initial_tag_period(it).indexmax();
    for (int ir=rmin;ir<=rmax;ir++)
    {
      if (initial_tag_year(it) !=  year(ir,initial_tag_period(it,ir)))
      {
        cerr << "sanity error" << endl
             << "initial_tag_year(" << it << ") = " << initial_tag_year(it)
             << endl
             << "while year(" << ir << "," << initial_tag_period(it,ir)
             << ") = " << year(ir,initial_tag_period(it,ir)) << endl;
        //ad_exit(1);
      }
    }
  }
}

extern ofstream * clogf;

void dvar_len_fish_stock_history::read_tagging_data(adstring& root,
   pmulti_species_data & pmsd) 
{
  int ir;
  adstring tmpstring=root+adstring(".tag");
  cifstream cifs((char *)(tmpstring));

  if (!cifs)
  {
    *clogf << "Error trying to open file " << tmpstring << endl;
    cerr << "Error trying to open file " << tmpstring << endl;
    exit(1);
  }
  //*clogf << "BD" << endl;
  int tmpint=num_tag_releases;
  cifs >> num_tag_releases >> tag_shlen >> tag_nlint >> tag_filen;
  if (num_tag_releases != tmpint)
  {
    cerr << "Error number of tag releases in frq file is not the "
         << "same as in tag file " << endl
         << " frq file is " << tmpint << " while tag file is "
         << num_tag_releases << endl;
    ad_exit(1);
  } 
  ivector tmp_itr(1,num_tag_releases);
  cifs >> tmp_itr;
  if (!cifs)
  {
    *clogf << "Error reading tag recovery vector file " << tmpstring << endl;
    cerr << "Error reading tag recovery vector file " << tmpstring << endl;
    ad_exit(1);
  }
  dmatrix tmp_initial_tag_release_by_length(1,num_tag_releases,1,tag_nlint);
  ivector tmp_tag_region(1,num_tag_releases);
  ivector tmp_tag_year(1,num_tag_releases);
  ivector tmp_tag_month(1,num_tag_releases);

  rep_rate.allocate(1,num_regions,1,num_fish_periods,1,num_fish_incidents);
  rep_dev_coffs.allocate(1,num_fisheries,2,num_fish_times);
  itind.allocate(1,num_fisheries);

  ivector tmpitr(1,num_tag_releases);
  tmpitr=tmp_itr;
  {
    for (int i=1;i<=num_tag_releases;i++) if (tmpitr(i)==0) tmpitr(i)=1;
  } 
  i3_array tmp_tag_recaptures_by_length(1,num_tag_releases,1,index_type(tmpitr),1,5);
  int it;
  int isp=1;
  int tagsum=0;
  if (pmsd)
  {
    pmsd->num_tag_release_by_species.allocate(1,pmsd->num_species);
    pmsd->tag_species_flag.allocate(1,num_tag_releases,1,pmsd->num_species);
    pmsd->num_tag_release_by_species.initialize();
    pmsd->tag_species_flag.initialize();
    pmsd->tag_index=0;
  }
  for (it=1;it<=num_tag_releases;it++)
  {
    cifs >> tmp_tag_region(it) >> tmp_tag_year(it) >> tmp_tag_month(it);
    if (pmsd)
    {
      cifs >> pmsd->tag_species_flag(it);
  
      pmsd->num_tag_release_by_species+=pmsd->tag_species_flag(it);
    }
    cifs >> tmp_initial_tag_release_by_length(it);
    if (tmp_itr(it)>0)
    {
      cifs >> tmp_tag_recaptures_by_length(it);
    } 
      // reset the time for tag returns
  }
  ivector iind;
  old_num_tag_releases=num_tag_releases;
  //ivector offset;
  if (pmsd)
  {
    offset.allocate(1,pmsd->num_species);
    offset.initialize();
    for (int is=2;is<=pmsd->num_species;is++)
    {
      offset(is)=offset(is-1)+pmsd->num_tag_release_by_species(is-1);
    }
    num_tag_releases=sum(pmsd->num_tag_release_by_species);
    itr.allocate(1,num_tag_releases);
    iind.allocate(1,pmsd->num_species);
    iind.initialize();
    for (it=1;it<=old_num_tag_releases;it++)
    {
      for (int is=1;is<=pmsd->num_species;is++)
      {
        if (pmsd->tag_species_flag(it,is))
        {
          int ioff=offset(is)+iind(is)+1;
          itr(ioff)=tmp_itr(it);
          iind(is)++;
        }
      }
    }
  }
  else
  {
    itr=tmp_itr;
  }
  ivector   itr2(1,num_tag_releases);
  itr2=itr;
  {
    for (int i=1;i<=num_tag_releases;i++) if (itr2(i)==0) itr2(i)=1;
  } 
  initial_tag_recruitment_period.allocate(1,num_tag_releases,1,num_regions);
  tag_region.allocate(1,num_tag_releases);
  tag_year.allocate(1,num_tag_releases);
  tag_month.allocate(1,num_tag_releases);
  true_tag_year.allocate(1,num_tag_releases);
  true_tag_month.allocate(1,num_tag_releases);
  tag_recaptures_by_length.allocate(1,num_tag_releases,1,index_type(itr2),1,5);
  initial_tag_release_by_length.allocate(1,num_tag_releases,1,tag_nlint);
  if (pmsd)
  {
    //offset.allocate(1,pmsd->num_species);
    //offset.initialize();
    iind.initialize();
    pmsd->num_tag_groups=pmsd->num_tag_release_by_species;
    int mmin=pmsd->tag_species_pointer.indexmin();
    int mmax=pmsd->tag_species_pointer.indexmax();
    ivector old_tag_species_pointer(mmin,mmax);
    old_tag_species_pointer=pmsd->tag_species_pointer;
    pmsd->tag_species_pointer.deallocate();
    pmsd->tag_species_pointer.allocate(1,num_tag_releases);
    //tag_flags.deallocate();
    //tag_flags.allocate(1,num_tag_releases,1,10);
  
    pmsd->tag_species_index.allocate(1,num_tag_releases);
    for (it=1;it<=old_num_tag_releases;it++)
    {
      if (!allocated(pmsd->tag_species_index(it)))
      {
        pmsd->tag_species_index(it).allocate(1,sum(pmsd->tag_species_flag(it)));
        pmsd->tag_species_index(it,1)=it;
      }

      for (int is=1;is<=pmsd->num_species;is++)
      {
        if (pmsd->tag_species_flag(it,is))
        {
          int ioff=offset(is)+iind(is)+1;
              
          if (ioff>it)
          {
            if (!allocated(pmsd->tag_species_index(ioff)))
              pmsd->tag_species_index(ioff).allocate(1,sum(pmsd->tag_species_flag(it)));
            if(sum(pmsd->tag_species_flag(it))>1)
            {
              pmsd->tag_species_index(it,2)=ioff;
              pmsd->combined_tags_flag=1;
            }
            pmsd->tag_species_index(ioff,1)=ioff;
            if (sum(pmsd->tag_species_flag(it))>1)
              pmsd->tag_species_index(ioff,2)=it;
          }
          tag_region(ioff)=tmp_tag_region(it);
          tag_year(ioff)=tmp_tag_year(it);
          tag_month(ioff)=tmp_tag_month(it);
          initial_tag_release_by_length(ioff)=tmp_initial_tag_release_by_length(it);
          tag_recaptures_by_length(ioff)=tmp_tag_recaptures_by_length(it);
          //pmsd->tag_species_pointer(ioff)=old_tag_species_pointer(it);
          pmsd->tag_species_pointer(ioff)=is;
          //tag_flags(ioff)=true_tag_flags(it);
          iind(is)++;
        }
      }
    }
    for (it=1;it<=num_tag_releases;it++)
    {
      pmsd->tag_species_index(it)=sort(pmsd->tag_species_index(it));
    }
    //num_tag_releases=sum(pmsd->num_tag_release_by_species);
  }
  else   //NMD 03Nov13
  {
    tag_region=tmp_tag_region;
    tag_year=tmp_tag_year;
    tag_month=tmp_tag_month;
    initial_tag_release_by_length=tmp_initial_tag_release_by_length;
    tag_recaptures_by_length=tmp_tag_recaptures_by_length;
  }   //NMD 03Nov13

  for (it=1;it<=num_tag_releases;it++)
  {
    if (pmsd)
    {
      if (it>tagsum+pmsd->num_tag_groups(isp))
      {
        tagsum+=pmsd->num_tag_groups(isp);
        isp++;
      }
    }
    date_struc new_date;  


    if (pmsd)
    {
      tag_region(it)=tag_region(it)+(isp-1)*pmsd->num_real_regions;
    }

    check_region_number(tag_region(it),num_regions,it,tag_year(it),
      tag_month(it),num_regions);

    true_tag_year(it)=tag_year(it);

    true_tag_month(it)=tag_month(it);

    // normalize years
    {
      normalize_tag_dates(tag_year,tag_month,true_tag_year,
        month_1,year1,direction_flag,month_factor,first_time,it);
    }
    // readin in the number of recaptures by length interval
    if (itr(it)>0)
    {
      //cifs >> tag_recaptures_by_length(it);
      // reset the time for tag returns
      int mmin = tag_recaptures_by_length(it).indexmin();
      int mmax = tag_recaptures_by_length(it).indexmax();
      for (int i=mmin;i<=mmax;i++)
      {
        if (pmsd)
        {   
          tag_recaptures_by_length(it,i,2)=tag_recaptures_by_length(it,i,2)
            +(isp-1)*pmsd->num_real_fisheries;
        }
        int yr = tag_recaptures_by_length(it,i,3)-year1+1;
        int mn = tag_recaptures_by_length(it,i,4);
        // ***************************************************
        // get right month 1

        if (month_1 !=1)
        {
          if (direction_flag==-1)
          {
            if (tag_month(it) >= month_1)
            {
              mn-=(month_1-1);
            }
            else
            {
              yr-=1;
              mn+=12;
              mn-=(month_1-1);
            }
          }  
          else if (direction_flag==1)
          {
            mn+=(13-month_1);
            if (mn>12)
            {
              yr+=1;
              mn-=12;
            }
          }
        }
        // ***************************************************
        date_struc newdate=get_new_time(yr,mn,1,month_factor,first_time);
        tag_recaptures_by_length(it,i,3)=newdate.year;
        tag_recaptures_by_length(it,i,4)=newdate.month;
      }
    }
    if (!cifs)
    {
      cerr << "Error reading tag release dat in file "
        << tmpstring << " for release number " << it << endl;
      ad_exit(1);
    }
  }
  *clogf << "BG" << endl;
  check_for_early_returns(tag_recaptures_by_length,tag_year,
    tag_month,num_tag_releases,itr);
    
  if (min(tag_region)<1 || max(tag_region)>num_regions || !cifs)
  {
    *clogf << "Error reading tag data from file  " << tmpstring << endl;
    cerr << "Error reading tag data from file  " << tmpstring << endl;
    ad_exit(1);
  }
  // get the fishing period corresponding to each tag release;

  tag_region_bounds.allocate(1,2,1,num_tag_releases);
  if (pmsd)
  {
    pmsd->tag_region_bounds.allocate(1,2,1,num_tag_releases);
    pmsd->tag_region_bounds.initialize();
    for (it=1;it<=num_tag_releases;it++)
    {
      int cs=pmsd->tag_species_pointer(it);
      pmsd->tag_region_bounds(1,it)=pmsd->region_bounds(cs,1);
      pmsd->tag_region_bounds(2,it)=pmsd->region_bounds(cs,2);
      tag_region_bounds(1,it)=pmsd->region_bounds(cs,1);
      tag_region_bounds(2,it)=pmsd->region_bounds(cs,2);
    }
  }
  else
  {
    for (it=1;it<=num_tag_releases;it++)
    {
      tag_region_bounds(1,it)=1;
      tag_region_bounds(2,it)=num_regions;
    }
  }
 
  initial_tag_period.allocate(1,num_tag_releases,tag_region_bounds(1),
    tag_region_bounds(2));


  for (it=1;it<=num_tag_releases;it++)
  {
    //for (int ir=1;ir<=num_regions;ir++)
    for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
    {
      int ip;
      int no_match_flag=1;
      for (ip=1;ip<=num_fish_periods(ir);ip++) 
      {
        if (tag_year(it) < year(ir,ip)) 
        {
          no_match_flag=0;
          break;
        }
        if ( tag_year(it) == year(ir,ip) 
           && tag_month(it) <= month(ir,ip))
        {
          //cout << tag_year(it) << " " << year(ir,ip) << endl;
          //cout << tag_month(it) << " " << month(ir,ip) << endl;
          no_match_flag=0;
          break;
        }
      }
      if (no_match_flag)
      {
        *clogf << endl << "Tag group " << it << " released after last fishery "
             << "remove it from tag file " << endl << endl;
        cerr << endl << "Tag group " << it << " released after last fishery "
             << "remove it from tag file " << endl << endl;
        ad_exit(1);
      }
      initial_tag_period(it,ir)=ip;
    }
  }
  if (age_flags(96) && age_flags(121))
  {
    *clogf << "At the moment can not have both age_flags(96) and"
      << endl << " age flags(121) active " << endl;
    cerr << "At the moment can not have both age_flags(96) and"
      << endl << " age flags(121) active " << endl;
    exit(1);
  }
  if (age_flags(96))
  {
//    terminal_tag_period.allocate(1,num_tag_releases,1,num_regions);
    terminal_tag_period.allocate(1,num_tag_releases,tag_region_bounds(1),
      tag_region_bounds(2));     //NMD_jan25-19

    for (it=1;it<=num_tag_releases;it++)
    {
//     for (int ir=1;ir<=num_regions;ir++)
      for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
      {
        int ip;
        // dave july 3 2018
        for (ip=1;ip<=num_fish_periods(ir);ip++) 
        {
          int ip1=min(num_fish_periods(ir),ip+1);
          if ((tag_year(it)+age_flags(96) < year(ir,ip))
            && move_flags(ir,ip1)>0)  break;
          if (( tag_year(it)+age_flags(96) == year(ir,ip) 
           && tag_month(it) <= month(ir,ip)) 
            && move_flags(ir,ip1)>0)  break;
        }
        terminal_tag_period(it,ir)=min(ip,num_fish_periods(ir));
      }
    }
  }
  *clogf << "BI" << endl;
  if (age_flags(121))
  {
    terminal_tag_period.allocate(1,num_tag_releases,1,num_regions);

    for (it=1;it<=num_tag_releases;it++)
    {
      ivector& trl=tag_recaptures_by_length(it,1);
      int mmonth=trl(4);
      int yyear= trl(3);
      for (int ir=1;ir<=num_regions;ir++)
      {
        int ip;
        for (ip=1;ip<=num_fish_periods(ir);ip++) 
        {
          if (yyear < year(ir,ip)) break;
          if ( yyear == year(ir,ip) 
           && mmonth <= month(ir,ip)) break; 
        }
        terminal_tag_period(it,ir)=min(ip,num_fish_periods(ir));
      }
    }
  }
  minimum_initial_tag_period.allocate(1,num_regions);
  initial_tag_year.allocate(1,num_tag_releases);
  initial_tag_year=tag_year;
  int bad_tag_flag=0;
  for (it=1;it<=num_tag_releases;it++)
  {
    if (initial_tag_year(it)> nyears)
    {
      cout << "tags released after nyears for tag release " << it 
           << " tags released in year "
             << initial_tag_year(it) << "total number of years is "
             << nyears << endl;
      bad_tag_flag=1;
    }
  }
  if (bad_tag_flag) exit(1);
    
  minimum_initial_tag_period=100000;
  for (ir=1;ir<=num_regions;ir++)
  {
    //minimum_initial_tag_period(ir)=initial_tag_period(1,ir);
    for (int it=1;it<=num_tag_releases;it++)
    {
      if (!pmsd)
      {
        if (minimum_initial_tag_period(ir)>initial_tag_period(it,ir)) 
          minimum_initial_tag_period(ir)=initial_tag_period(it,ir);
      }
      else
      {
        int cs=pmsd->tag_species_pointer(it);
        if (pmsd->region_bounds(cs,1)<=ir 
          && pmsd->region_bounds(cs,2)>=ir)
        {
          if (minimum_initial_tag_period(ir)>initial_tag_period(it,ir)) 
            minimum_initial_tag_period(ir)=initial_tag_period(it,ir);
        }
      }
    }
  }

  //min_tag_age.allocate(1,num_tag_releases,1,num_regions,initial_tag_period,
  //    imatrix(1,num_tag_releases,num_fish_periods));
    imatrix tmp=imatrix(1,num_tag_releases,
      tag_region_bounds(1),tag_region_bounds(2),
       num_fish_periods);

    

  min_tag_age.allocate(1,num_tag_releases,
    tag_region_bounds(1),tag_region_bounds(2),initial_tag_period,
    imatrix(1,num_tag_releases,
      tag_region_bounds(1),tag_region_bounds(2),
       num_fish_periods));

    
  for (it=1;it<=num_tag_releases;it++)
  {
    int ng=get_nage_tag(it);
    for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
    {
      int min_age=1;
      min_tag_age(it,ir,initial_tag_period(it,ir))=min_age;
      for (int ip=initial_tag_period(it,ir);ip<=num_fish_periods(ir);ip++)
      {
        min_tag_age(it,ir,ip)=
          min(ng,min_age+year(ir,ip)-initial_tag_year(it));
      }
    }
  }
  nage_by_tag_release.allocate(1,num_tag_releases);
  for (int it=1;it<=num_tag_releases;it++)
  {
    nage_by_tag_release(it)=nage_by_region(tag_region(it));
  }
    

  min_tag_age1.allocate(1,num_tag_releases,1,num_regions,
      index_type(initial_tag_year),nyears);

  for (it=1;it<=num_tag_releases;it++)
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      int min_age=1;
      min_tag_age1(it,ir,initial_tag_year(it))=min_age;
      for (int iy=initial_tag_year(it)+1;iy<=nyears;iy++)
      {
        if (min_age<nage_by_tag_release(it)) min_age++;
        min_tag_age1(it,ir,iy)=min_age;
      }
    }
  }
    
  sanity_check_2(initial_tag_year,year,initial_tag_period);
  if (!age_flags(96)) 
  {
    cout << " initial_tag_year(2) "  << initial_tag_year(2)  << endl;

    terminal_tag_period.allocate(1,num_tag_releases,
      tag_region_bounds(1),tag_region_bounds(2));
    //   tagN.allocate(1,num_tag_releases,1,num_regions,
    tag_num_fish_periods=imatrix(1,num_tag_releases,
      tag_region_bounds(1),tag_region_bounds(2),num_fish_periods);

    tagN.allocate(1,num_tag_releases,
      tag_region_bounds(1),
      tag_region_bounds(2),
      index_type(initial_tag_year),nyears,min_tag_age1,
      nage_by_tag_release);

    tagnum_fish.allocate(1,num_tag_releases,
      tag_region_bounds(1),
      tag_region_bounds(2),
      initial_tag_period,
      imatrix(1,num_tag_releases,tag_region_bounds(1),tag_region_bounds(2),
      num_fish_periods),
      min_tag_age,
      nage_by_tag_release);
    if (parest_flags(360))
    {
      if (!allocated(tag_tot_mort))
      {
        tag_tot_mort.allocate(1,num_tag_releases,
          tag_region_bounds(1),tag_region_bounds(2),initial_tag_period,
          imatrix(1,num_tag_releases,tag_region_bounds(1),
          tag_region_bounds(2),num_fish_periods),min_tag_age,
          nage_by_tag_release);
      }
      tag_tot_mort.initialize();
    }
    minttp.allocate(1,num_regions);   //NMD_jan30-19
  }
  else
  {
    int it;
    terminal_tag_year.allocate(1,num_tag_releases);
    terminal_tag_year=initial_tag_year+age_flags(96);
    for (it=1;it<=num_tag_releases;it++)
    {
      if (terminal_tag_year(it)>nyears)
        terminal_tag_year(it)=nyears;
      if (terminal_tag_year(it)<initial_tag_year(it))
      {
        cout << "Error terminal tag year for tag group " << it
             << " is less than the initial tag year " << endl;
        ad_exit(1);
      }
    }


  *clogf << "BL" << endl;
    min_tag_age5.allocate(1,num_tag_releases,tag_region_bounds(1),
      tag_region_bounds(2),initial_tag_period,terminal_tag_period);

    for (it=1;it<=num_tag_releases;it++)
    {
      for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
      {
        int min_age=1;
        min_tag_age5(it,ir,initial_tag_period(it,ir))=min_age;
        for (int ip=initial_tag_period(it,ir)+1;
          ip<=terminal_tag_period(it,ir);ip++)
        {
          if (year(ir,ip)>year(ir,ip-1)) 
          {
            if (min_age<nage) min_age++;
          }
          min_tag_age5(it,ir,ip)=min_age;
        }
      }
    }
    min_tag_age6.allocate(1,num_tag_releases,1,num_regions,initial_tag_year,
      terminal_tag_year);

    for (it=1;it<=num_tag_releases;it++)
    {
      for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
      {
        int min_age=1;
        min_tag_age6(it,ir,initial_tag_year(it))=min_age;
        for (int ip=initial_tag_year(it)+1;
          ip<=terminal_tag_year(it);ip++)
        {
          if (year(ir,ip)>year(ir,ip-1)) 
          {
            if (min_age<nage) min_age++;
          }
          min_tag_age6(it,ir,ip)=min_age;
        }
      }
    }

  *clogf << "BM" << endl;

    tagN.allocate(1,num_tag_releases,
      tag_region_bounds(1),tag_region_bounds(2),
      index_type(initial_tag_year),index_type(terminal_tag_year),
      min_tag_age6,nage);
    tagnum_fish.allocate(1,num_tag_releases,
      tag_region_bounds(1),tag_region_bounds(2),
      initial_tag_period,terminal_tag_period,min_tag_age5,nage);
    if (parest_flags(360))
    {
      if (!allocated(tag_tot_mort))
      {
        tag_tot_mort.allocate(1,num_tag_releases,
          tag_region_bounds(1),tag_region_bounds(2),
          initial_tag_period,terminal_tag_period,min_tag_age5,nage);
      }
      tag_tot_mort.initialize();
    }
    minttp.allocate(1,num_regions);
    int ir;
    for (ir=1;ir<=num_regions;ir++)
    {
      minttp(ir)=1000000;
      for (int it=1;it<=num_tag_releases;it++)
      {
        if (within_interval(ir,tag_region_bounds(1,it),tag_region_bounds(2,it)))
        { 
          if (minttp(ir)>terminal_tag_period(it,ir))
            minttp(ir)=terminal_tag_period(it,ir);
        }
      }
    }
    {
      cout << "Minttp report" << endl;
      //minttp(3)+=2;
      for (int ir=1;ir<=num_regions;ir++)
      {
        int mp=minttp(ir);
        //cout << ir << " " << mp << " " << year(ir,mp) 
        //      << " " << month(ir,mp) << endl;
        mp++;
        if (mp>month(ir).indexmax())
          mp=month(ir).indexmax();
          
        cout << ir << " " << mp << " " << year(ir,mp) 
              << " " << month(ir,mp) << endl;
      }
      //ad_exit(1);
    }
    min_tag_year=min(tag_year);
    pooledtagN.allocate(1,num_regions,min_tag_year,nyears,1,nage);
      
     pooled_tagnum_fish.allocate(1,num_regions,minttp+1,num_fish_periods,
       1,nage);
     epooled_tagnum_fish_recr.allocate(1,num_regions,minttp+1,num_fish_periods,
       1,nage);
  } 

  min_init_tag_period.allocate(1,num_regions);
  min_init_tag_period.initialize();
  for (ir=1;ir<=num_regions;ir++)
  {
    minttp(ir)=1000000;
    for (int it=1;it<=num_tag_releases;it++)
    {
      if (within_interval(ir,tag_region_bounds(1,it),tag_region_bounds(2,it)))
      {
        if (minttp(ir)>terminal_tag_period(it,ir))
          minttp(ir)=terminal_tag_period(it,ir);
      }
    }
  }


  for (int it=1;it<=num_tag_releases;it++)
  {
    for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
    {
      if (min_init_tag_period(ir)>initial_tag_period(it,ir)||
        min_init_tag_period(ir)==0)
        min_init_tag_period(ir)=initial_tag_period(it,ir);
    }
  }

  get_initial_tag_fishery_realization_index();

  if (!age_flags(96)) 
  {
    num_alltagfish_incidents.allocate(1,num_tag_releases, // 1,num_regions,
      tag_region_bounds(1),tag_region_bounds(2),
      initial_tag_period,tag_num_fish_periods);
    num_tagfish_incidents.allocate(1,num_tag_releases,  // 1,num_regions,
      tag_region_bounds(1),tag_region_bounds(2),
      initial_tag_period,tag_num_fish_periods);

  }
  else
  {
    num_alltagfish_incidents.allocate(1,num_tag_releases,
      tag_region_bounds(1),tag_region_bounds(2),
      initial_tag_period,terminal_tag_period);
    num_tagfish_incidents.allocate(1,num_tag_releases,
      tag_region_bounds(1),tag_region_bounds(2),
      initial_tag_period,terminal_tag_period);
    num_pooledtagfish_incidents.allocate(1,num_regions,
      minttp+1,num_fish_periods);
  }

  *clogf << "BO" << endl;
  for (it=1;it<=num_tag_releases;it++)
  {
    for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
    {
      int ip;
         
      if (initial_tag_period(it,ir))
      {
        if (!age_flags(96)) 
          for (ip=initial_tag_period(it,ir);ip<=num_fish_periods(ir);ip++)
          {
            num_alltagfish_incidents(it,ir,ip)=num_fish_incidents(ir,ip);
            num_tagfish_incidents(it,ir,ip)=num_fish_incidents(ir,ip);
          }
          else
            for (ip=initial_tag_period(it,ir);ip<=terminal_tag_period(it,ir);ip++)
            {
              num_alltagfish_incidents(it,ir,ip)=num_fish_incidents(ir,ip);
              num_tagfish_incidents(it,ir,ip)=num_fish_incidents(ir,ip);
            }
        }
      }
    }
    if (age_flags(96)) 
      for (int ir=1;ir<=num_regions;ir++)
        for (int ip=minttp(ir)+1;ip<=num_fish_periods(ir);ip++)
          num_pooledtagfish_incidents(ir,ip)=num_fish_incidents(ir,ip);

    get_initial_tag_recruitment_period();

  *clogf << "BP" << endl;

  if (!age_flags(96)) 
  {
    tot_tag_catch.allocate(1,num_tag_releases,
      //1,num_regions,
      tag_region_bounds(1),tag_region_bounds(2),
      initial_tag_period,
      //imatrix(1,num_tag_releases,num_fish_periods),
      tag_num_fish_periods,
      1,num_alltagfish_incidents);

    {
      ofstream ofs("obstag");
      ofs << "num_tag_releases" << endl;
      ofs << num_tag_releases << endl;
      ofs << "num_regions" << endl;
      ofs << num_regions << endl;
      ofs << "num_fish_periods" << endl;
      ofs << num_fish_periods << endl;
      ofs << "num_alltagfish_incidents" << endl;
      ofs << num_alltagfish_incidents << endl;
    }


    obstagcatch_by_length.allocate(1,num_tag_releases,
        //1,num_regions,
      tag_region_bounds(1),tag_region_bounds(2),
      initial_tag_period,tag_num_fish_periods,
      //initial_tag_period,imatrix(1,num_tag_releases,num_fish_periods),
      1,num_alltagfish_incidents,1,tag_nlint);

  }
  else
  {
    tot_tag_catch.allocate(1,num_tag_releases,
      tag_region_bounds(1),tag_region_bounds(2),
      initial_tag_period,
      terminal_tag_period,1,num_alltagfish_incidents);

    if (allocated(obstagcatch_by_length))
      obstagcatch_by_length.deallocate();
     
    obstagcatch_by_length.allocate(1,num_tag_releases,
      tag_region_bounds(1),tag_region_bounds(2),
      initial_tag_period,terminal_tag_period,1,num_alltagfish_incidents,
      1,tag_nlint);
  }


  initial_tag_release_by_age.allocate(1,num_tag_releases,1,nage_by_tag_release);


  if (age_flags(96)) 
  {
    min_tag_age4.allocate(1,num_tag_releases,
       tag_region_bounds(1),tag_region_bounds(2),
       initial_tag_period,terminal_tag_period,1,num_tagfish_incidents);
  
    for (it=1;it<=num_tag_releases;it++)
    {
      for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
      {
        int min_age=1;
        min_tag_age4(it,ir,initial_tag_period(it,ir))=min_age;
        for (int ip=initial_tag_period(it,ir)+1;
          ip<=terminal_tag_period(it,ir);ip++)
        {
          if (year(ir,ip)>year(ir,ip-1)) 
          {
            if (min_age<nage) min_age++;
          }
          for (int fi=1;fi<=num_tagfish_incidents(it,ir,ip);fi++)
          {
            min_tag_age4(it,ir,ip,fi)=min_age;
          }
        }
      }
    }
    
  *clogf << "BQ" << endl;

    ivector rmin(1,num_tag_releases);
    ivector rmax(1,num_tag_releases);
    if (pmsd)
    {
      for (it=1;it<=num_tag_releases;it++)
      {
        int cs=pmsd->tag_species_pointer(it);
        rmin(it)=pmsd->region_bounds(cs,1);
        rmax(it)=pmsd->region_bounds(cs,2);
      }
    }
    else
    {
      rmin=1;
      rmax=num_regions;
    }
    i3_array itp=i3_array(1,nage,initial_tag_period);
    i3_array ttp=i3_array(1,nage,terminal_tag_period);
    i4_array nti=i4_array(1,nage,num_tagfish_incidents);

    //probtagcatch.allocate(1,nage,1,num_tag_releases,
    //  1,num_regions,itp,ttp,1,nti);

    probtagcatch.allocate(1,num_tag_releases,
        // 1,num_regions,
      tag_region_bounds(1),tag_region_bounds(2),
      initial_tag_period,terminal_tag_period,
      1,num_tagfish_incidents,1,nage_by_tag_release);

    tagcatch.allocate(1,num_tag_releases,
      //1,num_regions,
      tag_region_bounds(1),tag_region_bounds(2),
      initial_tag_period,terminal_tag_period,
      1,num_tagfish_incidents,min_tag_age4,nage_by_tag_release);

    pooled_tagcatch.allocate(1,num_regions,minttp+1,
      num_fish_periods,
      1,num_pooledtagfish_incidents,1,nage_by_region);

    pooledtot_tag_catch.allocate(1,num_regions,minttp+1,
      num_fish_periods,
      1,num_pooledtagfish_incidents);

    obstagcatch.allocate(1,num_tag_releases,
      //  1,num_regions,
      tag_region_bounds(1),tag_region_bounds(2),
      initial_tag_period,
      terminal_tag_period,
      1,num_tagfish_incidents,min_tag_age4,nage_by_tag_release);

    pooledobstagcatch_by_length.allocate(1,num_regions,
      minttp+1,num_fish_periods,
      1,num_pooledtagfish_incidents,1,tag_nlint);

    pooledobstagcatch.allocate(1,num_regions,
      minttp+1,num_fish_periods,
      1,num_pooledtagfish_incidents,1,nage_by_region);

    obstagcatch1.allocate(1,num_tag_releases,
         //1,num_regions,
      tag_region_bounds(1),tag_region_bounds(2),
      initial_tag_period, terminal_tag_period,
      1,num_tagfish_incidents,min_tag_age4,nage_by_tag_release);

    tagcatch.initialize();

    pooled_tagcatch.initialize();

    pooledtot_tag_catch.initialize();

    obstagcatch.initialize();

    pooledobstagcatch_by_length.initialize();

    pooledobstagcatch.initialize();

    obstagcatch1.initialize();

  }
  *clogf << "BR" << endl;
}

