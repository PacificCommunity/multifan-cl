/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

int add_first_node(fishery_header_record * min_record,
  const fishery_header_record& _current_record,
   const  movement_info&  _mi,int num_regions);

void breakup_big_function(ivector& nfp,ivector& nft,
  dvar_fish_stock_history& fsh,
  imatrix& nfi,fishery_header_record_array& fhra,imatrix& temp,imatrix& ty,
  imatrix& tm,imatrix& tw,imatrix& test_move_index,
  ivector& header_movement_index);

    header_record& header_record::operator = (header_record& hr)
    {
      year=hr.year;
      month=hr.month;
      week=hr.week;
      return *this;
    }
  
  istream& operator >> (istream& ifs,fishery_header_record_array& fhra)
  {
    for (int i=fhra.indexmin();i<=fhra.indexmax();i++)
    {
      ifs >> fhra(i);
    }
    return ifs;
  }

  cifstream& operator >> (cifstream& ifs,fishery_header_record_array& fhra)
  {
    for (int i=fhra.indexmin();i<=fhra.indexmax();i++)
    {
      ifs >> fhra(i);
    }
    return ifs;
  }

  int movement_info::num_periods(void)
  {
    if (!allocated(move_weeks))
      return 0;
    else
      return move_weeks.indexmax();
  }

  void movement_info::allocate(const ivector & mw, int _year)
  {
    move_weeks.allocate(mw.indexmin(),mw.indexmax());
    move_weeks=mw;
    year=_year;
    current=1;
  }

  int movement_info::month(void)
  {
    return (move_weeks(current)-1)/4+1;
  }

  int movement_info::week(void)
  {
    return (move_weeks(current)-1)%4+1;
  }

  void movement_info::operator ++ (void) 
  { 
    if (current < move_weeks.indexmax())
    {
      current++;
    }
    else
    {
      current = move_weeks.indexmin();
      year++;
    }
  }
  void movement_info::initialize(const fishery_header_record& t1)
  {
    // sets year, period to be >= the fishery_header_record
    int weeks=4*(t1.month-1)+t1.week;
    year=t1.year;
    int i;
    for (i=move_weeks.indexmin();i<=move_weeks.indexmax();i++)
    {
      if (weeks<=move_weeks(i))
        break;
    }
    if (i>move_weeks.indexmax())
    {
      current=1;
      year++;
    }
    else
    {
      current=i;
    }
  }
  void movement_info::initialize_terminal(const fishery_header_record& t1)
  {
    // sets year, period to be >= the fishery_header_record
    int weeks=4*(t1.month-1)+t1.week;
    year=t1.year;
    int i;
    for (i=move_weeks.indexmin();i<=move_weeks.indexmax();i++)
    {
      if (weeks<=move_weeks(i))
        break;
    }
    if (i<=move_weeks.indexmax())
    {
    }
    if (i>move_weeks.indexmax())
    {
      current=1;
      year++;
    }
    else
    {
      current=i;
      if (i<=move_weeks.indexmax())
      {
        ++(*this);
      }
    }
  }

  int get_movement_index(int om,int ow,movement_info& mo)
  {
    int m=4*(om-1)+ow;
    for (int i=mo.move_weeks.indexmin();i<=mo.move_weeks.indexmax();i++)
    {
      if (m == mo.move_weeks(i)) return i;
    }
    return 0;
  }

int is_movement_period(const movement_info& _mi,const fishery_header_record& t1)
{
  movement_info& mi=(movement_info&) _mi;
  int match=0;
  int week=4*(t1.month-1)+t1.week;
  ivector& weeks=mi.move_weeks;
  int mmin=weeks.indexmin();
  int mmax=weeks.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    if (week == weeks(i))
    {
      match=1;
      break;
    }
  }
  return match;
}

  int operator < (const movement_info& mi,const fishery_header_record& t1)
  {
    int flag=0;
    if (mi.year<t1.year)
      flag=1;
    else if (mi.year==t1.year)
      if (mi.move_weeks(mi.current) < 4*(t1.month-1)+t1.week) 
        flag=1;
    return flag;
  }

  int operator <= (const movement_info& mi,const fishery_header_record& t1)
  {
    int flag=0;
    if (mi.year<t1.year)
      flag=1;
    else if (mi.year==t1.year)
      if (mi.move_weeks(mi.current) <= 4*(t1.month-1)+t1.week) 
        flag=1;
    return flag;
  }

  int operator == (const movement_info& _mi,const fishery_header_record& t1)
  {
    movement_info& mi=(movement_info&) _mi;
    int flag=0;
    if ( mi.year==t1.year && 
      mi.move_weeks(mi.current) == 4*(t1.month-1)+t1.week ) 
        flag=1;
    return flag;
  }

  int check_for_missing_node(const fishery_header_record& _t1,
     const fishery_header_record& _t2,const movement_info& _mi,int num_regions)
  {
    fishery_header_record& t1=(fishery_header_record&) _t1;
    fishery_header_record& t2=(fishery_header_record&) _t2;
    movement_info& mi=(movement_info&) _mi;
    int nnodes=0;
    if (mi==t1) ++mi; 
    if (t1<t2)
    {
      if (mi<t2)
      {                        // movement
        while (mi<t2)
        {
          nnodes++;
          ++mi;
        }          
      }
      if (mi==t2) // we are at a movement node but there is
      {
        ++mi;    // a fishing incident
      }
    }
    if (num_regions >1)
      return nnodes;
    else
      return 0;
  }

  int check_for_missing_node(const fishery_header_record& _t1,
     const fishery_header_record& _t2,const movement_info& _mi,
     int ip, ivector& year,ivector& month,ivector& week,
     const ivector& _move_flags,int num_regions)
  {
    ivector& move_flags= (ivector&) _move_flags;

    fishery_header_record& t1=(fishery_header_record&) _t1;
    fishery_header_record& t2=(fishery_header_record&) _t2;
    movement_info& mi=(movement_info&) _mi;
    int nnodes=0;
    if (mi<t1)
    {                        // movement
      while (mi<t1)
      {
        nnodes++;
        year(ip+nnodes)=mi.year;
        month(ip+nnodes)=mi.month();
        week(ip+nnodes)=mi.week();
        move_flags(ip+nnodes)=1;
        ++mi;
      }          
    }
    if (mi==t1) 
    {
      ++mi; 
    }
    if (t1<t2)
    {
      if (mi<t2)
      {                        // movement
        while (mi<t2)
        {
          nnodes++;
          year(ip+nnodes)=mi.year;
          month(ip+nnodes)=mi.month();
          week(ip+nnodes)=mi.week();
          move_flags(ip+nnodes)=1;
          ++mi;
        }          
      }
      if (mi==t2) // we are at a movement node but there is
      {
        move_flags(ip+1)=3;
        ++mi;    // a fishing incident
      }
    }
    if (num_regions >1)
      return nnodes;
    else
      return 0;
  }

  int add_end_nodes(const fishery_header_record& _max_record,
     const  movement_info&  _mi,int num_regions)
  {
    fishery_header_record& max_record=(fishery_header_record&) _max_record;
    movement_info& mi=(movement_info&) _mi;

    int nnodes=0;
    {                        // movement
      while (mi<=max_record)
      {
        //if (!(mi<=max_record))
        //  break;
        nnodes++;
        ++mi;
      }          
    }
    if (num_regions >1)
      return nnodes;
    else
      return 0;
  }

  int add_start_nodes(const fishery_header_record& _max_record,
     const  movement_info&  _mi,int num_regions)
  {
    fishery_header_record& max_record=(fishery_header_record&) _max_record;
    movement_info& mi=(movement_info&) _mi;

    int nnodes=0;
    {                        // movement
      while (mi<max_record)
      {
        nnodes++;
        ++mi;
      }          
    }
    if (num_regions >1)
      return nnodes;
    else
      return 0;
  }

  
  int add_first_node(fishery_header_record * min_record,
    const fishery_header_record& _current_record,
     const  movement_info&  _mi,int num_regions)
  {
    fishery_header_record& current_record=
      (fishery_header_record&) _current_record;

    int nnodes=0;
    {                        // movement
      if ((*min_record)<current_record)
      {
        nnodes++;
      }          
    }
    if (num_regions >1)
      return nnodes;
    else
      return 0;
  }

  int add_initial_nodes(const fishery_header_record& _max_record,
     const  movement_info&  _mi,int num_regions)
  {
    fishery_header_record& max_record=(fishery_header_record&) _max_record;
    movement_info& mi=(movement_info&) _mi;

    int nnodes=0;
    {                        // movement
      while (mi<=max_record)
      {
        nnodes++;
        ++mi;
      }          
    }
    if (num_regions >1)
      return nnodes;
    else
      return 0;
  }

  int add_end_nodes(const fishery_header_record& _max_record,
     const  movement_info&  _mi,
     int ip, ivector& year,ivector& month,ivector& week,
     const ivector& _move_flags,int num_regions)
  {
    ivector& move_flags= (ivector&) _move_flags;
    fishery_header_record& max_record=(fishery_header_record&) _max_record;
    movement_info& mi=(movement_info&) _mi;

    int nnodes=0;
    {                        // movement
      while (mi<=max_record)
      {
        //if (!(mi<=max_record))
        //  break;
        nnodes++;
        year(ip+nnodes)=mi.year;
        month(ip+nnodes)=mi.month();
        week(ip+nnodes)=mi.week();
        move_flags(ip+nnodes)=1;
        ++mi;
      }          
    }
    if (num_regions >1)
      return nnodes;
    else
      return 0;
  }

  int add_start_nodes(const fishery_header_record& _max_record,
     const  movement_info&  _mi,
     int ip, ivector& year,ivector& month,ivector& week,
     const ivector& _move_flags,int num_regions)
  {
    ivector& move_flags= (ivector&) _move_flags;
    fishery_header_record& max_record=(fishery_header_record&) _max_record;
    movement_info& mi=(movement_info&) _mi;

    int nnodes=0;
    {                        // movement
      while (mi<max_record)
      {
        year(ip+nnodes)=mi.year;
        month(ip+nnodes)=mi.month();
        week(ip+nnodes)=mi.week();
        move_flags(ip+nnodes)=1;
        nnodes++;
        ++mi;
      }          
    }
    if (num_regions >1)
      return nnodes;
    else
      return 0;
  }


  int add_first_node(fishery_header_record * mr,
     const fishery_header_record& _max_record,
     const  movement_info&  _mi,
     int ip, ivector& year,ivector& month,ivector& week,
     const ivector& _move_flags,int num_regions)
  {
    ivector& move_flags= (ivector&) _move_flags;
    fishery_header_record& current_record=(fishery_header_record&) _max_record;
    movement_info& mi=(movement_info&) _mi;

    int nnodes=0;
    {                        // movement
      if ((*mr)<current_record)
      {
        year(ip+nnodes)=mr->year;
        month(ip+nnodes)=mr->month;
        week(ip+nnodes)=mr->week;
        // How do we check if this is a movment period.
        if (mi==*mr)
        {
          move_flags(ip+nnodes)=1;
          nnodes++;
          ++mi;
        }
      }          
    }
    if (num_regions >1)
      return nnodes;
    else
      return 0;
  }


  void donothing(void){;}

dvar_fish_stock_history fishery_header_record_array::get_history_data(int ntg,
  int num_regions,int nage,ivector& _parest_flags,ivector& regmin,
  ivector& regmax,ivector& dataswitch,imatrix& Dflags,int month_1,
  int _mfactor,int& _first_time,ivector& mw,int _direction_flag,
  imatrix& _season_region_flags,pmulti_species_data & pmsd)
{
  int i;
  //check that the data are sorted
  ivector fishing_period(1,size());
  ivector fishing_region(1,size());
  ivector fishing_incident(1,size());

  // use num_real_regions for num_regions where necessary
  int num_real_regions=num_regions;
 
  if (pmsd)
  {
    num_real_regions=pmsd->num_real_regions;
  }
 


 // !!!!!!!!! temporary movement info
  if (!allocated(mw))
  {
    mw.allocate(1,1);
    mw(1)=1;
  }
  // renumber for month_1
  int nmove=mw.indexmax();
  if (month_1 !=1)
  {
    for (i=1;i<=nmove;i++)
    {
      int wk=mw(i)%4;
      int mn=(mw(i)-wk)/4+1;
      if (mn>= month_1)
        mw(i)-=4*(month_1-1);
      else
        mw(i)+=4*(12-month_1+1);
    }
  }    
  mw=sort(mw);
  movement_info mo; 
  mo.allocate(mw,1); 
  ivector header_movement_index(indexmin(),indexmax());
  header_movement_index.initialize();
  for (i=indexmin();i<=indexmax();i++)
  {
   header_movement_index(i)=get_movement_index(elem(i).month,elem(i).week,mo); 
  }

  int nfsh=1;
  for (i=indexmin();i<=indexmax();i++)
  {
    nfsh=max(nfsh,elem(i).fishery);
  }

  ivector nft(1,nfsh);  // the number of times a fishery occurred
  nft.initialize();

  i=indexmin();
  nft(elem(i).fishery)+=1;
  for (i=indexmin()+1;i<=indexmax();i++)
  {
    if (elem(i)>elem(i-1) || elem(i).fishery != elem(i-1).fishery)
    {
      int fi2=elem(i).fishery;
      nft(fi2)+=1;
    }
    else
    {
      //cout << "Duplicate header" << endl;
      if (!pmsd)
      {
        cerr << "This can't happen" << endl;
        ad_exit(1);
      }
    }
  }

  if (pmsd)
  {
    imatrix nfst(1,nfsh,1,pmsd->num_species);
    nfst.initialize();
    for (i=indexmin();i<=indexmax();i++)
    {
      nfst(elem(i).fishery)+=elem(i).species;
    }
    cout << nfst << " " << sum(nfst) << endl;
  }

  //int nyrs=elem(indexmax()).year-elem(indexmin()).year+1;

  int ir;
  ivector reg_nyrs_min(1,num_regions);
  ivector reg_nyrs_max(1,num_regions);
  int year_counter;
  // the records are sorted by region then by time --
  // regmin and regmax are ivectors for the min and max header
  // for each region
//***********************************************************************
//***********************************************************************
//***********************************************************************
  ivector nfp(1,num_regions);
  ivector nfp_start(1,num_regions);
  ivector nfp_end(1,num_regions);
  fishery_header_record * min_record = &elem(regmin(1));
  
  for (ir=2;ir<=num_regions;ir++)
  {
    // this should be time for earliest record so that other regions i
    // can be padded
    if (*min_record > elem(regmin(ir)))
      min_record =&elem(regmin(ir));
  }
//  imatrix nmiss1(1,num_regions,1,10000);
  imatrix nmiss1(1,num_regions,1,20000);   //NMD_27Nov2018
//  imatrix nmiss2(1,num_regions,1,10000);
  imatrix nmiss2(1,num_regions,1,20000);   //NMD_27Nov2018
  nmiss1.initialize();
  nmiss2.initialize();
//  imatrix xnmiss1(1,num_regions,1,10000);
//  imatrix xnmiss2(1,num_regions,1,10000);
  imatrix xnmiss1(1,num_regions,1,20000);   //NMD_27Nov2018
  imatrix xnmiss2(1,num_regions,1,20000);   //NMD_27Nov2018
  xnmiss1.initialize();
  xnmiss2.initialize();
  //for (ir=1;ir<=num_real_regions;ir++)
  for (ir=1;ir<=num_regions;ir++)
  {
    nfp(ir)=1;
    nfp_start(ir)=0;
    mo.initialize(*min_record);
    //nfp_start(ir)+=add_first_node(min_record,elem(regmin(ir)),mo,num_regions);
    nfp_start(ir)+=add_start_nodes(elem(regmin(ir)),mo,num_regions);
    nfp(ir)+=nfp_start(ir);
    mo.initialize(elem(regmin(ir)));
    for (i=regmin(ir);i<=regmax(ir)-1;i++)
    {
      fishery_header_record& t1=elem(i);
      fishery_header_record& t2=elem(i+1);
      if (elem(i)<elem(i+1))
      {
/*
	cout << " i=  " << i << endl;
	cout << " elem(i)=  " << elem(i) << endl;
	cout << " elem(i+1)=  " << elem(i+1) << endl;
	cout << " t1=  " << t1 << endl;
	cout << " t2=  " << t1 << endl;
*/
        //int nmiss=check_for_missing_node(t1,t2,mo);
        nmiss1(ir,i)=check_for_missing_node(t1,t2,mo,num_regions);
        xnmiss1(ir,i)=mo.current;
        if (nmiss1(ir,i))
          nfp(ir)+=nmiss1(ir,i);
        else
          donothing();
        // New fishing period
        nfp(ir)+=1;
      }
    }
  }
  fishery_header_record * max_record = &elem(regmax(1));
  //for (ir=2;ir<=num_real_regions;ir++)
  for (ir=2;ir<=num_regions;ir++)
  {
    if (*max_record < elem(regmax(ir)))
      max_record =&elem(regmax(ir));
  }
  //for (ir=1;ir<=num_real_regions;ir++)
  for (ir=1;ir<=num_regions;ir++)
  {
    mo.initialize_terminal(elem(regmax(ir)));
    nfp_end(ir)=add_end_nodes(*max_record,mo,num_regions);
    nfp(ir)+=nfp_end(ir);
  }
  //fishery_header_record& s1=elem(regmin(1));
  //fishery_header_record& s2=elem(regmin(2));
  //fishery_header_record& t1=elem(regmax(1));
  //fishery_header_record& t2=elem(regmax(2));

  int nr=num_regions;
  ivector knfp;
  imatrix nfi(1,num_regions,1,nfp);
  imatrix test_years(1,num_regions,1,nfp);
  imatrix test_months(1,num_regions,1,nfp);
  imatrix test_weeks(1,num_regions,1,nfp);
  imatrix test_move_flags(1,num_regions,1,nfp);
  imatrix test_move_index(1,num_regions,1,nfp);
  test_years.initialize();
  test_months.initialize();
  test_weeks.initialize();
  test_move_flags.initialize();
  nfi.initialize();
  for (ir=1;ir<=num_regions;ir++)
  {
    nfp(ir)=1;
    mo.initialize(*min_record);
    //nfp(ir)+=add_first_node(min_record,elem(regmin(ir)),mo,nfp(ir),
    //  test_years(ir),test_months(ir),test_weeks(ir),
    //  test_move_flags(ir),num_regions);
    nfp(ir)+=add_start_nodes(elem(regmin(ir)),mo,nfp(ir),
      test_years(ir),test_months(ir),test_weeks(ir),
      test_move_flags(ir),num_regions);
    nfi(ir,nfp(ir))=1;
    mo.initialize(elem(regmin(ir)));
   if (is_movement_period(mo,elem(regmin(ir))))
      test_move_flags(ir,nfp(ir))=4;
    for (i=regmin(ir);i<=regmax(ir)-1;i++)
    {
      if (elem(i)<elem(i+1))
      {
        //nfp(ir)+=check_for_missing_node(elem(i),elem(i+1),mo,nfp(ir),
        nmiss2(ir,i)= check_for_missing_node(elem(i),elem(i+1),mo,nfp(ir),
          test_years(ir),test_months(ir),test_weeks(ir),
          test_move_flags(ir),num_regions);
        xnmiss2(ir,i)=mo.current;
        if (nmiss2(ir,i))
          nfp(ir)+=nmiss2(ir,i);
        if (nmiss1(ir,i) != nmiss2(ir,i))
        {
          cout << "ir = " << ir << " i = " << i 
               << "  nmiss1 = " << nmiss1(ir,i)
               << "  nmiss2 = " << nmiss2(ir,i)
               << endl;
        }
        if (xnmiss1(ir,i) != xnmiss2(ir,i))
        {
          cout << "ir = " << ir << " i = " << i 
               << "  current1 = " << xnmiss1(ir,i)
               << "  current2 = " << xnmiss2(ir,i)
               << endl;
        }
        // New fishing period
        nfp(ir)+=1;
        // check if elem(i+1) is a movement period
        if (is_movement_period(mo,elem(i+1)))
          test_move_flags(ir,nfp(ir))=2;
      }
      nfi(ir,nfp(ir))+=1;
    }
  }
  for (ir=1;ir<=num_regions;ir++)
  {
    mo.initialize_terminal(elem(regmax(ir)));
    nfp(ir)+=add_end_nodes(*max_record,mo,nfp(ir),test_years(ir),
      test_months(ir),test_weeks(ir),test_move_flags(ir),num_regions);
  }
  //for (ir=1;ir<=num_regions;ir++)
  //{
  //  cout <<"rnaux2.cpp " << sum(nfi(ir)) << endl;
  //}


  int ii;
 /*
  imatrix temp(1,num_regions,1,3000);
  temp.initialize();
  
  for (ir=1;ir<=num_regions;ir++)
  {
    for (ii=1;ii<=nfp_start(ir);ii++)
    {
      year_counter = test_years(ir,ii);
      temp(ir,year_counter)+=1;
    }
    year_counter=elem(regmin(ir)).year;
    temp(ir,year_counter)+=1;

    for (i=regmin(ir);i<=regmax(ir)-1;i++)
    {
      if (elem(i).year<elem(i+1).year)
      {
        // new year
        year_counter+=elem(i+1).year-elem(i).year;
        // since there is at least one fishery in this year !!!
        // what if there isn't ?????
        if (temp(ir).indexmax()< year_counter)
        {
          ad_exit(1);
        }
        temp(ir,year_counter)+=1;
      }
      else if (elem(i).month < elem(i+1).month)
      {
        // So there is another fishing period in this year
        if (temp(ir).indexmax()< year_counter)
        {
          ad_exit(1);
        }
        temp(ir,year_counter)+=1;
      }
    }
    for (ii=0;ii<nfp_end(ir);ii++)
    {
      year_counter = test_years(ir,nfp(ir)-ii);
      temp(ir,year_counter)+=1;
    }
  }
  */
  // get the number of fishing periods for each region

  ivector nfp1(1,num_regions);
  nfp1=1;

  fishing_period(1)=1;
  fishing_region(1)=1;
  for (ir=1;ir<=num_regions;ir++)
  {
    nfp1(ir)=nfp_start(ir);
    fishing_period(ir)=nfp1(ir);
    fishing_region(ir)=ir;
    mo.initialize(elem(regmin(ir)));
    for (i=regmin(ir);i<=regmax(ir)-1;i++)
    {
      if (elem(i)<elem(i+1))
      {
        nfp1(ir)+=check_for_missing_node(elem(i),elem(i+1),mo,num_regions);
        // New fishing period
        nfp1(ir)+=1;
      }
      fishing_period(i+1)=nfp1(ir);
      fishing_region(i+1)=ir;
    }
  }
  imatrix nfi1(1,num_regions,1,nfp);
  nfi1.initialize();
  nfp1=1+nfp_start;;
  fishing_incident(1)=1;
  for (ir=1;ir<=num_regions;ir++)
  {
    nfi1(ir,nfp1(ir))=1;
    mo.initialize(elem(regmin(ir)));
    for (i=regmin(ir);i<=regmax(ir)-1;i++)
    {
      if ((elem(i)<elem(i+1)))
      {
        nfp1(ir)+=check_for_missing_node(elem(i),elem(i+1),mo,num_regions);
        // New fishing period
        nfp1(ir)+=1;
      }
      // one more fishing incident during this fishing period
      nfi1(ir,nfp1(ir))+=1;
      //nfi1(ir,nfp(ir))+=1;
      fishing_incident(i+1)=nfi1(ir,nfp1(ir));
      //fishing_incident(i+1)=nfi1(ir,nfp(ir));
    }
  }
  ivector fl1(1,200);  // These are the PAREST flags
  int ierr=0;
  for (int uu=nft.indexmin();uu<=nft.indexmax();uu++)
  {
    if ( nft(uu)==0) 
    {
      cerr << " No fishing incidents for fishery " << uu << endl
       << "this is not supported at present" << endl;
      ierr=1;
    }
  }
  if (ierr==1) 
  {
     cerr << "stopping because of this problem" << endl;
     exit(1);
  }  
     
  test_move_index.initialize();
  if (_mfactor!=0 && _mfactor!=1)
  {
    month_doubling_kludge(*this,_mfactor,regmin,_first_time,pmsd);
    int mmin=test_years.indexmin();
    int mmax=test_years.indexmax();
    for (int ir=mmin;ir<=mmax;ir++)
    {
      int mmin1=test_years(ir).indexmin();
      int mmax1=test_years(ir).indexmax();
      for (int ip=mmin1;ip<=mmax1;ip++) 
      {
        int oy=test_years(ir,ip);
        int om=test_months(ir,ip);
        int ow=test_weeks(ir,ip);
        test_move_index(ir,ip)=get_movement_index(om,ow,mo);
        if (oy)
        {
          date_struc newdate=get_new_time(oy,om,ow,_mfactor,_first_time);
          test_years(ir,ip)=newdate.year;
          test_months(ir,ip)=newdate.month;
          test_weeks(ir,ip)=newdate.week;
        }
      }
    }
  }
  else
  {
    int mmin=test_years.indexmin();
    int mmax=test_years.indexmax();
    for (int ir=mmin;ir<=mmax;ir++)
    {
      int mmin1=test_years(ir).indexmin();
      int mmax1=test_years(ir).indexmax();
      for (int ip=mmin1;ip<=mmax1;ip++) 
      {
        int oy=test_years(ir,ip);
        int om=test_months(ir,ip);
        int ow=test_weeks(ir,ip);
        test_move_index(ir,ip)=get_movement_index(om,ow,mo);
      }
    }
  }
  imatrix temp(1,num_regions,1,3000);
  temp.initialize();
  
  for (ir=1;ir<=num_regions;ir++)
  {
    for (ii=1;ii<=nfp_start(ir);ii++)
    {
      year_counter = test_years(ir,ii);
      temp(ir,year_counter)+=1;
    }
    year_counter=elem(regmin(ir)).year;
    temp(ir,year_counter)+=1;

    for (i=regmin(ir);i<=regmax(ir)-1;i++)
    {
      if (elem(i).year<elem(i+1).year)
      {
        // new year
        year_counter+=elem(i+1).year-elem(i).year;
        // since there is at least one fishery in this year !!!
        // what if there isn't ?????
        if (temp(ir).indexmax()< year_counter)
        {
          cerr << "this better not happen" << endl;
          ad_exit(1);
        }
        temp(ir,year_counter)+=1;
      }
      else if (elem(i).month < elem(i+1).month)
      {
        // So there is another fishing period in this year
        if (temp(ir).indexmax()< year_counter)
        {
          ad_exit(1);
        }
        temp(ir,year_counter)+=1;
      }
    }
    for (ii=0;ii<nfp_end(ir);ii++)
    {
      year_counter = test_years(ir,nfp(ir)-ii);
      temp(ir,year_counter)+=1;
    }
  }
  int min_year=elem(regmin(1)).year;
  int max_year=elem(regmax(1)).year;
  cout << regmax(1) << endl;
  for (ir=2;ir<=num_regions;ir++)
  {
    if (min_year>elem(regmin(ir)).year)
      min_year=elem(regmin(ir)).year;
    if (max_year<elem(regmax(ir)).year)
      max_year=elem(regmax(ir)).year;
  }
  int nyrs=max_year-min_year+1;
  ivector nage_by_region(1,num_regions);
  ivector nage_by_fishery(1,nfsh);
  if (!pmsd)
  {
    nage_by_region=nage;
    nage_by_fishery=nage;
  }
  else
  {
    if (nfsh / pmsd->num_real_fisheries != pmsd->num_species) 
    {
      cerr << "num_fisheries error" << endl;
      ad_exit(1);
    }
    int ii=0;
    for (int is=1;is<=pmsd->num_species;is++) 
    {
      for (int i=1;i<=pmsd->num_real_fisheries;i++) 
      {
        ii++;
        if (is==1)
          nage_by_fishery(ii)=nage;
        else
          nage_by_fishery(ii)=pmsd->nage(is);
      }
    }
    for (int ir=1;ir<=num_regions;ir++)
    {
      int is=pmsd->region_species_pointer(ir);
      if (is==1)
        nage_by_region(ir)=nage;
      else
        nage_by_region(ir)=pmsd->nage(is);
    }
  }
  dvar_fish_stock_history fsh(ntg,num_regions,nage,nfp,nfi,nfsh,nyrs,
    _parest_flags,fl1,nft,regmin,regmax,dataswitch,Dflags,mo,_direction_flag,
    _mfactor,_season_region_flags,pmsd,nage_by_region,nage_by_fishery);

//NMD 03Nov2011
   if (pmsd)
   {
     pmsd->nage_by_region=nage_by_region;
     int tmp;
     tmp = max(pmsd->nage);
     if (allocated(pmsd->length_month_yr)) pmsd->length_month_yr.deallocate();
     pmsd->length_month_yr.allocate(1,pmsd->num_species,1,_mfactor,1,tmp);
     if (allocated(pmsd->sdevs_yr)) pmsd->sdevs_yr.deallocate();
     pmsd->sdevs_yr.allocate(1,pmsd->num_species,1,nyrs,1,tmp);
     if (allocated(pmsd->nat_mort)) pmsd->nat_mort.deallocate();
     pmsd->nat_mort.allocate(2,pmsd->num_species,1,nyrs,1,tmp);
     if (allocated(pmsd->region_rec_diff_sums)) 
       pmsd->region_rec_diff_sums.deallocate();
     pmsd->region_rec_diff_sums.allocate(2,pmsd->num_species,1,nyrs);

     if (allocated(pmsd->fisc)) 
       pmsd->fisc.deallocate();
     pmsd->fisc.allocate(1,num_regions,1,nfp,1,nfi,1,pmsd->num_species);

     if (allocated(pmsd->fisc_lf)) 
       pmsd->fisc_lf.deallocate();
     pmsd->fisc_lf.allocate(1,num_regions,1,nfp,1,nfi,1,pmsd->num_species);

     if (allocated(pmsd->fisc_wf)) 
       pmsd->fisc_wf.deallocate();
     pmsd->fisc_wf.allocate(1,num_regions,1,nfp,1,nfi,1,pmsd->num_species);

     if (allocated(pmsd->fisn)) 
       pmsd->fisn.deallocate();
     pmsd->fisn.allocate(1,num_regions,1,nfp,1,nfi);

     if (allocated(pmsd->fisn_lf)) 
       pmsd->fisn_lf.deallocate();
     pmsd->fisn_lf.allocate(1,num_regions,1,nfp,1,nfi);

     if (allocated(pmsd->fisn_wf)) 
       pmsd->fisn_wf.deallocate();
     pmsd->fisn_wf.allocate(1,num_regions,1,nfp,1,nfi);

     if (allocated(pmsd->rec_delta)) 
       pmsd->rec_delta.deallocate(); 
     pmsd->rec_delta.allocate(2,pmsd->num_species,1,nyrs); 

     if (allocated(pmsd->sp_in_catch))
       pmsd->sp_in_catch.deallocate();
     pmsd->sp_in_catch.
       allocate(1,num_regions,1,nfp,1,nfi,1,pmsd->num_species);

     if (allocated(pmsd->reg_in_catch))
       pmsd->reg_in_catch.deallocate();
     pmsd->reg_in_catch.
       allocate(1,num_regions,1,nfp,1,nfi,1,pmsd->num_species);

     if (allocated(pmsd->sp_in_lf))
       pmsd->sp_in_lf.deallocate();
     pmsd->sp_in_lf.
       allocate(1,num_regions,1,nfp,1,nfi,1,pmsd->num_species);

     if (allocated(pmsd->reg_in_lf))
       pmsd->reg_in_lf.deallocate();
     pmsd->reg_in_lf.
       allocate(1,num_regions,1,nfp,1,nfi,1,pmsd->num_species);

     if (allocated(pmsd->sp_in_wf))
       pmsd->sp_in_wf.deallocate();
     pmsd->sp_in_wf.
       allocate(1,num_regions,1,nfp,1,nfi,1,pmsd->num_species);

     if (allocated(pmsd->reg_in_wf))
       pmsd->reg_in_wf.deallocate();
     pmsd->reg_in_wf.
       allocate(1,num_regions,1,nfp,1,nfi,1,pmsd->num_species);




     //NMD 03Nov2011
   } 

  fsh.age_nage=fishery_freq_record::age_nage;
  fsh.age_age1=fishery_freq_record::age_age1;

  breakup_big_function(nfp,nft,fsh,nfi,*this,temp,test_years,
    test_months,test_weeks,test_move_index,header_movement_index);

  if (!allocated(fsh.move_flags)) fsh.move_flags.allocate(test_move_flags);
  fsh.move_flags=test_move_flags;
  //fsh.do_recruitment_period_calculations();
  //fsh.recr.allocate(1,fsh.num_recruitment_periods);
  //  fsh.recr.allocate(1,fsh.nyears);
  fsh.recr.allocate(2,fsh.nyears);  //NMD_19May2016
  if (fsh.pmsd && fsh.pmsd->num_species>1)
  {
    for (int is=2;is<=fsh.pmsd->num_species;is++)
    {
      if (!allocated(fsh.pmsd->recr(is)))
//        fsh.pmsd->recr(is).allocate(1,fsh.nyears);
        fsh.pmsd->recr(is).allocate(2,fsh.nyears);  //NMD_19May2016
    }
  }
   
  // Need to modify this for different numbers of age classes
  // in each species  // DF 21 Aug 2018
  //fsh.N.allocate(1,num_regions,1,fsh.nyears,1,nage_by_region);
  //fsh.Nsave.allocate(1,num_regions,1,fsh.nyears,1,nage_by_region);
  //fsh.exp_N.allocate(1,num_regions,1,fsh.nyears,1,nage_by_region);

  for (int i=1;i<=num_regions;i++)
  {
    fsh.N(i)=1.e+225;
  }
  fsh.Rsave.allocate(1,num_regions,1,fsh.nyears);
  fsh.num_fish_data_recs=size();
  if (fsh.fishing_period.indexmax() != fishing_period.indexmax())
  {
    fsh.fishing_period.deallocate();
  }
  fsh.fishing_period=fishing_period;

  if (fsh.fishing_incident.indexmax() != fishing_incident.indexmax())
  {
    fsh.fishing_incident.deallocate();
  }
  fsh.fishing_incident=fishing_incident;
  if (fsh.num_tag_releases)
  {
    fsh.tag_flags.allocate(1,fsh.num_tag_releases,1,10);
    fsh.true_tag_flags.allocate(1,fsh.num_tag_releases,1,10);
  }
  fsh.true_num_tag_releases=fsh.num_tag_releases;

  //fsh.make_fishing_period_report();
  fsh.get_fishery_realization_index();

  return fsh;
}

void dvar_fish_stock_history::
  do_recruitment_period_calculations(void)
{
  if(!age_flags(93))
  {
    age_flags(93)=1;
    rec_times(1)=1;
  }
    
  int minregion=1;
  int mintime=48*(year(1,1)-1)+12.*(month(1,1)-1)
    +4.*(week(1,1)-1)+1;
  int ir;
  for (ir=2;ir<=num_regions;ir++)
  {
    int tm1=48.*(year(ir,1)-1)+12.*(month(ir,1)-1)+
      4.*(week(ir,1)-1)+1;
    if (tm1 < mintime)
    {
      mintime=tm1;
      minregion=ir;
    }
  }

  int i=0;
  for (i=1;i<=age_flags(93);i++)
  {
    if (month(minregion,1)<rec_times(i))
      break;
  }
  int irc=age_flags(93)+2-i;
  initial_recruitment_count=irc;
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++) 
    {
      int i=0;
      for (i=1;i<=age_flags(93);i++)
      {
        if (month(ir,ip)< rec_times(i))
           break;
      }
      recruitment_period(ir,ip)=irc+age_flags(93)*(year(ir,ip)-2.)
        +i-1;
    }
  }
  num_recruitment_periods=max(recruitment_period(1));
  int tm=0;
  for (ir=2;ir<=num_regions;ir++)
  {
    tm=max(recruitment_period(ir));
    if (num_recruitment_periods<tm) num_recruitment_periods=tm;
  }
}

void dvar_fish_stock_history::get_initial_tag_recruitment_period(void)
{
  int irc=initial_recruitment_count;
  for (int it=1;it<=num_tag_releases;it++)
  {  
    for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
    {
      int ip=initial_tag_period(it,ir);
      int i=0;
      for (i=1;i<=age_flags(93);i++)
      {
        if (month(ir,ip)< 1000)
           break;
        if (rec_times(i)<1000)
           break;
        if (month(ir,ip)< rec_times(i))
           break;
      }
      initial_tag_recruitment_period(it,ir)=irc+age_flags(93)*(year(ir,ip)-2.)
           +i-1;
    }
  }
}

void dvar_fish_stock_history::get_fishery_realization_index(void)
{
  fishery_realization_index.initialize();
  for (int fi=1;fi<=num_fisheries;fi++)
  {
    for (int nt=1;nt<=num_fish_times(fi);nt++)
    {
      int rr=realization_region(fi,nt);
      int rp=realization_period(fi,nt);
      int badflag=1;
      for (int i=1;i<=num_fish_incidents(rr,rp);i++)
      {
        if (parent(rr,rp,i)==fi)
        {
          fishery_realization_index(rr,rp,i)=nt;
          badflag=0;
          break;
        }
      }
      if (badflag==1)
      {
        cerr << "This can't happen" << endl;
        ad_exit(1);
      }
    }
  }
}

void do_between_time_calculations(dvar_fish_stock_history& fsh)
{
  MY_DOUBLE_TYPE weekdays=7.0;
  MY_DOUBLE_TYPE monthdays=7.0*4.0;
  MY_DOUBLE_TYPE yeardays=12.0*7.0*4.0;
  int tmult=1;
  if (fsh.age_flags(57)) tmult= fsh.age_flags(57);
  for (int i=1;i<=fsh.num_fisheries;i++)
  {
    if (fsh.fish_flags(i,23)==0)
    {
      for (int nt=2;nt<=fsh.num_fish_times(i);nt++)
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
        // check for time errors DF Feb22 05
        if (fsh.between_times(i,nt)<=0)
        {
          cerr << "Error in betwee times calculations generated by"
            << endl << "fishery data for fishery " << i << " in time "
            << nt << endl 
            << " whihc is year pair " << fsh.year(rr1,rp1)
            << "  " << fsh.year(rr2,rp2)
            << endl << " month pair " << fsh.month(rr1,rp1)
            << "  " << fsh.month(rr2,rp2)
            << endl << " month pair " << fsh.week(rr1,rp1)
            << "  " << fsh.week(rr2,rp2) << endl;
          ad_exit(1);
        }

        // inverse of the time period for use in penalty function
        // results are in 1/yr
        fsh.between_times(i,nt)=(tmult*yeardays)/fsh.between_times(i,nt);
      }
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


        int year=fsh.data_fish_flags(2,i);
        int month=fsh.data_fish_flags(3,i);
        int ty=fsh.really_true_year(rr2,rp2)+fsh.year1-1;
        int tm=fsh.really_true_month(rr2,rp2);

        if (tmp>t_interval && ( (year==0) || (ty<year) || 
           (ty==year && tm < month) ) )
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
}

void do_robust_mean_spread_calcs(dvar_fish_stock_history& fsh)
{
  MY_DOUBLE_TYPE weekdays=7.0;
  MY_DOUBLE_TYPE monthdays=7.0*4.0;
  MY_DOUBLE_TYPE yeardays=12.0*7.0*4.0;
  int tmult=1;
  if (fsh.age_flags(57)) tmult= fsh.age_flags(57);
  for (int i=1;i<=fsh.num_fisheries;i++)
  {
    MY_DOUBLE_TYPE tmp=0.0;
    MY_DOUBLE_TYPE t_interval = fsh.fish_flags(i,23)*tmult*monthdays;
    for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
    {
      if (nt==1 || fsh.between_times(i,nt))
      {
        tmp=0.0;
        if (fsh.fish_flags(i,25))
        {
          int ii;
          for (ii=nt;ii>=1;ii--)
          {
            int rp2=fsh.realization_period(i,nt);
            int rp1=fsh.realization_period(i,ii);
            int rr2=fsh.realization_region(i,nt);
            int rr1=fsh.realization_region(i,ii);
            tmp=yeardays*(fsh.year(rr2,rp2)-fsh.year(rr1,rp1));
            tmp+=monthdays*(fsh.month(rr2,rp2)-fsh.month(rr1,rp1));
            tmp+=weekdays*(fsh.week(rr2,rp2)-fsh.week(rr1,rp1));
	    if (tmp/monthdays >fsh.fish_flags(i,25)) break;
          }
	  if (ii<1) ii=1;
          fsh.imp_back(i,nt)=ii;
	  tmp=0.0;
          for (ii=nt;ii<=fsh.num_fish_times(i);ii++)
	  {
            int rp1=fsh.realization_period(i,nt);
            int rp2=fsh.realization_period(i,ii);
            int rr1=fsh.realization_region(i,nt);
            int rr2=fsh.realization_region(i,ii);
            tmp=yeardays*(fsh.year(rr2,rp2)-fsh.year(rr1,rp1));
            tmp+=monthdays*(fsh.month(rr2,rp2)-fsh.month(rr1,rp1));
            tmp+=weekdays*(fsh.week(rr2,rp2)-fsh.week(rr1,rp1));
	    if (tmp/monthdays >fsh.fish_flags(i,25)) break;
          }
	  if (ii>fsh.num_fish_times(i)) ii=fsh.num_fish_times(i);
          fsh.imp_forward(i,nt)=ii;
	  tmp=0.0;
        }
      }
    }
  }
}



#if !defined(__GNUDOS__) && !defined(__MSVC32__)
int max(int x,int y)
{
  if (x>=y) return x;
  return y;
}

int min(int x,int y)
{
  if (x<=y) return x;
  return y;
}
#endif

void breakup_big_function(ivector& nfp,ivector& nft,
  dvar_fish_stock_history& fsh,
  imatrix& nfi,fishery_header_record_array& fhra,imatrix& temp,
  imatrix& test_year,imatrix& test_month,imatrix& test_week,
  imatrix& test_move_index,ivector& header_movement_index)
{
  nft.initialize();
  fsh.obs_region_tot_catch.initialize();
  int i=1;
  //ofstream ofss("testyear");
  fsh.year.initialize();
  fsh.month.initialize();
  fsh.week.initialize();
  fsh.move_index.initialize();
  imatrix temp_test(1,fsh.num_regions,1,fsh.nyears);
  temp_test.initialize();
  int ir;
  fsh.fish_times.initialize();
  for (ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)  // Loop over fishing periods
    {
      if (fsh.num_fish_incidents(ir,ip)) // Loop over fishing
      {
        fsh.year(ir,ip)=fhra.elem(i).year;
        fsh.movement_period(ir,ip)=fhra.elem(i).movement_period;
        fsh.move_index(ir,ip)=header_movement_index(i);
        fsh.month(ir,ip)=fhra.elem(i).month;
        fsh.week(ir,ip)=fhra.elem(i).week;

        if (fsh.month(ir,ip) >12 || fsh.month(ir,ip) <1 || 
           fsh.week(ir,ip) >4 || fsh.week(ir,ip) <1  )
        {
          cerr << " Error reading in fsh.month or fsh.week field from record " << i << endl;
          cerr << " fsh.month value is " << fsh.month(ir,ip) << endl;
          cerr << " fsh.week value is " << fsh.week(ir,ip) << endl;
          exit(1);
        }
        fsh.fraction(ir,ip)=-log(double(temp(ir,fsh.year(ir,ip))));  // apportion natural mortality
      }
      else
      {
        fsh.year(ir,ip)=test_year(ir,ip);
        fsh.month(ir,ip)=test_month(ir,ip);
        fsh.week(ir,ip)=test_week(ir,ip);
        fsh.move_index(ir,ip)=test_move_index(ir,ip);
      }

      if (ip>1)
      {
        if ( ( fsh.year(ir,ip-1)== fsh.year(ir,ip))
          && (fsh.month(ir,ip-1)== fsh.month(ir,ip))
          && (fsh.week(ir,ip-1)== fsh.week(ir,ip)) )
        {
          cerr << "Ideintical record times at record " << i << endl
          << " for fishing period " << ip << " in region " << ir << endl;
          ad_exit(1);
        }
      }
      
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        fsh.header_record_index(ir,ip,fi)=i;
        nft(fhra.elem(i).fishery)+=1; // count the number of realizations 
                                 // of each fishery
        int fi1=fhra.elem(i).fishery;
        int nft1=nft(fi1);

        fsh.realization_region(fi1,nft1);

        fsh.realization_region(
          fhra.elem(i).fishery,nft(fhra.elem(i).fishery))=ir;
        fsh.realization_period(
          fhra.elem(i).fishery,nft(fhra.elem(i).fishery))=ip;
        fsh.realization_incident(
          fhra.elem(i).fishery,nft(fhra.elem(i).fishery))=fi;
        double ff=fhra.elem(i).fishery;
        fsh.parent(ir,ip,fi)=ff;
        
        fsh.fish_times(ir,ip,fi)=nft(fhra.elem(i).fishery);
        fsh.effort(ir,ip,fi)=fhra.elem(i).effort;

        fsh.effort_weight(ir,ip,fi)=fhra.elem(i).effort_weight;

        fsh.obs_tot_catch(ir,ip,fi)=fhra.elem(i).obs_tot_catch;
        fsh.obs_region_tot_catch(ir,ip)+=fsh.obs_tot_catch(ir,ip,fi);
        if (fsh.pmsd)
        {
          fsh.pmsd->fisc(ir,ip,fi)=fhra.elem(i).species1;
          fsh.pmsd->fisn(ir,ip,fi)=sum(fhra.elem(i).species1);

          fsh.pmsd->fisc_lf(ir,ip,fi)=fhra.elem(i).species2;
          fsh.pmsd->fisn_lf(ir,ip,fi)=sum(fhra.elem(i).species2);

          fsh.pmsd->fisc_wf(ir,ip,fi)=fhra.elem(i).species3;
          fsh.pmsd->fisn_wf(ir,ip,fi)=sum(fhra.elem(i).species3);

          if (allocated(fsh.pmsd->sp_in_catch(ir,ip,fi)))
            fsh.pmsd->sp_in_catch(ir,ip,fi).deallocate();
          fsh.pmsd
            ->sp_in_catch(ir,ip,fi).allocate(1,fsh.pmsd->fisn(ir,ip,fi));
          if (allocated(fsh.pmsd->reg_in_catch(ir,ip,fi)))
            fsh.pmsd->reg_in_catch(ir,ip,fi).deallocate();
          fsh.pmsd->
            reg_in_catch(ir,ip,fi).allocate(1,fsh.pmsd->fisn(ir,ip,fi));

          if (allocated(fsh.pmsd->sp_in_lf(ir,ip,fi)))
            fsh.pmsd->sp_in_lf(ir,ip,fi).deallocate();
          fsh.pmsd
            ->sp_in_lf(ir,ip,fi).allocate(1,fsh.pmsd->fisn_lf(ir,ip,fi));
          if (allocated(fsh.pmsd->reg_in_lf(ir,ip,fi)))
            fsh.pmsd->reg_in_lf(ir,ip,fi).deallocate();
          fsh.pmsd->
            reg_in_lf(ir,ip,fi).allocate(1,fsh.pmsd->fisn_lf(ir,ip,fi));

          if (allocated(fsh.pmsd->sp_in_wf(ir,ip,fi)))
            fsh.pmsd->sp_in_wf(ir,ip,fi).deallocate();
          fsh.pmsd
            ->sp_in_wf(ir,ip,fi).allocate(1,fsh.pmsd->fisn_wf(ir,ip,fi));
          if (allocated(fsh.pmsd->reg_in_wf(ir,ip,fi)))
            fsh.pmsd->reg_in_wf(ir,ip,fi).deallocate();
          fsh.pmsd->
            reg_in_wf(ir,ip,fi).allocate(1,fsh.pmsd->fisn_wf(ir,ip,fi));

          int jj=1;
          int jj1=1;
          int jj2=1;
          for (int ii=1;ii<=fsh.pmsd->num_species;ii++)
          {
            if (fhra.elem(i).species1(ii)>0)
            {
              fsh.pmsd->sp_in_catch(ir,ip,fi,jj)=ii;

              int ir1=(ir-1)%fsh.pmsd->num_real_regions+1;
              fsh.pmsd->reg_in_catch(ir,ip,fi,jj)
                =ir1+(ii-1)*fsh.pmsd->num_real_regions;
              jj++;
            }
            if (fhra.elem(i).species2(ii)>0)
            {
              fsh.pmsd->sp_in_lf(ir,ip,fi,jj1)=ii;

              int ir1=(ir-1)%fsh.pmsd->num_real_regions+1;
              fsh.pmsd->reg_in_lf(ir,ip,fi,jj1)
                =ir1+(ii-1)*fsh.pmsd->num_real_regions;
              jj1++;
            }
            if (fhra.elem(i).species3(ii)>0)
            {
              fsh.pmsd->sp_in_wf(ir,ip,fi,jj2)=ii;

              int ir1=(ir-1)%fsh.pmsd->num_real_regions+1;
              fsh.pmsd->reg_in_wf(ir,ip,fi,jj2)
                =ir1+(ii-1)*fsh.pmsd->num_real_regions;
              jj2++;
            }
          }
        }

        if (fsh.effort(ir,ip,fi)>-0.5L) fsh.total_num_obs+=1;
        if (fsh.obs_tot_catch(ir,ip,fi)>-0.5L) fsh.total_num_obs+=1;

        i++;
      }
    }
  }
  for (ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)  
    {
      temp_test(ir,fsh.year(ir,ip))+=1.0;
    }
  }
  for (ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)  
    {
      fsh.fraction(ir,ip)=-log(double(temp_test(ir,fsh.year(ir,ip))));  
     // apportion natural mortality
    }
  }
  //ofss << setw(10) << fsh.year(1)(1,40) << endl << setw(10) << test_year(1)(1,40) << endl;
  //ofss << setw(10) << fsh.month(1)(1,40) << endl << setw(10) << test_month(1)(1,40) << endl;
  //ofss << setw(10) << fsh.week(1)(1,40) << endl << setw(10) << test_week(1)(1,40) << endl;
}

void set_true_months(dvar_len_fish_stock_history& fsh)
{
  int ir;
  for (ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)  // Loop over fishing periods
    {
      if (fsh.month_factor!=0 && fsh.month_factor!=1)
      {         
        date_struc newdate=get_old_time(fsh.year(ir,ip),
          fsh.month(ir,ip),fsh.week(ir,ip),fsh.month_factor,fsh.first_time);
        fsh.true_month(ir,ip)=newdate.month;
        fsh.true_year(ir,ip)=newdate.year;
        fsh.true_week(ir,ip)=newdate.week;
      }
      else
      {
        fsh.true_month(ir,ip)=fsh.month(ir,ip);
        fsh.true_year(ir,ip)=fsh.year(ir,ip);
      }
    }
  }
  fsh.un_normalize_year();
}   


void  get_seasonal_catchability_pars_index(dvar_len_fish_stock_history& fsh)
{	
  ivector ff47=column(fsh.fish_flags,47); 
   
  dvector chk=dvector(ff47);
  if (norm2(chk-mean(chk))>1.e-6)
  {
    cerr << " Option for different number of coffs in explicit seasonal catchability"
            " not yet implemented" << endl;
    exit(1);	    
  }	    

  int month1=fsh.month_1;
  int num_pars=ff47(1);
  MY_DOUBLE_TYPE delta=12.0/num_pars;
  int im;
  for (int i=1;i<=12;i++)
  {
    if ( (i-month1)<0) 
      im=i-month1+13;
    else  
      im=i-month1+1;

     MY_DOUBLE_TYPE d1=(im-1)/delta+1.0;
     int id1=int(d1);
     MY_DOUBLE_TYPE diff1=d1-id1;
      
     fsh.seasonal_catchability_pars_index(i,1)=id1;
     if (int(d1)==num_pars)
       fsh.seasonal_catchability_pars_index(i,2)=1;
     else
       fsh.seasonal_catchability_pars_index(i,2)=id1+1;
	       
    fsh.seasonal_catchability_pars_mix(i,1)=1.0-diff1;
    fsh.seasonal_catchability_pars_mix(i,2)=diff1;
     
    
  }
    cout <<"rnaux2.cpp " << fsh.seasonal_catchability_pars_index << endl;
    cout <<"rnaux2.cpp " << fsh.seasonal_catchability_pars_mix << endl;
}  

dvar_fish_stock_history fishery_header_record_array::get_history_data(int ntg,
  int num_regions,int nage,ivector& _parest_flags,ivector& regmin,
  ivector& regmax,ivector& dataswitch,imatrix& Dflags,int month_1,
  int _mfactor,int& _first_time,ivector& mw,int _direction_flag,
  imatrix& _season_region_flags)
{
  cout << "this get history data got called" << endl;
  ad_exit(1);
  int i;
  //check that the data are sorted
  ivector fishing_period(1,size());
  ivector fishing_region(1,size());
  ivector fishing_incident(1,size());


 // !!!!!!!!! temporary movement info
  if (!allocated(mw))
  {
    mw.allocate(1,1);
    mw(1)=1;
  }
  // renumber for month_1
  int nmove=mw.indexmax();
  if (month_1 !=1)
  {
    for (i=1;i<=nmove;i++)
    {
      int wk=mw(i)%4;
      int mn=(mw(i)-wk)/4+1;
      if (mn>= month_1)
        mw(i)-=4*(month_1-1);
      else
        mw(i)+=4*(12-month_1+1);
    }
  }    
  mw=sort(mw);
  movement_info mo; 
  mo.allocate(mw,1); 
  ivector header_movement_index(indexmin(),indexmax());
  header_movement_index.initialize();
  for (i=indexmin();i<=indexmax();i++)
  {
   header_movement_index(i)=get_movement_index(elem(i).month,elem(i).week,mo); 
  }

  int nfsh=1;
  for (i=indexmin();i<=indexmax();i++)
  {
    nfsh=max(nfsh,elem(i).fishery);
  }

  ivector nft(1,nfsh);  // the number of times a fishery occurred
  nft.initialize();

  for (i=indexmin();i<=indexmax();i++)
  {
    nft(elem(i).fishery)+=1;
  }

  //int nyrs=elem(indexmax()).year-elem(indexmin()).year+1;

  int ir;
  ivector reg_nyrs_min(1,num_regions);
  ivector reg_nyrs_max(1,num_regions);
  int year_counter;
  // the records are sorted by region then by time --
  // regmin and regmax are ivectors for the min and max header
  // for each region
//***********************************************************************
//***********************************************************************
//***********************************************************************
  ivector nfp(1,num_regions);
  ivector nfp_start(1,num_regions);
  ivector nfp_end(1,num_regions);
  fishery_header_record * min_record = &elem(regmin(1));
  for (ir=2;ir<=num_regions;ir++)
  {
    // this should be time for earliest record so that other regions i
    // can be padded
    if (*min_record > elem(regmin(ir)))
      min_record =&elem(regmin(ir));
  }
  imatrix nmiss1(1,num_regions,1,10000);
  imatrix nmiss2(1,num_regions,1,10000);
  nmiss1.initialize();
  nmiss2.initialize();
  imatrix xnmiss1(1,num_regions,1,10000);
  imatrix xnmiss2(1,num_regions,1,10000);
  xnmiss1.initialize();
  xnmiss2.initialize();
  for (ir=1;ir<=num_regions;ir++)
  {
    nfp(ir)=1;
    nfp_start(ir)=0;
    mo.initialize(*min_record);
    //nfp_start(ir)+=add_first_node(min_record,elem(regmin(ir)),mo,num_regions);
    nfp_start(ir)+=add_start_nodes(elem(regmin(ir)),mo,num_regions);
    nfp(ir)+=nfp_start(ir);
    mo.initialize(elem(regmin(ir)));
    for (i=regmin(ir);i<=regmax(ir)-1;i++)
    {
      fishery_header_record& t1=elem(i);
      fishery_header_record& t2=elem(i+1);
      if (elem(i)<elem(i+1))
      {
        //int nmiss=check_for_missing_node(t1,t2,mo);
        nmiss1(ir,i)=check_for_missing_node(t1,t2,mo,num_regions);
        xnmiss1(ir,i)=mo.current;
        if (nmiss1(ir,i))
          nfp(ir)+=nmiss1(ir,i);
        else
          donothing();
        // New fishing period
        nfp(ir)+=1;
      }
    }
  }
  fishery_header_record * max_record = &elem(regmax(1));
  for (ir=2;ir<=num_regions;ir++)
  {
    if (*max_record < elem(regmax(ir)))
      max_record =&elem(regmax(ir));
  }
  for (ir=1;ir<=num_regions;ir++)
  {
    mo.initialize_terminal(elem(regmax(ir)));
    nfp_end(ir)=add_end_nodes(*max_record,mo,num_regions);
    nfp(ir)+=nfp_end(ir);
  }
  //fishery_header_record& s1=elem(regmin(1));
  //fishery_header_record& s2=elem(regmin(2));
  //fishery_header_record& t1=elem(regmax(1));
  //fishery_header_record& t2=elem(regmax(2));

  imatrix nfi(1,num_regions,1,nfp);
  imatrix test_years(1,num_regions,1,nfp);
  imatrix test_months(1,num_regions,1,nfp);
  imatrix test_weeks(1,num_regions,1,nfp);
  imatrix test_move_flags(1,num_regions,1,nfp);
  imatrix test_move_index(1,num_regions,1,nfp);
  test_years.initialize();
  test_months.initialize();
  test_weeks.initialize();
  test_move_flags.initialize();
  nfi.initialize();
  for (ir=1;ir<=num_regions;ir++)
  {
    nfp(ir)=1;
    mo.initialize(*min_record);
    //nfp(ir)+=add_first_node(min_record,elem(regmin(ir)),mo,nfp(ir),
    //  test_years(ir),test_months(ir),test_weeks(ir),
    //  test_move_flags(ir),num_regions);
    nfp(ir)+=add_start_nodes(elem(regmin(ir)),mo,nfp(ir),
      test_years(ir),test_months(ir),test_weeks(ir),
      test_move_flags(ir),num_regions);
    nfi(ir,nfp(ir))=1;
    mo.initialize(elem(regmin(ir)));
   if (is_movement_period(mo,elem(regmin(ir))))
      test_move_flags(ir,nfp(ir))=4;
    for (i=regmin(ir);i<=regmax(ir)-1;i++)
    {
      if (elem(i)<elem(i+1))
      {
        //nfp(ir)+=check_for_missing_node(elem(i),elem(i+1),mo,nfp(ir),
        nmiss2(ir,i)= check_for_missing_node(elem(i),elem(i+1),mo,nfp(ir),
          test_years(ir),test_months(ir),test_weeks(ir),
          test_move_flags(ir),num_regions);
        xnmiss2(ir,i)=mo.current;
        if (nmiss2(ir,i))
          nfp(ir)+=nmiss2(ir,i);
        if (nmiss1(ir,i) != nmiss2(ir,i))
        {
          cout << "ir = " << ir << " i = " << i 
               << "  nmiss1 = " << nmiss1(ir,i)
               << "  nmiss2 = " << nmiss2(ir,i)
               << endl;
        }
        if (xnmiss1(ir,i) != xnmiss2(ir,i))
        {
          cout << "ir = " << ir << " i = " << i 
               << "  current1 = " << xnmiss1(ir,i)
               << "  current2 = " << xnmiss2(ir,i)
               << endl;
        }
        // New fishing period
        nfp(ir)+=1;
        // check if elem(i+1) is a movement period
        if (is_movement_period(mo,elem(i+1)))
          test_move_flags(ir,nfp(ir))=2;
      }
      nfi(ir,nfp(ir))+=1;
    }
  }
  for (ir=1;ir<=num_regions;ir++)
  {
    mo.initialize_terminal(elem(regmax(ir)));
    nfp(ir)+=add_end_nodes(*max_record,mo,nfp(ir),test_years(ir),
      test_months(ir),test_weeks(ir),test_move_flags(ir),num_regions);
  }
  //for (ir=1;ir<=num_regions;ir++)
  //{
  //  cout <<"rnaux2.cpp " << sum(nfi(ir)) << endl;
  //}


  int ii;
 /*
  imatrix temp(1,num_regions,1,3000);
  temp.initialize();
  
  for (ir=1;ir<=num_regions;ir++)
  {
    for (ii=1;ii<=nfp_start(ir);ii++)
    {
      year_counter = test_years(ir,ii);
      temp(ir,year_counter)+=1;
    }
    year_counter=elem(regmin(ir)).year;
    temp(ir,year_counter)+=1;

    for (i=regmin(ir);i<=regmax(ir)-1;i++)
    {
      if (elem(i).year<elem(i+1).year)
      {
        // new year
        year_counter+=elem(i+1).year-elem(i).year;
        // since there is at least one fishery in this year !!!
        // what if there isn't ?????
        if (temp(ir).indexmax()< year_counter)
        {
          ad_exit(1);
        }
        temp(ir,year_counter)+=1;
      }
      else if (elem(i).month < elem(i+1).month)
      {
        // So there is another fishing period in this year
        if (temp(ir).indexmax()< year_counter)
        {
          ad_exit(1);
        }
        temp(ir,year_counter)+=1;
      }
    }
    for (ii=0;ii<nfp_end(ir);ii++)
    {
      year_counter = test_years(ir,nfp(ir)-ii);
      temp(ir,year_counter)+=1;
    }
  }
  */
  // get the number of fishing periods for each region

  ivector nfp1(1,num_regions);
  nfp1=1;

  fishing_period(1)=1;
  fishing_region(1)=1;
  for (ir=1;ir<=num_regions;ir++)
  {
    nfp1(ir)=nfp_start(ir);
    fishing_period(ir)=nfp1(ir);
    fishing_region(ir)=ir;
    mo.initialize(elem(regmin(ir)));
    for (i=regmin(ir);i<=regmax(ir)-1;i++)
    {
      if (elem(i)<elem(i+1))
      {
        nfp1(ir)+=check_for_missing_node(elem(i),elem(i+1),mo,num_regions);
        // New fishing period
        nfp1(ir)+=1;
      }
      fishing_period(i+1)=nfp1(ir);
      fishing_region(i+1)=ir;
    }
  }
  imatrix nfi1(1,num_regions,1,nfp);
  nfi1.initialize();
  nfp1=1+nfp_start;;
  fishing_incident(1)=1;
  for (ir=1;ir<=num_regions;ir++)
  {
    nfi1(ir,nfp1(ir))=1;
    mo.initialize(elem(regmin(ir)));
    for (i=regmin(ir);i<=regmax(ir)-1;i++)
    {
      if ((elem(i)<elem(i+1)))
      {
        nfp1(ir)+=check_for_missing_node(elem(i),elem(i+1),mo,num_regions);
        // New fishing period
        nfp1(ir)+=1;
      }
      // one more fishing incident during this fishing period
      nfi1(ir,nfp1(ir))+=1;
      //nfi1(ir,nfp(ir))+=1;
      fishing_incident(i+1)=nfi1(ir,nfp1(ir));
      //fishing_incident(i+1)=nfi1(ir,nfp(ir));
    }
  }
  ivector fl1(1,200);  // These are the PAREST flags
  int ierr=0;
  for (int uu=nft.indexmin();uu<=nft.indexmax();uu++)
  {
    if ( nft(uu)==0) 
    {
      cerr << " No fishing incidents for fishery " << uu << endl
       << "this is not supported at present" << endl;
      ierr=1;
    }
  }
  if (ierr==1) 
  {
     cerr << "stopping because of this problem" << endl;
     exit(1);
  }  
     
  test_move_index.initialize();
  if (_mfactor!=0 && _mfactor!=1)
  {
    cerr << "Need to fix month odubling kludge to get pmsd" << endl;
    ad_exit(1);
   // month_doubling_kludge(*this,_mfactor,regmin,_first_time);
    int mmin=test_years.indexmin();
    int mmax=test_years.indexmax();
    for (int ir=mmin;ir<=mmax;ir++)
    {
      int mmin1=test_years(ir).indexmin();
      int mmax1=test_years(ir).indexmax();
      for (int ip=mmin1;ip<=mmax1;ip++) 
      {
        int oy=test_years(ir,ip);
        int om=test_months(ir,ip);
        int ow=test_weeks(ir,ip);
        test_move_index(ir,ip)=get_movement_index(om,ow,mo);
        if (oy)
        {
          date_struc newdate=get_new_time(oy,om,ow,_mfactor,_first_time);
          test_years(ir,ip)=newdate.year;
          test_months(ir,ip)=newdate.month;
          test_weeks(ir,ip)=newdate.week;
        }
      }
    }
  }
  else
  {
    int mmin=test_years.indexmin();
    int mmax=test_years.indexmax();
    for (int ir=mmin;ir<=mmax;ir++)
    {
      int mmin1=test_years(ir).indexmin();
      int mmax1=test_years(ir).indexmax();
      for (int ip=mmin1;ip<=mmax1;ip++) 
      {
        int oy=test_years(ir,ip);
        int om=test_months(ir,ip);
        int ow=test_weeks(ir,ip);
        test_move_index(ir,ip)=get_movement_index(om,ow,mo);
      }
    }
  }
  imatrix temp(1,num_regions,1,3000);
  temp.initialize();
  
  for (ir=1;ir<=num_regions;ir++)
  {
    for (ii=1;ii<=nfp_start(ir);ii++)
    {
      year_counter = test_years(ir,ii);
      temp(ir,year_counter)+=1;
    }
    year_counter=elem(regmin(ir)).year;
    temp(ir,year_counter)+=1;

    for (i=regmin(ir);i<=regmax(ir)-1;i++)
    {
      if (elem(i).year<elem(i+1).year)
      {
        // new year
        year_counter+=elem(i+1).year-elem(i).year;
        // since there is at least one fishery in this year !!!
        // what if there isn't ?????
        if (temp(ir).indexmax()< year_counter)
        {
          ad_exit(1);
        }
        temp(ir,year_counter)+=1;
      }
      else if (elem(i).month < elem(i+1).month)
      {
        // So there is another fishing period in this year
        if (temp(ir).indexmax()< year_counter)
        {
          ad_exit(1);
        }
        temp(ir,year_counter)+=1;
      }
    }
    for (ii=0;ii<nfp_end(ir);ii++)
    {
      year_counter = test_years(ir,nfp(ir)-ii);
      temp(ir,year_counter)+=1;
    }
  }
  int min_year=elem(regmin(1)).year;
  int max_year=elem(regmax(1)).year;
  cout << regmax(1) << endl;
  for (ir=2;ir<=num_regions;ir++)
  {
    if (min_year>elem(regmin(ir)).year)
      min_year=elem(regmin(ir)).year;
    if (max_year<elem(regmax(ir)).year)
      max_year=elem(regmax(ir)).year;
  }
  int nyrs=max_year-min_year+1;
  pmulti_species_data  pmsd=0;
  ivector nage_by_region(1,num_regions);
  ivector nage_by_fishery(1,nfsh);
  nage_by_region=nage;
  nage_by_fishery=nage;
  dvar_fish_stock_history fsh(ntg,num_regions,nage,nfp,nfi,nfsh,nyrs,
    _parest_flags,fl1,nft,regmin,regmax,dataswitch,Dflags,mo,_direction_flag,
    _mfactor,_season_region_flags,pmsd,nage_by_region,nage_by_fishery);

  fsh.age_nage=fishery_freq_record::age_nage;
  fsh.age_age1=fishery_freq_record::age_age1;

  breakup_big_function(nfp,nft,fsh,nfi,*this,temp,test_years,
    test_months,test_weeks,test_move_index,header_movement_index);

  if (!allocated(fsh.move_flags)) fsh.move_flags.allocate(test_move_flags);
  fsh.move_flags=test_move_flags;
  //fsh.do_recruitment_period_calculations();
  //fsh.recr.allocate(1,fsh.num_recruitment_periods);
//  fsh.recr.allocate(1,fsh.nyears);
  fsh.recr.allocate(2,fsh.nyears);   //NMD_19May2016
  fsh.N.allocate(1,num_regions,1,fsh.nyears,1,nage_by_region);
  for (int i=1;i<=num_regions;i++)
  {
    fsh.N(i)=1.e+225;
  }
  fsh.Nsave.allocate(1,num_regions,1,fsh.nyears,1,nage_by_region);
  fsh.Rsave.allocate(1,num_regions,1,fsh.nyears);
  fsh.exp_N.allocate(1,num_regions,1,fsh.nyears,1,nage_by_region);
  fsh.num_fish_data_recs=size();
  fsh.fishing_period=fishing_period;
  fsh.fishing_incident=fishing_incident;
  if (fsh.num_tag_releases)
    fsh.tag_flags.allocate(1,fsh.num_tag_releases,1,10);

  //fsh.make_fishing_period_report();
  fsh.get_fishery_realization_index();

  return fsh;
}

  
#undef HOME_VERSION

