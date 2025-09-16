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

#include "newmprot.hpp"

            void xxxxxyyyy(){ ;}

void mytrap(const prevariable& x)
{
  if (value(x)>5.0)
    cout <<"htotcafi.cpp " << x << endl;
}

dvariable total_catch_fit(dvar_fish_stock_history& fsh,int print_switch,
  int avg_calc_flag)
{
  cout << "VVV" << endl;
  int no_pool_flag=fsh.check_total_catch_pooling_flags();
  dvariable xy=0.;
  dvariable tmp1=0.;
  int ntimes;
  int catflg=sum(column(fsh.fish_flags,45));
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    if (fsh.parest_flags(142)==0)
    {
      ntimes=fsh.num_fish_periods(ir);
    }
    else
    {
      ntimes=min(fsh.parest_flags(142),fsh.num_fish_periods(ir));
    }
    for (int ip=1;ip<=ntimes;ip++)
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
      {
        if (fsh.obs_tot_catch(ir,ip,fi)>-0.5L)
        {
          if (catflg)
          {
            if (avg_calc_flag)
            {
              fsh.totalcatch_by_numbers(ir)+=fsh.tot_catch(ir,ip,fi);
              fsh.obstotalcatch_by_numbers(ir)+=fsh.obs_tot_catch(ir,ip,fi);
              fsh.numtotalcatch_by_numbers(ir)+=1;
            }
            else
            {
              tmp1+=.01*fsh.fish_flags(fsh.parent(ir,ip,fi),45)
                *square(log(1.+fsh.tot_catch(ir,ip,fi))
                 -log(1.+fsh.obs_tot_catch(ir,ip,fi)));
            }
          }
          else
          {
            if (avg_calc_flag)
            {
              fsh.totalcatch_by_numbers(ir)+=fsh.tot_catch(ir,ip,fi);
              fsh.obstotalcatch_by_numbers(ir)+=fsh.obs_tot_catch(ir,ip,fi);
              fsh.numtotalcatch_by_numbers(ir)+=1;
            }
            else
            {
              MY_DOUBLE_TYPE tmpll=0.0; // YT 20170626
              if (no_pool_flag || !fsh.pmsd || fsh.pmsd->fisn(ir,ip,fi)==1)
              {
                dvariable tmp2=.01*fsh.age_flags(144)*
                  square(log(1.+fsh.tot_catch(ir,ip,fi))
                  -log(1.+fsh.obs_tot_catch(ir,ip,fi)));
                tmp1+=tmp2;
                tmpll=value(tmp2); // YT 20170626
              }
              else
              {
                if (ir==fsh.pmsd->reg_in_catch(ir,ip,fi,1))
                {
                  int n=fsh.pmsd->fisn(ir,ip,fi); 
                  // assume the total catch is in fishery corresponding
                  // to fishery realization incident 1
                  MY_DOUBLE_TYPE tmp_obs_tot_catch=fsh.obs_tot_catch(ir,ip,fi);
                  dvariable tmp_tot_catch=fsh.tot_catch(ir,ip,fi);
                  for (int i=2;i<=n;i++)
                  {
                     int rr=fsh.pmsd->reg_in_catch(ir,ip,fi,i);
                     tmp_tot_catch+=fsh.tot_catch(rr,ip,fi);
                  }
                  dvariable tmp2=.01*fsh.age_flags(144)*
                    square(log(1.+tmp_tot_catch)
                    -log(1.+tmp_obs_tot_catch));
                  tmp1+=tmp2;
                  tmpll=value(tmp2); // YT 20170626
                }
              }
//              if(fsh.ppstf) // YT Added to include tot catch likelihood in the case of 2 sex model with only catch in number
	      if(fsh.ppstf && !sum(fsh.q_flag)) // NMD_13jul2017
              {
                fsh.ppstf->tot_catch_like_by_realization(ir,ip,fi)
                        =tmpll;
              }
            }
          }
        }
      }
    }
  }
  xy+=tmp1;
  if (print_switch)
  {
    cout << " after total catch fit = " << xy << endl;
    cout << "exp(q0)= " << setprecision(10) <<exp(fsh.q0) << "  "
       << setprecision(10) << exp(fsh.q0(1)) << endl;
    cout << "total catch contribution = "<< tmp1 << endl;
  }
  return xy;
}

void dvar_fish_stock_history::set_zero_effdev_flag(void)
{
  zero_effdev_flag.initialize();
  for (int fi=1;fi<=num_fisheries;fi++) {
    for (int nt=1;nt<=num_real_fish_times(fi);nt++) 
    {
      int ir=realization_region(fi,nt);
      int ip=realization_period(fi,nt);
      int ri=realization_incident(fi,nt);
      int year=data_fish_flags(2,fi);
      int month=data_fish_flags(3,fi);
      int ty=really_true_year(ir,ip)+year1-1;
      int tm=really_true_month(ir,ip);
      if (obs_tot_catch(ir,ip,ri)>=-0.5L) 
      {  
          zero_effdev_flag(fi,nt)=1;
      }
     /*
      else
      {
        // period is after real end of fishing periods i.e. part 
        // of projection
        if ( (ty>year) || (ty==year && tm >= month) )
        {
          zero_effdev_flag(fi,nt)=1;
        }
        else
        {
          cout << "zero catch for fishery " << fi 
             << " incident " << nt << endl; 
        }
      }
      */
    }
  }
}
void dvar_fish_stock_history::allocate_fishery_projection_flags(void)
{
  if (allocated(fishery_projection_flag))
    fishery_projection_flag.deallocate();
  fishery_projection_flag.allocate(1,num_regions,1,num_fish_periods,
    1,num_fish_incidents);
  fishery_projection_flag.initialize();
}

void dvar_fish_stock_history::set_fishery_projection_flags(void)
{

  for (int fi=1;fi<=num_fisheries;fi++) 
  {
    for (int nt=1;nt<=num_fish_times(fi);nt++) 
    {
      int ir=realization_region(fi,nt);
      int ip=realization_period(fi,nt);
      int ri=realization_incident(fi,nt);
      int year=data_fish_flags(2,fi);
      int month=data_fish_flags(3,fi);
      int ty=really_true_year(ir,ip)+year1-1;
      int tm=really_true_month(ir,ip);
      if ( (ty>year) || (ty==year && tm >= month) )
      {
        fishery_projection_flag(ir,ip,ri)=1;
      }
    }
  }
}
void dvar_fish_stock_history::allocate_grouped_fishery_projection_flags(void)
{
  ivector group_flags32=column(fish_flags,32);
  int ngroups=max(group_flags32);
  if (allocated(grouped_fishery_projection_flag))
    grouped_fishery_projection_flag.deallocate();
  grouped_fishery_projection_flag.allocate(1,num_regions,1,num_fish_periods,
    1,ngroups);
  grouped_fishery_projection_flag.initialize();
}

void dvar_fish_stock_history::set_grouped_fishery_projection_flags(void)
{
  ivector group_flags32=column(fish_flags,32);
  int ngroups=max(group_flags32);
  
  for (int fi=1;fi<=num_fisheries;fi++) 
  {
    for (int nt=1;nt<=num_fish_times(fi);nt++) 
    {
      int ir=realization_region(fi,nt);
      int ip=realization_period(fi,nt);
      int ri=realization_incident(fi,nt);
      int year=data_fish_flags(2,fi);
      int month=data_fish_flags(3,fi);
      int ty=really_true_year(ir,ip)+year1-1;
      int tm=really_true_month(ir,ip);
      if ( (ty>year) || (ty==year && tm >= month) )
      {
        int gp=group_flags32(parent(ir,ip,ri));
        grouped_fishery_projection_flag(ir,ip,gp)+=
          fishery_projection_flag(ir,ip,ri);
      }
    }
  }
}

void dvar_fish_stock_history::set_missing_totcatch_flags(void)
{
  zero_catch_for_period_flag.initialize();
  zero_catch_for_incident_flag.initialize();
  missing_catch_by_fishery_flag.initialize();
  missing_catch_for_period_flag.initialize();
  missing_catch_flag=0;
  missing_catch_for_incident_flag.initialize();

  for (int fi=1;fi<=num_fisheries;fi++) 
  {
    for (int nt=1;nt<=num_fish_times(fi);nt++) 
    {
      int ir=realization_region(fi,nt);
      int ip=realization_period(fi,nt);
      int ri=realization_incident(fi,nt);
      if (obs_tot_catch(ir,ip,ri)<-0.5L) 
      {  
        missing_catch_by_fishery_flag(fi)+=1;
        missing_catch_for_period_flag(ir,ip)+=1; 
        missing_catch_for_incident_flag(ir,ip,ri)=1; 
        missing_catch_flag+=1;
      }
      if (obs_tot_catch(ir,ip,ri)==0.0) 
      {  
        zero_catch_for_period_flag(ir,ip)+=1; 
        zero_catch_for_incident_flag(ir,ip,ri)=1; 
      }
    }
  }
  if (missing_catch_flag>0)
  {
    fm_level_devs.allocate(1,num_fisheries,2,missing_catch_by_fishery_flag);
  }
}

void dvar_fish_stock_history::check_for_no_catch_or_effort(void)
{
  for (int fi=1;fi<=num_fisheries;fi++) 
  {
    for (int nt=1;nt<=num_fish_times(fi);nt++) 
    {
      int ir=realization_region(fi,nt);
      int ip=realization_period(fi,nt);
      int ri=realization_incident(fi,nt);
      {  
        if (missing_catch_for_incident_flag(ir,ip,ri) && missing_effort_by_region_flag(ir,ip,ri))
        {
          cerr << "Missing both catch and effort for fishery " << fi << " realization " << nt 
               << " year " << really_true_year(ir,ip)+year1 -1  << " month " << really_true_month(ir,ip) << endl;
          ad_exit(1);
        }
      }
    }
  }
}
#undef HOME_VERSION
