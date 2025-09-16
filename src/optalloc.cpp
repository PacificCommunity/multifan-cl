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
  dvariable mfexp(BOR_CONST prevariable& );
#endif

void dvar_fish_stock_history::set_fish_mort_q0(double u)
{
  if (allocated(fish_mort_q0))
  {
    for (int i=fish_mort_q0.indexmin();i<=fish_mort_q0.indexmax();
      i++)
    {
      fish_mort_q0(i)=u;
    }
  }
}

void dvar_len_fish_stock_history::allocate_optional_stuff(void)
{
  if (parest_flags(141)==9 || parest_flags(141)==10
     || parest_flags(141)==12
     || parest_flags(139)==9 || parest_flags(139)==10
     || parest_flags(139)==12)
  {
    if (allocated(m_estimator_coffs))
    {
      m_estimator_coffs.deallocate();
    }
    /*
    ifstream ifs("parests.5_200_new");
    if (!ifs)
    {
      cerr << "Error opening file parests.5_200_new" << endl;
      ad_exit(1);
    }
    */
   /*
    int nc,mmin,mmax;
    ifs >> nc >> mmin >> mmax;
    if (!ifs)
    {
      cerr << "Error reading from parests.out" << endl;
      ad_exit(1);
    }
    
    m_estimator_coffs.allocate(mmin,mmax,1,25);
    ifs >> m_estimator_coffs;
    if (!ifs)
    {
      cerr << "Error reading from parests.out" << endl;
      ad_exit(1);
    }
   */
    m_estimator_coffs=get_learner_code_coffs_low();
  }
    
  int sum_qflag=sum(column(fish_flags,55));
  if (age_flags(92)==3) ss2_flag=1;

  if (sum_qflag || age_flags(170))
  {
    if (allocated(fish_mort_q0)) fish_mort_q0.deallocate();
    fish_mort_q0.allocate(1,num_regions,1,num_fish_periods,
      1,num_fish_incidents,1,nage);

    if (allocated(tot_mort_q0)) tot_mort_q0.deallocate();
    tot_mort_q0.allocate(1,num_regions,1,num_fish_periods,1,nage);

    if (allocated(num_fish_q0))
    {
      num_fish_q0.deallocate();
    }
    num_fish_q0.allocate(1,num_regions,1,num_fish_periods,1,nage);
    if (allocated(survival_q0))
    {
      survival_q0.deallocate();
    }
    survival_q0.allocate(1,num_regions,1,num_fish_periods,1,nage);
    survival_q0.initialize();
    if (allocated(catchability_q0))
    {
      catchability_q0.deallocate();
    }
    catchability_q0.allocate(1,num_regions,1,num_fish_periods,
      1,num_fish_incidents);

    if (allocated(N_q0))
    {
      N_q0.deallocate();
    }
    N_q0.allocate(1,num_regions,1,nyears,1,nage);

    if (allocated(fish_mort_calcs_q0))
    {
      fish_mort_calcs_q0.deallocate();
    }
    fish_mort_calcs_q0.allocate(1,num_regions,1,num_fish_periods,
      1,num_fish_incidents,1,nage);

    if (allocated(catch_q0))
    {
      catch_q0.deallocate();
    }
    catch_q0.allocate(1,num_regions,1,num_fish_periods,
      1,num_fish_incidents,1,nage);
    N_q0.initialize();
    catch_q0.initialize();
    set_fish_mort_q0(-70.0);
  }
  if (age_flags(92))
  {
    qmu.allocate(1,num_fisheries,1,num_fish_times);
    qvar.allocate(1,num_fisheries,1,num_fish_times);
    qstudent.allocate(1,num_fisheries,1,num_fish_times);
  }
  if(check(column(fish_flags,26),3))
  {
    //age_weight_incident_sel.allocate(1,num_regions,1,num_fish_periods,
     // 1,num_fish_incidents,1,nage,1,nwint);
    //age_len_incident_sel.allocate(1,num_regions,1,num_fish_periods,
     // 1,num_fish_incidents,1,nage,1,nlint);
    //age_len_incident_sel_yr.allocate(1,age_flags(57),1,num_fisheries,
     // 1,nage,1,nlint);
    //age_weight_sel.allocate(1,12,1,4,1,num_fisheries,1,nage,1,nwint);
    //wlengths.allocate(1,nwint);
  }
  if (parest_flags(240))
  {
    qij.allocate(1,num_regions);
    for (int i=1;i<=num_regions;i++)
    {
      qij(i).allocate(1,num_fish_periods(i));
      for (int j=1;j<=num_fish_periods(i);j++)
      {
        qij(i,j).allocate(1,num_fish_incidents(i,j));
      }
    }
  }
  int ns;
  if (age_flags(57)==0)
  {
    age_flags(57)=1;
  }
  if (parest_flags(197)==1)
  {
    ns=1;
  }
  else
  {
    ns=age_flags(57);
  }
  int rem=last_real_year%ns;
  na=last_real_year/ns;
  if (rem) na++;

  if (!af170q0)  //NMD_9Feb2014
  {
    if (age_flags(94)==3  || age_flags(182) )
    {
      __SAFE_DEALLOCATE__(bh_recr_devs)
      bh_recr_devs.allocate(1,na);
      __SAFE_DEALLOCATE__(bh_bio)
      bh_bio.allocate(1,na);
      __SAFE_DEALLOCATE__(bh_predicted_recruits)
      bh_predicted_recruits.allocate(1,na);
      if (pmsd)  //NMD_12Apr2021
      {
        __SAFE_DEALLOCATE__(pmsd->bh_predicted_recruits)	  
	  pmsd->bh_predicted_recruits.allocate(1,pmsd->num_species,1,na);
      }  //NMD_12Apr2021
      __SAFE_DEALLOCATE__(bh_reproductive_biomass)
      bh_reproductive_biomass.allocate(1,na);
    }
    else
    {
      __SAFE_DEALLOCATE__(bh_recr_devs)
      bh_recr_devs.allocate(1,nyears);
      __SAFE_DEALLOCATE__(bh_bio)
      bh_bio.allocate(1,nyears);
      __SAFE_DEALLOCATE__(bh_predicted_recruits)
      bh_predicted_recruits.allocate(1,nyears);
      if (pmsd)  //NMD_12Apr2021
      {
        __SAFE_DEALLOCATE__(pmsd->bh_predicted_recruits)
        pmsd->bh_predicted_recruits.allocate(1,pmsd->num_species,1,nyears);
      }  //NMD_12Apr2021
      __SAFE_DEALLOCATE__(bh_reproductive_biomass)
      bh_reproductive_biomass.allocate(1,nyears);
    }
  }  //NMD_9Feb2014
  __SAFE_DEALLOCATE__(bh_numbers)
  bh_numbers.allocate(1,nyears);

  if (parest_flags(301)>0 && af170q0!=1) //exclude zero fishing case NMD_7jun-19
  {
    if (allocated(tc_wtprob))
    {
      tc_wtprob.deallocate();
    }
    if (allocated(tc_wght_freq))
    {
      tc_wght_freq.deallocate();
    }
    if (allocated(max_wght_obs))
    {
      max_wght_obs.deallocate();
    }
    tc_wtprob.allocate(1,num_regions,
      1,num_fish_periods,1,num_fish_incidents);
    tc_wght_freq.allocate(1,num_regions,
      1,num_fish_periods,1,num_fish_incidents);
    max_wght_obs.allocate(1,num_regions,
      1,num_fish_periods,1,num_fish_incidents);
  }     //  NMD_7jun-19
  if (parest_flags(311)>0 && af170q0!=1) //exclude zero fishing case NMD_7jun-19
  {
    if (allocated(tc_tprob))
    {
      tc_tprob.deallocate();
    }
    if (allocated(tc_len_freq))
    {
      tc_len_freq.deallocate();
    }
    if (allocated(max_len_obs))
    {
      max_len_obs.deallocate();
    }
    tc_tprob.allocate(1,num_regions,
      1,num_fish_periods,1,num_fish_incidents);
    tc_len_freq.allocate(1,num_regions,
      1,num_fish_periods,1,num_fish_incidents);
    max_len_obs.allocate(1,num_regions,
      1,num_fish_periods,1,num_fish_incidents);
  }     //  NMD_7jun-19
}

