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

void  dvar_fish_stock_history::copy_recruitment_to_population(dvar_matrix& Ninit,
  dvar3_array *pN)
{
  if (fabs(value(Ninit(1,1))-15.0)>1.e-3)
  {
    cout << " VV Here " << Ninit(1,1) << endl;
  }
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int iy=1;iy<=last_real_year;iy++)
    {
      if (pN)
      {
        (*pN)(ir,iy,1)=Ninit(ir,iy);
      }
      else
      {
        N(ir,iy,1)=Ninit(ir,iy);
      }
    }
  }
}

void  dvar_fish_stock_history::recinpop_standard(ivector * pq_flag,dvar3_array *pN,
  int override)
{
  //if (initial_orthp_estimate_flag || pN==0)
  {
    if (!allocated(Ninit_standard))
    {
      Ninit_standard.allocate(1,num_regions,1,last_real_year);
    }
    // _Ninit_standard.allocate(1,num_regions,1,last_real_year);
    Ninit_standard.initialize();
    // _Ninit_standard.initialize();
  }

  if (initial_orthp_estimate_flag==0 && pN != 0)
  {
    (*pN).initialize();
  }

  dvar_vector recr1=recr(2,last_real_year);  //NMD_19May2016
  if (parest_flags(400)==0)
  {
    recmean=mean(recr1);
  }
  else
  {
    if (norm2(recr(last_real_year-parest_flags(400)+1,last_real_year))
      > 1.e-10 && !parest_flags(398))
    {
      cerr << "This better not happen" << endl;
      // ad_exit(1);
    }
    recmean=mean(recr1(2,last_real_year-parest_flags(400)));  //NMD_19May2016
    if(parest_flags(398))
    {
      dvariable tmp=mean(exp(recr1(2,last_real_year-parest_flags(400))));  //NMD_19May2016
      recr1(last_real_year-parest_flags(400)+1,last_real_year)=log(tmp);
      recr(last_real_year-parest_flags(400)+1,last_real_year)=log(tmp);
    }
  }

  initmean=mean(initpop);
  if (pmsd && pmsd->num_species>1)
  {
    int ns=pmsd->num_species;
    for (int is=2;is<=ns;is++)
    {
      if (allocated(pmsd->recr1(is)))
        pmsd->recr1(is).deallocate();
      pmsd->recr1(is)=pmsd->recr(is)(2,last_real_year); //NMD_19May2016
      if (parest_flags(400)==0)  //NMD_19May2016
      {
        pmsd->recmean(is)=mean(pmsd->recr1(is));
      }
      else
      {
        if (norm2(pmsd->recr(is)(last_real_year-parest_flags(400)+1,last_real_year))
          > 1.e-10 && !parest_flags(398))
        {
          cerr << "This better not happen" << endl;
          ad_exit(1);
        }
        pmsd->recmean(is)=mean(pmsd->recr1(is)(2,last_real_year-parest_flags(400)));
        if(parest_flags(398))
        {
          dvariable tmp=mean(exp(pmsd->recr1(is)(2,last_real_year-parest_flags(400))));
          pmsd->recr1(is)(last_real_year-parest_flags(400)+1,last_real_year)=log(tmp);
          pmsd->recr(is)(last_real_year-parest_flags(400)+1,last_real_year)=log(tmp);
        }  //NMD_19May2016
      }
    }
  }
 
  // Here is where we distribute the recruitment over the regions
  // we make alpha year dependent to vary the recruitment distribution by time
  if (age_flags(71))
  {
    dvar_matrix tmp=mfexp(region_rec_diff_coffs);
    int iy;
    if (!pmsd)
    {
      for (iy=1;iy<last_real_year;iy++)
      {
        region_rec_diff_sums(iy)=mean(region_rec_diff_coffs(iy));
        if (!age_flags(178))
        {
          region_rec_diffs(iy)=region_rec_diff_coffs(iy)
                  -region_rec_diff_sums(iy);
        }
        else
        {
          region_rec_diffs(iy)=region_rec_diff_coffs(iy);
          rec_delta(iy)=
                log(sum(exp(region_rec_diffs(iy)+pop_delta)));
        }
      }
    }
    else
    {
      for (int is=1;is<=pmsd->num_species;is++)
      {
        int lb=pmsd->region_bounds(is,1);
        int ub=pmsd->region_bounds(is,2);
        for (iy=1;iy<last_real_year;iy++)
        {
          if(is==1)
          {
            region_rec_diff_sums(iy)=
              mean(region_rec_diff_coffs(iy)(lb,ub));
            if (!age_flags(178))
            {
              region_rec_diffs(iy)(lb,ub)=region_rec_diff_coffs(iy)(lb,ub)
                -region_rec_diff_sums(iy);
            }
            else
            {
              region_rec_diffs(iy)(lb,ub)=region_rec_diff_coffs(iy)(lb,ub);
              
              rec_delta(iy)=
              log(sum(exp(region_rec_diffs(iy)(lb,ub)+pop_delta(lb,ub))));
            }
          }
          else
          {
            pmsd->region_rec_diff_sums(is,iy)=
              mean(region_rec_diff_coffs(iy)(lb,ub));
            
            if (!age_flags(178))
            {
              region_rec_diffs(iy)(lb,ub)=region_rec_diff_coffs(iy)(lb,ub)
                -pmsd->region_rec_diff_sums(is,iy);
            }
            else
            {
              region_rec_diffs(iy)(lb,ub)=region_rec_diff_coffs(iy)(lb,ub);

              pmsd->rec_delta(is,iy)=
              log(sum(exp(region_rec_diffs(iy)(lb,ub)+pop_delta(lb,ub))));
            }
          }
        }
      }
    }

    for (iy=last_real_year;iy<=nyears;iy++)
    {
      region_rec_diff_sums(iy)=1.2348e-6;
      region_rec_diffs(iy)=0.0;
    }
    if (pmsd)
    {
      for (int is=2;is<=pmsd->num_species;is++)
      {
        for (iy=last_real_year;iy<=nyears;iy++)
        {
          pmsd->region_rec_diff_sums(is,iy)=1.2348e-6;
        }
      }
    }
  }

  if (age_flags(101)) 
    cout << "Recruitment correlate parameter = " << setprecision(4) 
    << sv(25) << endl;
  int ir;
  for (ir=1;ir<=num_regions;ir++)
  {
    if (!age_flags(71))
    {
      for (int iy=1;iy<=last_real_year;iy++)
      {
        if (age_flags(5)==0)
        {
          dvariable tmprcr=0.0;  //NMD_21Jun2016
          if (!pmsd || pmsd->region_species_pointer(ir)==1)
          {
            if(iy==1)
              tmprcr=0.0;
            else
              tmprcr=recr1(iy);  //NMD_21Jun2016
            Ninit_standard(ir,iy)=pop_delta(ir)+tmprcr-recmean+totpop;  //NMD_21Jun2016
          }
          else
          {
            int is=pmsd->region_species_pointer(ir);
            if(iy==1)
              tmprcr=0.0;
            else
              tmprcr=pmsd->recr1(is,iy);  //NMD_21Jun2016
            Ninit_standard(ir,iy)=pop_delta(ir)+tmprcr  //NMD_21Jun2016
                -pmsd->recmean(is)+pmsd->totpop(is);
            }
          }
        else
        {
          dvariable tmprcr=0.0;  //NMD_21Jun2016
          if (!pmsd || pmsd->region_species_pointer(ir)==1)
          {
            if(iy==1)
              tmprcr=0.0;
            else
              tmprcr=recr1(iy);  //NMD_21Jun2016

            Ninit_standard(ir,iy)=pop_delta(ir)+tmprcr-recmean+rec_init_diff+totpop;
          }
          else
          {
            int is=pmsd->region_species_pointer(ir);
            if(iy==1)
              tmprcr=0.0;
            else
              tmprcr=pmsd->recr1(is,iy);  //NMD_21Jun2016
            Ninit_standard(ir,iy)=pop_delta(ir)+tmprcr  //NMD_21Jun2016
                -pmsd->recmean(is)+pmsd->rec_init_diff(is)+pmsd->totpop(is);
          }
        }
        if (age_flags(101))
        {
          Ninit_standard(ir,iy)+=sv(25)*rec_covars(iy);
        }
        if (age_flags(101))
        {
          Ninit_standard(ir,iy)+=sv(25)*rec_covars(iy);
        }
      }
    }
    else
    {
      for (int iy=1;iy<=last_real_year;iy++)
      {
        if (age_flags(5)==0)
        {
          if (!pmsd || pmsd->region_species_pointer(ir)==1)
          {
            dvariable tmprcr;  //NMD_19May2016
            if(iy==1)
              tmprcr=0.0;
            else
              tmprcr=recr1(iy);  //NMD_19May2016
            if (age_flags(178))
            {
              Ninit_standard(ir,iy)=pop_delta(ir)+tmprcr  //NMD_19May2016
                  -recmean+totpop+region_rec_diffs(iy,ir)-rec_delta(iy);
            }
            else
            {
              Ninit_standard(ir,iy)=pop_delta(ir)+tmprcr  //NMD_19May2016
                  -recmean+totpop+region_rec_diffs(iy,ir);
            }
          }
          else
          {
            int is=pmsd->region_species_pointer(ir);
            dvariable tmprcr;  //NMD_19May2016
            if(iy==1)
              tmprcr=0.0;
            else
              tmprcr=pmsd->recr1(is,iy);  //NMD_19May2016
            if (age_flags(178))
            {
              Ninit_standard(ir,iy)=pop_delta(ir)+tmprcr  //NMD_19May2016
                  -pmsd->recmean(is)+pmsd->totpop(is)
                  +region_rec_diffs(iy,ir)-pmsd->rec_delta(is,iy);
            }
            else
            {
              Ninit_standard(ir,iy)=pop_delta(ir)+tmprcr  //NMD_19May2016
                  -pmsd->recmean(is)+pmsd->totpop(is)
                  +region_rec_diffs(iy,ir);
            }
          }
        }
        else
        {
          if (!pmsd || pmsd->region_species_pointer(ir)==1)
          {
            dvariable tmprcr;  //NMD_19May2016
            if(iy==1)
              tmprcr=0.0;
            else
              tmprcr=recr1(iy);  //NMD_19May2016
            Ninit_standard(ir,iy)=pop_delta(ir)+tmprcr
                -recmean
                +rec_init_diff
                +totpop+region_rec_diffs(iy,ir);
            }
          else
          {
            int is=pmsd->region_species_pointer(ir);
            dvariable tmprcr;  //NMD_19May2016
            if(iy==1)
              tmprcr=0.0;
            else
              tmprcr=pmsd->recr1(is,iy);  //NMD_19May2016

            Ninit_standard(ir,iy)=pop_delta(ir)+tmprcr
                -pmsd->recmean(is) +pmsd->rec_init_diff(is)
                +pmsd->totpop(is)+region_rec_diffs(iy,ir);
            }
          }
        if (age_flags(101) && initial_orthp_estimate_flag==0)
        {
          Ninit_standard(ir,iy)+=sv(25)*rec_covars(iy);
        }
      }
    }
  }
  copy_recruitment_to_population(Ninit_standard,pN);
  {
    ofstream ofs("std-rec" + str(override));
    if (override==6)
    {
      ofs << "Early call from nnewlan.cpp" << endl;
    }
    ofs << "norm(Ninit_standard)" << endl;
    ofs << norm(Ninit_standard) << endl;
    ofs << "pop_delta" << endl;
    ofs << pop_delta << endl;
    ofs << "number of years= " << last_real_year << endl;
    ofs << setscientific() << setw(7) << setprecision(3) << Ninit_standard
        << endl << endl;;
  }
}
