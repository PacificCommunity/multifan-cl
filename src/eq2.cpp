/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#define HOME_VERSION
#include "all.hpp"

dmatrix submatrix(const dmatrix& _M,int lr,int ur,int lc,int uc)
{
  dmatrix & M=(dmatrix&)(_M);
  int nr=ur-lr+1;
  int nc=uc-lc+1;
  dmatrix tmp(1,nr);
  for (int i=1;i<=nr;i++)
  {
    tmp(i)=M(lr+i-1)(lc,uc).shift(1);
  }
  return tmp;
}

dvar_matrix submatrix(const dvar_matrix& _M,int lr,int ur,int lc,int uc)
{
  dvar_matrix & M=(dvar_matrix&)(_M);
  int nr=ur-lr+1;
  int nc=uc-lc+1;
  dvar_matrix tmp(1,nr);
  for (int i=1;i<=nr;i++)
  {
    tmp(i)=M(lr+i-1)(lc,uc).shift(1);
  }
  return tmp;
}

dvar_matrix dvar_fish_stock_history::get_initial_equilibrium_plus_group
(const dvar_matrix & eqrec, ivector * pq_flag)
{
//  ivector * pq_flag=0;
  int ns=age_flags(57);
  int nya=age_flags(95); // num_years_for_average
  //d3_array Z(1,ns,1,num_regions,1,nage); // Total mortality
  dvar3_array S(1,ns,1,nage,1,num_regions); // Survival
  dmatrix eps(1,ns,1,num_regions);

  // get the number of movmement periods for each season
  ivector  mps = get_equilibrium_movements_per_season();
  //  check for 1 movement per season
  for (int i=mps.indexmin();i<=mps.indexmax();i++)
  {
    if (mps(i) !=1)
    {
      cerr << "Only implemented for one movement per season" << endl;
      ad_exit(1);
    }
  }
  if (parest_flags(373))
  {
    pmsd_error();
    if (!pmsd)
    {
      if (!pq_flag)
      {
        get_kludged_initial_equilibrium_survival(mps,pq_flag);
      }
      else
      {
        //surv=get_initial_equilibrium_survival(mps,nya,pq_flag);
        if (!generate_report)
        {
          cout << "ERROR: internal sequence during F0eval should " << endl;
          cout << "be generating the report. Debug is required" << endl;
          cout << "Exiting... " << endl;
          ad_exit(1);
        }
        surv=dvar4_array(csurv);
      }
      kludged_surv=surv;
    }
  }
  else
  {
    if (!pmsd)
    {
      surv=get_initial_equilibrium_survival(mps,nya,pq_flag);
    }
    else
    {
      surv=get_initial_equilibrium_survival_ms(mps,nya,pq_flag);
    }
  }
  for (int i=1;i<=ns;i++)
  {
    S(i)=surv(i,1);   // for the 1 mp per season
  }

  // At this point we need at least the relative recruitment for the
  // first year for each region and season 

  // matrix for ageing the plus group  (age nage)
  // top vector element is season 1 region 1  This goes to season 2 region 1
  // killed off via the survival rate for 
  dmatrix MM(1,ns*num_regions,1,ns*num_regions);
  MM.initialize();
  dvar_matrix tmp(1,num_regions,1,num_regions);

  for (int is=1;is<=ns;is++)
  {
    tmp.initialize();
    for (int ir=1;ir<=num_regions;ir++)
    {
      tmp(ir,ir)=S(is,nage,ir); // at this point tmp is the diagonal matrix
                                // for survival
    }
    //tmp=Minv(is)*tmp;  // Now move the fish
    tmp=Dad(is,nage)*tmp;  // Now move the fish
    
    if (is<ns)
    {
      // The next line picks out the part of the "ageing-survival" matrix 
      // which corresponds to tmp
      dvar_matrix M1=submatrix(MM,is*num_regions+1,(is+1)*num_regions,
        (is-1)*num_regions+1,is*num_regions);
      M1=tmp;
    }
    else
    {
      dvar_matrix M1=submatrix(MM,1,num_regions,(is-1)*num_regions+1,
        is*num_regions);
      M1=tmp;
    }
  }

  int ii=0;
  /*
  cout << "Ageing or transition matrix for plus group"<< endl;
  for (int is=1;is<=ns;is++)
  {
    for (int i=1;i<=num_regions;i++)
    {
      cout << setfixed() << setprecision(3) << setw(5) << MM(++ii) << endl;
    }
    cout << endl;
  }
  */
  // get "recruitment" to plus group
  dvector x(1,num_regions*ns);
  random_number_generator rng(981);
  x.fill_randn(rng);
  x=10.0*exp(x);   // Generate the "recruitment" to the plus group

  dmatrix Id=identity_matrix(1,num_regions*ns);

  dvector eq= solve(Id-MM,x);
  dmatrix eqmat(1,ns,1,num_regions);
  int offset=0;
  for (int is=1;is<=ns;is++)
  {
    eqmat(is)=eq(offset+1,offset+num_regions).shift(1);
    offset+=num_regions;
  }

  return eqmat;
}

    

