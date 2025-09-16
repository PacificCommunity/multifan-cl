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

  extern dvar_vector * psv;

void setup_Dyf_age_dep(dvar3_array& Dad,int nage,int num_regions,
  imatrix& Dflags, dvar_vector& diff_coffs,dvar_vector& diff_coffs2,
  dvar_vector& diff_coffs3,ivector& age_flags);


dvar_matrix dvar_fish_stock_history::get_equilibrium_survival_rate(void)
{
  dvar3_array _Dad = Dad(1);
  int ir;
  int ny=0;
//  MY_DOUBLE_TYPE nmult=age_flags(95)/10.;
  MY_DOUBLE_TYPE nmult=1.0;  //NMD_21jun2016
//  if (age_flags(128)>0)
  if (age_flags(128)>0 || !af170q0)   //NMD_6jun2022
  {
    nmult=age_flags(128)/100.;  //NMD_4jun2025
  }
  dvar_matrix es(1,num_regions,1,nage);
  es.initialize();
  switch(age_flags(94))
  {
  case 1:    
    if (!pmsd)
    {
      for (ir=1;ir<=num_regions;ir++)
      {
        es(ir)=-nmult*exp(nat_mort(1));
      }
    }
    else
    {
      for (ir=1;ir<=num_regions;ir++)
      {
        es(ir)=-nmult*exp(get_nat_mort_region(ir)(1));
      }
    }
    break;

  case 2:    
    if (age_flags(92))
    {
      cerr << endl << "Error -- can't use beginning fishing  mortalities to" << endl;
      cerr << " calculate the initial equilbrium age structure when" << endl;
      cerr << " employing catch-conditioned (total catches) option" << endl;
      cerr << " i.e. age_flags(92) !=0. Instead You must set age flags(94) to 1" << endl;
//      cerr << "  use nat mort with age_flags(95)/10. as multiplier" << endl;
      cerr << "  use nat mort with age_flags(128)/100. as multiplier" << endl; //NMD_21jun2016
      ad_exit(1);
    }
    ny=age_flags(95);
    for (ir=1;ir<=num_regions;ir++)
    {
      for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
      {
        if (ny<year(ir,ip)) break;
        dvar_matrix* pfm;
        if (af170q0)
        {
          pfm=&(fish_mort_q0(ir,ip));
        }
        else
        {
          pfm=&(fish_mort(ir,ip));
        }
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
        {                                         // incidents for this period
          es(ir)-=mfexp((*pfm)(fi));
        }
      }
      es(ir)/=double(ny);;
      if (!pmsd)
        es(ir)-=mfexp(nat_mort(year(ir,1)));
      else
        es(ir)-=mfexp(get_nat_mort_region(ir)(year(ir,1)));
    }
    //cout <<"equilib.cpp " << endl << es << endl;
    break;
  case 3:    
    if (!pmsd)
    {
      for (ir=1;ir<=num_regions;ir++)
      {
        es(ir)=-exp(nat_mort(1));
      }
    }
    else
    {
      for (ir=1;ir<=num_regions;ir++)
      {
        es(ir)=-exp(get_nat_mort_region(ir)(1));
      }
    }
    break;
  default:
    cerr << "illegal value for age_flags(94) "
      "in get_equilibrium_survival_rate" << endl;
  }
  return es;
}

void dvar_fish_stock_history::get_initial_age_structure_equilibrium(void)
{
  // get the initial age structure assuming equilibrium conditions
  // the inputs are the (log) equilibrium recruitment vector
  // and the (log) equilibrium survivial rate matrix
  dvar3_array _Dad = Dad(1);
  dvar_matrix les=get_equilibrium_survival_rate();

  dvar_matrix tmpN(1,num_regions,1,nage);
  tmpN.initialize();
  if (num_regions==1) 
  {
    tmpN(1,1)=N(1,1,1);
    for (int j=1;j<nage-1;j++)
    {
      tmpN(1,j+1)=les(1,j)+tmpN(1,j);
    }
    tmpN(1,nage)=les(1,nage-1)+tmpN(1,nage-1);
    dvariable tt= log(0.99999-mfexp(les(1,nage)));
    tmpN(1,nage)-= tt;
  }
  else // do equilbrium calculations for multi region situation
  {
    int ir;
    dvar_matrix es=mfexp(les);
    dvar_vector tt(1,num_regions);
    tt.initialize();
    for (ir=1;ir<=num_regions;ir++)
    {
      int ny=age_flags(95);
      for (int iy=2;iy<=ny;iy++)
        tt(ir)+=mfexp(N(ir,iy,1));
      tt(ir)/=(ny-1);
      //tt(ir)=mfexp(N(ir,1,1));
    }
    //cout <<"equilib.cpp " << tt << endl;
    tt=_Dad(1)*tt;
    //cout <<"equilib.cpp " << tt << endl;
    for (ir=1;ir<=num_regions;ir++)
    {
      tmpN(ir,1)=tt(ir);
    }
    for (int j=2;j<=nage;j++)
    {
      for (ir=1;ir<=num_regions;ir++)
      {
        tt(ir)=es(ir,j-1)*tmpN(ir,j-1);
      }
      tt=_Dad(j)*tt;
      for (ir=1;ir<=num_regions;ir++)
      {
        tmpN(ir,j)=tt(ir);
      }
    }
    // now for the last age class
    dvar_matrix meq(1,num_regions,1,num_regions);
    meq.initialize();
    for (ir=1;ir<=num_regions;ir++)
    {
      for (int ir1=1;ir1<=num_regions;ir1++)
      {
        meq(ir,ir1)=-es(ir1,nage)*_Dad(nage)(ir,ir1);
      }
    }
    for (ir=1;ir<=num_regions;ir++)
    {
      meq(ir,ir)+=1;
    }
    tt=inv(meq)*tt;
    for (ir=1;ir<=num_regions;ir++)
    {
      tmpN(ir,nage)=tt(ir);
    }
    tmpN=log(1.e-10+tmpN);
  }
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int j=1;j<=nage;j++)
    {
      N(ir,1,j)=tmpN(ir,j);
    }
  }

}

void dvar_len_fish_stock_history::get_equilibrium_structure_for_yield(void)
{
  // get the initial age structure assuming equilibrium conditions
  // the inputs are the (log) equilibrium recruitment vector
  // and the (log) equilibrium survivial rate matrix
  //dvar3_array _Dad = Dad(1);

  get_fishing_mortality_by_age_by_year_by_region();
  dvar3_array& Fayr=F_by_age_by_year_by_region;
  dvar_vector n1(1,nage);
  const int nsteps=500;
  dvector tmp1(0,nsteps);
  dvar_vector tmp2(0,nsteps);
  int i,ir;
  int ii=1;
  tmp1(0)=0.0;
  tmp2(0)=0.0;
  int navg=age_flags(148);
  int tmult=age_flags(57);
  if (!navg)navg=tmult;
  if (!navg)navg=1;

  dvar_matrix F(1,num_regions,1,nage);
  dvar_matrix les(1,num_regions,1,nage);
  dvar_matrix es(1,num_regions,1,nage);
  dvar_matrix Z(1,num_regions,1,nage);
  F.initialize();
  es.initialize();
  les.initialize();
  Z.initialize();
  for (ir=1;ir<=num_regions;ir++)
  {
    for (i=1;i<=navg;i++)
      F(ir)+=Fayr(ir,nyears-i+1);

    F(ir)/=double(navg);
    //cout <<"equilib.cpp " << F(ir) << endl;
  }

  //ofstream ofs("testequ");
  //ofstream ofs1("catchequ");
  for (i=1;i<=nsteps;i++)
  {
    MY_DOUBLE_TYPE lambda=i/10.0;
    if (!pmsd)
    {
      for (ir=1;ir<=num_regions;ir++)
      {
         Z(ir)=lambda*F(ir)+exp(nat_mort(nyears));
         es(ir)=exp(-Z(ir));
      }
    }
    else
    {
      for (ir=1;ir<=num_regions;ir++)
      {
         Z(ir)=lambda*F(ir)+exp(get_nat_mort_region(ir)(nyears));
         es(ir)=exp(-Z(ir));
      }
    }
    les=log(1.e-10+es);

    dvar_matrix tmpN(1,num_regions,1,nage);
    tmpN.initialize();
    if (num_regions==1) 
    {
      tmpN(1,1)=N(1,1,1);
      for (int j=1;j<nage-1;j++)
      {
        tmpN(1,j+1)=les(1,j)+tmpN(1,j);
      }
      tmpN(1,nage)=les(1,nage-1)+tmpN(1,nage-1);
      dvariable tt= log(0.99999-mfexp(les(1,nage)));
      tmpN(1,nage)-= tt;
    }
    else if (!pmsd|| pmsd->num_species==1) // do equilbrium calculations for multi region situation
    {
      
      int ir;
      dvar_vector tt=epop_delta;

      for (ir=1;ir<=num_regions;ir++)
      {
        tmpN(ir,1)=tt(ir);
      }
      int nmp=mo.num_periods();
      dvar_matrix es=mfexp(les/double(nmp));
      for (int j=1;j<nage;j++)
      {
        for (int mp=1;mp<=nmp;mp++)
        {
          if (num_regions>1)
          {
            tt=Dad(move_map(mp),j)*tt;
          }
          for (ir=1;ir<=num_regions;ir++)
          {
            tt(ir)=es(ir,j)*tt(ir);
          }
        }
        for (ir=1;ir<=num_regions;ir++)
        {
          tmpN(ir,j+1)=tt(ir);
        }
      }
      // now for the last age class
      {
        dvar_matrix meq(1,num_regions,1,num_regions);
        meq.initialize();
        for (ir=1;ir<=num_regions;ir++) meq(ir,ir)=1.0;
        for (int mp=1;mp<=nmp;mp++)
        {
          meq=Dad(move_map(mp),nage)*meq;
          for (ir=1;ir<=num_regions;ir++)
          {
            for (int ir1=1;ir1<=num_regions;ir1++)
            {
              meq(ir,ir1)=es(ir1,nage)*Dad(move_map(mp),nage)(ir,ir1);
            }
          }
        }
        for (ir=1;ir<=num_regions;ir++)
        {
          meq(ir,ir)-=1;
        }
        tt=inv(meq)*tt;
        tt=-tt; // cause we used minus the matrix
        for (ir=1;ir<=num_regions;ir++)
        {
          tmpN(ir,nage)=tt(ir);
        }
      }
    }
    else
    {
      
      int ir;

      dvar_vector tt1=epop_delta;

      int nrr=pmsd->num_real_regions;
      int ns=pmsd->num_species;

      for (int isp=1;isp<=ns;isp++)
      {
        if (pmsd->num_real_regions==1)    // - no movement NMD_jan28-19
        {
          ivector isb=get_region_bounds(isp);
	  int irr=max(isb);
          tmpN(irr,1)=tt1(irr);
          for (int j=1;j<nage-1;j++)
          {
            tmpN(irr,j+1)=les(irr,j)+tmpN(irr,j);
          }
          tmpN(irr,nage)=les(irr,nage-1)+tmpN(irr,nage-1);
	  
          dvariable ttt= log(0.99999-mfexp(les(irr,nage)));
          tmpN(irr,nage)-= ttt;
        }
        else    //NMD_jan28-19
        {	
          int offset=(isp-1)*nrr;
          dvar_vector tt=tt1(1+offset,nrr+offset).shift(1);
          dvar_matrix tmpN1=tmpN.sub(1+offset,nrr+offset).shift(1);
          dvar_matrix les1=les.sub(1+offset,nrr+offset).shift(1);
          for (ir=1;ir<=nrr;ir++)
          {
            tmpN1(ir,1)=tt(ir);
          }
          int nmp=mo.num_periods();
          dvar_matrix es=mfexp(les1/double(nmp));	  
          for (int j=1;j<nage;j++)
          {
            for (int mp=1;mp<=nmp;mp++)
            {
              dvar_matrix Dd=
                Dad(move_map(mp),j).sub(1+offset,nrr+offset).shift(1);
              if (nrr>1)
              {
                tt=Dd*tt;
              }
              for (ir=1;ir<=nrr;ir++)
              {
                tt(ir)=es(ir,j)*tt(ir);
              }
            }
            for (ir=1;ir<=nrr;ir++)
            {
              tmpN1(ir,j+1)=tt(ir);
            }
          }
          // now for the last age class
          {
            dvar_matrix meq(1,nrr,1,nrr);
            meq.initialize();
            for (ir=1;ir<=nrr;ir++) meq(ir,ir)=1.0;
            for (int mp=1;mp<=nmp;mp++)
            {
              dvar_matrix Dd=
                Dad(move_map(mp),nage).sub(1+offset,nrr+offset).shift(1);
              meq=Dd*meq;
              //meq=Dad(move_map(mp),nage)*meq;
              for (ir=1;ir<=nrr;ir++)
              {
                for (int ir1=1;ir1<=nrr;ir1++)
                {
                  meq(ir,ir1)=es(ir1,nage)*Dd(ir,ir1);
                  //meq(ir,ir1)=es(ir1,nage)*Dad(move_map(mp),nage)(ir,ir1);
                }
              }
            }
            for (ir=1;ir<=nrr;ir++)
            {
              meq(ir,ir)-=1;
            }
            tt=inv(meq)*tt;
            tt=-tt; // cause we used minus the matrix
            for (ir=1;ir<=nrr;ir++)
            {
              tmpN1(ir,nage)=tt(ir);
            }
          }
        }  //end of movement case      //NMD_jan28-19
      } 
    }
  }
}
