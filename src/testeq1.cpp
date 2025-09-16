/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#define HOME_VERSION
#include "all.hpp"
      void check_sanity(imatrix& mpy);
      void check_sanity_trans(imatrix& mpy);

void check_code1(void)
{
  double x=1.0;
  double y=sqrt(x);
}

void check_code(void)
{
  check_code1();
}




  void print_P(dvar_matrix& P)
  {
    dvar_matrix TP=trans(P);
    ofstream ofs("testn");
    int imin=TP.indexmin();
    int imax=TP.indexmax();
    for (int i=imin;i<=imax;i++)
    {
      ofs << setprecision(3) << setscientific() << setw(11) 
          << TP(i) << endl;
    }
  }

 
  dvar_matrix mycolumn(const dvar_matrix& _M,int i,int j)
  {
    ADUNCONST(dvar_matrix,M)
    int mmin=M.indexmin();
    int mmax=M.indexmax();
    dvar_matrix tmp(mmin,mmax);
    for (int ii=mmin;ii<=mmax;ii++)
    {
      tmp(ii)=M(ii)(i,j);
    }
    return tmp;
  }
 

void dvar_fish_stock_history::test_initial_equilibrium_code(ivector * pq_flag)
{
  int nya=age_flags(95); // num_years_for_average
  if (nya==0 && age_flags(94)==2)
  {
    cerr << "Error age_flags(95) must be >0 " << endl;
    ad_exit(1);
  }

  // get the number of movmement periods for each season
  ivector  mps = get_equilibrium_movements_per_season();

  imatrix equilib_move_index = get_initial_equilibrium_movement(mps);

  //dvar4_array  surv;
  
  if (parest_flags(373) && parest_flags(374)>0)
  {
    pmsd_error();
    if (!pmsd)
    {
      // just to have some reasonable survival numbers to work with for now
      //surv=get_initial_equilibrium_survival_for_catch_conditioned(mps,
      //  nya,pq_flag);

      //This will be the real code
      if (!pq_flag)
      {
        switch(age_flags(94))
        {
        case 0:
//          get_kludged_initial_equilibrium_survival(mps,pq_flag);
            cerr << " Incorrect option N_init survival in catch-conditioned "
              << " case - exiting" << endl;
            ad_exit(1);
          break;
        case 1:
          kludged_surv=get_initial_equilibrium_survival_natmort
            (mps,pq_flag);
          break;
        case 2:
          if (age_flags(94)==2 && nya==0)
          {
            cout << "Initial exploited equilibrium population requires"
                 << "parest_flags(95)>0 " << endl;
            ad_exit(1);
          }
          kludged_surv=get_kludged_initial_equilibrium_survival
            (mps,pq_flag);
          break;
        default:
          cerr << "bad switch age_flags(94) " << endl;
          ad_exit(1);
        }
      }
      else
      {
        //surv=get_initial_equilibrium_survival(mps,nya,pq_flag);
        if (!generate_report && parest_flags(145)!=2 && parest_flags(145)!=3
            && parest_flags(145)!=7)
        {
          cout << "Catch-conditioned initial survival - " << endl;
          cout << "Fzero evaluation for Hessian calcs - incorrect " << endl;
          cout << "evaluation for report & parest_flags(145) settings " << endl;
          ad_exit(1);
        }
//        surv=dvar4_array(csurv);
        kludged_surv=dvar4_array(csurv);  //NMD_6May2022
      }
      //kludged_surv=surv;
    }
  }
  else
  {
    if (!pmsd)
    {
      if (age_flags(94)==1)
      {
        if (age_flags(92)==0) //NMD_25jan2022 - ccond or cerrs
        {
          surv=get_initial_equilibrium_survival_natmort
              (mps,pq_flag);
        }
        else
        {
          kludged_surv=get_initial_equilibrium_survival_natmort
            (mps,pq_flag);
        }
      }
      else
      {
//        kludged_surv=get_initial_equilibrium_survival
//          (mps,nya,pq_flag);
        if (age_flags(92)==0)
        {
          surv=get_initial_equilibrium_survival
            (mps,nya,pq_flag);
        }
        else
        {
          cerr << " Invalid oprtion for survival calc" << endl;
          ad_exit(1);
        }
      }
    }
    else
    {                 //NMD_7feb2025
      if (age_flags(92)==2 && age_flags(94)==1)
      {
        kludged_surv=get_initial_equilibrium_survival_natmort
              (mps,pq_flag);
      }
      else if(age_flags(92)==0)
      {
        surv=get_initial_equilibrium_survival_ms(mps,nya,pq_flag);	
      }
      else
      {
        pmsd_error();
      }
    }
  }

  // new code for initial equilibrium population  just here to see 
  //what it needs and then to make a variable object
  //int ns=age_flags(57);
  dvar_matrix  eqrec=get_initial_equilibrium_recruitment(nya);
  //cout << "BBX " << sum(eqrec) << endl;
  // Next call is ignored - in development for as part of catch-conditioning
  // NMD_2dec2020
  //  dvar_matrix eqmat=get_initial_equilibrium_plus_group(eqrec, pq_flag);
  //

  dvar3_array EN=get_equilibrium_cohorts(mps,eqrec,equilib_move_index);

  // P is the fish up to age class nage-1
  dvar_matrix P=get_initial_equilibrium_population1(EN);

  dvar_vector A=get_initial_equilibrium_A(mps,P,equilib_move_index);

  dvar_matrix B=get_initial_equilibrium_B(mps,equilib_move_index);

  // X should contain the equilbrium number of age classs nage fish
  dvar_vector X(1,num_regions);
  if (!pmsd)
  {
    dmatrix Id=identity_matrix(1,num_regions);
    X=solve(Id-B,A);
    Xchk(iloop)=value(X);
  }
  else
  {
    int nrr=pmsd->num_real_regions;
    int numsp=pmsd->num_species;
    dmatrix Id=identity_matrix(1,nrr);
    int isp;
    for (isp=1;isp<=numsp;isp++)  // loop over species
    {
      int offset=(isp-1)*nrr;
      X(1+offset,nrr+offset).shift(1)
        =solve(Id-mycolumn(B,1+offset,nrr+offset).colshift(1),
          A(1+offset,nrr+offset).shift(1));
//      cout << "II" << isp << endl 
//           << mycolumn(B,1+offset,nrr+offset).colshift(1) << endl
//           << A(1+offset,nrr+offset).shift(1) << endl;
    }
    for (isp=1;isp<=numsp;isp++)  // loop over species
    {
      int offset=(isp-1)*nrr;
//      cout << "LL" << isp << endl << X(1+offset,nrr+offset) << endl;
    }
  }
  // put N into the initial equilbriumn population matrix
  put_X_in_P(X,P);

  // print it out to see what it looks like
  //print_P(P);
  //ad_exit(1);
  // put EN into N
  put_P_in_N(P);
}

void dvar_fish_stock_history::put_X_in_P(dvar_vector& X,dvar_matrix& P)
{
  P(nage)=X;
}

void dvar_fish_stock_history::put_P_in_N(dvar_matrix& P)
{
  dvar_matrix TP=trans(P);
  int ir;
  if (af170q0==0)
  {
    for (ir=1;ir<=num_regions;ir++)
    {
      N(ir,1)=log(1.e-10+TP(ir));
    }
  }
  else
  {
    for (ir=1;ir<=num_regions;ir++)
    {
      for (int j=1;j<=nage;j++)
      {
        N_q0(ir,1)=log(1.e-10+TP(ir));
      }
    }
  }
}

dvar3_array dvar_fish_stock_history::get_equilibrium_cohorts
  (ivector& mps,dvar_matrix& eqrec,imatrix& emi)
{
  int ns=age_flags(57);
  dvar3_array EN(1,nage,1,ns,1,num_regions);
  EN.initialize();
  int is,im;
  int ir=0;
  if (!pmsd)
  {
    int ng1=nage-1;
    if (age_flags(94)==3) ng1=nage;
    for (is=1;is<=ns;is++)
    {
      EN(1,is)=eqrec(is);
      dvar_vector tmp(1,num_regions);
      for (int j=1;j<ng1;j++)
      {
        int si=get_season(is,j);
        tmp=EN(j,is);
        for (im=1;im<=mps(si);im++)
        {
          if (!parest_flags(394))
          {
            if (age_flags(92)==0)
            {
              tmp = elem_prod(tmp,surv(si,im,j));
            }
            else
            {
              tmp = elem_prod(tmp,kludged_surv(si,im,j));
            }
          }
          else
          {
            tmp = 0.8*tmp;
          }
          // now do the movement 
          if (num_regions>1)
          {
            tmp=Dad(emi(si,im),j)*tmp;
          }
        } 
        EN(j+1,is)=tmp;
      }
    }
  }
  else   // mult species code
  {
    int nrr=pmsd->num_real_regions;
    int numsp=pmsd->num_species;

    dvar3_array TEN(1,nage,1,ns,1,nrr);
    TEN.initialize();
    dvar_vector tmp(1,nrr);
    for (int isp=1;isp<=numsp;isp++)  // loop over species
    {
      int offset=(isp-1)*nrr;
      for (is=1;is<=ns;is++)
      {
        EN(1,is)(1+offset,nrr+offset).shift(1)
          =eqrec(is)(1+offset,nrr+offset).shift(1);
        TEN(1,is)=eqrec(is)(1+offset,nrr+offset).shift(1);
        dvar_vector tmp(1,nrr);
        for (int j=1;j<nage-1;j++)
        {
          int si=get_season(is,j);
          tmp=TEN(j,is);
          for (im=1;im<=mps(si);im++)
          {
            if (!parest_flags(394))
            {
              if (age_flags(92)==0)
              {
                tmp = 
                elem_prod(tmp,surv(si,im,j).sub(1+offset,nrr+offset).shift(1));
              }
              else
              {
                tmp = 
                elem_prod(tmp,
                kludged_surv(si,im,j).sub(1+offset,nrr+offset).shift(1));
              }
            }
            else
            {
              tmp = 0.8*tmp;
            }
//            tmp = 
//              elem_prod(tmp,surv(si,im,j).sub(1+offset,nrr+offset).shift(1));
            // now do themovement 
            if (num_regions>1)
            {
              tmp=
                Dad(emi(si,im),j).sub(1+offset,nrr+offset).shift(1)*tmp;
            }
          } 
          //if (j>=nage-3)
          //  cout <<"testeq1.cpp " << j << endl;
          TEN(j+1,is)=tmp;
          EN(j+1,is).sub(1+offset,nrr+offset).shift(1)=tmp;
            //cout << tmp << endl;
        }
      }
    }
  }

  /*
  {
    ofstream ofs("EN");
    int rmin=1;
    int rmax=6;
    if (pmsd)
    {
      rmin=7;
      rmax=12;
    }
    for (is=1;is<=ns;is++)
    {
     for (int ir=rmin;ir<=rmax;ir++)
     {
       for (int j=1;j<=nage;j++)
       {
         ofs <<setw(10) << EN(j,is,ir) << " ";
       }
       ofs << endl;
     }
     ofs << endl;
    }
  }
  //cout << EN << endl;
  */
  return EN;
}



dvar_matrix dvar_fish_stock_history::get_initial_equilibrium_population1
  (dvar3_array& EN)
{
  // gets up to the nage -1 fish
  int ns=age_flags(57);
  int is;
  dvar_matrix ipop(1,nage,1,num_regions);
  ipop.initialize();
  for (int j=1;j<nage;j++)
  {
    ipop(j)=EN(j,get_season2(j));
  }
//  cout << "DD " << endl; 
//  if (num_regions>7)
//    cout <<  ipop(nage)(7,12) << endl;
//  else
//    cout <<  ipop(nage) << endl;

  return ipop;

}

dvar_vector dvar_fish_stock_history::get_initial_equilibrium_A
  (ivector& mps,dvar_matrix& ipop,imatrix& emi)
{
  // calculate the number of age nage fish when we start out with 0.
  int ns=age_flags(57);
  int is,im;
  int ir=0;
  dvar_vector ntmp1(1,num_regions);
  ntmp1.initialize();

  
  // we only need the last few cohorts
  if (!pmsd) 
  {
    dvar_matrix tmp1(1,nage,1,num_regions);
    tmp1=ipop;
    int jmin=max(1,nage-ns);
    for (is=1;is<=ns;is++)
    {
      for (int j=jmin;j<=nage;j++)
      {
        dvar_vector tmp=tmp1(j);
        for (im=1;im<=mps(is);im++)
        {
          if (!parest_flags(394))
          {
            if (age_flags(92)==0)
            {
              tmp = elem_prod(tmp,surv(is,im,j));
            }
            else
            {
              tmp = elem_prod(tmp,kludged_surv(is,im,j));
            }
          }
          else
          {
            tmp = 0.8*tmp;
          }
          if (num_regions>1)
          {
            tmp=Dad(emi(is,im),j)*tmp;
          }
        } 
        if (j<nage)
        {
          tmp1(j+1)=tmp;
        }
        else
        {
          tmp1(nage)+=tmp;
        }
      }
    }
    ntmp1=tmp1(nage);
  }
  else
  {
    int nrr=pmsd->num_real_regions;
    int numsp=pmsd->num_species;

    int jmin=max(1,nage-ns);
    for (int isp=1;isp<=numsp;isp++)
    {
      dvar_matrix tmp1(1,nage,1,nrr);
      tmp1.initialize();
      int offset=(isp-1)*nrr;
      tmp1=mycolumn(ipop,1+offset,nrr+offset).colshift(1);
      for (is=1;is<=ns;is++)
      {
        for (int j=jmin;j<=nage;j++)
        {
          dvar_vector tmp=tmp1(j);
          for (im=1;im<=mps(is);im++)
          {
            if (!parest_flags(394))
            {
              if (age_flags(92)==0)
              {
                tmp = 
                elem_prod(tmp,surv(is,im,j).sub(1+offset,nrr+offset).shift(1));
              }
              else
              {
                tmp = 
                elem_prod(tmp,
                  kludged_surv(is,im,j).sub(1+offset,nrr+offset).shift(1));
              }
            }
            else
            {
              tmp = 0.8*tmp;
            }
//            tmp = 
//              elem_prod(tmp,surv(is,im,j).sub(1+offset,nrr+offset).shift(1));
            if (num_regions>1)
            {
              tmp=Dad(emi(is,im),j).sub(1+offset,nrr+offset).shift(1)*tmp;
            }
          } //NMD_YT_jun26-17 
          if (j<nage)
          {
            tmp1(j+1)=tmp;
          }  
          else
          {
            tmp1(nage)+=tmp;
          }
        }
      }
      ntmp1(1+offset,nrr+offset).shift(1)=tmp1(nage);
    }
  }
  return ntmp1;
}

dvar_matrix dvar_fish_stock_history::get_initial_equilibrium_B
  (ivector& mps,imatrix& emi)
{
  // calculate the number of age nage fish when we start out with 0.
  int ns=age_flags(57);
  int is,im;
  int ir=0;
  dvar_vector tmp1(1,num_regions);
  dvar_matrix B;

  int pflag=0;
  if (!pmsd)
  {
    B.allocate(1,num_regions,1,num_regions);
    for (ir=1;ir<=num_regions;ir++)
    {
      tmp1.initialize();
      tmp1(ir)=1.0;
      for (is=1;is<=ns;is++)
      {
        dvar_vector tmp=tmp1;
        for (im=1;im<=mps(is);im++)
        {
          if (!parest_flags(394))
          {
            if (age_flags(92)==0)
            {
              tmp = elem_prod(tmp,surv(is,im,nage));
            }
            else
            {
              tmp = elem_prod(tmp,kludged_surv(is,im,nage));
            }
          }
          else
          {
            tmp = 0.8*tmp;
          }
          if (pflag)
          {
            //cout <<"testeq1.cpp " << surv(is,im,nage) << endl;
          }
          // now do themovement 
          if (num_regions>1)
          {
            tmp=Dad(emi(is,im),nage)*tmp;
          }
        } 
        tmp1=tmp;
      }
      B(ir)=tmp1;
    }
  }
  else
  {
    int nrr=pmsd->num_real_regions;
    int numsp=pmsd->num_species;
    B.allocate(1,num_regions,1,nrr);
    dvar_vector tmp1(1,nrr);
    for (int isp=1;isp<=numsp;isp++)
    {
      int offset=(isp-1)*nrr;
      for (ir=1;ir<=nrr;ir++)
      {
        tmp1.initialize();
        tmp1(ir)=1.0;
        for (is=1;is<=ns;is++)
        {
          dvar_vector tmp=tmp1;
          for (im=1;im<=mps(is);im++)
          {
            if (!parest_flags(394))
            {
              if (age_flags(92)==0)
              {
                tmp = 
                  elem_prod(
                  tmp,surv(is,im,nage).sub(1+offset,nrr+offset).shift(1));
              }
              else
              {
                tmp = 
                  elem_prod(tmp,
                  kludged_surv(is,im,nage).sub(1+offset,nrr+offset).shift(1));
              }
            }
            else
            {
              tmp = 0.8*tmp;
            }
//            tmp = 
//              elem_prod(
//               tmp,surv(is,im,nage).sub(1+offset,nrr+offset).shift(1));
            if (pflag)
            {
              //cout <<"testeq1.cpp " << surv(is,im,nage) << endl;
            }
            // now do themovement 
            if (num_regions>1)
            {
              tmp=Dad(emi(is,im),nage).sub(1+offset,nrr+offset).shift(1)*tmp;
            }
          } 
          tmp1=tmp;
        }
        B(ir+offset)=tmp1;
      }
    }
  }

  return trans(B);
}
dvar_matrix  dvar_fish_stock_history::get_initial_R(void)
{
  int ns=age_flags(57);
  dvar_matrix R(1,ns,1,num_regions);
  R.initialize();
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int iy=1;iy<=ns;iy++)
    {
      R(iy,ir)=exp(N(ir,iy,1));
    }
  }
  return R;
}

dvar_matrix  dvar_fish_stock_history::get_initial_equilibrium_recruitment
  (int _nya)
{
  int nya=_nya;
  if (nya==0)
  {
    cerr << "putting in kludge for nya=0 setting it equal to 1" << endl;
    nya=max(2,age_flags(57)+1);
  }
  
  dvar_matrix rec(1,age_flags(57),1,num_regions);
  imatrix ny(1,num_regions,1,age_flags(57));
  rec.initialize();
  ny.initialize();
  int ir;
  dvar3_array * pN=0;
  if (af170q0==0)
  {
    pN=&N;
  }
  else
  {
    pN=&N_q0;
  }
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int iy=2;iy<=nya;iy++)
    {
      int is=(iy-1)%age_flags(57)+1;
      rec(is,ir)+=exp((*pN)(ir,iy,1));  
      ny(ir,is)++;
    }
  }
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int is=1;is<=age_flags(57);is++)
    {
      rec(is,ir)/=double(ny(ir,is));  
    }
  }
  return rec;
}


ivector  dvar_fish_stock_history::get_equilibrium_movements_per_season(void)
{
  int i,j;
  // this is the number of years to use in the average

  // get the number of movement periods in each season;
  // number of seasons is suposed to be equal to the number of recruitments
  // i.e.  age_flags(57)
  ivector mps(1,age_flags(57));  //NMD_13dec2021
  if (num_regions>1)  //NMD_13dec2021
  {
    int ir;
    imatrix moves_per_year(1,num_regions,1,last_real_year);
    moves_per_year.initialize();
    for (ir=1;ir<=num_regions;ir++)
    {
      int oldyr=year(ir,1);
      for (int i=1;i<=oldyr;i++)
      {
        moves_per_year(ir,i)=1;
      }
      for (int ip=2;ip<=num_real_fish_periods(ir);ip++)  
      {
        int yr=year(ir,ip);
        if (yr>oldyr && move_flags(ir,ip)==0)
        {
          cerr << "first fishing period of year is not a movement period"
               << endl;
//        cerr << "enter a character to continue " << endl;
//        char ch;
//        cin >> ch;
        }
        oldyr=year(ir,ip);
        if (move_flags(ir,ip))
        {
          moves_per_year(ir,yr)++;
        }
      }
    }
    check_sanity_trans(moves_per_year);

//  ivector mps(1,age_flags(57));
    mps=moves_per_year(1)(1,age_flags(57));
  }
  else  //NMD_13dec2021
  {
    mps=1;
  }

  return mps;
}

ivector  dvar_fish_stock_history::kludged_get_equilibrium_movements_per_season(void)
{
  int i,j;
  // this is the number of years to use in the average

  // get the number of movement periods in each season;
  // number of seasons is suposed to be equal to the number of recruitments
  // i.e.  age_flags(57)
  ivector mps(1,age_flags(57));  //NMD_13dec2021
  if (num_regions>1)  //NMD_13dec2021
  {
    int ir;
    imatrix moves_per_year(1,num_regions,1,1000);
    moves_per_year.initialize();
    moves_per_year(1,1)=1;
    for (ir=1;ir<=num_regions;ir++)
    {
      int oldyr=year(ir,1);
      for (int i=1;i<=oldyr;i++)
      {
        moves_per_year(ir,i)=1;
      }
      for (int ip=2;ip<=num_real_fish_periods(ir);ip++)  
      {
        int yr=year(ir,ip);
        if (yr>oldyr && move_flags(ir,ip)==0)  //NMD_3dec2021
        {
          cerr << "first fishing period of year is not a movement period"
               << endl;
//        cerr << "enter a character to continue " << endl;
//        char ch;
//        cin >> ch;
        }
        oldyr=year(ir,ip);
        if (move_flags(ir,ip))
        {
          moves_per_year(ir,yr)++;
        }
      }
    }
    check_sanity_trans(moves_per_year);

//  ivector mps(1,age_flags(57)); //NMD_13dec2021
    mps=moves_per_year(1)(1,age_flags(57));
  }
  else  //NMD_13dec2021
  {
    mps=1;
  }
  
  return mps;
}


imatrix  dvar_fish_stock_history::get_initial_equilibrium_movement(ivector& mps)
{

  int ir;
  imatrix movement(1,age_flags(57),1,mps);
  movement(1,1)=1;

  ivector moves_per_year(1,age_flags(57));
  moves_per_year.initialize();
  int oldyr=year(1,1);
  for (int ip=1;ip<=num_real_fish_periods(1);ip++)  
  {
    int yr=year(1,ip);
    if (yr>age_flags(57)) break;
    if (yr>oldyr && move_flags(1,ip)==0)
    {
//      cerr << "first fishing period of year is not a movement period"
//           << endl;
//      cerr << "enter a character to continue " << endl;
//      char ch;
//      cin >> ch;
    }
    oldyr=year(1,ip);  //NMD 28Mar2012
    if (move_flags(1,ip))
    {
      moves_per_year(yr)++;
      movement(yr,moves_per_year(yr))=move_index(1,ip);
      if (parest_flags(357)>0)
//      if (fish_nonmv_flag>0)  //NMD_14Sep2018
      {
        movement(yr,moves_per_year(yr))= 
          (movement(yr,moves_per_year(yr)))%age_flags(57)+1;
      }
    }
  }
  return movement;
}


dvar4_array  dvar_fish_stock_history::get_initial_equilibrium_survival
  (ivector & mps,int nya,ivector * pq_flag)
{
  // Now we need to determine the total morality in each
  // between movement period

  int af57=age_flags(57);
  dvar4_array tmort(1,af57,
    1,mps,1,num_regions,1,nage);  
//  dvar4_array surv(1,af57,
//    1,mps,1,nage,1,num_regions);  
  dvar4_array local_surv(1,af57,
    1,mps,1,nage,1,num_regions);  
  if (!allocated(csurv))
  {
    csurv.allocate(1,af57,1,mps,1,nage,1,num_regions);  
    csurv_chk.allocate(0,10);
    for (int i=0;i<=10;i++)
    {
      csurv_chk(i).allocate(1,af57,1,mps,1,nage,1,num_regions);  
    }
    csurv_chk.initialize();
  }

  tmort.initialize();
  local_surv.initialize();  //NMD_9May2022
  i3_array ny(1,num_regions,1,age_flags(57),1,imatrix(1,num_regions,mps));
  ny.initialize();

  int ir;
  int newyear_flag=0;
  for (ir=1;ir<=num_regions;ir++)
  {
    int oldyr=year(ir,1);
    int mpp=1;
    int is=(oldyr-1)%age_flags(57)+1;
    switch(age_flags(94))
    {
    case 1:
//      if (age_flags(95)==0)
//      if (age_flags(128)==0)  //NMD_21jun2016
      if (age_flags(128)==0 || pq_flag)  //NMD_6jun2022
      {
        tmort(is,mpp,ir)+=mfexp(nat_mort(year(ir,1))+fraction(ir,1));
      }
      else
      {
//        tmort(is,mpp,ir)+=age_flags(95)/10.*
        tmort(is,mpp,ir)+=age_flags(128)/100.*  //NMD_21jun2016
          mfexp(get_nat_mort_region(ir)(year(ir,1))+fraction(ir,1));
      }
      break;
    case 2:
      check_code();
      /*
      if (parest_flags(373) && parest_flags(374))
      {
        cerr << "Initial population option for average totmort without"
             << "initial coefficients estimation is not possible for "
             << "the catch conditioned case - quitting" << endl;
        ad_exit(1);
      }
      */
      if (pq_flag)
      {
        tmort(is,mpp,ir)+=tot_mort_q0(ir,1);
      }
      else
      {
        tmort(is,mpp,ir)+=tot_mort(ir,1);
      }
      break;
    case 3:
        tmort(is,mpp,ir)+=mfexp(nat_mort(year(ir,1))+fraction(ir,1));
    break;
    default:
//      cerr << "Illegal value for age_flags(95) = " << age_flags(95)
      cerr << "Illegal value for age_flags(128) = " << age_flags(128)
// NMD_21jun2016 
           << endl;
      ad_exit(1);
    }

    for (int ip=2;ip<=num_real_fish_periods(ir);ip++)  
    {
      int yr=year(ir,ip);
      if (yr>nya) break;
      int is=(yr-1)%age_flags(57)+1;
//      if (yr>1)  //NMD_3May2022 - causes error in mpp for year=1 & mps>1
//      {
        if (yr>oldyr)  // new year
        {
          newyear_flag=1;
          oldyr=yr;
          mpp=1;
        }
        else if (move_flags(ir,ip))
        {
          mpp++;
        }
//      }

      switch(age_flags(94))
      {
      case 1:
//        if (age_flags(95)==0) // NMD_21jun2016 
//        if (age_flags(128)==0)
        if (age_flags(128)==0 || pq_flag)  //NMD_6jun2022
          tmort(is,mpp,ir)+=mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
        else
          tmort(is,mpp,ir)+=
//            age_flags(95)/10.*mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
            age_flags(128)/100.*mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
// NMD_21jun2016 
        break;
      case 2:
      check_code();
         /*
        if (parest_flags(373) && parest_flags(374))
        {
          cerr << "Initial population option for average totmort without"
               << "initial coefficients estimation is not possible for "
               << "the catch conditioned case - quitting" << endl;
          ad_exit(1);
        }
        */
        if (pq_flag)
        {
          tmort(is,mpp,ir)+=tot_mort_q0(ir,1);
        }
        else
        {
          tmort(is,mpp,ir)+=tot_mort(ir,ip);
        }
        break;
      case 3:
          tmort(is,mpp,ir)+=mfexp(nat_mort(year(ir,1))+fraction(ir,1));
        break;
      default:
//        cerr << "Illegal value for gage_flags(95) = " << age_flags(95)
        cerr << "Illegal value for gage_flags(128) = " << age_flags(128)
//NMD_21jun2016
             << endl;
        ad_exit(1);
      }
    }
  }
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int is=1;is<=age_flags(57);is++)
    {
      for (int pi=1;pi<=mps(is);pi++)
      {
        //tmort(is,pi,ir)/=ny(ir,is,pi);
        tmort(is,pi,ir)/=(nya/double(month_factor));
      }
    }
  }
  for (int is=1;is<=age_flags(57);is++)
  {
    for (int pi=1;pi<=mps(is);pi++)
    {
      if (pq_flag)
      {
//        surv(is,pi)=mfexp(-trans(tmort(is,pi))-0.0001);
        local_surv(is,pi)=mfexp(-trans(tmort(is,pi))-0.0001);
      }
      else
      {
//        surv(is,pi)=mfexp(-trans(tmort(is,pi))-0.0001);
        local_surv(is,pi)=mfexp(-trans(tmort(is,pi))-0.0001);
      }
//      csurv(is,pi)=value(surv(is,pi));
//      csurv_chk(iloop,is,pi)=value(surv(is,pi));
      csurv(is,pi)=value(local_surv(is,pi));
      csurv_chk(iloop,is,pi)=value(local_surv(is,pi));
    }
  }
//  return surv;
  return local_surv;
}


dvar4_array  dvar_fish_stock_history::get_initial_equilibrium_survival_ms
  (ivector & mps,int nya,ivector * pq_flag)
{
  // Now we need to determine the total morality in each
  // between movement period

  int af57=age_flags(57);
  dvar4_array tmort(1,af57,
    1,mps,1,num_regions,1,nage);  
  dvar4_array surv(1,af57,
    1,mps,1,nage,1,num_regions);  
  if (!allocated(csurv))
  {
    csurv.allocate(1,af57,1,mps,1,nage,1,num_regions);  
    csurv_chk.allocate(0,10);
    for (int i=0;i<=10;i++)
    {
      csurv_chk(i).allocate(1,af57,1,mps,1,nage,1,num_regions);  
    }
    csurv_chk.initialize();
  }
  tmort.initialize();  
  i3_array ny(1,num_regions,1,af57,1,imatrix(1,num_regions,mps));
  ny.initialize();

  int ir;
  int newyear_flag=0;
  for (ir=1;ir<=num_regions;ir++)
  {
    int af94=get_age_flags_region(ir,94);
    int af95=get_age_flags_region(ir,95);
    int oldyr=year(ir,1);
    int mpp=1;
    int is=(oldyr-1)%af57+1;
    switch(af94)
    {
    case 1:
      if (af95==0)
      {
        tmort(is,mpp,ir)+=mfexp(get_nat_mort_region(ir)(year(ir,1))
          +fraction(ir,1));
      }
      else
      {
        tmort(is,mpp,ir)+=af95/10.*
          mfexp(get_nat_mort_region(ir)(year(ir,1))+fraction(ir,1));
      }
      break;
    case 2:
      if (pq_flag)
      {
        tmort(is,mpp,ir)+=tot_mort_q0(ir,1);
      }
      else
      {
        tmort(is,mpp,ir)+=tot_mort(ir,1);
      }
      break;
    default:
      cerr << "Illegal value for gaf95 = " << af95
           << endl;
      ad_exit(1);
    }

    for (int ip=2;ip<=num_real_fish_periods(ir);ip++)  
    {
      int yr=year(ir,ip);
      if (yr>nya) break;
      int is=(yr-1)%af57+1;
//      if (yr>1)  //NMD_3May2022 - causes error in mpp for year=1 & mps>1
//      {
        if (yr>oldyr)  // new year
        {
          newyear_flag=1;
          oldyr=yr;
          mpp=1;
        }
        else if (move_flags(ir,ip))
        {
          mpp++;
        }
//      }

      switch(af94)
      {
      case 1:
        if (af95==0)
          tmort(is,mpp,ir)+=mfexp(get_nat_mort_region(ir)(year(ir,ip))+fraction(ir,ip));
        else
          tmort(is,mpp,ir)+=af95/10.*
            mfexp(get_nat_mort_region(ir)(year(ir,ip))+fraction(ir,ip));
        break;
      case 2:
         if (pq_flag)
         {
           tmort(is,mpp,ir)+=tot_mort_q0(ir,1);
         }
         else
         {
           tmort(is,mpp,ir)+=tot_mort(ir,ip);
         }
         break;
      default:
        cerr << "Illegal value for gaf95 = " << af95
             << endl;
        ad_exit(1);
      }
    }
  }
  for (ir=1;ir<=num_regions;ir++)
  {
    int af94=get_age_flags_region(ir,94);
    int af95=get_age_flags_region(ir,95);
    for (int is=1;is<=af57;is++)
    {
      for (int pi=1;pi<=mps(is);pi++)
      {
        //tmort(is,pi,ir)/=ny(ir,is,pi);
        tmort(is,pi,ir)/=(nya/double(month_factor));
      }
    }
  }
  for (int is=1;is<=af57;is++)
  {
    for (int pi=1;pi<=mps(is);pi++)
    {
      if (pq_flag)
      {
        surv(is,pi)=mfexp(-trans(tmort(is,pi))-0.0001);
      }
      else
      {
        surv(is,pi)=mfexp(-trans(tmort(is,pi))-0.0001);
      }
      csurv(is,pi)=value(surv(is,pi));
      csurv_chk(iloop,is,pi)=value(surv(is,pi));
    }
  }
  {
    ofstream ofs("SURV");
    ofs << setw(10) << surv << endl;
  }
//  {
//    ofstream ofs("TOTMS");
//    for (int is=7;is<=12;is++)
//    {
//      ofs << "region " << is << endl;
//      ofs << setw(10) << tot_mort(is) << endl;
//    }
//  }
  return surv;
}


void check_sanity(imatrix& m)
{
  int cmin=m.indexmin();
  int cmax=m.indexmax();

  for (int i=cmin;i<=cmax;i++)
  {
    int rmin = m(i).indexmin();
    int rmax = m(i).indexmax();

    int m1=m(i,rmin);
    for (int j=rmin+1;j<=rmax;j++)
    {
      if (m1 != m(i,j))
      {
        cerr << "Sanity error" << endl;
        ad_exit(1);
      }
    }
  }
}


void check_sanity_trans(imatrix& _m)
{
  imatrix m=trans(_m);
  int cmin=m.indexmin();
  int cmax=m.indexmax();

  for (int i=cmin;i<=cmax;i++)
  {
    int rmin = m(i).indexmin();
    int rmax = m(i).indexmax();

    int m1=m(i,rmin);
    for (int j=rmin+1;j<=rmax;j++)
    {
      if (m1 != m(i,j))
      {
        cerr << "Sanity error" << endl;
        ad_exit(1);
      }
    }
  }
}



int dvar_fish_stock_history::get_season2(int j)
{
  int ns=age_flags(57);
  int i=(ns-(j-1))%ns;
  if(i<0) i+=ns;
  return (i+1);
}

int dvar_fish_stock_history::get_season(int i,int j)
{
  //int n_s=age_flags(57);
  return ( (i+j-2) % age_flags(57) )  + 1;
}
