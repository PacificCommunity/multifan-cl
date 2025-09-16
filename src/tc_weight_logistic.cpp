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


void dvar_len_fish_stock_history::
  tail_compress_predicted_weight_frequencies(void)
{
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)
      {
        if (wght_sample_size(ir,ip,fi)>0)
        {
          int mmin=tc_wght_freq(ir,ip,fi).indexmin();
          int mmax=tc_wght_freq(ir,ip,fi).indexmax();
          int nlint1=mmax-mmin+1;
          dvar_vector & tp =wtprob(ir,ip,fi);
          int mintp=tp.indexmin();
          int maxtp=tp.indexmax();
          if (allocated(tc_wtprob(ir,ip,fi)))
            tc_wtprob(ir,ip,fi).deallocate();
          tc_wtprob(ir,ip,fi).allocate(mmin,mmax);
          dvar_vector & tc=tc_wtprob(ir,ip,fi);
          tc=tp(mmin,mmax);
          if (mintp<mmin)
          {
            tc(mmin)+=sum(tp(mintp,mmin-1));
          }
          if (maxtp>mmax)
          {
            tc(mmax)+=sum(tp(mmax+1,maxtp));
          }
        }
      }
    }
  }
}

void dvar_len_fish_stock_history::make_tail_compressed_weight_samples(void)
{
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)  
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {     
        if (wght_sample_size(ir,ip,fi))
        {
          dvector & wf=wght_freq(ir,ip,fi);
          int mmin=wf.indexmin();
          int mmax=wf.indexmax();
          int i=0;
          int j=0;
          // for now just use one flag for lef and right bounds
          MY_DOUBLE_TYPE leftbound=parest_flags(303)/100.;
          MY_DOUBLE_TYPE rightbound=parest_flags(303)/100.;
          MY_DOUBLE_TYPE ssum=0.0;
          ssum=0.0;
          for (i=mmin;i<=mmax;i++)
          {
            ssum+=wf(i);
            if (ssum>leftbound)
             break;
          }
          ssum=0.0;
          for (j=mmax;j>=mmin;j--)
          {
            ssum+=wf(j);
            if (ssum>rightbound)
             break;
          }

          if (i>mmax || j<mmin)
          {
            cerr << "this cant happen" << endl;
            ad_exit(1);
          }
          if (i==j)
          {
            cout << "Only one non zero slot in this sample"
                 << endl << "setting weight sample size to zero" << endl;
            wght_sample_size(ir,ip,fi)=0.0;
            max_wght_obs(ir,ip,fi)=i;
          }
          else if (parest_flags(302)>0 && 
            wght_sample_size(ir,ip,fi)< parest_flags(302))
          {
            cout << "sample size of " <<  wght_sample_size(ir,ip,fi)
             << " is less than minimum of " <<  parest_flags(302)
                 << endl << "setting weight sample size to zero" << endl;
            wght_sample_size(ir,ip,fi)=0.0;
            max_wght_obs(ir,ip,fi)=mmax;
          }
          else
          {
            tc_wght_freq(ir,ip,fi).allocate(i,j);
            tc_wght_freq(ir,ip,fi)=wf(i,j);
            if (i>mmin)
              tc_wght_freq(ir,ip,fi)(i)+=sum(wf(mmin,i-1));
            if (j<mmax)
              tc_wght_freq(ir,ip,fi)(j)+=sum(wf(j+1,mmax));
            MY_DOUBLE_TYPE tmp=wf(i);
            int imax=i;
            for (int ii=i+1;ii<=j;ii++)
            {
              if (wf(ii)>tmp)
              {
                tmp=wf(ii);
                imax=ii;
              }
            }
            max_wght_obs(ir,ip,fi)=imax;
          }
        }
      }
    }
  }
}


dvariable tail_compressed_weight_logistic_normal_fit_heteroscedastic
  (dvar_len_fish_stock_history& fsh,d3_array& total_freq)
{
  dvariable rho=fsh.weight_rho;
  dvariable psi=fsh.weight_psi;
  dvariable var=exp(fsh.log_weight_variance);
  int nwint=fsh.nwint;
  
  dvar_matrix S(1,nwint,1,nwint);

  dvar_vector rp(0,nwint-1);
  rp(0)=1.0;
  switch (fsh.parest_flags(287))
  {
  case 0:
    rp(1)=rho;  //LN1  LN2
    break;
  case 1:   // LN3m
    rp(1)=rho+psi/(1.0+square(rho+psi)/(1.0-square(rho)));
    break;
  case 2:   // LN3
    {
      dvariable phi2=-1+(2.0-sfabs(rho))*psi;
      rp(1)=rho/(1.0-phi2);
    }
    break;
  default:
    cerr << "illegal value for fsh.parest_flags(287)" << endl;
    ad_exit(1);
  }
  switch (fsh.parest_flags(287))
  {
  case 0:
  case 1:   // LN3m
    for (int i=2;i<nwint;i++)
    {
      rp(i)=rho*rp(i-1);
    }
    break;
  case 2:   // LN3
    {
      dvariable phi2=-1+(2.0-sfabs(rho))*psi;
      for (int i=2;i<nwint;i++)
      {
        rp(i)=rho*rp(i-1)+phi2*rp(i-2);
      }
    }
    break;
  default:
    cerr << "illegal value for fsh.parest_flags(287)" << endl;
    ad_exit(1);
  }

  for (int i=1;i<=nwint;i++)
  {
    S(i,i)=1.0;
    for (int j=1;j<i;j++)
    {
      S(i,j)=rp(i-j);
      S(j,i)=rp(i-j);
    }
  }
  dvariable wght_tmp=0.;
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    int ntimes;
    if (fsh.parest_flags(142)==0)
    {
      ntimes=fsh.num_fish_periods(ir);
    }
    else 
    {
      ntimes=min(fsh.parest_flags(142),fsh.num_fish_periods(ir));
    }
    funnel_dvariable ftmp;
    //dvariable ftmp;
    dvariable end_tmp=0.0;
    int blocksize=40;
    //for (int ipp=1;ipp<=10;ipp+=blocksize)
    for (int ipp=1;ipp<=ntimes;ipp+=blocksize)
    {
      dvariable loc_tmp=0.0;
      ad_begin_funnel();
      for (int ip=ipp;ip<=ipp+blocksize-1;ip++)
      {
        if (ip>ntimes) break;
        for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
        {
          if (fsh.wght_sample_size(ir,ip,fi)>0)
          {
            int mmin=fsh.tc_wght_freq(ir,ip,fi).indexmin();
            int mmax=fsh.tc_wght_freq(ir,ip,fi).indexmax();
            int nwint1=mmax-mmin+1;
            // lmax will be the index of the slot for the 
            // largest proportion
            int lmax=static_cast<int>(fsh.max_wght_obs(ir,ip,fi));
            lmax=mmax;
            MY_DOUBLE_TYPE epsilon=1.e-6;
            if (fsh.parest_flags(284))
              epsilon=fsh.parest_flags(284)/1.e+8;
              
            dvector ole=fsh.tc_wght_freq(ir,ip,fi)+epsilon;
            ole/=sum(ole);

            dvar_vector & tp =fsh.tc_wtprob(ir,ip,fi);
            dvar_vector pole=(tp+epsilon);
           /*
            dvar_vector pole=(tp+epsilon)(mmin,mmax);
            int maxtp=tp.indexmax();
            int mintp=tp.indexmin();
            if (mintp<mmin)
            {
              pole(mmin)+=sum(tp(mintp,mmin-1));
            }
            if (maxtp>mmax)
            {
              pole(mmax)+=sum(tp(mmax+1,maxtp));
            }
            pole/=sum(pole);
           */
            dvar_vector lole1=log(ole);
            dvar_vector lpole1=log(pole);
            dvar_vector diff1=lole1-lpole1;
            
            dvar_vector vp;
            if (!fsh.parest_flags(304))
            {
              vp=pow(elem_prod(pole,(1.0-pole)),fsh.weight_tot_exp);
            }
            else
            {
              vp=pow(pole,fsh.weight_tot_exp);
            }

            dvar_matrix SUB(mmin,mmax-1,mmin,mmax-1);
            for (int i=mmin;i<=mmax-1;i++)
            {
              SUB(i,i)=vp(i)*vp(i)*S(i,i)+vp(mmax)*vp(mmax)*S(mmax,mmax)
                -2.0*vp(i)*vp(mmax)*S(i,mmax);
              for (int j=mmin;j<i;j++)
              {
                SUB(i,j)=vp(i)*vp(j)*S(i,j)+vp(mmax)*vp(mmax)*S(mmax,mmax)
                  -vp(i)*vp(mmax)*S(i,mmax)-vp(j)*vp(mmax)*S(j,mmax);
                SUB(j,i)=SUB(i,j);
              }
            }
            SUB*=var;
            dvar_vector diff(mmin,mmax-1);
            diff=diff1(mmin,mmax-1)-diff1(mmax);
           
            dvariable ldet=0.0;
            int sgn=1;
            dvariable sgn1=0.0;
            {
              diff.shift(1);
              for (int i=mmin;i<=mmax-1;i++)
              {
                SUB(i).shift(1);
              }
              SUB.shift(1);
            }
            dvar_vector e=choleski_factor_solve(SUB,diff,ldet,sgn);
            dvariable ln_det_choleski=ldet;

            if (fsh.parest_flags(295))
            {
               e*=pow(fsh.wght_sample_size(ir,ip,fi),fsh.weight_exp);
            }
  
            int pdftype=fsh.parest_flags(283);
            switch(pdftype)
            {
            case 0: // normal
              {
                const MY_DOUBLE_TYPE ltpi=0.5*log(2.0*3.1415926535);
                wght_tmp+=(nwint1-1)*ltpi + 0.5*(e*e) + ln_det_choleski;
                if (fsh.parest_flags(295))
                  wght_tmp-=(nwint1-1)*fsh.weight_exp
                    *log(fsh.wght_sample_size(ir,ip,fi));
                break;
              }
            case 1: // student
              {
                //dvariable v=50.0;  // this will be the degrees of freedoom for
                                 // the students t
                dvariable v=exp(fsh.log_weight_dof);  // this will be 
                                            // the degrees of freedoom for
                MY_DOUBLE_TYPE p=nwint1-1.0;
                const MY_DOUBLE_TYPE lppi2=0.5*p*log(3.1415926535);
                loc_tmp+= -gammln(0.5*(v+p)) + gammln(0.5*v) 
                        + lppi2 +(0.5*p)*log(v)
                        +0.5*(v+p)*log(1.0+e*e/v) 
                        +ln_det_choleski;
  
                //cout << "EEE loc_tmp " << loc_tmp << endl;
                if (fsh.parest_flags(295))
                  loc_tmp-=p*fsh.weight_exp*log(fsh.wght_sample_size(ir,ip,fi));
                break;
              }
            default: // error
              cerr << "illegal value for parest_flags(283)" << endl;
              ad_exit(1);
            }
          }
        }
      }
      ftmp=loc_tmp;
      end_tmp=ftmp;
      wght_tmp+=end_tmp;
    }
  }
  return wght_tmp;
}

