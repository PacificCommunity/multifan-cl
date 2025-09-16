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
  tail_compress_predicted_length_frequencies(void)
{
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)
      {
        if (len_sample_size(ir,ip,fi)>0)
        {
          int mmin=tc_len_freq(ir,ip,fi).indexmin();
          int mmax=tc_len_freq(ir,ip,fi).indexmax();
          int nlint1=mmax-mmin+1;
          dvar_vector & tp =tprob(ir,ip,fi);
          int mintp=tp.indexmin();
          int maxtp=tp.indexmax();
          if (allocated(tc_tprob(ir,ip,fi)))
            tc_tprob(ir,ip,fi).deallocate();
          tc_tprob(ir,ip,fi).allocate(mmin,mmax);
          dvar_vector & tc=tc_tprob(ir,ip,fi);
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

void dvar_len_fish_stock_history::make_tail_compressed_samples(void)
{
  ofstream * pofs=0;
  if (parest_flags(313))
    pofs= new ofstream("lfreqs");
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)  
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {     
        if (len_sample_size(ir,ip,fi))
        {
          dvector & lf=len_freq(ir,ip,fi);
          int mmin=lf.indexmin();
          int mmax=lf.indexmax();
          int i=0;
          int j=0;
          // for now just use one flag for lef and right bounds
          MY_DOUBLE_TYPE leftbound=parest_flags(313)/100.;
          MY_DOUBLE_TYPE rightbound=parest_flags(313)/100.;
          MY_DOUBLE_TYPE ssum=0.0;
          ssum=0.0;
          for (i=mmin;i<=mmax;i++)
          {
            ssum+=lf(i);
            if (ssum>leftbound)
             break;
          }
          ssum=0.0;
          for (j=mmax;j>=mmin;j--)
          {
            ssum+=lf(j);
            if (ssum>rightbound)
             break;
          }
          if (i>mmax || j<mmin)
          {
            cerr << "make_tail_compressed_samples - implausible range" << endl;
            ad_exit(1);
          }
          if (i==j)
          {
            cout << "Only one non zero slot in this sample"
                 << endl << "setting length sample size to zero" << endl;
            len_sample_size(ir,ip,fi)=0.0;
            max_len_obs(ir,ip,fi)=i;
          }
          else if (parest_flags(312)>0 &&
            len_sample_size(ir,ip,fi)< parest_flags(312))
          {
            cout << "sample size of " <<  len_sample_size(ir,ip,fi)
             << " is less than minimum of " <<  parest_flags(312)
                 << endl << "setting length sample size to zero" << endl;
            len_sample_size(ir,ip,fi)=0.0;
            max_len_obs(ir,ip,fi)=mmax;
          }
          else
          {
            tc_len_freq(ir,ip,fi).allocate(i,j);
            tc_len_freq(ir,ip,fi)=lf(i,j);
            if (i>mmin)
              tc_len_freq(ir,ip,fi)(i)+=sum(lf(mmin,i-1));
            if (j<mmax)
              tc_len_freq(ir,ip,fi)(j)+=sum(lf(j+1,mmax));
            if (pofs)
              (*pofs) << tc_len_freq(ir,ip,fi) << endl;
             
            MY_DOUBLE_TYPE tmp=lf(i);
            int imax=i;
            for (int ii=i+1;ii<=j;ii++)
            {
              if (lf(ii)>tmp)
              {
                tmp=lf(ii);
                imax=ii;
              }
            }
            max_len_obs(ir,ip,fi)=imax;
          }
        }
      }
    }
  }
  if (pofs)
    delete pofs;
  pofs=0;
}

dvariable tail_compressed_logistic_normal_fit_heteroscedastic
  (dvar_len_fish_stock_history& fsh,d3_array& total_freq)
{
  dvariable rho=fsh.length_rho;
  dvariable psi=fsh.length_psi;
  dvariable var=exp(fsh.log_length_variance);
  int nlint=fsh.nlint;
  
  dvar_matrix S(1,nlint,1,nlint);
  dvar_matrix CS(1,nlint,1,nlint);

  dvar_vector rp(0,nlint-1);
  rp(0)=1.0;
  switch (fsh.parest_flags(297))
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
    cerr << "illegal value for fsh.parest_flags(297)" << endl;
    ad_exit(1);
  }
  switch (fsh.parest_flags(297))
  {
  case 0:
  case 1:   // LN3m
    for (int i=2;i<nlint;i++)
    {
      rp(i)=rho*rp(i-1);
    }
    break;
  case 2:   // LN3
    {
      dvariable phi2=-1+(2.0-sfabs(rho))*psi;
      for (int i=2;i<nlint;i++)
      {
        rp(i)=rho*rp(i-1)+phi2*rp(i-2);
      }
    }
    break;
  default:
    cerr << "illegal value for fsh.parest_flags(297)" << endl;
    ad_exit(1);
  }

  for (int i=1;i<=nlint;i++)
  {
    S(i,i)=1.0;
    for (int j=1;j<i;j++)
    {
      S(i,j)=rp(i-j);
      S(j,i)=rp(i-j);
    }
  }
  CS=value(S);
  dvariable len_tmp=0.;
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
          if (fsh.len_sample_size(ir,ip,fi)>0)
          {
            int mmin=fsh.tc_len_freq(ir,ip,fi).indexmin();
            int mmax=fsh.tc_len_freq(ir,ip,fi).indexmax();
            int nlint1=mmax-mmin+1;
            // lmax will be the index of the slot for the 
            // largest proportion
            int lmax=static_cast<int>(fsh.max_len_obs(ir,ip,fi));
            lmax=mmax;
            MY_DOUBLE_TYPE epsilon=1.e-6;
            if (fsh.parest_flags(294))
              epsilon=fsh.parest_flags(294)/1.e+8;
              
            dvector ole=fsh.tc_len_freq(ir,ip,fi)+epsilon;
            ole/=sum(ole);
            dvar_vector & tp =fsh.tc_tprob(ir,ip,fi);
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
            if (!fsh.parest_flags(314))
            {
              vp=pow(elem_prod(pole,(1.0-pole)),fsh.length_tot_exp);
            }
            else
            {
              vp=pow(pole,fsh.length_tot_exp);
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
               e*=pow(fsh.len_sample_size(ir,ip,fi),fsh.length_exp);
            }
  
            int pdftype=fsh.parest_flags(293);
            switch(pdftype)
            {
            case 0: // normal
              {
                const MY_DOUBLE_TYPE ltpi=0.5*log(2.0*3.1415926535);
                len_tmp+=(nlint1-1)*ltpi + 0.5*(e*e) + ln_det_choleski;
                if (fsh.parest_flags(295))
                  len_tmp-=(nlint1-1)*fsh.length_exp
                    *log(fsh.len_sample_size(ir,ip,fi));
                break;
              }
            case 1: // student
              {
                dvariable v=exp(fsh.log_length_dof);  // this will be 
                                            // the degrees of freedoom for
                MY_DOUBLE_TYPE p=nlint1-1.0;
                const MY_DOUBLE_TYPE lppi2=0.5*p*log(3.1415926535);
                loc_tmp+= -gammln(0.5*(v+p)) + gammln(0.5*v) 
                        + lppi2 +(0.5*p)*log(v)
                        +0.5*(v+p)*log(1.0+e*e/v) 
                        +ln_det_choleski;
                 if ((e.indexmax()-e.indexmin()+1) != p)
                 {
                   cerr << "Logistic-normal heteroscedastic PDF error " << endl;
                   cerr << "integral of student-t != proportion" << endl;
                   ad_exit(1);
                 }
  
                if (fsh.parest_flags(295))
                  loc_tmp-=p*fsh.length_exp*log(fsh.len_sample_size(ir,ip,fi));
                break;
              }
            default: // error
              cerr << "illegal value for parest_flags(293)" << endl;
              ad_exit(1);
            }
          }
        }
      }
      ftmp=loc_tmp;
      end_tmp=ftmp;
      len_tmp+=end_tmp;
    }
  }
  return len_tmp;
}

