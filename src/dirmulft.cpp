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

// dirichlet multinomail mixture

dvariable dirichlet_multinomial_fit(dvar_len_fish_stock_history& fsh, 
  d3_array& total_freq)
{
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

    dvar_vector A=mfexp(fsh.fish_pars(12));
    for (int ip=1;ip<=ntimes;ip++)
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
      {
        if (fsh.len_sample_size(ir,ip,fi)>0)
        {
          int fishery=fsh.parent(ir,ip,fi);

          dvar_vector alpha=
            A(fishery)*(0.05/fsh.nlint+fsh.tprob(ir,ip,fi));

          len_tmp-=gammln(A(fishery))
           + sum(gammln(alpha+fsh.len_freq(ir,ip,fi)))
           -gammln(A(fishery)+fsh.len_sample_size(ir,ip,fi))
           -sum(gammln(alpha));

        }
      }
    }
  }
  return len_tmp;
}

dvariable dirichlet_multinomial_mixture_fit(dvar_len_fish_stock_history& fsh, 
  d3_array& total_freq)
{
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

    dvar_vector A=mfexp(fsh.fish_pars(12));
    MY_DOUBLE_TYPE q=0.05;
    for (int ip=1;ip<=ntimes;ip++)
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
      {
        if (fsh.len_sample_size(ir,ip,fi)>0)
        {
          int fishery=fsh.parent(ir,ip,fi);

          dvar_vector alpha=
            A(fishery)*(0.05/fsh.nlint+fsh.tprob(ir,ip,fi));

          dvariable tmp=gammln(A(fishery))
           -gammln(A(fishery)+fsh.len_sample_size(ir,ip,fi))
           + sum(gammln(alpha+fsh.len_freq(ir,ip,fi)))
           -sum(gammln(alpha));

          dvariable A1=A(fishery)/5.0;
          dvar_vector alpha1=
            A1*(0.05/fsh.nlint+fsh.tprob(ir,ip,fi));

          dvariable tmp1=gammln(A1)
           -gammln(A1+fsh.len_sample_size(ir,ip,fi))
           + sum(gammln(alpha1+fsh.len_freq(ir,ip,fi)))
           -sum(gammln(alpha1));

          dvariable xtmp;
          if (tmp>tmp1)
            xtmp=tmp+log(q*exp(tmp1-tmp)+(1.0-q));
          else
            xtmp=tmp1+log(q+(1.0-q)*exp(tmp-tmp1));
          xtmp+=gammln(1.0+fsh.len_sample_size(ir,ip,fi))
                -sum(gammln(1.0+fsh.len_freq(ir,ip,fi)));
          len_tmp-=xtmp;

        }
      }
    }
  }
  return len_tmp;
}

dvariable multinomial_fit(dvar_len_fish_stock_history& fsh, 
  d3_array& total_freq)
{
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

    for (int ip=1;ip<=ntimes;ip++)
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
      {
        if (fsh.len_sample_size(ir,ip,fi)>0)
        {
          len_tmp-=fsh.len_sample_size(ir,ip,fi)*
            (fsh.len_freq(ir,ip,fi)*log(1.e-10+fsh.tprob(ir,ip,fi)));
        }
      }
    }
  }
  return len_tmp;
}

dvariable logistic_normal_fit(dvar_len_fish_stock_history& fsh, 
  d3_array& total_freq)
{
  // *********************************************************
  // *********************************************************
  // need to calculate these. for the moment they are just
  // here so that this mess will compile.
  dvariable rho=fsh.length_rho;
  dvariable psi=fsh.length_psi;
  dvariable var=exp(fsh.log_length_variance);
  //dvariable rho=0.9;
  //dvariable var=4.;
  // *********************************************************
  // *********************************************************
  int nlint=fsh.nlint;
  
  dvar_matrix S(1,nlint,1,nlint);
  dvar_matrix SUB(1,nlint-1,1,nlint-1);

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
  dmatrix M(1,nlint-1,1,nlint);
  dvar_matrix SM(1,nlint-1,1,nlint);
  M.initialize();
  for (int i=1;i<=nlint-1;i++)
  {
    M(i,i)=1.0;
    M(i,nlint)=-1;
  }
  for (int i=1;i<=nlint-1;i++)
  {
    SM(i)=S*M(i);
  }
  for (int i=1;i<=nlint-1;i++)
  {
    SUB(i,i)=M(i)*SM(i);
    for (int j=1;j<i;j++)
    {
      SUB(i,j)=M(i)*SM(j);
      SUB(j,i)=SUB(i,j);
    }
  }
  SUB*=var;
  
  dvar_matrix chinv=inv(choleski_decomp(SUB));
  //cout << chinv << endl;
  //exit(1);
  dvariable ln_det_choleski_inv=0.0;
  int sgn=1;
  for (int i=1;i<=nlint-1;i++)
  {
    if (value(chinv(i,i))>0.0)
    {
      ln_det_choleski_inv+=log(chinv(i,i));
    }
    else
    {
      sgn*=-1;
      ln_det_choleski_inv+=log(-chinv(i,i));
    }
  }
  // *********************************************************
  // *********************************************************

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

    for (int ip=1;ip<=ntimes;ip++)
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
      {
        // if (fsh.len_sample_size(ir,ip,fi)=0 there is no
        // sample for these data

        if (fsh.len_sample_size(ir,ip,fi)>0)
        {
          // fsh.len_freq(ir,ip,fi) are the observed length frequencies
          // normalized to sum to 1. for now we add small epsilon to
          // them here
          MY_DOUBLE_TYPE epsilon=1.e-6;
          if (fsh.parest_flags(294))
            epsilon=fsh.parest_flags(294)/1.e+8;
            
          dvector ole=fsh.len_freq(ir,ip,fi)+epsilon;
          ole/=sum(ole);
          dvar_vector pole=fsh.tprob(ir,ip,fi)+epsilon;
          pole/=sum(pole);
          dvar_vector lole1=log(ole);
          dvar_vector lpole1=log(pole);
          if (fsh.parest_flags(299))
          {
             lole1=elem_prod(lole1,sqrt(pole));
             lpole1=elem_prod(lpole1,sqrt(pole));
          }

          dvar_vector lole=lole1(1,nlint-1)-lole1(nlint);
          dvar_vector lpole=lpole1(1,nlint-1)-lpole1(nlint);

          dvar_vector diff=lole-lpole;

          // chinv is the inverse of the Choleski decomp of 
          // the covariance matrix

          dvar_vector e=chinv*diff;
          if (fsh.parest_flags(295))
          {
             //e*=sqrt(fsh.len_sample_size(ir,ip,fi));
             e*=pow(fsh.len_sample_size(ir,ip,fi),fsh.length_exp);
          }
          // ln_det_sigma is the log of the det of the covariance matrix

          int pdftype=fsh.parest_flags(293);
          switch(pdftype)
          {
          case 0: // normal
            {
              const MY_DOUBLE_TYPE ltpi=0.5*log(2.0*3.1415926535);
              len_tmp+=(nlint-1)*ltpi + 0.5*(e*e) - ln_det_choleski_inv;
              if (fsh.parest_flags(299))
              {
                cerr << "not implemented yet" << endl;
                ad_exit(1);
              }
              if (fsh.parest_flags(295))
                len_tmp-=(nlint-1)*fsh.length_exp
                  *log(fsh.len_sample_size(ir,ip,fi));
              break;
            }
          case 1: // student
            {
              //dvariable v=50.0;  // this will be the degrees of freedoom for
                               // the students t
              dvariable v=exp(fsh.log_length_dof);  // this will be 
                                          // the degrees of freedoom for
              MY_DOUBLE_TYPE p=nlint-1.0;
              const MY_DOUBLE_TYPE lppi2=0.5*p*log(3.1415926535);
              len_tmp+= -gammln(0.5*(v+p)) + gammln(0.5*v) 
                      + lppi2 +(0.5*p)*log(v)
                      +0.5*(v+p)*log(1.0+e*e/v) 
                      - ln_det_choleski_inv;
              if (fsh.parest_flags(299))
                len_tmp-=0.5*sum(log(pole));
              if (fsh.parest_flags(295))
                len_tmp-=p*fsh.length_exp*log(fsh.len_sample_size(ir,ip,fi));
              break;
            }
          default: // error
            cerr << "illegal value for parest_flags(293)" << endl;
            ad_exit(1);
          }
        }
      }
    }
  }
  return len_tmp;
}

dvariable logistic_normal_weight_fit(dvar_len_fish_stock_history& fsh,
  d3_array& total_wght)
{
  // *********************************************************
  // *********************************************************
  // need to calculate these. for the moment they are just
  // here so that this mess will compile.
  dvariable rho=fsh.weight_rho;
  dvariable psi=fsh.weight_psi;
  dvariable var=exp(fsh.log_weight_variance);
  // *********************************************************
  // *********************************************************
  int nwint=fsh.nwint;
  
  dvar_matrix S(1,nwint,1,nwint);
  dvar_matrix SUB(1,nwint-1,1,nwint-1);

  dvar_vector rp(0,nwint-1);
  rp(0)=1;
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
  dmatrix M(1,nwint-1,1,nwint);
  dvar_matrix SM(1,nwint-1,1,nwint);
  M.initialize();
  for (int i=1;i<=nwint-1;i++)
  {
    M(i,i)=1.0;
    M(i,nwint)=-1;
  }
  for (int i=1;i<=nwint-1;i++)
  {
    SM(i)=S*M(i);
  }
  for (int i=1;i<=nwint-1;i++)
  {
    SUB(i,i)=M(i)*SM(i);
    for (int j=1;j<i;j++)
    {
      SUB(i,j)=M(i)*SM(j);
      SUB(j,i)=SUB(i,j);
    }
  }
  SUB*=var;
  
  dvar_matrix chinv=inv(choleski_decomp(SUB));
  //cout << chinv << endl;
  //exit(1);
  dvariable ln_det_choleski_inv=0.0;
  int sgn=1;
  for (int i=1;i<=nwint-1;i++)
  {
    if (value(chinv(i,i))>0.0)
    {
      ln_det_choleski_inv+=log(chinv(i,i));
    }
    else
    {
      sgn*=-1;
      ln_det_choleski_inv+=log(-chinv(i,i));
    }
  }
  // *********************************************************
  // *********************************************************

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

    for (int ip=1;ip<=ntimes;ip++)
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
      {
        // if (fsh.wght_sample_size(ir,ip,fi)=0 there is no
        // sample for these data

        if (fsh.wght_sample_size(ir,ip,fi)>0)
        {
          // fsh.wght_freq(ir,ip,fi) are the observed weight frequencies
          // normalized to sum to 1. for now we add small epsilon to
          // them here
          MY_DOUBLE_TYPE epsilon=1.e-6;
          if (fsh.parest_flags(284))
            epsilon=fsh.parest_flags(284)/1.e+8;
            

          dvector ole=fsh.wght_freq(ir,ip,fi)+epsilon;
          ole/=sum(ole);

          dvar_vector pole=fsh.wtprob(ir,ip,fi)+epsilon;
          pole/=sum(pole);
          dvar_vector lole1=log(ole);
          dvar_vector lpole1=log(pole);
          if (fsh.parest_flags(289))
          {
             lole1=elem_prod(lole1,sqrt(pole));
             lpole1=elem_prod(lpole1,sqrt(pole));
          }

          dvar_vector lole=lole1(1,nwint-1)-lole1(nwint);
          dvar_vector lpole=lpole1(1,nwint-1)-lpole1(nwint);

          dvar_vector diff=lole-lpole;

          // chinv is the inverse of the Choleski decomp of 
          // the covariance matrix

          dvar_vector e=chinv*diff;
          if (fsh.parest_flags(285))
          {
             //e*=sqrt(fsh.wght_sample_size(ir,ip,fi));
             e*=pow(fsh.wght_sample_size(ir,ip,fi),fsh.weight_exp);
          }
          // ln_det_sigma is the log of the det of the covariance matrix

          int pdftype=fsh.parest_flags(283);
          switch(pdftype)
          {
          case 0: // normal
            {
              if (fsh.parest_flags(289))
              {
                cerr << "not implemented yet" << endl;
                ad_exit(1);
              }

              const MY_DOUBLE_TYPE ltpi=0.5*log(2.0*3.1415926535);
              wght_tmp+=(nwint-1)*ltpi + 0.5*(e*e) - ln_det_choleski_inv;
              if (fsh.parest_flags(285))
#if !defined(NO_MY_DOUBLE_TYPE)
                wght_tmp-=(nwint-1.0L)*fsh.weight_exp
#else
                wght_tmp-=(nwint-1.0)*fsh.weight_exp
#endif
                  *log(fsh.wght_sample_size(ir,ip,fi));
              break;
            }
          case 1: // student
            {
              //dvariable v=50.0;  // this will be the degrees of freedoom for
                               // the students t
              dvariable v=exp(fsh.log_weight_dof);  // this will 
                                          // be the degrees of freedoom for
              MY_DOUBLE_TYPE p=nwint-1.0;
              const MY_DOUBLE_TYPE lppi2=0.5*p*log(3.1415926535);
              wght_tmp+= -gammln(0.5*(v+p)) + gammln(0.5*v) 
                      + lppi2 +(0.5*p)*log(v)
                      +0.5*(v+p)*log(1.0+e*e/v) 
                      - ln_det_choleski_inv;
              if (fsh.parest_flags(289))
                wght_tmp-=0.5*sum(log(pole));
              if (fsh.parest_flags(285))
                wght_tmp-=p*fsh.weight_exp*log(fsh.wght_sample_size(ir,ip,fi));
              break;
            }
          default: // error
            cerr << "illegal value for parest_flags(283)" << endl;
            ad_exit(1);
          }
        }
      }
    }
  }
  return wght_tmp;
}









