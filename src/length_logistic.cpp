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

void get_grouping_vector(dvar_len_fish_stock_history& fsh) 
{
  MY_DOUBLE_TYPE cutoff=0.01;
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
      {
        if (fsh.len_sample_size(ir,ip,fi)>0)
        {
          dvector &lf=fsh.len_freq(ir,ip,fi);
          dvector tmp(1,lf.indexmax());
          tmp.initialize();
          int ic=1;
          for (int i=1;i<=lf.indexmax();i++)
          {
            tmp(ic)+=lf(i);
            if (tmp(ic)>cutoff)
            {
              ic++;
            }
          }            
          dmatrix M(1,ic-1,1,lf.indexmax());
          M.initialize();
          tmp.initialize();
          for (int i=1;i<=lf.indexmax();i++)
          {
            M(ic,i)=1.0; 
            tmp(ic)+=lf(i);
            if (tmp(ic)>cutoff)
            {
              ic++;
            }
          }            
          cout << M << endl; 
          ad_exit(1);
        }
      }
    } 
  }
}


dvariable logistic_normal_fit_heteroscedastic(dvar_len_fish_stock_history& fsh, 
  d3_array& total_freq)
{
  //get_grouping_vector(fsh);
 
  dvariable rho=fsh.length_rho;
  dvariable psi=fsh.length_psi;
  dvariable var=exp(fsh.log_length_variance);
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
      cout << " ipp= " << ipp << endl;
      for (int ip=ipp;ip<=ipp+blocksize-1;ip++)
      {
        if (ip>ntimes) break;
        for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
        {
          if (fsh.len_sample_size(ir,ip,fi)>0)
          {
            // lmax will be the index of the slot for the 
            // largest proportion
            int lmax=nlint;
            // fsh.len_freq(ir,ip,fi) are the observed length frequencies
            // normalized to sum to 1. for now we add small epsilon to
            // them here
            MY_DOUBLE_TYPE epsilon=1.e-6;
            if (fsh.parest_flags(294))
              epsilon=fsh.parest_flags(294)/1.e+8;
              
            //get_grouping_vector();
            
            dvector ole=fsh.len_freq(ir,ip,fi)+epsilon;
            ole/=sum(ole);
            dvar_vector pole=fsh.tprob(ir,ip,fi)+epsilon;
            pole/=sum(pole);
            dvar_vector lole1=log(ole);
            dvar_vector lpole1=log(pole);
            dvar_vector diff1=lole1-lpole1;
            dvar_matrix M(1,nlint-1,1,nlint);
            dvar_matrix SM(1,nlint-1,1,nlint);
            M.initialize();
            if (lmax !=nlint)
            {
              cerr << "will need to fix followoing code for this" << endl;
              ad_exit(1);
            }
            for (int i=1;i<=nlint-1;i++)
            {
              M(i,i)=pow(pole(i),fsh.length_tot_exp);
              M(i,lmax)=-pow(pole(lmax),fsh.length_tot_exp);
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
           /*
            ofstream ofs("SUB");
            ofs << "S" << endl<<endl;
            ofs << S << endl<<endl;
            ofs << "M" << endl<<endl;
            ofs << M << endl<<endl;
            ofs << SUB << endl<<endl;
           */
            
            SUB*=var;
            
            dvar_vector diff(1,nlint-1);
            if (lmax==nlint)
            {
              diff=diff1(1,nlint-1)-diff1(nlint);
            }
            else
            {
              diff(1,lmax-1)=diff(1,lmax-1)-diff(lmax);
              diff(lmax,nlint-1)=diff(lmax+1,nlint).shift(lmax)
                -diff(lmax);
            }
  
            dvariable ldet=0.0;
            int sgn=1;
            dvariable sgn1=0.0;
            dvar_vector e=choleski_factor_solve(SUB,diff,ldet,sgn);
            dvariable ln_det_choleski_inv=-ldet;

           /*
            ofs << SUB << endl << endl;
            ofs << eigenvalues(SUB) << endl;
            cout << "MMMM" << endl;
            ad_exit(1);
           */

           /*
            {
              dvariable ldet=0.0;
              dvar_vector e1=choleski_solve(SUB,diff,ldet,sgn);
              cout << "norm2(e-e1)" << endl;
              cout << norm2(e-e1) << endl;
              cout << "ldet" << endl;
              cout << ldet << endl;
              cout << "ln_det_choleski_inv" << endl;
              cout << ln_det_choleski_inv << endl;
              ad_exit(1);
            }
           */


            if (fsh.parest_flags(295))
            {
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
                loc_tmp+= -gammln(0.5*(v+p)) + gammln(0.5*v) 
                        + lppi2 +(0.5*p)*log(v)
                        +0.5*(v+p)*log(1.0+e*e/v) 
                        -ln_det_choleski_inv;
  
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

