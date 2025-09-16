/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#define HOME_VERSION
#include "all.hpp"
//#include "here.h"
#include "variable.hpp"
#include "newmprot.hpp"
extern int TURNOFF;

dvariable recrpen(dvar_len_fish_stock_history& fsh,int is,int print_switch,
		  MY_DOUBLE_TYPE recr_weight) //NMD 14May 2012
{
  dvariable xy=0.0;
  dvariable ppf_tmp=0.0;   //NMD_22nov2023
  MY_DOUBLE_TYPE pen=0.01;
  if (fsh.age_flags(111)>0) pen=fsh.age_flags(111);
  dvar_vector trecr;
  if (is==1)
    trecr=fsh.recr(2,fsh.last_real_year);  //NMD_19May_2016
  else
    trecr=fsh.pmsd->recr(is)(2,fsh.last_real_year);  //NMD_19May2016

  dvariable tmp4=pen*norm2(first_difference(trecr));
  ppf_tmp=value(tmp4);  //NMD_22nov2023
  if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_22nov2023
  {
    fsh.ppstf->temporal_recr_dev_pen+=ppf_tmp;
  }

  if (fsh.parest_flags(149)>=0)
  {
    dvariable tmp_stuff;
    tmp_stuff=recr_weight*norm2(trecr-mean(trecr));
    tmp4+=tmp_stuff;
  }
  else
  {
    dvariable tmp_stuff;
    tmp_stuff=.5*fsh.last_real_year*log(1.e-6+norm2(trecr-mean(trecr)));
    tmp4+=tmp_stuff;
  }
  ppf_tmp=value(tmp4)-ppf_tmp;  //NMD_22nov2023
  if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_22nov2023
  {
    fsh.ppstf->norm_temporal_recr_pen+=ppf_tmp;
    ppf_tmp=0.0;
  }

  dvector trend(trecr.indexmin(),trecr.indexmax());
  trend.fill_seqadd(-1.,2./(trend.size()-1.));
  trend=trend/norm(trend);

  dvariable tmp6=0.0;
  if (fsh.parest_flags(155)==0)
  {
    MY_DOUBLE_TYPE trend_wt=.0;
    if (fsh.parest_flags(153)>0)
    {
      trend_wt=fsh.parest_flags(153)/10.;
      tmp6=trend_wt*square(trecr*trend);
    }
    else if (fsh.parest_flags(153)==0)
    {
      trend_wt=.01;
      tmp6=trend_wt*square(trecr*trend);
    }
  }
  else
  {
    MY_DOUBLE_TYPE trend_wt=.0;
    if (fsh.parest_flags(153)>0)
    {
      trend_wt=fsh.parest_flags(153)/10.;
    }
    else if (fsh.parest_flags(153)==0)
    {
      trend_wt=.01;
    }
    dvar_vector lintrend=column(fsh.orth_recr,1);
    tmp6=trend_wt*norm2(lintrend);
  }
  ppf_tmp=value(tmp6);  //NMD_22nov2023
  if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_22nov2023
  {
    fsh.ppstf->recr_trend_pen+=ppf_tmp;
    ppf_tmp=0.0;
  }

  if (TURNOFF==1)
  {
   // tmp4=0.0;
   // tmp6=0.0;
  }
  //  cout << "#### Dbug callpen: before recr_trend pen= " << setprecision(12) << xy << endl; //NMD
    cout << "fishery recruitment trend penalty = "<< tmp6 << endl;
    cout << "tmp4 = "<< tmp4 << endl;
  if (print_switch)
  {
    cout << "fishery recruitment trend penalty = "<< tmp6 << endl;
  }
  xy+=tmp6;
  xy+=tmp4;
  if (print_switch)
  {
    cout << " after recruitment trend penalty  = " << xy << endl;
  }

  //  cout << "#### Dbug callpen: species " << is << " tmp4= " 
  //       << setprecision(12) << tmp4 << endl; //NMD
  //  cout << "#### Dbug callpen: species " << is << " tmp6= " 
  //       << setprecision(12) << tmp6 << endl; //NMD
  //  cout << "#### Dbug callpen: after recr_trend pen= " << setprecision(12) << xy << endl; //NMD

  dvariable r2pen=0.0;
  MY_DOUBLE_TYPE pen_wght=.1;
  
  if (fsh.age_flags(76)>=0) 
  {
    if (fsh.age_flags(76)>0) 
    {
      pen_wght=fsh.age_flags(76);
    }
    if (!fsh.age_flags(57))
    {
      if (fsh.nyears>=3)
      {
        r2pen+=pen_wght*norm2(first_difference(first_difference(trecr)));
      }
    }
    else // split up recruitment curverature penalty by period
    {
      int nyears=fsh.last_real_year;
      int tmult=fsh.age_flags(57);
      if (nyears%tmult)
      {
        cerr << "Number of years is not a multiple of kludge factor" << endl;
      }
      else
      {
        int ryears=nyears/tmult;
        int icount=0;
        int i;
        for (i=1;i<=nyears;i+=tmult)
        {
          icount++;
        } 
        dvar_matrix recr_per(1,tmult,1,icount);
        int ii=1;
        for (i=1;i<=nyears;i+=tmult)
        {
          for (int j=1;j<=tmult;j++)
          {
            if(i==1)  //NMD_19May2016
              recr_per(j,ii)=0.0;
            else
              recr_per(j,ii)=trecr(i+j-1);
          }
          ii++;
        } 
        for (int j=1;j<=tmult;j++)
        {
          r2pen+=.01*norm2(first_difference(first_difference(recr_per(j))));
        }
      }
    }
  }
  
  //if (TURNOFF==0)
  {
    if (print_switch)
    {
      cout << "fishery recruitment curv penalty = "<< r2pen << endl;
    }
    xy+=r2pen;
  }
  ppf_tmp=value(r2pen);  //NMD_22nov2023
  if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_22nov2023
  {
    fsh.ppstf->recr_curve_pen+=ppf_tmp;
    ppf_tmp=0.0;
  }

  if (fsh.parest_flags(148))
  {
    dvariable last_wt=fsh.parest_flags(148)/10.0;
    dvariable rlast_pen=last_wt*square(trecr(fsh.last_real_year)-
                        fsh.recr(fsh.last_real_year-1));
    xy+=rlast_pen;
    cout << "Last recruit. penalty = " << rlast_pen << endl; 
  }
  //  cout << "#### Dbug callpen: after last recr pen= " << setprecision(12) << xy << endl; //NMD
  return xy;
}
dvariable recruitment_autocorrelation(dvar_len_fish_stock_history& fsh,int is,
  int print_switch, MY_DOUBLE_TYPE recr_weight) 
{
  dvar_vector trecr(2,fsh.last_real_year);
  if (is==1)
    trecr=fsh.recr(2,fsh.last_real_year); 
  else
    trecr=fsh.pmsd->recr(is)(2,fsh.last_real_year); 
  trecr.shift(1);
  trecr-=mean(trecr);
  
  int mmin=1;
  int mmax=trecr.indexmax();
  //dvar_matrix M(mmin,mmax,mmin,mmax);
  dvariable rho=fsh.species_pars(1,is);

  dvar_vector vtmp=chinv_multiply(rho,trecr);
  //dvariable ld2=(mmax-1)*log(1-rho*rho);

  //dvariable tmp1=recr_weight*(vtmp*vtmp);
  MY_DOUBLE_TYPE var=0.5/recr_weight;
  dvariable tmp1=0.5*mmax*log(var)-0.5*log(1.0-rho*rho)
      +0.5*((1-rho*rho)/var)*(vtmp*vtmp);


  return tmp1;
  //return tmp1+0.5*ld2;
}

/*
dvariable recruitment_autocorrelation_moment_estimator
  (dvar_len_fish_stock_history& fsh,int is,int print_switch,
   MY_DOUBLE_TYPE recr_weight)
{
  int lry=fsh.last_real_year;
  dvar_vector trecr0(2,lry);
  //because fsh.recr has dimensions (2,last_real_year)
  if (is==1)
    trecr0=fsh.recr;
  else
    trecr0=fsh.pmsd->recr(is);
  trecr0=trecr0.shift(1);
  trecr0-=mean(trecr0);
  int n=trecr0.indexmax();
  dvar_vector trecr=trecr0(2,n);
  trecr.shift(1);
  dvar_vector trecr1=trecr0(1,n-1);
  trecr.shift(1);
  // should we set rho=rhoest?
#if !defined(NO_MY_DOUBLE_TYPE)
  dvariable rhoest=(trecr*trecr1)/norm2(trecr0)*(n/(n-1.0L));
#else
  dvariable rhoest=(trecr*trecr1)/norm2(trecr0)*(n/(n-1.0));
#endif
//  dvariable rho=fsh.species_pars(1,is);

  int mmin=1;
  int mmax=trecr.indexmax();
  dvar_vector vtmp=chinv_multiply(rhoest,trecr);
  dvariable ld2=(mmax-1)*log(1-rhoest*rhoest);
  dvariable tmp1=recr_weight*(vtmp*vtmp);
  return tmp1+0.5*ld2;
}
*/
dvariable recruitment_autocorrelation_moment_estimator
  (dvar_len_fish_stock_history& fsh,int is,int print_switch,
   MY_DOUBLE_TYPE recr_weight)
{
  cout << "XXX test using species pars for moment estimator" << endl;
  int lry=fsh.last_real_year;
  dvar_vector trecr0(2,lry);
  //because fsh.recr has dimensions (2,last_real_year)
  if (is==1)
    trecr0=fsh.recr;
  else
    trecr0=fsh.pmsd->recr(is);
  trecr0=trecr0.shift(1);
  trecr0-=mean(trecr0);
  int n=trecr0.indexmax();
  dvar_vector trecr=trecr0(2,n);
  trecr.shift(1);
  dvar_vector trecr1=trecr0(1,n-1);
  trecr.shift(1);
  // should we set rho=rhoest?
#if !defined(NO_MY_DOUBLE_TYPE)
  dvariable rhoest=(trecr*trecr1)/norm2(trecr0)*(n/(n-1.0L));
#else
  dvariable rhoest=(trecr*trecr1)/norm2(trecr0)*(n/(n-1.0));
#endif
  // save the rhoest in species pars
  fsh.species_pars(6,is)=rhoest;

  int mmin=1;
  int mmax=trecr.indexmax();
  dvar_vector vtmp=chinv_multiply(rhoest,trecr);
  dvariable ld2=(mmax-1)*log(1-rhoest*rhoest);
  dvariable tmp1=recr_weight*(vtmp*vtmp);
  return tmp1+0.5*ld2;
}


dvariable bh_recruitment_autocorrelation_moment_estimator
  (dvar_len_fish_stock_history& fsh,int is,int print_switch,
   MY_DOUBLE_TYPE wght)
{
  fsh.pmsd_error();
  dvariable tmp=0;
  int ibh;
  int jbh;
  int nsize;
  dvar_vector w;
  if (fsh.age_flags(199)<=0) // Use residuals over full model period ND 2Mar2011
  {
    if (is==1)
    {
      ibh=fsh.bh_recr_devs.indexmin();
      jbh=fsh.bh_recr_devs.indexmax();
      w=fsh.bh_recr_devs;
    }
    else
    {
      ibh=fsh.pmsd->bh_recr_devs.indexmin();
      jbh=fsh.pmsd->bh_recr_devs.indexmax();
      w=fsh.pmsd->bh_recr_devs(is);
    }
  }
  else     // Use residuals for the user-specified period
  {
    ibh=fsh.last_real_year-fsh.age_flags(199)+1;    //ND 8Mar2011
    jbh=fsh.last_real_year-fsh.age_flags(200);    //ND 8Mar2011
    if (fsh.age_flags(94)==3  || fsh.age_flags(182) ) //annualised case
    {
      ibh=(ibh-1)/fsh.age_flags(57)+1;
      jbh=(jbh-1)/fsh.age_flags(57)+1;
    }
    w=fsh.pmsd->bh_recr_devs(ibh,jbh);
  }
  nsize=jbh-ibh+1;
  // estimate rho with moment estimator
  dvar_vector mw=w-mean(w);
  int switch_empirical=1;
  dvariable num=mw(ibh,jbh-1)*mw(ibh+1,jbh).shift(ibh)/(jbh-ibh);
  dvariable den=norm2(mw)/(jbh-ibh+1);
  dvariable t=num/den;
  dvariable fpen=0.0;
  dvariable rhoest =  0.99-posfun(0.99-t,.01,fpen);
  dvariable rho;
  if (switch_empirical)
  {
    rho=rhoest;
    fsh.species_pars(6,is)=rhoest;
  }
  else
  {
    rho=fsh.species_pars(2,is);
  }
  dvar_vector v=chinv_multiply(rho,w);
  dvariable ld2=(nsize-1)*log(1-rho*rho);
  if (fsh.age_flags(179))
  {
    tmp=fsh.age_flags(179)*v*v+0.5*ld2;
  }
  else
  {
    tmp=wght*v*v+0.5*ld2;
  }
  return tmp;
}

void check_gramschmidt(dvector & w,dvector& v,int mmax,MY_DOUBLE_TYPE rho);

dvar_vector chinv_multiply(const prevariable rho,const dvar_vector& v)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();

  dvar_vector w(mmin,mmax);

  dvariable a=1.0/sqrt(1-rho*rho);
  dvariable b=-rho*a;

  for (int i=mmin;i<mmax;i++)
  {
    w(i)=a*v(i)+b*v(i+1);
  }
  w(mmax)=v(mmax);

 /*
  dvector cw=value(w);
  dvector cv=value(v);
  MY_DOUBLE_TYPE crho=value(rho);
  check_gramschmidt(cw,vv,mmax-mmin+1,crho);
 */

  return w;
}

void check_gramschmidt(dvector & w,dvector& v,int mmax,MY_DOUBLE_TYPE rho)
{
  w.shift(1);
  v.shift(1);
  dmatrix M(1,mmax,1,mmax);
  dvector ttmp(1,mmax);
  ttmp(1)=rho;
  for (int i=2;i<=mmax;i++)
  {
    ttmp(i)=ttmp(i-1)*rho;
  }
  for (int i=1;i<=mmax;i++)
  {
    M(i,i)=1.0;
    for (int j=1;j<i;j++)
    {
       M(i,j)= ttmp(i-j);
       M(j,i)=M(i,j);
    }
  }

  dmatrix N=inv(M);
  cout << setfixed() << setprecision(4) << setw(8) << M << endl << endl;
  cout << setfixed() << setprecision(4) << setw(8) << N << endl << endl;
  cout << setfixed() << setprecision(4) << setw(8) << choleski_decomp(N) 
       << endl;
  dvector u=v*choleski_decomp(N);
  cout <<  v  << endl;
  cout <<  u  << endl;
  cout <<  w  << endl;
  cout <<  v * N * v - w*w  << endl;
  ad_exit(1);
}  
#undef HOME_VERSION

