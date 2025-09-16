/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
#include "phi_stuff.h"
#include "tridiagonal_dmatrix.h"
void check_grouping_flags(ivector v);
ivector get_inv_group_ptr(ivector v);

extern int stupid_print_switch;
extern "C" void adfloat_except(int k);


dvariable dirichlet_multinomial(const dvector& _q,const dvar_vector& _p,
  MY_DOUBLE_TYPE Nmax,const prevariable& dmlambda)
{
  dvector q=_q/(1.e-20+sum(_q));
  dvar_vector p=_p/(1.e-20+sum(_p));
#if !defined(NO_MY_DOUBLE_TYPE)
  dvariable tmp=gammln(Nmax+1.0L);
#else
  dvariable tmp=gammln(Nmax+1.0);
#endif
  dvar_vector lp=dmlambda*p;
  dvector Nq=Nmax*q;
#if !defined(NO_MY_DOUBLE_TYPE)
  tmp-=sum(gammln(Nq+1.0L));
#else
  tmp-=sum(gammln(Nq+1.0));
#endif
  tmp+=gammln(dmlambda);
  tmp-=gammln(Nmax+dmlambda);
  tmp+=sum(gammln(Nq+lp));
  tmp-=sum(gammln(lp));
  return tmp;
}
  

dvariable len_dm_nore(dvar_len_fish_stock_history& fsh,
  d3_array& total_freq,int print_switch)
{
 
  ofstream * pofs=0;
  ofstream * pofs1=0;
  if (stupid_print_switch)
  {
    pofs=new ofstream("tc_len_ssmult_fits");
    pofs1=new ofstream("len_ssmult_fits");
  }
  if (print_switch && !fsh.af170q0)
  {
    fsh.ppstf->lencontrib.initialize();
    fsh.ppstf->lencontrib_by_realization.initialize();
  }
  ivector & sizes =fsh.tcs->sizes;

  dvar_vector fp22=fsh.fish_pars(22);
  dvar_vector fp23=fsh.fish_pars(23);

  int blocksize=5000;
  dvariable len_tmp_dm=0.;
  // !!!!!  need to get fish_pars for this
  ofstream * ppofs = 0;
  ofstream * ppofs1 = 0;
  if (print_switch)
  {
    ppofs = new ofstream("sizemult");
    ppofs1 = new ofstream("dmsizemult");
    int maxgp=max(fsh.tcs->len_rho_group_ptr);
    (*ppofs1) << " average sample sizes " <<endl;
    for (int i=1;i<=maxgp;i++)
    {
      (*ppofs1) << "  " << fsh.tcs->average_sample_size(i);
    }
    (*ppofs1) << endl;
    (*ppofs1) << " Group   Max   Effective lambda  Covariate Average  Relative          Moment      Observed" <<endl;
    (*ppofs1) << "         size    size                       size    samp size        estimator    size" <<endl;
  }

  for (int ii=1;ii<=fsh.tcs->tot_nfi;ii+=blocksize)
  {
    //funnel_dvariable ftmp;
    dvariable ftmp;
    //ad_begin_funnel();
    int ic=0;
    for (ic=ii;ic<=ii+blocksize-1;ic++)
    {
      if (ic>fsh.tcs->tot_nfi) break;
      int sz=sizes(ic);
      if (sz<fsh.tcs->minsize)
      {
        cout << "this cant happen" << endl;
        ad_exit(1);
      }
      if (sz>1)
      {
        int ir=fsh.tcs->reg(ic);
        int ip=fsh.tcs->per(ic);
        int fi=fsh.tcs->finc(ic);
        int pp=fsh.parent(ir,ip,fi);
        int gp=fsh.tcs->len_rho_group_ptr(pp);  //group this fishery is in   
        int lbic=fsh.tcs->left_bound(ic);
        int rbic=fsh.tcs->right_bound(ic);
        dvector lf=fsh.len_freq(ir,ip,fi);
        dvector OP(lbic,rbic);
        OP=lf(lbic,rbic);
        dvar_vector P(lbic,rbic);
        dvar_vector plf=fsh.tprob(ir,ip,fi);
        // XXXXXXXXXXXXXXXXX
        P=plf(lbic,rbic);
        OP.shift(1);
        P.shift(1);
        P+=1.e-10;

        // *************************************************************
        // *************************************************************
        //  size stuff to be moved to function
        // *************************************************************
        // *******************************************************
        // *******************************************************
        MY_DOUBLE_TYPE lss=fsh.len_sample_size(ir,ip,fi);

        //  add a small amount to relative sample size in case it is almost 0
        // which screws up the differentiabity
        const MY_DOUBLE_TYPE min_sample_size=0.001;
        MY_DOUBLE_TYPE relative_sample_size=min_sample_size+
           lss/(min_sample_size+fsh.tcs->average_sample_size(gp));
        dvariable size_multiplier_dm=
          pow(relative_sample_size,fp23(pp));

        dvariable dmlambda=exp(fp22(pp))*size_multiplier_dm;
        // !!!   XXXX need to deal with truncation in likelihood in
        // some way

        dvariable ttmp0=0.;
        MY_DOUBLE_TYPE Nmax=1000.;
        if (fsh.parest_flags(342)>0) Nmax=fsh.parest_flags(342);
        {
          dvariable cvN=Nmax*(1+dmlambda)/(Nmax+dmlambda);
          dvariable S=sum(P);
          MY_DOUBLE_TYPE OS=sum(OP);
  
          if (OS>1.e-10)
          {
            len_tmp_dm-=cvN*log(S+1.e-6);
            ttmp0=cvN*log(S+1.e-6);
          }
        }

        dvariable ttmp=dirichlet_multinomial(OP,P,Nmax,dmlambda);
        len_tmp_dm-=ttmp;

        if (ppofs1)
        {
          int m=P.indexmax()-P.indexmin()+1;
          MY_DOUBLE_TYPE moment_estimator=m/sum(elem_div(square(value(OP-P)),
            value(1.e-6+P)));
          (*ppofs1) << "   " << fsh.tcs->len_rho_group_ptr(pp) 
                    << std::right << setfixed() << setw(9) << setprecision(2)
                    << Nmax  
                    << std::right << setfixed() << setw(9) << setprecision(2)
                    <<  Nmax*(1+dmlambda)/(Nmax+dmlambda) 
                    << std::right << setfixed() << setw(9) << setprecision(2)
                    <<  dmlambda  
                    << std::right << setfixed() << setw(9) << setprecision(2)
                    <<  fp23(pp)  
                    << std::right << setfixed() << setw(9) << setprecision(2)
                    <<  exp(fp22(pp))  
                    << setfixed() << setw(9) << setprecision(2)
                    <<  relative_sample_size
                    << setfixed() << setw(9) << setprecision(2)
                    <<  size_multiplier_dm
                    << setfixed() << setw(9) << setprecision(2)
                    <<  moment_estimator
                    << setfixed() << setw(9) << setprecision(2)
                    <<  lss << endl;
        }

        if (pofs1)
        {
          int realmin=fsh.tprob(ir,ip,fi).indexmin();
          int realmax=fsh.tprob(ir,ip,fi).indexmax();
          dvector tmp1(realmin,realmax);
          tmp1.initialize();
          tmp1(fsh.tcs->left_bound(ic),fsh.tcs->right_bound(ic)).shift(1)=OP;
//          (*pofs1) << ir << " " << ip << " " << fi << endl;
          (*pofs1) << ir << " " << ip << " " << fi << " " << pp << endl;
          (*pofs1) << tmp1 << endl;
          tmp1(fsh.tcs->left_bound(ic),fsh.tcs->right_bound(ic)).shift(1)=value(P); //NMD_31May2016
          (*pofs1) << tmp1/sum(tmp1) << endl;
          tmp1(fsh.tcs->left_bound(ic),fsh.tcs->right_bound(ic)).shift(1)=value(P); //NMD_31May2016
          (*pofs1) << tmp1/sum(tmp1) <<   endl;
        }
        if (print_switch && !fsh.af170q0)
        {
          //contrib by fishery ....NMD_29May2017
          fsh.ppstf->lencontrib(fsh.parent(ir,ip,fi)) -= value(ttmp0);
          fsh.ppstf->lencontrib(fsh.parent(ir,ip,fi)) -= value(ttmp);
//          fsh.ppstf->lencontrib_by_realization(ir,ip,fi) = value(ttmp0);
          fsh.ppstf->lencontrib_by_realization(ir,ip,fi) -= value(ttmp0); //NMD_25Apr2022
          fsh.ppstf->lencontrib_by_realization(ir,ip,fi) -= value(ttmp);
        }
      }
    }
  }
  if (ppofs)
  {
    delete ppofs;
    ppofs=0;
  }
  if (ppofs1)
  {
    delete ppofs1;
    ppofs1=0;
  }
  if (pofs)
  {
    delete pofs;
    pofs=0;
  }
  if (pofs1)
  {
    delete pofs1;
    pofs1=0;
  }
  return len_tmp_dm;
}

dvariable len_dm_nore_notc(dvar_len_fish_stock_history& fsh,
  d3_array& total_freq,int print_switch)
{
  dvar_vector fp22=fsh.fish_pars(22);
  dvar_vector fp23=fsh.fish_pars(23);

  dvariable len_tmp_dm=0.;
  if (print_switch && !fsh.af170q0)
  {
    fsh.ppstf->lencontrib.initialize();
    fsh.ppstf->lencontrib_by_realization.initialize();
  }
  int mmin=fsh.fish_flags.indexmin();	  
  int mmax=fsh.fish_flags.indexmax();	  
  ivector gpvector(mmin,mmax);
  if(allocated(fsh.tcs->len_rho_group_ptr))
  {
    gpvector=fsh.tcs->len_rho_group_ptr;
  }
  else if(sum(column(fsh.fish_flags,68))>0)
  {	  
    gpvector=column(fsh.fish_flags,68);
  }  
  else
  {
    gpvector.fill_seqadd(1,1);
  }

// NMD_30May2017 - ensure average_sample_size calculation
  dvector average_sample_size;
  ivector num_samples;
  if(allocated(fsh.tcs->average_sample_size))
  {
    average_sample_size=fsh.tcs->average_sample_size;
  }
  else
  {
    int maxg=max(gpvector);
    average_sample_size.allocate(1,maxg);
    average_sample_size.initialize();
    num_samples.allocate(1,maxg);
    num_samples.initialize();
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
          int pp=fsh.parent(ir,ip,fi);
          int gp=gpvector(pp);  //group this fishery is in
          if (fsh.len_sample_size(ir,ip,fi)>0)
          {
            average_sample_size(gp)+=fsh.len_sample_size(ir,ip,fi);
            num_samples(gp)+=1;
          }
        }
      }
    }
    for (int i=1;i<=maxg;i++)
    {
      if (num_samples(i)>0)
        average_sample_size(i)/=num_samples(i);
    }
  }
// NMD_30May2017 - ensure average_sample_size calculation

  ofstream * ppofs1 = 0;
  if (print_switch)
  {
    ppofs1 = new ofstream("dmsizemult");
    int maxgp=max(gpvector);
    (*ppofs1) << " average sample sizes " <<endl;
    for (int i=1;i<=maxgp;i++)
    { 
      (*ppofs1) << "  " << average_sample_size(i);
    }
    (*ppofs1) << endl;
    (*ppofs1) << " Group   Max   Effective lambda  Covariate Average  Relative          Moment      Observed" <<endl;
    (*ppofs1) << "         size    size                       size    samp size        estimator    size" <<endl;
  }

  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++)
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
      {
        if(fsh.len_sample_size(ir,ip,fi)>0)
        {
          int pp=fsh.parent(ir,ip,fi);
          dvector OP=fsh.len_freq(ir,ip,fi);
          dvar_vector P=fsh.tprob(ir,ip,fi)+.0001;
          P/=sum(P);
          MY_DOUBLE_TYPE lss=fsh.len_sample_size(ir,ip,fi);
          //  add a small amount to relative sample size in case it is almost 0
          // which screws up the differentiabity
          const MY_DOUBLE_TYPE min_sample_size=0.001;
          //cerr << "Need to deal with average sample size before"
          //        " we can use this" << endl;
	  //
	  int gp=gpvector(pp);
//          MY_DOUBLE_TYPE relative_sample_size=min_sample_size+
//             lss/(min_sample_size+fsh.tcs->average_sample_size(gp));
          MY_DOUBLE_TYPE relative_sample_size=min_sample_size+
            lss/(min_sample_size+average_sample_size(gp));  //NMD_30May2017
          dvariable size_multiplier_dm=
            pow(relative_sample_size,fp23(pp));

          dvariable dmlambda=exp(fp22(pp))*size_multiplier_dm;

          dvariable ttmp0=0.;
          MY_DOUBLE_TYPE Nmax=1000.;
          if (fsh.parest_flags(342)>0) Nmax=fsh.parest_flags(342);
          {
            dvariable cvN=Nmax*(1+dmlambda)/(Nmax+dmlambda);
            dvariable S=sum(P);
            MY_DOUBLE_TYPE OS=sum(OP);
  
            if (OS>1.e-10)
            {
              len_tmp_dm-=cvN*log(S+1.e-6);
              ttmp0=cvN*log(S+1.e-6);
            }
          }
          dvariable ttmp=dirichlet_multinomial(OP,P,Nmax,dmlambda);
          len_tmp_dm-=ttmp;


          if (ppofs1)
          {
            int m=P.indexmax()-P.indexmin()+1;
            MY_DOUBLE_TYPE moment_estimator=m/sum(elem_div(square(value(OP-P)),
              value(1.e-6+P)));
//            (*ppofs1) << "   " << fsh.tcs->len_rho_group_ptr(pp) 
            (*ppofs1) << "   " << gpvector(pp)  //NMD_30May2017 
                    << std::right << setfixed() << setw(9) << setprecision(2)
                    << Nmax  
                    << std::right << setfixed() << setw(9) << setprecision(2)
                    <<  Nmax*(1+dmlambda)/(Nmax+dmlambda) 
                    << std::right << setfixed() << setw(9) << setprecision(2)
                    <<  dmlambda  
                    << std::right << setfixed() << setw(9) << setprecision(2)
                    <<  fp23(pp)  
                    << std::right << setfixed() << setw(9) << setprecision(2)
                    <<  exp(fp22(pp))  
                    << setfixed() << setw(9) << setprecision(2)
                    <<  relative_sample_size
                    << setfixed() << setw(9) << setprecision(2)
                    <<  size_multiplier_dm
                    << setfixed() << setw(9) << setprecision(2)
                    <<  moment_estimator
                    << setfixed() << setw(9) << setprecision(2)
                    <<  lss << endl;
          }
          if (print_switch && !fsh.af170q0)
          {
            //contrib by fishery ....NMD_29May2017
            fsh.ppstf->lencontrib(fsh.parent(ir,ip,fi)) -= value(ttmp0);
            fsh.ppstf->lencontrib(fsh.parent(ir,ip,fi)) -= value(ttmp);
//            fsh.ppstf->lencontrib_by_realization(ir,ip,fi) = value(ttmp0);
            fsh.ppstf->lencontrib_by_realization(ir,ip,fi) -= value(ttmp0);  //NMD_25apr2022
            fsh.ppstf->lencontrib_by_realization(ir,ip,fi) -= value(ttmp);
          }
        }
      }
    }
  }
  if (ppofs1)
  {
    delete ppofs1;
    ppofs1=0;
  }
  return len_tmp_dm;
}
