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

dvariable old_len_self_scaling_multinomial_nore
  (dvar_len_fish_stock_history& fsh,d3_array& total_freq,int print_switch)
{
  //adfloat_except(1);
  ofstream * pofs=0;
  ofstream * pofs1=0;
  if (stupid_print_switch)
  {
    pofs=new ofstream("tc_len_ssmult_fits");
    pofs1=new ofstream("len_ssmult_fits");
  }
  int minsize=1;
  if (fsh.parest_flags(320)>0)
    minsize=fsh.parest_flags(320);
  int nlint=fsh.nlint;
  // ***************************************************************
  // ***************************************************************
  //  III   will move this later to near beginning of program

  ivector& len_rho_group_ptr=fsh.tcs->len_rho_group_ptr;
  int& tot_nfi=fsh.tcs->tot_nfi;
  dvector & average_sample_size=fsh.tcs->average_sample_size;
  ivector & reg=fsh.tcs->reg;
  ivector & per=fsh.tcs->per;
  ivector & finc=fsh.tcs->finc;
  ivector & left_bound=fsh.tcs->left_bound;
  ivector & right_bound=fsh.tcs->right_bound;
  ivector & sizes =fsh.tcs->sizes;
  int maxsize=max(sizes);
  ivector ss=sort(sizes);
  int mmin=sizes.indexmin();
  int mmax=sizes.indexmax();
  ivector tmp(mmin,mmax);
  tmp(mmin)=ss(mmin);
  //ivector& is=fsh.tcs->is;
 
  int ic=mmin;
  for (int i=sizes.indexmin()+1;i<=sizes.indexmax();i++)
  {
    if (ss(i)>ss(i-1))
    {
      tmp(++ic)=ss(i);
    }
  }
  ivector is=tmp(mmin,ic);
 
  dvar_vector fp20=fsh.fish_pars(20);

  int blocksize=5000;
  dvariable len_tmp=0.;
  // !!!!!  need to get fish_pars for this
  ofstream * ppofs = 0;
  if (print_switch)
    ppofs = new ofstream("sizemult");

  for (int ii=1;ii<=tot_nfi;ii+=blocksize)
  {
    //funnel_dvariable ftmp;
    dvariable ftmp;
    //ad_begin_funnel();
    int ic=0;
    for (ic=ii;ic<=ii+blocksize-1;ic++)
    {
      if (ic>tot_nfi) break;
      int sz=sizes(ic);
      if (sz<minsize)
      {
        cout << "this cant happen" << endl;
        ad_exit(1);
      }
      if (sz>1)
      {
        int ir=reg(ic);
        int ip=per(ic);
        int fi=finc(ic);
        int pp=fsh.parent(ir,ip,fi);
        int gp=len_rho_group_ptr(pp);  //group this fishery is in   
        int lbic=left_bound(ic);
        int rbic=right_bound(ic);
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
        int old_lssmu_flag=0;
        dvariable cvN;

        MY_DOUBLE_TYPE ss=fsh.len_sample_size(ir,ip,fi);
        dvariable vN=exp(fsh.fish_pars(14,pp));
        dvariable size_multiplier=
          pow(ss/average_sample_size(gp),fp20(pp));
        dvariable xx=vN*size_multiplier;
        dvariable fpen=0.0;
        xx=posfun(xx-0.04,.01,fpen)+0.04;
        MY_DOUBLE_TYPE ub=1000.;
        if (fsh.parest_flags(336))
          ub=fsh.parest_flags(336);
        cvN=-posfun(-xx+ub,.01,fpen)+ub;

        if (ppofs)
        {
          *ppofs << cvN << " " <<  ss/average_sample_size(gp)   
                 << " " << ss  << " " << gp << " " << pp << endl;
        }

        dvariable ln=log(cvN);
        const int nslots=sz;

        dvector pp1=fsh.m_estimator_coffs(nslots);
        int ioff=1;
        int nss=static_cast<int>(pp1(1+ioff));
        ioff++;
        dvector sample_sizes(1,nss);
        sample_sizes=pp1(1+ioff,nss+ioff).shift(1);
        ioff+=nss;
        MY_DOUBLE_TYPE ftoff=pp1(1+ioff);
        ioff++;
        dvector aa(1,7);
        aa=pp1(1+ioff,7+ioff).shift(1);

        
        //if (fsh.parest_flags(320)>0 && fsh.parest_flags(322)>0)
        {
          dvariable S=sum(P);
          MY_DOUBLE_TYPE OS=sum(OP);
          len_tmp-=cvN*log(S+1.e-6);
        }
        MY_DOUBLE_TYPE mvv=mean(log(sample_sizes));
        const MY_DOUBLE_TYPE l1000=log(1000.);
        dvariable t1=(ln-l1000)/l1000;
        dvariable tmp=exp((aa(6)+aa(7)*t1)*ln);
        len_tmp-=
          gammln(tmp+exp(aa(1)+t1*(aa(2)+t1*(aa(3)+t1*(aa(4)+t1*aa(5))))));
        dvar_vector sobs=cvN*OP;
        for (int ii=1;ii<=sobs.indexmax();ii++)
        {
          if (value(sobs(ii))>0.0)
          len_tmp+=gammln(sobs(ii)+ftoff);
        }
        len_tmp-=sobs*log(P);
        if (pofs1)
        {
          int realmin=fsh.tprob(ir,ip,fi).indexmin();
          int realmax=fsh.tprob(ir,ip,fi).indexmax();
          dvector tmp1(realmin,realmax);
          tmp1.initialize();
          tmp1(left_bound(ic),right_bound(ic)).shift(1)=OP;
//          (*pofs1) << ir << " " << ip << " " << fi << endl;
          (*pofs1) << ir << " " << ip << " " << fi << " " << pp << endl;
          (*pofs1) << tmp1 << endl;
          tmp1(left_bound(ic),right_bound(ic)).shift(1)=value(P); //NMD_31May2016
          (*pofs1) << tmp1/sum(tmp1) << endl;
          tmp1(left_bound(ic),right_bound(ic)).shift(1)=value(P); //NMD_31May2016
          (*pofs1) << tmp1/sum(tmp1) <<   endl;
        }
      }
    }
  }
  if (ppofs)
  {
    delete ppofs;
    ppofs=0;
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
  return len_tmp;
}
