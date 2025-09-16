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

dvariable new_wght_self_scaling_multinomial_nore
  (dvar_len_fish_stock_history& fsh,d3_array& total_freq,int print_switch)
{
  //adfloat_except(1);
  ofstream * pofs=0;
  ofstream * pofs1=0;



  if (stupid_print_switch)
  {
    pofs=new ofstream("tc_wght_ssmult_fits");
    pofs1=new ofstream("wghr_ssmult_fits");
  }
  int minsize=1;
  if (fsh.parest_flags(330)>0)
    minsize=fsh.parest_flags(330);
  // ***************************************************************
  // ***************************************************************
  //  III   will move this later to near beginning of program

  ivector& wght_rho_group_ptr=fsh.wtcs->wght_rho_group_ptr;
  int& tot_nfi=fsh.wtcs->tot_nfi;
  dvector & average_sample_size=fsh.wtcs->average_sample_size;
  ivector & reg=fsh.wtcs->reg;
  ivector & per=fsh.wtcs->per;
  ivector & finc=fsh.wtcs->finc;
  ivector & left_bound=fsh.wtcs->left_bound;
  ivector & right_bound=fsh.wtcs->right_bound;
  ivector & sizes =fsh.wtcs->sizes;
  int maxsize=max(sizes);
  ivector ss=sort(sizes);
  int mmin=sizes.indexmin();
  int mmax=sizes.indexmax();
  ivector tmp(mmin,mmax);
  tmp(mmin)=ss(mmin);
 
  int ic=mmin;
  for (int i=sizes.indexmin()+1;i<=sizes.indexmax();i++)
  {
    if (ss(i)>ss(i-1))
    {
      tmp(++ic)=ss(i);
    }
  }
  ivector is=tmp(mmin,ic);
 
  dvar_vector fp21=fsh.fish_pars(21);

  int blocksize=5000;
  dvariable wght_tmp=0.;
  // !!!!!  need to get fish_pars for this
  ofstream * ppofs = 0;
  if (print_switch)
    ppofs = new ofstream("wsizemult");

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
        ftmp=0.0;
        int ir=reg(ic);
        int ip=per(ic);
        int fi=finc(ic);
        int pp=fsh.parent(ir,ip,fi);
        int gp=wght_rho_group_ptr(pp);  //group this fishery is in   
        int lbic=left_bound(ic);
        int rbic=right_bound(ic);
        dvector obsc(lbic,rbic);
        obsc=fsh.wght_freq(ir,ip,fi)(lbic,rbic);
        dvar_vector predc(lbic,rbic);
        predc=fsh.wtprob(ir,ip,fi)(lbic,rbic);
        obsc.shift(1);
        predc.shift(1);
        predc+=1.e-10;
        int nslots=rbic-lbic+1;

        int old_lssmu_flag=0;
        dvariable vN;
//        {     //NMD_8jul2020
          // this outputs vN scaled by the len_sample_size
          MY_DOUBLE_TYPE ss=fsh.wght_sample_size(ir,ip,fi);
          dvariable tvN=exp(fsh.fish_pars(15,pp));
          dvariable size_multiplier=
            pow(ss/average_sample_size(gp),fp21(pp));
          dvariable xx=tvN*size_multiplier;
          dvariable fpen=0.0;
          xx=posfun(xx-0.04,.01,fpen)+0.04;
          MY_DOUBLE_TYPE ub=1000.;
          if (fsh.parest_flags(337))
            ub=fsh.parest_flags(337);
          vN=-posfun(-xx+ub,.01,fpen)+ub;
//        }

        if (ppofs)
        {
          *ppofs << vN << " " <<  ss/average_sample_size(gp)   
                 << " " << ss  << " " << gp << " " << pp << endl;
        }
        dvariable Phi=0.0;
        dvariable tmp=0.0;
        int np=0;
        int nz=0;
        for (int ii=1;ii<=nslots;ii++)
        {
          if (obsc(ii)>0.0)
          {
            np++;
          }
          else
          {
            nz++;
          }
        }

        MY_DOUBLE_TYPE eps=0.001;
        obsc+=eps/nslots;
        obsc/=sum(obsc);
        predc/=sum(predc);
        for (int ii=1;ii<=nslots;ii++)
        {
          tmp+=obsc(ii)*log(obsc(ii));
        }
        if (np<nslots) np++;
        MY_DOUBLE_TYPE multiplier=1.0;
#if !defined(NO_MY_DOUBLE_TYPE)
        Phi=multiplier*0.5*(np-1.0L)*log(vN)-vN*tmp;
#else
        Phi=multiplier*0.5*(np-1.0)*log(vN)-vN*tmp;
#endif
        wght_tmp-=Phi;
        ftmp-=Phi; //NMD_8dec2023
        //wght_tmp-=vN*(obsc*log(predc)-sum(predc));
        wght_tmp-=vN*(obsc*log(predc));
        ftmp-=vN*(obsc*log(predc)); //NMD_8dec2023

        if (fsh.ppstf && print_switch && !fsh.af170q0)
        {
          fsh.ppstf->wghtcontrib(pp) += value(ftmp);
          fsh.ppstf->wghtcontrib_by_realization(ir,ip,fi) = value(ftmp);
        }

        if (pofs1)
        {
          int realmin=fsh.wtprob(ir,ip,fi).indexmin();
          int realmax=fsh.wtprob(ir,ip,fi).indexmax();
          dvector tmp1(realmin,realmax);
          tmp1.initialize();
          tmp1(left_bound(ic),right_bound(ic)).shift(1)=obsc;
          (*pofs1) << ir << " " << ip << " " << fi << " " << pp << endl;
          (*pofs1) << tmp1 << endl;
          tmp1(left_bound(ic),right_bound(ic)).shift(1)=value(predc); 
          (*pofs1) << tmp1/sum(tmp1) << endl;
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
  return wght_tmp;
}
