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
//extern int wght_print_flag;

void check_grouping_flags(ivector v);
ivector get_inv_group_ptr(ivector v);
symmetric_tridiagonal_dmatrix get_tridag_BH(MY_DOUBLE_TYPE crho,
  MY_DOUBLE_TYPE cvar,int mmin,int mmax);

extern int stupid_print_switch;

dvariable wght_self_scaling_multinomial_re_multi_rho_multi_var
  (dvar_len_fish_stock_history& fsh,d3_array& total_wght,int print_switch)
{
  ofstream * pofs=0;
  ofstream * pofs1=0;
  if (stupid_print_switch)
  {
    pofs=new ofstream("tc_wght_ssmult_fits");
    pofs1=new ofstream("wght_ssmult_fits");
  }
  int minsize=1;
  if (fsh.parest_flags(330)>0)
    minsize=fsh.parest_flags(330);
  // we assume that rho and var are grouped the same way
  // Don't have this extended to the parameter psi yet
 
  dvar_vector wrho=fsh.fish_pars(17);
  dvariable psi=fsh.weight_psi;
  dvar_vector wght_var=exp(fsh.fish_pars(19));
  int nwint=fsh.nwint;
  // ***************************************************************
  // ***************************************************************
  //  III   will move this later to near beginning of program
  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ivector& wght_rho_group_ptr=fsh.wtcs->wght_rho_group_ptr;
  imatrix& group_size_flag=fsh.wtcs->group_size_flag;
  ivector inv_wght_rho_group_ptr=get_inv_group_ptr(wght_rho_group_ptr);
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
  ivector & is=fsh.wtcs->is;
 
  int maxg=max(wght_rho_group_ptr);

  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ofstream * ppofs = 0;
  if (print_switch)
    ppofs = new ofstream("wsizemult");


  dvar_matrix wght_ld(minsize,maxsize,1,maxg);
  wght_ld.initialize();

  // correlation for random effects
  for (int ig=1;ig<=maxg;ig++)
  {
    dvariable rho=wrho(inv_wght_rho_group_ptr(ig));
    // multiply by variance
    dvariable var=wght_var(inv_wght_rho_group_ptr(ig));
    //rp*=var;
    for (int ii=minsize;ii<=maxsize;ii++)
    {
      if (group_size_flag(ig,ii))
      {
        wght_ld(ii,ig)=0.5*((ii-1)*log(1-rho*rho)+ii*log(var));
      }
    }
  }

  int blocksize=5000;
  dvariable wght_tmp=0.;
  const MY_DOUBLE_TYPE lstpi=0.91893853320467274177;
  // !!!!!  need to get fish_pars for this
  dvar_vector fp21=fsh.fish_pars(21);
  for (int ii=1;ii<=tot_nfi;ii+=blocksize)
  {
    //funnel_dvariable ftmp;
    dvariable ftmp;
    dvariable end_tmp=0.0;
    dvariable loc_tmp=0.0;
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
        int gp=wght_rho_group_ptr(pp);  //group this fishery is in   

        int lbic=left_bound(ic);
        int rbic=right_bound(ic);
        dvector wf=fsh.wght_freq(ir,ip,fi);
        dvector OP(lbic,rbic);
        OP=wf(lbic,rbic);
        //OP(lbic)=sum(wf(wf.indexmin(),lbic));
        //OP(rbic)=sum(wf(rbic,wf.indexmax()));
        dvector& eps=fsh.wght_eps(ir,ip,fi);
        if (!allocated(eps))
        {
          eps.allocate(1,rbic-lbic+1);
          eps.initialize();
        }
        dvar_vector P(lbic,rbic);
        dvar_vector pwf=fsh.wtprob(ir,ip,fi);
        P=pwf(lbic,rbic);
        //P(lbic)=sum(pwf(pwf.indexmin(),lbic));
        //P(rbic)=sum(pwf(rbic,pwf.indexmax()));
        OP.shift(1);
        P.shift(1);
        P+=1.e-10;
        // scale the predicted probabilities after RE's to sum
        // to sum(P) in the log-likelihood function
        //dvariable size_multiplier=
        //  pow(wght_sample_size(ir,ip,fi)/average_sample_size(gp),fp21(pp));

        // *************************************************************
        // *************************************************************
        //  size stuff to be moved to function
        // *************************************************************
        // *******************************************************
        // *******************************************************
        int old_lssmu_flag=0;
        dvariable cvN;
        if (old_lssmu_flag)
        {
          dvariable vN=exp(fsh.fish_pars(15,pp));
          const MY_DOUBLE_TYPE a=0.75;
          const MY_DOUBLE_TYPE beta1=0.99500;
          const MY_DOUBLE_TYPE beta2=0.9990;
          MY_DOUBLE_TYPE ss=fsh.wght_sample_size(ir,ip,fi);
          dvariable size_multiplier=
            pow(ss/average_sample_size(gp),fp21(pp));
          cvN=vN*size_multiplier;
          MY_DOUBLE_TYPE offset=0.0;
          const int n=sz;
          offset=(n-1)/10.;
          if (n>20) offset+=(n-20)/10.;
          if (n>40) offset+=(n-40)/20.;
          if (n>60) offset+=(n-60)/20.;
  
          dvariable avn;
          if (n<50)
            avn=pow(cvN,a);
          else if (n<75)
            avn=pow(cvN,0.9*a);
          else
            avn=pow(cvN,0.65*a);

#if !defined(NO_MY_DOUBLE_TYPE)
          dvariable alpha=beta1+(beta2-beta1)*(avn-1.0L)/avn;
#else
          dvariable alpha=beta1+(beta2-beta1)*(avn-1.0)/avn;
#endif
//          dvariable alpha=beta1+(beta2-beta1)*(avn-1.0L)/avn;
          dvariable svn=pow(cvN,alpha);
  
          if (old_lssmu_flag)
          {
            // this informs the model that the data were end compressed
            // via the old method
            // scale the predicted proababilities after RE's to sum
            // to sum(P) in the log-likelihood function
            if (fsh.parest_flags(330)>0 && fsh.parest_flags(331)>0)
            {
              dvariable S=sum(P);
              loc_tmp+=fsh.parest_flags(331)*square(log(S));
            }
            if (fsh.parest_flags(330)>0 && fsh.parest_flags(332)>0)
            {
              dvariable S=sum(P);
              loc_tmp-=cvN*log(S);
            }
          }
          loc_tmp+=wght_ld(sz,gp);
          loc_tmp+=sz*lstpi;
          loc_tmp-=gammln(svn+1.0+offset);
          dvar_vector sobs=cvN*OP;
          for (int ii=1;ii<=sobs.indexmax();ii++)
          {
            if (value(sobs(ii))>0.0)
            loc_tmp+=gammln(sobs(ii));
          }
        }
        else
        {
          dvariable vN=exp(fsh.fish_pars(15,pp));
          MY_DOUBLE_TYPE ss=fsh.wght_sample_size(ir,ip,fi);
          dvariable size_multiplier=
            pow(ss/average_sample_size(gp),fp21(pp));
          dvariable xx=vN*size_multiplier;
          dvariable fpen=0.0;
          xx=posfun(xx-0.04,.01,fpen)+0.04;
          MY_DOUBLE_TYPE ub=1000.;
          if (fsh.parest_flags(337))
            ub=fsh.parest_flags(337);
          cvN=-posfun(-xx+ub,.01,fpen)+ub;

          if (ppofs)
          {
            *ppofs << cvN << " " <<  ss/average_sample_size(gp)   
                   << " " << ss << " " << gp << " " << pp << endl;
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

          // this informs the model that the data were end compressed
          // scale the predicted proababilities after RE's to sum
          // to sum(P) in the log-likelihood function
          {
            dvariable S=sum(P);
            MY_DOUBLE_TYPE OS=sum(OP);
  
            loc_tmp-=cvN*log(S+1.e-6);
          }
  
          loc_tmp+=wght_ld(sz,gp);
          loc_tmp+=sz*lstpi;
          
          MY_DOUBLE_TYPE mvv=mean(log(sample_sizes));
          const MY_DOUBLE_TYPE l1000=log(1000.);
          dvariable t1=(ln-l1000)/l1000;
          dvariable tmp=exp((aa(6)+aa(7)*t1)*ln);
          loc_tmp-=
            gammln(tmp+exp(aa(1)+t1*(aa(2)+t1*(aa(3)+t1*(aa(4)+t1*aa(5))))));
          dvar_vector sobs=cvN*OP;
          for (int ii=1;ii<=sobs.indexmax();ii++)
          {
            if (value(sobs(ii))>0.0)
            loc_tmp+=gammln(sobs(ii)+ftoff);
          }
        }   
 
        dvariable vrho=wrho(pp);
        dvariable  vvar=wght_var(pp);
        MY_DOUBLE_TYPE tmp=0.5*((sz-1)*log(value(1-vrho*vrho))+sz*log(value(vvar)));
        if ( (fabs(tmp-value(wght_ld(sz,gp)))> 1.e-15))
        {
           cerr << tmp << " " << wght_ld(sz,gp) << endl;
           ad_exit(1);
        }

        loc_tmp+=lognormal_multinomial_log_likelihood(cvN,OP,P,
           vrho,vvar,eps,ic,fsh.pccsa);
        if (pofs)
        {
          int mmin=eps.indexmin();
          int mmax=eps.indexmax();
          dvector tmp(mmin,mmax);
          for (int i=mmin;i<=mmax;i++)
          {
            tmp(i)=phi(value(P(i)),eps(i));
          }
          int realmin=fsh.wtprob(ir,ip,fi).indexmin();
          int realmax=fsh.wtprob(ir,ip,fi).indexmax();
          dvector tmp1(realmin,realmax);
          tmp1.initialize();
          tmp1(left_bound(ic),right_bound(ic)).shift(1)=OP;
//          (*pofs) << ir << " " << ip << " " << fi << endl;
//          (*pofs1) << ir << " " << ip << " " << fi << endl;
          (*pofs) << ir << " " << ip << " " << fi << " " << pp << endl;
          (*pofs1) << ir << " " << ip << " " << fi << " " << pp << endl;
          (*pofs) << OP << endl;
          (*pofs1) << tmp1 << endl;
          (*pofs) << P/sum(P) <<   endl;
          tmp1(left_bound(ic),right_bound(ic)).shift(1)=value(P);
          (*pofs1) << tmp1/sum(tmp1) << endl;
          (*pofs) << tmp/sum(tmp) <<   endl;
          tmp1(left_bound(ic),right_bound(ic)).shift(1)=tmp;
          (*pofs1) << tmp1/sum(tmp1) <<   endl;
          (*pofs1) << "sum(eps) " <<  sum(eps) <<   endl;
          (*pofs) << "sum(eps) " <<  sum(eps) <<   endl;
        }
      }
    }
    ftmp=loc_tmp;
    end_tmp=ftmp;
    wght_tmp+=end_tmp;
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
