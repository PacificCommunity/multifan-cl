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
//#include "makebig2.cpp"


extern int my_debug_flag;
void check_grouping_flags(ivector v);
extern unsigned char _fsb[];
extern unsigned char _fsb6[];
extern unsigned char _fsb7[];

void  add_test_par(const prevariable & _y,const prevariable& _x);
void  add_test_par1(const prevariable & _y,const prevariable& _x);
void  check_test_par(const dvar_vector & _y);
void  check_test_par(const dvar_vector & _y);
void  check_test_par(const prevariable & _y);

void print_covariance_matrices(const banded_symmetric_dvar_matrix& vbsd2,
  ofstream & ofcv,int ir,int ip,int fi);
  

ivector get_inv_group_ptr(ivector v);
symmetric_tridiagonal_dmatrix get_tridag_BH(MY_DOUBLE_TYPE crho,
  MY_DOUBLE_TYPE cvar,int mmin,int mmax);

extern int stupid_print_switch;


dvariable new_wght_self_scaling_multinomial_re_multi_rho_multi_var
  (dvar_len_fish_stock_history& fsh,d3_array& total_freq,int print_switch)
{
  int one_time=0.0;
  ofstream * pofcv=0;
  ofstream * pofs=0;
  ofstream * pofs1=0;
  int kludge_flag=0;
  if (stupid_print_switch)
  {
    pofcv=new ofstream("ssmult_re_covariance");
    pofs=new ofstream("tc_wght_ssmult_fits");
    pofs1=new ofstream("wght_ssmult_fits");
  }
  int minsize=1;
  if (fsh.parest_flags(330)>0)
    minsize=fsh.parest_flags(330);
 
  dvar_vector wrho=fsh.fish_pars(17);
  dvariable psi=fsh.weight_psi;
  dvar_vector wght_var=exp(fsh.fish_pars(19));
  dvar_vector wght_vexp=fsh.fish_pars(27);
  dvariable lvpen=norm2(wght_var);

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

  //dvar_matrix wght_ld(minsize,maxsize,1,maxg);
  dvar_matrix wght_ld2(minsize,maxsize,1,maxg);
  wght_ld2.initialize();

  if(!fsh.wpcsa)
  {
    fsh.wpcsa = new pcubic_spline_array[maxg];
    fsh.wpcsa--;
    for (int ig=1;ig<=maxg;ig++)
    {
      fsh.wpcsa[ig]=0;
    } 
    for (int ig=1;ig<=maxg;ig++)
    {
      switch (fsh.parest_flags(355))
      {
      case 0:
        fsh.wpcsa[ig] = new cubic_spline_array(_fsb);
        break;
      case 1:
        fsh.wpcsa[ig] = new cubic_spline_array(_fsb7);
        break;
      case 2:
        fsh.wpcsa[ig] = new cubic_spline_array(_fsb6);
        break;
      default:
        cerr << "Illegal value for fsh.parest_flags(355)"
          " value is " << fsh.parest_flags(355) << endl;
        ad_exit(1);
      }
    }
  }
  for (int ig=1;ig<=maxg;ig++)
  {
    dvariable rho=wrho(inv_wght_rho_group_ptr(ig));
    fsh.wpcsa[ig]->set_z(rho);
    // !!!!!!!!!!!!!  new au-au stuff
    // divide by variance
    dvariable var=wght_var(inv_wght_rho_group_ptr(ig));
    for (int ii=minsize;ii<=maxsize;ii++)
    {
      if (group_size_flag(ig,ii))
      {
        wght_ld2(ii,ig)=-fsh.wpcsa[ig]->get_ln_det_choleski_inv(ii)
          +0.5*ii*log(var);
      }
    }
  }

  int blocksize=5000;
  dvariable wght_tmp=0.;
  dvariable wght_tmp1=0.;
  const MY_DOUBLE_TYPE lstpi=0.91893853320467274177;
  // !!!!!  need to get fish_pars for this
  dvar_vector fp21=fsh.fish_pars(21);
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

        int new_switch=1;
        dvariable vN;
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
        if (ppofs)
        {
          *ppofs << vN << " " <<  ss/average_sample_size(gp)   
                 << " " << ss  << " " << gp << " " << pp << endl;
        }

        dvector& eps=fsh.wght_eps(ir,ip,fi);
        if (!allocated(eps))
        {
          eps.allocate(1,rbic-lbic+1);
          eps.initialize();
        }
        if (fsh.parest_flags(137)==0)   //NMD_9Dec2019
        {
          dvariable Phi=0.0;
          dvariable tmp=0.0;
          int np=0;
          int nz=0;
          for (int ii=1;ii<=nslots;ii++)
          {
            if (obsc(ii)>0.0)
            {
              np++;
              tmp+=obsc(ii)*log(obsc(ii));
            }
            else
            {
              nz++;
            }
          }
  
          obsc/=sum(obsc);
          predc/=sum(predc);
          if (np<nslots) np++;
          MY_DOUBLE_TYPE multiplier=1.0;
#if !defined(NO_MY_DOUBLE_TYPE)
          Phi=multiplier*0.5*(np-1.0L)*log(vN)-vN*tmp;
#else
          Phi=multiplier*0.5*(np-1.0)*log(vN)-vN*tmp;
#endif
          wght_tmp-=Phi;
          //wght_tmp-=vN*(obsc*log(predc)-sum(predc));
          //wght_tmp-=vN*(obsc*log(predc));  //NMD_9Dec2019
        }
        else
        {
          dvariable lt2=wght_multiplier_calcs (fsh,ir,ip,fi,pp,gp,sz,
            average_sample_size,obsc);

          dvariable ln=log(vN);
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
          // to sum(predc) in the log-likelihood function
          {
            dvariable S=sum(predc);
            MY_DOUBLE_TYPE OS=sum(obsc);
  
            if (OS>1.e-10)
            {
              wght_tmp-=vN*log(S+1.e-6);
            }
          }
          
          MY_DOUBLE_TYPE mvv=mean(log(sample_sizes));
          const MY_DOUBLE_TYPE l1000=log(1000.);
          dvariable t1=(ln-l1000)/l1000;
          dvariable tmp=exp((aa(6)+aa(7)*t1)*ln);
          // KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
          dvariable lt3=0.0;
          lt3-=
            gammln(tmp+exp(aa(1)+t1*(aa(2)+t1*(aa(3)+t1*(aa(4)+t1*aa(5))))));
          // KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
          dvar_vector sobs=vN*obsc;
          for (int ii=1;ii<=sobs.indexmax();ii++)
          {
            if (value(sobs(ii))>0.0)
            lt3+=gammln(sobs(ii)+ftoff);
          }
          wght_tmp+=lt3;
        }
        dvariable vrho=wrho(pp);
        dvariable vvar=wght_var(pp);
        dvariable vexp=wght_vexp(pp);
        
        int ararflag=1;
        // KKKKKKKKKKKKKKKK
        //if (one_time<3)
        {
          one_time++;
          
          if (fsh.wpcsa[gp]==0)
            fsh.wpcsa[gp] = new cubic_spline_array(_fsb);
          fsh.wpcsa[gp]->set_z(vrho);
          int nn=predc.indexmax()-predc.indexmin()+1;
          banded_lower_triangular_dvar_matrix vbltd=
            fsh.wpcsa[gp]->get_ltcholeski_inv(nn);
          vbltd=(1.0/sqrt(vvar))*vbltd;
          int bw=vbltd.bandwidth();
          
          MY_DOUBLE_TYPE e=0.0001;
          for (int j=1;j<=nn;j++)
          {
            dvariable svinv=1.0/sqrt(vvar);
            for (int i=j;i<=min(j+bw-1,nn);i++)
            {
              vbltd(i,j)*=pow((e+1.0/predc(i))/(e+1.0/nn),vexp);
            }
          }
          
          
          banded_symmetric_dvar_matrix vbsd2=mult_trans_mult(vbltd);
          if (pofcv)
          {
            print_covariance_matrices(vbsd2,*pofcv,ir,ip,fi);
          }
          MY_DOUBLE_TYPE crho=value(vrho);
          MY_DOUBLE_TYPE cvar=value(vvar);
          MY_DOUBLE_TYPE N=value(vN);
          dvector p=value(predc);
          dvector eta_hat=lognormal_ss_multinomial_newton_raphson
              (N,obsc,p,value(vbltd),value(vbsd2),eps);
          wght_tmp1+=lognormal_multinomial_log_likelihood_cubic_spline
            (vN,obsc,predc,eps,ic,vbsd2,eta_hat);
          
          dvariable the_ld=wght_ld2(sz,gp);
          dvariable tmp=vexp*log(e+1.0/nn);
          int mmin=predc.indexmin();
          int mmax=predc.indexmax();
          for (int i=mmin;i<=mmax;i++)
          {
            the_ld-=vexp*(log(e+1.0/predc(i)));
          }
          the_ld+=(mmax-mmin+1)*tmp;
          //cout << the_ld << ln_det(make_dmatrix(value(vbsd2)) 
          wght_tmp1+=the_ld;
          wght_tmp1+=sz*lstpi;
        }
       
        kludge_flag++;
        if (pofs)
        {
          int mmin=eps.indexmin();
          int mmax=eps.indexmax();
          dvector tmp(mmin,mmax);
          for (int i=mmin;i<=mmax;i++)
          {
            tmp(i)=phi(value(predc(i)),eps(i));
          }
          int realmin=fsh.wtprob(ir,ip,fi).indexmin();
          int realmax=fsh.wtprob(ir,ip,fi).indexmax();
          dvector tmp1(realmin,realmax);
          tmp1.initialize();
          tmp1(left_bound(ic),right_bound(ic)).shift(1)=obsc;
//          (*pofs) << ir << " " << ip << " " << fi << endl;
//          (*pofs1) << ir << " " << ip << " " << fi << endl;
          (*pofs) << ir << " " << ip << " " << fi << " " << pp << endl;
          (*pofs1) << ir << " " << ip << " " << fi << " " << pp << endl;
          (*pofs) << obsc << endl;
          (*pofs1) << tmp1 << endl;
          (*pofs) << predc/sum(predc) <<   endl;
          tmp1(left_bound(ic),right_bound(ic)).shift(1)=value(predc);
          (*pofs1) << tmp1/sum(tmp1) << endl;
          (*pofs) << tmp/sum(tmp) <<   endl;
          tmp1(left_bound(ic),right_bound(ic)).shift(1)=tmp;
          (*pofs1) << tmp1/sum(tmp1) <<   endl;
          (*pofs1) << "sum(eps) " <<  sum(eps) <<   endl;
          (*pofs) << "sum(eps) " <<  sum(eps) <<   endl;
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
  // KKKKKKKKKKKKKKKK
  one_time=0;
  cout << "WWWWWWWWWW " << wght_tmp1 << endl;
  wght_tmp+=wght_tmp1;
  return wght_tmp;
}

