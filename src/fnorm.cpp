/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


#define USE_DD_NOT
#if defined(USE_DD)
#  if !defined(_MAC) && !defined(DAVES_VERSION)
#    include <qd/qd.h>
#    include <qd/fpu.h>
#  else
#    include <qd.h>
#  endif
#endif

#include "all.hpp"
//extern int ctlc_flag;
void test_fpu(void);
void dvar_fish_stock_history::get_f_normal(MY_DOUBLE_TYPE & f,dmatrix& Y,dvector& x,
  dvector& g,banded_symmetric_dmatrix& H,int Hessflag,int robflag,
  dmatrix& resids,dmatrix & wts,int wtflag,
  ivector& iswitch,imatrix& indx,MY_DOUBLE_TYPE _c)
{
  MY_DOUBLE_TYPE c=_c/3.0;
  int jj=1;
  dvector fvec(1,10000);
  //ofstream ofs("imp"); 
  MY_DOUBLE_TYPE s_eps=0.3;
  MY_DOUBLE_TYPE pen1=1.0/square(s_eps);
  MY_DOUBLE_TYPE thpen1=1.0/square(3.0*s_eps);
  MY_DOUBLE_TYPE sqpen1=square(pen1);
  MY_DOUBLE_TYPE sqthpen1=square(thpen1);
  //double s_eta=1.0;
  MY_DOUBLE_TYPE s_eta=0.1;
  MY_DOUBLE_TYPE pen2=1.0/(2.0*square(s_eta));
  int n=x.indexmax();
  //double c=0.05;
  f=0;
  g.initialize();
  H.initialize();
  int i;
  int ii=1;
  //dmatrix grouped_catchability_coffs(1,ngroups,1,num_real_grouped_fish_times);
  int j;
  
  int offset=0;
  for (j=1;j<=ngroups;j++)
  {
    if (iswitch(j))
    {
      grouped_catchability_coffs(j)=
        x(1+offset,num_real_grouped_fish_times(j)+offset).shift(1);
      offset+=num_real_grouped_fish_times(j);
    }
  }
  
  ivector grouping=column(fish_flags,29);
  MY_DOUBLE_TYPE fold=0.0;
  
  for (j=1;j<=num_fisheries;j++)
  {
    int nrft=num_real_fish_times(j);
    dvector eff_wt(1,nrft);
    if (fish_flags(j,13)==0)
    {
      eff_wt=10.0;
    }
    else if (fish_flags(j,13)<0)
    {
      eff_wt=-fish_flags(j,13)*
        sqrt(.01+effort_by_fishery(j)(1,nrft));
    }
    else
    {
      eff_wt=fish_flags(j,13);
    }
   
    if (fish_flags(j,10)==1)
    {
      if (robflag==1)
      {
        for (i=1;i<=num_real_fish_times(j);i++)
        {
          ii=indx(grouping(j),gfish_ptr(j,i)); 
    
          // only do this if we have fishing effort
          if (true_effort_by_fishery(j,i)>-0.5L)
          {
            MY_DOUBLE_TYPE r=Y(j,i)
              -grouped_catchability_coffs(grouping(j),gfish_ptr(j,i));
  
            resids(j,i)=r;
            MY_DOUBLE_TYPE r2=square(r);
            /*********************************************************/
            MY_DOUBLE_TYPE u1=exp(-0.5*pen1*r2);
            MY_DOUBLE_TYPE v1=exp(-0.5*thpen1*r2);
            MY_DOUBLE_TYPE u=1.e-30+(u1+c*v1);
            /*********************************************************/
  
            //ofs << setprecision(14) << u << " " << r2 << endl;
            f-=log(u);
            fvec(jj++)=-log(u);
            /*********************************************************/
            //double tmp1= exp(-0.5*pen1*r2)+c/square(1.0+0.5*pen1*r2);
  
            MY_DOUBLE_TYPE tmp1= u1*pen1+c*v1*thpen1;
  
            /*********************************************************/
            MY_DOUBLE_TYPE tmp= tmp1*r;
            
            g(ii)-=1.0/u*tmp;
           
            
            if (Hessflag)
            {
              MY_DOUBLE_TYPE a1= 1.0/square(u)*square(tmp);
              MY_DOUBLE_TYPE a2= 1.0/u*tmp1;
             /*********************************************************/
              MY_DOUBLE_TYPE a3= 1.0/u*(u1*sqpen1+c*v1*sqthpen1)*r2;
              H(ii,ii)+=a1+a2-a3;
            }
          }
        }
      }
      else
      {
        for (i=1;i<=num_real_fish_times(j);i++)
        {
          ii=indx(grouping(j),gfish_ptr(j,i)); 
    
          grouped_catchability_coffs(grouping(j),gfish_ptr(j,i));
    
          // only do this if we have fishing effort
          if (true_effort_by_fishery(j,i)>-0.5L)
          {
            MY_DOUBLE_TYPE r=Y(j,i)
              -grouped_catchability_coffs(grouping(j),gfish_ptr(j,i));
            resids(j,i)=r;
            MY_DOUBLE_TYPE r2=square(r);
            if (wtflag==0)
            {
              f+=0.5*pen1*r2;
              fvec(jj++)=0.5*pen1*r2;
              g(ii)-=pen1*r;
              if (Hessflag)
              {
                H(ii,ii)+=pen1;
              }
            }
            else
            {
              f+=0.5*wts(j,i)*pen1*r2;
              fvec(jj++)=0.5*wts(j,i)*pen1*r2;
              g(ii)-=wts(j,i)*pen1*r;
              if (Hessflag)
              {
                H(ii,ii)+=wts(j,i)*pen1;
              }
            }
          }
        }
      }
    }
  }

  const ivector& xcatflags=(column(fish_flags,15));  
  ivector& catflags=(ivector&) xcatflags;  
  ii=1;
  for (j=1;j<=ngroups;j++)
  {
    if (fish_flags(gfish_index(j,1),10)==1)
    {
      ii++;
      for (i=2;i<=num_real_grouped_fish_times(j);i++)
      {
        //double r=(x(i)-x(i-1));
        
        int ff15=catflags(gfish_index(j,1));
        MY_DOUBLE_TYPE pp2;
        if (ff15)
        {
          pp2=ff15*grouped_between_times(j,i);
        }
        else
        {
          pp2=pen2*grouped_between_times(j,i);
        }
        MY_DOUBLE_TYPE r=grouped_catchability_coffs(j,i)
          -grouped_catchability_coffs(j,i-1);
        f+=0.5*pp2*square(r);
        fvec(jj++)=0.5*pp2*square(r);
        g(ii)+=pp2*r;
        g(ii-1)-=pp2*r;
        if (Hessflag)
        {
          H(ii,ii-1)=-pp2;
          H(ii-1,ii-1)+=pp2;
          H(ii,ii)+=pp2;
          if (H(ii,ii-1)==0)
          {
            cerr << "this cant happen" << endl;
            ad_exit(1);
          }
          if (pp2==0)
          {
            cerr << "this cant happen" << endl;
            ad_exit(1);
          }
        }
        ii++;
      }
    }
  }
  /*
  ofstream ofs("resids");
  ofs << resids << endl << endl;
  for (j=1;j<=num_fisheries;j++)
  {
    ofs << sort(resids(j)) << endl;
  }
  ofs << endl;
  for (j=1;j<=num_fisheries;j++)
  {
    ofs << sort(fabs(resids(j))) << endl;
  }
  */
  test_fpu();
  MY_DOUBLE_TYPE df=0.0;
  for (j=1;j<=jj-1;j++)
  {
    df+=fvec(j);
  }
#if (defined(USE_DD) && ( ((defined(__MSVC32__) &&  __MSVC32__ >=8)) || defined(linux) || defined(__ADMING__) || defined(_MAC)) )
  f=to_double(df);
#else
  f=df;
#endif
  //cout <<"fnorm.cpp " << df-f << endl;
}
