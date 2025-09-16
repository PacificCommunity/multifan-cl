/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#define HOME_VERSION
#include "all.hpp"
//uostream& operator<<(uostream& ostr,const dvector& z);
extern adstring * _current_path;

void getpaths(adstring& current_path,adstring& root);
#include "variable.hpp"
  //double fcomp(const dvar_len_fish_stock_history& fsh,dvar_vector& x,int nvar,
    //int print_switch,const dvector& gbest,int grad_switch);
//  MY_DOUBLE_TYPE fcomp(const dvar_len_fish_stock_history& _fsh,const dvar_vector& _x, int nvar,int print_switch,const dvector& _gbest,int gradient_switch);
MY_DOUBLE_TYPE useless(MY_DOUBLE_TYPE u){return 0.0;}
#if defined(_WIN32)
#  if defined(close)
#    undef close
#  endif
#endif


void hess_routine(independent_variables& x,dvar_len_fish_stock_history& fsh,
  int nvar)
{
  gradient_structure::set_YES_DERIVATIVES();
  //fmc.ireturn=fsh.parest_flags(1);
  int extra_precision=1;  
  MY_DOUBLE_TYPE f;
  MY_DOUBLE_TYPE eps=0.5;
  MY_DOUBLE_TYPE delta=1.0e-6;
  if (extra_precision) delta=1.0e-5;   // smaller step size
  MY_DOUBLE_TYPE sdelta1;
  MY_DOUBLE_TYPE sdelta2;
  //MY_DOUBLE_TYPE volatile sdelta1;
  //MY_DOUBLE_TYPE volatile sdelta2;
  dvector g1(1,nvar);
  dvector g2(1,nvar);
  dvector gbest(1,nvar);
  dvector hess(1,nvar);
  dvector hess1(1,nvar);
  gbest.fill_seqadd(1.e+50,0.);

  adstring current_path;
  adstring root;
  //getpaths(current_path,root);
/*
#if !defined(linux)
  adstring hessfile=current_path + adstring("\\") + root + adstring(".hes"); 
#else
  adstring hessfile=current_path + adstring("/") + root + adstring(".hes"); 
#endif
*/
  root=ad_root;
  int hswitch=fsh.parest_flags(223);
  int hmin=1;
  int hmax=nvar;
  if (hswitch)
  {
    hmin=fsh.parest_flags(223);
    hmax=fsh.parest_flags(224);
  } 
  adstring hessfile=root + adstring(".hes"); 
  uostream ofs(hessfile);
  if (!ofs)
  {
    cerr << "Error opening file " << (char*) hessfile << endl;
    exit(1);
  }
  f=fcomp(fsh,dvar_vector(x),nvar,0,gbest,0);
  gradcalc(nvar,g1);
 // ofs << setfixed()<< setprecision(18) << g1  << endl << endl;
  MY_DOUBLE_TYPE xsave;
  //ofs << "  " << nvar << endl;
  if (hswitch) 
    ofs << nvar << "  " << hmin << "  " << hmax;
  else
    ofs <<  "  " << nvar;

  {
    for (int i=hmin;i<=hmax;i++)
    {
      cout << "Estimating row " << i << " out of " << nvar 
           << " for hessian" << endl;

      xsave=x(i);
      sdelta1=xsave+delta;
      useless(sdelta1);
      sdelta1-=xsave;
      x(i)=xsave+sdelta1;
      f=fcomp(fsh,dvar_vector(x),nvar,0,gbest,0);
      cout << " Calling gradcalc " << endl;
      gradcalc(nvar,g1);
      cout << " Finished gradcalc " << endl;

      sdelta2=xsave-delta;
      useless(sdelta2);
      sdelta2-=xsave;
      x(i)=xsave+sdelta2;
      f=fcomp(fsh,dvar_vector(x),nvar,0,gbest,0);
      cout << " Calling gradcalc " << endl;
      gradcalc(nvar,g2);
      cout << " Finished gradcalc " << endl;
      x(i)=xsave;

      hess=(g1-g2)/(sdelta1-sdelta2);

      if (extra_precision)
      {

        MY_DOUBLE_TYPE delta1=eps*delta;
        xsave=x(i);
        sdelta1=xsave+delta1;
        useless(sdelta1);
        sdelta1-=xsave;
        x(i)=xsave+sdelta1;
        f=fcomp(fsh,dvar_vector(x),nvar,0,gbest,0);
        cout << " Calling gradcalc " << endl;
        gradcalc(nvar,g1);
        cout << " Finished gradcalc " << endl;
  
        sdelta2=xsave-delta1;
        useless(sdelta2);
        sdelta2-=xsave;
        x(i)=xsave+sdelta2;
        f=fcomp(fsh,dvar_vector(x),nvar,0,gbest,0);
        cout << " Calling gradcalc " << endl;
        gradcalc(nvar,g2);
        cout << " Finished gradcalc " << endl;
        x(i)=xsave;
  
        hess1=(g1-g2)/(sdelta1-sdelta2);
        MY_DOUBLE_TYPE eps2=eps*eps;
        hess=(eps2*hess-hess1) /(eps2-1.);

      }

      ofs <<  hess ;
    }
  } 
  ofs.CLOSE();
  {
    uistream ifs(hessfile);
    int nvar;
    ifs >> nvar;
    dmatrix h(1,nvar,1,nvar);
    for (int i=1;i<=nvar;i++)
    {
      ifs >> h(i); 
    }
    MY_DOUBLE_TYPE fmax=0.0;
    for (int i=1;i<=nvar;i++)
    {
      for (int j=1;j<=nvar;j++)
      {
        if (fabs(h(i,j)-h(j,i))> fmax)
          fmax=fabs(h(i,j)-h(j,i));
      }
    }
    cout << "fmax = " << setprecision(15) << fmax << endl;
  }
  gradient_structure::set_NO_DERIVATIVES();
}

#undef HOME_VERSION
