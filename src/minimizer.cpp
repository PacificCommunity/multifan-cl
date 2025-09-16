/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

/**
 * \file
 * Description not yet available.
 */
#include <sstream>
using std::istringstream;

#  include <admodel.h>
#include "phi_stuff.h"
#include "tridiagonal_dmatrix.h"
MY_DOUBLE_TYPE SUMPEN=0.0;
MY_DOUBLE_TYPE get_f(int mmin,int mmax,const dvector& p,const dvector& eta,
  const banded_lower_triangular_dmatrix & chinv,const MY_DOUBLE_TYPE N,
  const symmetric_tridiagonal_dmatrix& stdsinv,
  const dvector& q);

MY_DOUBLE_TYPE get_f(int mmin,int mmax,const dvector& p,const dvector& eta,
  const banded_lower_triangular_dmatrix & chinv,const MY_DOUBLE_TYPE N,
  const banded_symmetric_dmatrix& stdsinv,
  const dvector& q);

dvector get_g(int mmin,int mmax,const dvector& p,const dvector& eta,
  const banded_lower_triangular_dmatrix & chinv,const MY_DOUBLE_TYPE N,
  const symmetric_tridiagonal_dmatrix& stdsinv,
  const dvector& q);

dvector get_g(int mmin,int mmax,const dvector& p,const dvector& eta,
  const banded_lower_triangular_dmatrix & chinv,const MY_DOUBLE_TYPE N,
  const banded_symmetric_dmatrix& stdsinv,
  const dvector& q);

dvector get_fd_ders(int mmin,int mmax,const dvector& p,const dvector& eta,
  const banded_lower_triangular_dmatrix & chinv,const MY_DOUBLE_TYPE N,
  const symmetric_tridiagonal_dmatrix& stdsinv,
  const dvector& q,MY_DOUBLE_TYPE delta);

dvector get_fd_ders(int mmin,int mmax,const dvector& p,const dvector& eta,
  const banded_lower_triangular_dmatrix & chinv,const MY_DOUBLE_TYPE N,
  const banded_symmetric_dmatrix& stdsinv,
  const dvector& q,MY_DOUBLE_TYPE delta);


dvector operator * (const dvector&,const banded_lower_triangular_dmatrix&);
dmatrix make_dmatrix(const banded_lower_triangular_dmatrix&);

dvector inner_minimize(int n,MY_DOUBLE_TYPE N,const dvector & q,const dvector& p,
  const dmatrix& Sinv)
{
  fmmt1 fmc1(n,7);
  MY_DOUBLE_TYPE f=0.0;
  MY_DOUBLE_TYPE fb=1.e+100;
  dvector g(1,n);
  dvector ub(1,n);
  fmc1.itn=0;
  fmc1.ifn=0;
  fmc1.iprint=0;
  fmc1.ireturn=0;
  fmc1.ialph=0;
  fmc1.ihang=0;
  fmc1.ihflag=0;
  fmc1.use_control_c=0;
  
  fmc1.dfn=1.e-2;
  dvariable pen=0.0;
  dvector Nq=N*q;
  dvector eta(1,n);
  dvector t(1,n);
  eta.initialize();
  int mmin=1;
  int mmax=n;
  while (fmc1.ireturn>=0)
  {
    fmc1.fmin(f,eta,g);
    if (fmc1.ireturn>0)
    {
      MY_DOUBLE_TYPE v=0.0;
      for (int i=mmin;i<=mmax;i++)
      {
        // old lognormal version
        t(i)=phi(p(i),eta(i));
        v+=t(i);
      }
      f= N*log(v)-Nq*log(t)+0.5*(eta*(Sinv*eta));
      for (int i=mmin;i<=mmax;i++)
      {
        g(i)=(-N*q(i)/t(i)+N/v)*dphi_deta(p(i),eta(i));
      }
      g+=Sinv*eta;

      if (f<fb) 
      {
        fb=f;
        ub=eta;
      }
    }
  }
  eta=ub;
  return eta;
}
dvector inner_minimize(int n,MY_DOUBLE_TYPE N,const dvector & q,const dvector& p,
//  const dmatrix& Sinv,const dmatrix& cch,
  //const dmatrix& cch,
  const banded_lower_triangular_dmatrix & chinv)
{
  fmmt1 fmc1(n,7);
  MY_DOUBLE_TYPE f=0.0;
  MY_DOUBLE_TYPE fb=1.e+100;
  dvector g(1,n);
  dvector ub(1,n);
  fmc1.itn=0;
  fmc1.ifn=0;
  fmc1.iprint=0;
  fmc1.ireturn=0;
  fmc1.ialph=0;
  fmc1.ihang=0;
  fmc1.ihflag=0;
  fmc1.use_control_c=0;
  
  fmc1.dfn=1.e-2;
  dvariable pen=0.0;
  dvector Nq=N*q;
  dvector eps(1,n);
  dvector eta(1,n);
  dvector t(1,n);
  eps.initialize();
  eta.initialize();
  int mmin=1;
  int mmax=n;
  while (fmc1.ireturn>=0)
  {
    fmc1.fmin(f,eps,g);
    if (fmc1.ireturn>0)
    {
      //eta=cch*eps;
      eta=solve(chinv,eps);
      //cout << "RRR " << norm(eta-eta2) << " " << norm(eta) << endl;
      //ad_exit(81);
      MY_DOUBLE_TYPE v=0.0;
      for (int i=mmin;i<=mmax;i++)
      {
        // old lognormal version
        t(i)=phi(p(i),eta(i));
        v+=t(i);
      }
      f= N*log(v)-Nq*log(t)+0.5*norm2(eps);
      for (int i=mmin;i<=mmax;i++)
      {
        g(i)=(-N*q(i)/t(i)+N/v)*dphi_deta(p(i),eta(i));
      }
      g=solve_trans(chinv,g);
      //g=g*cch;
      //cout << "BBB " << norm2(g1-g) << " " << norm2(g1) << endl;
      //ad_exit(11);
      g+=eps;

      if (f<fb) 
      {
        fb=f;
        ub=eta;
      }
    }
  }
  eta=ub;
  return eta;
}
dvector inner_minimize(int n,MY_DOUBLE_TYPE N,const dvector & q,const dvector& p,
  const banded_lower_triangular_dmatrix & chinv, MY_DOUBLE_TYPE & gmax)
{
  fmmt1 fmc1(n,7);
  MY_DOUBLE_TYPE f=0.0;
  MY_DOUBLE_TYPE fb=1.e+100;
  dvector g(1,n);
  dvector ub(1,n);
  fmc1.itn=0;
  fmc1.ifn=0;
  fmc1.iprint=0;
  fmc1.ireturn=0;
  fmc1.ialph=0;
  fmc1.ihang=0;
  fmc1.ihflag=0;
  fmc1.use_control_c=0;
  
  fmc1.dfn=1.e-2;
  dvariable pen=0.0;
  dvector Nq=N*q;
  dvector eps(1,n);
  dvector eta(1,n);
  dvector t(1,n);
  eps.initialize();
  eta.initialize();
  int mmin=1;
  int mmax=n;
  while (fmc1.ireturn>=0)
  {
    fmc1.fmin(f,eps,g);
    if (fmc1.ireturn>0)
    {
      //eta=cch*eps;
      eta=solve(chinv,eps);
      //cout << "RRR " << norm(eta-eta2) << " " << norm(eta) << endl;
      //ad_exit(81);
      MY_DOUBLE_TYPE v=0.0;
      for (int i=mmin;i<=mmax;i++)
      {
        // old lognormal version
        t(i)=phi(p(i),eta(i));
        v+=t(i);
      }
      f= N*log(v)-Nq*log(t)+0.5*norm2(eps);
      for (int i=mmin;i<=mmax;i++)
      {
        g(i)=(-N*q(i)/t(i)+N/v)*dphi_deta(p(i),eta(i));
      }
      g=solve_trans(chinv,g);
      //g=g*cch;
      //cout << "BBB " << norm2(g1-g) << " " << norm2(g1) << endl;
      //ad_exit(11);
      g+=eps;

      if (f<fb) 
      {
        fb=f;
        ub=eta;
      }
    }
  }
  eta=ub;
  gmax=fmc1.gmax;
  return eta;
}
dvector inner_minimize(int n,MY_DOUBLE_TYPE N,const dvector & q,const dvector& p,
  const banded_lower_triangular_dmatrix & chinv, MY_DOUBLE_TYPE & gmax,dvector& eta,
  const symmetric_tridiagonal_dmatrix& stdsinv)
{
  fmmt1 fmc1(n,7);
  MY_DOUBLE_TYPE f=0.0;
  MY_DOUBLE_TYPE fb=1.e+100;
  dvector g(1,n);
  dvector ub(1,n);
  fmc1.itn=0;
  fmc1.ifn=0;
  fmc1.iprint=0;
  fmc1.ireturn=0;
  fmc1.ialph=0;
  fmc1.ihang=0;
  fmc1.ihflag=0;
  fmc1.use_control_c=0;
  
  fmc1.dfn=1.e-2;
  dvariable pen=0.0;
  dvector eps(1,n);
  // !!!!
  //eps=eta*chinv;
  eps=chinv*eta;
  int mmin=1;
  int mmax=n;
  while (fmc1.ireturn>=0)
  {
    fmc1.fmin(f,eps,g);
    // !!
    //eta=solve_trans(chinv,eps);
    eta=solve(chinv,eps);
    if (fmc1.ireturn>0)
    {
      f=get_f(mmin,mmax,p,eta,chinv,N,stdsinv,q);
      //if (fmc1.ifn==0) cout << f << endl;
      g=get_g(mmin,mmax,p,eta,chinv,N,stdsinv,q);
      //  cout << 
      //  norm(g-get_fd_ders(mmin,mmax,p,eta,chinv,N,stdsinv,q,1.e-6))<< endl;
      g=solve_trans(chinv,g);
      if (f<fb) 
      {
        fb=f;
        ub=eta;
      }
    }
  }
  eta=ub;
  gmax=fmc1.gmax;
  //cout << f << endl;
  return eta;
}
dvector inner_minimize(int n,MY_DOUBLE_TYPE N,const dvector & q,const dvector& p,
  const banded_lower_triangular_dmatrix & chinv, MY_DOUBLE_TYPE & gmax,dvector& eta,
  const banded_symmetric_dmatrix& stdsinv)
{
  fmmt1 fmc1(n,7);
  MY_DOUBLE_TYPE f=0.0;
  MY_DOUBLE_TYPE fb=1.e+100;
  dvector g(1,n);
  dvector ub(1,n);
  fmc1.itn=0;
  fmc1.ifn=0;
  fmc1.iprint=0;
  fmc1.ireturn=0;
  fmc1.ialph=0;
  fmc1.ihang=0;
  fmc1.ihflag=0;
  fmc1.use_control_c=0;
  
  fmc1.dfn=1.e-2;
  dvariable pen=0.0;
  dvector eps(1,n);
  // !!!!
  //eps=eta*chinv;
  eps=chinv*eta;
  int mmin=1;
  int mmax=n;
  while (fmc1.ireturn>=0)
  {
    fmc1.fmin(f,eps,g);
    // !!
    //eta=solve_trans(chinv,eps);
    eta=solve(chinv,eps);
    if (fmc1.ireturn>0)
    {
      f=get_f(mmin,mmax,p,eta,chinv,N,stdsinv,q);
      //if (fmc1.ifn==0) cout << f << endl;
      g=get_g(mmin,mmax,p,eta,chinv,N,stdsinv,q);
        //cout << 
        //norm(g-get_fd_ders(mmin,mmax,p,eta,chinv,N,stdsinv,q,1.e-6))<< endl;
      g=solve_trans(chinv,g);
      if (f<fb) 
      {
        fb=f;
        ub=eta;
      }
    }
  }
  eta=ub;
  gmax=fmc1.gmax;
  //cout << f << endl;
  return eta;
}
