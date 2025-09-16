/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#define USE_DD_NOT
//#define USE_DD_NOTOUBLE
#include "all.hpp"
//#include <admodel.h>
#include "df22fun.h"
#include "phi_stuff.h"
extern int new_flag;
#include "tridiagonal_dmatrix.h"
df2_two_variable d2phi(MY_DOUBLE_TYPE _p,MY_DOUBLE_TYPE _eta);
dvector operator * (const dvector&,const banded_lower_triangular_dmatrix&);
dvector solve(const symmetric_tridiagonal_dmatrix& _STD1,const dvector& _w,
  const dvector& _v);

dvector inner_minimize(int,MY_DOUBLE_TYPE, const dvector&, const dvector&, const dmatrix&, const banded_lower_triangular_dmatrix&);

dvector inner_minimize(int,MY_DOUBLE_TYPE, const dvector&, const dvector&,const banded_lower_triangular_dmatrix&);

dvector inner_minimize(int n,MY_DOUBLE_TYPE N,const dvector & q,const dvector& p,
  const banded_lower_triangular_dmatrix & chinv, MY_DOUBLE_TYPE & gmax,dvector& eta,
  const symmetric_tridiagonal_dmatrix& stdsinv);

dvector inner_minimize(int n,MY_DOUBLE_TYPE N,const dvector & q,const dvector& p,
  const banded_lower_triangular_dmatrix & chinv, MY_DOUBLE_TYPE & gmax,dvector& eta,
  const banded_symmetric_dmatrix& stdsinv);

//dvector inner_minimize(int n,MY_DOUBLE_TYPE N,const dvector & q,const dvector& p,
//  const dmatrix& Sinv);
//dvector inner_minimize(int n,MY_DOUBLE_TYPE N,const dvector & q,const dvector& p,
//  const dmatrix& Sinv,const dmatrix& S);
extern int my_debug_flag;

dvector inner_minimize(int n,const MY_DOUBLE_TYPE & N,const dvector & q,
  const dvector& p,
  const dmatrix& Sinv,const dmatrix& cch,const banded_lower_triangular_dmatrix & chinv);

dvector inner_minimize(int n,MY_DOUBLE_TYPE N,const dvector & q,const dvector& p,
  const banded_lower_triangular_dmatrix & chinv, MY_DOUBLE_TYPE & gmax,dvector& eta);

dvector inner_minimize(int n,MY_DOUBLE_TYPE N,const dvector & q,const dvector& p,
  const banded_lower_triangular_dmatrix & chinv, MY_DOUBLE_TYPE & gmax);

dvector sherman_morrison_woodbury(symmetric_tridiagonal_dmatrix& BHess,
  dvector& grad,dvector & w)
{
  // ********************************************
  // rank one index update
  // ********************************************
  dvector h=solve(BHess,grad);
  dvector b=solve(BHess,w);
  h = h + (w*h)*b/(1.0 - w*b);
  return h;
}
dvector solve(const banded_symmetric_dmatrix& _S,const dmatrix& u,
  const dmatrix& v,const dvector& x);


dvector solve(const symmetric_tridiagonal_dmatrix& _STD1,const dvector& _w,
  const dvector& _v)
{
  ADUNCONST(symmetric_tridiagonal_dmatrix,STD1)
  ADUNCONST(dvector,w)
  ADUNCONST(dvector,v)
  dvector u1=solve(STD1,v);
  dvector x=solve(STD1,w);
  //cout << make_dmatrix(STD)*u1-v << endl;
  //cout << make_dmatrix(STD)*x-w << endl;
  //ad_exit(1);
  //dvector ans = u1 - (w*u1)*x / ( 1 + w*x);
  dvector ans = u1 + (w*u1)*x / ( 1 - w*x);
  return ans;
}

dvector solve(const banded_symmetric_dmatrix& _STD1,const dvector& _w,
  const dvector& _v)
{
  ADUNCONST(banded_symmetric_dmatrix,STD1)
  ADUNCONST(dvector,w)
  ADUNCONST(dvector,v)
  dvector u1=solve(STD1,v);
  dvector x=solve(STD1,w);
  //cout << make_dmatrix(STD)*u1-v << endl;
  //cout << make_dmatrix(STD)*x-w << endl;
  //ad_exit(1);
  //dvector ans = u1 - (w*u1)*x / ( 1 + w*x);
  dvector ans = u1 + (w*u1)*x / ( 1 - w*x);
  return ans;
}


extern MY_DOUBLE_TYPE SUMPEN;
MY_DOUBLE_TYPE get_f(int mmin,int mmax,const dvector& p,const dvector& eta,
  const banded_lower_triangular_dmatrix & chinv,const MY_DOUBLE_TYPE N,
  const symmetric_tridiagonal_dmatrix& stdsinv,
  const dvector& q)
{
  dvector Nq=N*q;
  MY_DOUBLE_TYPE v=0.0;
  dvector t(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    // old lognormal version
    t(i)=phi(p(i),eta(i));
    v+=t(i);
  }
  int n=mmax-mmin+1;
  MY_DOUBLE_TYPE f= N*log(v)-Nq*log(t);
  f+=0.5*(eta*(stdsinv*eta));
  return f;
}

MY_DOUBLE_TYPE get_f(int mmin,int mmax,const dvector& p,const dvector& eta,
  const banded_lower_triangular_dmatrix & chinv,const MY_DOUBLE_TYPE N,
  const banded_symmetric_dmatrix& stdsinv,const dvector& q)
{
  dvector Nq=N*q;
  MY_DOUBLE_TYPE v=0.0;
  dvector t(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    // old lognormal version
    t(i)=phi(p(i),eta(i));
    v+=t(i);
  }
  int n=mmax-mmin+1;
  MY_DOUBLE_TYPE f= N*log(v)-Nq*log(t);
  f+=0.5*(eta*(stdsinv*eta));
  return f;
}

dvector get_fd_ders(int mmin,int mmax,const dvector& p,const dvector& eta,
  const banded_lower_triangular_dmatrix & chinv,const MY_DOUBLE_TYPE N,
  const symmetric_tridiagonal_dmatrix& stdsinv,
  const dvector& q,MY_DOUBLE_TYPE delta)
{
  
  dvector ee(mmin,mmax);
  ee=eta;
  dvector gg(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    MY_DOUBLE_TYPE esave=ee(i);
    ee(i)=esave+delta;
    MY_DOUBLE_TYPE fp=get_f(mmin,mmax,p,ee,chinv,N,stdsinv,q);
    ee(i)=esave-delta;
    MY_DOUBLE_TYPE fm=get_f(mmin,mmax,p,ee,chinv,N,stdsinv,q);
    gg(i)=(fp-fm)/(2.0*delta);
    ee(i)=esave;
  }
  return gg;
}

dvector get_fd_ders(int mmin,int mmax,const dvector& p,const dvector& eta,
  const banded_lower_triangular_dmatrix & chinv,const MY_DOUBLE_TYPE N,
  const banded_symmetric_dmatrix& stdsinv,
  const dvector& q,MY_DOUBLE_TYPE delta)
{
  
  dvector ee(mmin,mmax);
  ee=eta;
  dvector gg(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    MY_DOUBLE_TYPE esave=ee(i);
    ee(i)=esave+delta;
    MY_DOUBLE_TYPE fp=get_f(mmin,mmax,p,ee,chinv,N,stdsinv,q);
    ee(i)=esave-delta;
    MY_DOUBLE_TYPE fm=get_f(mmin,mmax,p,ee,chinv,N,stdsinv,q);
    gg(i)=(fp-fm)/(2.0*delta);
    ee(i)=esave;
  }
  return gg;
}

dmatrix get_fd_H(int mmin,int mmax,const dvector& p,const dvector& eta,
  const banded_lower_triangular_dmatrix & chinv,const MY_DOUBLE_TYPE N,
  const symmetric_tridiagonal_dmatrix& stdsinv,
  const dvector& q)
{
  
  dvector ee(mmin,mmax);
  ee=eta;
  dvector gg(mmin,mmax);
  dmatrix H(mmin,mmax,mmin,mmax);
  H.initialize();
  //double delta=1.e-5;
  MY_DOUBLE_TYPE delta=5.e-6;
  //double delta=1.e-6;
  for (int i=mmin;i<=mmax;i++)
  {
    MY_DOUBLE_TYPE esave=ee(i);
    ee(i)=esave+delta;
    dvector Hu=get_fd_ders(mmin,mmax,p,ee,chinv,N,stdsinv,q,delta);
    ee(i)=esave-delta;
    dvector Hm=get_fd_ders(mmin,mmax,p,ee,chinv,N,stdsinv,q,delta);
    H(i)=(Hu-Hm)/(2.0*delta);
    ee(i)=esave;
  }
  return H;
}

dmatrix get_fd_H(int mmin,int mmax,const dvector& p,const dvector& eta,
  const banded_lower_triangular_dmatrix & chinv,const MY_DOUBLE_TYPE N,
  const banded_symmetric_dmatrix& stdsinv,
  const dvector& q)
{
  
  dvector ee(mmin,mmax);
  ee=eta;
  dvector gg(mmin,mmax);
  dmatrix H(mmin,mmax,mmin,mmax);
  H.initialize();
  MY_DOUBLE_TYPE delta=2.e-5;
  for (int i=mmin;i<=mmax;i++)
  {
    MY_DOUBLE_TYPE esave=ee(i);
    ee(i)=esave+delta;
    dvector Hu=get_fd_ders(mmin,mmax,p,ee,chinv,N,stdsinv,q,delta);
    ee(i)=esave-delta;
    dvector Hm=get_fd_ders(mmin,mmax,p,ee,chinv,N,stdsinv,q,delta);
    H(i)=(Hu-Hm)/(2.0*delta);
    ee(i)=esave;
  }
  return H;
}


dvector get_g(int mmin,int mmax,const dvector& p,const dvector& eta,
  const banded_lower_triangular_dmatrix & chinv,const MY_DOUBLE_TYPE N,
  const symmetric_tridiagonal_dmatrix& stdsinv,
  const dvector& q)
{
  dvector Nq=N*q;
  MY_DOUBLE_TYPE v=0.0;
  dvector t(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    // old lognormal version
    t(i)=phi(p(i),eta(i));
    v+=t(i);
  }
  int n=mmax-mmin+1;
  MY_DOUBLE_TYPE f= N*log(v)-Nq*log(t);
  dvector grad(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    grad(i)=(-N*q(i)/t(i)+N/v)*dphi_deta(p(i),eta(i));
  }
  grad+=stdsinv*eta;
  return grad;
}

dvector get_g(int mmin,int mmax,const dvector& p,const dvector& eta,
  const banded_lower_triangular_dmatrix & chinv,const MY_DOUBLE_TYPE N,
  const banded_symmetric_dmatrix& stdsinv,
  const dvector& q)
{
  dvector Nq=N*q;
  MY_DOUBLE_TYPE v=0.0;
  dvector t(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    // old lognormal version
    t(i)=phi(p(i),eta(i));
    v+=t(i);
  }
  int n=mmax-mmin+1;
  MY_DOUBLE_TYPE f= N*log(v)-Nq*log(t);
  dvector grad(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    grad(i)=(-N*q(i)/t(i)+N/v)*dphi_deta(p(i),eta(i));
  }
  grad+=stdsinv*eta;
  return grad;
}

symmetric_tridiagonal_dmatrix get_H(int mmin,int mmax,const dvector& p,
  const dvector& eta,
  const banded_lower_triangular_dmatrix & chinv,const MY_DOUBLE_TYPE N,
  const symmetric_tridiagonal_dmatrix& _stdsinv,
  const dvector& q,dvector& w)
{
  symmetric_tridiagonal_dmatrix BHess(mmin,mmax);
  ADUNCONST(symmetric_tridiagonal_dmatrix,stdsinv)
  BHess.initialize();
  
  dvector Nq=N*q;
  MY_DOUBLE_TYPE v=0.0;
  dvector t(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    // old lognormal version
    t(i)=phi(p(i),eta(i));
    v+=t(i);
  }
  int n=mmax-mmin+1;
  dvector eps=eta*chinv;
 /*
  dvector e(1,n);
  e=1.0/sqrt(n);
  dvector ff=e*chinv;
  dvector h=solve(chinv,e);
  dvector se=stdsinv*e;
  dvector v2=stdsinv*eta;
 */
  dvector d1(mmin,mmax);
  dvector d2(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    df2_two_variable tmp=d2phi(p(i),eta(i));
    d1(i)=tmp.get_u_y();
    d2(i)=tmp.get_u_yy();
  }
  MY_DOUBLE_TYPE Nv=N/v;
  MY_DOUBLE_TYPE Nvv=Nv/v;
  dvector Nqt(mmin,mmax);
  //dvector w(mmin,mmax);
  dvector Nqtt(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    Nq(i)=N*q(i);
  }
  for (int i=mmin;i<=mmax;i++)
  {
    Nqt(i)=Nq(i)/t(i);
  }
  for (int i=mmin;i<=mmax;i++)
  {
    Nqtt(i)=Nqt(i)/t(i);
  }
     
  for (int i=mmin;i<=mmax;i++)
  {
    BHess(i,i)+=(-Nqt(i)+Nv)*d2(i);
    BHess(i,i)+=Nqtt(i)*square(d1(i));
    // for rank one update
    w(i)=sqrt(Nvv)*d1(i);
  }
  
  for (int i=mmin;i<mmax;i++)
  {
    BHess(i,i)+=stdsinv(i,i);
    BHess(i+1,i)=stdsinv(i+1,i);
  }
  BHess(mmax,mmax)+=stdsinv(mmax,mmax);
  
  return BHess;
}

symmetric_tridiagonal_dmatrix get_H(int mmin,int mmax,const dvector& p,
  const dvector& eta,
  const banded_lower_triangular_dmatrix & chinv,const MY_DOUBLE_TYPE N,
  const banded_symmetric_dmatrix& _stdsinv,
  const dvector& q,dvector& w)
{
  symmetric_tridiagonal_dmatrix BHess(mmin,mmax);
  ADUNCONST(symmetric_tridiagonal_dmatrix,stdsinv)
  BHess.initialize();
  
  dvector Nq=N*q;
  MY_DOUBLE_TYPE v=0.0;
  dvector t(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    // old lognormal version
    t(i)=phi(p(i),eta(i));
    v+=t(i);
  }
  int n=mmax-mmin+1;
  dvector eps=eta*chinv;
 /*
  dvector e(1,n);
  e=1.0/sqrt(n);
  dvector ff=e*chinv;
  dvector h=solve(chinv,e);
  dvector se=stdsinv*e;
  dvector v2=stdsinv*eta;
 */
  dvector d1(mmin,mmax);
  dvector d2(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    df2_two_variable tmp=d2phi(p(i),eta(i));
    d1(i)=tmp.get_u_y();
    d2(i)=tmp.get_u_yy();
  }
  MY_DOUBLE_TYPE Nv=N/v;
  MY_DOUBLE_TYPE Nvv=Nv/v;
  dvector Nqt(mmin,mmax);
  //dvector w(mmin,mmax);
  dvector Nqtt(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    Nq(i)=N*q(i);
  }
  for (int i=mmin;i<=mmax;i++)
  {
    Nqt(i)=Nq(i)/t(i);
  }
  for (int i=mmin;i<=mmax;i++)
  {
    Nqtt(i)=Nqt(i)/t(i);
  }
     
  for (int i=mmin;i<=mmax;i++)
  {
    BHess(i,i)+=(-Nqt(i)+Nv)*d2(i);
    BHess(i,i)+=Nqtt(i)*square(d1(i));
    // for rank one update
    w(i)=sqrt(Nvv)*d1(i);
  }
  
  for (int i=mmin;i<mmax;i++)
  {
    BHess(i,i)+=stdsinv(i,i);
    BHess(i+1,i)=stdsinv(i+1,i);
  }
  BHess(mmax,mmax)+=stdsinv(mmax,mmax);
  
  return BHess;
}

symmetric_tridiagonal_dmatrix get_H(int mmin,int mmax,const dvector& p,
  const dvector& eta,
  const banded_lower_triangular_dmatrix & chinv,const MY_DOUBLE_TYPE N,
  const symmetric_tridiagonal_dmatrix& _stdsinv,
  const dvector& q,dmatrix& U,dmatrix& V)
{
  symmetric_tridiagonal_dmatrix BHess(mmin,mmax);
  ADUNCONST(symmetric_tridiagonal_dmatrix,stdsinv)
  BHess.initialize();
  
  
  dvector Nq=N*q;
  MY_DOUBLE_TYPE v=0.0;
  dvector t(mmin,mmax);
  dvector w(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    // old lognormal version
    t(i)=phi(p(i),eta(i));
    v+=t(i);
  }
  int n=mmax-mmin+1;
  dvector eps=eta*chinv;
  dvector e(1,n);
  e=1.0/sqrt(n);
  dvector se=stdsinv*e;
 /*
  dvector ff=e*chinv;
  dvector h=solve(chinv,e);
  dvector v2=stdsinv*eta;
 */
  dvector d1(mmin,mmax);
  dvector d2(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    df2_two_variable tmp=d2phi(p(i),eta(i));
    d1(i)=tmp.get_u_y();
    d2(i)=tmp.get_u_yy();
  }
  MY_DOUBLE_TYPE Nv=N/v;
  MY_DOUBLE_TYPE Nvv=Nv/v;
  dvector Nqt(mmin,mmax);
  //dvector w(mmin,mmax);
  dvector Nqtt(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    Nq(i)=N*q(i);
  }
  for (int i=mmin;i<=mmax;i++)
  {
    Nqt(i)=Nq(i)/t(i);
  }
  for (int i=mmin;i<=mmax;i++)
  {
    Nqtt(i)=Nqt(i)/t(i);
  }
     
  for (int i=mmin;i<=mmax;i++)
  {
    BHess(i,i)+=(-Nqt(i)+Nv)*d2(i);
    BHess(i,i)+=Nqtt(i)*square(d1(i));
    // for rank one update
    w(i)=sqrt(Nvv)*d1(i);
  }
  U(1)=-w;
  V(1)=w;
  //U(2)=-e;
  //V(2)=se;
  //U(3)=e;
  //V(3)=e*(e*se);
  
  for (int i=mmin;i<mmax;i++)
  {
    BHess(i,i)+=stdsinv(i,i);
    BHess(i+1,i)=stdsinv(i+1,i);
  }
  BHess(mmax,mmax)+=stdsinv(mmax,mmax);
  
  return BHess;
}

banded_symmetric_dmatrix get_H(int mmin,int mmax,const dvector& p,
  const dvector& eta,
  const banded_lower_triangular_dmatrix & chinv,const MY_DOUBLE_TYPE N,
  const banded_symmetric_dmatrix& _stdsinv,
  const dvector& q,dmatrix& U,dmatrix& V)
{
  ADUNCONST(banded_symmetric_dmatrix,stdsinv)
  int bw=stdsinv.bandwidth();
  banded_symmetric_dmatrix BHess(mmin,mmax,bw);
  BHess.initialize();
  
  
  dvector Nq=N*q;
  MY_DOUBLE_TYPE v=0.0;
  dvector t(mmin,mmax);
  dvector w(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    // old lognormal version
    t(i)=phi(p(i),eta(i));
    v+=t(i);
  }
  int n=mmax-mmin+1;
  dvector eps=eta*chinv;
  dvector e(1,n);
  e=1.0/sqrt(n);
  dvector se=stdsinv*e;
 /*
  dvector ff=e*chinv;
  dvector h=solve(chinv,e);
  dvector v2=stdsinv*eta;
 */
  dvector d1(mmin,mmax);
  dvector d2(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    df2_two_variable tmp=d2phi(p(i),eta(i));
    d1(i)=tmp.get_u_y();
    d2(i)=tmp.get_u_yy();
  }
  MY_DOUBLE_TYPE Nv=N/v;
  MY_DOUBLE_TYPE Nvv=Nv/v;
  dvector Nqt(mmin,mmax);
  //dvector w(mmin,mmax);
  dvector Nqtt(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    Nq(i)=N*q(i);
  }
  for (int i=mmin;i<=mmax;i++)
  {
    Nqt(i)=Nq(i)/t(i);
  }
  for (int i=mmin;i<=mmax;i++)
  {
    Nqtt(i)=Nqt(i)/t(i);
  }
     
  for (int i=mmin;i<=mmax;i++)
  {
    BHess(i,i)+=(-Nqt(i)+Nv)*d2(i);
    BHess(i,i)+=Nqtt(i)*square(d1(i));
    // for rank one update
    w(i)=sqrt(Nvv)*d1(i);
  }
  U(1)=-w;
  V(1)=w;
  //U(2)=-e;
  //V(2)=se;
  //U(3)=e;
  //V(3)=e*(e*se);
  
  for (int i=mmin;i<=mmax;i++)
  {
    for (int j=i;j<=min(mmax,i+bw-1);j++)
    {
      BHess(j,i)+=stdsinv(j,i);
    }
  }
  //BHess(mmax,mmax)+=stdsinv(mmax,mmax);
  
  return BHess;
}


dvector mynr(int & badflag,
  int & negflag,
  const MY_DOUBLE_TYPE N,
  const dvector& q,
  const dvector& p,
  const dvector& eta,
  const banded_lower_triangular_dmatrix & chinv,
  const symmetric_tridiagonal_dmatrix& stdsinv);

dvector lognormal_ss_multinomial_newton_raphson(const MY_DOUBLE_TYPE N,
  const dvector& q,const dvector& p,
  //const dmatrix & cch,
  const banded_lower_triangular_dmatrix & chinv,
  const symmetric_tridiagonal_dmatrix& _stdsinv,const dvector& _ee)
{
  ADUNCONST(dvector,ee)
  ADUNCONST(symmetric_tridiagonal_dmatrix,stdsinv)
  int mmin=q.indexmin();
  int mmax=q.indexmax();
  // Newton Raphson
  
  dvector t(mmin,mmax);
  dvector Nq=N*q;
  dvector grad(mmin,mmax);
  symmetric_tridiagonal_dmatrix BHess(mmin,mmax);
  dvector w(mmin,mmax);
  // get initial norm grad

  MY_DOUBLE_TYPE oldnc=1.e+100;
  int goodflag=1;
  int itcount=0;
  dvector h(mmin,mmax);
  dvector h3(mmin,mmax);
  dvector e(mmin,mmax);
  dvector eta(mmin,mmax);
  e=1.0;
  h.initialize();
  h3.initialize();
  MY_DOUBLE_TYPE mult=1.0;
  MY_DOUBLE_TYPE gmax=0.0;
  dvector eta_inner=inner_minimize(mmax,N,q,p,chinv,gmax,ee,stdsinv);
  if (fabs(gmax)>1.e-3)
  {
     eta_inner=inner_minimize(mmax,N,q,p,chinv,gmax,eta_inner,stdsinv);
  }
  if (fabs(gmax)>1.e-1)
  {
     eta_inner=inner_minimize(mmax,N,q,p,chinv,gmax,eta_inner,stdsinv);
  }
  eta=eta_inner;
 

  //cout << "norm2(eta)" <<  norm2(eta) << endl;
  //exit(1);
  do
  {
    MY_DOUBLE_TYPE f=get_f(mmin,mmax,p,eta,chinv,N,stdsinv,q);
    dvector grad=get_g(mmin,mmax,p,eta,chinv,N,stdsinv,q);
    //dvector gg1=get_fd_ders(mmin,mmax,p,eta,chinv,N,stdsinv,q,1.e-6);
    //cout << norm(grad-gg1) << endl;
    //ad_exit(1);

    MY_DOUBLE_TYPE nc=norm(grad);
    if (nc <1.e-12 || itcount > 8 ||
      nc <1.e-8 && itcount > 10)
    {
      if (itcount >15)
      {
        cout << "A convergence difficulty in newtraph " << nc << endl;
#if !defined(NO_MY_DOUBLE_TYPE)
        if (nc>1.0L)
#else
        if (nc>1.0)
#endif
        {
          cout << "trapped" << endl;
        }
      }
      // need U V  for rank three update
      dmatrix U(1,1,mmin,mmax);
      dmatrix V(1,1,mmin,mmax);
      BHess=get_H(mmin,mmax,p,eta,chinv,N,stdsinv,q,U,V);
      break;
    }
    dvector vvv(mmin,mmax);
    if (nc<oldnc)
    {
      // need U V  for rank three update
      dmatrix U(1,1,mmin,mmax);
      dmatrix V(1,1,mmin,mmax);
      BHess=get_H(mmin,mmax,p,eta,chinv,N,stdsinv,q,U,V);
      //dmatrix BHess1=make_dmatrix(BHess);
      //int m1=U.indexmax();
      //for (int i=1;i<=m1;i++)
      //{
      //  BHess1+=outer_prod(U(i),V(i));
      //}
      
      //dmatrix Hess=get_fd_H(mmin,mmax,p,eta,chinv,N,stdsinv,q);
      //cout << norm2(BHess1-Hess) << endl;
      //ad_exit(1);
      
      // ********************************************
      // rank one index update
      // ********************************************
      dvector h3=solve(BHess,U,V,grad);
      //dvector hk=solve(BHess1,grad);
      //dvector hh1=solve(Hess,grad);
      //cout << norm(h3-hh1) << endl;
      //cout << norm(hk-hh1) << endl;
      //cout << norm(hk-h3) << endl;
      eta-=h3;
      itcount++;
      oldnc=nc;
      mult=1.0;
    }
    else
    {
      mult*=0.5;
      eta+=mult*h3;
      //eta+=mult*vvv;
      itcount++;
    }
  }
  while(1);
  //cout << endl;
  /*
  int badflag=0;
  int negflag=0;
  eta=mynr(badflag,negflag,N,q,p,eta,chinv,stdsinv);
  while(badflag || negflag)
  {
    eta=mynr(badflag,negflag,N,q,p,eta,chinv,stdsinv);
  }
  //cout << "XX " << norm2(eta1-eta) << endl;
  */
  return eta;
}

dvector lognormal_ss_multinomial_newton_raphson(const MY_DOUBLE_TYPE N,
  const dvector& q,const dvector& p,
  const banded_lower_triangular_dmatrix & chinv,
  const banded_symmetric_dmatrix& _stdsinv,const dvector& _ee)
{
  ADUNCONST(banded_symmetric_dmatrix,stdsinv)
  int bw=stdsinv.bandwidth();
  ADUNCONST(dvector,ee)
  int mmin=q.indexmin();
  int mmax=q.indexmax();
  // Newton Raphson
  
  dvector t(mmin,mmax);
  dvector Nq=N*q;
  dvector grad(mmin,mmax);
  banded_symmetric_dmatrix BHess(mmin,mmax,bw);
  dvector w(mmin,mmax);
  // get initial norm grad

  MY_DOUBLE_TYPE oldnc=1.e+100;
  int goodflag=1;
  int itcount=0;
  dvector h(mmin,mmax);
  dvector h3(mmin,mmax);
  dvector e(mmin,mmax);
  dvector eta(mmin,mmax);
  e=1.0;
  h.initialize();
  h3.initialize();
  MY_DOUBLE_TYPE mult=1.0;
  MY_DOUBLE_TYPE gmax=0.0;
  dvector eta_inner=inner_minimize(mmax,N,q,p,chinv,gmax,ee,stdsinv);
  if (fabs(gmax)>1.e-3)
  {
     eta_inner=inner_minimize(mmax,N,q,p,chinv,gmax,eta_inner,stdsinv);
  }
  if (fabs(gmax)>1.e-1)
  {
     eta_inner=inner_minimize(mmax,N,q,p,chinv,gmax,eta_inner,stdsinv);
  }
  eta=eta_inner;
 

  //cout << "norm2(eta)" <<  norm2(eta) << endl;
  //exit(1);
  do
  {
    MY_DOUBLE_TYPE f=get_f(mmin,mmax,p,eta,chinv,N,stdsinv,q);
    dvector grad=get_g(mmin,mmax,p,eta,chinv,N,stdsinv,q);
    //dvector gg1=get_fd_ders(mmin,mmax,p,eta,chinv,N,stdsinv,q,1.e-6);
    //cout << norm(grad-gg1) << endl;
    //ad_exit(1);

    MY_DOUBLE_TYPE nc=norm(grad);
    if (nc <1.e-12 || itcount > 12 ||
      (nc <1.e-8 && itcount > 5) ||
      (nc <1.e-7 && itcount > 7)  ||
      (nc <1.e-6 && itcount > 8)) 
    {
      if (itcount >10)
      {
        cout << "B convergence difficulty in newtraph " << nc << endl;
#if !defined(NO_MY_DOUBLE_TYPE)
        if (nc>1.0L)
#else
        if (nc>1.0)
#endif
        {
          cout << "trapped" << endl;
        }
      }
      // need U V  for rank three update
      dmatrix U(1,1,mmin,mmax);
      dmatrix V(1,1,mmin,mmax);
      BHess=get_H(mmin,mmax,p,eta,chinv,N,stdsinv,q,U,V);
      break;
    }
    dvector vvv(mmin,mmax);
    if (nc<oldnc)
    {
      // need U V  for rank three update
      dmatrix U(1,1,mmin,mmax);
      dmatrix V(1,1,mmin,mmax);
      BHess=get_H(mmin,mmax,p,eta,chinv,N,stdsinv,q,U,V);
      /*
      dmatrix BHess1=make_dmatrix(BHess);
      int m1=U.indexmax();
      for (int i=1;i<=m1;i++)
      {
        BHess1+=outer_prod(U(i),V(i));
      }
      */
      
      //dmatrix Hess=get_fd_H(mmin,mmax,p,eta,chinv,N,stdsinv,q);
      //cout << "VVV   " << norm2(BHess1-Hess) << endl;
      //cout << "WWW   " << max(BHess1-Hess) << endl;
      //cout << "WWW   " << BHess1-Hess << endl;
      //ad_exit(1);
      
      // ********************************************
      // rank one index update
      // ********************************************
      //dvector h3=solve(BHess,U,V,grad);
      banded_symmetric_dmatrix ddBH=banded_symmetric_dmatrix(BHess);
      dvector ddg=dvector(grad);
      dmatrix ddU=dmatrix(U);
      dmatrix ddV=dmatrix(V);
      dvector ddh3=solve(ddBH,ddU,ddV,ddg);
#if defined(USE_DD)
      dvector h3=make_dvector(ddh3);
#else
      dvector h3=ddh3;
#endif
      /*
      dmatrix BHess1=make_dmatrix(BHess);
      int m1=U.indexmax();
      for (int i=1;i<=m1;i++)
      {
        BHess1+=outer_prod(U(i),V(i));
      }
      dvector hh=solve(BHess1,grad);
      cout << norm(BHess*hh+U(1)*(V(1)*hh)-grad) << endl;
      cout << norm(BHess1*hh-grad) << endl;
      cout << norm(BHess*h3+U(1)*(V(1)*h3)-grad) << endl;
      MY_DOUBLE_TYPE ddiff=norm(BHess*h3+U(1)*(V(1)*h3)-grad);
      if (ddiff>1.e-10)
      {
        cout << "LLL  " << ddiff << endl;
        dmatrix BHess1=make_dmatrix(BHess);
        dvector hhu=make_dvector(ddh3);
        cout << "DDDD   " << norm(BHess*hhu+U(1)*(V(1)*hhu)-grad) << endl;
        cout << "DDDD2   " 
             << norm(ddBH*ddh3+ddU(1)*(ddV(1)*ddh3)-ddvector(grad)) << endl;
        int m1=U.indexmax();
        for (int i=1;i<=m1;i++)
        {
          BHess1+=outer_prod(U(i),V(i));
        }
        dmatrix ddBH1=dmatrix(BHess1);
        dvector ddg=dvector(grad);
        dvector hh=solve(BHess1,grad);
        dvector ddhh=solve(ddBH1,ddg);
        cout << norm(BHess*hh+U(1)*(V(1)*hh)-grad) << endl;
        cout << norm(BHess*hh+U(1)*(V(1)*hh)-grad) << endl;
        cout << norm(BHess1*hh-grad) << endl;
        cout << "DDDD   " << norm(ddBH1*ddhh-ddg) << endl;
     }
     */
      
      //dvector hk=solve(BHess1,grad);
      //dvector hh1=solve(Hess,grad);
      //cout << norm(h3-hh1) << endl;
      //cout << norm(hk-hh1) << endl;
      //cout << norm(hk-h3) << endl;
      eta-=h3;
      itcount++;
      oldnc=nc;
      mult=1.0;
    }
    else
    {
      mult*=0.5;
      eta+=mult*h3;
      //eta+=mult*vvv;
      itcount++;
    }
  }
  while(1);
  //cout << endl;
  /*
  int badflag=0;
  int negflag=0;
  eta=mynr(badflag,negflag,N,q,p,eta,chinv,stdsinv);
  while(badflag || negflag)
  {
    eta=mynr(badflag,negflag,N,q,p,eta,chinv,stdsinv);
  }
  //cout << "XX " << norm2(eta1-eta) << endl;
  */
  return eta;
}
dvector lognormal_ss_multinomial_newton_raphson(const MY_DOUBLE_TYPE N,
  const dvector& q,const dvector& p,
  //const dmatrix & cch,
  const banded_lower_triangular_dmatrix & chinv,
  const banded_symmetric_dmatrix& _stdsinv)
{
  ADUNCONST(banded_symmetric_dmatrix,stdsinv)
  int mmin=q.indexmin();
  int mmax=q.indexmax();
  // Newton Raphson
  
  dvector t(mmin,mmax);
  dvector Nq=N*q;
  dvector grad(mmin,mmax);
  int bw=stdsinv.bandwidth();
  banded_symmetric_dmatrix BHess(mmin,mmax,bw);
  dvector w(mmin,mmax);
  // get initial norm grad

  MY_DOUBLE_TYPE oldnc=1.e+100;
  int goodflag=1;
  int itcount=0;
  dvector h(mmin,mmax);
  dvector h3(mmin,mmax);
  dvector e(mmin,mmax);
  dvector ee(mmin,mmax);
  dvector eta(mmin,mmax);
  e=1.0;
  h.initialize();
  h3.initialize();
  MY_DOUBLE_TYPE mult=1.0;
  MY_DOUBLE_TYPE gmax=0.0;
  dvector eta_inner=inner_minimize(mmax,N,q,p,chinv,gmax,ee,stdsinv);
  if (fabs(gmax)>1.e-3)
  {
     eta_inner=inner_minimize(mmax,N,q,p,chinv,gmax,eta_inner,stdsinv);
  }
  if (fabs(gmax)>1.e-1)
  {
     eta_inner=inner_minimize(mmax,N,q,p,chinv,gmax,eta_inner,stdsinv);
  }
  eta=eta_inner;
 

  //cout << "norm2(eta)" <<  norm2(eta) << endl;
  //exit(1);
  do
  {
    MY_DOUBLE_TYPE v=0.0;
    for (int i=mmin;i<=mmax;i++)
    {
      // old lognormal version
      t(i)=phi(p(i),eta(i));
      v+=t(i);
    }
    //cout << eta << endl;
    //  f= -N*q*log(t)+N*log(v)
    for (int i=mmin;i<=mmax;i++)
    {
      grad(i)=(-N*q(i)/t(i)+N/v)*dphi_deta(p(i),eta(i));
    }
    dvector v2=stdsinv*eta;
    grad+=v2;
    
    MY_DOUBLE_TYPE nc=norm(grad);
    if (nc <1.e-12 || itcount > 18 ||
      nc <1.e-8 && itcount > 10)
    {
      if (itcount >15)
      {
        cout << "C convergence difficulty in newtraph " << nc << endl;
      }
      // get final hessian and quit
      BHess.initialize();
      dvector d1(mmin,mmax);
      dvector d2(mmin,mmax);
      for (int i=mmin;i<=mmax;i++)
      {
        df2_two_variable tmp=d2phi(p(i),eta(i));
        d1(i)=tmp.get_u_y();
        d2(i)=tmp.get_u_yy();

        //d1(i)=dphi_deta(p(i),eta(i));
        //d2(i)=d2phi_deta2(p(i),eta(i));
      }
      MY_DOUBLE_TYPE Nv=N/v;
      MY_DOUBLE_TYPE Nvv=Nv/v;
      dvector Nq(mmin,mmax);
      dvector Nqt(mmin,mmax);
      dvector Nqtt(mmin,mmax);
      for (int i=mmin;i<=mmax;i++)
      {
        Nq(i)=N*q(i);
      }
      for (int i=mmin;i<=mmax;i++)
      {
        Nqt(i)=Nq(i)/t(i);
      }
      for (int i=mmin;i<=mmax;i++)
      {
        Nqtt(i)=Nqt(i)/t(i);
      }
        
      for (int i=mmin;i<=mmax;i++)
      {
        BHess(i,i)+=(-Nqt(i)+Nv)*d2(i);
        BHess(i,i)+=Nqtt(i)*square(d1(i));
        w(i)=sqrt(Nvv)*d1(i);
      }
      for (int i=mmin;i<=mmax;i++)
      { 
        for (int j=i;j<=min(mmax,i+bw-1);j++)
        {
          BHess(j,i)+=stdsinv(j,i);
        }
      }
      /*
      for (int i=mmin;i<mmax;i++)
      {
        BHess(i,i)+=stdsinv(i,i);
        BHess(i+1,i)=stdsinv(i+1,i);
      }
      BHess(mmax,mmax)+=stdsinv(mmax,mmax);
      */
      //badglag=0.0;
      break;
    }
    if (nc<oldnc)
    {
      BHess.initialize();
      dvector d1(mmin,mmax);
      dvector d2(mmin,mmax);
      for (int i=mmin;i<=mmax;i++)
      {
        df2_two_variable tmp=d2phi(p(i),eta(i));
        d1(i)=tmp.get_u_y();
        d2(i)=tmp.get_u_yy();
        //d1(i)=dphi_deta(p(i),eta(i));
        //d2(i)=d2phi_deta2(p(i),eta(i));
      }
      MY_DOUBLE_TYPE Nv=N/v;
      MY_DOUBLE_TYPE Nvv=Nv/v;
      dvector Nq(mmin,mmax);
      dvector Nqt(mmin,mmax);
      dvector Nqtt(mmin,mmax);
      for (int i=mmin;i<=mmax;i++)
      {
        Nq(i)=N*q(i);
      }
      for (int i=mmin;i<=mmax;i++)
      {
        Nqt(i)=Nq(i)/t(i);
      }
      for (int i=mmin;i<=mmax;i++)
      {
        Nqtt(i)=Nqt(i)/t(i);
      }
      w.initialize();
        
      for (int i=mmin;i<=mmax;i++)
      {
       
        BHess(i,i)+=(-Nqt(i)+Nv)*d2(i);
        // this next line is part of the rank one update
        BHess(i,i)+=Nqtt(i)*square(d1(i));
        w(i)=sqrt(Nvv)*d1(i);
      }
      for (int i=mmin;i<mmax;i++)
      {
        BHess(i,i)+=stdsinv(i,i);
        BHess(i+1,i)=stdsinv(i+1,i);
      }
      BHess(mmax,mmax)+=stdsinv(mmax,mmax);

      // ********************************************
      // rank one index update
      // ********************************************
      dvector hk=solve(BHess,grad);
      dvector bx=solve(BHess,w);
      h3 = hk + (w*hk)*bx / ( 1 - w*bx);

      eta-=h3;
      itcount++;
      oldnc=nc;
      mult=1.0;
    }
    else
    {
      mult*=0.5;
      eta+=mult*h3;
      itcount++;
    }
  }
  while(1);
  //cout << endl;
  /*
  int badflag=0;
  int negflag=0;
  eta=mynr(badflag,negflag,N,q,p,eta,chinv,stdsinv);
  while(badflag || negflag)
  {
    eta=mynr(badflag,negflag,N,q,p,eta,chinv,stdsinv);
  }
  //cout << "XX " << norm2(eta1-eta) << endl;
  */
  return eta;
}


dvector lognormal_ss_multinomial_newton_raphson(const MY_DOUBLE_TYPE N,
  const dvector& q,const dvector& p,
  //const dmatrix & cch,
  const banded_lower_triangular_dmatrix & chinv,
  const symmetric_tridiagonal_dmatrix& _stdsinv)
{
  ADUNCONST(symmetric_tridiagonal_dmatrix,stdsinv)
  int mmin=q.indexmin();
  int mmax=q.indexmax();
  // Newton Raphson
  
  dvector t(mmin,mmax);
  dvector Nq=N*q;
  dvector grad(mmin,mmax);
  symmetric_tridiagonal_dmatrix BHess(mmin,mmax);
  dvector w(mmin,mmax);
  // get initial norm grad

  MY_DOUBLE_TYPE oldnc=1.e+100;
  int goodflag=1;
  int itcount=0;
  dvector h(mmin,mmax);
  dvector h3(mmin,mmax);
  dvector e(mmin,mmax);
  dvector ee(mmin,mmax);
  dvector eta(mmin,mmax);
  e=1.0;
  h.initialize();
  h3.initialize();
  MY_DOUBLE_TYPE mult=1.0;
  MY_DOUBLE_TYPE gmax=0.0;
  dvector eta_inner=inner_minimize(mmax,N,q,p,chinv,gmax,ee,stdsinv);
  if (fabs(gmax)>1.e-3)
  {
     eta_inner=inner_minimize(mmax,N,q,p,chinv,gmax,eta_inner,stdsinv);
  }
  if (fabs(gmax)>1.e-1)
  {
     eta_inner=inner_minimize(mmax,N,q,p,chinv,gmax,eta_inner,stdsinv);
  }
  eta=eta_inner;
 

  //cout << "norm2(eta)" <<  norm2(eta) << endl;
  //exit(1);
  do
  {
    MY_DOUBLE_TYPE v=0.0;
    for (int i=mmin;i<=mmax;i++)
    {
      // old lognormal version
      t(i)=phi(p(i),eta(i));
      v+=t(i);
    }
    //cout << eta << endl;
    //  f= -N*q*log(t)+N*log(v)
    for (int i=mmin;i<=mmax;i++)
    {
      grad(i)=(-N*q(i)/t(i)+N/v)*dphi_deta(p(i),eta(i));
    }
    dvector v2=stdsinv*eta;
    grad+=v2;
    
    MY_DOUBLE_TYPE nc=norm(grad);
    if (nc <1.e-12 || itcount > 18 ||
      nc <1.e-8 && itcount > 10)
    {
      if (itcount >15)
      {
        cout << "D convergence difficulty in newtraph " << nc << endl;
      }
      // get final hessian and quit
      BHess.initialize();
      dvector d1(mmin,mmax);
      dvector d2(mmin,mmax);
      for (int i=mmin;i<=mmax;i++)
      {
        df2_two_variable tmp=d2phi(p(i),eta(i));
        d1(i)=tmp.get_u_y();
        d2(i)=tmp.get_u_yy();

        //d1(i)=dphi_deta(p(i),eta(i));
        //d2(i)=d2phi_deta2(p(i),eta(i));
      }
      MY_DOUBLE_TYPE Nv=N/v;
      MY_DOUBLE_TYPE Nvv=Nv/v;
      dvector Nq(mmin,mmax);
      dvector Nqt(mmin,mmax);
      dvector Nqtt(mmin,mmax);
      for (int i=mmin;i<=mmax;i++)
      {
        Nq(i)=N*q(i);
      }
      for (int i=mmin;i<=mmax;i++)
      {
        Nqt(i)=Nq(i)/t(i);
      }
      for (int i=mmin;i<=mmax;i++)
      {
        Nqtt(i)=Nqt(i)/t(i);
      }
        
      for (int i=mmin;i<=mmax;i++)
      {
        BHess(i,i)+=(-Nqt(i)+Nv)*d2(i);
        BHess(i,i)+=Nqtt(i)*square(d1(i));
        w(i)=sqrt(Nvv)*d1(i);
      }
      for (int i=mmin;i<mmax;i++)
      {
        BHess(i,i)+=stdsinv(i,i);
        BHess(i+1,i)=stdsinv(i+1,i);
      }
      BHess(mmax,mmax)+=stdsinv(mmax,mmax);
      //badglag=0.0;
      break;
    }
    if (nc<oldnc)
    {
      BHess.initialize();
      dvector d1(mmin,mmax);
      dvector d2(mmin,mmax);
      for (int i=mmin;i<=mmax;i++)
      {
        df2_two_variable tmp=d2phi(p(i),eta(i));
        d1(i)=tmp.get_u_y();
        d2(i)=tmp.get_u_yy();
        //d1(i)=dphi_deta(p(i),eta(i));
        //d2(i)=d2phi_deta2(p(i),eta(i));
      }
      MY_DOUBLE_TYPE Nv=N/v;
      MY_DOUBLE_TYPE Nvv=Nv/v;
      dvector Nq(mmin,mmax);
      dvector Nqt(mmin,mmax);
      dvector Nqtt(mmin,mmax);
      for (int i=mmin;i<=mmax;i++)
      {
        Nq(i)=N*q(i);
      }
      for (int i=mmin;i<=mmax;i++)
      {
        Nqt(i)=Nq(i)/t(i);
      }
      for (int i=mmin;i<=mmax;i++)
      {
        Nqtt(i)=Nqt(i)/t(i);
      }
      w.initialize();
        
      for (int i=mmin;i<=mmax;i++)
      {
       
        BHess(i,i)+=(-Nqt(i)+Nv)*d2(i);
        // this next line is part of the rank one update
        BHess(i,i)+=Nqtt(i)*square(d1(i));
        w(i)=sqrt(Nvv)*d1(i);
      }
      for (int i=mmin;i<mmax;i++)
      {
        BHess(i,i)+=stdsinv(i,i);
        BHess(i+1,i)=stdsinv(i+1,i);
      }
      BHess(mmax,mmax)+=stdsinv(mmax,mmax);

      // ********************************************
      // rank one index update
      // ********************************************
      dvector hk=solve(BHess,grad);
      dvector bx=solve(BHess,w);
      h3 = hk + (w*hk)*bx / ( 1 - w*bx);

      eta-=h3;
      itcount++;
      oldnc=nc;
      mult=1.0;
    }
    else
    {
      mult*=0.5;
      eta+=mult*h3;
      itcount++;
    }
  }
  while(1);
  //cout << endl;
  /*
  int badflag=0;
  int negflag=0;
  eta=mynr(badflag,negflag,N,q,p,eta,chinv,stdsinv);
  while(badflag || negflag)
  {
    eta=mynr(badflag,negflag,N,q,p,eta,chinv,stdsinv);
  }
  //cout << "XX " << norm2(eta1-eta) << endl;
  */
  return eta;
}

dvector mynr(int & badflag,
  int & negflag,
  const MY_DOUBLE_TYPE N,
  const dvector& q,
  const dvector& p,
  const dvector& _eta,
  const banded_lower_triangular_dmatrix & chinv,
  const symmetric_tridiagonal_dmatrix& _stdsinv)
{ 
  ADUNCONST(symmetric_tridiagonal_dmatrix,stdsinv)
  ADUNCONST(dvector,eta)
  int mmin=q.indexmin();
  int mmax=q.indexmax();
  dvector grad(mmin,mmax);
  int itcount=0;
  MY_DOUBLE_TYPE mult=1.0;
  symmetric_tridiagonal_dmatrix BHess(mmin,mmax);
  dvector w(mmin,mmax);
  MY_DOUBLE_TYPE oldnc=1.e+100;
  dvector h3(mmin,mmax);
  dvector t(mmin,mmax);
  do
  {
    MY_DOUBLE_TYPE v=0.0;
    for (int i=mmin;i<=mmax;i++)
    {
      // old lognormal version
      t(i)=phi(p(i),eta(i));
      v+=t(i);
    }
    //cout << eta << endl;
    //  f= -N*q*log(t)+N*log(v)
    for (int i=mmin;i<=mmax;i++)
    {
      grad(i)=(-N*q(i)/t(i)+N/v)*dphi_deta(p(i),eta(i));
    }
    dvector v2=stdsinv*eta;
    grad+=v2;
    
    MY_DOUBLE_TYPE nc=norm(grad);
    if (nc <1.e-12 || itcount > 18 ||
      nc <1.e-8 && itcount > 10)
    {
      if (itcount >18)
      {
        cout << "E convergence difficulty in newtraph " << nc << endl;
        badflag++;
      }
      else 
      {
        badflag=0;
        negflag=0;
      }
      // get final hessian and quit
      BHess.initialize();
      dvector d1(mmin,mmax);
      dvector d2(mmin,mmax);
      for (int i=mmin;i<=mmax;i++)
      {
        d1(i)=dphi_deta(p(i),eta(i));
        d2(i)=d2phi_deta2(p(i),eta(i));
      }
      MY_DOUBLE_TYPE Nv=N/v;
      MY_DOUBLE_TYPE Nvv=Nv/v;
      dvector Nq(mmin,mmax);
      dvector Nqt(mmin,mmax);
      dvector Nqtt(mmin,mmax);
      for (int i=mmin;i<=mmax;i++)
      {
        Nq(i)=N*q(i);
      }
      for (int i=mmin;i<=mmax;i++)
      {
        Nqt(i)=Nq(i)/t(i);
      }
      for (int i=mmin;i<=mmax;i++)
      {
        Nqtt(i)=Nqt(i)/t(i);
      }
        
      for (int i=mmin;i<=mmax;i++)
      {
        BHess(i,i)+=(-Nqt(i)+Nv)*d2(i);
        BHess(i,i)+=Nqtt(i)*square(d1(i));
        w(i)=sqrt(Nvv)*d1(i);
      }
      for (int i=mmin;i<mmax;i++)
      {
        BHess(i,i)+=stdsinv(i,i);
        BHess(i+1,i)=stdsinv(i+1,i);
      }
      BHess(mmax,mmax)+=stdsinv(mmax,mmax);
      break;
    }
    if (nc<oldnc)
    {
      BHess.initialize();
      dvector d1(mmin,mmax);
      dvector d2(mmin,mmax);
      for (int i=mmin;i<=mmax;i++)
      {
        d1(i)=dphi_deta(p(i),eta(i));
        d2(i)=d2phi_deta2(p(i),eta(i));
      }
      MY_DOUBLE_TYPE Nv=N/v;
      MY_DOUBLE_TYPE Nvv=Nv/v;
      dvector Nq(mmin,mmax);
      dvector Nqt(mmin,mmax);
      dvector Nqtt(mmin,mmax);
      for (int i=mmin;i<=mmax;i++)
      {
        Nq(i)=N*q(i);
      }
      for (int i=mmin;i<=mmax;i++)
      {
        Nqt(i)=Nq(i)/t(i);
      }
      for (int i=mmin;i<=mmax;i++)
      {
        Nqtt(i)=Nqt(i)/t(i);
      }
      w.initialize();
        
      for (int i=mmin;i<=mmax;i++)
      {
       
        BHess(i,i)+=(-Nqt(i)+Nv)*d2(i);
        // this next line is part of the rank one update
        BHess(i,i)+=Nqtt(i)*square(d1(i));
        w(i)=sqrt(Nvv)*d1(i);
      }
      for (int i=mmin;i<mmax;i++)
      {
        BHess(i,i)+=stdsinv(i,i);
        BHess(i+1,i)=stdsinv(i+1,i);
      }
      BHess(mmax,mmax)+=stdsinv(mmax,mmax);
      if (badflag)
      {
        for (int i=mmin;i<=mmax;i++)
        {
          BHess(i,i)+=0.5*badflag;
        }
      }
      if (negflag)
      {
        for (int i=mmin;i<=mmax;i++)
        {
          BHess(i,i)+=0.1*negflag;
        }
      }
        
      // ********************************************
      // rank one index update
      // ********************************************
      dvector hk=solve(BHess,grad);
      dvector bx=solve(BHess,w);
      h3 = hk + (w*hk)*bx / ( 1 - w*bx);
      if (grad*h3 < 0)
      {
        dmatrix tmpmat(mmin,mmax,mmin,mmax);
        tmpmat.initialize();
        for (int i=mmin;i<mmax;i++)
        {
          tmpmat(i,i)=BHess(i,i);
          tmpmat(i+1,i)=BHess(i+1,i);
          tmpmat(i,i+1)=tmpmat(i+1,i);
        }
        tmpmat(mmax,mmax)=BHess(mmax,mmax);
        for (int i=mmin;i<=mmax;i++)
        {
          for (int j=mmin;j<=mmax;j++)
          {
            tmpmat(i,j)-=w(i)*w(j);
          }
        }
        if (badflag)
        {
          for (int i=mmin;i<=mmax;i++)
          {
            BHess(i,i)+=0.5*badflag;
          }
        }
        cerr << "infeasible direction" 
             << grad*hk/(1.e-50+norm(grad)*(norm(hk))) << "  "
             << grad*h3/(1.e-50+norm(grad)*norm(h3))<< endl;
        cerr << 1-w*bx << endl;
        cerr << norm2(solve(tmpmat,grad) - h3) << endl;
        negflag++;
        break;
      }

      eta-=h3;
      itcount++;
      oldnc=nc;
      mult=1.0;
    }
    else
    {
      mult*=0.5;
      eta+=mult*h3;
      itcount++;
    }
  }
  while(1);
  //badflag=0;
 
  return eta;
}
