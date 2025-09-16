/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include <admodel.h>
#include "phi_stuff.h"
#include "tridiagonal_dmatrix.h"

dvector inner_minimize(int n,MY_DOUBLE_TYPE N,const dvector & q,const dvector& p,
  const symmetric_tridiagonal_dmatrix & Sinv);


dvector lognormal_ss_multinomial_newton_raphson(const MY_DOUBLE_TYPE N,
  const dvector& q,const dvector& p,const MY_DOUBLE_TYPE rho)
{
  int mmin=q.indexmin();
  int mmax=q.indexmax();
  symmetric_tridiagonal_dmatrix Sinv(mmin,mmax);
  dvector w(mmin,mmax);
  // Newton Raphson
  
  dvector t(mmin,mmax);
  dvector Nq=N*q;
  dvector grad(mmin,mmax);
  dmatrix Hess(mmin,mmax,mmin,mmax);
  // get initial norm grad

  MY_DOUBLE_TYPE oldnc=1.e+100;
  int goodflag=1;
  int itcount=0;
  dvector h(mmin,mmax);
  dvector e(mmin,mmax);
  e=1.0;
  h.initialize();
  MY_DOUBLE_TYPE mult=1.0;
  dvector eta=inner_minimize(mmax,N,q,p,Sinv);
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
    grad+=Sinv*eta;
    MY_DOUBLE_TYPE nc=norm(grad);
    //cout << "iteration " << itcount << " norm grad " << nc << endl;
    if (itcount==12)
    {
     // cout << "HERE" << endl;
    }
    if (nc <1.e-12 || itcount > 15 ||
      nc <1.e-6 && itcount > 10)
    {
      // get final hessian and quit
      Hess.initialize();
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
        Hess(i,i)+=(-Nqt(i)+Nv)*d2(i);
        Hess(i,i)+=(Nqtt(i)-Nvv)*square(d1(i));
        for (int j=mmin;j<i;j++)
        {
          Hess(i,j)-=Nvv*d1(i)*d1(j);
          Hess(j,i)= Hess(i,j);
        }
      }
      for (int i=mmin;i<=mmax;i++)
      {
        Sinv(i,i)+=(-Nqt(i)+Nv)*d2(i);
        Sinv(i,i)+=(Nqtt(i)-Nvv)*square(d1(i));
        w(i)=sqrt(Nvv)*d1(i);
      }
      Hess+=make_dmatrix(Sinv);
      choleski_decomp(Hess);
      break;
    }
    if (nc<oldnc)
    {
      Hess.initialize();
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
        Hess(i,i)+=(-Nqt(i)+Nv)*d2(i);
        Hess(i,i)+=(Nqtt(i)-Nvv)*square(d1(i));
        for (int j=mmin;j<i;j++)
        {
          Hess(i,j)-=Nvv*d1(i)*d1(j);
          Hess(j,i)= Hess(i,j);
        }
      }
      for (int i=mmin;i<=mmax;i++)
      {
        Sinv(i,i)+=(-Nqt(i)+Nv)*d2(i);
        Sinv(i,i)+=(Nqtt(i)-Nvv)*square(d1(i));
        w(i)=sqrt(Nvv)*d1(i);
      }
      Hess+=make_dmatrix(Sinv);
      h=solve(Hess,grad);
      eta-=h;
      itcount++;
      oldnc=nc;
      mult=1.0;
    }
    else
    {
      mult*=0.5;
      eta+=mult*h;
      itcount++;
    }
  }
  while(1);
  //cout << endl;
  return eta;
}

dvector inner_minimize(int n,MY_DOUBLE_TYPE N,const dvector & q,const dvector& p,
  const symmetric_tridiagonal_dmatrix & Sinv)
{
  dmatrix cch;
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
      eta=cch*eps;
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
      g=g*cch;
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
