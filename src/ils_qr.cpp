/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
//#include <admodel.h>
#include "ils_qr.hpp"

/*
dvar_vector max(double m,const dvar_vector& _v)
{
  const double delta=0.0005;
  const double a2=3.0/(delta*delta);
  const double a3=-2.0/(delta*delta*delta);
  ADUNCONST(dvar_vector,v)
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  dvar_vector tmp(mmin,mmax);
  double bound=m+delta;
  for (int i=mmin;i<=mmax;i++)
  {
    if (v(i)<=bound)
    {
      if (v(i)<=m)
      {
        tmp(i)=m;
      }
      else
      {
        dvariable x=bound-v(i);
        dvariable p=square(x)*(3.0-2.0*x);
        tmp(i)=(1.0-p)*m + p*v(i);
      }
    }
  }
} 
*/
    
dvar_vector max(double m,const dvar_vector& _v)
{
  const double delta=0.001;
  const double a2=4.0;
  const double a3=-1.0;
  ADUNCONST(dvar_vector,v)
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  dvar_vector tmp(mmin,mmax);
  tmp=v;
  double bound=m+delta;
  for (int i=mmin;i<=mmax;i++)
  {
    if (v(i)<=bound)
    {
      if (v(i)<=m)
      {
        tmp(i)=m;
      }
      else
      {
        dvariable x=(v(i)-m)/delta;
        dvariable p=square(x)*(3.0-2.0*x);
        tmp(i)=(1.0-p)*m + p*v(i);
        //cout << x << " ";
      }
    }
  }
  return tmp;
} 
  

ils_manager::ils_manager(dmatrix& _M,dvector& _obs,int _niter) : 
  obs(_obs), niter(_niter) 
{ 
  if (_M.indexmin() !=1 || _obs.indexmin() !=1)
  {  
    cerr << "shape error in ils_manager  minindex must equal 1" << endl;
    ad_exit(1);
  }
  if (_M.indexmax() != _obs.indexmax())
  {  
    cerr << "shape error in ils_manager  M and obs must be the" 
            " same shape" << endl;
    ad_exit(1);
  }
  nobs=_M.indexmax();
  obs1.allocate(1,nobs);
  int mmin=_M.indexmin();
  int mmax=_M.indexmax();
  M.allocate(mmin,mmax,_M(mmin).indexmin(),_M(mmin).indexmax());
  M=_M;
  M1.allocate(1,niter,mmin,mmax,_M(mmin).indexmin(),
    _M(mmin).indexmax());
}

double ils_manager::fit_data(void)
{
  cerr << "this is not properly implemented" << endl;
  ad_exit(1);
  double r2=0.0;
  gram_schmidt_qr(M,Q,R);
  theta_hat= solve(trans(R)*R,trans(M)*obs);
  pred_obs=M*theta_hat;
  vhat=norm2(obs-pred_obs)/nobs;
  weights=1.0/(0.8*(0.2+fabs(obs-pred_obs)/sqrt(vhat)));
  //dmatrix M1(1,nobs,1,npar);

  int it=1;
  M1(it)=M;
  for (int i=1;i<=nobs;i++)
  {
    M1(it,i)*=sqrt(weights(i));
  }
  obs1=elem_prod(weights,obs);
  gram_schmidt_qr(M1(it),Q,R);
  theta_hat= solve(trans(R)*R,trans(M1(it))*obs1);
  pred_obs=M1(it)*theta_hat;
  r2=norm2(obs1-pred_obs);
  vhat=norm2(obs1-pred_obs)/nobs;
  //weights=1.0/(0.8*(0.2+fabs(obs1-pred_obs)/sqrt(vhat)));

  /*
  it=2;
  M1(it)=M1(it-1);
  for (int i=1;i<=nobs;i++)
  {
    M1(it,i)*=sqrt(weights(i));
  }
  obs1=elem_prod(weights,obs1);
  gram_schmidt_qr(M1(it),Q,R);
  theta_hat= solve(trans(R)*R,trans(M1(it))*obs1);
  pred_obs=M1(it)*theta_hat;
  vhat=norm2(obs1-pred_obs)/nobs;
  weights=1.0/(0.8*(0.2+fabs(obs1-pred_obs)/sqrt(vhat)));
  for (int it=2;it<=niter;it++)
  {
    M1(it)=M1(it-1);
    for (int i=1;i<=nobs;i++)
    {
      M1(it,i)*=sqrt(weights(i));
    }
    obs1=elem_prod(weights,obs1);
    gram_schmidt_qr(M1(it),Q,R);
    theta_hat= solve(trans(R)*R,trans(M1(it))*obs1);
    pred_obs=M1(it)*theta_hat;
    r2=norm2(obs1-pred_obs);
    vhat=r2/nobs;
    weights=1.0/(0.8*(0.2+fabs(obs1-pred_obs)/sqrt(vhat)));
  }
  */
  return r2;
}

double ils(dmatrix& M,dvector& obs,dvector& weights,int _iter)
{
  int nobs=obs.indexmax()-obs.indexmin()+1;
  dmatrix Q;
  dmatrix R;
  gram_schmidt_qr(M,Q,R);
  dvector  theta_hat= solve(trans(R)*R,trans(M)*obs);
  dvector pred_obs=M*theta_hat;
  double vhat=norm2(obs-pred_obs)/nobs;
  weights=1.0/(0.8*(0.2+fabs(obs-pred_obs)/sqrt(vhat)));
  return 1.0;
}


/*
//this is already in newrshimp_experiment.cpp
void gram_schmidt_qr(dmatrix& M,dmatrix& Q,dmatrix& R)
{
  int m=M.indexmax();
  int n=M(m).indexmax();
  dmatrix TQ(1,n,1,m);
  TQ=trans(M);

  dmatrix ID(1,n,1,n);
  ID.initialize();
  for (int i=1;i<=n;i++)
  {
    ID(i,i)=1.0;
  }

  dmatrix TR(1,n,1,n);
  TR=ID;
  for (int i=1;i<=n;i++)
  {
    double a=norm(TQ(i));
    TQ(i)/=a;
    TR(i)/=a;
    for (int j=i+1;j<=n;j++)
    {
      double a=TQ(i)*TQ(j);
      TQ(j)-=a*TQ(i);
      TR(j)-=a*TR(i);
    }
  }

  Q=trans(TQ);
  R=inv(trans(TR));

}
*/

void gram_schmidt_qr(dvar_matrix& M,dvar_matrix& Q,dvar_matrix& R)
{
  int m=M.indexmax();
  int n=M(m).indexmax();
  dvar_matrix TQ(1,n,1,m);
  TQ=trans(M);

  dmatrix ID(1,n,1,n);
  ID.initialize();
  for (int i=1;i<=n;i++)
  {
    ID(i,i)=1.0;
  }

  dvar_matrix TR(1,n,1,n);
  TR=ID;
  for (int i=1;i<=n;i++)
  {
    dvariable a=norm(TQ(i));
    TQ(i)/=a;
    TR(i)/=a;
    for (int j=i+1;j<=n;j++)
    {
      dvariable a=TQ(i)*TQ(j);
      TQ(j)-=a*TQ(i);
      TR(j)-=a*TR(i);
    }
  }

  Q=trans(TQ);
  R=inv(trans(TR));

}

dfils_manager::dfils_manager(dmatrix& _M,dvar_vector& _obs,int _niter) : 
  obs(_obs), niter(_niter) , theta_hat(0,_niter)  
{ 
  if (_M.indexmin() !=1 || _obs.indexmin() !=1)
  {  
    cerr << "shape error in ils_manager  minindex must equal 1" << endl;
    ad_exit(1);
  }
  if (_M.indexmax() != _obs.indexmax())
  {  
    cerr << "shape error in ils_manager  M and obs must be the" 
            " same shape" << endl;
    ad_exit(1);
  }
  nobs=_M.indexmax();
  obs1.allocate(1,nobs);
  //minval.allocate(1,nobs);
  //minval=0.001;
  int mmin=_M.indexmin();
  int mmax=_M.indexmax();
  int cmin=_M(M.indexmin()).indexmin();
  int cmax=_M(M.indexmin()).indexmax();
  M.allocate(mmin,mmax,_M(mmin).indexmin(),_M(mmin).indexmax());
  M=_M;
  M1.allocate(1,niter,mmin,mmax,_M(mmin).indexmin(),
    _M(mmin).indexmax());
  theta_hat.initialize();
  silly_pred.allocate(0,_niter,mmin,mmax);
}
void dfils_manager::fit_data(void)
{
  dvariable r2=0.0;
  gram_schmidt_qr(M,Q,R);
  theta_hat(0)= solve(trans(R)*R,trans(M)*obs);
  pred_obs=M*theta_hat(0);
  //silly_pred(0)=value(pred_obs);
  //r2=norm2(obs-pred_obs);
  weights=sqrt(1.0/max(0.001,fabs(obs-pred_obs)));

  int it=1;
  M1(it)=M;
  for (int i=1;i<=nobs;i++)
  {
    M1(it,i)*=weights(i);
  }
  obs1=elem_prod(weights,obs);
  gram_schmidt_qr(M1(it),Q1,R1);
  theta_hat(it)= solve(trans(R1)*R1,trans(M1(it))*obs1);
  pred_obs=M1(it)*theta_hat(it);
  //silly_pred(it)=value(pred_obs);
  //weights=sqrt(1.0/max(0.001,fabs(obs-pred_obs)));

  /*
  it=2;
  M1(it)=M;
  for (int i=1;i<=nobs;i++)
  {
    M1(it,i)*=weights(i);
  }
  obs1=elem_prod(weights,obs);
  gram_schmidt_qr(M1(it),Q1,R1);
  theta_hat(it)= solve(trans(R1)*R1,trans(M1(it))*obs1);
  pred_obs=M1(it)*theta_hat(it);
  silly_pred(it)=value(pred_obs);
  weights=sqrt(1.0/max(0.001,fabs(obs-pred_obs)));
  for (int it=2;it<=niter;it++)
  {
    M1(it)=M;
    for (int i=1;i<=nobs;i++)
    {
      M1(it,i)*=weights(i);
    }
    obs1=elem_prod(weights,obs);
    gram_schmidt_qr(M1(it),Q1,R1);
    theta_hat(it)= solve(trans(R1)*R1,trans(M1(it))*obs1);
    pred_obs=M1(it)*theta_hat(it);
    silly_pred(it)=value(pred_obs);
    r2=norm2(obs1-pred_obs);
    weights=sqrt(1.0/max(0.001,fabs(obs-pred_obs)));
  }
  return r2;
  */
}
