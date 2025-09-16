/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include <admodel.h>

dvariable normal_length_to_weight(MY_DOUBLE_TYPE delta,MY_DOUBLE_TYPE xmin,
  const dvariable& _mean,const dvariable & _sigma,const dvariable& _a,
  const dvariable & _b)
{
  MY_DOUBLE_TYPE eps=.01;
  MY_DOUBLE_TYPE mean=value(_mean);
  MY_DOUBLE_TYPE sigma=value(_sigma);
  MY_DOUBLE_TYPE b=value(_b);
  MY_DOUBLE_TYPE a=value(_a);
  int n= int((-2*xmin)/delta)+1;
  dvector x(1,n);
  dvector w(1,n);
  dvector y(1,n);
  dvector z(1,n);
  x[1]=xmin-delta;
  MY_DOUBLE_TYPE f=0;
  int i;
  for (i=1;i<=n;i++)
  {
    w[i]=exp(-0.5*x(i)*x(i));
    y[i]=sigma*x(i)+mean;
    if (y[i]<eps) {
      z[i]=eps/(2-y[i]/eps);
      f+=pow(z[i],b)*w[i];
    }
    else
    {
      f+=pow(y[i],b)*w[i];
    }
    if (i<n) x[i+1]=x[i]+delta;
  }
  MY_DOUBLE_TYPE f1=a*f*0.39894*delta;
  dvariable tmp;
  value(tmp)=f1;
  dvector dfy(1,n);
  dvector dfz(1,n);
  dfy.initialize();
  dfz.initialize();
  MY_DOUBLE_TYPE dfa,dfb,dfsigma,dfmean,df1,dff;
  dfsigma=0.0;
  dfmean=0.0;
  // adjoint code
  df1=1.0;
  dfb=0.0;   // DF corrected 19NOV03
  //double f1=f*0.39894*delta;
  dff=df1*a*0.39894*delta;
  dfa=df1*f*0.39894*delta;
  for (i=n;i>=1;i--)
  {
    if (y[i]<eps) {
      //f+=pow(z[i],b)*w[i];
#if !defined(NO_MY_DOUBLE_TYPE)
      dfz[i]+=dff*b*pow(z[i],b-1.0L)*w[i];
#else
      dfz[i]+=dff*b*pow(z[i],b-1.0)*w[i];
#endif
      dfb+=dff*pow(z[i],b)*log(z[i])*w[i];
      //z[i]=eps/(2-y[i]/eps);
      dfy[i]+=dfz[i]/square(2-y[i]/eps);
    }
    else
    {
      //f+=pow(y[i],b)*w[i];
#if !defined(NO_MY_DOUBLE_TYPE)
      dfy[i]+=dff*b*pow(y[i],b-1.0L)*w[i];
#else
      dfy[i]+=dff*b*pow(y[i],b-1.0)*w[i];
#endif
      dfb+=dff*pow(y[i],b)*log(y[i])*w[i];
    }
    //y[i]=sigma*x(i)+mean;
    dfsigma+=dfy[i]*x(i);
    dfmean+=dfy[i];
  }
  gradient_structure::GRAD_STACK1->set_gradient_stack(default_evaluation4ind,
    address(tmp),
    address(_mean),dfmean ,
    address(_sigma),dfsigma ,
    address(_a),dfa,
    address(_b),dfb);

  return tmp;
}
