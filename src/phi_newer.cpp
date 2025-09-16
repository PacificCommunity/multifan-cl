/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include <admodel.h>
#include "df22fun.h"
#include "df32fun.h"


static const MY_DOUBLE_TYPE  eps_p=0.50;
static const MY_DOUBLE_TYPE  p_exp=0.01;
static const  int xxplain_flag=1;

MY_DOUBLE_TYPE psiprime(MY_DOUBLE_TYPE p)
{
  cerr << "cant happen" << endl;
  ad_exit(1);
  MY_DOUBLE_TYPE tmp=1/p + 1/(1-p);
  return tmp;
}
MY_DOUBLE_TYPE psiprime(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE e)
{
  cerr << "cant happen" << endl;
  ad_exit(1);
  MY_DOUBLE_TYPE tmp=1/p + 1/(1-p);
  return tmp;
  //return psiprime(p)-e/(square(eps_p+sqrt(p))*(2.0*sqrt(p)));
}
 
MY_DOUBLE_TYPE psi(MY_DOUBLE_TYPE p)
{
  cerr << "cant happen" << endl;
  ad_exit(1);
  MY_DOUBLE_TYPE tmp=log(p) - log(1-p);
  return tmp;
}
 
MY_DOUBLE_TYPE d2phi_deta2_1(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE eta)
{
  cerr << "cant happen" << endl;
  ad_exit(1);
  //double z = eta + psi(p);
  MY_DOUBLE_TYPE q=eps_p+sqrt(p);
  MY_DOUBLE_TYPE z = eta/q + psi(p);
  MY_DOUBLE_TYPE u=mfexp(-z);
  MY_DOUBLE_TYPE v=1/(1+u);
  MY_DOUBLE_TYPE r=2.0*u*v-1.0;
  MY_DOUBLE_TYPE s=u*square(v);
  MY_DOUBLE_TYPE f=s*r/square(q);

  MY_DOUBLE_TYPE dfu=0.0;
  MY_DOUBLE_TYPE dfq=0.0;
  MY_DOUBLE_TYPE dfv=0.0;
  MY_DOUBLE_TYPE dfz=0.0;
  MY_DOUBLE_TYPE dfs=0.0;
  MY_DOUBLE_TYPE dfr=0.0;
  MY_DOUBLE_TYPE dfp=0.0;
  MY_DOUBLE_TYPE dfeta=0.0;

  //double f=u*square(v);

  MY_DOUBLE_TYPE dff=1.0;

  //double f=s*r/square(q);
  dfs+=dff*r/square(q);
  dfr+=dff*s/square(q);
  dfq-=2.0*dff*s*r/cube(q);

  //double s=u*square(v);
  dfu+=dfs*square(v);
  dfv+=dfs*2*u*v;

  //double r=2.0*u*v-1.0;
  dfu+=dfr*2.0*v;
  dfv+=dfr*2.0*u;

  //double v=1/(1+u);
  dfu-=dfv*square(v);

  //double u=exp(-z);
  dfz-=dfu*u;

  //double z = eta + psi(p);
  //dfeta+=dfz;
  //dfp+=dfz*psiprime(p);

  //double z = eta/q + psi(p);
  dfeta+=dfz/q; 
  dfp+=dfz*psiprime(p);
  dfq-=dfz*eta/square(q);
  dfz=0.0;

  //double q=eps_p+sqrt(p);
  dfp+=dfq/(2.0*sqrt(p));
  dfq=0.0;
  
  return dfp;
}
  
MY_DOUBLE_TYPE d2phi_deta2_2(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE eta)
{
  cerr << "cant happen" << endl;
  ad_exit(1);
  MY_DOUBLE_TYPE q=eps_p+sqrt(p);
  MY_DOUBLE_TYPE z = eta/q + psi(p);
  //double z = eta + psi(p);
  MY_DOUBLE_TYPE u=mfexp(-z);
  MY_DOUBLE_TYPE v=1/(1+u);
  MY_DOUBLE_TYPE r=2.0*u*v-1.0;
  MY_DOUBLE_TYPE s=u*square(v);
  MY_DOUBLE_TYPE f=s*r/square(q);

  MY_DOUBLE_TYPE dfu=0.0;
  MY_DOUBLE_TYPE dfq=0.0;
  MY_DOUBLE_TYPE dfv=0.0;
  MY_DOUBLE_TYPE dfz=0.0;
  MY_DOUBLE_TYPE dfs=0.0;
  MY_DOUBLE_TYPE dfr=0.0;
  MY_DOUBLE_TYPE dfp=0.0;
  MY_DOUBLE_TYPE dfeta=0.0;

  MY_DOUBLE_TYPE dff=1.0;

  //double f=s*r/square(q);
  dfs+=dff*r/square(q);
  dfr+=dff*s/square(q);
  dfq-=2.0*dff*s*r/cube(q);
  dff=0.0;

  //double s=u*square(v);
  dfu+=dfs*square(v);
  dfv+=dfs*2*u*v;
  dfs=0.0;

  //double r=2.0*u*v-1.0;
  dfu+=dfr*2.0*v;
  dfv+=dfr*2.0*u;
  dfr=0.0;

  //double v=1/(1+u);
  dfu-=dfv*square(v);
  dfv=0.0;

  //double u=exp(-z);
  dfz-=dfu*u;
  dfu=0.0;

  //double z = eta + psiprime(p);
  //dfeta+=dfz;
  //dfp+=dfz*psiprime(p);

  //double z = eta/q + psi(p);
  dfeta+=dfz/q; 
  dfp+=dfz*psiprime(p);
  dfq-=dfz*eta/square(q);
  dfz=0.0;

  //double q=eps_p+sqrt(p);
  dfp+=dfq/(2.0*sqrt(p));
  dfq=0.0;

  return dfeta;
}
  
MY_DOUBLE_TYPE d2phi_deta2(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE eta)
{
  cerr << "cant happen" << endl;
  ad_exit(1);
  //double z = eta + psi(p);
  MY_DOUBLE_TYPE q=eps_p+sqrt(p);
  MY_DOUBLE_TYPE z = eta/q + psi(p);
  MY_DOUBLE_TYPE u=mfexp(-z);
  MY_DOUBLE_TYPE v=1/(1+u);
  MY_DOUBLE_TYPE r=2.0*u*v-1.0;
  MY_DOUBLE_TYPE s=u*square(v);
  MY_DOUBLE_TYPE f=s*r/square(q);
  return f;
}
  
MY_DOUBLE_TYPE dphi_deta_2(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE eta)
{
  cerr << "cant happen" << endl;
  ad_exit(1);
  //double z = eta + psi(p);
  MY_DOUBLE_TYPE q=eps_p+sqrt(p);
  MY_DOUBLE_TYPE z = eta/q + psi(p);
  MY_DOUBLE_TYPE u=mfexp(-z);
  MY_DOUBLE_TYPE v=1/(1+u);
  MY_DOUBLE_TYPE f=u*square(v)/q;

  MY_DOUBLE_TYPE dff=1.0;
  MY_DOUBLE_TYPE dfq=0.0;
  MY_DOUBLE_TYPE dfu=0.0;
  MY_DOUBLE_TYPE dfv=0.0;
  MY_DOUBLE_TYPE dfz=0.0;
  MY_DOUBLE_TYPE dfp=0.0;
  MY_DOUBLE_TYPE dfeta=0.0;

  //double f=u*square(v)/q;
  dfu+=dff*square(v)/q;
  dfv+=dff*2*u*v/q;
  dfq-=dff*u*square(v)/square(q);
  

  //double v=1/(1+u);
  dfu-=dfv*square(v);

  //double u=exp(-z);
  dfz-=dfu*u;

  //double z = eta + psiprime(p);
  //dfeta+=dfz;
  //dfp+=dfz*psiprime(p);
 
  //double z = eta/q + psi(p);
  dfeta+=dfz/q; 
  dfp+=dfz*psiprime(p);
  dfq-=dfz*eta/square(q);
  dfz=0.0;

  //double q=eps_p+sqrt(p);
  dfp+=dfq/(2.0*sqrt(p));
  dfq=0.0;

  return dfeta;
}
MY_DOUBLE_TYPE dphi_deta_1(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE eta)
{
  cerr << "cant happen" << endl;
  ad_exit(1);
  //double z = eta + psi(p);
  MY_DOUBLE_TYPE q=eps_p+sqrt(p);
  MY_DOUBLE_TYPE z = eta/q + psi(p);
  MY_DOUBLE_TYPE u=mfexp(-z);
  MY_DOUBLE_TYPE v=1/(1+u);
  MY_DOUBLE_TYPE f=u*square(v)/q;

  MY_DOUBLE_TYPE dff=1.0;
  MY_DOUBLE_TYPE dfq=0.0;
  MY_DOUBLE_TYPE dfu=0.0;
  MY_DOUBLE_TYPE dfv=0.0;
  MY_DOUBLE_TYPE dfz=0.0;
  MY_DOUBLE_TYPE dfp=0.0;
  MY_DOUBLE_TYPE dfeta=0.0;

  //double f=u*square(v)/q;
  dfq-=dff*u*square(v)/square(q);
  dfu+=dff*square(v)/q;
  dfv+=dff*2*u*v/q;

  //double v=1/(1+u);
  dfu-=dfv*square(v);

  //double u=exp(-z);
  dfz-=dfu*u;

  //double z = eta + psiprime(p);
  //dfeta+=dfz;
  //dfp+=dfz*psiprime(p);
  //double z = eta/q + psi(p);
  dfeta+=dfz/q; 
  dfp+=dfz*psiprime(p);
  dfq-=dfz*eta/square(q);
  dfz=0.0;

  //double q=eps_p+sqrt(p);
  dfp+=dfq/(2.0*sqrt(p));
  dfq=0.0;
  return dfp;

}
  
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


df2_two_variable d2phi(MY_DOUBLE_TYPE _p,MY_DOUBLE_TYPE _eta)
{
  df2_two_variable p;
  df2_two_variable eta;
  df2_two_variable pexp;
  p=_p;
  eta=_eta;
  p.set_independent_1();
  eta.set_independent_2();
  pexp=p_exp;
  if (xxplain_flag)
  {
    df2_two_variable f=p*exp(eta);
    return f;
  }
  else
  {
    df2_two_variable q=eps_p+pow(p,pexp);
    df2_two_variable eq=eta/q;
    df2_two_variable f=p*exp(eq);
    return f;
  }
}

df3_two_variable d3phi(MY_DOUBLE_TYPE _p,MY_DOUBLE_TYPE _eta)
{
  df3_two_variable p;
  df3_two_variable eta;
  df3_two_variable pexp;
  p=_p;
  eta=_eta;
  p.set_independent_1();
  eta.set_independent_2();
  pexp=p_exp;
  if (xxplain_flag)
  {
    df3_two_variable f=p*exp(eta);
    return f;
  }
  else
  {
    df3_two_variable q=eps_p+pow(p,pexp);
    df3_two_variable eq=eta/q;
    df3_two_variable f=p*exp(eq);
    return f;
  }
}



MY_DOUBLE_TYPE phi(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE eta)
{
  if (xxplain_flag)
  {
    MY_DOUBLE_TYPE f=p*exp(eta);
    return f;
  }
  else
  {
    MY_DOUBLE_TYPE q=eps_p+pow(p,p_exp);
    MY_DOUBLE_TYPE eq=eta/q;
    MY_DOUBLE_TYPE f=p*exp(eq);
    return f;
  }
}


MY_DOUBLE_TYPE dphi_deta(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE eta)
{
 /* old version
  if (xxplain_flag)
  {
    //double f=p*exp(eta);
    MY_DOUBLE_TYPE df=1.0;
    MY_DOUBLE_TYPE dfeta=p*exp(eta);
    return dfeta;
  }
  else
  {
    MY_DOUBLE_TYPE q=eps_p+pow(p,p_exp);
    MY_DOUBLE_TYPE eq=eta/q;
    MY_DOUBLE_TYPE eeq=exp(eq);
    MY_DOUBLE_TYPE f=p*eeq;
    MY_DOUBLE_TYPE df=1.0;
    //double f=p*eeq;
    MY_DOUBLE_TYPE dfp=df*eeq;
    MY_DOUBLE_TYPE dfeeq=df*p;
    df=0;
    //double eeq=exp(eq);
    MY_DOUBLE_TYPE dfeq=dfeeq*eeq;
    dfeeq=0.0;
    //double eq=eta/q;
    MY_DOUBLE_TYPE dfeta=dfeq*1.0/q;
    dfeq=0.0;
    return dfeta;
  }
  */
  df2_two_variable tmp=d2phi(p,eta);
  return tmp.get_u_y();
}

MY_DOUBLE_TYPE dphi_dp(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE eta)
{
  if (xxplain_flag)
  {
    //double f=p*exp(eta);
    MY_DOUBLE_TYPE df=1.0;
    MY_DOUBLE_TYPE dfp=exp(eta);
    return dfp;
  }
  else
  {
    MY_DOUBLE_TYPE q=eps_p+pow(p,p_exp);
    MY_DOUBLE_TYPE eq=eta/q;
    MY_DOUBLE_TYPE eeq=exp(eq);
    MY_DOUBLE_TYPE f=p*eeq;
    MY_DOUBLE_TYPE df=1.0;
    //double f=p*eeq;
    MY_DOUBLE_TYPE dfp=df*eeq;
    MY_DOUBLE_TYPE dfeeq=df*p;
    df=0;
    //double eeq=exp(eq);
    MY_DOUBLE_TYPE dfeq=dfeeq*eeq;
    dfeeq=0.0;
    //double eq=eta/q;
    MY_DOUBLE_TYPE dfq=-dfeq*eta/(q*q);
    dfeq=0.0;
    //double q=eps_p+pow(p,p_exp);
#if !defined(NO_MY_DOUBLE_TYPE)
    dfp+=dfq*p_exp*pow(p,p_exp-1.0L);
#else
    dfp+=dfq*p_exp*pow(p,p_exp-1.0);
#endif
    dfq=0.0;
    return dfp;
  }
}

// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  


/*

int main()
{
   MY_DOUBLE_TYPE p=.3;
   MY_DOUBLE_TYPE eta=.1;
   cout << phi(p,eta) << endl;
   cout << phi(p,0.0) << endl;
   cout << phi(.6,0.0) << endl;
   df2_two_variable z=d2phi(p,eta);
   df3_two_variable w=d3phi(p,eta);

   MY_DOUBLE_TYPE delta=1.e-5;
   MY_DOUBLE_TYPE f11=
     (phi(p+delta,eta)-2.0*phi(p,eta) +phi(p-delta,eta))/(delta*delta);
   MY_DOUBLE_TYPE f22=
     (phi(p,eta+delta)-2.0*phi(p,eta) +phi(p,eta-delta))/(delta*delta);
   MY_DOUBLE_TYPE f1u=(phi(p+delta,eta+delta)- phi(p+delta,eta-delta))/(2.*delta);
   MY_DOUBLE_TYPE f1d=(phi(p-delta,eta+delta)- phi(p-delta,eta-delta))/(2.*delta);
   cout << "A " << (f1u-f1d)/(2.0*delta) << " " << z.get_u_xy() << endl;
   cout << "B " << f11 << " " << z.get_u_xx() << endl;
   cout << "C " << f22 << " " << z.get_u_yy() << endl;
   cout << "****************************" << endl;
   cout << "dphi_deta(p,eta)" <<endl;
   //cout << dphi_deta(p,eta)<<endl;
   
   cout << z.get_u_y() << " " 
        << (phi(p,eta+delta)- phi(p,eta-delta))/(2.*delta) << endl;
   cout << "****************************" << endl;
   cout << "dphi_dp(p,eta)" <<endl;
   //cout << dphi_dp(p,eta) <<endl;
   cout << (phi(p+delta,eta)- phi(p-delta,eta))/(2.*delta) << endl;
   cout << "u_x " << " " << z.get_u_x() << endl;
   cout << "u_y " << " " << z.get_u_y() << endl;
   cout << "u_xx " << " " << z.get_u_xx() << endl;
   cout << "u_xy " << " " << z.get_u_xy() << endl;
   cout << "u_yy " << " " << z.get_u_yy() << endl;

   cout << "****************************" << endl;
   cout << "d2phi_deta2(p,eta)" << endl << endl;
   //cout << d2phi_deta2(p,eta) <<  endl;
   cout << (phi(p,eta+delta)-2.0*phi(p,eta)
       + phi(p,eta-delta))/square(delta) << endl;
   cout << "****************************" << endl;

    delta=1.e-4;
    MY_DOUBLE_TYPE x1=(phi(p+delta,eta+delta)-2.0*phi(p+delta,eta)
       + phi(p+delta,eta-delta))/square(delta);
    MY_DOUBLE_TYPE x2=(phi(p-delta,eta+delta)-2.0*phi(p-delta,eta)
       + phi(p-delta,eta-delta))/square(delta);
   cout << "d2phi_deta2_1(p,eta)" << endl;
   cout << (x1-x2)/(2.0*delta) << endl;;
    delta=1.e-5;
   cout << "d2phi_deta2_2(p,eta)" << endl;
    delta=1.e-4;
    MY_DOUBLE_TYPE y1=(phi(p,eta+2.0*delta)-2.0*phi(p,eta+delta)
       + phi(p,eta))/square(delta);
    MY_DOUBLE_TYPE y2=(phi(p,eta)-2.0*phi(p,eta-delta)
       + phi(p,eta-2.0*delta))/square(delta);
    cout << (y1-y2)/(2.0*delta) << endl;
    delta=1.e-5;
   cout << d2phi_deta2_2(p,eta) << endl;
   cout << (d2phi_deta2(p,eta+delta) - 
          d2phi_deta2(p,eta-delta))/(2.*delta) << endl;
   cout << "u_x " << " " << *w.get_u_x() << endl;
   cout << "u_y " << " " << *w.get_u_y() << endl;
   cout << "u_xx " << " " << *w.get_u_xx() << endl;
   cout << "u_xy " << " " << *w.get_u_xy() << endl;
   cout << "u_yy " << " " << *w.get_u_yy() << endl;
   cout << "u_xxx " << " " << *w.get_u_xxx() << endl;
   cout << "u_xxy " << " " << *w.get_u_xxy() << endl;
   cout << "u_xyy " << " " << *w.get_u_xyy() << endl;
   cout << "u_yyy " << " " << *w.get_u_yyy() << endl;


   cout << "****************************" << endl;
   cout << "d2phi_deta2_1(p,eta)" << endl;
   cout << d2phi_deta2_1(p,eta) << endl;
   cout << (d2phi_deta2(p+delta,eta) - 
          d2phi_deta2(p-delta,eta))/(2.*delta) << endl;
   cout << "****************************" << endl;
   cout << "dphi_deta_1(p,eta)" << endl;
   cout << dphi_deta_1(p,eta) << endl;
   cout << (dphi_deta(p+delta,eta) - 
          dphi_deta(p-delta,eta))/(2.*delta) << endl;
   cout << "****************************" << endl;
   cout << "dphi_deta_2(p,eta)" << endl;
   cout << dphi_deta_2(p,eta) << endl;
   cout << (dphi_deta(p,eta+delta) - 
          dphi_deta(p,eta-delta))/(2.*delta) << endl;
   cout << "****************************" << endl;


}

*/
