/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

  //#include <df1b2fun.h>
  #include <admodel.h>
  #include "df12fun.h"
  //#include <df32fun.h>
  //#include <df3fun.h>
  #define EPS 3.0e-7
  #define FPMIN 1.0e-30
  #define MAXIT 100
  

  df1_two_variable betacf(df1_two_variable& a, df1_two_variable& b, 
    MY_DOUBLE_TYPE x);
  df1_two_variable betacf(df1_two_variable& a, MY_DOUBLE_TYPE b, df1_two_variable& x);
    
  df1_two_variable betacf(MY_DOUBLE_TYPE a,df1_two_variable& b,df1_two_variable&  x);
  df1_two_variable xgammln(const df1_two_variable& xx);
  dvariable xgammln(const prevariable& xx);
  MY_DOUBLE_TYPE xgammln(MY_DOUBLE_TYPE xx);

  static df1_two_variable xbetai(df1_two_variable& a,df1_two_variable& b,
    MY_DOUBLE_TYPE x, int maxit)
  {
    df1_two_variable bt;
  
#if !defined(NO_MY_DOUBLE_TYPE)
    if (x < 0.0 || x > 1.0L) cerr << "Bad x in routine xbetai" << endl;
#else
    if (x < 0.0 || x > 1.0) cerr << "Bad x in routine xbetai" << endl;
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
    if (x == 0.0 || x == 1.0L) bt=0.0;
#else
    if (x == 0.0 || x == 1.0) bt=0.0;
#endif
    else
      bt=exp(xgammln(a+b)-xgammln(a)-xgammln(b)+a*log(x)+b*log(1.0-x));
#if !defined(NO_MY_DOUBLE_TYPE)
    if (x < (value(a)+1.0L)/(value(a)+value(b)+2.0))
#else
    if (x < (value(a)+1.0)/(value(a)+value(b)+2.0))
#endif
      return bt*betacf(a,b,x)/a;
    else
    {
      MY_DOUBLE_TYPE xm1=1.0-x;
      return 1.0-bt*betacf(b,a,xm1)/b;
    }
  }
    
  static df1_two_variable xbetai(df1_two_variable& a,MY_DOUBLE_TYPE b,
    df1_two_variable& x,int maxit)
  {
    df1_two_variable bt;
  
#if !defined(NO_MY_DOUBLE_TYPE)
    if (value(x) < 0.0 || value(x) > 1.0L) cerr << "Bad x in routine xbetai" << endl;
#else
    if (value(x) < 0.0 || value(x) > 1.0) cerr << "Bad x in routine xbetai" << endl;
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
    if (value(x) == 0.0 || value(x) == 1.0L) bt=0.0;
#else
    if (value(x) == 0.0 || value(x) == 1.0) bt=0.0;
#endif
    else
      bt=exp(xgammln(a+b)-xgammln(a)-xgammln(b)+a*log(x)+b*log(1.0-x));
#if !defined(NO_MY_DOUBLE_TYPE)
    if (value(x) < (value(a)+1.0L)/(value(a)+b+2.0))
#else
    if (value(x) < (value(a)+1.0)/(value(a)+b+2.0))
#endif
      return bt*betacf(a,b,x)/a;
    else
    {
      df1_two_variable xm1=1.0-x;
      return 1.0-bt*betacf(b,a,xm1)/b;
    }
  }
    
  df1_two_variable xbetai(MY_DOUBLE_TYPE a,df1_two_variable& b,df1_two_variable& x,
    int maxit)
  {
    df1_two_variable bt;
  
#if !defined(NO_MY_DOUBLE_TYPE)
    if (value(x) < 0.0 || value(x) > 1.0L) cerr << "Bad x in routine xbetai" << endl;
#else
    if (value(x) < 0.0 || value(x) > 1.0) cerr << "Bad x in routine xbetai" << endl;
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
    if (value(x) == 0.0 || value(x) == 1.0L) bt=0.0;
#else
    if (value(x) == 0.0 || value(x) == 1.0) bt=0.0;
#endif
    else
      bt=exp(xgammln(a+b)-xgammln(a)-xgammln(b)+a*log(x)+b*log(1.0-x));
#if !defined(NO_MY_DOUBLE_TYPE)
    if (value(x) < (a+1.0L)/(a+value(b)+2.0))
#else
    if (value(x) < (a+1.0)/(a+value(b)+2.0))
#endif
      return bt*betacf(a,b,x)/a;
    else
    {
      df1_two_variable xm1=1.0-x;
      return 1.0-bt*betacf(b,a,xm1)/b;
    }
  }
    
  
  df1_two_variable betacf(df1_two_variable& a, df1_two_variable& b, 
    MY_DOUBLE_TYPE x)
  {
    int m,m2;
    df1_two_variable aa,c,d,del,h,qab,qam,qap;
     
    qab=a+b;
    qap=a+1.0;
    qam=a-1.0;
    c=1.0;
    d=1.0-qab*x/qap;
    if (fabs(value(d)) < FPMIN) d=FPMIN;
    d=1.0/d;
    h=d;
    for (m=1;m<=MAXIT;m++) 
    {
      m2=2*m;
      aa=m*(b-m)*x/((qam+m2)*(a+m2));
      d=1.0+aa*d;
      if (fabs(value(d)) < FPMIN) d=FPMIN;
      c=1.0+aa/c;
      if (fabs(value(c)) < FPMIN) c=FPMIN;
      d=1.0/d;
      h *= d*c;
      aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
      d=1.0+aa*d;
      if (fabs(value(d)) < FPMIN) d=FPMIN;
      c=1.0+aa/c;
      if (fabs(value(c)) < FPMIN) c=FPMIN;
      d=1.0/d;
  
      del=d*c;
      h *= del;
#if !defined(NO_MY_DOUBLE_TYPE)
      if (fabs(value(del)-1.0L) < EPS) break;
#else
      if (fabs(value(del)-1.0) < EPS) break;
#endif
    }
    if (m > MAXIT) 
    {
      cerr << "mum interations exceeded " << endl;
      ad_exit(1);
    }
    return h;
  }
  
  
  df1_two_variable betacf(MY_DOUBLE_TYPE a,df1_two_variable& b,df1_two_variable&  x)
  {
    int m,m2;
    df1_two_variable aa,c,d,del,h,qab,qam,qap;
     
    qab=a+b;
    qap=a+1.0;
    qam=a-1.0;
    c=1.0;
    d=1.0-qab*x/qap;
    if (fabs(value(d)) < FPMIN) d=FPMIN;
    d=1.0/d;
    h=d;
    for (m=1;m<=MAXIT;m++) 
    {
      m2=2*m;
      aa=m*(b-m)*x/((qam+m2)*(a+m2));
      d=1.0+aa*d;
      if (fabs(value(d)) < FPMIN) d=FPMIN;
      c=1.0+aa/c;
      if (fabs(value(c)) < FPMIN) c=FPMIN;
      d=1.0/d;
      h *= d*c;
      aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
      d=1.0+aa*d;
      if (fabs(value(d)) < FPMIN) d=FPMIN;
      c=1.0+aa/c;
      if (fabs(value(c)) < FPMIN) c=FPMIN;
      d=1.0/d;
  
      del=d*c;
      h *= del;
#if !defined(NO_MY_DOUBLE_TYPE)
      if (fabs(value(del)-1.0L) < EPS) break;
#else
      if (fabs(value(del)-1.0) < EPS) break;
#endif
    }
    if (m > MAXIT) 
    {
      cerr << "mum interations exceeded " << endl;
      ad_exit(1);
    }
    return h;
  }
  
  df1_two_variable betacf(df1_two_variable& a, MY_DOUBLE_TYPE b, df1_two_variable& x) 
  {
    int m,m2;
    df1_two_variable aa,c,d,del,h,qab,qam,qap;
     
    qab=a+b;
    qap=a+1.0;
    qam=a-1.0;
    c=1.0;
    d=1.0-qab*x/qap;
    if (fabs(value(d)) < FPMIN) d=FPMIN;
    d=1.0/d;
    h=d;
    for (m=1;m<=MAXIT;m++) 
    {
      m2=2*m;
      aa=m*(b-m)*x/((qam+m2)*(a+m2));
      d=1.0+aa*d;
      if (fabs(value(d)) < FPMIN) d=FPMIN;
      c=1.0+aa/c;
      if (fabs(value(c)) < FPMIN) c=FPMIN;
      d=1.0/d;
      h *= d*c;
      aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
      d=1.0+aa*d;
      if (fabs(value(d)) < FPMIN) d=FPMIN;
      c=1.0+aa/c;
      if (fabs(value(c)) < FPMIN) c=FPMIN;
      d=1.0/d;
  
      del=d*c;
      h *= del;
#if !defined(NO_MY_DOUBLE_TYPE)
      if (fabs(value(del)-1.0L) < EPS) break;
#else
      if (fabs(value(del)-1.0) < EPS) break;
#endif
    }
    if (m > MAXIT) 
    {
      cerr << "mum interations exceeded " << endl;
      ad_exit(1);
    }
    return h;
  }
  
static df1_two_variable xgammlnguts(const df1_two_variable& _z)
{
  //ADUNCONST(df1_two_variable,z)
  df1_two_variable z;
  z=_z;
  df1_two_variable x;
  const MY_DOUBLE_TYPE lpi =1.1447298858494001741434272;
  const MY_DOUBLE_TYPE pi =3.1415926535897932384626432;
  const MY_DOUBLE_TYPE lpp =0.9189385332046727417803297;
  int n=7;
  const MY_DOUBLE_TYPE c[9]={0.99999999999980993, 
    676.5203681218851, 
    -1259.1392167224028,
     771.32342877765313, 
    -176.61502916214059, 
    12.507343278686905,
     -0.13857109526572012, 
    9.9843695780195716e-6, 
    1.5056327351493116e-7};
  z-=0.0;
  z-=1.0;
  x=c[0];
  for (int i=1;i<=n+1;i++)
  {
    df1_two_variable zinv=1.0/(z+i);
    x+=c[i]*zinv;
  }    
  df1_two_variable t=z+n+0.5;
#if !defined(NO_MY_DOUBLE_TYPE)
  df1_two_variable ans= lpp + (z+0.5L)*log(t) -t + log(x);
#else
  df1_two_variable ans= lpp + (z+0.5)*log(t) -t + log(x);
#endif
  return(ans);
}

static MY_DOUBLE_TYPE xgammlnguts(MY_DOUBLE_TYPE z)
{
  MY_DOUBLE_TYPE x;
  const MY_DOUBLE_TYPE lpi =1.1447298858494001741434272;
  const MY_DOUBLE_TYPE pi =3.1415926535897932384626432;
  const MY_DOUBLE_TYPE lpp =0.9189385332046727417803297;
  int n=7;
  const MY_DOUBLE_TYPE c[9]={0.99999999999980993, 
    676.5203681218851, 
    -1259.1392167224028,
     771.32342877765313, 
    -176.61502916214059, 
    12.507343278686905,
     -0.13857109526572012, 
    9.9843695780195716e-6, 
    1.5056327351493116e-7};
  z-=0.0;
  z-=1.0;
  x=c[0];
  for (int i=1;i<=n+1;i++)
  {
    MY_DOUBLE_TYPE zinv=1.0/(z+i);
    x+=c[i]*zinv;
  }    
  MY_DOUBLE_TYPE t=z+n+0.5;
#if !defined(NO_MY_DOUBLE_TYPE)
  MY_DOUBLE_TYPE ans= lpp + (z+0.5L)*log(t) -t + log(x);
#else
  MY_DOUBLE_TYPE ans= lpp + (z+0.5)*log(t) -t + log(x);
#endif
  return(ans);
}

df1_two_variable xgammln(const df1_two_variable& z)
{
  const MY_DOUBLE_TYPE lpi =1.1447298858494001741434272;
  const MY_DOUBLE_TYPE pi =3.1415926535897932384626432;
#if !defined(NO_MY_DOUBLE_TYPE)
  if (value(z)<0.5L)
#else
  if (value(z)<0.5)
#endif
  {
    df1_two_variable zm1=1.0-z;
    return lpi - log(sin(pi*z)) - xgammlnguts(zm1);
  }
  else
  {
    return xgammlnguts(z);
  }
}
static MY_DOUBLE_TYPE xgammln(MY_DOUBLE_TYPE z)
{
  const MY_DOUBLE_TYPE lpi =1.1447298858494001741434272;
  const MY_DOUBLE_TYPE pi =3.1415926535897932384626432;
#if !defined(NO_MY_DOUBLE_TYPE)
  if (z<0.5L)
#else
  if (z<0.5)
#endif
  {
    return lpi - log(sin(pi*z)) - xgammlnguts(1.0-z);
  }
  else
  {
    return xgammlnguts(z);
  }
}
  
static dvariable xbetai(dvariable& _a,dvariable& _b,MY_DOUBLE_TYPE x,
    int maxit)
{
  init_df1_two_variable a(_a);
  init_df1_two_variable b(_b);
  df1_two_variable ret=xbetai(a,b,x,maxit);
  dvariable tmp;
  tmp=ret;
  return tmp;
}

dvariable xbetai(MY_DOUBLE_TYPE& a,dvariable& _b,dvariable& _x,
    int maxit)
{
  init_df1_two_variable b(_b);
  init_df1_two_variable x(_x);
  df1_two_variable ret=xbetai(a,b,x,maxit);
  dvariable tmp;
  tmp=ret;
  return tmp;
}
  
dvariable xbetai(dvariable& _a,MY_DOUBLE_TYPE& b,dvariable& _x,
    int maxit)
{
  init_df1_two_variable a(_a);
  init_df1_two_variable x(_x);
  df1_two_variable ret=xbetai(a,b,x,maxit);
  dvariable tmp;
  tmp=ret;
  return tmp;
}
  
//#include <df1b2fun.h>
#define ITMAX 100
#define EPS 1.0e-9
//#define EPS 3.0e-7
#define FPMIN 1.0e-30
MY_DOUBLE_TYPE get_values(MY_DOUBLE_TYPE x,MY_DOUBLE_TYPE y,int print_switch);


dvariable log_negbinomial_density(MY_DOUBLE_TYPE x,const dvariable& _xmu, 
  const dvariable& _xtau)
{
  ADUNCONST(dvariable,xmu)
  ADUNCONST(dvariable,xtau)
  init_df1_two_variable mu(xmu);
  init_df1_two_variable tau(xtau);
  if (value(tau)-1.0<0.0)
  {
    cerr << "tau <=1 in log_negbinomial_density " << endl;
    ad_exit(1);
  }
#if !defined(NO_MY_DOUBLE_TYPE)
  df1_two_variable r=mu/(1.e-120+(tau-1.0L));
#else
  df1_two_variable r=mu/(1.e-120+(tau-1.0));
#endif
  df1_two_variable tmp;
  tmp=xgammln(x+r)-xgammln(r) -xgammln(x+1)
    +r*log(r)+x*log(mu)-(r+x)*log(r+mu);
  dvariable tmp1;
  tmp1=tmp;
  return tmp1;
}



static void gcf(const df1_two_variable& _gammcf,const df1_two_variable& a,
  const df1_two_variable& x,const df1_two_variable& _gln)
{
  ADUNCONST(df1_two_variable,gln)
  ADUNCONST(df1_two_variable,gammcf)
  int i;
  df1_two_variable an,b,c,d,del,h;

  gln=xgammln(a);
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for (i=1;i<=ITMAX;i++) {
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(value(d)) < FPMIN) d=FPMIN;
    c=b+an/c;
    if (fabs(value(c)) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
#if !defined(NO_MY_DOUBLE_TYPE)
    if (fabs(value(del)-1.0L) < EPS) break;
#else
    if (fabs(value(del)-1.0) < EPS) break;
#endif
  }
  if (i > ITMAX) 
    cerr << "a too large, ITMAX too small in gcf" << endl;
  gammcf=exp(-x+a*log(x)-(gln))*h;
}

static void gser(const df1_two_variable& _gamser,const df1_two_variable& a,
  const df1_two_variable& x,const df1_two_variable& _gln)
{
  int n;
  ADUNCONST(df1_two_variable,gln)
  ADUNCONST(df1_two_variable,gamser)
  df1_two_variable sum,del,ap;

  gln=xgammln(a);

  if (value(x) <= 0.0) {
    if (value(x) < 0.0) 
      cerr << "x less than 0 in routine gser" << endl;
    gamser=0.0;
    return;
  } 
  else 
  {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ap+=1.0;
      del *= x/ap;
      sum += del;
      if (fabs(value(del)) < fabs(value(sum))*EPS) {
        gamser=sum*exp(-x+a*log(x)-(gln));
        return;
      }
    }
    cerr << "a too large, ITMAX too small in routine gser" << endl;
    return;
  }
}

dvariable xcumd_gamma(const prevariable&  x,const prevariable&  a)
{
  init_df1_two_variable xx(x);
  init_df1_two_variable aa(a);
  df1_two_variable zz=cumd_gamma(xx,aa);
  dvariable tmp;
  tmp=zz;
  return tmp;
}
 


df1_two_variable cumd_gamma(const df1_two_variable& x,
  const df1_two_variable& a)
{
  df1_two_variable gamser,gammcf,gln;

  if (value(x) < 0.0 || value(a) <= 0.0) 
    cerr << "Invalid arguments in routine gammp" << endl;
#if !defined(NO_MY_DOUBLE_TYPE)
  if (value(x) < (value(a)+1.0L)) {
#else
  if (value(x) < (value(a)+1.0)) {
#endif
    gser(gamser,a,x,gln);
    return gamser;
  } else {
    gcf(gammcf,a,x,gln);
    return 1.0-gammcf;
  }
}

dvariable& dvariable::operator = (const df1_two_variable& v)
{
  const prevariable * px=df1_two_variable::ind_var[0];
  const prevariable * py=df1_two_variable::ind_var[1];
  MY_DOUBLE_TYPE  dfx= *v.get_u_x();
  MY_DOUBLE_TYPE  dfy= *v.get_u_y();
  value(*this)=*v.get_u();

  gradient_structure::GRAD_STACK1->set_gradient_stack(default_evaluation3,
    &(value(*this)),&(value(*px)),dfx,&(value(*py)),dfy);

  return *this;
}

df1_two_variable betaln(const df1_two_variable & a,const df1_two_variable& b)
{
  return xgammln(a)+xgammln(b)-xgammln(a+b);
}

dvariable betaln(const prevariable & _xa,const prevariable& _xb)
{
  ADUNCONST(prevariable,xa)
  ADUNCONST(prevariable,xb)
  init_df1_two_variable a(xa);
  init_df1_two_variable b(xb);
  dvariable tmp;
  tmp=betaln(a,b);
  return tmp;
}

  dvariable ln_beta_density(MY_DOUBLE_TYPE y,const prevariable & xmu,
    const prevariable& xphi)
  {
    init_df1_two_variable mu(xmu);
    init_df1_two_variable phi(xphi);
    df1_two_variable omega=mu*phi;
    df1_two_variable tau=phi-mu*phi;
    df1_two_variable lb=betaln(omega,tau);
    df1_two_variable d=(omega-1)*log(y)+(tau-1)*log(1.0-y)-lb;
    dvariable tmp;
    tmp=d;
    return tmp;
  }



