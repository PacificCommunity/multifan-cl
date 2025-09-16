/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#define USE_DD_NOT
#if defined(USE_DD)
#  include <qd/fpu.h>
#endif
#include "all.hpp"

extern int ctlc_flag;
void test_fpu(void);
void positivize(const banded_symmetric_dmatrix& _m,MY_DOUBLE_TYPE id);

MY_DOUBLE_TYPE get_init_lambda_l(const banded_symmetric_dmatrix& m,dvector& g,
  MY_DOUBLE_TYPE delta);
MY_DOUBLE_TYPE get_init_lambda_u(const banded_symmetric_dmatrix& m,dvector& g,
  MY_DOUBLE_TYPE delta);
banded_symmetric_dmatrix positivize
  (MY_DOUBLE_TYPE id,const banded_symmetric_dmatrix& _m);

extern "C" /* _export*/  void dd_newton_raphson(int n,MY_DOUBLE_TYPE * v,MY_DOUBLE_TYPE * diag,
    MY_DOUBLE_TYPE * ldiag, MY_DOUBLE_TYPE * yy);
//extern "C" _export  void dd_newton_raphson(int n,MY_DOUBLE_TYPE * v,MY_DOUBLE_TYPE * diag,
 //   MY_DOUBLE_TYPE * ldiag, MY_DOUBLE_TYPE * yy);
dmatrix get_dmatrix_from_master(void);
void send_dvector_to_master(const dvector& v);

MY_DOUBLE_TYPE get_new_delta(MY_DOUBLE_TYPE ftmp,MY_DOUBLE_TYPE f0,MY_DOUBLE_TYPE mf,MY_DOUBLE_TYPE old_delta,
  MY_DOUBLE_TYPE ns)
{
  MY_DOUBLE_TYPE delta=old_delta;
  if (f0<ftmp)
  {
    delta*=0.5;
    if (ns<delta)
    {
      delta=0.5*ns;
    }
  }
  else
  {
    MY_DOUBLE_TYPE tmp=(ftmp-f0)/(1.e-100+mf);
    if (tmp>0.25)
    {
      delta*=2.0;
    }
  }
  return delta;
}

dvector solve_secular_equation(MY_DOUBLE_TYPE * pmf,dvector & g,
  banded_symmetric_dmatrix &H,MY_DOUBLE_TYPE delta,int * pfirst_try,
  MY_DOUBLE_TYPE & ns)
{
  int i;
  //bracket the solution  TRM pg 192
  MY_DOUBLE_TYPE lammin=get_init_lambda_l(H,g,delta);
  MY_DOUBLE_TYPE lammax=get_init_lambda_u(H,g,delta);
  //double lammin=-1.e+100;;
  //double lammax=+1.e+100;
  // need upper and lower bounds on lambda for bracketing
  // for now increase diagonal until choleski works
  int ierr=0;
  banded_symmetric_dmatrix S;
  MY_DOUBLE_TYPE lambda0=0;
  MY_DOUBLE_TYPE wt=pow(0.5,double(*pfirst_try));
  MY_DOUBLE_TYPE lambda=(lammin+wt*(lammax-lammin))/(1.0+wt);
  int posflag=0;
  ierr=-1;
  S=positivize(lambda,H);
  banded_lower_triangular_dmatrix L=choleski_decomp_trust_bound(S,ierr);
  while (ierr==1) 
  {
    posflag=1;
    lammin=lambda+L(1,1);
    (*pfirst_try)-=1;
     wt=pow(0.5,double(*pfirst_try));
    lambda=(lammin+wt*(lammax-lammin))/(1.0+wt);
    lambda=max(lammin,lambda);
    S=positivize(lambda,H);
    ierr=-1;
    L=choleski_decomp_trust_bound(S,ierr);
  }
  int check_alg=0;
  if (check_alg)
  {
    do
    {
      cin >> lambda;
      S=positivize(lambda,H);
      ierr=-1;
      if (ierr<=0)
      {
        L=choleski_decomp_trust_bound(S,ierr);
        dvector s=solve_trans(L,solve(L,-g));
        MY_DOUBLE_TYPE ns=norm(s);
        cout <<"mfclimplicit.cpp " << ns << endl;
      }
    }
    while(1);
  }
  
  dvector mg=-g;
  dvector s=solve_trans(L,solve(L,mg));
  ns=norm(s);
  // check if newton step succeeded
  if (posflag==0)
  {
    (*pfirst_try)+=1;
  }
  else
  {
    (*pfirst_try)-=1;
  }
  if (ns<delta && posflag==0)
  {
    *pmf=(g+H*(0.5*s))*s;
    return s;
  }
  int iconv=0;
  int mmax=S.indexmax();
  MY_DOUBLE_TYPE keasy=0.05;
  
  MY_DOUBLE_TYPE lamold=0;
  do
  { 
    if (iconv) 
      break;

    if (fabs(ns-delta)<keasy*delta)
    {
      // *pmf=g*s+0.5*(s*(H*s));
      *pmf=(g+H*(0.5*s))*s;
      iconv=1.0;
      return s;
    }
    else
    {
      if (ns>delta)  // lambda is too small so we are in L
      {
        lammin=lambda;
      }
      if (ns<delta)  // lambda is too large so we are in G
      {
        lammax=lambda;
      }
    }

    //int ierr=0;
    //banded_lower_triangular_dmatrix L=choleski_decomp(S,ierr);
    {
      dvector wv =solve(L,s);
      ns=norm(s);
      MY_DOUBLE_TYPE lamtmp=lambda+((ns-delta)/delta)* square(ns)/norm2(wv);
      if (lamtmp<lammin)
      {
        lambda=avg(lammin,lambda);
      }
      else if (lamtmp>lammax)
      {
        lambda=avg(lambda,lammax);
      }
      else
      {
        lamold=lambda;
        lambda=lamtmp;
      }
      
      do 
      {
        S=positivize(lambda,H);
        ierr=-1;
        
        L=choleski_decomp(S,ierr);
        if (ierr<=0)
        {
          s=solve_trans(L,solve(L,mg));
          ns=norm(s);
          break;
        }
        else
        {
          lammin=lambda;
          //lambda=0.5*(lammin+lamold);
          lambda=0.5*(lammin+lammax);
        }
      }
      while(1);
    }

  }
  while(1);
}

dvector safe_solve(const banded_symmetric_dmatrix& _m,const dvector&_v,
  MY_DOUBLE_TYPE id);

int minindex(MY_DOUBLE_TYPE f,MY_DOUBLE_TYPE f2,MY_DOUBLE_TYPE f3)
{
  if (f<f2)
  {
    if (f<f3)
      return 1;
  }
  else
  {
    if (f2<f3)
      return 2;
  }
  return 3;
}
    
static MY_DOUBLE_TYPE min(MY_DOUBLE_TYPE f,MY_DOUBLE_TYPE f2,MY_DOUBLE_TYPE f3)
{
  if (f<f2)
  {
    if (f<f3)
      return f;
  }
  else
  {
    if (f2<f3)
      return f2;
  }
  return f3;
}


MY_DOUBLE_TYPE poly(MY_DOUBLE_TYPE a,MY_DOUBLE_TYPE b,MY_DOUBLE_TYPE c,MY_DOUBLE_TYPE x)
{
  return c*x*x+b*x+a;
}

MY_DOUBLE_TYPE dpoly(MY_DOUBLE_TYPE a,MY_DOUBLE_TYPE b,MY_DOUBLE_TYPE c,MY_DOUBLE_TYPE x)
{
  return 2.0*c*x+b;
}

dvector get_coffs( MY_DOUBLE_TYPE f1, MY_DOUBLE_TYPE f2, MY_DOUBLE_TYPE pf1, MY_DOUBLE_TYPE pf2,MY_DOUBLE_TYPE h)
{
  MY_DOUBLE_TYPE ap=f1;
  MY_DOUBLE_TYPE bp=(f2-f1)/h+0.5*(pf1-pf2);
  MY_DOUBLE_TYPE cp=(f2-f1-bp*h)/(h*h);
  dvector tmp(1,3);
  tmp(1)=ap;
  tmp(2)=bp;
  tmp(3)=cp;
  return tmp; 
}

dvector get_coffs2( MY_DOUBLE_TYPE f1, MY_DOUBLE_TYPE f2, MY_DOUBLE_TYPE pf1,MY_DOUBLE_TYPE h)
{
  MY_DOUBLE_TYPE ap=f1;
  MY_DOUBLE_TYPE bp=pf1;
  MY_DOUBLE_TYPE cp=(f2-ap-bp*h)/(h*h);
  dvector tmp(1,3);
  tmp(1)=ap;
  tmp(2)=bp;
  tmp(3)=cp;
  return tmp; 
}

MY_DOUBLE_TYPE get_min(dvector& tmp)
{
  return -tmp(2)/(2.0*tmp(3));
}


class fmmtx : public fmm_control
{
private:
  dvector w;
  dvector funval;
  int xm;
  dmatrix xstep;
  dvector xrho;
  dvector rrr;
  dmatrix xy; 
  dvector xold; 
  dvector gold; 
public:
  MY_DOUBLE_TYPE dmin,fbest,df;

  long int llog,n1,ic,iconv,i1,link;
  MY_DOUBLE_TYPE z,zz,gys,gs,sig,gso,alpha,tot,fy,dgs;
  long int itn,icc,np,nn,is,iu,iv,ib;
  int i, j;
  MY_DOUBLE_TYPE gmax;
  MY_DOUBLE_TYPE fsave;
  dvector xx;
  dvector gbest;
  dvector xsave;
  dvector gsave;

  int n;

public:
  fmmtx(int nvar,int _xm=7);
 /*
  fmmt1(int nvar,_CONST lvector& ipar);
  MY_DOUBLE_TYPE minimize(BOR_CONST independent_variables & x,MY_DOUBLE_TYPE (*pf)(_CONST dvar_vector&));

  MY_DOUBLE_TYPE minimize(BOR_CONST independent_variables & x,BOR_CONST dvector& c,
        MY_DOUBLE_TYPE (*pf)(BOR_CONST dvar_vector&,BOR_CONST dvector&) );

  void fmin2(BOR_CONST MY_DOUBLE_TYPE& f, BOR_CONST independent_variables & x,BOR_CONST dvector& g, function_minimizer *);

  void fmin(BOR_CONST MY_DOUBLE_TYPE& f, const dvector & x,BOR_CONST dvector& g);

 */
//  dmatrix& hessian();
};



fmmtx::fmmtx(int nvar,int _xm)
: w(1,4*nvar), funval(1,10) , 
  xx(0,nvar) /* , gbest(0,nvar) ,xsave(0,nvar), gsave(0,nvar) ,
  xstep(0,_xm+1,1,nvar) , xy(0,_xm+1,1,nvar), xrho(0,_xm+1), 
  xold(1,nvar), gold(1,nvar), rrr(1,nvar)*/ 
{
/*
  ctlc_flag = 0;
  n = nvar;
  xm=_xm;
  xrho.initialize();
*/
//  cout << " In fmm::fmm(int nvar) nvar = " << nvar 
//       << " and n = " << n << "\n";
}


void get_initial_x(dvector& Y, dvector& x)
{
  int n=Y.indexmax();
  for (int i=1;i<=n;i++)
  {
    int mmin=max(i-6,1);
    int mmax=min(i+6,n);
    x(i)=mean(Y(mmin,mmax));
  }
}

void dvar_fish_stock_history::get_f(MY_DOUBLE_TYPE & f,dmatrix& Y,dvector& x,
  dvector& g,banded_symmetric_dmatrix& H,int Hessflag,int robflag,
  dmatrix& resids,dmatrix & wts,int wtflag,
  ivector& iswitch,imatrix& indx,MY_DOUBLE_TYPE c)
{
 
  get_f_normal(f,Y,x,g,H,Hessflag,robflag,resids,wts,wtflag,iswitch,
    indx,c);
  return;
 
  {
    int jj=1;
    dvector fvec(1,10000);
    //ofstream ofs("imp"); 
    MY_DOUBLE_TYPE s_eps=0.3;
    MY_DOUBLE_TYPE pen1=1.0/square(s_eps);
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
      else
      {
        grouped_catchability_coffs(j).initialize();
      }
    }
    
    ivector grouping=column(fish_flags,29);
    
    for (j=1;j<=num_fisheries;j++)
    {
      if (fish_flags(j,10)==1)
      {
        if (robflag==1)
        {
          for (i=1;i<=num_real_fish_times(j);i++)
          {
            //double r=Y(i)-x(i);
            ii=indx(grouping(j),gfish_ptr(j,i)); 
      
            grouped_catchability_coffs(grouping(j),gfish_ptr(j,i));
      
            MY_DOUBLE_TYPE r=Y(j,i)-grouped_catchability_coffs(grouping(j),gfish_ptr(j,i));
            resids(j,i)=r;
            MY_DOUBLE_TYPE r2=square(r);
            /*********************************************************/
            //double u=exp(-0.5*pen1*r2)+c/(1+0.5*pen1*r2);
            MY_DOUBLE_TYPE u1=exp(-0.5*pen1*r2);
            MY_DOUBLE_TYPE v1=1.0/(1+0.5*pen1*r2);
            MY_DOUBLE_TYPE u=u1+c*v1;
            /*********************************************************/
  
            //ofs << setprecision(14) << u << " " << r2 << endl;
            f-=log(u);
            fvec(jj++)=-log(u);
            /*********************************************************/
            //double tmp1= exp(-0.5*pen1*r2)+c/square(1.0+0.5*pen1*r2);
  
            MY_DOUBLE_TYPE tmp1= u1+c*square(v1);
  
            /*********************************************************/
            MY_DOUBLE_TYPE tmp= tmp1*pen1*r;
            
            g(ii)-=1.0/u*tmp;
           
            MY_DOUBLE_TYPE testh;
            MY_DOUBLE_TYPE eps=1.e-6;
            MY_DOUBLE_TYPE gt=0.0;
            MY_DOUBLE_TYPE rt=r+eps;;
            {
              MY_DOUBLE_TYPE r=rt;
              MY_DOUBLE_TYPE r2=square(r);
              MY_DOUBLE_TYPE u=exp(-0.5*pen1*r2)+c/(1+0.5*pen1*r2);
              MY_DOUBLE_TYPE tmp1= exp(-0.5*pen1*r2)+c/square(1.0+0.5*pen1*r2);
              MY_DOUBLE_TYPE tmp= tmp1*pen1*r;
              gt=-1.0/u*tmp;
            }
            
            if (Hessflag)
            {
              MY_DOUBLE_TYPE a1= 1.0/square(u)*square(tmp);
              MY_DOUBLE_TYPE a2= 1.0/u*tmp1*pen1;
              MY_DOUBLE_TYPE a3= 1.0/u*(u1+2.0*c*cube(v1))
                  *square(pen1*r);
             /*********************************************************/
             // MY_DOUBLE_TYPE a3= 1.0/u*(exp(-0.5*pen1*r2)+2.0*c/cube(1.0+0.5*pen1*r2))
             //    *square(pen1*r);
              /*
              H(ii,ii)+=a1;
              H(ii,ii)+=a2;
              H(ii,ii)-=a3;
              */
              H(ii,ii)+=a1+a2-a3;
              cout << H(ii,ii) << " " << (gt-g(ii))/eps << endl;
            }
          }
        }
        else
        {
          for (i=1;i<=num_real_fish_times(j);i++)
          {
            ii=indx(grouping(j),gfish_ptr(j,i)); 
      
            grouped_catchability_coffs(grouping(j),gfish_ptr(j,i));
      
            MY_DOUBLE_TYPE r=Y(j,i)-grouped_catchability_coffs(grouping(j),gfish_ptr(j,i));
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
  
    ii=1;
    for (j=1;j<=ngroups;j++)
    {
      if (fish_flags(gfish_index(j,1),10)==1)
      {
        ii++;
        for (i=2;i<=num_real_grouped_fish_times(j);i++)
        {
          //double r=(x(i)-x(i-1));
          MY_DOUBLE_TYPE pp2=pen2*grouped_between_times(j,i);
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
#if ( defined(USE_DD) && ( (defined(__MSVC32__) &&  __MSVC32__ >=8) || defined(linux) || defined(__ADMING__)))
    f=to_double(df);
#else
    f=df;
#endif
  }
  //cout <<"mfclimplicit.cpp " << df-f << endl;
}

dvariable dvar_fish_stock_history::
  grouped_implicit_catchability_deviations_calc(void)
{
  const MY_DOUBLE_TYPE tpi=2.e0*3.14159e0;
  int n=0;
  int mmin=gfish_index.indexmin();
  int mmax=gfish_index.indexmax();
  
  ivector iswitch(mmin,mmax);
  int i,j;
  ivector ff10=column(fish_flags,10);
  for (i=mmin;i<=mmax;i++)
  {
    int rmin=gfish_index(i).indexmin();
    int rmax=gfish_index(i).indexmax();
    int fv=fish_flags(gfish_index(i,rmin),10);
    int fv1=ff10(gfish_index(i,rmin));
    if (fv!=fv1) cerr << "error" << endl;
    for (int j=rmin+1;j<=rmax;j++)
    {
      if (fv != fish_flags(gfish_index(i,j),10))
      {
        cerr << "sanity error in ff(i,10)" << 
         " for group " << i << " and fishery " << gfish_index(i,j) 
         << endl;
         ad_exit(1);
      }
    }
    iswitch(i)=fv;
    if (fv==1)
    {
      n+=num_real_grouped_fish_times(i);
    }
  }
  if (n==0)
  {
    cerr << "Error can not have all fish flags 10 =0 with"
      " af(104)=1" << endl;
    ad_exit(1);
  }
  imatrix indx(1,ngroups,1,num_real_grouped_fish_times);
  indx.initialize();
  int ii=1;
  for (j=1;j<=ngroups;j++)
  {
    if (iswitch(j))
    {
      for (i=1;i<=num_real_grouped_fish_times(j);i++)
      {
        indx(j,i)=ii++;
      }
    }
  }
  dvector x(1,n);
  MY_DOUBLE_TYPE c=0.05;
  MY_DOUBLE_TYPE f=grouped_implicit_catchability_deviations_calc_part1(n,x,
    iswitch,indx,c);
  dvariable vf=grouped_implicit_catchability_deviations_calc_part2(x,
    iswitch,indx,c);
  cout << "fabs(c-v) = " << fabs(f-value(vf)) << endl;
  
  int offset=0;
  int i27=sum(column(fish_flags,27));
  for (j=1;j<=ngroups;j++)
  {
    if (iswitch(j))
    {
      grouped_catch_dev_coffs(j)=
       first_difference
        (x(1+offset,num_real_grouped_fish_times(j)+offset)).shift(2);
      if (i27)
      {
        i=gfish_index(j,1);
        dvar_vector seasonal=
          fish_pars(1,i)*sin(tpi*(grouped_true_month(j)/12.-fish_pars(2,i)));
        grouped_catchability_coffs(j)+=value(seasonal);
        grouped_catch_dev_coffs(j)+=first_difference(seasonal).shift(2);
      }
      int mmin=gfish_index(j).indexmin();
      int mmax=gfish_index(j).indexmax();
      for (i=mmin;i<=mmax;i++)
      { 
        q0(gfish_index(j,i))=x(1+offset);
      }
      offset+=num_real_grouped_fish_times(j);
    }
    else
    {
      //grouped_catchability_coffs(j).initialize();
    }
  }
  return vf;
}

void dvar_fish_stock_history::
  grouped_implicit_catchability_deviations_calc_part1_thread(void)
{
  int n=0;
  int mmin=gfish_index.indexmin();
  int mmax=gfish_index.indexmax();
  
  ivector iswitch(mmin,mmax);
  int i,j;
  for (i=mmin;i<=mmax;i++)
  {
    int rmin=gfish_index(i).indexmin();
    int rmax=gfish_index(i).indexmax();
    int fv=fish_flags(gfish_index(i,rmin),10);
    for (int j=rmin+1;j<=rmax;j++)
    {
      if (fv != fish_flags(gfish_index(i,j),10))
      {
        cerr << "sanity error in ff(i,10)" << 
         " for group " << i << " and fishery " << gfish_index(i,j) 
         << endl;
         ad_exit(1);
      }
    }
    iswitch(i)=fv;
    if (fv==1)
    {
      n+=num_real_grouped_fish_times(i);
    }
  }
  if (n==0)
  {
    cerr << "Error can not have all fish flags 10 =0 with"
      " af(104)=1" << endl;
    ad_exit(1);
  }
  imatrix indx(1,ngroups,1,num_real_grouped_fish_times);
  indx.initialize();
  int ii=1;
  for (j=1;j<=ngroups;j++)
  {
    if (iswitch(j))
    {
      for (i=1;i<=num_real_grouped_fish_times(j);i++)
      {
        indx(j,i)=ii++;
      }
    }
  }
  dvector x(1,n);
  MY_DOUBLE_TYPE c=0.05;
  MY_DOUBLE_TYPE f=grouped_implicit_catchability_deviations_calc_part1(n,x,
    iswitch,indx,c);
  if (!allocated(thread_xsave))
  {
    thread_xsave.allocate(1,n);
  }
  thread_xsave=x;
  thread_f=f;
  //dvariable vf=grouped_implicit_catchability_deviations_calc_part2(x,
  //  iswitch,indx,c);
  //cout <<"mfclimplicit.cpp " << fabs(f-value(vf)) << endl;
  //return vf;
}
dvariable dvar_fish_stock_history::
  grouped_implicit_catchability_deviations_calc_part2_thread(void)
{
  int n=0;
  int mmin=gfish_index.indexmin();
  int mmax=gfish_index.indexmax();
  
  ivector iswitch(mmin,mmax);
  int i,j;
  for (i=mmin;i<=mmax;i++)
  {
    int rmin=gfish_index(i).indexmin();
    int rmax=gfish_index(i).indexmax();
    int fv=fish_flags(gfish_index(i,rmin),10);
    for (int j=rmin+1;j<=rmax;j++)
    {
      if (fv != fish_flags(gfish_index(i,j),10))
      {
        cerr << "sanity error in ff(i,10)" << 
         " for group " << i << " and fishery " << gfish_index(i,j) 
         << endl;
         ad_exit(1);
      }
    }
    iswitch(i)=fv;
    if (fv==1)
    {
      n+=num_real_grouped_fish_times(i);
    }
  }
  if (n==0)
  {
    cerr << "Error can not have all fish flags 10 =0 with"
      " af(104)=1" << endl;
    ad_exit(1);
  }
  imatrix indx(1,ngroups,1,num_real_grouped_fish_times);
  indx.initialize();
  int ii=1;
  for (j=1;j<=ngroups;j++)
  {
    if (iswitch(j))
    {
      for (i=1;i<=num_real_grouped_fish_times(j);i++)
      {
        indx(j,i)=ii++;
      }
    }
  }
  //dvector x(1,n);
  MY_DOUBLE_TYPE c=0.05;
  //double f=grouped_implicit_catchability_deviations_calc_part1(n,x,
  //  iswitch,indx,c);
  //if (!allocated(thread_xsave))
  //{
  //  thread_xsave.allocate(1,n);
  //}
  //thread_xsave=x;
  dvariable vf=grouped_implicit_catchability_deviations_calc_part2
    (thread_xsave,iswitch,indx,c);
  cout <<"mfclimplicit.cpp " << fabs(thread_f-value(vf)) << endl;
  return vf;
}
#if defined(USE_ADPVM)
dvariable dvar_fish_stock_history::
  grouped_implicit_catchability_deviations_calc_part2_pvm(void)
{
  const MY_DOUBLE_TYPE tpi=2.e0*3.14159e0;
  int n=0;
  int mmin=gfish_index.indexmin();
  int mmax=gfish_index.indexmax();
  
  ivector iswitch(mmin,mmax);
  int i,j;
  for (i=mmin;i<=mmax;i++)
  {
    int rmin=gfish_index(i).indexmin();
    int rmax=gfish_index(i).indexmax();
    int fv=fish_flags(gfish_index(i,rmin),10);
    for (int j=rmin+1;j<=rmax;j++)
    {
      if (fv != fish_flags(gfish_index(i,j),10))
      {
        cerr << "sanity error in ff(i,10)" << 
         " for group " << i << " and fishery " << gfish_index(i,j) 
         << endl;
         ad_exit(1);
      }
    }
    iswitch(i)=fv;
    if (fv==1)
    {
      n+=num_real_grouped_fish_times(i);
    }
  }
  if (n==0)
  {
    cerr << "Error can not have all fish flags 10 =0 with"
      " af(104)=1" << endl;
    ad_exit(1);
  }
  imatrix indx(1,ngroups,1,num_real_grouped_fish_times);
  indx.initialize();
  int ii=1;
  for (j=1;j<=ngroups;j++)
  {
    if (iswitch(j))
    {
      for (i=1;i<=num_real_grouped_fish_times(j);i++)
      {
        indx(j,i)=ii++;
      }
    }
  }
  //dvector x(1,n);
  MY_DOUBLE_TYPE c=0.05;
  thread_xsave=mfget_dvector_from_slave(1);
  dvariable vf=grouped_implicit_catchability_deviations_calc_part2
    (thread_xsave,iswitch,indx,c);
  cout <<"mfclimplicit.cpp " << fabs(thread_f-value(vf)) << endl;
  int offset=0;
  dvector& x =thread_xsave;
  int i27=sum(column(fish_flags,27));
  for (j=1;j<=ngroups;j++)
  {
    if (iswitch(j))
    {
      grouped_catch_dev_coffs(j)=
       first_difference
        (x(1+offset,num_real_grouped_fish_times(j)+offset)).shift(2);
      if (i27)
      {
        i=gfish_index(j,1);
        dvar_vector seasonal=
          fish_pars(1,i)*sin(tpi*(grouped_true_month(j)/12.-fish_pars(2,i)));
        grouped_catchability_coffs(j)+=value(seasonal);
        grouped_catch_dev_coffs(j)+=first_difference(seasonal).shift(2);
      }
      int mmin=gfish_index(j).indexmin();
      int mmax=gfish_index(j).indexmax();
      for (i=mmin;i<=mmax;i++)
      { 
        q0(gfish_index(j,i))=x(1+offset);
      }
      offset+=num_real_grouped_fish_times(j);
    }
    else
    {
      //grouped_catchability_coffs(j).initialize();
    }
  }
  return vf;
}
#endif

dvariable dvar_fish_stock_history::get_vf(dvar_matrix& Y,dvector& x,
  ivector& iswitch,imatrix& indx)
{
  //ofstream ofs("vimp");
  cerr << "error we should not be here" << endl;
  ad_exit(1);
  dvariable vf;
  vf=0.0;
  
  MY_DOUBLE_TYPE s_eps=0.3;
  MY_DOUBLE_TYPE pen1=1.0/square(s_eps);
  MY_DOUBLE_TYPE s_eta=0.1;
  MY_DOUBLE_TYPE pen2=1.0/(2.0*square(s_eta));
  int n=x.indexmax();
  MY_DOUBLE_TYPE c=0.05;
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
  ii=1;
  for (j=1;j<=num_fisheries;j++)
  {
    if (iswitch(grouping(j)))
    {
      for (i=1;i<=num_real_fish_times(j);i++)
      {
        dvariable r=Y(j,i)-grouped_catchability_coffs(grouping(j),
          gfish_ptr(j,i));
        dvariable r2=square(r);
        dvariable u=exp(-0.5*pen1*r2)+c/(1+0.5*pen1*r2);
        vf-=log(u);
        //ofs << setprecision(14) << u << " " << r2 << endl;
        //cout << j << " " << i << log(u) << endl;
      }
    }
  }
  //cout << "AA " << vf << endl;
  //exit(1);

  ii=1;
  for (j=1;j<=ngroups;j++)
  {
    if (iswitch(j))
    {
      ii++;
      for (i=2;i<=num_real_grouped_fish_times(j);i++)
      {
        MY_DOUBLE_TYPE pp2=pen2*grouped_between_times(j,i);
        if (pp2==0)
        {
          cerr << "this can't happen" << endl;
          ad_exit(1);
        }
        dvariable r=grouped_catchability_coffs(j,i)
          -grouped_catchability_coffs(j,i-1);
        vf+=0.5*pp2*square(r);
      }
    }
  }
  //cout << "BB " << vf << endl;
  return vf;
}

/*
void dvar_fish_stock_history::grouped_catchability_calc(void)
{
  int i;
  ivector grouping=column(fish_flags,29);
  for (i=1;i<=ngroups;i++)
  {
    for (int nt=2;nt<=num_grouped_fish_times(i);nt++)
    {
      grouped_catch_dev_coffs(i,nt)=grouped_catchability_coffs(i,nt)
        -grouped_catchability_coffs(i,nt-1);
    }
  }
}
*/
banded_symmetric_dmatrix positivize
  (MY_DOUBLE_TYPE id,const banded_symmetric_dmatrix& _m)
{
  ADUNCONST(banded_symmetric_dmatrix,m)
  banded_symmetric_dmatrix tmp;

  tmp=m;

  int mmin=m.indexmin();
  int mmax=m.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    tmp(i,i)+=id;
  }
  return tmp;
}

void positivize(const banded_symmetric_dmatrix& _m,MY_DOUBLE_TYPE id)
{
  ADUNCONST(banded_symmetric_dmatrix,m)
  int mmin=m.indexmin();
  int mmax=m.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    m(i,i)+=id;
  }
}
banded_lower_triangular_dmatrix quiet_choleski_decomp(
  const banded_symmetric_dmatrix& _M,const int& _ierr)
{
  int & ierr = (int &) _ierr;
  ADUNCONST(banded_symmetric_dmatrix,M)
  int minsave=M.indexmin();
  M.shift(1);
  int n=M.indexmax();
  
  int bw=M.bandwidth();
  banded_lower_triangular_dmatrix L(1,n,bw);
#ifndef SAFE_INITIALIZE
    L.initialize();
#endif

  int i,j,k;
  MY_DOUBLE_TYPE tmp;
    if (M(1,1)<=0)
    {
      ierr=1;
      return L;
    }
  L(1,1)=sqrt(M(1,1));
  for (i=2;i<=bw;i++)
  {
    L(i,1)=M(i,1)/L(1,1);
  }

  for (i=2;i<=n;i++)
  {
    for (j=i-bw+1;j<=i-1;j++)
    {
      if (j>1)
      {	
        tmp=M(i,j);
        for (k=i-bw+1;k<=j-1;k++)
        {
	  if (k>0 && k>j-bw)
            tmp-=L(i,k)*L(j,k);
        }
        L(i,j)=tmp/L(j,j);
      }
    }
    tmp=M(i,i);
    for (k=i-bw+1;k<=i-1;k++)
    {
      if (k>0)	
        tmp-=L(i,k)*L(i,k);
    }
    if (tmp<=0)
    {
      ierr=1;
      return L;
    }
    L(i,i)=sqrt(tmp);
  }
  M.shift(minsave);
  L.shift(minsave);

  return L;
}

dvector safe_solve(const banded_symmetric_dmatrix& _m,const dvector&_v,
  MY_DOUBLE_TYPE id)
{
  int ierr=0;
  ADUNCONST(dvector,v)
  ADUNCONST(banded_symmetric_dmatrix,m)
  int mmin=m.indexmin();
  int mmax=m.indexmax();
  if (id>0.0)
  {
    positivize(m,id);
  }
  m.shift(1);
  v.shift(1);
  int ibreak=1;
  dvector w;
  do
  {
    const banded_lower_triangular_dmatrix& C=quiet_choleski_decomp(m,ierr);
    if (ierr==0)
    {
      w=solve_trans(C,solve(C,v));
      ibreak=0;
    }
    else
    {
      positivize(m,25.0);
      ierr=0;
    }
  }
  while(ibreak);
  m.shift(mmin);
  w.shift(mmin);
  v.shift(mmin);
  return w;
}


void dvar_fish_stock_history::cimplicit_seasonal_catchability_calc
  (dmatrix& Y)
{
  int i,it;
  const MY_DOUBLE_TYPE tpi=2.e0*3.14159e0;
  for (i=1;i<=num_fisheries;i++)
  {
    int rr=realization_region(i,1);
    if (fish_flags(i,27)==1)
    {
      for (it=1;it<=num_real_fish_times(i);it++)
      {                                         // incidents for this period
        int rp=realization_period(i,it);
        //int ri=realization_incident(i,it);
        //Y(i,it)-=value(fish_pars(1,i))*
         // sin(tpi*(true_month(rr,rp)/12.-value(fish_pars(2,i))));
      }
    }
  }
}
dvar_vector dvar_fish_stock_history::get_seasonal_catchability_effect(int i)
{
  const MY_DOUBLE_TYPE tpi=2.e0*3.14159e0;
  dvar_vector tmp(1,12); 
  for (int j=1;j<=12;j++)
  {
    tmp(j)=fish_pars(1,i)*sin(tpi*j/12.-fish_pars(2,i));
  }
  return tmp;
}

void dvar_fish_stock_history::vimplicit_seasonal_catchability_calc
  (dvar_matrix& Y)
{
  int i,it;
  const MY_DOUBLE_TYPE tpi=2.e0*3.14159e0;
  for (i=1;i<=num_fisheries;i++)
  {
    int rr=realization_region(i,1);
    if (fish_flags(i,27)==1)
    {
      dvar_vector  sce=get_seasonal_catchability_effect(i);
      for (it=1;it<=num_real_fish_times(i);it++)
      {                                         // incidents for this period
        int rp=realization_period(i,it);
        //Y(i,it)-=sce(true_month(rr,rp));
      }
    }
  }
}

MY_DOUBLE_TYPE dvar_fish_stock_history::
  grouped_implicit_catchability_deviations_calc_part1(int n,dvector& x,
    ivector& iswitch,imatrix& indx,MY_DOUBLE_TYPE c)
{
  MY_DOUBLE_TYPE f;
  x.initialize();
  //get_initial_x(Y,x);
  dvector g(1,n);
  banded_symmetric_dmatrix H(1,n,2);
  if (pfmin1)
  {
    if (pfmin1->n != n)
    {
# if !defined(__BORLANDC__)
      delete pfmin1;
# endif
      pfmin1=0;
    }
  }
    
  if (pfmin1==0)
  {
    pfmin1 = new fmmt1(n,3);
  }
  else
  {
    pfmin1->set_defaults();
  }
  pfmin1->use_control_c=0;
 
  pfmin1->noprintx=1;
  pfmin1->iprint=0;
  pfmin1->crit=1.e-3;
  //pfmin1->crit=.100;
  //dmatrix Y=value(effort_dev_coffs);
  dmatrix Y=value(implicit_catchability);
  int i27=sum(column(fish_flags,27));
  if (i27)
  {
    cimplicit_seasonal_catchability_calc(Y);
  }
  dmatrix resids(1,num_fisheries,1,num_real_fish_times);
  resids.initialize();
  dmatrix wts(1,num_fisheries,1,num_real_fish_times);
  wts.initialize();
  MY_DOUBLE_TYPE normg;
  // ****************************************************************
  // ****************************************************************
  //   !!!!!!!!!!!!!!!!  TEST  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // ****************************************************************
  //
  // least squares
  //
  get_f(f,Y,x,g,H,1,0,resids,wts,0,iswitch,indx);
  cout <<"mfclimplicit.cpp " << setprecision(15) << f << endl;
  dvector h=safe_solve(H,g,0.0); 
  x-= h;
  get_f(f,Y,x,g,H,1,0,resids,wts,0,iswitch,indx);
  //cout <<"mfclimplicit.cpp " << setprecision(15) << f << endl;
  MY_DOUBLE_TYPE ng=norm(g);
  MY_DOUBLE_TYPE nh=norm(h);
  /*
  cout << " g*h/(norm(g)*norm(h)) "  << endl;
  cout <<"mfclimplicit.cpp " <<  g*h/(ng*nh)  << endl;
  cout << "   norm(g) = " << ng << endl;
  */
  //double c=0.05;
  // start with truncated newton
  while (pfmin1->ireturn>=0)
  {
    pfmin1->fmin(f,x,g);
    if (pfmin1->ireturn>0)
    {
      get_f(f,Y,x,g,H,0,1,resids,wts,0,iswitch,indx,c);
      //cout <<"mfclimplicit.cpp " << setprecision(15) << f << endl;
    }
  }
  //then trust regionm to get good convergence 
  int trust_switch=1;
  if (trust_switch)
  {
  // ****************************************************************
  //   !!!!!!!!!!!!!!!!  TEST TRUST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // ****************************************************************
  //
  //  get g and  H
    MY_DOUBLE_TYPE ftmp;
    MY_DOUBLE_TYPE Delta=5.0;
    MY_DOUBLE_TYPE mf;
    dvector gbest=g;
    ng=norm(gbest);
    MY_DOUBLE_TYPE ngold=ng;
    get_f(ftmp,Y,x,g,H,1,1,resids,wts,0,iswitch,indx,c);
    MY_DOUBLE_TYPE f0=ftmp;
    int icount1=0;
    int maxcount1=20;
    int first_try=5;
    MY_DOUBLE_TYPE ns=0.0;
    do
    {
      
      MY_DOUBLE_TYPE change;
  //
  //  solve secular equation 
      dvector s=solve_secular_equation(&mf,g,H,Delta,&first_try,ns);
  //  compare f change to model
      dvector xtry=x+s;
      get_f(ftmp,Y,xtry,g,H,1,1,resids,wts,0,iswitch,indx,c);

      //cout << f0-ftmp << " " << -mf << " "
      //     << (ftmp-f0)/(1.e-100+mf);
      if (f0<=ftmp)
      {
       // cout << " " << Delta << " " << norm(s); 
      }
      //cout <<"mfclimplicit.cpp " << endl; 
      Delta=get_new_delta(ftmp,f0,mf,Delta,ns);
      int backswitch=1;
      if (ftmp<f0)
      {
        backswitch=0;
        // keep new g and H
        f0=ftmp;
        x=xtry;
        gbest=g;
        ngold=ng;
        ng=norm(gbest);
      }
      else if (ftmp==f0)
      {
        MY_DOUBLE_TYPE ng1=norm(g);
        if (ng1<ng) 
        {
          // keep new g and H
          f0=ftmp;
          x=xtry;
          gbest=g;
          ngold=ng;
          ng=ng1;
        }
      }
      else
      {
        // get back old g and H
        get_f(ftmp,Y,x,g,H,1,1,resids,wts,0,iswitch,indx,c);
      }
      icount1++;
       
      if (ng<1.e-10)
      {
        break;
      }
      else if (ng<1.e-8 && icount1 > 4)
      {
        break;
      }
      else if (ng<1.e-7 && icount1 > 7)
      {
        break;
      }
      else if (ng<1.e-6 && icount1 > 10)
      {
        break;
      }
      else if (icount1 >= maxcount1)
      {
        break;
      }
    }
    while(1);
  }
  //
  // robust c=0.001
  //
  /*
  while (pfmin1->ireturn>=0)
  {
    pfmin1->fmin(f,x,g);
    if (pfmin1->ireturn>0)
    {
      get_f(f,Y,x,g,H,0,1,resids,wts,0,iswitch,indx,c);
      //cout <<"mfclimplicit.cpp " << setprecision(15) << f << endl;
    }
  }
  */
  /*   for testing that Hessian etc is correct
  {
    MY_DOUBLE_TYPE f1,f2;
    int mmin=x.indexmin();
    int mmax=x.indexmax();
    dvector xt(mmin,mmax);
    dvector tg1(mmin,mmax);
    dvector tg2(mmin,mmax);
    dmatrix TH(mmin,mmax);
    x=xt;
    MY_DOUBLE_TYPE delta=1.e-5;
    for (int i=1;i<=n;i++)
    {
      xt(i)-=delta;
      get_f(f1,Y,xt,tg1,H,0,1,resids,wts,0,iswitch,indx,c);
      xt(i)+=2.0*delta;
      get_f(f2,Y,xt,tg2,H,0,1,resids,wts,0,iswitch,indx,c);
      xt(i)-=delta;
      TH(i)=(tg2-tg1)/(2.0*delta);
    }
    get_f(f2,Y,x,g,H,1,1,resids,wts,0,iswitch,indx,c);
    cout <<"mfclimplicit.cpp " << norm2(H-TH) << endl;
  }
  */
      
 /*
  {
    get_f(f,Y,x,g,H,1,1,resids,wts,0,iswitch,indx,c);
    MY_DOUBLE_TYPE oldf=f;
    //for (int ik=1;ik<=15;ik++)
    int ic1=0;
    int badflag=0;
    do
    {
      dvector h=safe_solve(H,g,0.0); 
      MY_DOUBLE_TYPE nh=norm(h);
      MY_DOUBLE_TYPE pf1=g*h/nh;
      if (pf1<0)
      {
        cerr << "infeasible direction" << endl;
        ad_exit(1);
      }
      int ic2=0;
      do
      {
        
        dvector xtry=x-h; 
        get_f(f,Y,xtry,g,H,1,1,resids,wts,0,iswitch,indx,c);
        MY_DOUBLE_TYPE pf2=g*h/nh;
  
        dvector tmp=get_coffs2(oldf,f,pf1,-nh);
        MY_DOUBLE_TYPE tmin=get_min(tmp);
        dvector x2=x+tmin/nh*h;
        MY_DOUBLE_TYPE f2;
        get_f(f2,Y,x2,g,H,1,1,resids,wts,0,iswitch,indx,c);
  
        dvector tmp2=get_coffs(oldf,f,pf1,pf2,-nh);
        MY_DOUBLE_TYPE tmin2=get_min(tmp2);
        dvector x3=x+tmin2/nh*h;
        MY_DOUBLE_TYPE f3;
        get_f(f3,Y,x3,g,H,1,1,resids,wts,0,iswitch,indx,c);
  
        MY_DOUBLE_TYPE pf3=g*h/nh;
        // cout <<"mfclimplicit.cpp " << pf3/norm(g) << endl;
        dvector tmp3=get_coffs(oldf,f3,pf1,pf3,tmin2/nh);
        MY_DOUBLE_TYPE tmin3=get_min(tmp3);
  
        dvector tmp4=get_coffs2(oldf,f3,pf1,tmin2/nh);
        MY_DOUBLE_TYPE tmin4=get_min(tmp4);
  
        if (min(f,f2,f3)<oldf)
        {
          switch (minindex(f,f2,f3))
          {
          case 1:
            x=xtry; 
            oldf=f;
            break;
          case 2:
            x=x2; 
            oldf=f2;
            break;
          case 3:
            x=x3; 
            oldf=f3;
            break;
          }
          break;
        }
        else
        {
          MY_DOUBLE_TYPE f4;
          dvector x4=x+tmin3*fabs(tmin2)/nh*h;
          get_f(f4,Y,x4,g,H,1,1,resids,wts,0,iswitch,indx,c);
  
          MY_DOUBLE_TYPE f5;
          dvector x5=x+tmin4*fabs(tmin2)/nh*h;
          get_f(f5,Y,x5,g,H,1,1,resids,wts,0,iswitch,indx,c);
  
  
          dvector tmp5=get_coffs2(oldf,f4,pf1,tmin3*fabs(tmin2)/nh);
          MY_DOUBLE_TYPE tmin5=get_min(tmp5);
          dvector x6=x+tmin5*fabs(tmin3*tmin2)/nh*h;
          MY_DOUBLE_TYPE f6;
          get_f(f6,Y,x6,g,H,1,1,resids,wts,0,iswitch,indx,c);
  
          dvector tmp6=get_coffs2(oldf,f6,pf1,tmin5*fabs(tmin3*tmin2)/nh);
          MY_DOUBLE_TYPE tmin6=get_min(tmp6);
          MY_DOUBLE_TYPE f7;
          dvector x7=x+tmin6*fabs(tmin5*tmin3*tmin2)/nh*h;
          get_f(f7,Y,x7,g,H,1,1,resids,wts,0,iswitch,indx,c);
  
          if (f5<oldf)
          {
            x=x5; 
            oldf=f5;
            break;
          }
          else if (f4<oldf)
          {
            x=x4; 
            oldf=f4;
            break;
          }
          else
          {
            //cerr << "need better minimzier" << endl;
            badflag=1;
            break;
          }
        }
      }
      while(1);
      get_f(f,Y,x,g,H,1,1,resids,wts,0,iswitch,indx,c);
      normg=norm(g);
      //cout << "norm(g) = " << normg << endl;
      ic1++;
      if (badflag)
        break;
      if (ic2>=10) break;
      if (normg<1.e-9) break;
      if (normg<1.e-8 && ic1>10) break;
      if (normg<1.e-6 && ic1>20) break;
      if (f!=oldf)
      {
         cerr << "error" << endl;
      }
    }
    while (ic1<100);
    //cout << " ic1 = " << ic1 << "  norm(g) = " << normg << endl;
  }
 */
  // ****************************************************************
  // ****************************************************************
  // ****************************************************************
  // ****************************************************************
 
 /*
  get_f(f,Y,x,g,H,1,0,resids,wts,0,iswitch,indx);
  x-=safe_solve(H,g,0.0); 
  get_f(f,Y,x,g,H,1,0,resids,wts,0,iswitch,indx);
  for (int i=1;i<=num_fisheries;i++)
  {
    for (int j=1;j<=num_real_fish_times(i);j++)
    {
      wts(i,j)=exp(-square(resids(i,j)))+.02/(1+square(resids(i,j)));
      //wts(i,j)=exp(-square(resids(i,j)));
      //wts(i,j)=1.0/(1+square(resids(i,j)));
    }
  }
  get_f(f,Y,x,g,H,1,0,resids,wts,1,iswitch,indx);
  x-=safe_solve(H,g,0.0); 
  for (int ic=1;ic<=50;ic++) 
  {
    for (i=1;i<=num_fisheries;i++)
    {
      for (int j=1;j<=num_real_fish_times(i);j++)
      {
        wts(i,j)=1.0/(1+square(resids(i,j)));
      }
    }
    get_f(f,Y,x,g,H,1,0,resids,wts,1,iswitch,indx);
    x-=safe_solve(H,g,0.0); 
  }
  
 */

  int icount=0;
 
 /* 
  dvector oldx(1,n);
  MY_DOUBLE_TYPE oldf;
  MY_DOUBLE_TYPE imult=5.0;
  dvector h(1,n);
  int aflag=0;
  int ic1=0;
  do
  {
    get_f(f,Y,x,g,H,1,1,resids,wts,0);
    if (icount>0 && oldf<f && aflag==0)
    {
      //cout <<"mfclimplicit.cpp " << f-oldf << endl;
      ic1++;
      x=oldx;
      if (ic1>5)
      {
        aflag=1;
      }
      else
      {
#if !defined(NO_MY_DOUBLE_TYPE)
        x-=(1.0/(imult+1.0L))*h; 
#else
        x-=(1.0/(imult+1.0))*h; 
#endif
        imult*=10.0;
      }
    }
    else
    {
      oldf=f; 
      oldx=x;
      MY_DOUBLE_TYPE id;
      aflag=0;
      id=ic1;
      ic1=0;
      h=safe_solve(H,g,id); 
      //cout <<"mfclimplicit.cpp " << h*g/(norm(h)*norm(g)) << endl; 
      x-=h;
      ng=norm(g);
      //cout << setscientific() << setprecision(5) << ng << " " << max(g) << endl;
      imult=1.0;
    }
    if (icount>30 && ng<1.e-8)
      break;
    if (++icount>50)
      break;
  }
  while (ng>1.e-13);

  
  if (ng>1.e-8)
  {
    cout << "Initial  " << setscientific() << setprecision(5) 
         << ng << " " << max(g) << endl;
  }
*/ 

  if (ng>1.e-4)
  {
    pfmin1->ireturn=0;
    while (pfmin1->ireturn>=0)
    {
      pfmin1->fmin(f,x,g);
      if (pfmin1->ireturn>0)
      {
         get_f(f,Y,x,g,H,0,1,resids,wts,0,iswitch,indx);
      }
    }
    icount=0;
    do 
    {
      get_f(f,Y,x,g,H,1,1,resids,wts,0,iswitch,indx);
      x-=safe_solve(H,g,0.0); 
      ng=norm(g);
      //cout <<"mfclimplicit.cpp " << setscientific() << setprecision(5) << ng 
      //       << " " << max(g) << endl;
      if (++icount>20)
        break;
    }
    while (ng>1.e-10);

    get_f(f,Y,x,g,H,1,1,resids,wts,0,iswitch,indx);
  // one more time with quad double
  }

  if (ng>1.e-11)
  {
    for (int ix=1;ix<=7;ix++)
    {
      //cout <<"mfclimplicit.cpp " << ng << endl;
      dvector h(1,n);
      test_fpu();
#  if ( defined(__MSVC32__) &&  __MSVC32__ >=8)
      dd_newton_raphson(n,&(g(1)),&H(1,1),&H(2,1),&(h(1)));
      //cerr << "not implemented" << endl;
     // ad_exit(1);
#  else
      dd_newton_raphson(n,&(g(1)),&H(1,1),&H(2,1),&(h(1)));
#  endif
      x-=h;
      get_f(f,Y,x,g,H,1,1,resids,wts,0,iswitch,indx,c);
      ng=norm(g);
      if (ng<1.e-12) break;
    }
  }
  
  //if (ng>1.e-8)
  {
    //cout << "Final  " << setscientific() << setprecision(5) 
    //     << normg << " " << max(g) << endl;
  }
  // put these here so that when we start without catch conditioning
  // the effort_dev_coffs will be correct
  ivector ff10=column(fish_flags,10);
  for (int j=1;j<=num_fisheries;j++)
  {
    if (ff10(j))
    {
      effort_dev_coffs(j)=resids(j);
    }
  }
  return f;
}

dvariable dvar_fish_stock_history::
  grouped_implicit_catchability_deviations_calc_part2(dvector& x,
  ivector& iswitch,imatrix& indx,MY_DOUBLE_TYPE _c)
{
  dvar_matrix vY(1,num_fisheries,1,num_fish_times);
  vY=implicit_catchability;
  int i27=sum(column(fish_flags,27));
  if (i27)
  {
    vimplicit_seasonal_catchability_calc(vY);
  }
  dvariable vf=get_vf_normal(vY,x,iswitch,indx,_c);
  return vf;
}

// routines to get lower bound on lambda TRM page 92

MY_DOUBLE_TYPE lb1(const banded_symmetric_dmatrix& m)
{
  MY_DOUBLE_TYPE tmp=min(m.get_dmatrix()(0));
  return -tmp;
}
  
MY_DOUBLE_TYPE lb2(const banded_symmetric_dmatrix& m)
{
  int mmin=m.indexmin();
  int mmax=m.indexmax();
  dvector tmp(mmin,mmax);
  tmp.initialize();
  int bw=m.bandwidth();
  for (int i=mmin;i<=mmax;i++)
  {
    for (int j=-bw+1;j<=bw-1;j++)
    {
      if (j==0) 
        tmp(i)+=m(i,i);
      else
      {
        if (i>i+j)
        {
          if (i+j>=mmin) 
            tmp(i)+=fabs(m(i,i+j));
        }
        else 
        {
          if (i+j<=mmax)
            tmp(i)+=fabs(m(i+j,i));
        }
      }
    }
  }
  return max(tmp);
}
  
MY_DOUBLE_TYPE lb2u(const banded_symmetric_dmatrix& m)
{
  int mmin=m.indexmin();
  int mmax=m.indexmax();
  dvector tmp(mmin,mmax);
  tmp.initialize();
  int bw=m.bandwidth();
  for (int i=mmin;i<=mmax;i++)
  {
    for (int j=-bw+1;j<=bw-1;j++)
    {
      if (j==0) 
        tmp(i)-=m(i,i);
      else
      {
        if (i>i+j)
        {
          if (i+j>=mmin) 
            tmp(i)+=fabs(m(i,i+j));
        }
        else 
        {
          if (i+j<=mmax)
            tmp(i)+=fabs(m(i+j,i));
        }
      }
    }
  }
  return max(tmp);
}
  
MY_DOUBLE_TYPE lb4(const banded_symmetric_dmatrix& m)
{
  // infinity norm
  MY_DOUBLE_TYPE fm=max(m.get_dmatrix()(0));
  int bw=m.bandwidth();
  for (int i=0;i<=bw-1;i++)
  {
    fm=max(fm,max(m.get_dmatrix()(i)));
  }
  return fm;
}
  
  
MY_DOUBLE_TYPE get_init_lambda_l(const banded_symmetric_dmatrix& m,dvector& g,
  MY_DOUBLE_TYPE delta)
{
  dvector tmp(1,3);
  tmp(1)=lb2(m);
  tmp(2)=norm(m);
  tmp(3)=lb4(m);
  dvector tmp2(1,3);
  tmp2(3)=norm(g)/delta-min(tmp);
  tmp2(2)=lb1(m);
  tmp2(1)=0;
  return max(tmp2);
}


MY_DOUBLE_TYPE get_init_lambda_u(const banded_symmetric_dmatrix& m,dvector& g,
  MY_DOUBLE_TYPE delta)
{
  dvector tmp(1,3);
  tmp(1)=lb2u(m);
  tmp(2)=norm(m);
  tmp(3)=lb4(m);
  dvector tmp2(1,2);
  tmp2(2)=norm(g)/delta+min(tmp);
  tmp2(1)=0;
  return max(tmp2);
}

dvariable dvar_fish_stock_history::get_vf_normal(dvar_matrix& Y,dvector& x,
  ivector& iswitch,imatrix& indx,MY_DOUBLE_TYPE _c)
{
  MY_DOUBLE_TYPE c=_c/3.0;
  //ofstream ofs("vimp");
  dvariable vf;
  vf=0.0;
  
  //  s_eps is the std dev of the effort devs
  MY_DOUBLE_TYPE s_eps=0.3;
  MY_DOUBLE_TYPE pen1=1.0/square(s_eps);
  // s_eta is the std dev of the time series changes in catchability
  MY_DOUBLE_TYPE s_eta=0.1;
  MY_DOUBLE_TYPE thpen1=1.0/square(3.0*s_eps);
  MY_DOUBLE_TYPE sqpen1=square(pen1);
  MY_DOUBLE_TYPE sqthpen1=square(thpen1);
  MY_DOUBLE_TYPE pen2=1.0/(2.0*square(s_eta));
  int n=x.indexmax();
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
  ii=1;
  for (j=1;j<=num_fisheries;j++)
  {
    if (iswitch(grouping(j)))
    {
      for (i=1;i<=num_real_fish_times(j);i++)
      {
        // only do this if we have fishing effort
#if !defined(NO_MY_DOUBLE_TYPE)
        if (true_effort_by_fishery(j,i)>-0.5L)
#else
        if (true_effort_by_fishery(j,i)>-0.5)
#endif
        {
          dvariable r=Y(j,i)-grouped_catchability_coffs(grouping(j),
            gfish_ptr(j,i));
          dvariable r2=square(r);
          dvariable u1=exp(-0.5*pen1*r2);
          dvariable v1=exp(-0.5*thpen1*r2);
          dvariable u=1.e-30+u1+c*v1;
          /*********************************************************/
          vf-=log(u);
        }
      }
    }
  }
  //cout << "AA " << vf << endl;
  //exit(1);

  const ivector& xcatflags=(column(fish_flags,15));  
  ivector& catflags=(ivector&) xcatflags;  
  ii=1;
  for (j=1;j<=ngroups;j++)
  {
    if (iswitch(j))
    {
      ii++;
      for (i=2;i<=num_real_grouped_fish_times(j);i++)
      {
        MY_DOUBLE_TYPE pp2;
        int ff15=catflags(gfish_index(j,1));
        if (ff15)
        {
          pp2=ff15*grouped_between_times(j,i);
        }
        else
        {
          pp2=pen2*grouped_between_times(j,i);
        }
        
        if (pp2==0)
        {
          cerr << "this can't happen" << endl;
          ad_exit(1);
        }
        dvariable r=grouped_catchability_coffs(j,i)
          -grouped_catchability_coffs(j,i-1);
        vf+=0.5*pp2*square(r);
      }
    }
  }
  //cout << "BB " << vf << endl;
  return vf;
}

#if defined(USE_ADPVM)
MY_DOUBLE_TYPE dvar_fish_stock_history::
  grouped_implicit_catchability_deviations_calc_part1_pvm(int n,dvector& x,
    ivector& iswitch,imatrix& indx,MY_DOUBLE_TYPE c)
{
  MY_DOUBLE_TYPE f;
  x.initialize();
  //get_initial_x(Y,x);
  dvector g(1,n);
  banded_symmetric_dmatrix H(1,n,2);
  if (pfmin1)
  {
    if (pfmin1->n != n)
    {
# if !defined(__BORLANDC__)
      delete pfmin1;
# endif
      pfmin1=0;
    }
  }
    
  if (pfmin1==0)
  {
    pfmin1 = new fmmt1(n,3);
  }
  else
  {
    pfmin1->set_defaults();
  }
  pfmin1->use_control_c=0;
 
  pfmin1->noprintx=1;
  pfmin1->iprint=0;
  pfmin1->crit=1.e-3;
  //pfmin1->crit=.100;
  //dmatrix Y=value(effort_dev_coffs);
  dmatrix Y=get_dmatrix_from_master();
  dmatrix resids(1,num_fisheries,1,num_real_fish_times);
  resids.initialize();
  dmatrix wts(1,num_fisheries,1,num_real_fish_times);
  wts.initialize();
  MY_DOUBLE_TYPE normg;
  // ****************************************************************
  // ****************************************************************
  //   !!!!!!!!!!!!!!!!  TEST  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // ****************************************************************
  //
  // least squares
  //
  get_f(f,Y,x,g,H,1,0,resids,wts,0,iswitch,indx);
  cout <<"mfclimplicit.cpp " << setprecision(15) << f << endl;
  dvector h=safe_solve(H,g,0.0); 
  x-= h;
  get_f(f,Y,x,g,H,1,0,resids,wts,0,iswitch,indx);
  //cout <<"mfclimplicit.cpp " << setprecision(15) << f << endl;
  MY_DOUBLE_TYPE ng=norm(g);
  MY_DOUBLE_TYPE nh=norm(h);
  /*
  cout << " g*h/(norm(g)*norm(h)) "  << endl;
  cout <<"mfclimplicit.cpp " <<  g*h/(ng*nh)  << endl;
  cout << "   norm(g) = " << ng << endl;
  */
  //double c=0.05;
  // start with truncated newton
  while (pfmin1->ireturn>=0)
  {
    pfmin1->fmin(f,x,g);
    if (pfmin1->ireturn>0)
    {
      get_f(f,Y,x,g,H,0,1,resids,wts,0,iswitch,indx,c);
      //cout <<"mfclimplicit.cpp " << setprecision(15) << f << endl;
    }
  }
  //then trust regionm to get good convergence 
  int trust_switch=1;
  if (trust_switch)
  {
  // ****************************************************************
  //   !!!!!!!!!!!!!!!!  TEST TRUST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // ****************************************************************
  //
  //  get g and  H
    MY_DOUBLE_TYPE ftmp;
    MY_DOUBLE_TYPE Delta=5.0;
    MY_DOUBLE_TYPE mf;
    dvector gbest=g;
    ng=norm(gbest);
    MY_DOUBLE_TYPE ngold=ng;
    get_f(ftmp,Y,x,g,H,1,1,resids,wts,0,iswitch,indx,c);
    MY_DOUBLE_TYPE f0=ftmp;
    int icount1=0;
    int maxcount1=20;
    int first_try=5;
    MY_DOUBLE_TYPE ns=0.0;
    do
    {
      
      MY_DOUBLE_TYPE change;
  //
  //  solve secular equation 
      dvector s=solve_secular_equation(&mf,g,H,Delta,&first_try,ns);
  //  compare f change to model
      dvector xtry=x+s;
      get_f(ftmp,Y,xtry,g,H,1,1,resids,wts,0,iswitch,indx,c);

      //cout << f0-ftmp << " " << -mf << " "
      //     << (ftmp-f0)/(1.e-100+mf);
      if (f0<=ftmp)
      {
       // cout << " " << Delta << " " << norm(s); 
      }
      //cout <<"mfclimplicit.cpp " << endl; 
      Delta=get_new_delta(ftmp,f0,mf,Delta,ns);
      int backswitch=1;
      if (ftmp<f0)
      {
        backswitch=0;
        // keep new g and H
        f0=ftmp;
        x=xtry;
        gbest=g;
        ngold=ng;
        ng=norm(gbest);
      }
      else if (ftmp==f0)
      {
        MY_DOUBLE_TYPE ng1=norm(g);
        if (ng1<ng) 
        {
          // keep new g and H
          f0=ftmp;
          x=xtry;
          gbest=g;
          ngold=ng;
          ng=ng1;
        }
      }
      else
      {
        // get back old g and H
        get_f(ftmp,Y,x,g,H,1,1,resids,wts,0,iswitch,indx,c);
      }
      icount1++;
       
      if (ng<1.e-10)
      {
        break;
      }
      else if (ng<1.e-8 && icount1 > 4)
      {
        break;
      }
      else if (ng<1.e-7 && icount1 > 7)
      {
        break;
      }
      else if (ng<1.e-6 && icount1 > 10)
      {
        break;
      }
      else if (icount1 >= maxcount1)
      {
        break;
      }
    }
    while(1);
  }
  //
  // robust c=0.001
  //
  /*
  while (pfmin1->ireturn>=0)
  {
    pfmin1->fmin(f,x,g);
    if (pfmin1->ireturn>0)
    {
      get_f(f,Y,x,g,H,0,1,resids,wts,0,iswitch,indx,c);
      //cout <<"mfclimplicit.cpp " << setprecision(15) << f << endl;
    }
  }
  */
  /*   for testing that Hessian etc is correct
  {
    MY_DOUBLE_TYPE f1,f2;
    int mmin=x.indexmin();
    int mmax=x.indexmax();
    dvector xt(mmin,mmax);
    dvector tg1(mmin,mmax);
    dvector tg2(mmin,mmax);
    dmatrix TH(mmin,mmax);
    x=xt;
    MY_DOUBLE_TYPE delta=1.e-5;
    for (int i=1;i<=n;i++)
    {
      xt(i)-=delta;
      get_f(f1,Y,xt,tg1,H,0,1,resids,wts,0,iswitch,indx,c);
      xt(i)+=2.0*delta;
      get_f(f2,Y,xt,tg2,H,0,1,resids,wts,0,iswitch,indx,c);
      xt(i)-=delta;
      TH(i)=(tg2-tg1)/(2.0*delta);
    }
    get_f(f2,Y,x,g,H,1,1,resids,wts,0,iswitch,indx,c);
    cout <<"mfclimplicit.cpp " << norm2(H-TH) << endl;
  }
  */
      
 /*
  {
    get_f(f,Y,x,g,H,1,1,resids,wts,0,iswitch,indx,c);
    MY_DOUBLE_TYPE oldf=f;
    //for (int ik=1;ik<=15;ik++)
    int ic1=0;
    int badflag=0;
    do
    {
      dvector h=safe_solve(H,g,0.0); 
      MY_DOUBLE_TYPE nh=norm(h);
      MY_DOUBLE_TYPE pf1=g*h/nh;
      if (pf1<0)
      {
        cerr << "infeasible direction" << endl;
        ad_exit(1);
      }
      int ic2=0;
      do
      {
        
        dvector xtry=x-h; 
        get_f(f,Y,xtry,g,H,1,1,resids,wts,0,iswitch,indx,c);
        MY_DOUBLE_TYPE pf2=g*h/nh;
  
        dvector tmp=get_coffs2(oldf,f,pf1,-nh);
        MY_DOUBLE_TYPE tmin=get_min(tmp);
        dvector x2=x+tmin/nh*h;
        MY_DOUBLE_TYPE f2;
        get_f(f2,Y,x2,g,H,1,1,resids,wts,0,iswitch,indx,c);
  
        dvector tmp2=get_coffs(oldf,f,pf1,pf2,-nh);
        MY_DOUBLE_TYPE tmin2=get_min(tmp2);
        dvector x3=x+tmin2/nh*h;
        MY_DOUBLE_TYPE f3;
        get_f(f3,Y,x3,g,H,1,1,resids,wts,0,iswitch,indx,c);
  
        MY_DOUBLE_TYPE pf3=g*h/nh;
        // cout <<"mfclimplicit.cpp " << pf3/norm(g) << endl;
        dvector tmp3=get_coffs(oldf,f3,pf1,pf3,tmin2/nh);
        MY_DOUBLE_TYPE tmin3=get_min(tmp3);
  
        dvector tmp4=get_coffs2(oldf,f3,pf1,tmin2/nh);
        MY_DOUBLE_TYPE tmin4=get_min(tmp4);
  
        if (min(f,f2,f3)<oldf)
        {
          switch (minindex(f,f2,f3))
          {
          case 1:
            x=xtry; 
            oldf=f;
            break;
          case 2:
            x=x2; 
            oldf=f2;
            break;
          case 3:
            x=x3; 
            oldf=f3;
            break;
          }
          break;
        }
        else
        {
          MY_DOUBLE_TYPE f4;
          dvector x4=x+tmin3*fabs(tmin2)/nh*h;
          get_f(f4,Y,x4,g,H,1,1,resids,wts,0,iswitch,indx,c);
  
          MY_DOUBLE_TYPE f5;
          dvector x5=x+tmin4*fabs(tmin2)/nh*h;
          get_f(f5,Y,x5,g,H,1,1,resids,wts,0,iswitch,indx,c);
  
  
          dvector tmp5=get_coffs2(oldf,f4,pf1,tmin3*fabs(tmin2)/nh);
          MY_DOUBLE_TYPE tmin5=get_min(tmp5);
          dvector x6=x+tmin5*fabs(tmin3*tmin2)/nh*h;
          MY_DOUBLE_TYPE f6;
          get_f(f6,Y,x6,g,H,1,1,resids,wts,0,iswitch,indx,c);
  
          dvector tmp6=get_coffs2(oldf,f6,pf1,tmin5*fabs(tmin3*tmin2)/nh);
          MY_DOUBLE_TYPE tmin6=get_min(tmp6);
          MY_DOUBLE_TYPE f7;
          dvector x7=x+tmin6*fabs(tmin5*tmin3*tmin2)/nh*h;
          get_f(f7,Y,x7,g,H,1,1,resids,wts,0,iswitch,indx,c);
  
          if (f5<oldf)
          {
            x=x5; 
            oldf=f5;
            break;
          }
          else if (f4<oldf)
          {
            x=x4; 
            oldf=f4;
            break;
          }
          else
          {
            //cerr << "need better minimzier" << endl;
            badflag=1;
            break;
          }
        }
      }
      while(1);
      get_f(f,Y,x,g,H,1,1,resids,wts,0,iswitch,indx,c);
      normg=norm(g);
      //cout << "norm(g) = " << normg << endl;
      ic1++;
      if (badflag)
        break;
      if (ic2>=10) break;
      if (normg<1.e-9) break;
      if (normg<1.e-8 && ic1>10) break;
      if (normg<1.e-6 && ic1>20) break;
      if (f!=oldf)
      {
         cerr << "error" << endl;
      }
    }
    while (ic1<100);
    //cout << " ic1 = " << ic1 << "  norm(g) = " << normg << endl;
  }
 */
  // ****************************************************************
  // ****************************************************************
  // ****************************************************************
  // ****************************************************************
 
 /*
  get_f(f,Y,x,g,H,1,0,resids,wts,0,iswitch,indx);
  x-=safe_solve(H,g,0.0); 
  get_f(f,Y,x,g,H,1,0,resids,wts,0,iswitch,indx);
  for (int i=1;i<=num_fisheries;i++)
  {
    for (int j=1;j<=num_real_fish_times(i);j++)
    {
      wts(i,j)=exp(-square(resids(i,j)))+.02/(1+square(resids(i,j)));
      //wts(i,j)=exp(-square(resids(i,j)));
      //wts(i,j)=1.0/(1+square(resids(i,j)));
    }
  }
  get_f(f,Y,x,g,H,1,0,resids,wts,1,iswitch,indx);
  x-=safe_solve(H,g,0.0); 
  for (int ic=1;ic<=50;ic++) 
  {
    for (i=1;i<=num_fisheries;i++)
    {
      for (int j=1;j<=num_real_fish_times(i);j++)
      {
        wts(i,j)=1.0/(1+square(resids(i,j)));
      }
    }
    get_f(f,Y,x,g,H,1,0,resids,wts,1,iswitch,indx);
    x-=safe_solve(H,g,0.0); 
  }
  
 */

  int icount=0;
 
 /* 
  dvector oldx(1,n);
  MY_DOUBLE_TYPE oldf;
  MY_DOUBLE_TYPE imult=5.0;
  dvector h(1,n);
  int aflag=0;
  int ic1=0;
  do
  {
    get_f(f,Y,x,g,H,1,1,resids,wts,0);
    if (icount>0 && oldf<f && aflag==0)
    {
      //cout <<"mfclimplicit.cpp " << f-oldf << endl;
      ic1++;
      x=oldx;
      if (ic1>5)
      {
        aflag=1;
      }
      else
      {
#if !defined(NO_MY_DOUBLE_TYPE)
        x-=(1.0/(imult+1.0L))*h; 
#else
        x-=(1.0/(imult+1.0))*h; 
#endif
        imult*=10.0;
      }
    }
    else
    {
      oldf=f; 
      oldx=x;
      MY_DOUBLE_TYPE id;
      aflag=0;
      id=ic1;
      ic1=0;
      h=safe_solve(H,g,id); 
      //cout <<"mfclimplicit.cpp " << h*g/(norm(h)*norm(g)) << endl; 
      x-=h;
      ng=norm(g);
      //cout << setscientific() << setprecision(5) << ng << " " << max(g) << endl;
      imult=1.0;
    }
    if (icount>30 && ng<1.e-8)
      break;
    if (++icount>50)
      break;
  }
  while (ng>1.e-13);

  
  if (ng>1.e-8)
  {
    cout << "Initial  " << setscientific() << setprecision(5) 
         << ng << " " << max(g) << endl;
  }
*/ 

  if (ng>1.e-4)
  {
    pfmin1->ireturn=0;
    while (pfmin1->ireturn>=0)
    {
      pfmin1->fmin(f,x,g);
      if (pfmin1->ireturn>0)
      {
         get_f(f,Y,x,g,H,0,1,resids,wts,0,iswitch,indx);
      }
    }
    icount=0;
    do 
    {
      get_f(f,Y,x,g,H,1,1,resids,wts,0,iswitch,indx);
      x-=safe_solve(H,g,0.0); 
      ng=norm(g);
      //cout <<"mfclimplicit.cpp " << setscientific() << setprecision(5) << ng 
      //       << " " << max(g) << endl;
      if (++icount>20)
        break;
    }
    while (ng>1.e-10);

    get_f(f,Y,x,g,H,1,1,resids,wts,0,iswitch,indx);
  // one more time with quad double
  }

  if (ng>1.e-11)
  {
    for (int ix=1;ix<=5;ix++)
    {
      //cout <<"mfclimplicit.cpp " << ng << endl;
      dvector h(1,n);
      test_fpu();
#  if ( defined(__MSVC32__) &&  __MSVC32__ >=8)
      dd_newton_raphson(n,&(g(1)),&H(1,1),&H(2,1),&(h(1)));
      //cerr << "not implemented" << endl;
     // ad_exit(1);
#  else
      dd_newton_raphson(n,&(g(1)),&H(1,1),&H(2,1),&(h(1)));
#  endif
      x-=h;
      get_f(f,Y,x,g,H,1,1,resids,wts,0,iswitch,indx,c);
      ng=norm(g);
      if (ng<1.e-12) break;
    }
  }
  
  //if (ng>1.e-8)
  {
    //cout << "Final  " << setscientific() << setprecision(5) 
    //     << normg << " " << max(g) << endl;
  }
 
  return f;
}
void dvar_fish_stock_history::
  grouped_implicit_catchability_deviations_calc_part1_pvm(void)
{
  int n=0;
  int mmin=gfish_index.indexmin();
  int mmax=gfish_index.indexmax();
  
  ivector iswitch(mmin,mmax);
  int i,j;
  for (i=mmin;i<=mmax;i++)
  {
    int rmin=gfish_index(i).indexmin();
    int rmax=gfish_index(i).indexmax();
    int fv=fish_flags(gfish_index(i,rmin),10);
    for (int j=rmin+1;j<=rmax;j++)
    {
      if (fv != fish_flags(gfish_index(i,j),10))
      {
        cerr << "sanity error in ff(i,10)" << 
         " for group " << i << " and fishery " << gfish_index(i,j) 
         << endl;
         ad_exit(1);
      }
    }
    iswitch(i)=fv;
    if (fv==1)
    {
      n+=num_real_grouped_fish_times(i);
    }
  }
  if (n==0)
  {
    cerr << "Error can not have all fish flags 10 =0 with"
      " af(104)=1" << endl;
    ad_exit(1);
  }
  imatrix indx(1,ngroups,1,num_real_grouped_fish_times);
  indx.initialize();
  int ii=1;
  for (j=1;j<=ngroups;j++)
  {
    if (iswitch(j))
    {
      for (i=1;i<=num_real_grouped_fish_times(j);i++)
      {
        indx(j,i)=ii++;
      }
    }
  }
  dvector x(1,n);
  MY_DOUBLE_TYPE c=0.05;

  do
  {
    int nopt=mfget_int_from_master();
    if (nopt<1)
      break;
    MY_DOUBLE_TYPE f=grouped_implicit_catchability_deviations_calc_part1_pvm(n,x,
      iswitch,indx,c);
    send_dvector_to_master(x);
  }
  while(1);

  if (!allocated(thread_xsave))
  {
    thread_xsave.allocate(1,n);
  }
  ad_exit(0);
}
#endif
