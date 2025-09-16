/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


#include "all.hpp"
#include <admodel.h>
#include "df22fun.h"
#include "df32fun.h"
#include "phi_stuff.h"
#include "tridiagonal_dmatrix.h"
//#include <newmult.hpp>
df2_two_variable d2phi(MY_DOUBLE_TYPE _p,MY_DOUBLE_TYPE _eta);
df3_two_variable d3phi(MY_DOUBLE_TYPE _p,MY_DOUBLE_TYPE _eta);
#include "makebig2.hpp"

extern MY_DOUBLE_TYPE SUMPEN;
static void df_lognormal_ss_multinomial_laplace_approximation_banded(void);

static void df_add_test_par(void)
{
  verify_identifier_string("u2");
  prevariable_position ypos=restore_prevariable_position();
  verify_identifier_string("u1");
  prevariable_position xpos=restore_prevariable_position();
  verify_identifier_string("u0");
  MY_DOUBLE_TYPE dfx=restore_prevariable_derivative(xpos);
  MY_DOUBLE_TYPE dfy=restore_prevariable_derivative(ypos);
  cout << " dfx " << dfx << " dfy = " << dfy;   
  dfx+=dfy;
  cout <<  " dfx+=dfy = " << dfx << endl;
  save_double_derivative(dfx,xpos);
}

void  add_test_par(const prevariable & _y,const prevariable& _x)
{
  ADUNCONST(prevariable,y)
  ADUNCONST(prevariable,x)
  cout << " x " << x << " y = " << y;  
  value(y)+=value(x);
  cout << " y+=x " << y << endl;
//    save_identifier_string("u0");
  const char * str1;
  str1="u0";
  char* strx1=const_cast <char*> (str1);
  save_identifier_string(strx1);
  x.save_prevariable_position();
//    save_identifier_string("u1");
  const char * str2;
  str2="u1";
  char* strx2=const_cast <char*> (str2);
  save_identifier_string(strx2);
  y.save_prevariable_position();
//    save_identifier_string("u2");
  const char * str3;
  str3="u2";
  char* strx3=const_cast <char*> (str3);
  save_identifier_string(strx3);
  
  gradient_structure::GRAD_STACK1->set_gradient_stack(df_add_test_par);
}

static void df_add_test_par1(void)
{
  verify_identifier_string("u2");
  prevariable_position ypos=restore_prevariable_position();
  verify_identifier_string("u1");
  prevariable_position xpos=restore_prevariable_position();
  verify_identifier_string("u0");
  MY_DOUBLE_TYPE dfx=restore_prevariable_derivative(xpos);
  MY_DOUBLE_TYPE dfy=restore_prevariable_derivative(ypos);
  cout << "11 dfx " << dfx << " dfy = " << dfy;   
  dfx+=dfy;
  cout <<  " dfx+=dfy = " << dfx << endl;
  save_double_derivative(dfx,xpos);
}

void  add_test_par1(const prevariable & _y,const prevariable& _x)
{
  ADUNCONST(prevariable,y)
  ADUNCONST(prevariable,x)
  cout << "11 x " << x << " y = " << y;  
  value(y)+=value(x);
  cout << " y+=x " << y << endl;
//    save_identifier_string("u0");
  const char * str4;
  str4="u0";
  char* strx4=const_cast <char*> (str4);
  save_identifier_string(strx4);
  x.save_prevariable_position();
//    save_identifier_string("u1");
  const char * str5;
  str5="u1";
  char* strx5=const_cast <char*> (str5);
  save_identifier_string(strx5);
  y.save_prevariable_position();
//    save_identifier_string("u2");
  const char * str6;
  str6="u2";
  char* strx6=const_cast <char*> (str6);
  save_identifier_string(strx6);
  
  gradient_structure::GRAD_STACK1->set_gradient_stack(df_add_test_par1);
}
static void df_check_test_par(void)
{
  verify_identifier_string("u2");
  prevariable_position ypos=restore_prevariable_position();
  verify_identifier_string("u1");
  MY_DOUBLE_TYPE dfy=restore_prevariable_derivative(ypos);
  cout << " dfy = " << dfy;   
}

void  check_test_par(const prevariable & _y)
{
  ADUNCONST(prevariable,y)
  cout << " y = " << y << endl;  
//    save_identifier_string("u1");
  const char * str7;
  str7="u1";
  char* strx7=const_cast <char*> (str7);
  save_identifier_string(strx7);
  y.save_prevariable_position();
//    save_identifier_string("u2");
  const char * str8;
  str8="u2";
  char* strx8=const_cast <char*> (str8);
  save_identifier_string(strx8);
  gradient_structure::GRAD_STACK1->set_gradient_stack(df_check_test_par);
}

static void df_check_test_par_vec(void)
{
  verify_identifier_string("w2");
  dvar_vector_position ypos=restore_dvar_vector_position();
  verify_identifier_string("w1");
  dvector dfy=restore_dvar_vector_derivatives(ypos);
  cout << " dfy = " << dfy << endl;   
}

void  check_test_par(const dvar_vector & _y)
{
  ADUNCONST(dvar_vector,y)
  cout << " y = " << y << endl;  
//    save_identifier_string("w1");
  const char * str9;
  str9="w1";
  char* strx9=const_cast <char*> (str9);
  save_identifier_string(strx9);
  y.save_dvar_vector_position();
//    save_identifier_string("w2");
  const char * str10;
  str10="w2";
  char* strx10=const_cast <char*> (str10);
  save_identifier_string(strx10);
  gradient_structure::GRAD_STACK1->set_gradient_stack(df_check_test_par_vec);
}

dvariable lognormal_multinomial_log_likelihood_cubic_spline
  (const dvariable& vN,const dvector& q,const dvar_vector& _vp,
  const dvector& _eps,int ic,
  banded_symmetric_dvar_matrix& vbsd2,const dvector& eta_hat)
{
  ADUNCONST(dvector,eps)
  ADUNCONST(dvar_vector,vp)


  MY_DOUBLE_TYPE N=value(vN);
  dvector p=value(vp);
  //banded_lower_triangular_dmatrix bltd = value(vbltd);
  
  //int mmin=p.indexmin();
  //int mmax=p.indexmax();
  //int nn=mmax-mmin+1;

  int pps=0;
  banded_symmetric_dmatrix bsd2=value(vbsd2);
  dvector eta_save;
  eta_save=eta_hat;
  eps=eta_hat;
  MY_DOUBLE_TYPE ln_det1=lognormal_ss_multinomial_laplace_approximation
    (N,q,p,eta_hat,bsd2);
  
  dvariable vlog_det=nograd_assign(ln_det1);
//    save_identifier_string("t4");
  const char * str11;
  str11="t4";
  char* strx11=const_cast <char*> (str11);
  save_identifier_string(strx11);
  eta_hat.save_dvector_value();
//    save_identifier_string("t3");
  const char * str12;
  str12="t3";
  char* strx12=const_cast <char*> (str12);
  save_identifier_string(strx12);
  eta_hat.save_dvector_position();
//    save_identifier_string("ps");
  const char * str13;
  str13="ps";
  char* strx13=const_cast <char*> (str13);
  save_identifier_string(strx13);
  vlog_det.save_prevariable_position();
//    save_identifier_string("rt");
  const char * str14;
  str14="rt";
  char* strx14=const_cast <char*> (str14);
  save_identifier_string(strx14);
  vN.save_prevariable_value();
//    save_identifier_string("s7");
  const char * str15;
  str15="s7";
  char* strx15=const_cast <char*> (str15);
  save_identifier_string(strx15);
  vN.save_prevariable_position();
//    save_identifier_string("s6");
  const char * str16;
  str16="s6";
  char* strx16=const_cast <char*> (str16);
  save_identifier_string(strx16);
  q.save_dvector_value();
//    save_identifier_string("s5");
  const char * str17;
  str17="s5";
  char* strx17=const_cast <char*> (str17);
  save_identifier_string(strx17);
  q.save_dvector_position();
//    save_identifier_string("s4");
  const char * str18;
  str18="s4";
  char* strx18=const_cast <char*> (str18);
  save_identifier_string(strx18);
  vp.save_dvar_vector_value();
//    save_identifier_string("s3");
  const char * str19;
  str19="s3";
  char* strx19=const_cast <char*> (str19);
  save_identifier_string(strx19);
  vp.save_dvar_vector_position();
//    save_identifier_string("s1");
  const char * str20;
  str20="s1";
  char* strx20=const_cast <char*> (str20);
  save_identifier_string(strx20);
  save_int_value(ic);
//    save_identifier_string("s14");
  const char * str21;
  str21="s14";
  char* strx21=const_cast <char*> (str21);
  save_identifier_string(strx21);
  bsd2.save_dmatrix_value();
//    save_identifier_string("s14a");
  const char * str22;
  str22="s14a";
  char* strx22=const_cast <char*> (str22);
  save_identifier_string(strx22);
  vbsd2.save_dvar_matrix_position();
  //save_pointer_value((void*)(pccsa));
//    save_identifier_string("s15");
  const char * str23;
  str23="s15";
  char* strx23=const_cast <char*> (str23);
  save_identifier_string(strx23);
  
  gradient_structure::GRAD_STACK1->
      set_gradient_stack
        (df_lognormal_ss_multinomial_laplace_approximation_banded);
  
  //cout << "BB log_det = " << vlog_det << endl;
  return vlog_det;
}
void check_deriv(banded_lower_triangular_dmatrix& L1,dvector& w);


MY_DOUBLE_TYPE lognormal_ss_multinomial_laplace_approximation(const MY_DOUBLE_TYPE& N,
  const dvector& q,const dvector& p,
  const dvector& eta,const banded_symmetric_dmatrix& stdsinv)
{
  int i,j,k;
  int mmin=q.indexmin();
  int mmax=q.indexmax();
  
  dvector t(mmin,mmax);
  dvector Nq=N*q;
  dvector grad(mmin,mmax);
  dvector w(mmin,mmax);
  int bw=stdsinv.bandwidth();
  banded_symmetric_dmatrix BHess(mmin,mmax,bw);
  banded_symmetric_dmatrix& BSinv=
    (banded_symmetric_dmatrix&) (stdsinv);
  BHess.initialize();
  MY_DOUBLE_TYPE v=0.0;
  for (int i=mmin;i<=mmax;i++)
  {
    t(i)=phi(p(i),eta(i));
    v+=t(i);
  }
  MY_DOUBLE_TYPE ml1=0.5*eta*(BSinv*eta);
  MY_DOUBLE_TYPE ml2=-Nq*log(t)+N*log(v); 
  MY_DOUBLE_TYPE mll= ml1+ml2;
  //double mll=-Nq*log(t)+N*log(v)+0.5*eta*(BSinv*eta); 
  int myflag=0;
  if (myflag)
  {
    int ierr=0;
    cout << 0.5*ln_det(make_dmatrix(BSinv),ierr) << endl;
  }
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

  for (int i=mmin;i<=mmax;i++)
  {
    for (int j=i;j<=min(mmax,i+bw-1);j++)
    {
      BHess(j,i)+=BSinv(j,i);
    }
  }
  // choleski decomp to get log det
  MY_DOUBLE_TYPE log_lik=mll;
  int ierr=0;
  MY_DOUBLE_TYPE ld=0.0;
  int minsave=BHess.indexmin();
  BHess.shift(1);
  int n=BHess.indexmax();
  banded_lower_triangular_dmatrix L(1,n,bw);
  L.initialize();

  MY_DOUBLE_TYPE tmp;
  if (BHess(1,1)<=0)
  {
    if (ierr==0)
      cerr << "Z Error matrix not positive definite in choleski_decomp"
        <<endl;
    ierr=1;
  }
  L(1,1)=sqrt(BHess(1,1));
  for (i=2;i<=bw;i++)
  {
    L(i,1)=BHess(i,1)/L(1,1);
  }

  for (i=2;i<=n;i++)
  {
    for (j=i-bw+1;j<=i-1;j++)
    {
      if (j>1)
      {	
        tmp=BHess(i,j);
        for (k=i-bw+1;k<=j-1;k++)
        {
	  if (k>0 && k>j-bw)
            tmp-=L(i,k)*L(j,k);
        }
        L(i,j)=tmp/L(j,j);
      }
    }
    tmp=BHess(i,i);
    for (k=i-bw+1;k<=i-1;k++)
    {
      if (k>0)	
        tmp-=L(i,k)*L(i,k);
    }
    if (tmp<=0)
    {
      if (ierr==0)
        cerr << "Z Error matrix not positive definite in choleski_decomp"
          <<endl;
      ierr=1;
    }
    L(i,i)=sqrt(tmp);
  }
  for (int i=mmin;i<=mmax;i++)
  {
    ld+=log(L(i,i));
  }
  int ii=0;
  //cout << "CCC  " << ld << " " << 0.5*ln_det(make_dmatrix(M),ii) << endl;

  // rank 1 update for determinant
  //dvector bx=solve(BHess,w);
  dvector bx=xlt1solve_solvet(L,w);
  tmp=1.0-(bx*w);
  if (tmp<=0.0)
  {
    cout << tmp << endl;
  }
  MY_DOUBLE_TYPE fpen=0.0;
  ld+=0.5*log(1.0-bx*w);
  int pps2=0;
  if (pps2)
  {
    dmatrix H=make_dmatrix(BHess);
    int sgn=0;
    cout << "A " << 0.5*ln_det(H,sgn) << endl;

    H-=outer_prod(w,w);
    cout << "B " << 0.5*ln_det(H,sgn) << endl;
    cout << "D " << ld << endl;
    ad_exit(1);
  }

  // KKKKKKKKKKKKKKKK
  log_lik=mll+ld;
  BHess.shift(minsave);
  L.shift(minsave);
  return log_lik;
}



static void df_lognormal_ss_multinomial_laplace_approximation_banded(void)
{
  int i,j,k;
  //double dflog_det=restore_prevariable_derivative();

  verify_identifier_string("s15");

  dvar_matrix_position BSinvpos=restore_dvar_matrix_position();
  verify_identifier_string("s14a");
  banded_symmetric_dmatrix BSinv=
    restore_banded_symmetric_dvar_matrix_value(BSinvpos);

  verify_identifier_string("s14");
  int ic=restore_int_value();
  verify_identifier_string("s1");
  dvar_vector_position vppos=restore_dvar_vector_position();
  verify_identifier_string("s3");
  dvector p=restore_dvar_vector_value(vppos);
  verify_identifier_string("s4");
  dvector_position qpos=restore_dvector_position();
  verify_identifier_string("s5");
  dvector q=restore_dvector_value(qpos);
  verify_identifier_string("s6");
  prevariable_position Npos=restore_prevariable_position();
  verify_identifier_string("s7");
  MY_DOUBLE_TYPE N=restore_prevariable_value();
  verify_identifier_string("rt");
  prevariable_position vcpos=restore_prevariable_position();
  MY_DOUBLE_TYPE dflog_det=restore_prevariable_derivative(vcpos);
  verify_identifier_string("ps");
  dvector_position etapos=restore_dvector_position();
  verify_identifier_string("t3");
  dvector eta=restore_dvector_value(etapos);
  verify_identifier_string("t4");
  int mmin=q.indexmin();
  int mmax=q.indexmax();
  int nn=mmax-mmin+1;

  banded_symmetric_dmatrix dfBSinv(BSinv.indexmin(),BSinv.indexmax(),
    BSinv.bandwidth());
  dfBSinv.initialize();
  banded_symmetric_dmatrix BHess(BSinv.indexmin(),BSinv.indexmax(),
    BSinv.bandwidth());
  BHess.initialize();
  banded_symmetric_dmatrix dfBHess(BSinv.indexmin(),BSinv.indexmax(),
    BSinv.bandwidth());
  dfBHess.initialize();

  
  dvector dfeta(mmin,mmax);
  dfeta.initialize();
  dvector t(mmin,mmax);
  dvector dft(mmin,mmax);
  dft.initialize();
  dvector Nq=N*q;
  MY_DOUBLE_TYPE v=0.0;
  MY_DOUBLE_TYPE dfv=0.0;
  for (i=mmin;i<=mmax;i++)
  {
    t(i)=phi(p(i),eta(i));
    //t(i)=p(i)*exp(eta(i));
    v+=t(i);
  }
  MY_DOUBLE_TYPE ml1=0.5*eta*(BSinv*eta);
  MY_DOUBLE_TYPE ml2=-Nq*log(t)+N*log(v);
  MY_DOUBLE_TYPE mll=ml1+ml2;
  //double mll=-Nq*log(t)+N*log(v)+0.5*eta*(BSinv*eta);
  // !!!!!!!!!GGGG
  MY_DOUBLE_TYPE dfmll=0.0;
  MY_DOUBLE_TYPE dfnv=0.0;
  MY_DOUBLE_TYPE dfnvv=0.0;
  //symmetric_tridiagonal_dmatrix dfBSinv(mmin,mmax);
  dfBSinv.initialize();
  //dmatrix dfSinv(mmin,mmax,mmin,mmax);
  //dfSinv.initialize();
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
  //dvector Nq(mmin,mmax);
  dvector Nqt(mmin,mmax);
  dvector Nqtt(mmin,mmax);
  dvector w(mmin,mmax);
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
  int bw=BHess.bandwidth();
  for (int i=mmin;i<=mmax;i++)
  {
    BHess(i,i)+=BSinv(i,i);
    for (int k=1;k<=min(mmax-i,bw-1);k++)
    {
      BHess(i+k,i)=BSinv(i+k,i);
    }
  }

  int n=mmax;

  dvector dfw(mmin,mmax);
  dvector dfx(mmin,mmax);
  banded_lower_triangular_dmatrix L1(1,n,bw);
  banded_lower_triangular_dmatrix dfL1(1,n,bw);
  dvector tmp(1,n);
  dmatrix tmp1(1,n,1,n);
  dmatrix dftmp1(1,n,1,n);
  dvector dftmp(1,n);
  tmp.initialize();
  tmp1.initialize();
  dftmp.initialize();
  dftmp1.initialize();
  dfL1.initialize();
  L1.initialize();
  int ierr=0;
  MY_DOUBLE_TYPE log_det1=mll;
  {
    if (BHess(1,1)<=0)
    {
      cerr << "Z Error matrix not positive definite in choleski_decomp"
          <<endl;
      ierr=1;
    }
    L1(1,1)=sqrt(BHess(1,1));
    for (i=2;i<=bw;i++)
    {
      L1(i,1)=BHess(i,1)/L1(1,1);
    }
  
    for (i=2;i<=n;i++)
    {
      for (j=i-bw+1;j<=i-1;j++)
      {
        if (j>1)
        {	
          tmp1(i,j)=BHess(i,j);
          for (k=i-bw+1;k<=j-1;k++)
          {
  	  if (k>0 && k>j-bw)
              tmp1(i,j)-=L1(i,k)*L1(j,k);
          }
          L1(i,j)=tmp1(i,j)/L1(j,j);
        }
      }
      tmp(i)=BHess(i,i);
      for (k=i-bw+1;k<=i-1;k++)
      {
        if (k>0)	
          tmp(i)-=L1(i,k)*L1(i,k);
      }
      if (tmp(i)<=0)
      {
        cerr << "Z Error matrix not positive definite in choleski_decomp"
              <<endl;
        ierr=1;
      }
      L1(i,i)=sqrt(tmp(i));
    }
    for (int i=mmin;i<=mmax;i++)
    {
      log_det1+=log(L1(i,i));
    }
  }

  dmatrix zz=make_dmatrix(L1);
  //cout << " VVV0 " << norm2(make_dmatrix(BHess)-zz*trans(zz)) << endl;
  dvector x1=solve(BHess,w);
  dvector x=xlt1solve_solvet(L1,w);
  //cout << " VVV " << norm2(x1-x) << endl;
  MY_DOUBLE_TYPE uu=1.0-x*w;
  MY_DOUBLE_TYPE fpen=0.0;
  //double puu=posfun(uu,.01,fpen);
  // XXXXXXXXXXXXx
  //log_det1+=0.5*log(puu);
  log_det1+=0.5*log(uu);
  MY_DOUBLE_TYPE log_det=log_det1;
  //cout << "MM log_det = " << log_det << " mll = " << mll 
    //   << " log_det - mll " << log_det-mll << endl;
  MY_DOUBLE_TYPE dfN=0;
  dvector dfp(mmin,mmax);
  dfp.initialize();

 //*******************************************************************8
  {
    dfx.initialize();
    dfw.initialize();
    //double log_det=ml+log_det1;
    MY_DOUBLE_TYPE dflog_det1=dflog_det;
    dfmll+=dflog_det;
    dflog_det=0.0;
    MY_DOUBLE_TYPE dfuu=0.0;
    //log_det1+=0.5*log(puu);
    dfuu+=dflog_det1*0.5/uu;
    //double puu=posfun(uu,.01);
    //dfuu+=dfpuu*dfposfun(uu,.01);
    //dfpuu=0.0;
    //double uu=1.0-x*w;
    dfx-=dfuu*w;
    dfw-=dfuu*x;
    dfuu=0.0;
    xdflt1solve_solvet(L1,w,dfw,dfx,dfL1);
    
    for (i=1;i<=n;i++)
    {
      // EEEEEEE
      //log_det1+=log(L1(i,i));
      dfL1(i,i)+=dflog_det1/L1(i,i);
    }
    //double log_det1=mll;
    //dfmll+=dflog_det1;
    dflog_det1=0.0;

  
    for (i=n;i>=2;i--)
    {
      //L1(i,i)=sqrt(tmp(i));
      dftmp(i)+=dfL1(i,i)/(2.0*L1(i,i));
      dfL1(i,i)=0.0;

      for (k=i-1;k>=i-bw+1;k--)
      {
        if (k>0)	
        {
          tmp(i)-=L1(i,k)*L1(i,k);
          dfL1(i,k)-=2.*dftmp(i)*L1(i,k);
        }
      }
      //tmp(i)=BHess(i,i);
      dfBHess(i,i)+=dftmp(i);
      dftmp(i)=0.0;
      for (j=i-1;j>=i-bw+1;j--)
      {
        if (j>1)
        {
          //L1(i,j)=tmp1(i,j)/L1(j,j);
          MY_DOUBLE_TYPE linv=1./L1(j,j);
          dftmp1(i,j)+=dfL1(i,j)*linv;
          dfL1(j,j)-=dfL1(i,j)*tmp1(i,j)*linv*linv;
          dfL1(i,j)=0.0;
          for (k=j-1;k>=i-bw+1;k--)
          {
            if (k>0 && k>j-bw)
            {
              // tmp1(i,j)-=L1(i,k)*L1(j,k);
              dfL1(i,k)-=dftmp1(i,j)*L1(j,k);
              dfL1(j,k)-=dftmp1(i,j)*L1(i,k);
            }
          }
  
          //tmp1(i,j)=BHess(i,j);
          dfBHess(i,j)+=dftmp1(i,j);
          dftmp1(i,j)=0.0;
        }
      }
    }
    MY_DOUBLE_TYPE xlinv=1./L1(1,1);
    MY_DOUBLE_TYPE xlinv2=xlinv*xlinv;
    for (i=bw;i>=2;i--)
    {
      //L1(i,1)=BHess(i,1)/L1(1,1);
      dfL1(1,1)-=dfL1(i,1)*BHess(i,1)*xlinv2;
      dfBHess(i,1)+=dfL1(i,1)/L1(1,1);
      dfL1(i,1)=0.0;
    }
    //L1(1,1)=sqrt(BHess(1,1));
    dfBHess(1,1)+=dfL1(1,1)/(2.*L1(1,1));
    dfL1(1,1)=0.0;
  }
  //BHess(mmax,mmax)+=BSinv(mmax,mmax);
  // why is this HERE
  //dfBHess.initialize();
  
  for (int i=mmin;i<=mmax;i++)
  {
    //BHess(i,i)+=BSinv(i,i);
    dfBSinv(i,i)+=dfBHess(i,i);
    for (int k=1;k<=min(mmax-i,bw-1);k++)
    {
      //BHess(i+k,i)=BSinv(i+k,i);
      dfBSinv(i+k,i)+=dfBHess(i+k,i);
      dfBHess(i+k,i)=0.0;
    }
  }
  dvector dfNq(mmin,mmax);
  dfNq.initialize();
  dvector dfd1(mmin,mmax);
  dvector dfd2(mmin,mmax);
  dfd1.initialize();
  dfd2.initialize();
  MY_DOUBLE_TYPE dfNvv=0.0;
  dvector dfNqt(mmin,mmax);
  dvector dfNqtt(mmin,mmax);
  dfNqtt.initialize();
  dfNqt.initialize();
  MY_DOUBLE_TYPE dfNv=0;

  // XXXXX
  //dfw.initialize();

  for (int i=mmax;i>=mmin;i--)
  {
    //w(i)=sqrt(Nvv)*d1(i);
    dfNvv+=0.5*dfw(i)/sqrt(Nvv)*d1(i);
    dfd1(i)+=dfw(i)*sqrt(Nvv);
    dfw(i)=0.0;
    //BHess(i,i)+=Nqtt(i)*square(d1(i));
    dfNqtt(i)+=dfBHess(i,i)*square(d1(i));
    dfd1(i)+=2.0*dfBHess(i,i)*Nqtt(i)*d1(i);
    
    //BHess(i,i)+=(-Nqt(i)+Nv)*d2(i);
    dfNqt(i)-=dfBHess(i,i)*d2(i);
    dfNv+=dfBHess(i,i)*d2(i);
    dfd2(i)+=dfBHess(i,i)*(-Nqt(i)+Nv);
  }

  for (int i=mmin;i<=mmax;i++)
  {
    //Nqtt(i)=Nqt(i)/t(i);
    dfNqt(i)+=dfNqtt(i)/t(i);
    dft(i)-=dfNqtt(i)*Nqt(i)/(t(i)*t(i));
    dfNqtt(i)=0.0;
  }
  for (int i=mmin;i<=mmax;i++)
  {
    //Nqt(i)=Nq(i)/t(i);
    dfNq(i)+=dfNqt(i)/t(i);
    dft(i)-=dfNqt(i)*Nq(i)/(t(i)*t(i));
    dfNqt(i)=0.0;
  }
  dfBHess.initialize();

  //double Nvv=Nv/v;
  dfNv+=dfNvv/v;
  dfv-=dfNvv*Nv/(v*v);
  dfNvv=0.0;

  // DDDDD
  
  //double Nv=N/v;
  dfN+=dfNv/v;
  dfv-=dfNv*N/(v*v);
  dfNv=0.0;
 
 
  //cout << "H1H1H1 " << dfeta(1) << " " << dfd1(i) << "  " << dfd2(i) << endl;

  for (int i=mmin;i<=mmax;i++)
  {
    df3_two_variable tmp=d3phi(p(i),eta(i));
    //d2(i)=d2phi_deta2(p(i),eta(i));
    //dfp(i)+=dfd2(i)*d2phi_deta2_1(p(i),eta(i));
    //dfeta(i)+=dfd2(i)*d2phi_deta2_2(p(i),eta(i));
    dfp(i)+=dfd2(i)* *tmp.get_u_xyy();
    dfeta(i)+=dfd2(i)* *tmp.get_u_yyy();
    dfd2(i)=0.0;
    //d1(i)=dphi_deta(p(i),eta(i));
    //dfp(i)+=dfd1(i)*dphi_deta_1(p(i),eta(i));
    //dfeta(i)+=dfd1(i)*dphi_deta_2(p(i),eta(i));
    dfp(i)+=dfd1(i)*  *tmp.get_u_xy();
    dfeta(i)+=dfd1(i)* *tmp.get_u_yy();
    dfd1(i)=0.0;
  }
  //cout << "HHH " << dfeta(1) << " " << dfd1(i) << "  " << dfd2(i) << endl;


  //double mll=-Nq*log(t)+N*log(v)+0.5*eta*(BSinv*eta); 

  dfNq-=dfmll*log(t);
  dft-=dfmll*elem_div(Nq,t);
  dfN+=dfmll*log(v);
  dfv+=dfmll*N/v;
  dfeta+=BSinv*eta;
  dfmll=0.0;
  for (int i=mmin;i<=mmax;i++)
  {
    dfBSinv(i,i)+=0.5*eta(i)*eta(i);
    for (int k=1;k<=min(mmax-i,bw-1);k++)
    {
      dfBSinv(i+k,i)+=eta(i+k)*eta(i);
    }
  }
  

  for (i=mmax;i>=mmin;i--)
  {
    df2_two_variable tmp=d2phi(p(i),eta(i));
    //v+=t(i);
    dft(i)+=dfv;
    //t(i)=phi(p(i),eta(i));
    //dfp(i)+=dft(i)*dphi_dp(p(i),eta(i));
    dfp(i)+=dft(i)* tmp.get_u_x();
    //dfp(i)+=dft(i)*dphi_deta(p(i),eta(i))*psiprime(p(i),eta(i));
    //dfeta(i)+=dft(i)*dphi_deta(p(i),eta(i));
    dfeta(i)+=dft(i)* tmp.get_u_y();
    dft(i)=0.0;
  }
  //dvector Nq=N*q;
  dfN+=dfNq*q;

  // v=0.0;
  dfv=0.0;
  // OOOOOOOOOOOOOOOOOOOOOOOOO
  //cout << "AA09 " << "dfBSinv(1,1) "  << dfBSinv(1,1) << endl;
  //df_mult_trans_mult(bltd,BSinv,dfBSinv,dfbltd);
  //cout << "AA08 " << "dfbltd(1,1) "  << dfbltd(1,1)  << endl;
  //dfbltd.initialize();
  
  
 
  dvector cd=get_cross_derivatives(L1,BHess,w,q,t,v,eta,N,dfeta,p);
  // add in the cross derivatives
  // KKKKKKKKKKKKKKKKKKKKKKKKK
  //cd.initialize();
 /*
  dfN+=cd(1);
  dfp+=cd(2,n+1).shift(1);
  for (int i=1;i<=n-1;i++)
  {
    //dfSinv(i,j)+=cd(1+n+(i-1)*n+j);
    dfBSinv(i,i)+=cd(1+n+(i-1)*n+i);
    dfBSinv(i+1,i)+=2.0*cd(1+n+i*n+i);
  }
  dfBSinv(n,n)+=cd(1+n+(n-1)*n+n);
  */
  
 
  dfN+=cd(1);
  dfp+=cd(2,n+1).shift(1);
  for (int i=1;i<=n;i++)
  {
    dfBSinv(i,i)+=cd(1+n+(i-1)*n+i);
    for (int k=1;k<=min(mmax-i,bw-1);k++)
    {
      //dfBSinv(i+1,i)+=2.0*cd(1+n+i*n+i);
      dfBSinv(i+k,i)+=2.0*cd(1+n+(i+k-1)*n+i);
    }
  }
  //dfBSinv(n,n)+=cd(1+n+(n-1)*n+n);
  // OOOOOOOOOOOOOOOOOOOOOOOOO
  //ad_exit(1);
  // OOOOOOOOOOOOOOOOOOOOOOOOO
 
  //cout << "AA11 " << "dfBSinv(1,1) "  << dfBSinv(1,1) << endl;
  //df_mult_trans_mult(bltd,BSinv,dfBSinv,dfbltd);
  //cout << "AA11 " << "dfbltd(1,1) "  << dfbltd(1,1)  << endl;
  //ad_exit(1);

  save_double_derivative(dfN,Npos);
  dfp.save_dvector_derivatives(vppos);
  dfBSinv.save_dmatrix_derivatives(BSinvpos);
}


dvector get_cross_derivatives(
  const banded_lower_triangular_dmatrix L1,
  const banded_symmetric_dmatrix& BH,
  const dvector& w,
  const dvector& q,
  const dvector& t,const MY_DOUBLE_TYPE v,const dvector& eta,const MY_DOUBLE_TYPE& N,
  const dvector dfeta,const dvector& p)
{
  //lower_triangular dmatrix LT(1,n);
  int n=q.indexmax();

  dmatrix NN(1,n,1,1+n);
  NN.initialize();
  //dmatrix NN(1,n,1+n+n*n);
 /*
  // take derivative of this wrt N

  for (int i=mmin;i<=mmax;i++)
  {
    grad(i)=(-N*q(i)/t(i)+N/v)*dphi_deta(p(i),eta(i));
  }
 */
 
  for (int i=1;i<=n;i++)
  {
    df2_two_variable tmp=d2phi(p(i),eta(i));
    //NN(i,1)=(-q(i)/t(i)+1.0/v)*dphi_deta(p(i),eta(i));
    NN(i,1)=(-q(i)/t(i)+1.0/v)* tmp.get_u_y();
  }
  df2_two_vector xxx(1,n);

  for (int i=1;i<=n;i++)
  {
    xxx(i)=d2phi(p(i),eta(i));
  }
  for (int i=1;i<=n;i++)
  {
    //NN(i,i+1)+=N*q(i)/square(t(i))*dphi_dp(p(i),eta(i))
    //  *dphi_deta(p(i),eta(i));
    NN(i,i+1)+=N*q(i)/square(t(i))* xxx(i).get_u_x()
      * xxx(i).get_u_y();

    //NN(i,i+1)+=N*(-q(i)/t(i)+1.0/v)*dphi_deta_1(p(i),eta(i));
    NN(i,i+1)+=N*(-q(i)/t(i)+1.0/v)* xxx(i).get_u_xy();
    for (int j=1;j<=n;j++)
    {
      //df2_two_variable tmpj=d2phi(p(j),eta(j));
      //NN(i,j+1)-= N/square(v)*dphi_dp(p(j),eta(j))*dphi_deta(p(i),eta(i));
      NN(i,j+1)-= N/square(v)*xxx(j).get_u_x()*xxx(i).get_u_y();
    }
  }
 

  dmatrix Hess(1,n,1,n);
  int mmin=1;
  int mmax=n;

  Hess.initialize();
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
    
 /*
  dvector u(mmin,mmax);
  for (int j=1;j<=n;j++)
  {
    //u(j)=psiprime(p(j),eta(j));
  }
 */

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
  dvector tmp(1,1+n+n*n);
  tmp.initialize();
  //w.initialize();
  //dmatrix Minv=inv(BH,w);
  dmatrix Minv=inv(L1,w);
  //dmatrix uhatp=-Minv*NN;
  dmatrix uhatp=-xlt1solve_solvet(L1,w,NN);
  tmp(1,n+1)=dfeta*uhatp;
 
 /*
  for (int i=1;i<=n-1;i++)
  {
    dvector t=2.0*Minv(i)*eta(i);
    tmp(1+n+(i-1)*n+i)=-0.5*(dfeta*t);
    t=Minv(i+1)*eta(i)+Minv(i)*eta(i+1);
    tmp(1+n+i*n+i)=-0.5*(dfeta*t);
  }
  {
    dvector t=2.0*Minv(n)*eta(n);
    tmp(1+n+(n-1)*n+n)=-0.5*(dfeta*t);
  }
  */
  int bw=L1.bandwidth();
  for (int i=1;i<=n;i++)
  {
    dvector t=2.0*Minv(i)*eta(i);
    tmp(1+n+(i-1)*n+i)=-0.5*(dfeta*t);
    for (int k=1;k<=min(mmax-i,bw-1);k++)
    {
      t=Minv(i+k)*eta(i)+Minv(i)*eta(i+k);
      tmp(1+n+(i+k-1)*n+i)=-0.5*(dfeta*t);
    }
  }
  return tmp;
}

void xdflt1solve_solvet(
  const banded_lower_triangular_dmatrix & L,
  const dvector& y,
  const dvector& dfx,const dvector& dfy,
  const banded_lower_triangular_dmatrix & dfL
);

dvector lt1solve_solvet(
  const banded_lower_triangular_dmatrix & L,
  const dvector& y
);

void xdflt1solve_solvet(
  const banded_lower_triangular_dmatrix & L,
  const dvector& w,
  const dvector& _dfw,const dvector& _dfx,
  const banded_lower_triangular_dmatrix & _dfL)
{
  ADUNCONST(dvector,dfx)
  ADUNCONST(dvector,dfw)
  ADUNCONST(banded_lower_triangular_dmatrix,dfL)
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  int bw=L.bandwidth();
  dvector v(mmin,mmax);
  dvector x(mmin,mmax);
  dvector asum(mmin,mmax);
  dvector bsum(mmin,mmax);
  v.initialize();
  asum.initialize();
  bsum.initialize();
  v(mmin)=w(mmin)/L(mmin,mmin);
  for (int i=mmin+1;i<=mmax;i++)
  {
    for (int k=1;k<=min(i-mmin,bw-1);k++)
    {
      asum(i)+=L(i,i-k)*v(i-k);
    }
    v(i)=(w(i)-asum(i))/L(i,i);
  }
 
  x(mmax)=v(mmax)/L(mmax,mmax);
  for (int i=mmax-1;i>=mmin;i--)
  {
    for (int k=1;k<=min(mmax-i,bw-1);k++)
    {
      bsum(i)+=L(i+k,i)*x(i+k);
    }
    x(i)=(v(i)-bsum(i))/L(i,i);
  }
  //return x;
  // get dfx;
  //dvector dfx(mmin,mmax);
  //dvector dfw(mmin,mmax);
  dvector dfv(mmin,mmax);
  dvector dfasum(mmin,mmax);
  dvector dfbsum(mmin,mmax);
  //banded_loxer_triangular_dmatrix dfL(mmin,mmax,2);
  dfv.initialize();
  //dfw.initialize();
  dfasum.initialize();
  dfbsum.initialize();
  //dfL.initialize();
  
  for (int i=mmin;i<=mmax-1;i++)
  {
    //x(i)=(v(i)-bsum(i))/L(i,i);
    dfv(i)+=dfx(i)/L(i,i);
    dfbsum(i)-=dfx(i)/L(i,i);
    dfL(i,i)-=dfx(i)*x(i)/L(i,i);
    dfx(i)=0.0;
    //bsum(i)=L(i+1,i)*x(i+1);
    //for (int k=1;k<=min(mmax-i,bw-1);k++)
    //{
    //  bsum(i)+=L(i+k,i)*x(i+k);
    //}
    for (int k=min(mmax-i,bw-1);k>=1;k--)
    {
      //  bsum(i)+=L(i+k,i)*x(i+k);
      dfL(i+k,i)+=dfbsum(i)*x(i+k);
      dfx(i+k)+=dfbsum(i)*L(i+k,i);
    }
    dfbsum(i)=0.0;
  }
  //x(mmax)=v(mmax)/L(mmax,mmax);
  dfv(mmax)+=dfx(mmax)/L(mmax,mmax);
  dfL(mmax,mmax)-=dfx(mmax)*x(mmax)/L(mmax,mmax);
  dfx(mmax)=0.0;
  for (int i=mmax;i>=mmin+1;i--)
  {
    //v(i)=(w(i)-asum(i))/L(i,i);
    dfw(i)+=dfv(i)/L(i,i);
    dfasum(i)-=dfv(i)/L(i,i);
    dfL(i,i)-=dfv(i)*v(i)/L(i,i);
    dfv(i)=0.0;
    //for (int k=1;k<=min(i-mmin,bw-1);k++)
    //{
    //  asum(i)+=L(i,i-k)*v(i-k);
    //}
    for (int k=min(i-mmin,bw-1);k>=1;k--)
    {
      dfL(i,i-k)+=dfasum(i)*v(i-k);
      dfv(i-k)+=dfasum(i)*L(i,i-k);
    }
    dfasum(i)=0.0;
  }
  //v(mmin)=w(mmin)/L(mmin,mmin);
  dfw(mmin)+=dfv(mmin)/L(mmin,mmin);
  dfL(mmin,mmin)-=dfv(mmin)*v(mmin)/L(mmin,mmin);
  dfv(mmin)=0.0;

}
void xdflt1solve_solvet_1(
  const banded_lower_triangular_dmatrix & L,
  const dvector& w,
  const dvector& _dfw,const dvector& _dfx,
  const banded_lower_triangular_dmatrix & _dfL)
{
  ADUNCONST(dvector,dfx)
  ADUNCONST(dvector,dfw)
  ADUNCONST(banded_lower_triangular_dmatrix,dfL)
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  int bw=L.bandwidth();
  dvector v(mmin,mmax);
  dvector x(mmin,mmax);
  dvector asum(mmin,mmax);
  dvector bsum(mmin,mmax);
  v.initialize();
  asum.initialize();
  bsum.initialize();
  v(mmin)=w(mmin)/L(mmin,mmin);
  for (int i=mmin+1;i<=mmax;i++)
  {
    for (int k=1;k<=min(i-mmin,bw-1);k++)
    {
      asum(i)+=L(i,i-k)*v(i-k);
    }
    v(i)=(w(i)-asum(i))/L(i,i);
  }
  x=v;
 /*
  x(mmax)=v(mmax)/L(mmax,mmax);
  for (int i=mmax-1;i>=mmin;i--)
  {
    for (int k=1;k<=min(mmax-i,bw-1);k++)
    {
      bsum(i)+=L(i+k,i)*x(i+k);
    }
    x(i)=(v(i)-bsum(i))/L(i,i);
  }
  */
  //return x;
  // get dfx;
  //dvector dfx(mmin,mmax);
  //dvector dfw(mmin,mmax);
  dvector dfv(mmin,mmax);
  dvector dfasum(mmin,mmax);
  dvector dfbsum(mmin,mmax);
  //banded_loxer_triangular_dmatrix dfL(mmin,mmax,2);
  dfv.initialize();
  //dfw.initialize();
  dfasum.initialize();
  dfbsum.initialize();
  //dfL.initialize();
  
  /*
  for (int i=mmin;i<=mmax-1;i++)
  {
    //x(i)=(v(i)-bsum(i))/L(i,i);
    dfv(i)+=dfx(i)/L(i,i);
    dfbsum(i)-=dfx(i)/L(i,i);
    dfL(i,i)-=dfx(i)*x(i)/L(i,i);
    dfx(i)=0.0;
    //bsum(i)=L(i+1,i)*x(i+1);
    //for (int k=1;k<=min(mmax-i,bw-1);k++)
    //{
    //  bsum(i)+=L(i+k,i)*x(i+k);
    //}
    for (int k=min(mmax-i,bw-1);k>=1;k--)
    {
      //  bsum(i)+=L(i+k,i)*x(i+k);
      dfL(i+k,i)+=dfbsum(i)*x(i+k);
      dfx(i+k)+=dfbsum(i)*L(i+k,i);
    }
    dfbsum(i)=0.0;
  }
  //x(mmax)=v(mmax)/L(mmax,mmax);
  dfv(mmax)+=dfx(mmax)/L(mmax,mmax);
  dfL(mmax,mmax)-=dfx(mmax)*x(mmax)/L(mmax,mmax);
  */
  dfv=dfx;
  dfx(mmax)=0.0;
  for (int i=mmax;i>=mmin+1;i--)
  {
    //v(i)=(w(i)-asum(i))/L(i,i);
    dfw(i)+=dfv(i)/L(i,i);
    dfasum(i)-=dfv(i)/L(i,i);
    dfL(i,i)-=dfv(i)*v(i)/L(i,i);
    dfv(i)=0.0;
    //for (int k=1;k<=min(i-mmin,bw-1);k++)
    //{
    //  asum(i)+=L(i,i-k)*v(i-k);
    //}
    for (int k=min(i-mmin,bw-1);k>=1;k--)
    {
      dfL(i,i-k)+=dfasum(i)*v(i-k);
      dfv(i-k)+=dfasum(i)*L(i,i-k);
    }
    dfasum(i)=0.0;
  }
  //v(mmin)=w(mmin)/L(mmin,mmin);
  dfw(mmin)+=dfv(mmin)/L(mmin,mmin);
  dfL(mmin,mmin)-=dfv(mmin)*v(mmin)/L(mmin,mmin);
  dfv(mmin)=0.0;

}

dvector xlt1solve_solvet(
  const banded_lower_triangular_dmatrix & L,
  const dvector& w)
{
  int bw=L.bandwidth();
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  dvector v(mmin,mmax);
  dvector x(mmin,mmax);
  dvector asum(mmin,mmax);
  dvector bsum(mmin,mmax);
  v.initialize();
  asum.initialize();
  bsum.initialize();
  v(mmin)=w(mmin)/L(mmin,mmin);
  for (int i=mmin+1;i<=mmax;i++)
  {
    for (int k=1;k<=min(i-mmin,bw-1);k++)
    {
      asum(i)+=L(i,i-k)*v(i-k);
    }
    v(i)=(w(i)-asum(i))/L(i,i);
  }
  //cout << " BB " << v << endl;
  x(mmax)=v(mmax)/L(mmax,mmax);
  for (int i=mmax-1;i>=mmin;i--)
  {
    for (int k=1;k<=min(mmax-i,bw-1);k++)
    {
      bsum(i)+=L(i+k,i)*x(i+k);
    }
    x(i)=(v(i)-bsum(i))/L(i,i);
  }
  return x;
}
dvector xlt1solve_solvet_1(
  const banded_lower_triangular_dmatrix & L,
  const dvector& w)
{
  int bw=L.bandwidth();
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  dvector v(mmin,mmax);
  dvector x(mmin,mmax);
  dvector asum(mmin,mmax);
  dvector bsum(mmin,mmax);
  v.initialize();
  asum.initialize();
  bsum.initialize();
  v(mmin)=w(mmin)/L(mmin,mmin);
  for (int i=mmin+1;i<=mmax;i++)
  {
    for (int k=1;k<=min(i-mmin,bw-1);k++)
    {
      asum(i)+=L(i,i-k)*v(i-k);
    }
    v(i)=(w(i)-asum(i))/L(i,i);
  }
  x=v;
  return x;
  //cout << " BB " << v << endl;
  x(mmax)=v(mmax)/L(mmax,mmax);
  for (int i=mmax-1;i>=mmin;i--)
  {
    for (int k=1;k<=min(mmax-i,bw-1);k++)
    {
      bsum(i)+=L(i+k,i)*x(i+k);
    }
    x(i)=(v(i)-bsum(i))/L(i,i);
  }
  return x;
}

void check_deriv(banded_lower_triangular_dmatrix& L1,dvector& w)
{
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  int bw=L1.bandwidth();
  banded_lower_triangular_dmatrix  dfL1(mmin,mmax,bw);
  MY_DOUBLE_TYPE eta=1.e-4;
  MY_DOUBLE_TYPE ds=L1(1,1);
  dvector x=xlt1solve_solvet_1(L1,w);
  dmatrix M=make_dmatrix(L1);
  dvector y=solve(trans(M),solve(M,w));
  cout << norm2(x-y) << endl;
  L1(1,1)+=eta;
  dvector xu=xlt1solve_solvet_1(L1,w);
  MY_DOUBLE_TYPE fu=norm2(xu);
  L1(1,1)=ds;
  L1(1,1)-=eta;
  dvector xl=xlt1solve_solvet_1(L1,w);
  MY_DOUBLE_TYPE fl=norm2(xl);
  L1(1,1)=ds;
  cout << (fu-fl)/(2.0*eta)<< endl;
  dvector dfw(mmin,mmax);
  dfw.initialize();
  dfL1.initialize();
  dvector dfx(mmin,mmax);
  dfx=2.0*x;
  
  xdflt1solve_solvet_1(L1,w,dfw,dfx,dfL1);
  
  cout << dfL1(1,1) << endl;
}
    //dvector x=lt1solve_solvet(L1,w);
