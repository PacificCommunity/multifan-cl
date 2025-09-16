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
#include "phi_stuff.h"
#include "tridiagonal_dmatrix.h"
//#include <newmult.hpp>
df2_two_variable d2phi(MY_DOUBLE_TYPE _p,MY_DOUBLE_TYPE _eta);
df3_two_variable d3phi(MY_DOUBLE_TYPE _p,MY_DOUBLE_TYPE _eta);
  // dmatrix make_dmatrix(banded_lower_triangular_dmatrix B);
  // 
   extern MY_DOUBLE_TYPE SUMPEN;
  // 
  // dmatrix make_dmatrix(banded_lower_triangular_dmatrix B)
  // {
  //   int mmin=B.indexmin();
  //   if (mmin !=1)
  //   {
  //     ad_exit(1);
  //   }
  //   int mmax=B.indexmax();
  //   int bw=B.bandwidth();
  //   dmatrix M(1,mmax,1,mmax);
  //   M.initialize();
  //   for (int i=1;i<=mmax;i++)
  //   {
  //     for (int j=max(i-bw+1,1);j<=i;j++)
  //     {
  //       M(i,j)=B(i,j);
  //     }
  //   }
  //   return M;
  // }
  // 
  // dmatrix make_ar1matrix(int mmax,MY_DOUBLE_TYPE rho,MY_DOUBLE_TYPE var)
  // {
  //   dmatrix M(1,mmax,1,mmax);
  //   dvector a(0,mmax-1);
  //   a(0)=1.0;
  //   for (int i=1;i<mmax;i++)
  //   {
  //     a(i)=a(i-1)*rho;
  //   }
  //   for (int i=1;i<=mmax;i++)
  //   {
  //     for (int j=1;j<=i;j++)
  //     {
  //       M(i,j)=a(i-j);
  //       M(j,i)=M(i,j);
  //     }
  //   }
  //   M*=var/(1-rho*rho);
  //   return M;
  // }
  //     
  // 
  // static void  xxx(int ic) {ic++;}
   void  xdf_get_tridag_BH(MY_DOUBLE_TYPE crho,
     MY_DOUBLE_TYPE cvar,int mmin,int mmax,
     const symmetric_tridiagonal_dmatrix&  dfstsinv,const MY_DOUBLE_TYPE& _crho,
     const MY_DOUBLE_TYPE& _cvar);
   
   dmatrix xinv(const banded_lower_triangular_dmatrix & L,
     const dvector& ww);
   
   void yxdflt1solve_solvet(
     const banded_lower_triangular_dmatrix & L,
     const dvector& y,
     const dvector& dfx,const dvector& dfy,
     const banded_lower_triangular_dmatrix & dfL
     );
  // 
   dvector yxlt1solve_solvet(
     const banded_lower_triangular_dmatrix & L,
     const dvector& y);
  // 
   dmatrix yxlt1solve_solvet(
     const banded_lower_triangular_dmatrix & L,
     const dvector& ww,
     const dmatrix& N);
  // 
  // dvector lt1solve_solvet(const dmatrix& L,const dvector& y);
  // void dflt1solve_solvet(const dmatrix& L,const dvector& y,
  //   const dvector& dfx,const dvector& dfy,const dmatrix& dfL);
   dvector lognormal_ss_multinomial_newton_raphson(const MY_DOUBLE_TYPE N,
     const dvector& q,const dvector& p,
     const banded_lower_triangular_dmatrix & _chinv,
     const symmetric_tridiagonal_dmatrix& lensinv,const dvector& ee);
  // 
   dvector lognormal_ss_multinomial_newton_raphson(const MY_DOUBLE_TYPE N,
     const dvector& q,const dvector& p,
     const banded_lower_triangular_dmatrix & _chinv,
     const symmetric_tridiagonal_dmatrix& lensinv);
  // 
   MY_DOUBLE_TYPE xlognormal_ss_multinomial_laplace_approximation(const MY_DOUBLE_TYPE& N,
     const dvector& q,const dvector& p,
    const dvector& eta,const symmetric_tridiagonal_dmatrix& stdsinv);
  // 
  // dvector get_cross_derivatives(const dmatrix& M,const dvector& q,
  //   const dvector& t,const MY_DOUBLE_TYPE v,const dvector& eta,const MY_DOUBLE_TYPE& N,
  //   const dvector dfeta,const dvector& p);
   dvector xget_cross_derivatives(
     const banded_lower_triangular_dmatrix L1,
     const symmetric_tridiagonal_dmatrix& BH,
     const dvector& w,
     const dvector& q,
     const dvector& t,const MY_DOUBLE_TYPE v,const dvector& eta,const MY_DOUBLE_TYPE& N,
     const dvector dfeta,const dvector& p);
  // 
  // 
  // //extern int global_except_switch;
  // 
  // dmatrix inv(const symmetric_tridiagonal_dmatrix& _BH,const dvector&w)
  // {
  //   ADUNCONST(symmetric_tridiagonal_dmatrix,BH)
  //   int mmin=BH.indexmin();
  //   int mmax=BH.indexmax();
  //   dmatrix M(mmin,mmax,mmin,mmax);
  //   M.initialize();
  //   for (int i=mmin;i<=mmax-1;i++)
  //   {
  //     M(i,i)=BH(i,i);
  //     M(i+1,i)=BH(i+1,i);
  //     M(i,i+1)=BH(i+1,i);
  //   }
  //   M(mmax,mmax)=BH(mmax,mmax);
  //   for (int i=mmin;i<=mmax;i++)
  //   {
  //     for (int j=mmin;j<=mmax;j++)
  //     {
  //       M(i,j)-=w(i)*w(j);
  //     }
  //   }
  //  /*
  //   ofstream ofs("newhess");
  //   ofs << M.indexmin() << " " << M.indexmax()  << endl;
  //   ofs << M << endl;
  //   ad_exit(1);
  //  */
  //   return inv(M);
  // }
  // extern "C" void adfloat_except(int k);
  // 
  // /*
  // banded_lower_triangular_dmatrix get_ltdchinv(MY_DOUBLE_TYPE rho,
  //   MY_DOUBLE_TYPE cvar,int mmin,int mmax)
  // {
  //   banded_lower_triangular_dmatrix tmp(mmin,mmax,2);
  //   MY_DOUBLE_TYPE sinv=1.0/sqrt(cvar);
  //   MY_DOUBLE_TYPE rho1=sqrt(1-rho*rho);
  //   tmp(mmin,mmin)=sinv;
  //   tmp(mmin+1,mmin)=-sinv*rho/rho1;
  //   for (int i=mmin+1;i<mmax;i++)
  //   {
  //     tmp(i,i)=sinv/rho1;
  //     tmp(i+1,i)=-sinv*rho/rho1;
  //   }
  //   tmp(mmax,mmax)=sinv/rho1;
  //   return tmp;
  // }
  // */
   static banded_lower_triangular_dmatrix get_ltdchinv(MY_DOUBLE_TYPE rho,
     MY_DOUBLE_TYPE cvar,int mmin,int mmax)
   {
     banded_lower_triangular_dmatrix tmp(mmin,mmax,2);
     MY_DOUBLE_TYPE rho1=sqrt(1-rho*rho);
     MY_DOUBLE_TYPE rho1i=1.0/rho1;
     MY_DOUBLE_TYPE sinv=1.0/sqrt(cvar);
     MY_DOUBLE_TYPE tt=rho1i*sinv;
     MY_DOUBLE_TYPE tt1=-rho*rho1i*sinv;
     tmp(1,1)=tt;
     for (int i=mmin+1;i<=mmax;i++)
     {
       tmp(i,i)=tt;
       tmp(i,i-1)=tt1;
     }
     tmp(mmax,mmax)=sinv;
     return tmp;
   }
  // 
   static void xdf_lognormal_ss_multinomial_laplace_approximation_1(void);
     symmetric_tridiagonal_dmatrix xget_tridag_BH(MY_DOUBLE_TYPE crho,
     MY_DOUBLE_TYPE cvar,int mmin,int mmax);
  // 
  // dvariable lognormal_multinomial_log_likelihood
  //   (const dvariable& vN,const dvector& q,const dvar_vector& vp,
  //   const dvariable&  _vrho,const dvariable&  _vvar,int ic)
  // {
  //   ADUNCONST(dvariable,vrho)
  //   ADUNCONST(dvariable,vvar)
  // 
  //   MY_DOUBLE_TYPE crho=value(vrho);
  //   MY_DOUBLE_TYPE cvar=value(vvar);
  // 
  //   MY_DOUBLE_TYPE N=value(vN);
  //   const dvector p=value(vp);
  //   int mmin=p.indexmin();
  //   int mmax=p.indexmax();
  // 
  //   symmetric_tridiagonal_dmatrix stdsinv=xget_tridag_BH(crho,cvar,mmin,
  //      mmax);
  // 
  //   banded_lower_triangular_dmatrix chinv=
  //     get_ltdchinv(crho,cvar,mmin,mmax);
  // 
  //   dvector eta_hat=lognormal_ss_multinomial_newton_raphson(N,q,p,chinv,stdsinv);
  //   
  //   MY_DOUBLE_TYPE ln_det=lognormal_ss_multinomial_laplace_approximation
  //     (N,q,p,eta_hat,stdsinv);
  //   
  // 
  //   dvariable vlog_det=nograd_assign(ln_det);
  //   save_identifier_string("t5");
  //   //vSinv.save_dvar_matrix_position();
  //   save_identifier_string("t4");
  //   eta_hat.save_dvector_value();
  //   save_identifier_string("t3");
  //   eta_hat.save_dvector_position();
  //   save_identifier_string("ps");
  //   vlog_det.save_prevariable_position();
  //   save_identifier_string("rt");
  //   vN.save_prevariable_value();
  //   save_identifier_string("s7");
  //   vN.save_prevariable_position();
  //   save_identifier_string("s6");
  //   q.save_dvector_value();
  //   save_identifier_string("s5");
  //   q.save_dvector_position();
  //   save_identifier_string("s4");
  //   vp.save_dvar_vector_value();
  //   save_identifier_string("s3");
  //   vp.save_dvar_vector_position();
  //   //save_identifier_string("s2");
  //   //save_pointer_value((void*)(&Sinv));
  //   save_identifier_string("s1");
  //   vrho.save_prevariable_value();
  //   save_identifier_string("s10");
  //   vvar.save_prevariable_value();
  //   save_identifier_string("s11");
  //   vrho.save_prevariable_position();
  //   save_identifier_string("s12");
  //   vvar.save_prevariable_position();
  //   save_identifier_string("s13");
  //   save_int_value(ic);
  //   save_identifier_string("s14");
  //   
  //   gradient_structure::GRAD_STACK1->
  //       set_gradient_stack(df_lognormal_ss_multinomial_laplace_approximation_1);
  //   return vlog_det;
  // }
  // 
dvariable lognormal_multinomial_log_likelihood
  (const dvariable& vN,const dvector& q,const dvar_vector& vp,
  const dvariable&  _vrho,const dvariable&  _vvar,const dvector& _eps,int ic)
{
  ADUNCONST(dvariable,vrho)
  ADUNCONST(dvariable,vvar)
  ADUNCONST(dvector,eps)

  MY_DOUBLE_TYPE crho=value(vrho);
  MY_DOUBLE_TYPE cvar=value(vvar);

  MY_DOUBLE_TYPE N=value(vN);
  const dvector p=value(vp);
  int mmin=p.indexmin();
  int mmax=p.indexmax();

  // XXXXXXXXXXXXX
  //crho=0.0;
  // XXXXXXXXXXXXX
  symmetric_tridiagonal_dmatrix stdsinv=xget_tridag_BH(crho,cvar,mmin,
     mmax);

  banded_lower_triangular_dmatrix chinv=
    get_ltdchinv(crho,cvar,mmin,mmax);

  int pps=0;

  if (pps)
  {
    dmatrix X=make_dmatrix(get_ltdchinv(0.9,1,1,6));
    cout << setfixed() << setprecision(3) << setw(6) << endl;
    cout  << endl << X << endl;
    cout << setfixed() << setprecision(3) << setw(6) << endl;
    cout  << endl << "inv(X*trans(X))"  << endl;
    cout  <<  inv(X*trans(X))  << endl;
    dmatrix Y=make_dmatrix(xget_tridag_BH(0.9,1,1,6));
    cout  << endl << "Y"  << endl;
    cout  <<  Y  << endl;
    ad_exit(1);
   /*
    int sgn;
    MY_DOUBLE_TYPE rho=0.8;
    MY_DOUBLE_TYPE v=1.0;
    MY_DOUBLE_TYPE rho1=sqrt(1.0-rho*rho);
    banded_lower_triangular_dmatrix xchinv=get_ltdchinv1(rho,v,1,8);
    dmatrix XM=make_dmatrix(xchinv);
    banded_lower_triangular_dmatrix xchinv1=get_ltdchinv1(rho,v,1,8);
    dmatrix XM1=make_dmatrix(xchinv1);
    ofstream ofs("ccc");
    ofs << "XM" << endl;
    ofs << setprecision(3) << setw(6) << setfixed() << endl;
    ofs << XM << endl;

    ofs << "XM1" << endl;
    ofs << setprecision(3) << setw(6) << setfixed() << endl;
    ofs << XM1 << endl;

    ofs << endl << "choleski_decomp(trans(XM)*XM)" << endl;
    ofs << setprecision(3) << setw(6) << setfixed() << endl;
    ofs << endl << choleski_decomp(trans(XM)*XM) << endl;

    ofs << endl << "choleski_decomp(XM1*trans(XM1))" << endl;
    ofs << setprecision(3) << setw(6) << setfixed() << endl;
    ofs << endl << choleski_decomp(XM1*trans(XM1)) << endl;
    ofs << "choleski_decomp(inv(XM)*trans(inv(XM)))" << endl;
    ofs << setprecision(3) << setw(6) << setfixed() << endl;
    ofs << choleski_decomp(inv(XM)*trans(inv(XM))) << endl;
    ofs << "inv(XM)*trans(inv(XM))" << endl;
    ofs << setprecision(3) << setw(6) << setfixed() << endl;
    ofs << inv(XM)*trans(inv(XM)) << endl;
    ofs << "trans(inv(XM))*inv(XM)" << endl;
    ofs << setprecision(3) << setw(6) << setfixed() << endl;
    ofs << trans(inv(XM))*inv(XM) << endl;
    ofs << "inv(XM)" << endl;
    ofs << setprecision(3) << setw(6) << setfixed() << endl;
    ofs << inv(XM) << endl;
    ofs << endl << "choleski_decomp(trans(XM)*XM)" << endl;
    ofs << setprecision(3) << setw(6) << setfixed() << endl;
    ofs << endl << choleski_decomp(trans(XM)*XM) << endl;
    ofs << endl << "rho1 " << rho1 << endl;
    ofs << endl << "rho1*choleski_decomp(trans(XM)*XM)" << endl;
    ofs << setprecision(3) << setw(6) << setfixed() << endl;
    ofs << endl << rho1*choleski_decomp(trans(XM)*XM) << endl;
    ofs << endl << "trans(XM)*XM" << endl;
    ofs << setprecision(3) << setw(6) << setfixed() << endl;
    ofs << endl << trans(XM)*XM << endl;
    ofs << endl << "inv(trans(XM)*XM)" << endl;
    ofs << setprecision(3) << setw(6) << setfixed() << endl;
    ofs << endl << inv(trans(XM)*XM) << endl;
    dmatrix YM=choleski_decomp(trans(XM)*XM);
    ofs << endl << "YM*trans(YM)" << endl;
    ofs << setprecision(3) << setw(6) << setfixed() << endl;
    ofs << endl << YM*trans(YM) << endl;
    ofs << endl << "trans(YM)*YM" << endl;
    ofs << setprecision(3) << setw(6) << setfixed() << endl;
    ofs << endl << trans(YM)*YM << endl;
    ofs << endl << "inv(trans(YM)*YM)" << endl;
    ofs << setprecision(3) << setw(6) << setfixed() << endl;
    ofs << endl << inv(trans(YM)*YM) << endl;
    ofs << endl << "XM1" << endl;
    ofs << setprecision(3) << setw(6) << setfixed() << endl;
    ofs << endl << XM1 << endl;
    ofs << endl << "XM1*trans(XM1)" << endl;
    ofs << setprecision(3) << setw(6) << setfixed() << endl;
    ofs << endl << XM1*trans(XM1) << endl;
    ofs << endl << "trans(XM1)*XM1" << endl;
    ofs << setprecision(3) << setw(6) << setfixed() << endl;
    ofs << endl << trans(XM1)*XM1 << endl;
    ofs << "NN " << norm2(choleski_decomp(trans(XM)*XM)-XM1) << endl;
  */
   

    ad_exit(1);
    
   /* 
#if !defined(NO_MY_DOUBLE_TYPE)
    dmatrix M1=make_ar1matrix(2,0.8,1.0L);
#else
    dmatrix M1=make_ar1matrix(2,0.8,1.0);
#endif
    cout << M1 << endl;
    cout << exp(ln_det(M1,sgn)) << endl;
#if !defined(NO_MY_DOUBLE_TYPE)
    dmatrix M2=make_ar1matrix(3,0.8,1.0L);
#else
    dmatrix M2=make_ar1matrix(3,0.8,1.0);
#endif
    cout << M2 << endl;
    cout << exp(ln_det(M2,sgn)) << endl;
#if !defined(NO_MY_DOUBLE_TYPE)
    dmatrix M3=make_ar1matrix(4,0.8,1.0L);
#else
    dmatrix M3=make_ar1matrix(4,0.8,1.0);
#endif
    cout << M3 << endl;
    cout << exp(ln_det(M3,sgn)) << endl;
#if !defined(NO_MY_DOUBLE_TYPE)
    dmatrix M4=make_ar1matrix(5,0.8,1.0L);
#else
    dmatrix M4=make_ar1matrix(5,0.8,1.0);
#endif
    cout << M4 << endl;
    cout << exp(ln_det(M4,sgn)) << endl;
    ad_exit(1);
    
    
    dmatrix M=make_ar1matrix(5,0.8,2);
    cout << endl << "A" << endl;
    cout << setprecision(4) << setfixed() << inv(M) << endl;
    symmetric_tridiagonal_dmatrix std=xget_tridag_BH(0.8,2.0,1,
      5);
    dmatrix SM=make_dmatrix(std);
    cout << endl << "B" << endl;
    cout << SM << endl;
    banded_lower_triangular_dmatrix ch1=get_ltdchinv(0.8,2.0,1,5);
    dmatrix cm=make_dmatrix(ch1);
    cout << endl << "C" << endl;
    cout << cm*trans(cm) << endl;
    ad_exit(1);
   */

    //dmatrix LDM=make_dmatrix(chinv);
    //cout << LDM << endl;
  }

  dvector eta_hat=lognormal_ss_multinomial_newton_raphson(N,q,p,chinv,stdsinv,
    eps);
  
  eps=eta_hat;

  MY_DOUBLE_TYPE ln_det=xlognormal_ss_multinomial_laplace_approximation
    (N,q,p,eta_hat,stdsinv);
  

  dvariable vlog_det=nograd_assign(ln_det);
//    save_identifier_string("t5");
  const char * str18;
  str18="t5";
  char* strx18=const_cast <char*> (str18);
  save_identifier_string(strx18);
  //vSinv.save_dvar_matrix_position();
//    save_identifier_string("t4");
  const char * str19;
  str19="t4";
  char* strx19=const_cast <char*> (str19);
  save_identifier_string(strx19);
  eta_hat.save_dvector_value();
//    save_identifier_string("t3");
  const char * str20;
  str20="t3";
  char* strx20=const_cast <char*> (str20);
  save_identifier_string(strx20);
  eta_hat.save_dvector_position();
//    save_identifier_string("ps");
  const char * str21;
  str21="ps";
  char* strx21=const_cast <char*> (str21);
  save_identifier_string(strx21);
  vlog_det.save_prevariable_position();
//    save_identifier_string("rt");
  const char * str22;
  str22="rt";
  char* strx22=const_cast <char*> (str22);
  save_identifier_string(strx22);
  vN.save_prevariable_value();
//    save_identifier_string("s7");
  const char * str23;
  str23="s7";
  char* strx23=const_cast <char*> (str23);
  save_identifier_string(strx23);
  vN.save_prevariable_position();
//    save_identifier_string("s6");
  const char * str24;
  str24="s6";
  char* strx24=const_cast <char*> (str24);
  save_identifier_string(strx24);
  q.save_dvector_value();
//    save_identifier_string("s5");
  const char * str25;
  str25="s5";
  char* strx25=const_cast <char*> (str25);
  save_identifier_string(strx25);
  q.save_dvector_position();
//    save_identifier_string("s4");
  const char * str26;
  str26="s4";
  char* strx26=const_cast <char*> (str26);
  save_identifier_string(strx26);
  vp.save_dvar_vector_value();
//    save_identifier_string("s3");
  const char * str27;
  str27="s3";
  char* strx27=const_cast <char*> (str27);
  save_identifier_string(strx27);
  vp.save_dvar_vector_position();
  //save_identifier_string("s2");
  //save_pointer_value((void*)(&Sinv));
//    save_identifier_string("s1");
  const char * str29;
  str29="s1";
  char* strx29=const_cast <char*> (str29);
  save_identifier_string(strx29);
  vrho.save_prevariable_value();
//    save_identifier_string("s10");
  const char * str30;
  str30="s10";
  char* strx30=const_cast <char*> (str30);
  save_identifier_string(strx30);
  vvar.save_prevariable_value();
//    save_identifier_string("s11");
  const char * str31;
  str31="s11";
  char* strx31=const_cast <char*> (str31);
  save_identifier_string(strx31);
  vrho.save_prevariable_position();
//    save_identifier_string("s12");
  const char * str32;
  str32="s12";
  char* strx32=const_cast <char*> (str32);
  save_identifier_string(strx32);
  vvar.save_prevariable_position();
//    save_identifier_string("s13");
  const char * str33;
  str33="s13";
  char* strx33=const_cast <char*> (str33);
  save_identifier_string(strx33);
  save_int_value(ic);
//    save_identifier_string("s14");
  const char * str34;
  str34="s14";
  char* strx34=const_cast <char*> (str34);
  save_identifier_string(strx34);
  
  gradient_structure::GRAD_STACK1->
      set_gradient_stack(xdf_lognormal_ss_multinomial_laplace_approximation_1);
  return vlog_det;
}

  MY_DOUBLE_TYPE xlognormal_ss_multinomial_laplace_approximation(const MY_DOUBLE_TYPE& N,
    const dvector& q,const dvector& p,
    const dvector& eta,const symmetric_tridiagonal_dmatrix& stdsinv)
  {
    int i,j,k;
    int mmin=q.indexmin();
    int mmax=q.indexmax();
    
    dvector t(mmin,mmax);
    dvector Nq=N*q;
    dvector grad(mmin,mmax);
    dvector w(mmin,mmax);
    symmetric_tridiagonal_dmatrix BHess(mmin,mmax);
    symmetric_tridiagonal_dmatrix& BSinv=
      (symmetric_tridiagonal_dmatrix&) (stdsinv);
    BHess.initialize();
    MY_DOUBLE_TYPE v=0.0;
    for (int i=mmin;i<=mmax;i++)
    {
      t(i)=phi(p(i),eta(i));
      v+=t(i);
    }
  
    MY_DOUBLE_TYPE mll=-Nq*log(t)+N*log(v)+0.5*eta*(BSinv*eta); 
    int myflag=0;
    if (myflag)
    {
      int ierr=0;
      cout << 0.5*ln_det(make_dmatrix(BSinv),ierr) << endl;
    }
   /*
    dmatrix MM(mmin,mmax,mmin,mmax);
    MM.initialize();
    for (int i=mmin;i<mmax;i++)
    {
      MM(i,i)=BSinv(i,i);
      MM(i+1,i)=BSinv(i+1,i);
      MM(i,i+1)=BSinv(i,i+1);
    }
    MM(mmax,mmax)=BSinv(mmax,mmax);
    int sgn=0;
    cout << 0.5*ln_det(MM,sgn) << endl;
   */
    
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
    for (int i=mmin;i<mmax;i++)
    {
      BHess(i,i)+=BSinv(i,i);
      BHess(i+1,i)=BSinv(i+1,i);
    }
    BHess(mmax,mmax)+=BSinv(mmax,mmax);
  
    // choleski decomp to get log det
    MY_DOUBLE_TYPE log_lik=mll;
    int ierr=0;
    MY_DOUBLE_TYPE ld=0.0;
    {
      {
        //int & ierr = (int &) _ierr;
        int ierr=0;
        int mmin=BHess.indexmin();
        int mmax=BHess.indexmax();
        
        symmetric_tridiagonal_dmatrix L(mmin,mmax);
        L.initialize();
      
        int i,j,k;
        MY_DOUBLE_TYPE tmp;
        if (BHess(mmin,mmin)<=0)
        {
          if (ierr==0)
            cerr << "U Error matrix not positive definite in choleski_decomp"
              <<endl;
          ierr=1;
          return 0;
        }
        L(mmin,mmin)=sqrt(BHess(mmin,mmin));
        for (i=mmin;i<=mmin+1;i++)
        {
          L(i,mmin)=BHess(i,mmin)/L(mmin,mmin);
        }
      
        for (i=mmin+1;i<=mmax;i++)
        {
          if (i>2)
          {	
            tmp=BHess(i,i-1);
            L(i,i-1)=tmp/L(i-1,i-1);
          }
          tmp=BHess(i,i);
          if (i>1)	
            tmp-=L(i,i-1)*L(i,i-1);
          if (tmp<=0)
          {
            if (ierr==0)
              cerr << "Z Error matrix not positive definite in choleski_decomp"
                <<endl;
            ierr=1;
            return 0;
          }
          L(i,i)=sqrt(tmp);
        }
        for (int i=mmin;i<=mmax;i++)
        {
          ld+=log(L(i,i));
        }
      }
    }
  
    dvector vv(mmin,mmax);
    vv=sqrt(SUMPEN);
    //dvector hk=solve(BHess,grad);
    dvector bx=solve(BHess,w);
    //h3 = hk + (w*hk)*bx / ( 1 - w*bx);
    dvector v1=solve(BHess,vv);
    dvector u = v1 + (w*v1)*bx / ( 1 - w*bx);
  
    // rank 1 update for determinant
    //dvector x=solve(BHess,w);
    MY_DOUBLE_TYPE tmp=1.0-(bx*w);
    if (tmp<=0.0)
    {
      cout << tmp << endl;
    }
    MY_DOUBLE_TYPE fpen=0.0;
    //ln_det+=0.5*log(posfun(tmp,.01,fpen));
    MY_DOUBLE_TYPE ldd=ld;
    ld+=0.5*log(1.0-bx*w);
    ld+=0.5*log(1.0+u*vv);
    int pps2=0;
    if (pps2)
    {
      dmatrix H=make_dmatrix(BHess);
      int sgn=0;
      cout << "A " << 0.5*ln_det(H,sgn) << endl;
      H-=outer_prod(w,w);
      cout << "B " << 0.5*ln_det(H,sgn) << endl;
      cout << "C " << ldd << endl;
      cout << "D " << ld << endl;
      ad_exit(1);
    }
  
    log_lik=ld+mll;
    return log_lik;
  }
  
static void xdf_lognormal_ss_multinomial_laplace_approximation_1(void)
{
  int i,j,k;
  //double dflog_det=restore_prevariable_derivative();

  verify_identifier_string("s14");
  int ic=restore_int_value();
  verify_identifier_string("s13");
  prevariable_position varpos=restore_prevariable_position();
  verify_identifier_string("s12");
  prevariable_position rhopos=restore_prevariable_position();
  verify_identifier_string("s11");
  MY_DOUBLE_TYPE cvar=restore_prevariable_value();
  verify_identifier_string("s10");
  MY_DOUBLE_TYPE crho=restore_prevariable_value();

  verify_identifier_string("s1");
  //dmatrix& Sinv=*(dmatrix*)(restore_pointer_value());
  //verify_identifier_string("s2");
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
  //dvar_matrix_position Sinvpos=restore_dvar_matrix_position();
  verify_identifier_string("t5");
  int mmin=q.indexmin();
  int mmax=q.indexmax();
  
  dvector dfeta(mmin,mmax);
  dfeta.initialize();
  dvector t(mmin,mmax);
  dvector dft(mmin,mmax);
  dft.initialize();
  symmetric_tridiagonal_dmatrix BHess(mmin,mmax);
  BHess.initialize();
  symmetric_tridiagonal_dmatrix BSinv=xget_tridag_BH(crho,cvar,
    mmin,mmax);
  symmetric_tridiagonal_dmatrix dfBHess(mmin,mmax);
  dfBHess.initialize();
  MY_DOUBLE_TYPE dfN=0;
  dvector dfp(mmin,mmax);
  dfp.initialize();

  dvector Nq=N*q;
  dvector dfNq(mmin,mmax);
  dfNq.initialize();
  MY_DOUBLE_TYPE v=0.0;
  MY_DOUBLE_TYPE dfv=0.0;
  for (i=mmin;i<=mmax;i++)
  {
    t(i)=phi(p(i),eta(i));
    //t(i)=p(i)*exp(eta(i));
    v+=t(i);
  }
  MY_DOUBLE_TYPE mll=-Nq*log(t)+N*v+0.5*eta*(BSinv*eta);
  // !!!!!!!!!GGGG
  MY_DOUBLE_TYPE dfmll=0.0;
  MY_DOUBLE_TYPE dfnv=0.0;
  MY_DOUBLE_TYPE dfnvv=0.0;
  symmetric_tridiagonal_dmatrix dfBSinv(mmin,mmax);
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
  for (int i=mmin;i<mmax;i++)
  {
    BHess(i,i)+=BSinv(i,i);
    BHess(i+1,i)=BSinv(i+1,i);
  }
  BHess(mmax,mmax)+=BSinv(mmax,mmax);

  int n=mmax;

  dvector dfw(mmin,mmax);
  dvector dfx(mmin,mmax);
  banded_lower_triangular_dmatrix L1(1,n,2);
  banded_lower_triangular_dmatrix dfL1(1,n,2);
  dvector tmp(1,n);
  dmatrix tmp1(1,n,1,n);
  dmatrix dftmp1(1,n,1,n);
  dvector dftmp(1,n);
  tmp.initialize();
  tmp1.initialize();
  dftmp.initialize();
  dftmp1.initialize();
  dfL1.initialize();

  //if (global_except_switch)
  //{
  //  adfloat_except(1);
  //}
  MY_DOUBLE_TYPE log_det1=mll;
  {
    if (BHess(1,1)<=0)
    {
      cerr << "X Error matrix not positive definite in choleski_decomp"
        <<endl;
      //adfloat_except(1);
      ad_exit(1);
    }
    L1(1,1)=sqrt(BHess(1,1));
    for (i=2;i<=2;i++)
    {
      L1(i,1)=BHess(i,1)/L1(1,1);
    }
  
    for (i=2;i<=n;i++)
    {
      for (j=i-1;j<=i-1;j++)
      {
        tmp1(i,j)=BHess(i,j);
        L1(i,j)=tmp1(i,j)/L1(j,j);
      }
      tmp(i)=BHess(i,i);
      tmp(i)-=L1(i,i-1)*L1(i,i-1);
      if (tmp(i)<=0)
      {
        cerr << "Y Error matrix not positive definite in choleski_decomp"
          <<endl;
        //adfloat_except(1);
        ad_exit(1);
      }
      L1(i,i)=sqrt(tmp(i));
    }
    //cout << norm(L-L1) << endl;
    for (i=1;i<=n;i++)
    {
      // EEEEEEE
      log_det1+=log(L1(i,i));
    }
  }
  dvector x=yxlt1solve_solvet(L1,w);
  dfx.initialize();
  dfw.initialize();
  MY_DOUBLE_TYPE uu=1.0-x*w;
  MY_DOUBLE_TYPE fpen=0.0;
  //double puu=posfun(uu,.01,fpen);
  //double dfpuu=0.0;
  MY_DOUBLE_TYPE dfuu=0.0;
  // XXXXXXXXXXXXx
  //log_det1+=0.5*log(puu);
  log_det1+=0.5*log(uu);
  MY_DOUBLE_TYPE log_det=log_det1;
 //*******************************************************************8
  {
    //double log_det=log_det1;
    MY_DOUBLE_TYPE dflog_det1=dflog_det;
    dflog_det=0.0;
    //log_det1+=0.5*log(puu);
  // XXXXXXXXXXXXx
    dfuu+=dflog_det1*0.5/uu;
    //double puu=posfun(uu,.01);
    //dfuu+=dfpuu*dfposfun(uu,.01);
    //dfpuu=0.0;
    //double uu=1.0-x*w;
    dfx-=dfuu*w;
    dfw-=dfuu*x;
    dfuu=0.0;
    //dvector x=lt1solve_solvet(L1,w);
    yxdflt1solve_solvet(L1,w,dfw,dfx,dfL1);
    
    for (i=1;i<=n;i++)
    {
      // EEEEEEE
      //log_det1+=log(L1(i,i));
      dfL1(i,i)+=dflog_det1/L1(i,i);
    }
    //double log_det1=mll;
    dfmll+=dflog_det1;
    dflog_det1=0.0;

  
    for (i=n;i>=2;i--)
    {
      //L1(i,i)=sqrt(tmp(i));
      dftmp(i)+=dfL1(i,i)/(2.0*L1(i,i));
      dfL1(i,i)=0.0;
      //tmp(i)-=L1(i,i-1)*L1(i,i-1);
      dfL1(i,i-1)-=2.*dftmp(i)*L1(i,i-1);
      //tmp(i)=BHess(i,i);
      dfBHess(i,i)+=dftmp(i);
      dftmp(i)=0.0;
      for (j=i-1;j>=i-1;j--)
      {
        //L1(i,j)=tmp1(i,j)/L1(j,j);
        MY_DOUBLE_TYPE linv=1./L1(j,j);
        dftmp1(i,j)+=dfL1(i,j)*linv;
        dfL1(j,j)-=dfL1(i,j)*tmp1(i,j)*linv*linv;
        dfL1(i,j)=0.0;
        //tmp1(i,j)=BHess(i,j);
        dfBHess(i,j)+=dftmp1(i,j);
        dftmp1(i,j)=0.0;
      }
    }
    MY_DOUBLE_TYPE xlinv=1./L1(1,1);
    MY_DOUBLE_TYPE xlinv2=xlinv*xlinv;
    for (i=2;i>=2;i--)
    {
      //L1(i,1)=BHess(i,1)/L1(1,1);
      dfBHess(i,1)+=dfL1(i,1)*xlinv;
      dfL1(1,1)-=dfL1(i,1)*BHess(i,1)*xlinv2;
      dfL1(i,1)=0.0;
    }
    //L1(1,1)=sqrt(BHess(1,1));
    dfBHess(1,1)+=dfL1(1,1)/(2.*L1(1,1));
  }
  //BHess(mmax,mmax)+=BSinv(mmax,mmax);
  dfBSinv(mmax,mmax)+=dfBHess(mmax,mmax);
  for (int i=mmin;i<mmax;i++)
  {
    //BHess(i+1,i)=BSinv(i+1,i);
    dfBSinv(i+1,i)+=dfBHess(i+1,i);
    dfBHess(i+1,i)=0.0;
    //BHess(i,i)+=BSinv(i,i);
    dfBSinv(i,i)+=dfBHess(i,i);
  }
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
  dfBHess.initialize();

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
  for (int i=mmin;i<mmax;i++)
  {
    dfBSinv(i,i)+=0.5*eta(i)*eta(i);
    dfBSinv(i+1,i)+=eta(i+1)*eta(i);
  }
  dfBSinv(mmax,mmax)+=0.5*eta(mmax)*eta(mmax);
  

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
  
  dvector cd=xget_cross_derivatives(L1,BHess,w,q,t,v,eta,N,dfeta,p);

  // add in the cross derivatives
 
  dfN+=cd(1);
  dfp+=cd(2,n+1).shift(1);
 
  for (int i=1;i<=n-1;i++)
  {
    //dfSinv(i,j)+=cd(1+n+(i-1)*n+j);
    dfBSinv(i,i)+=cd(1+n+(i-1)*n+i);
    dfBSinv(i+1,i)+=2.0*cd(1+n+i*n+i);
  }
  dfBSinv(n,n)+=cd(1+n+(n-1)*n+n);
 
  MY_DOUBLE_TYPE dfrho;
  MY_DOUBLE_TYPE dfvar;

  xdf_get_tridag_BH(crho,cvar,mmin,mmax,dfBSinv,dfrho,dfvar);
 
  //dfSinv.save_dmatrix_derivatives(Sinvpos);
  save_double_derivative(dfrho,rhopos);
  save_double_derivative(dfvar,varpos);
  save_double_derivative(dfN,Npos);
  dfp.save_dvector_derivatives(vppos);
}
  

symmetric_tridiagonal_dmatrix xget_tridag_BH(MY_DOUBLE_TYPE crho,
  MY_DOUBLE_TYPE cvar,int mmin,int mmax)
{
  symmetric_tridiagonal_dmatrix stsinv(mmin,mmax);
  MY_DOUBLE_TYPE cvinv=1.0/((1.0-crho*crho)*cvar);
  MY_DOUBLE_TYPE tt1=-cvinv*crho;
  MY_DOUBLE_TYPE tt2=(1+crho*crho)*cvinv;
  stsinv(1,1)=cvinv;
  for (int i=mmin+1;i<mmax;i++)
  {
    stsinv(i,i)=tt2;
  }
  stsinv(mmax,mmax)=cvinv;

  for (int i=mmin+1;i<=mmax;i++)
  {
    stsinv(i,i-1)=tt1;
  }
  return stsinv;
}
dvector xget_cross_derivatives(
  const banded_lower_triangular_dmatrix L1,
  const symmetric_tridiagonal_dmatrix& BH,
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
  dmatrix Minv=xinv(L1,w);
  //dmatrix uhatp=-Minv*NN;
  dmatrix uhatp=-yxlt1solve_solvet(L1,w,NN);
  tmp(1,n+1)=dfeta*uhatp;
 
 /*
  for (int i=1;i<=n;i++)
  {
    for (int j=1;j<=n;j++)
    {
      dvector t=Minv(i)*eta(j)+Minv(j)*eta(i);
      tmp(1+n+(i-1)*n+j)=-0.5*(dfeta*t);
    }
  }
 */
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
  
  return tmp;
}

dvector yxlt1solve_solvet(
  const banded_lower_triangular_dmatrix & L,
  const dvector& w)
{
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  dvector v(mmin,mmax);
  dvector x(mmin,mmax);
  dvector asum(mmin,mmax);
  dvector bsum(mmin,mmax);
  v.initialize();
  v(mmin)=w(mmin)/L(mmin,mmin);
  for (int i=mmin+1;i<=mmax;i++)
  {
    asum(i)=L(i,i-1)*v(i-1);
    v(i)=(w(i)-asum(i))/L(i,i);
  }
  //cout << " BB " << v << endl;
  x(mmax)=v(mmax)/L(mmax,mmax);
  for (int i=mmax-1;i>=mmin;i--)
  {
    bsum(i)=L(i+1,i)*x(i+1);
    x(i)=(v(i)-bsum(i))/L(i,i);
  }
  //cout << " ZZ " << x << endl;
  return x;
}

void yxdflt1solve_solvet(
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
  dvector v(mmin,mmax);
  dvector x(mmin,mmax);
  dvector asum(mmin,mmax);
  dvector bsum(mmin,mmax);
  v.initialize();
  v(mmin)=w(mmin)/L(mmin,mmin);
  for (int i=mmin+1;i<=mmax;i++)
  {
    asum(i)=L(i,i-1)*v(i-1);
    v(i)=(w(i)-asum(i))/L(i,i);
  }
  x(mmax)=v(mmax)/L(mmax,mmax);
  for (int i=mmax-1;i>=mmin;i--)
  {
    bsum(i)=L(i+1,i)*x(i+1);
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
    dfL(i+1,i)+=dfbsum(i)*x(i+1);
    dfx(i+1)+=dfbsum(i)*L(i+1,i);
    dfbsum(i)=0.0;
  }
  //x(mmax)=v(mmax)/L(mmax,mmax);
  dfv(mmax)+=dfx(mmax)/L(mmax,mmax);
  dfL(mmax,mmax)-=dfx(mmax)*x(mmax)/L(mmax,mmax);
  for (int i=mmax;i>=mmin+1;i--)
  {
    //v(i)=(w(i)-asum(i))/L(i,i);
    dfw(i)+=dfv(i)/L(i,i);
    dfasum(i)-=dfv(i)/L(i,i);
    dfL(i,i)-=dfv(i)*v(i)/L(i,i);
    dfv(i)=0.0;
    //asum(i)=L(i,i-1)*v(i-1);
    dfL(i,i-1)+=dfasum(i)*v(i-1);
    dfv(i-1)+=dfasum(i)*L(i,i-1);
    dfasum(i)=0.0;
  }
  //v(mmin)=w(mmin)/L(mmin,mmin);
  dfw(mmin)+=dfv(mmin)/L(mmin,mmin);
  dfL(mmin,mmin)-=dfv(mmin)*v(mmin)/L(mmin,mmin);
  dfv(mmin)=0.0;

}
dmatrix yxlt1solve_solvet(
  const banded_lower_triangular_dmatrix & L,
  const dvector& ww,
  const dmatrix& N)
{
  dmatrix TN=trans(N);
  int mmin=ww.indexmin();
  int mmax=ww.indexmax();
  int cmin=TN.indexmin();
  int cmax=TN.indexmax();
  dmatrix tmp(cmin,cmax,mmin,mmax);
  dvector v(mmin,mmax);
  dvector x(mmin,mmax);
  dvector asum(mmin,mmax);
  dvector bsum(mmin,mmax);
  dvector aw(mmin,mmax);

  const dvector& w=ww;
  v.initialize();
  asum.initialize();
  bsum.initialize();
  v(mmin)=w(mmin)/L(mmin,mmin);
  for (int i=mmin+1;i<=mmax;i++)
  {
    asum(i)=L(i,i-1)*v(i-1);
    v(i)=(w(i)-asum(i))/L(i,i);
  }
  //cout << " BB " << v << endl;
  x(mmax)=v(mmax)/L(mmax,mmax);
  for (int i=mmax-1;i>=mmin;i--)
  {
    bsum(i)=L(i+1,i)*x(i+1);
    x(i)=(v(i)-bsum(i))/L(i,i);
  }
  aw=x;
  for (int ii=cmin;ii<=cmax;ii++)
  {
    dvector& w=TN(ii);
    v.initialize();
    asum.initialize();
    bsum.initialize();
    
    v(mmin)=w(mmin)/L(mmin,mmin);
    for (int i=mmin+1;i<=mmax;i++)
    {
      asum(i)=L(i,i-1)*v(i-1);
      v(i)=(w(i)-asum(i))/L(i,i);
    }
    //cout << " BB " << v << endl;
    x(mmax)=v(mmax)/L(mmax,mmax);
    for (int i=mmax-1;i>=mmin;i--)
    {
      bsum(i)=L(i+1,i)*x(i+1);
      x(i)=(v(i)-bsum(i))/L(i,i);
    }
    tmp(ii)=x + (aw*(ww*x))/(1.0-(ww*aw));
              
  }
  return trans(tmp);
} 
dmatrix xinv(const banded_lower_triangular_dmatrix & L,
  const dvector& ww)
{
  int mmin=ww.indexmin();
  int mmax=ww.indexmax();
  dmatrix tmp(mmin,mmax,mmin,mmax);
  dvector v(mmin,mmax);
  dvector x(mmin,mmax);
  dvector asum(mmin,mmax);
  dvector bsum(mmin,mmax);
  dvector aw(mmin,mmax);

  const dvector& w=ww;
  v.initialize();
  asum.initialize();
  bsum.initialize();
  v(mmin)=w(mmin)/L(mmin,mmin);
  for (int i=mmin+1;i<=mmax;i++)
  {
    asum(i)=L(i,i-1)*v(i-1);
    v(i)=(w(i)-asum(i))/L(i,i);
  }
  //cout << " BB " << v << endl;
  x(mmax)=v(mmax)/L(mmax,mmax);
  for (int i=mmax-1;i>=mmin;i--)
  {
    bsum(i)=L(i+1,i)*x(i+1);
    x(i)=(v(i)-bsum(i))/L(i,i);
  }
  aw=x;
  dvector ee(mmin,mmax);
  ee.initialize();
  for (int ii=mmin;ii<=mmax;ii++)
  {
    if (ii>mmin) ee(ii-1)=0.0;
    ee(ii)=1.0;
    dvector& w=ee;
    v.initialize();
    asum.initialize();
    bsum.initialize();
    
    v(mmin)=w(mmin)/L(mmin,mmin);
    for (int i=mmin+1;i<=mmax;i++)
    {
      asum(i)=L(i,i-1)*v(i-1);
      v(i)=(w(i)-asum(i))/L(i,i);
    }
    //cout << " BB " << v << endl;
    x(mmax)=v(mmax)/L(mmax,mmax);
    for (int i=mmax-1;i>=mmin;i--)
    {
      bsum(i)=L(i+1,i)*x(i+1);
      x(i)=(v(i)-bsum(i))/L(i,i);
    }
    tmp(ii)=x + (aw*(ww*x))/(1.0-(ww*aw));
              
  }
  return trans(tmp);
} 
/*  old code
void  df_get_tridag_BH(MY_DOUBLE_TYPE crho,
  MY_DOUBLE_TYPE cvar,int mmin,int mmax,
  const symmetric_tridiagonal_dmatrix&  _dfstsinv,const MY_DOUBLE_TYPE& _dfcrho,
  const MY_DOUBLE_TYPE& _dfcvar)
{
  ADUNCONST(double,dfcrho)
  ADUNCONST(double,dfcvar)
  ADUNCONST(symmetric_tridiagonal_dmatrix,dfstsinv)
  symmetric_tridiagonal_dmatrix stsinv(mmin,mmax);
  MY_DOUBLE_TYPE ttm=(1.0-crho*crho)*cvar;
  //double tt0=1.0/(1.0-crho*crho)*cvar);
  //double tt1=1.0/sqrt((1.0-crho*crho)*cvar);
  MY_DOUBLE_TYPE tt0=1.0/ttm;
  MY_DOUBLE_TYPE tt1=1.0/sqrt(ttm);
  MY_DOUBLE_TYPE tt2=-crho*tt1;
  MY_DOUBLE_TYPE tt3=1.0/sqrt(cvar);
  MY_DOUBLE_TYPE tt4=-crho*tt0;
  MY_DOUBLE_TYPE tt5=1.0+crho*crho;
  MY_DOUBLE_TYPE tt6=tt0*tt5;
  stsinv(1,1)=tt0;
  for (int i=mmin+1;i<mmax;i++)
  {
    stsinv(i,i)=tt6;
  }
  stsinv(mmax,mmax)=tt0;

  for (int i=mmin+1;i<=mmax;i++)
  {
    stsinv(i,i-1)=tt4;
  }
  MY_DOUBLE_TYPE dftt0=0.0;
  MY_DOUBLE_TYPE dftt1=0.0;
  MY_DOUBLE_TYPE dftt2=0.0;
  MY_DOUBLE_TYPE dftt3=0.0;
  MY_DOUBLE_TYPE dftt4=0.0;
  MY_DOUBLE_TYPE dftt5=0.0;
  MY_DOUBLE_TYPE dftt6=0.0;
  MY_DOUBLE_TYPE dfttm=0.0;
  MY_DOUBLE_TYPE dfcrho1=0.0;
  dfcvar=0.0;
  dfcrho=0.0;
  for (int i=mmin+1;i<=mmax;i++)
  {
    //stsinv(i,i-1)=tt4;
    dftt4+=dfstsinv(i,i-1);
    dfstsinv(i,i-1)=0.0;
  }
  //stsinv(mmax,mmax)=tt0;
  dftt0+=dfstsinv(mmax,mmax);
  dfstsinv(mmax,mmax)=0.0;
  for (int i=mmin+1;i<mmax;i++)
  {
    //stsinv(i,i)=tt6;
    dftt6+=dfstsinv(i,i);
    dfstsinv(i,i)=0.0;
  }
  //stsinv(1,1)=tt0;
  dftt0+=dfstsinv(1,1);
  dfstsinv(1,1)=0.0;
  //double tt6=tt0*tt5;
  dftt0+=dftt6*tt5;
  dftt5+=dftt6*tt0;
  dftt6=0.0;
  //double tt5=1.0+crho*crho;
  dfcrho+=dftt5*2.0*crho;
  dftt5=0.0;
  //double tt4=-crho*tt0;
  dfcrho-=dftt4*tt0;
  dftt0-=dftt4*crho;
  dftt4=0.0;
  //double tt3=1.0/sqrt(cvar);
  dfcvar-=dftt3*tt3/cvar*1.5;
  dftt3=0.0;
  //double tt2=-crho*tt1;
  dftt1-=dftt2*crho;
  dfcrho-=dftt2*tt1;
  dftt2=0.0;
  //double tt1=1.0/sqrt(ttm);
  dfttm-=dftt1*1.5*tt1/ttm;
  dftt1=0.0;
  //double tt0=1.0/ttm
  dfttm-=dftt0*tt0/ttm;
  dftt0=0;
  //double ttm=(1.0-crho*crho)*cvar;
  dfcrho-=dfttm*2.0*crho*cvar;
  dfcvar+=dfttm*(1.0-crho*crho);
  dfttm=0.0;
}
*/
//  new code
static void  xdf_get_tridag_BH(MY_DOUBLE_TYPE crho,
  MY_DOUBLE_TYPE cvar,int mmin,int mmax,
  const symmetric_tridiagonal_dmatrix&  _dfstsinv,const MY_DOUBLE_TYPE& _dfcrho,
  const MY_DOUBLE_TYPE& _dfcvar)
{
  ADUNCONST(MY_DOUBLE_TYPE,dfcrho)
  ADUNCONST(MY_DOUBLE_TYPE,dfcvar)
  ADUNCONST(symmetric_tridiagonal_dmatrix,dfstsinv)
  symmetric_tridiagonal_dmatrix stsinv(mmin,mmax);
  MY_DOUBLE_TYPE crho1=1-crho*crho;
  MY_DOUBLE_TYPE cvinv=1/(crho1*cvar);
  MY_DOUBLE_TYPE tt1=-cvinv*crho;
  MY_DOUBLE_TYPE tt2=(1+crho*crho)*cvinv;
  stsinv(1,1)=cvinv;
  for (int i=mmin+1;i<mmax;i++)
  {
    stsinv(i,i)=tt2;
  }
  stsinv(mmax,mmax)=cvinv;

  for (int i=mmin+1;i<=mmax;i++)
  {
    stsinv(i,i-1)=tt1;
  }
  MY_DOUBLE_TYPE dfcvinv=0.0;
  MY_DOUBLE_TYPE dftt1=0.0;
  MY_DOUBLE_TYPE dftt2=0.0;
  MY_DOUBLE_TYPE dfcrho1=0.0;
  dfcvar=0.0;
  dfcrho=0.0;
  for (int i=mmin+1;i<=mmax;i++)
  {
    //stsinv(i,i-1)=tt1;
    dftt1+=dfstsinv(i,i-1);
  }
  //stsinv(mmax,mmax)=cvinv;
  dfcvinv+=dfstsinv(mmax,mmax);
  for (int i=mmin+1;i<mmax;i++)
  {
    //stsinv(i,i)=tt2;
    dftt2+=dfstsinv(i,i);
  }
  //stsinv(1,1)=cvinv;
  dfcvinv+=dfstsinv(1,1);
  //double tt2=(1+crho*crho)*cvinv;
  dfcrho+=dftt2*2*crho*cvinv;
  dfcvinv+=dftt2*(1+crho*crho);
  //double tt1=-cvinv*crho;
  dfcvinv-=dftt1*crho;
  dfcrho-=dftt1*cvinv;
  //double cvinv=1/(crho1*cvar);
  dfcvar-=dfcvinv/(crho1*square(cvar));
  dfcrho1-=dfcvinv/(square(crho1)*cvar);
  //double crho1=1-crho*crho;
  dfcrho-=2*dfcrho1*crho;
}
