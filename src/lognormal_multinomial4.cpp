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

  dvector lt1solve_solvet(const dmatrix& L,const dvector& y);
void dflt1solve_solvet(const dmatrix& L,const dvector& y,
  const dvector& dfx,const dvector& dfy,const dmatrix& dfL);
dvector lognormal_ss_multinomial_newton_raphson(const MY_DOUBLE_TYPE N,
  const dvector& q,const dvector& p,const dmatrix & Sinv);
dvector lognormal_ss_multinomial_newton_raphson(const MY_DOUBLE_TYPE N,
  const dvector& q,const dvector& p,const dmatrix & Sinv,const dmatrix & S);
MY_DOUBLE_TYPE lognormal_ss_multinomial_laplace_approximation(const MY_DOUBLE_TYPE& N,
  const dvector& q,const dvector& p,const dmatrix & Sinv,
  const dvector& eta);
MY_DOUBLE_TYPE lognormal_ss_multinomial_laplace_approximation(const MY_DOUBLE_TYPE& N,
  const dvector& q,const dvector& p,
 const dvector& eta,const symmetric_tridiagonal_dmatrix& stdsinv);
static void df_lognormal_ss_multinomial_laplace_approximation(void);
static void df_lognormal_ss_multinomial_laplace_approximation_1(void);
dvector lognormal_ss_multinomial_newton_raphson(const MY_DOUBLE_TYPE N,
  const dvector& q,const dvector& p,
  //const dmatrix & cch,
  const banded_lower_triangular_dmatrix & _chinv,
  const symmetric_tridiagonal_dmatrix& lensinv);

dvector get_cross_derivatives(const dmatrix& M,const dvector& q,
  const dvector& t,const MY_DOUBLE_TYPE v,const dvector& eta,const MY_DOUBLE_TYPE& N,
  const dvector dfeta,const dvector& p);
extern int new_flag;

MY_DOUBLE_TYPE phi(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE e);
MY_DOUBLE_TYPE dphi_deta(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE e);
MY_DOUBLE_TYPE d2phi_deta2(MY_DOUBLE_TYPE p,MY_DOUBLE_TYPE e);

int stupid_print_switch=0;
int wght_print_flag=0;

dvector ltsolve(const dmatrix& L,const dvector& eps)
{
  int mmin=eps.indexmin();
  int mmax=eps.indexmax();
  dvector v(mmin,mmax);
  dvector ssum(mmin,mmax);
  v.initialize();
  for (int i=mmin;i<=mmax;i++)
  {
    ssum(i)=0.0;
    for (int j=mmin;j<i;j++)
    {
      ssum(i)+=L(i,j)*v(j);
    }
    v(i)=(eps(i)-ssum(i))/L(i,i);
  }
  return v;
}
ofstream  stupid_print("ssmult.dat"); 


dmatrix importance_sampling(const dvector& Nq,const dvector& p,
  const dvector& eta_hat,const dmatrix& L,int nsamp,const dmatrix& Sinv,
  MY_DOUBLE_TYPE N)
{
  random_number_generator rng(9813);
  int mmin=eta_hat.indexmin();
  int mmax=eta_hat.indexmax();
  dmatrix fvalues(1,2,0,nsamp);
  for (int ii=0;ii<=nsamp;ii++)
  {
    dvector eps(mmin,mmax);
    if (ii==0)
    {
      eps.initialize();
    }
    else
    {
      eps.fill_randn(rng);
    }
    MY_DOUBLE_TYPE v=0.0;
    dvector chi=ltsolve(L,eps);
    //dvector nu=eta_hat+L*eps;
    dvector nu=eta_hat+solve(L,eps);
    dvector t(mmin,mmax);
    for (int i=mmin;i<=mmax;i++)
    {
      t(i)=phi(p(i),nu(i));
      v+=t(i);
    }
    fvalues(2,ii)=-(Nq*log(t))+N*log(v)+0.5*nu*(Sinv*nu)-0.5*norm2(eps); 
    fvalues(1,ii)=0.5*norm2(eps); 
  }
  fvalues(2)-=fvalues(2,0); 
  return fvalues;
}

/*
dvariable lognormal_multinomial_log_likelihood
  (const dvariable& vN,const dvector& q,const dvar_vector& vp,
  const dvar_matrix& vSinv,const dmatrix& Sinv)
{

  MY_DOUBLE_TYPE N=value(vN);
  const dvector p=value(vp);

  dvector eta_hat=lognormal_ss_multinomial_newton_raphson(N,q,p,Sinv);
  
  MY_DOUBLE_TYPE ln_det=lognormal_ss_multinomial_laplace_approximation
    (N,q,p,Sinv,eta_hat);
  
  int mmin=p.indexmin();
  int mmax=p.indexmax();

  dvector tt(mmin,mmax);
  MY_DOUBLE_TYPE v=0.0;
  for (int i=mmin;i<=mmax;i++)
  {
    tt(i)=phi(p(i),eta_hat(i));
    v+=tt(i);
  }

  dvariable vlog_det=nograd_assign(ln_det);
//    save_identifier_string("t5");
  const char * str1;
  str1="t5";
  char* strx1=const_cast <char*> (str1);
  save_identifier_string(strx1);
  vSinv.save_dvar_matrix_position();
//    save_identifier_string("t4");
  const char * str2;
  str2="t4";
  char* strx2=const_cast <char*> (str2);
  save_identifier_string(strx2);
  eta_hat.save_dvector_value();
//    save_identifier_string("t3");
  const char * str3;
  str3="t3";
  char* strx3=const_cast <char*> (str3);
  save_identifier_string(strx3);
  eta_hat.save_dvector_position();
//    save_identifier_string("ps");
  const char * str4;
  str4="ps";
  char* strx4=const_cast <char*> (str4);
  save_identifier_string(strx4);
  vlog_det.save_prevariable_position();
//    save_identifier_string("rt");
  const char * str5;
  str5="rt";
  char* strx5=const_cast <char*> (str5);
  save_identifier_string(strx5);
  vN.save_prevariable_value();
//    save_identifier_string("s7");
  const char * str6;
  str6="s7";
  char* strx6=const_cast <char*> (str6);
  save_identifier_string(strx6);
  vN.save_prevariable_position();
//    save_identifier_string("s6");
  const char * str7;
  str7="s6";
  char* strx7=const_cast <char*> (str7);
  save_identifier_string(strx7);
  q.save_dvector_value();
//    save_identifier_string("s5");
  const char * str8;
  str8="s5";
  char* strx8=const_cast <char*> (str8);
  save_identifier_string(strx8);
  q.save_dvector_position();
//    save_identifier_string("s4");
  const char * str9;
  str9="s4";
  char* strx9=const_cast <char*> (str9);
  save_identifier_string(strx9);
  vp.save_dvar_vector_value();
//    save_identifier_string("s3");
  const char * str10;
  str10="s3";
  char* strx10=const_cast <char*> (str10);
  save_identifier_string(strx10);
  vp.save_dvar_vector_position();
//    save_identifier_string("s2");
  const char * str11;
  str11="s2";
  char* strx11=const_cast <char*> (str11);
  save_identifier_string(strx11);
  save_pointer_value((void*)(&Sinv));
//    save_identifier_string("s1");
  const char * str12;
  str12="s1";
  char* strx12=const_cast <char*> (str12);
  save_identifier_string(strx12);
  
  gradient_structure::GRAD_STACK1->
      set_gradient_stack(df_lognormal_ss_multinomial_laplace_approximation);
  return vlog_det;
}
*/

// *********************************************************************
// *********************************************************************
/*
dvariable lognormal_multinomial_log_likelihood
  (const dvariable& vN,const dvector& q,const dvar_vector& vp,
  const dvar_matrix& vSinv,const dmatrix& Sinv,const dmatrix & ch,
  const banded_lower_triangular_dmatrix & chinv,
  const symmetric_tridiagonal_dmatrix& lensinv)
{

  MY_DOUBLE_TYPE N=value(vN);
  const dvector p=value(vp);

  dvector eta_hat=lognormal_ss_multinomial_newton_raphson(N,q,p,Sinv,ch,chinv,
    lensinv);
  
  MY_DOUBLE_TYPE ln_det=lognormal_ss_multinomial_laplace_approximation
    (N,q,p,Sinv,eta_hat);
  
  int mmin=p.indexmin();
  int mmax=p.indexmax();

  dvector tt(mmin,mmax);
  MY_DOUBLE_TYPE v=0.0;
  for (int i=mmin;i<=mmax;i++)
  {
    tt(i)=phi(p(i),eta_hat(i));
    v+=tt(i);
  }

  dvariable vlog_det=nograd_assign(ln_det);
//    save_identifier_string("t5");
  const char * str13;
  str13="t5";
  char* strx13=const_cast <char*> (str13);
  save_identifier_string(strx13);
  vSinv.save_dvar_matrix_position();
//    save_identifier_string("t4");
  const char * str14;
  str14="t4";
  char* strx14=const_cast <char*> (str14);
  save_identifier_string(strx14);
  eta_hat.save_dvector_value();
//    save_identifier_string("t3");
  const char * str15;
  str15="t3";
  char* strx15=const_cast <char*> (str15);
  save_identifier_string(strx15);
  eta_hat.save_dvector_position();
//    save_identifier_string("ps");
  const char * str16;
  str16="ps";
  char* strx16=const_cast <char*> (str16);
  save_identifier_string(strx16);
  vlog_det.save_prevariable_position();
//    save_identifier_string("rt");
  const char * str17;
  str17="rt";
  char* strx17=const_cast <char*> (str17);
  save_identifier_string(strx17);
  vN.save_prevariable_value();
//    save_identifier_string("s7");
  const char * str18;
  str18="s7";
  char* strx18=const_cast <char*> (str18);
  save_identifier_string(strx18);
  vN.save_prevariable_position();
//    save_identifier_string("s6");
  const char * str19;
  str19="s6";
  char* strx19=const_cast <char*> (str19);
  save_identifier_string(strx19);
  q.save_dvector_value();
//    save_identifier_string("s5");
  const char * str20;
  str20="s5";
  char* strx20=const_cast <char*> (str20);
  save_identifier_string(strx20);
  q.save_dvector_position();
//    save_identifier_string("s4");
  const char * str21;
  str21="s4";
  char* strx21=const_cast <char*> (str21);
  save_identifier_string(strx21);
  vp.save_dvar_vector_value();
//    save_identifier_string("s3");
  const char * str22;
  str22="s3";
  char* strx22=const_cast <char*> (str22);
  save_identifier_string(strx22);
  vp.save_dvar_vector_position();
//    save_identifier_string("s2");
  const char * str23;
  str23="s2";
  char* strx23=const_cast <char*> (str23);
  save_identifier_string(strx23);
  save_pointer_value((void*)(&Sinv));
//    save_identifier_string("s1");
  const char * str24;
  str24="s1";
  char* strx24=const_cast <char*> (str24);
  save_identifier_string(strx24);
  
  gradient_structure::GRAD_STACK1->
      set_gradient_stack(df_lognormal_ss_multinomial_laplace_approximation);
  return vlog_det;
}
*/
// *********************************************************************
// *********************************************************************
// *********************************************************************
dvariable lognormal_multinomial_log_likelihood
  (const dvariable& vN,const dvector& q,const dvar_vector& vp,
  const dvar_matrix& vSinv,const dmatrix& Sinv,const dmatrix & ch,
  const banded_lower_triangular_dmatrix & chinv,
  const symmetric_tridiagonal_dmatrix& stdsinv)
{

  MY_DOUBLE_TYPE N=value(vN);
  const dvector p=value(vp);

  dvector eta_hat=lognormal_ss_multinomial_newton_raphson(N,q,p,
    //ch,
    chinv,
    stdsinv);
  
  MY_DOUBLE_TYPE ln_det=lognormal_ss_multinomial_laplace_approximation
    (N,q,p,eta_hat,stdsinv);
  
  int mmin=p.indexmin();
  int mmax=p.indexmax();

  dvector tt(mmin,mmax);
  MY_DOUBLE_TYPE v=0.0;
  for (int i=mmin;i<=mmax;i++)
  {
    tt(i)=phi(p(i),eta_hat(i));
    v+=tt(i);
  }

  dvariable vlog_det=nograd_assign(ln_det);
//    save_identifier_string("t5");
  const char * str25;
  str25="t5";
  char* strx25=const_cast <char*> (str25);
  save_identifier_string(strx25);
  vSinv.save_dvar_matrix_position();
//    save_identifier_string("t4");
  const char * str26;
  str26="t4";
  char* strx26=const_cast <char*> (str26);
  save_identifier_string(strx26);
  eta_hat.save_dvector_value();
//    save_identifier_string("t3");
  const char * str27;
  str27="t3";
  char* strx27=const_cast <char*> (str27);
  save_identifier_string(strx27);
  eta_hat.save_dvector_position();
//    save_identifier_string("ps");
  const char * str28;
  str28="ps";
  char* strx28=const_cast <char*> (str28);
  save_identifier_string(strx28);
  vlog_det.save_prevariable_position();
//    save_identifier_string("rt");
  const char * str29;
  str29="rt";
  char* strx29=const_cast <char*> (str29);
  save_identifier_string(strx29);
  vN.save_prevariable_value();
//    save_identifier_string("s7");
  const char * str30;
  str30="s7";
  char* strx30=const_cast <char*> (str30);
  save_identifier_string(strx30);
  vN.save_prevariable_position();
//    save_identifier_string("s6");
  const char * str31;
  str31="s6";
  char* strx31=const_cast <char*> (str31);
  save_identifier_string(strx31);
  q.save_dvector_value();
//    save_identifier_string("s5");
  const char * str32;
  str32="s5";
  char* strx32=const_cast <char*> (str32);
  save_identifier_string(strx32);
  q.save_dvector_position();
//    save_identifier_string("s4");
  const char * str33;
  str33="s4";
  char* strx33=const_cast <char*> (str33);
  save_identifier_string(strx33);
  vp.save_dvar_vector_value();
//    save_identifier_string("s3");
  const char * str34;
  str34="s3";
  char* strx34=const_cast <char*> (str34);
  save_identifier_string(strx34);
  vp.save_dvar_vector_position();
//    save_identifier_string("s2");
  const char * str35;
  str35="s2";
  char* strx35=const_cast <char*> (str35);
  save_identifier_string(strx35);
  save_pointer_value((void*)(&Sinv));
//    save_identifier_string("s1");
  const char * str36;
  str36="s1";
  char* strx36=const_cast <char*> (str36);
  save_identifier_string(strx36);
  
  gradient_structure::GRAD_STACK1->
      set_gradient_stack(df_lognormal_ss_multinomial_laplace_approximation);
  return vlog_det;
}
/*
dvariable lognormal_multinomial_log_likelihood
  (const dvariable& vN,const dvector& q,const dvar_vector& vp,
  const dvar_matrix& vSinv,const dmatrix& Sinv,const dmatrix & ch,
  const banded_lower_triangular_dmatrix & chinv,
  const symmetric_tridiagonal_dmatrix& stdsinv,
  MY_DOUBLE_TYPE crho,MY_DOUBLE_TYPE cvar)
{

  MY_DOUBLE_TYPE N=value(vN);
  const dvector p=value(vp);

  dvector eta_hat=lognormal_ss_multinomial_newton_raphson(N,q,p,Sinv,ch,chinv,
    stdsinv);
  
  MY_DOUBLE_TYPE ln_det=lognormal_ss_multinomial_laplace_approximation
    (N,q,p,Sinv,eta_hat,stdsinv);
  
  int mmin=p.indexmin();
  int mmax=p.indexmax();

  dvector tt(mmin,mmax);
  MY_DOUBLE_TYPE v=0.0;
  for (int i=mmin;i<=mmax;i++)
  {
    tt(i)=phi(p(i),eta_hat(i));
    v+=tt(i);
  }

  dvariable vlog_det=nograd_assign(ln_det);
//    save_identifier_string("t5");
  const char * str37;
  str37="t5";
  char* strx37=const_cast <char*> (str37);
  save_identifier_string(strx37);
  vSinv.save_dvar_matrix_position();
//    save_identifier_string("t4");
  const char * str38;
  str38="t4";
  char* strx38=const_cast <char*> (str38);
  save_identifier_string(strx38);
  eta_hat.save_dvector_value();
//    save_identifier_string("t3");
  const char * str39;
  str39="t3";
  char* strx39=const_cast <char*> (str39);
  save_identifier_string(strx39);
  eta_hat.save_dvector_position();
//    save_identifier_string("ps");
  const char * str40;
  str40="ps";
  char* strx40=const_cast <char*> (str40);
  save_identifier_string(strx40);
  vlog_det.save_prevariable_position();
//    save_identifier_string("rt");
  const char * str41;
  str41="rt";
  char* strx41=const_cast <char*> (str41);
  save_identifier_string(strx41);
  vN.save_prevariable_value();
//    save_identifier_string("s7");
  const char * str42;
  str42="s7";
  char* strx42=const_cast <char*> (str42);
  save_identifier_string(strx42);
  vN.save_prevariable_position();
//    save_identifier_string("s6");
  const char * str43;
  str43="s6";
  char* strx43=const_cast <char*> (str43);
  save_identifier_string(strx43);
  q.save_dvector_value();
//    save_identifier_string("s5");
  const char * str44;
  str44="s5";
  char* strx44=const_cast <char*> (str44);
  save_identifier_string(strx44);
  q.save_dvector_position();
//    save_identifier_string("s4");
  const char * str45;
  str45="s4";
  char* strx45=const_cast <char*> (str45);
  save_identifier_string(strx45);
  vp.save_dvar_vector_value();
//    save_identifier_string("s3");
  const char * str46;
  str46="s3";
  char* strx46=const_cast <char*> (str46);
  save_identifier_string(strx46);
  vp.save_dvar_vector_position();
//    save_identifier_string("s2");
  const char * str47;
  str47="s2";
  char* strx47=const_cast <char*> (str47);
  save_identifier_string(strx47);
  save_pointer_value((void*)(&Sinv));
//    save_identifier_string("s1");
  const char * str48;
  str48="s1";
  char* strx48=const_cast <char*> (str48);
  save_identifier_string(strx48);
  save_double_value(crho);
//    save_identifier_string("s10");
  const char * str49;
  str49="s10";
  char* strx49=const_cast <char*> (str49);
  save_identifier_string(strx49);
  save_double_value(cvar);
//    save_identifier_string("s11");
  const char * str50;
  str50="s11";
  char* strx50=const_cast <char*> (str50);
  save_identifier_string(strx50);
  
  gradient_structure::GRAD_STACK1->
      set_gradient_stack(df_lognormal_ss_multinomial_laplace_approximation_1);
  return vlog_det;
}
*/



/*
dvariable lognormal_multinomial_log_likelihood
  (const dvariable& vN,const dvector& q,const dvar_vector& vp,
  const dvar_matrix& vSinv,const dmatrix& Sinv,const dmatrix & S)
{

  MY_DOUBLE_TYPE N=value(vN);
  const dvector p=value(vp);

  dvector eta_hat=lognormal_ss_multinomial_newton_raphson(N,q,p,Sinv,S);
  
  MY_DOUBLE_TYPE ln_det=lognormal_ss_multinomial_laplace_approximation
    (N,q,p,Sinv,eta_hat);
  
  int mmin=p.indexmin();
  int mmax=p.indexmax();

  dvector tt(mmin,mmax);
  MY_DOUBLE_TYPE v=0.0;
  for (int i=mmin;i<=mmax;i++)
  {
    tt(i)=phi(p(i),eta_hat(i));
    v+=tt(i);
  }

  dvariable vlog_det=nograd_assign(ln_det);
//    save_identifier_string("t5");
  const char * str51;
  str51="t5";
  char* strx51=const_cast <char*> (str51);
  save_identifier_string(strx51);
  vSinv.save_dvar_matrix_position();
//    save_identifier_string("t4");
  const char * str52;
  str52="t4";
  char* strx52=const_cast <char*> (str52);
  save_identifier_string(strx52);
  eta_hat.save_dvector_value();
//    save_identifier_string("t3");
  const char * str53;
  str53="t3";
  char* strx53=const_cast <char*> (str53);
  save_identifier_string(strx53);
  eta_hat.save_dvector_position();
//    save_identifier_string("ps");
  const char * str54;
  str54="ps";
  char* strx54=const_cast <char*> (str54);
  save_identifier_string(strx54);
  vlog_det.save_prevariable_position();
//    save_identifier_string("rt");
  const char * str55;
  str55="rt";
  char* strx55=const_cast <char*> (str55);
  save_identifier_string(strx55);
  vN.save_prevariable_value();
//    save_identifier_string("s7");
  const char * str56;
  str56="s7";
  char* strx56=const_cast <char*> (str56);
  save_identifier_string(strx56);
  vN.save_prevariable_position();
//    save_identifier_string("s6");
  const char * str57;
  str57="s6";
  char* strx57=const_cast <char*> (str57);
  save_identifier_string(strx57);
  q.save_dvector_value();
//    save_identifier_string("s5");
  const char * str58;
  str58="s5";
  char* strx58=const_cast <char*> (str58);
  save_identifier_string(strx58);
  q.save_dvector_position();
//    save_identifier_string("s4");
  const char * str59;
  str59="s4";
  char* strx59=const_cast <char*> (str59);
  save_identifier_string(strx59);
  vp.save_dvar_vector_value();
//    save_identifier_string("s3");
  const char * str60;
  str60="s3";
  char* strx60=const_cast <char*> (str60);
  save_identifier_string(strx60);
  vp.save_dvar_vector_position();
//    save_identifier_string("s2");
  const char * str61;
  str61="s2";
  char* strx61=const_cast <char*> (str61);
  save_identifier_string(strx61);
  save_pointer_value((void*)(&Sinv));
//    save_identifier_string("s1");
  const char * str62;
  str62="s1";
  char* strx62=const_cast <char*> (str62);
  save_identifier_string(strx62);
  
  gradient_structure::GRAD_STACK1->
      set_gradient_stack(df_lognormal_ss_multinomial_laplace_approximation);
  return vlog_det;
}
*/
MY_DOUBLE_TYPE lognormal_ss_multinomial_laplace_approximation(const MY_DOUBLE_TYPE& N,
  const dvector& q,const dvector& p,const dmatrix & Sinv,
 const dvector& eta)
{
  int i,j,k;
  int mmin=q.indexmin();
  int mmax=q.indexmax();
  
  dvector t(mmin,mmax);
  dvector Nq=N*q;
  dvector grad(mmin,mmax);
  dvector w(mmin,mmax);
  symmetric_tridiagonal_dmatrix BHess(mmin,mmax);
  banded_symmetric_dmatrix BSinv(mmin,mmax,2);
  for (int i=mmin;i<mmax;i++)
  {
    BSinv(i,i)=Sinv(i,i);
    BSinv(i+1,i)=Sinv(i+1,i);
  }
  BSinv(mmax,mmax)=Sinv(mmax,mmax);

  dmatrix M(mmin,mmax,mmin,mmax);
  M.initialize();
  BHess.initialize();

  MY_DOUBLE_TYPE v=0.0;
  for (int i=mmin;i<=mmax;i++)
  {
    t(i)=phi(p(i),eta(i));
    v+=t(i);
  }
 /*
  MY_DOUBLE_TYPE mll=0.0; 
  for (int i=mmin;i<=mmax;i++)
  {
    mll-=Nq(i)*log(t(i));
  }
 */
  if (stupid_print_switch)
  {
    dvector pred=t/v;
    dvector obs=q;
    dvector r2=elem_div(square(obs-pred),elem_prod(pred,1.0-pred));
    stupid_print << setscientific() << setprecision(4) << setw(11) << t/v << endl;
    stupid_print << setscientific() << setprecision(4) << setw(11) << q << endl;
    stupid_print << setscientific() << setprecision(4) << setw(11) << p << endl;
    stupid_print << endl;
    stupid_print << setfixed() << setprecision(3) << setw(8) << eta << endl;
    if (wght_print_flag)
        stupid_print << "wght" << endl;
    else
        stupid_print << "length" << endl;
    stupid_print << setscientific() << setprecision(4) << setw(11) << r2 << endl;
    stupid_print << endl;
    stupid_print << endl;
  }
  MY_DOUBLE_TYPE mll=-Nq*log(t)+N*log(v)+0.5*eta*(Sinv*eta); 
  dvector d1(mmin,mmax);
  dvector d2(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    d1(i)=dphi_deta(p(i),eta(i));
    d2(i)=d2phi_deta2(p(i),eta(i));
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
    M(i,i)+=(-Nqt(i)+Nv)*d2(i);
    BHess(i,i)+=(-Nqt(i)+Nv)*d2(i);
    M(i,i)+=(Nqtt(i)-Nvv)*square(d1(i));
    // this next line is part of the rank one update
    BHess(i,i)+=Nqtt(i)*square(d1(i));
    for (int j=mmin;j<i;j++)
    {
      M(i,j)-=Nvv*d1(i)*d1(j);
      M(j,i)= M(i,j);
    }
    w(i)=sqrt(Nvv)*d1(i);
  }
  M+=Sinv;
  for (int i=mmin;i<mmax;i++)
  {
    BHess(i,i)+=BSinv(i,i);
    BHess(i+1,i)=BSinv(i+1,i);
  }
  BHess(mmax,mmax)+=BSinv(mmax,mmax);

  // choleski decompa to get log det
  if (M(1,1)<=0)
  {
    cerr << "Error matrix not positive definite in choleski_decomp"
      <<endl;
    ad_exit(1);
  }

  int rowsave=M.rowmin();
  int colsave=M.colmin();
  M.rowshift(1);
  M.colshift(1);
  int n=M.rowmax();

  dmatrix L(1,n,1,n);
  L.initialize();

  MY_DOUBLE_TYPE tmp;
  L(1,1)=sqrt(M(1,1));
  for (i=2;i<=n;i++)
  {
    L(i,1)=M(i,1)/L(1,1);
  }

  for (i=2;i<=n;i++)
  {
    for (j=2;j<=i-1;j++)
    {
      tmp=M(i,j);
      for (k=1;k<=j-1;k++)
      {
        tmp-=L(i,k)*L(j,k);
      }
      L(i,j)=tmp/L(j,j);
    }
    tmp=M(i,i);
    for (k=1;k<=i-1;k++)
    {
      tmp-=L(i,k)*L(i,k);
    }
    if (tmp<=0)
    {
      cerr << "Error matrix not positive definite in ln_det_choleski"
        <<endl;
      ad_exit(1);
    }
    L(i,i)=sqrt(tmp);
  }
  /*
  cout << norm2(Hess-L*trans(L)) << endl;
  dmatrix fvals=importance_sampling(Nq,p,eta,L,200,Sinv,N);
  ofstream ofs("import_sample")e
  ofs << trans(fvals) << endl;
  ad_exit(1);
  */

  MY_DOUBLE_TYPE log_lik=mll;
  // DDDDDD
  
  MY_DOUBLE_TYPE ld1=0.0;
  for (i=1;i<=n;i++)
  {
   // EEEEEEEE
    log_lik+=log(L(i,i));
    ld1+=log(L(i,i));
  }
  int ierr=0;
  //double ld=ln_det_choleski(BHess,ierr);
  // put the code in here for derivatives later
  MY_DOUBLE_TYPE ln_det=0.0;
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
          cerr << "Error matrix not positive definite in choleski_decomp"
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
            cerr << "Error matrix not positive definite in choleski_decomp"
              <<endl;
          ierr=1;
          return 0;
        }
        L(i,i)=sqrt(tmp);
      }
      for (int i=mmin;i<=mmax;i++)
      {
        if (L(i,i)>0.0)
          ln_det+=log(L(i,i));
      }
    }
  }


  dvector x=solve(BHess,w);
  ln_det+=0.5*log(1.0-(x*w));
 
  //cout << " ld1   ln_det  " << endl;
  //cout <<  ld1 << "  "  <<   ln_det << endl;
  //cout << " log_lik   ln_det+mll  " << endl;
  //cout <<  log_lik << "  "  <<   ln_det+mll << endl;
 
  log_lik=ln_det+mll;
  return log_lik;
}


static void df_lognormal_ss_multinomial_laplace_approximation(void)
{
  int i,j,k;
  //double dflog_det=restore_prevariable_derivative();
  verify_identifier_string("s1");
  dmatrix& Sinv=*(dmatrix*)(restore_pointer_value());
  verify_identifier_string("s2");
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
  dvar_matrix_position Sinvpos=restore_dvar_matrix_position();
  verify_identifier_string("t5");
  int mmin=q.indexmin();
  int mmax=q.indexmax();
  
  dvector dfeta(mmin,mmax);
  dfeta.initialize();
  dvector t(mmin,mmax);
  dvector dft(mmin,mmax);
  dft.initialize();
  dmatrix M(mmin,mmax,mmin,mmax);
  M.initialize();
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
  MY_DOUBLE_TYPE mll=-Nq*log(t)+N*v+0.5*eta*(Sinv*eta);
  // !!!!!!!!!GGGG
  MY_DOUBLE_TYPE dfmll=0.0;
  MY_DOUBLE_TYPE dfnv=0.0;
  MY_DOUBLE_TYPE dfnvv=0.0;
  dmatrix dfSinv(mmin,mmax,mmin,mmax);
  dfSinv.initialize();
  dvector d1(mmin,mmax);
  dvector d2(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    d1(i)=dphi_deta(p(i),eta(i));
    d2(i)=d2phi_deta2(p(i),eta(i));
  }
  MY_DOUBLE_TYPE Nv=N/v;
  MY_DOUBLE_TYPE Nvv=Nv/v;
  //dvector Nq(mmin,mmax);
  dvector Nqt(mmin,mmax);
  dvector Nqtt(mmin,mmax);
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
    M(i,i)+=(-Nqt(i)+Nv)*d2(i);
    M(i,i)+=(Nqtt(i)-Nvv)*square(d1(i));
    for (int j=mmin;j<i;j++)
    {
      M(i,j)-=Nvv*d1(i)*d1(j);
      M(j,i)= M(i,j);
    }
  }
  M+=Sinv;

  int rowsave=M.rowmin();
  int colsave=M.colmin();
  M.rowshift(1);
  M.colshift(1);
  int n=M.rowmax();

  if (M.colsize() != M.rowsize())
  {
    cerr << "Error in chol_decomp. Matrix not square" << endl;
    ad_exit(1);
  }
  M.rowshift(1);
  M.colshift(1);
  //int n=M.rowmax();

  dmatrix L(1,n,1,n);
  dmatrix dfL(1,n,1,n);
  dvector tmp(1,n);
  dmatrix tmp1(1,n,1,n);
  dmatrix dftmp1(1,n,1,n);
  dmatrix dfM(1,n,1,n);
  dvector dftmp(1,n);
  tmp.initialize();
  tmp1.initialize();
  dftmp.initialize();
  dftmp1.initialize();
  dfM.initialize();
  dfL.initialize();
#ifndef SAFE_INITIALIZE
    L.initialize();
#endif

  if (M(1,1)<=0)
  {
    cerr << "Error matrix not positive definite in choleski_decomp"
      <<endl;
    ad_exit(1);
  }
  L(1,1)=sqrt(M(1,1));
  for (i=2;i<=n;i++)
  {
    L(i,1)=M(i,1)/L(1,1);
  }

  for (i=2;i<=n;i++)
  {
    for (j=2;j<=i-1;j++)
    {
      tmp1(i,j)=M(i,j);
      for (k=1;k<=j-1;k++)
      {
        tmp1(i,j)-=L(i,k)*L(j,k);
      }
      L(i,j)=tmp1(i,j)/L(j,j);
    }
    tmp(i)=M(i,i);
    for (k=1;k<=i-1;k++)
    {
      tmp(i)-=L(i,k)*L(i,k);
    }
    if (tmp(i)<=0)
    {
      cerr << "Error matrix not positive definite in choleski_decomp"
        <<endl;
      ad_exit(1);
    }
    L(i,i)=sqrt(tmp(i));
  }
  MY_DOUBLE_TYPE log_det1=mll;
  for (i=1;i<=n;i++)
  {
    // EEEEEEE
    log_det1+=log(L(i,i));
  }
  MY_DOUBLE_TYPE log_det=log_det1;
 //*******************************************************************8
  //double log_det=log_det1;
  MY_DOUBLE_TYPE dflog_det1=dflog_det;
  for (i=1;i<=n;i++)
  {
    // EEEEEEE
    //log_det1+=log(L(i,i));
    dfL(i,i)+=dflog_det1/L(i,i);
  }
  //double log_det1=mll;
  dfmll+=dflog_det1;
  dflog_det1=0.0;

  for (i=n;i>=2;i--)
  {
    //L(i,i)=sqrt(tmp(i));
    dftmp(i)+=dfL(i,i)/(2.0*L(i,i));
    dfL(i,i)=0.0;
    for (k=i-1;k>=1;k--)
    {
      //tmp(i)-=L(i,k)*L(i,k);
      dfL(i,k)-=2.*dftmp(i)*L(i,k);
    }
    //tmp(i)=M(i,i);
    dfM(i,i)+=dftmp(i);
    dftmp(i)=0.0;
    for (j=i-1;j>=2;j--)
    {
      //L(i,j)=tmp1(i,j)/L(j,j);
      MY_DOUBLE_TYPE linv=1./L(j,j);
      dftmp1(i,j)+=dfL(i,j)*linv;
      dfL(j,j)-=dfL(i,j)*tmp1(i,j)*linv*linv;
      dfL(i,j)=0.0;
      for (k=j-1;k>=1;k--)
      {
        //tmp(i,j)-=L(i,k)*L(j,k);
        dfL(i,k)-=dftmp1(i,j)*L(j,k);
        dfL(j,k)-=dftmp1(i,j)*L(i,k);
      }
      //tmp(i,j)=M(i,j);
      dfM(i,j)+=dftmp1(i,j);
      dftmp1(i,j)=0.0;
    }
  }
  MY_DOUBLE_TYPE linv=1./L(1,1);
  MY_DOUBLE_TYPE linv2=linv*linv;
  for (i=n;i>=2;i--)
  {
    //L(i,1)=M(i,1)/L(1,1);
    dfM(i,1)+=dfL(i,1)*linv;
    dfL(1,1)-=dfL(i,1)*M(i,1)*linv2;
    dfL(i,1)=0.0;
  }
  //L(1,1)=sqrt(M(1,1));
  dfM(1,1)+=dfL(1,1)/(2.*L(1,1));
  //M+=Sinv;
  dfSinv+=dfM;
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

  for (int i=mmax;i>=mmin;i--)
  {
    for (int j=i-1;j>=mmin;j--)
    {
      //M(j,i)= M(i,j);
      dfM(i,j)+=dfM(j,i);
      dfM(j,i)=0.0;
      //M(i,j)-=Nvv*d1(i)*d1(j);
      dfNvv-=dfM(i,j)*d1(i)*d1(j);
      dfd1(i)-=dfM(i,j)*Nvv*d1(j);
      dfd1(j)-=dfM(i,j)*Nvv*d1(i);
    }
    //M(i,i)+=(Nqtt(i)-Nvv)*square(d1(i));
    dfNqtt(i)+=dfM(i,i)*square(d1(i));
    dfNvv-=dfM(i,i)*square(d1(i));
    dfd1(i)+=2.0*dfM(i,i)*(Nqtt(i)-Nvv)*d1(i);
    
    //M(i,i)+=(-Nqt(i)+Nv)*d2(i);
    dfNqt(i)-=dfM(i,i)*d2(i);
    dfNv+=dfM(i,i)*d2(i);
    dfd2(i)+=dfM(i,i)*(-Nqt(i)+Nv);
  }
  dfM.initialize();

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
    //d2(i)=d2phi_deta2(p(i),eta(i));
    dfp(i)+=dfd2(i)*d2phi_deta2_1(p(i),eta(i));
    dfeta(i)+=dfd2(i)*d2phi_deta2_2(p(i),eta(i));
    dfd2(i)=0.0;
    //d1(i)=dphi_deta(p(i),eta(i));
    dfp(i)+=dfd1(i)*dphi_deta_1(p(i),eta(i));
    dfeta(i)+=dfd1(i)*dphi_deta_2(p(i),eta(i));
    dfd1(i)=0.0;
  }
  //cout << "HHH " << dfeta(1) << " " << dfd1(i) << "  " << dfd2(i) << endl;


  //double mll=-Nq*log(t)+N*log(v)+0.5*eta*(Sinv*eta); 

  dfNq-=dfmll*log(t);
  dft-=dfmll*elem_div(Nq,t);
  dfN+=dfmll*log(v);
  dfv+=dfmll*N/v;
  dfSinv+=outer_prod(0.5*eta,eta);
  dfeta+=Sinv*eta;
  dfmll=0.0;
  

  for (i=mmax;i>=mmin;i--)
  {
    //v+=t(i);
    dft(i)+=dfv;
    //t(i)=phi(p(i),eta(i));
    //dfp(i)+=dft(i)*dphi_deta(p(i),eta(i))*psiprime(p(i),eta(i));
    dfp(i)+=dft(i)*dphi_dp(p(i),eta(i));
    dfeta(i)+=dft(i)*dphi_deta(p(i),eta(i));
    dft(i)=0.0;
  }
  //dvector Nq=N*q;
  dfN+=dfNq*q;

  // v=0.0;
  dfv=0.0;

  // DDDDD
   //cout << "PPP deriv of lap part wrt eta(1)" << dfeta(1) << endl;

  
  dvector cd=get_cross_derivatives(M,q,t,v,eta,N,dfeta,p);

  // add in the cross derivatives
 
  dfN+=cd(1);
  dfp+=cd(2,n+1).shift(1);
  for (int i=1;i<=n;i++)
  {
    for (int j=1;j<=n;j++)
    {
      dfSinv(i,j)+=cd(1+n+(i-1)*n+j);
    }
  }
 
 

  dfM.rowshift(rowsave);
  dfM.colshift(colsave);

  dfSinv.save_dmatrix_derivatives(Sinvpos);
  save_double_derivative(dfN,Npos);
  dfp.save_dvector_derivatives(vppos);
  //dfeta.save_dvector_derivatives(vetapos);
  //dfM.save_dmatrix_derivatives(MMpos);
}

  
dvector get_cross_derivatives(const dmatrix& M,const dvector& q,
  const dvector& t,const MY_DOUBLE_TYPE v,const dvector& eta,const MY_DOUBLE_TYPE& N,
  const dvector dfeta,const dvector& p)
{
  //lower_triangular dmatrix LT(1,n);
  int n=q.indexmax();

  dmatrix NN(1,n,1,1+n);
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
    NN(i,1)=(-q(i)/t(i)+1.0/v)*dphi_deta(p(i),eta(i));
  }

  dmatrix Hess(1,n,1,n);
  int mmin=1;
  int mmax=n;

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
  for (int i=1;i<=n;i++)
  {
    for (int j=1;j<=n;j++)
    {
      NN(i,j+1)=Hess(i,j)*psiprime(p(j),eta(j));
      cerr << "ERROR - SSMULT cross-derivatives invalid" << endl;
      ad_exit(1);
    }
  }
  dvector tmp(1,1+n+n*n);
  tmp.initialize();
 
  dmatrix Minv=inv(M);
  dmatrix uhatp=-Minv*NN;
  tmp(1,n+1)=dfeta*uhatp;
  for (int i=1;i<=n;i++)
  {
    for (int j=1;j<=n;j++)
    {
      dvector t=Minv(i)*eta(j)+Minv(j)*eta(i);
      tmp(1+n+(i-1)*n+j)=-0.5*(dfeta*t);
    }
  }
  
  return tmp;
}


  dvector lt1solve_solvet(const dmatrix& L,const dvector& y)
  {
    int mmin=y.indexmin();
    int mmax=y.indexmax();
    dvector v(mmin,mmax);
    dvector w(mmin,mmax);
    dvector asum(mmin,mmax);
    dvector bsum(mmin,mmax);
    v.initialize();
    v(mmin)=y(mmin)/L(mmin,mmin);
    for (int i=mmin+1;i<=mmax;i++)
    {
      asum(i)=L(i,i-1)*v(i-1);
      v(i)=(y(i)-asum(i))/L(i,i);
    }
    //cout << " BB " << v << endl;
    w(mmax)=v(mmax)/L(mmax,mmax);
    for (int i=mmax-1;i>=mmin;i--)
    {
      bsum(i)=L(i+1,i)*w(i+1);
      w(i)=(v(i)-bsum(i))/L(i,i);
    }
    //cout << " ZZ " << w << endl;
    return w;
  }
  
void dflt1solve_solvet(const dmatrix& L,const dvector& y,
  const dvector& _dfw,const dvector& _dfy,const dmatrix& _dfL)
{
  ADUNCONST(dvector,dfw)
  ADUNCONST(dvector,dfy)
  ADUNCONST(dmatrix,dfL)
  int mmin=y.indexmin();
  int mmax=y.indexmax();
  dvector v(mmin,mmax);
  dvector w(mmin,mmax);
  dvector asum(mmin,mmax);
  dvector bsum(mmin,mmax);
  v.initialize();
  v(mmin)=y(mmin)/L(mmin,mmin);
  for (int i=mmin+1;i<=mmax;i++)
  {
    asum(i)=L(i,i-1)*v(i-1);
    v(i)=(y(i)-asum(i))/L(i,i);
  }
  w(mmax)=v(mmax)/L(mmax,mmax);
  for (int i=mmax-1;i>=mmin;i--)
  {
    bsum(i)=L(i+1,i)*w(i+1);
    w(i)=(v(i)-bsum(i))/L(i,i);
  }
  //return w;
  // get dfw;
  //dvector dfw(mmin,mmax);
  //dvector dfy(mmin,mmax);
  dvector dfv(mmin,mmax);
  dvector dfasum(mmin,mmax);
  dvector dfbsum(mmin,mmax);
  //banded_lower_triangular_dmatrix dfL(mmin,mmax,2);
  dfv.initialize();
  dfy.initialize();
  dfasum.initialize();
  dfbsum.initialize();
  dfL.initialize();
  
  for (int i=mmin;i<=mmax-1;i++)
  {
    //w(i)=(v(i)-bsum(i))/L(i,i);
    dfv(i)+=dfw(i)/L(i,i);
    dfbsum(i)-=dfw(i)/L(i,i);
    dfL(i,i)-=dfw(i)*w(i)/L(i,i);
    dfw(i)=0.0;
    //bsum(i)=L(i+1,i)*w(i+1);
    dfL(i+1,i)+=dfbsum(i)*w(i+1);
    dfw(i+1)+=dfbsum(i)*L(i+1,i);
    dfbsum(i)=0.0;
  }
  //w(mmax)=v(mmax)/L(mmax,mmax);
  dfv(mmax)+=dfw(mmax)/L(mmax,mmax);
  dfL(mmax,mmax)-=dfw(mmax)*w(mmax)/L(mmax,mmax);
  for (int i=mmax;i>=mmin+1;i--)
  {
    //v(i)=(y(i)-asum(i))/L(i,i);
    dfy(i)+=dfv(i)/L(i,i);
    dfasum(i)-=dfv(i)/L(i,i);
    dfL(i,i)-=dfv(i)*v(i)/L(i,i);
    dfv(i)=0.0;
    //asum(i)=L(i,i-1)*v(i-1);
    dfL(i,i-1)+=dfasum(i)*v(i-1);
    dfv(i-1)+=dfasum(i)*L(i,i-1);
    dfasum(i)=0.0;
  }
  //v(mmin)=y(mmin)/L(mmin,mmin);
  dfy(mmin)+=dfv(mmin)/L(mmin,mmin);
  dfL(mmin)-=dfv(mmin)*v(mmin)/L(mmin,mmin);
  dfv(mmin)=0.0;

}


  dvector xlt1solve_solvet(const dmatrix& L,const dvector& w)
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
  
void xdflt1solve_solvet(const dmatrix& L,const dvector& w,
  const dvector& _dfw,const dvector& _dfx,const dmatrix& _dfL)
{
  ADUNCONST(dvector,dfx)
  ADUNCONST(dvector,dfw)
  ADUNCONST(dmatrix,dfL)
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
  dfw.initialize();
  dfasum.initialize();
  dfbsum.initialize();
  dfL.initialize();
  
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
#undef HOME_VERSIONN


