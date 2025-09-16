/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


#include "all.hpp"
#ifdef __MSC__
  MY_DOUBLE_TYPE mfexp(const MY_DOUBLE_TYPE& );
  dvariable mfexp(const prevariable& );
  dvar_vector mfexp(const dvar_vector& );
  dvector mfexp(const dvector& );
#else
  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
  dvariable mfexp(prevariable& );
  //dvar_vector mfexp(dvar_vector& );
  //dvector mfexp(dvector& );
#endif

vnrfss2::vnrfss2(dvar_matrix& _fm,int _nf,int _ng,dvar_matrix & _s,dvar_vector& _M,dvector& _Cobs,
 dvar_vector& _N,dvar_vector& _w,ivector& _wtflag,MY_DOUBLE_TYPE _beta,
 dvar_vector& _q, dvar_len_fish_stock_history * _pfsh,int _ir,int _ip,
  MY_DOUBLE_TYPE _rmax) : 
 nf(_nf), ng(_ng), s(_s),M(_M),Cobs(_Cobs), N(_N),
 q(_q),C(1,nf,1,ng),Fraw(1,nf,1,ng),F(_fm), Z(1,ng), Sraw(1,ng), S(1,ng),Chat(1,nf),
 wtflag(_wtflag),beta(_beta),phi(1.0),logq(1,nf),simflag(1),
 pfsh(_pfsh),ir(_ir),ip(_ip), sN(1,_nf), swN(1,_nf),rmax(_rmax), Zraw(1,ng) 
{
  w=_w/1000.;
  Nbefore=elem_prod(N,exp(-0.5*M));
  for (int fi=1;fi<=nf;fi++)
  {
    sN(fi)=elem_prod(s(fi),Nbefore);
    if (wtflag(fi))
    {
      swN(fi)=elem_prod(w,sN(fi));
    }
  }
}

vnrfss2::vnrfss2(int _nf,int _ng,dvar_matrix _s,dvar_vector& _M,dvector _Cobs,
 dvar_vector _N,dvar_vector _w,ivector _wtflag,MY_DOUBLE_TYPE _beta,
 dvar_vector _q, dvar_len_fish_stock_history * _pfsh,int _ir,int _ip,
 ivector _fishin) : 
 nf(_nf), ng(_ng), s(_s),M(_M),Cobs(_Cobs), N(_N),
 q(_q),C(1,nf,1,ng),F(1,nf,1,ng),Fraw(1,nf,1,ng),Zraw(1,ng), Z(1,ng), Sraw(1,ng), S(1,ng),
 Chat(1,nf),
 wtflag(_wtflag),beta(_beta),phi(1.0),logq(1,nf),simflag(1),
 pfsh(_pfsh),ir(_ir),ip(_ip),sN(1,_nf),swN(1,_nf),fishin(_fishin),rmax(0.7) 
{
  cerr << "This needs fixing" << endl;
  ad_exit(1);
  w=_w/1000.;
  Nbefore=elem_prod(N,M);
  for (int fi=1;fi<=nf;fi++)
  {
    sN(fi)=elem_prod(s(fi),Nbefore);
    if (wtflag(fi))
    {
      swN(fi)=elem_prod(w,sN(fi));
    }
  }
}
  
 dvariable vnrfss2::testnr(void)
 {
   dvariable fpen=0.0;
   int fi;
   int log_flag=0;
    
   if (!log_flag)
   {
     get_initial_q();
     calculate_Fraw(); 
     calculate_Zraw(); 
     fpen=calculate_S(); 
     calculate_F(); 
     int mmin=Z.indexmin();
     int mmax=Z.indexmax();
   }
   dvar_vector tmp=elem_prod(Nbefore,S);
   Nafter=elem_prod(tmp,exp(-0.5*M));
   return fpen;
 }
  
  
 void vnrfss2::get_initial_q(void)
 {
   int fi;
   int j;
    // get initial q values
    for (fi=1;fi<=nf;fi++)
    {
      dvariable ssum=0.0;
      if (!wtflag(fi))
      {
        ssum=sum(sN(fi));
        q(fi)=Cobs(fi)/ssum;
      }
      else
      {
        ssum=sum(swN(fi));
        q(fi)=Cobs(fi)/ssum;
      }
    }
    logq=log(1.e-10+q);
  }
      

  void vnrfss2::calculate_F(void)
  {
    int fi;
    int j;
    dvar_vector mult=elem_div(Z,Zraw);
    for (fi=1;fi<=nf;fi++)
    {
      F(fi)=log(1.e-20+elem_prod(Fraw(fi),mult));
    }
  }


  void vnrfss2::calculate_Fraw(void)
  {
    int fi;
    int j;
    for (fi=1;fi<=nf;fi++)
    {
      Fraw(fi)=q(fi)*s(fi);
    }
  }

  void vnrfss2::calculate_Zraw(void)
  {
    int fi;
    int j;
    Zraw.initialize();
    for (fi=1;fi<=nf;fi++)
    {
      for (j=1;j<=ng;j++)
      {
        Zraw(j)+=Fraw(fi,j);
      }
    }
  }

  dvariable vnrfss2::calculate_S()
  {
    int mmin=S.indexmin();
    int mmax=S.indexmax();
    dvariable tpen=0.0;
    Sraw=1.0-Zraw;
    for (int i=1;i<=mmax;i++)
    { 
      dvariable ppen=0.0;
      int rmin=static_cast<int>(1.0-rmax);
      S(i)=posfun(0.9-Zraw(i),0.1,ppen)+0.1;
      Z(i)=1.0-S(i);
      tpen+=ppen;
    }
    return tpen;
  }

