/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
int dumflag=1;

# define USE_DD

#include "all.hpp"
extern int compare_switch;

double my_double(const dd_real& x)
{
#if (defined(USE_DD) && (__MSVC32__>=8 || defined(linux) || defined(__ADMING__) || defined(_MAC)))
  return to_double(x);
#else
  return double(x);
#endif
}
int nrcounter=0;
ddnrf::ddnrf(int _nf,int _ng,const ddmatrix& _s,const ddvector& _M,
 const ddvector& _N,dvector& _Cobs,const ddvector& _w,
 ivector& _wtflag,dd_real _beta, dvar_fish_stock_history * _pfsh,
 int _ir,int _ip, dd_real _rmax) : 
 F(1,_nf,1,ng),nf(_nf), ng(_ng), s(_s),M(_M),Cobs(_Cobs), N(_N),
 q(1,_nf),C(1,nf,1,ng),Z(1,_ng), S(1,_ng),Chat(1,nf),w(_w),
 wtflag(_wtflag),beta(_beta),phi(1.0),logq(1,nf),simflag(1),
 pfsh(_pfsh),ir(_ir),ip(_ip), sN(1,_nf), swN(1,_nf),rmax(_rmax),
 have_weights(0)
{
  w/=1000.;
  
  dfF.allocate(1,nf,1,ng,1,nf);
  dfsN.allocate(1,nf,1,ng,1,nf);
  dfswN.allocate(1,nf,1,ng,1,nf);
  dfs.allocate(1,nf,1,ng,1,nf);
  dfC.allocate(1,nf,1,ng,1,nf);
  dfS.allocate(1,ng,1,nf);
  dfq.allocate(1,nf,1,nf);
  dfZ.allocate(1,ng,1,nf);
  dfChat.allocate(1,nf,1,nf);
  dfw.allocate(1,ng,1,nf);
  dfM.allocate(1,ng,1,nf);
  dfN.allocate(1,ng,1,nf);

  dFdq.allocate(1,nf,1,ng,1,nf);
  dsdq.allocate(1,nf,1,ng,1,nf);
  //dCdq.allocate(1,nf,1,ng,1,nf);
  dSdq.allocate(1,ng,1,nf);
  dZdq.allocate(1,ng,1,nf);
  dwdq.allocate(1,ng,1,nf);
  dMdq.allocate(1,ng,1,nf);
  dNdq.allocate(1,ng,1,nf);

  dfF.initialize();
  dfsN.initialize();
  dfswN.initialize();
  dfs.initialize();
  dfC.initialize();
  dfS.initialize();
  dfq.initialize();
  dfZ.initialize();
  dfChat.initialize();
  dfw.initialize();
  dfM.initialize();
  dfN.initialize();
  
  dFdq.initialize();
  dsdq.initialize();
  //dCdq.initialize();
  dSdq.initialize();
  dZdq.initialize();
  dwdq.initialize();
  dMdq.initialize();
  dNdq.initialize();


  int fi;
  for (fi=1;fi<=nf;fi++)
  {
    have_weights=1;
    switch (pfsh->age_flags(92))
    {
    case 0:   // added for tags  DF May 08 08
    case 2:
      sN(fi)=elem_prod(s(fi),N);
      break;
    case 3:
    case 4:
      sN(fi)=elem_prod(s(fi),elem_prod(N,exp(-0.5*M))); // ss2 option
      break;
    default:
      cerr << "illegal age_flags(92) value" << endl;
      ad_exit(1);
    }
    if (wtflag(fi))
    {
      swN(fi)=elem_prod(w,sN(fi));
    }
  }
}

ddnrf::ddnrf(int _nf,int _ng,ddmatrix _s,ddvector& _M,dvector _Cobs,
 ddvector _N,ddvector _w,ivector _wtflag,dd_real _beta,
 ddvector _q, dvar_fish_stock_history * _pfsh,int _ir,int _ip,
 ivector _fishin, double _rmax) : 
 nf(_nf), ng(_ng), s(_s),M(_M),Cobs(_Cobs), N(_N),
 q(_q),C(1,nf,1,ng),F(1,nf,1,ng), Z(1,ng),S(1,ng),
 Chat(1,nf),w(_w),
 wtflag(_wtflag),beta(_beta),phi(1.0),logq(1,nf),simflag(1),
 pfsh(_pfsh),ir(_ir),ip(_ip),sN(1,_nf),swN(1,_nf),fishin(_fishin),rmax(_rmax) 
{
  w/=1000.;
}

 dd_real ddnrf::testnr(void)
 {
   dd_real fpen=0.0;
   ddvector d1(1,nf);
   ddvector d2(1,nf);
   ddvector q1(1,nf);
   ddvector q2(1,nf);
   ddvector qders(1,nf);
   ddvector dders(1,nf);
   dd_real delta=1.e-6;
   int fi;
   int log_flag=0;
   
   if (!log_flag)
   {
     get_initial_q_ss2();
     icount=1;
     if(pfsh->age_flags(92) !=3)  // check for ss2 option
     {
       newton_raphson_loop();
     } 
     else
     {
       calculate_sN();
       calculate_F(); 
       calculate_Z(); 
       calculate_S(); 
     }
   }
   else
   {
      get_initial_q();
      return fpen;
      get_initial_q();
      icount=1;
      do
      {
        // newton raphson loop
        ddvector d=calculate_ld();
        dd_real nd=norm(d);
        if (nd<1.e-8) 
        {
          cout << ir << " " << ip << " " << icount  << endl;
          break;
        }
        if (icount > 5)
        {
          cout << ir << " " << ip << " " << icount  << endl;
          break;
        }
        if (icount > 15 && nd < 1.e-6) 
        {
          cout << ir << " " << ip << " " << icount  << endl;
          break;
        }
        if (icount > 25 && nd < 1.e-4) 
        {
          cout << ir << " " << ip << " " << icount  << endl;
          break;
        }
        if (icount > 40 ) 
        {
          cout << ir << " " << ip << " " << icount  << endl;
          break;
        }
        ddmatrix J=calculate_lJ();

        ddvector h=solve(trans(J),d);
        logq-=h;
        
        icount++;
      }
      while(1);
    }


    if (ir==1 && ip==1)
    {
      //ofstream ofs("catch_pred_obs",ios::app);
      //ofs << " Chat: " << Chat << " Cobs: " << Cobs << endl;
      /*
      ofs << "q(6)" << endl << q(6) << " " << endl
          << "s(6)" << endl << s(6) << endl
          << endl;
      ofs << "make_dmatrix(F)" << endl;
      ofs << make_dmatrix(F) << endl << endl;
      ofs << "Z" << endl;
      ofs << Z << endl;
      ofs << "Chat" << endl;
      ofs << Chat << endl;
      ofs << "Cobs" << endl;
      ofs << Cobs << endl;
      ofs << "Chat-Cobs" << endl;
      ofs << Chat-Cobs << endl;
      ofs << "N" << endl;
      ofs << N << endl;
      */
    }

    calculate_derivatives();

    implict_derivative_calculations();

    return fpen;
  }
  
 void ddnrf::implict_derivative_calculations(void)
 {
   dFdq.initialize();
   dsdq.initialize();
    //dCdq.initialize();
   dSdq.initialize();
   dZdq.initialize();
   dwdq.initialize();
   dMdq.initialize();
   dNdq.initialize();
   int sgn=0;
#if defined(USE_DD)
   double lnd=ln_det(make_dmatrix(dfq),sgn);
   if (lnd<-10.0)
   {
     //cout << " sum(enf) + " << sum(N) << endl;
     //cerr << "possible underlfow in" 
       // " implict_derivative_calculations" << lnd << endl
       //  << setw(11) << make_dmatrix(dfq) << endl;
   }
#else
   double lnd=ln_det(dfq,sgn);
   if (lnd<-10.0)
   {
     //cout << " sum(enf) + " << sum(N) << endl;
     //cerr << "possible underlfow in" 
     //   " implict_derivative_calculations" << lnd << endl
     //    << setw(11) << dfq << endl;
   }
#endif
   ddmatrix dfqinv=inv(trans(dfq));
   int fi,j;
   for (fi=1;fi<=nf;fi++)
   {
     for (j=1;j<=ng;j++)
     {
       dFdq(fi,j)=-(dfqinv*dfF(fi,j));
       dsdq(fi,j)=-(dfqinv*dfs(fi,j));
       //dCdq(fi,j)=dfqinv*dfC(fi,j);
     }
   }
 
   for (j=1;j<=ng;j++)
   {
     dSdq(j)=-(dfqinv*dfS(j));
     dZdq(j)=-(dfqinv*dfZ(j));
     dwdq(j)=-(dfqinv*dfw(j));
     dMdq(j)=-(dfqinv*dfM(j));
     dNdq(j)=-(dfqinv*dfN(j));
   }
 }

 void ddnrf::get_initial_q(void)
 {
   int fi;
   int j;
    // get initial q values
    for (fi=1;fi<=nf;fi++)
    {
      dd_real ssum=0.0;
      if (!wtflag(fi))
      {
        //sum=s(fi)*N;
        ssum=sum(sN(fi));
        q(fi)=Cobs(fi)/ssum;
        if (q(fi)<0.0)
        {
          cout << " q("<<fi <<") < 0 " << endl;
        }
      }
      else
      {
        //sum=s(fi)*elem_prod(w,N);
        ssum=sum(swN(fi));
        q(fi)=Cobs(fi)/ssum;
        if (q(fi)<0.0)
        {
          cout << " q("<<fi <<") < 0 " << endl;
        }
      }
    }
    logq=log(1.e-15+q);
  }
 void ddnrf::get_initial_q_ss2(void)
 {
   int fi;
   int j;
    // get initial q values
    //ddvector N1=elem_prod(N,exp(-0.5*M));
    ddvector N1=N;
    ddvector sN1(1,nf);
    ddvector wN1=elem_prod(w,N1);
    for (fi=1;fi<=nf;fi++)
    {
      dd_real ssum=0.0;
      if (!wtflag(fi))
      {
        //sum=s(fi)*N;
        ssum=s(fi)*N1;
        //ssum=sum(sN(fi));
        q(fi)=Cobs(fi)/(1.e-10+ssum);
        if (q(fi)<0.0)
        {
          cout << " q("<<fi <<") < 0 " << endl;
        }
      }
      else
      {
        //sum=s(fi)*elem_prod(w,N);
        ssum=s(fi)*wN1;
        q(fi)=Cobs(fi)/(1.e-8+ssum);
        if (q(fi)<0.0)
        {
          cout << " q("<<fi <<") < 0 " << endl;
        }
      }
    }
    logq=log(1.e-15+q);
  }

  ddvector ddnrf::calculate_ld()
  {
    q=exp(logq);
    return calculate_d();
  }
  
  ddvector ddnrf::calculate_d()
  {
   //calculate_sN();
   calculate_F(); 
   calculate_Z(); 
   calculate_S(); 
   //calculate_T(); 
   //calculate_That(); 
   calculate_Chat(); 
   //ddvector d=Chat-phi*Cobs;
   ddvector d=Chat-Cobs;
   return d;
  }
  void ddnrf::calculate_derivatives()
  {
    dfF.initialize();
    dfsN.initialize();
    dfswN.initialize();
    dfs.initialize();
    dfC.initialize();
    dfS.initialize();
    dfq.initialize();
    dfZ.initialize();
    dfChat.initialize();
    dfw.initialize();
    dfM.initialize();
    dfN.initialize();

   df_calculate_Chat(); 
   df_calculate_S(); 
   df_calculate_Z(); 
   df_calculate_F(); 
   //df_calculate_sN();
  }
 // 
  ddmatrix ddnrf::calculate_J()
  {
   //1/280 	−4/105 	1/5 	−4/5 	0 	4/5 	−1/5 	4/105 	−1/280 	
  
    int fi;
    dd_real diff1=1.e-12;
    dd_real diff2=3.e-12;
    ddmatrix J(1,nf,1,nf);
    //cout << " enter diff " << endl;
    //cin >> diff; 
    for (fi=1;fi<=nf;fi++)
    {
      dd_real sd1=(diff1-q(fi))+q(fi);
      dd_real tmp1=q(fi);
      q(fi)+=sd1;
      ddvector fp1=calculate_d();
      q(fi)=tmp1;
      q(fi)-=sd1;
      ddvector fm1=calculate_d();
      q(fi)=tmp1;
      ddvector t1=(fp1-fm1)/(2.0*sd1);
      //J(fi)=(fp1-fm1)/(2.0*sd1);

     
      dd_real sd2=(diff2-q(fi))+q(fi);
      dd_real tmp2=q(fi);
      q(fi)+=sd2;
      ddvector fp2=calculate_d();
      q(fi)=tmp2;
      q(fi)-=sd2;
      ddvector fm2=calculate_d();
      q(fi)=tmp2;
      ddvector t2=(fp2-fm2)/(2.0*sd2);
      ddvector tmp=calculate_d();

      J(fi)=(9.0*t1-t2)/8.0;
      
    }
    return trans(J);
  }

  dd3_array ddnrf::calculate_alpha()
  {
    int fi,j;
    dd_real diff=1.e-6;
    dd3_array u_alpha(1,nf,1,ng,1,nf);
    ddmatrix J(1,nf,1,nf);
    //cout << " enter diff " << endl;
    //cin >> diff; 
    for (fi=1;fi<=nf;fi++)
    {
      for (j=1;j<=ng;j++)
      {
        dd_real sd=(diff-s(fi,j))+s(fi,j);
        dd_real tmp=s(fi,j);
        s(fi,j)+=sd;
        ddvector fp=calculate_d();
        //cout << " fp = " << fp << endl;
        s(fi,j)=tmp;
        s(fi,j)-=sd;
        ddvector fm=calculate_d();
        s(fi,j)=tmp;
        u_alpha(fi,j)=(fp-fm)/(2.0*sd);
      }
    }
    return u_alpha;
  }

  ddmatrix ddnrf::calculate_lJ()
  {
  
    int fi;
    ddmatrix J(1,nf,1,nf);
    dd_real diff=1.e-5;
    for (fi=1;fi<=nf;fi++)
    {
      dd_real sd=(diff-logq(fi))+logq(fi);
      dd_real tmp=logq(fi);
      logq(fi)+=sd;
      ddvector fp=calculate_ld();
      logq(fi)=tmp;
      logq(fi)-=sd;
      ddvector fm=calculate_ld();
      logq(fi)=tmp;
      J(fi)=(fp-fm)/(2.0*sd);
    }
    return J;
  }
 // 
 // 
 //     
  void ddnrf::calculate_F(void)
  {
    int fi;
    int j;
    for (fi=1;fi<=nf;fi++)
    {
      F(fi)=q(fi)*s(fi);
    }
  }

  void ddnrf::calculate_Z(void)
  {
    int fi;
    int j;
    Z.initialize();
    for (fi=1;fi<=nf;fi++)
    {
      for (j=1;j<=ng;j++)
      {
        Z(j)+=F(fi,j);
      }
    }
    if (pfsh->age_flags(92)!=3 && pfsh->age_flags(92)!=4)
    {
      Z+=M;
    }
  }

  void ddnrf::calculate_sN()
  {
    for (int fi=1;fi<=nf;fi++)
    {
      sN(fi)=elem_prod(s(fi),N);
      if (wtflag(fi))
      {
        swN(fi)=elem_prod(w,sN(fi));
      }
    }
  }

  void ddnrf::df_calculate_sN()
  {
    int j;
    for (int fi=1;fi<=nf;fi++)
    {
      if (wtflag(fi))
      {
        //swN(fi)=elem_prod(w,sN(fi));
        for (j=1;j<=ng;j++)
        {
          dfw(j)+=dfswN(fi,j)*sN(fi,j);
          dfsN(fi,j)+=dfswN(fi,j)*w(j);
        }
        dfswN(fi)=0.0;
      }
      //sN(fi)=elem_prod(s(fi),N);
      for (j=1;j<=ng;j++)
      {
        dfs(fi,j)+=dfsN(fi,j)*N(j);
        dfN(j)+=dfsN(fi,j)*s(fi,j);
      }
      dfsN(fi)=0.0;
    }
  }

  void ddnrf::calculate_S()
  {
    int mmin=S.indexmin();
    int mmax=S.indexmax();
    for (int i=1;i<=mmax;i++)
    { 
      if (Z(i)<=rmax)
      {
          S(i)=exp(-Z(i));
      }
      else
      {
        dd_real dd=Z(i)-rmax;
        S(i)=exp(-rmax)-exp(-rmax)*dd;
        dd_real ppen=0.0;
        //Z(i)=rmax-posfun(0.2-dd,0.1,ppen)+0.2;
      }
    }
  }

  void ddnrf::calculate_Chat(void)
  {
    int fi,j;
    Chat.initialize();
    for (j=1;j<=ng;j++)
    {
      dd_real tmp=1.0/Z(j)*(1.0-S(j))*N(j);
      if (pfsh->age_flags(92)==4)  // Popes aprox
      {
        tmp*=exp(-0.5*M(j));
      }
      for (fi=1;fi<=nf;fi++)
      {
        C(fi,j)=F(fi,j)*tmp;
      }
    }
    for (fi=1;fi<=nf;fi++)
    {
      if (!wtflag(fi))
      {
        Chat(fi)=sum(C(fi));
      }
      else
      {
        Chat(fi)=w*C(fi);
      }
    }
  }

        
  ddmatrix ddnrf::calculate_J_vanal(void)
  {
    int i,j,k;
    ddmatrix Jac(1,nf,1,nf);
    Jac.initialize();
    for (i=1;i<=nf;i++)
    {
      for (k=1;k<=nf;k++)
      {
        for (j=1;j<=ng;j++)
        {
          if (i==k)  
          {
            if (Z(j)<=rmax)
            {
              dd_real tmp=((1.0-S(j))*(1.0-F(i,j)/Z(j))
                +F(i,j)*S(j))*s(k,j)/Z(j)*N(j);
              if (isnan(tmp))
              {
                cout << "got a nan" << endl;
              }
              if (pfsh->age_flags(92)==4)
              {
                tmp*=exp(-0.5*M(j));
              }
              if (wtflag(i)==0)
              {
                Jac(i,k)+=tmp;
              }
              else
              {
                Jac(i,k)+=w(j)*tmp;
              }
            }
            else
            {
              dd_real tmp=((1.0-S(j))*(1.0-F(i,j)/Z(j))
                +F(i,j)*exp(-rmax))*s(k,j)/Z(j)*N(j);
              if (isnan(tmp))
              {
                cout << "got a nan" << endl;
              }
              if (pfsh->age_flags(92)==4)
              {
                tmp*=exp(-0.5*M(j));
              }
              if (wtflag(i)==0)
              {
                Jac(i,k)+=tmp;
              }
              else
              {
                Jac(i,k)+=w(j)*tmp;
              }
            }
          }
          else
          {
            if (Z(j)<=rmax)
            {
              dd_real tmp=( -(1.0-S(j))/Z(j)+S(j))*s(k,j)*F(i,j)/Z(j)*N(j);
              if (isnan(tmp))
              {
                cout << "got a nan" << endl;
              }
              if (pfsh->age_flags(92)==4)
              {
                tmp*=exp(-0.5*M(j));
              }
              if (wtflag(i)==0)
              {
                Jac(i,k)+=tmp;
              }
              else
              {
                Jac(i,k)+=w(j)*tmp;
              }
            }
            else
            {
              dd_real tmp=
                ( -(1.0-S(j))/Z(j)+exp(-rmax))*s(k,j)*F(i,j)/Z(j)*N(j);
              if (isnan(tmp))
              {
                cout << "got a nan" << endl;
              }
              if (pfsh->age_flags(92)==4)
              {
                tmp*=exp(-0.5*M(j));
              }
              if (wtflag(i)==0)
              {
                Jac(i,k)+=tmp;
              }
              else
              {
                Jac(i,k)+=w(j)*tmp;
              }
            }
          }
        }
      }
      if (fabs(Jac(i,i))<1.e-6)
      {
        Jac(i,i)=1.e-6;
      }
    }
    return Jac;
  }
  ddmatrix ddnrf::calculate_J_anal(void)
  {
    int i,j,k;
    ddmatrix Jac(1,nf,1,nf);
    Jac.initialize();
    for (i=1;i<=nf;i++)
    {
      for (k=1;k<=nf;k++)
      {
        for (j=1;j<=ng;j++)
        {
          if (i==k)  
          {
            if (Z(j)<=rmax)
            {
              Jac(i,k)+=
                ((1.0-S(j))*(1.0-F(i,j)/Z(j))
                +F(i,j)*S(j))*s(k,j)/Z(j)*N(j);
            }
            else
            {
              Jac(i,k)+=
                ((1.0-S(j))*(1.0-F(i,j)/Z(j))
                +F(i,j)*exp(-rmax))*s(k,j)/Z(j)*N(j);
            }
          }
          else
          {
            if (Z(j)<=rmax)
            {
              Jac(i,k)+=
                ( -(1.0-S(j))/Z(j)+S(j))*s(k,j)*F(i,j)/Z(j)*N(j);
            }
            else
            {
              Jac(i,k)+=
                ( -(1.0-S(j))/Z(j)+exp(-rmax))*s(k,j)*F(i,j)/Z(j)*N(j);
            }
          }
        }
      }
    }
    return Jac;
  }
        
  dd3_array ddnrf::calculate_alpha_ders(void)
  {
    int i,j,k;
    dd3_array u_alpha(1,nf,1,ng,1,nf);
    for (i=1;i<=nf;i++)
    {
      for (k=1;k<=nf;k++)
      {
        for (j=1;j<=ng;j++)
        {
          if (i==k)  
          {
            u_alpha(k,j,i)=
              ((1.0-S(j))*(1.0-F(i,j)/Z(j))+F(i,j)*S(j))*q(k)/Z(j)*N(j);
          }
          else
          {
            u_alpha(k,j,i)=
              ( -(1.0-S(j))/Z(j)+S(j))*q(k)*F(i,j)/Z(j)*N(j);
              //(-S(j)*(1.0-F(i,j)/Z(j))+F(i,j)*S(j))*q(k)/Z(j)*N(j);
          }
        }
      }
    }
    return u_alpha;
  }
        
  ddvector ddnrf::calculate_N_ders(void)
  {
    int i,j,k;
    ddvector N_alpha(1,ng);
    N_alpha.initialize();
    for (i=1;i<=nf;i++)
    {
      for (j=1;j<=ng;j++)
      {
        N_alpha(j)+=F(i,j)/Z(j)*(1.0-S(j));
      }
    }
    return N_alpha;
  }

  void ddnrf::df_calculate_Chat(void)
  {
  /*
    int fi,j;
    Chat.initialize();
    for (j=1;j<=ng;j++)
    {
      dd_real tmp=1.0/Z(j)*(1.0-S(j))*N(j);
      if (pfsh->age_flags(92)==4)  // Popes aprox
      {
        tmp*=exp(-0.5*M(j));
      }
      for (fi=1;fi<=nf;fi++)
      {
        C(fi,j)=F(fi,j)*tmp;
      }
    }
    for (fi=1;fi<=nf;fi++)
    {
      if (!wtflag(fi))
      {
        Chat(fi)=sum(C(fi));
      }
      else
      {
        Chat(fi)=w*C(fi);
      }
    }
   */
 
    int fi,j;
    dfChat=0.0; 
    for (fi=1;fi<=nf;fi++)
    {
      dfChat(fi,fi)=1.0;
    }
    for (fi=1;fi<=nf;fi++)
    {
      if (!wtflag(fi))
      {
        //Chat(fi)=sum(C(fi));
        for (j=1;j<=ng;j++)
        {
          dfC(fi,j)+=dfChat(fi);
        }
        dfChat(fi)=0.0; 
      }
      else
      {
        //Chat(fi)=w*C(fi);
        for (j=1;j<=ng;j++)
        {
          dfw(j)+=C(fi,j)*dfChat(fi);
          dfC(fi,j)+=w(j)*dfChat(fi);
        }
        dfChat(fi)=0.0;
      }
    }
    for (j=1;j<=ng;j++)
    {
      dd_real tmp=1.0/Z(j)*(1.0-S(j))*N(j);
      ddvector dftmp(1,nf);
      dftmp.initialize();
      for (fi=1;fi<=nf;fi++)
      {
        //C(fi,j)=F(fi,j)*tmp;
        dfF(fi,j)+=dfC(fi,j)*tmp;
        dftmp+=dfC(fi,j)*F(fi,j);
        dfC(fi,j)=0.0;
      }
      if (pfsh->age_flags(92)==4)  // Popes aprox
      {
         cerr << "Not done yet" << endl;
         ad_exit(1);
        //tmp*=exp(-0.5*M(j));
        //dfM(j)
      }
      //dd_real tmp=1.0/Z(j)*(1.0-S(j))*N(j);
      dfN(j)+=dftmp/Z(j)*(1.0-S(j));
      dfS(j)-=dftmp/Z(j)*N(j);
      dfZ(j)-=dftmp*tmp/Z(j);
      dftmp=0.0;
    }
  }

  void ddnrf::df_calculate_S()
  {
    int mmin=S.indexmin();
    int mmax=S.indexmax();
    for (int i=1;i<=mmax;i++)
    { 
      if (Z(i)<=rmax)
      {
        //S(i)=exp(-Z(i));
        dfZ(i)-=dfS(i)*S(i);
        dfS(i)=0.0;
      }
      else
      {
        //S(i)=exp(-rmax)-exp(-rmax)*dd;
        ddvector dfdd=-exp(-rmax)*dfS(i);
        dfS(i)=0.0;
        //dd_real dd=Z(i)-rmax;
        dfZ(i)+=dfdd;
        dfdd=0.0;
      }
    }
  }

  void ddnrf::df_calculate_Z(void)
  {
    int fi;
    int j;
    //Z.initialize();
    if (pfsh->age_flags(92)!=3 && pfsh->age_flags(92)!=4)
    {
      //Z+=M;
      dfM+=dfZ;
    }
    for (fi=1;fi<=nf;fi++)
    {
      for (j=1;j<=ng;j++)
      {
        //Z(j)+=F(fi,j);
        dfF(fi,j)+=dfZ(j);
      }
    }
    //Z.initialize();
    dfZ=0.0;
  }

  void ddnrf::df_calculate_F(void)
  {
    int fi;
    int j;
    for (fi=1;fi<=nf;fi++)
    {
      for (j=1;j<=ng;j++)
      {
        //F(fi,j)=q(fi)*s(fi,j);
        dfq(fi)+=dfF(fi,j)*s(fi,j);
        dfs(fi,j)+=dfF(fi,j)*q(fi);
        dfF(fi,j)=0.0;
      }
      if (fabs(dfq(fi,fi))<1.e-8) dfq(fi,fi)=1.e-8;
    }
  }

  void ddnrf::newton_raphson_loop(void)
  {
    int conv_flag=0;
    dd_real nd;
    // here we use the ss2 parameterization to get vector of initial q's
    // but this is already done in containing function so comment out
    // get_inital_q_ss2();
    
    do
    {
      nrcounter++;
      // newton raphson loop
      ddvector d=calculate_d();
       
      if (norm(q)> 1.e+4)
      {
        cerr << "real_ddnrc3.cpp Not enough fish for projection in region "
             << endl << ir <<  "  period  " << ip << endl;
        //break;
      }
      nd=norm(d);

      if (nd>1.e+05)
      {
        cout << ir << " " << ip << " " << icount  << " " 
             << setscientific() << setprecision(4) << my_double(nd) << endl;
      }
      
#if !defined(MY_REAL_DOUBLE) 
      if (nd<1.e-22) 
#else
      if (nd<1.e-13) 
#endif
      {
        //cout << ir << " " << ip << " " << icount  << endl;
        conv_flag=1;
        break;
      }
      if (icount == 8 ) 
      {
        //break;
      }
 
      if (icount > 40 ) 
      {
        cout << ir << " " << ip << " " << icount  << " " 
             << setscientific() << setprecision(4) << my_double(nd) << endl;
        break;
      }
 
      
      //if (pfsh->age_flags(92)==4)
      //  J=calculate_J();
      //else
      ddmatrix J=calculate_J_vanal();
      //cout << eigenvalues(make_dmatrix(J)) << endl;
      //cout << endl << eigenvectors(make_dmatrix(J)) << endl;
#if defined(MY_REAL_DOUBLE) 
#endif
      
      if (compare_switch)
      {
        ddmatrix J1=calculate_J_vanal();
      }
      ddvector h=solve(J,d);
      dd_real cd=norm(J*h-d);
//      cout << " d = " << d << " h = " << h << " cd = " << cd << endl;
      if (cd>1.e-12)
      {
        int icount=0;
        do 
        {
          ddvector delta=J*h-d;
          h-=solve(J,delta);
          cd=norm(J*h-d);
          //cout <<  " cd = " << setscientific() << setprecision(5) 
          //   << cd << endl;
          icount++;
        }
        while(cd>1.e-12 && icount<3);
      }
      int mflag=0;
      int miter=0;
      do
      {
        q-=h;
        miter++;
        
        int fi;
        for (fi=1;fi<=nf;fi++)
        {
          mflag=0;
          if (q(fi)<0.0) 
          {
            q+=h;
            h/=2.0;
            mflag=1;
            break;
          }
        }
      }
      while(mflag && miter<20);
      
      icount++;
    }
    while(1);
    if (conv_flag==0)
    {
      cerr << "conv_flg error" << endl;
      ad_exit(1);
    }
  }

#if !defined(USE_DD_NOT)
void dbgprint(const ddvector & v)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    cout << my_double(v(i)) << " " << endl;
  }
}

void dbgprint(dd_real v)
{
  cout << my_double(v) << " " << endl;
}

void ddbgprint(const ddvector & v)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    cout << my_double(v(i)) << " ";
  }
  cout <<"ddnrc3.cpp " << endl;
}
void dbgprint(const ddmatrix & v)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    dbgprint(v(i));
  }
}

void ddbgprint(const ddmatrix & v)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    ddbgprint(v(i));
  }
}

void mp(const ddvector & v)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    cout << my_double(v(i)) << " " << endl;
  }
}

void mp(const ddvector & v,int i)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  cout << my_double(v(i)) << " " << endl;
}

void mp(dd_real v)
{
  cout << my_double(v) << " " << endl;
}

void mp(const ddmatrix & v,int i)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  mp(v(i));
}


void mp(const ddmatrix & v)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    ddbgprint(v(i));
  }
}
#endif
