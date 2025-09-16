/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#define HOME_VERSION
#include <fvar.hpp>
  void mytestcin(void);
  void insert_identifier_string_10(const char * _s);

// this is to get UNIX systems to use getchar

// #define UNIXKLUDGE

#if defined(__SPDLL__)
#  if !defined(linux)
#    include <windows.h>
#  endif
#include <admodel.h>
#include <s.h>
#include <newredef.h>
#endif

#if defined(__BORLANDC__)
  #include <signal.h>
#endif

int need_this_int=0;   // to ensure this gets linked in

#ifdef __ZTC__
  #include <conio.h>
#endif

#ifdef __GNUDOS__
  #define ADGETCH getch
#endif

#include <fvar.hpp>
extern int ctlc_flag;

#if defined(__TURBOC__) && defined(__linux__)
  void clrscr(void);
#endif

#if defined(__TURBOC__) && !defined(__linux__)
  #pragma hdrstop
  #include <iostream.h>
  #include <conio.h>
#endif

#if defined (__WAT32__) || defined (__MSVC32__) 
  #include <conio.h>
#endif

#if defined(__CYGWIN__)
  int kbhit(void)
  {
    return 0;
  }
#endif

#ifdef __ZTC__
  #include <iostream.hpp>
  #include <disp.h>
  #define endl "\n"
  //#define if (ad_printf) (*ad_printf) disp_if (ad_printf) (*ad_printf)
  void clrscr(void);
#endif

#ifdef __SUN__
  #define getch getchar
  void clrscr(void); //{ if (ad_printf) (*ad_printf)("\n"); }
#endif

#if defined(__GNU__) || defined(UNIXKLUDGE)
  #define getch getchar
#endif
  #include <signal.h>
  
  extern "C" void onintr(int k)
  {
    signal(SIGINT, exit_handler);
    ctlc_flag = 1;
    if (ad_printf) (*ad_printf)("\npress q to quit or c to invoke derivative checker: ");
    //return(1);
  }


#ifdef __NDPX__
  #include <iostream.hxx>
  extern "C" {
    void clrscr();
  };
#endif


#if defined (__MSVC32__)
  void __cdecl clrscr(void){;}
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
//void ijtest(int i, int j, int ij, int place)
//{
//  cout << "At location " << place
//       << ", i = " << i
//       << ", j = " << j
//       << ", ij = " << ij << "\n";
//}

MY_DOUBLE_TYPE dafsqrt( MY_DOUBLE_TYPE x );


  void tracing_message(int traceflag,const char *s);
  void tracing_message(int traceflag,const char *s,int *pn);
  void tracing_message(int traceflag,const char *s,MY_DOUBLE_TYPE *pd);
  void tracing_message(int traceflag,const char *s,MY_DOUBLE_TYPE d);

  void tracing_message(int traceflag,const char *s)
  {
    if (traceflag)
    {
      ofstream ofs;
      ofs.open("adtrace",ios::app);
      ofs << s << endl;
      ofs.close();
    }
  }
    
  void tracing_message(int traceflag,const char *s,int *pn)
  {
    if (traceflag)
    {
      ofstream ofs;
      ofs.open("adtrace",ios::app);
      ofs << s << " " << *pn << endl;
      ofs.close();
    }
  }
    
  void tracing_message(int traceflag,const char *s,MY_DOUBLE_TYPE *pd)
  {
    if (traceflag)
    {
      ofstream ofs;
      ofs.open("adtrace",ios::app);
      ofs << s << " " << *pd << endl;
      ofs.close();
    }
  }
    
  void tracing_message(int traceflag,const char *s,MY_DOUBLE_TYPE d)
  {
    if (traceflag)
    {
      ofstream ofs;
      ofs.open("adtrace",ios::app);
      ofs << s << " " << d << endl;
      ofs.close();
    }
  }
    
int log_values_switch=0;
ofstream logstream("fmin.log");

void print_values(const MY_DOUBLE_TYPE& f, const dvector & x,const dvector& g)
{
  logstream << setprecision(13) << f << endl;
  logstream << setprecision(13) << x << endl;
  logstream << setprecision(13) << g << endl;
}

extern adtimer * pfmintime=0;
void set_grad_profiler_10(const char * s);

extern int traceflag;
#pragma warn -sig
void fmm::fmin(BOR_CONST MY_DOUBLE_TYPE& _f, const dvector & _x,BOR_CONST dvector& _g)
{
  mytestcin();
  if (log_values_switch)
  {
    print_values(_f,_x,_g);
  }
  if (pfmintime==0) pfmintime=new adtimer;

  tracing_message(traceflag,"A3");
  //{
  //  ofstream ofs("log2");
   // ofs << (void*)ad_printf << endl;
   // ofs << (void*)printf << endl;
  //}
  dvector& g=(dvector&) _g;
  MY_DOUBLE_TYPE& f=(MY_DOUBLE_TYPE&) _f;
  independent_variables& x= (independent_variables&) _x;
    #ifdef DIAG
      cout << "On entry to fmin: " << *this << endl;
    #endif

  tracing_message(traceflag,"A4");
  //{
#if !defined (__MSVC32__)
    #if defined( __SUN__) && !(defined __GNU__)
      #if defined( __HP__)
        if (ireturn <= 0 )
        {
	   signal(SIGINT, &onintr);
        }
      #else
        if (ireturn <= 0 )
        {
	   signal(SIGINT, (SIG_PF)&onintr);
        }
      #endif
    #endif
 #endif

  #if defined (__MSVC32__)           // INTR
      if (ireturn <= 0 )     // INTR
      {
	 signal(SIGINT, &onintr);    // INTR
      }    // INTR
  #endif    // INTR
      

    #if defined( __GNU__) || defined (__BORLANDC__)
      if (ireturn <= 0 ) 
      {
	 signal(SIGINT, &onintr);
      }
    #endif

    #ifdef __ZTC__
      if (ireturn <= 0 ) // change to if (disp_inited == 0)
      {
        if (disp_inited == 0)
        {
          disp_open();
          disp_usebios();
        }
      }
    #endif
  mytestcin();

  tracing_message(traceflag,"A5");
      if (ireturn ==1 && dcheck_flag ==0) 
      {
        ireturn = 3;
      }
      
      if (ireturn >= 3)
      {
         derch(f, x, g, n, ireturn);
         return;
      }

  mytestcin();
      if (ireturn == 1) goto call1;
      if (ireturn == 2) goto call2;

      fbest=1.e+100;
  tracing_message(traceflag,"A6");
      if (!h) h.allocate(n);
      w.initialize();
      alpha=1.0;
      cout << " 1 alpha = " << setscientific() << setprecision(16) << alpha  << endl;

      ihflag=0;

     if (n==0) 
     { 
       cerr << "Error -- the number of active parameters"
         " fmin must be > 0\n";
       ad_exit(1);
     } 

  tracing_message(traceflag,"A7");
     if (x.indexmin() !=1) 
     { 
       cerr << "Error -- minimum valid index"
         " for independent_variables in fmin must be 1\n"
        << " it is " << x.indexmin() << "\n";
        ad_exit(1);
     } 

     if (x.size() <n) 
     { 
       cerr << "Error -- the size of the independent_variables"
        " which is " << x.size() << " must be >= " << n << "\n"
        << " the number of independent variables in fmin\n";
        ad_exit(1);
     } 


  tracing_message(traceflag,"A8");
     if (g.indexmin() !=1) 
     { 
       cerr << "Error -- minimum valid index"
         " for the gradient vector in fmin must be 1\n"
        << " it is " << g.indexmin() << "\n";
        ad_exit(1);
     } 

     if (g.size() <n) 
     { 
       cerr << "Error -- the size of the gradient vector"
        " which is " << g.size() << " must be >=\n"
        << " the number of independent variables in fmin\n";
        ad_exit(1);
     } 

  tracing_message(traceflag,"A9");

     for (i=1; i<=n; i++)
           xx.elem(i)=x.elem(i); /* save current feasible estimates in xx */

/*      do 45200 i=1,n
         xx(i)=x(i)
45200 continue   */

  tracing_message(traceflag,"A10");
      itn=0;
      icc=0;

  /*    itn=0
        icc=0  */

       for (i=1; i< 11; i++)
          funval.elem(i)=0.;


   /*  do 835 i=1,10
         funval(i)= 0.
  835 continue  */

  tracing_message(traceflag,"A11");
      ihang = 0;
      llog=1;
      np=n+1;
      n1=n-1;
      nn=n*np/2;
      is=n;
      iu=n;
      iv=n+n;
      ib=iv+n;
      iexit=0;
//      ij=nn+1;

 /*   ihang = 0
      llog=.true.
      np=n+1
      n1=n-1
      nn=n*np/2
      is=n
      iu=n
      iv=n+n
      ib=iv+n
      iexit=0
      ij=nn+1 */

  tracing_message(traceflag,"A12");
      h.elem(1,1) = 1;
      for (i=2; i<=n; i++)
      {
        for ( j=1; j<i; j++)
        {
           h.elem(i,j)=0;
        }
        h.elem(i,i)=1;
      }

  tracing_message(traceflag,"A13");
/*
      for (i=1; i<n+1; i++)
      {
        for ( j=1; j<i+1; j++)
        {
//         ij=ij-1;
//         ijtest(n-j+1,n-i+1,ij,-1);
//         h.elem(ij)=0;
           h.elem(n-j+1,n-i+1)=0;
        }
//      ijtest(n-i+1,n-i+1,ij,0);
        h.elem(n-i+1,n-i+1)=1;
      }
*/
/*    do 5 i=1,n
      do 6 j=1,i
      ij=ij-1
6     h.elem(ij)=0.
5     h.elem(ij)=1. */

//      ij=np;
      dmin=h.elem(1,1);

 /*   ij=np
      dmin=h(1)  */


      for ( i=2; i<=n; i++)
      {
         if(h.elem(i,i)<dmin)
            dmin=h.elem(i,i);
//       ijtest(i,i,ij,1);
//       ij=ij+np-i;
      }


 /*   do 16 i=2,n
      if(h(ij).lt.dmin)  then
         dmin=h(ij)
      end if
16    ij=ij+np-i   */

      if (dmin <= 0.)
         goto label7020;
      if(dfn == 0.)
         z = 0.0; //f;

  /*  if(dmin.le.0.d0)go to 7020
      if(dfn.eq.0.d0) z=f  */

  tracing_message(traceflag,"A14");
      for (i=1; i<=n; i++)
      {
        xsave.elem(i)=x.elem(i);
        x.elem(i)=xx.elem(i);
      }
      ireturn=1;
  tracing_message(traceflag,"A15");
      if (h.disk_save())
      {
        cout << "starting hessian save" << endl;
        h.save();
        cout << "finished hessian save" << endl;
      }
  tracing_message(traceflag,"A16");
      return;

  call1:
  mytestcin();
  tracing_message(traceflag,"A17");
      if (h.disk_save())
      {
        cout << "starting hessian restore" << endl;
        h.restore();
        cout << "finished hessian restore" << endl;
      }
  tracing_message(traceflag,"A18");


 /*   fcomp(f, xx, g, n);*/

      for (i=1; i<=n; i++)
      {
        x.elem(i)=xsave.elem(i);
      }
      ireturn=3;

  tracing_message(traceflag,"A19");

      {
       // cout << " newfbest = " << fbest << endl;
      }
 /*   fbest=f  */

      for ( i=1; i<=n; i++)
         gbest.elem(i)=g.elem(i);
  tracing_message(traceflag,"A20");

 /*   do 5175 i=1,n
         gbest(i)=g(i)
 5175 continue   */

      funval.elem(10) = f;
      df=dfn;
      if(dfn == 0.0)
         df = f - z;
      if(dfn < 0.0)
         df=fabs(df * f);
      if(df <= 0.0)
         df=1.;
label20:
/*  20 continue  */

      ic=0;
      iconv = 1;

 /*   funval(10)=f
      df=dfn
      if(dfn.eq.0.0)df=f-z
      if(dfn.lt.0.0)df=fabs(df*f)
      if(df.le.0.0)df=1.
20    continue
      ic=0
      iconv = 1  */

      for ( i=1; i<=9; i++)
         funval.elem(i)= funval.elem(i+1);

   /* do 840 i=1,9
         funval(i)= funval(i+1)
  840 continue */

      funval.elem(10) = f;

   /* funval(10)=f  */

      if ( itn>15 && fabs( funval.elem(1)-funval.elem(10))< min_improve )
         ihang = 1;

   /* if ( itn.gt.15 .and. fabs( funval(1)-funval(10)).lt. .0001 ) then
         ihang = 1
      end if   */

      gmax = 0;
      for ( i=1; i<=n; i++)
      {
        if(fabs(g.elem(i)) > crit) iconv = 2;
        if(fabs(g.elem(i)) > fabs(gmax) ) gmax = g.elem(i);
      }

 /*   do 1021 i=1,n
      if( fabs(g(i)) .gt. crit) iconv=2
 1021 continue   */

      if( iconv == 1 || ihang == 1)
      {
         ireturn=-1;
         goto label92;
      }

      if(iprint == 0)
         goto label21;

      if(itn == 0)
         goto label7000;

#if defined(USE_DDOUBLE)
#undef double
      if(fmod(double(itn),double(iprint)) != 0)
         goto label21;
      
#define MY_DOUBLE_TYPE MY_DOUBLE_TYPE
#else
      if(fmod(double(itn),double(iprint)) != 0)
         goto label21;
     
#endif

      if (llog) goto label7010;
#     if   !defined (__MSVC32__)  && !defined (__WAT32__) && !defined(linux) && !defined(__ADMING__)
        if (!scroll_flag) clrscr();
#     endif


label7003:
      if (iprint>0)
      {
	if (ad_printf) (*ad_printf)("%d variables; iteration %ld; function evaluation %ld\n",
              n, itn, ifn);
	if (ad_printf) (*ad_printf)("Function value %15.7le; maximum gradient component mag %13.5le\n",
#if defined(USE_DDOUBLE)
#undef double
             (MY_DOUBLE_TYPE)( f),(MY_DOUBLE_TYPE)( gmax));
#define MY_DOUBLE_TYPE MY_DOUBLE_TYPE
#elif defined(MY_REAL_DOUBLE)
              (MY_REAL_DOUBLE)f, (MY_REAL_DOUBLE)gmax);
#else
              f, gmax);
#endif
      }
label7002:
      if(iprint>0) 
      {
#if defined(MY_REAL_DOUBLE)
        {
          int mynp=min(65,n);
          for (int i=1;i<=mynp;i++)
          {
            cout << setw(5) << i <<  "  " 
                 << setscientific() << setprecision(5) << (MY_REAL_DOUBLE)(x(i)) << " " 
                 << setscientific() << setprecision(5) << (MY_REAL_DOUBLE)(g(i)); 
            if (i%4) 
            {
              cout << " : ";
            }
            else  
            {
              cout << endl;
            }
          }
        }
#else
        fmmdisp(x, g, n, this->scroll_flag,noprintx);
#endif
      }

label21 :
      itn=itn+1;
      for (i=1; i<=n; i++)
         x.elem(i)=xx.elem(i);

      w.elem(1)=-g.elem(1);

      pfmintime->get_elapsed_time_and_reset();
      for (i=2; i<=n; i++)
      {
	 i1=i-1;
         z=-g.elem(i);
         MY_DOUBLE_TYPE * pd=&(h.elem(i,1));
         MY_DOUBLE_TYPE * pw=&(w.elem(1));
         for (j=1; j<=i1; j++)
         {
            z-=*pd++ * *pw++;
         }
         w.elem(i)=z;
      }

      w.elem(is+n)=w.elem(n)/h.elem(n,n);

      {
        dvector tmp(1,n);
        tmp.initialize();
        for (i=1; i<=n1; i++)
        {
          j=i;
          MY_DOUBLE_TYPE * pd=&(h.elem(n-j+1,n-1));
          //double qd=tmp1(np-j);
          MY_DOUBLE_TYPE qd=w.elem(is+np-j);
          MY_DOUBLE_TYPE * pt=&(tmp(1));
          for (int ii=1; ii<=n1; ii++)
          {
            //tmp(ii)+=*pd-- * qd;
            *pt++ +=*pd-- * qd;
          }
          //tmp1(n-i)=w.elem(n-i)/h.elem(n-i,n-i)-tmp(i);
          w.elem(is+n-i)=w.elem(n-i)/h.elem(n-i,n-i)-tmp(i);
        }

       /*
        for (i=1; i<=n1; i++)
        {
           z=0.0;
           for (j=1; j<=i; j++)
           {
              z+=h.elem(n-j+1,n-i)*w.elem(is+np-j);
           }
           w.elem(is+n-i)=w.elem(n-i)/h.elem(n-i,n-i)-z;
           tmp2(n-i)=w.elem(n-i)/h.elem(n-i,n-i)-z;
        }
        cout << norm2(tmp1-tmp2) << endl;
      */
      }

 /*   do 25 i=1,n1
      ij=ij-1
      z=0.0
      do 26 j=1,i
      z=z+h(ij)*w(is+np-j)
26    ij=ij-1
25    w(is+n-i)=w(n-i)/h(ij)-z */

      gs=0.0;

   /* gs=0.0 */

      for (i=1; i<=n; i++)
         gs+=w.elem(is+i)*g.elem(i);

/*    do 29 i=1,n
29    gs=gs+w(is+i)*g(i) */

      iexit=2;
      if(gs >= 0.0)
         goto label92;
      gso=gs;
 
      alpha=-2.0*df/gs;
      cout << " 2 alpha = " << setscientific() << setprecision(16) << alpha  << endl;
#if !defined(NO_MY_DOUBLE_TYPE)
      if(alpha > 1.0L)
        alpha=1.0;
#else
      if(alpha > 1.0)
        alpha=1.0;
#endif
      df=f;
      tot=0.0;

 /*   iexit=2
      if(gs.ge.0.0)goto92
      gso=gs
      alpha=-2.d0*df/gs
      if(alpha.gt.1.d0)alpha=1.d0
      df=f
      tot=0.d0 */

label30:
/*30    continue   */

      iexit=3;
      if (ialph) { goto label92; }
      if( ifn >= maxfn)
      {
         maxfn_flag=1;
         goto label92;
      }
      else
      {
         maxfn_flag=0;
         iexit=1;
      }

      if(quit_flag) goto label92;
   /* iexit=3
      if(ifn.ge.maxfn)goto92
      iexit=1 */

      for (i=1; i<=n; i++)
         {
         z=alpha*w.elem(is+i);
         xx.elem(i)+=z;
         }

  /*  do 31 i=1,n
         z=alpha*w(is+i)
         xx(i)=xx(i)+z
   31 continue */



      for (i=1; i<=n; i++)
      {
        xsave.elem(i)=x.elem(i);
        gsave.elem(i)=g.elem(i);
        x.elem(i)=xx.elem(i);
        fsave = f;
      }
      fsave = f;
      ireturn=2;
      if (h.disk_save())
      {
        cout << "starting hessian save" << endl;
        h.save();
        cout << "finished hessian save" << endl;
      }
      return;

  call2:
  mytestcin();
      if (h.disk_save())
      {
        cout << "starting hessian restore" << endl;
        h.restore();
        cout << "finished hessian restore" << endl;
      }

 /*  fcomp(&fy, xx, w, n );*/

      for (i=1; i<=n; i++)
      {
        x.elem(i)=xsave.elem(i);
        w.elem(i)=g.elem(i);
        g.elem(i)=gsave.elem(i);
      }
      fy = f;
      f = fsave;
      ireturn=-1;

/*c
c   ***************************************************************
c   experimental change to deal with situation where prgram refuses
c   to make iteration even though function value has been reduced
c   */
     // cout << setprecision(17) << fbest << " " << fy << endl;
      if (fy <= fbest)
      {
        fbest=fy;
        for (i=1; i<=n; i++)
        {
          x.elem(i)=xx.elem(i);
          gbest.elem(i)=w.elem(i);
        }
      }
      

#if (defined( __SUN__) && !defined(__GNU__)) || defined(UNIXKLUDGE) || defined(linux) || defined(__MSVC32__)
       if(ctlc_flag || ifn == dcheck_flag )
#elif defined(__BORLANDC__)
       if ( kbhit() || ctlc_flag|| ifn == dcheck_flag )
#else
       if ( kbhit() || ifn == dcheck_flag )
#endif
       {
          int c=0;
          if (ifn != dcheck_flag)
          {
          #if !defined(__GNUDOS__)  || defined(UNIXKLUDGE)  || defined(linux) \
	      || defined(__CYGWIN32__) || defined(__MINGW32__) || defined (__MSVC32__)
            c = toupper(getch());
          #else
            c = toupper(getxkey());
          #endif
          }
          else
            c='C';

          if ( c == 'C')
          {
            for (i=1; i<=n; i++)
            {
              x.elem(i)=xx.elem(i);
            }
            ireturn = 3;
            dvector w1n=w(1,n);
            derch(f, x , w1n, n, ireturn);
            return;
          }
          else 
          {
            if ( c == 'Q'|| c == 'N') 
            {
              quit_flag=c;
              goto label92;
            }
            else
            {
              quit_flag=0;
            }
          }
       }

       if (quit_flag)
       {
	 if (quit_flag==1) quit_flag='Q';
	 if (quit_flag==2) quit_flag='N';
         goto label92;
       }

       icc+=1;
       if( icc >= 5)
          icc=0;

/*     icc=icc+1
#if !defined(NO_MY_DOUBLE_TYPE)
       if( icc .ge.5L) then
#else
       if( icc .ge.5) then
#endif
          icc=0
       end if*/

/* 4648 format(' value of f is',g15.6) */


      ic++;

/*    ic = ic+1 */

      if( ic >imax)
      {
         if (iprint>0)
         {
           if (ad_printf) (*ad_printf)("  ic > imax  in fminim is answer attained ?\n" );
           fmmdisp(x, g, n, this->scroll_flag,noprintx);
         }
         ihflag=1;
         ihang=1;
         goto label92;
      }

      ifn++;
      gys=0.0;

      for (i=1; i<= n; i++)
         gys+=w.elem(i)*w.elem(is+i);

      if(fy>f+fringe)
      {
         //cout << "old f = " << setprecision(12) << f << endl;
         //cout << "fy    = " << setprecision(12) << fy << endl;
         goto label40;
      }

      {
        MY_DOUBLE_TYPE ratio=gys/gso;
        if(fabs(ratio)<=0.25)
           goto label50;
        if(fabs(ratio)<=0.95 && ic > 4)
           goto label50;
      }


      if(gys>0.0)
         goto  label40;

      tot+=alpha;
      z=10.0;

      if(gs<gys)
         z=gys/(gs-gys);

      if(z>10.0)
         z=10.0;
      // XXXXXXXXXXXXXXXX
      //if (z>0.9) z=0.9;
      alpha=alpha*z;
      cout << " 3 alpha = " << setscientific() << setprecision(16) << alpha << endl;

      if (alpha == 0.)
      {
         ialph=1;
         #ifdef __ZTC__
         if (ireturn <= 0)
         {
           disp_close();
         }
         #endif
         return;
      }

      f=fy;
      gs=gys;
      goto label30;

label40:

      for (i=1;i<=n;i++)
         xx.elem(i)-=alpha*w.elem(is+i);

     /*
      cout << setprecision(15) << " f =  " <<  f  
           << " fy =  " <<  fy << endl;  
      cout << " f-fy " << f-fy
           << " (f-fy)/alpha " << (f-fy)/alpha 
           << " gs = " << gs << endl;;
      cout << " 3.0*(f-fy)/alpha+gys+gs " << 3.0*(f-fy)/alpha+gys+gs
           << " gs = " << gs << endl;;
      */
      if (alpha == 0.)
      {
        ialph=1;
        return;
      }

      z=3.0*(f-fy)/alpha+gys+gs;
      //cout << setprecision(15) << " z*z-gs*gys = " << z*z-gs*gys << endl;
      zz=dafsqrt(z*z-gs*gys);
      z=1.0-(gys+zz-z)/(2.0*zz+gys-gs);
      if (fabs(fy-1.e+95) < 1.e-66)
      {
        alpha*=.001;
        //cout << "after reduction alpha = " << setscientific
        //  << alpha << "  " << alpha << endl;
      }
      else
      {
       /*
        cout << " alpha =  " << alpha << " z = " << z 
             << " gys = "  << gys <<  " zz = " << zz  
             << " 2.0*zz+gys-gs = " << 2.0*zz+gys-gs  << endl;
       */ 
        if (z>10.0) 
        {
          cout << "large z" << z << endl;
          z=10.0;
        }
        alpha=alpha*z;
      }
      cout << " 4 alpha = " << setscientific() << setprecision(16) << alpha  << endl;

/*c
c   trap for alpha = 0
c  */

      if (alpha == 0.)
      {
         ialph=1;
        if (ialph)
        {
	   if (ad_printf) (*ad_printf)("\nFunction minimizer: Step size"
            "  too small -- ialph=1");
        }
         //#ifdef __ZTC__
         //if (ireturn <= 0)
         //{
         //  disp_close();
         //}
         //#endif
         return;
      }
      goto label30;

/*    if (alpha.eq.0.) then
         ialph=1
         return
      end if
      go to 30  */

label50:

/*50    continue  */

      alpha+=tot;
      cout << " 5 alpha = " << setscientific() << setprecision(16) << alpha  << endl;
      f=fy;
      df-=f;
      dgs=gys-gso;
      xxlink=1;
      if(dgs+alpha*gso>0.0)
         goto label52;
/*    alpha=tot+alpha
      f=fy
      df=df-f
      dgs=gys-gso
      xxlink=1
      if(dgs+alpha*gso.gt.0.d0)goto52  */

      for (i=1;i<=n;i++)
         w.elem(iu+i)=w.elem(i)-g.elem(i);

/*    do 51 i=1,n
51    w(iu+i)=w(i)-g(i)  */

      sig=1.0/(alpha*dgs);
      goto label70;

/*    sig=1.d0/(alpha*dgs)
      goto70  */

label52:

/*52    continue  */
      zz=alpha/(dgs-alpha*gso);
      z=dgs*zz-1.0;

/*    zz=alpha/(dgs-alpha*gso)
      z=dgs*zz-1.d0  */

      for (i=1;i<=n;i++)
         w.elem(iu+i)=z*g.elem(i)+w.elem(i);
/*    do 53 i=1,n
53    w(iu+i)=z*g(i)+w(i)  */

      sig=1.0/(zz*dgs*dgs);
      goto label70;
label60:
      xxlink=2;

/*    sig=1.d0/(zz*dgs*dgs)
      goto70
60    continue
      xxlink=2  */

      for (i=1;i<=n;i++)
         w.elem(iu+i)=g.elem(i);

/*    do 61 i=1,n
61    w(iu+i)=g(i)  */

      if(dgs+alpha*gso>0.0)
         goto label62;
      sig=1.0/gso;
      goto  label70;

/*    if(dgs+alpha*gso.gt.0.d0)goto62
      sig=1.d0/gso
      goto70  */

label62:

/*62    continue  */

      sig=-zz;
      goto label70;

/*    sig=-zz
      goto70  */

label65:

/*65    continue  */

      for (i=1;i<=n;i++)
         g.elem(i)=w.elem(i);

/*    do 66 i=1,n
66    g(i)=w(i)  */

      goto  label20;

/*    goto20  */

label70:

/*70    continue  */

      w.elem(iv+1)=w.elem(iu+1);

/*    w(iv+1)=w(iu+1)  */

      pfmintime->get_elapsed_time_and_reset();
      for (i=2;i<=n;i++)
      {
         i1=i-1;
         z=w.elem(iu+i);
         MY_DOUBLE_TYPE * pd=&(h.elem(i,1));
         MY_DOUBLE_TYPE * pw=&(w.elem(iv+1));
         for (j=1;j<=i1;j++)
         {
           z-=*pd++ * *pw++;
           //z-=h.elem(i,j)*w.elem(iv+j);
         }
         w.elem(iv+i)=z;
      }

/*    do 71 i=2,n
      ij=i
      i1=i-1
      z=w(iu+i)
      do 72 j=1,i1
      z-=h(ij)*w(iv+j)
72    ij=ij+n-j
71    w(iv+i)=z  */

//      ij=1;

/*    ij=1  */

      pfmintime->get_elapsed_time_and_reset();
      for (i=1;i<=n;i++)
      {
//      ijtest(i,i,ij,6);
         z=h.elem(i,i)+sig*w.elem(iv+i)*w.elem(iv+i);
//       z=h.elem(ij)+sig*w.elem(iv+i)*w.elem(iv+i);
         if(z <= 0.0)
            z=dmin;
         if(z<dmin)
            dmin=z;
         h.elem(i,i)=z;
//       h.elem(ij)=z;
         w.elem(ib+i)=w.elem(iv+i)*sig/z;
         sig-=w.elem(ib+i)*w.elem(ib+i)*z;
//       ij=ij+np-i;
       }

     /*
      for (i=1;i<=n1;i++)
      {
         i1=i+1;
         for (j=i1;j<=n;j++)
         {
            w.elem(iu+j)=w.elem(iu+j)-h.elem(j,i)*w.elem(iv+i);
            h.elem(j,i)+=w.elem(ib+i)*w.elem(iu+j);
         }
      }
      */

      cout << "START" << endl;
      for (j=2;j<=n;j++)
      {
         MY_DOUBLE_TYPE * pd=&(h.elem(j,1));
         MY_DOUBLE_TYPE * qd=&(w.elem(iu+j));
         MY_DOUBLE_TYPE * rd=&(w.elem(iv+1));
         for (i=1;i<j;i++)
         {
            *qd-=*pd * *rd++;
            *pd++ +=w.elem(ib+i)* *qd;
         }
      }
      cout << "END" << endl;



      if (xxlink == 1) goto label60;
      if (xxlink == 2) goto label65;

/*    goto(60,65),xxlink  */

label90:

/* 90    continue  */

      for (i=1;i<=n;i++)
         g.elem(i)=w.elem(i);

/*    do 91 i=1,n
91    g(i)=w(i)  */

label92:

/*92    continue

c
c   check for error and convergence conditions
c                                    */
      if (iprint>0)
      {
        if (ialph)
        {
	   if (ad_printf) (*ad_printf)("\nFunction minimizer: Step size too small -- ialph=1");
        }
        if (ihang == 1)
        {
	   if (ad_printf) (*ad_printf)("Function minimizer not making progress ... is minimum attained?\n");
#if defined(USE_DDOUBLE)
#undef double
           if (ad_printf) (*ad_printf)("Minimprove criterion = %13.5le\n",double(min_improve));
#define MY_DOUBLE_TYPE MY_DOUBLE_TYPE
#else
           if (ad_printf) (*ad_printf)("Minimprove criterion = %13.5le\n",min_improve);
#endif
        }
      }
      if(iexit == 2)
      {
        if (iprint>0)
        {
          if (ad_printf) (*ad_printf)("*** grad transpose times delta x greater >= 0\n"
           " --- convergence critera may be too strict\n");
          ireturn=-1;
        }
      }
#     if defined (__MSVC32__)  && !defined (__WAT32__)
        if (scroll_flag == 0) clrscr();
#     endif

      if (maxfn_flag == 1)
      {
        if (iprint>0)
        {
	  if (ad_printf) (*ad_printf)("Maximum number of function evaluations exceeded");
        }
      }

      if (iprint>0)
      {
        if (quit_flag == 'Q')
          if (ad_printf) (*ad_printf)("User initiated interrupt");
      }

      if(iprint == 0) goto label777;

      if (ad_printf) (*ad_printf)("\n - Final statistics:\n");
      if (ad_printf) (*ad_printf)("%d variables; iteration %ld; function evaluation %ld\n",
	      n, itn, ifn);
#if defined(USE_DDOUBLE)
#undef double
      if (ad_printf) (*ad_printf)("Function value %15.7le; maximum gradient component mag %13.5le\n",
             (MY_DOUBLE_TYPE)( f),(MY_DOUBLE_TYPE)( gmax));
      if (ad_printf) (*ad_printf)("Exit code = %ld;  converg criter %13.5le\n",iexit,double(crit));

      
#define MY_DOUBLE_TYPE MY_DOUBLE_TYPE
#else
      if (ad_printf) (*ad_printf)("Function value %15.7le; maximum gradient component mag %13.5le\n",
              f, gmax);
      if (ad_printf) (*ad_printf)("Exit code = %ld;  converg criter %13.5le\n",iexit,crit);

     
#endif

      fmmdisp(x, g, n, this->scroll_flag,noprintx);

label777:
         if (ireturn <= 0)
         #ifdef DIAG
           if (ad_printf) (*ad_printf)("Final values of h in fmin:\n");
           cout << h << "\n";
         #endif
         #ifdef __ZTC__
         {
           disp_close();
         }
         #endif
         return;
/*
c *** print headings **
c                       */
label7000:
      if (iprint>0)
      {
#     if defined (__MSVC32__)  && !defined (__WAT32__)
        if (!scroll_flag) clrscr();
#endif
        if (ad_printf) (*ad_printf)("\nInitial statistics: ");
      }
      goto label7003;

label7010:
   if (iprint>0)
   {
#     if defined (__MSVC32__)  && !defined (__WAT32__)
     if (!scroll_flag)  clrscr();
#endif
     if (ad_printf) (*ad_printf)("\nIntermediate statistics: ");
   }

   llog=0;
   goto label7003;

label7020:
   if (iprint>0)
   {
     if (ad_printf) (*ad_printf)("*** hessian not positive definite\n");
   } 
         #ifdef __ZTC__
         if (ireturn <= 0)
         {
           disp_close();
         }
         #endif
         return;
   }

   MY_DOUBLE_TYPE dafsqrt( MY_DOUBLE_TYPE x )
   {
     if (x>0)
       return(sqrt(x));
     else
       cout << "bad value in dafsqrt " << endl;
     return(0.);
   }

#undef HOME_VERSION
