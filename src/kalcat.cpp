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
  dvariable mfexp(BOR_CONST prevariable& );
  //dvar_vector mfexp(dvar_vector& );
  //dvector mfexp(dvector& );
#endif

dvariable dvar_fish_stock_history::robust_kalman_filter_for_catchability(void)
{

  dvariable var_E=0.01;
  //dvariable var_E=0.01;
  dvariable f=0.0;
  funnel_dvariable fd;
  
  
  //dvar_matrix qmu(1,num_fisheries,1,num_fish_times);
  //dvar_matrix qvar(1,num_fisheries,1,num_fish_times);
  dvar_matrix lmax(1,num_fisheries,1,num_fish_times);
  //dmatrix student(1,num_fisheries,1,num_fish_times);

  dvariable d2f;
  dvariable d0f;
  dvariable h;
  dvariable kf_contrib=0.0;

  ad_begin_funnel();
  for (int i=1;i<=num_fisheries;i++)
  {

    if(fish_flags(i,38))
      var_E=fish_flags(i,38)/100.;
    dvar_vector& q = implicit_catchability(i);
    dvariable d;
    dvariable var_q;
    dvariable  a_t;
    const MY_DOUBLE_TYPE cont=0.05;   // amount of contamination by outliers
    const MY_DOUBLE_TYPE cont1 = 1.0-cont;
    const MY_DOUBLE_TYPE var_random_walk = .003; // variance of random walk 
                                 // component of catchability

    dvariable qmu0;
    for (int nt=1;nt<=num_fish_times(i);nt++)
    {
      if (nt == 1)   // average over first 5 q's
      {
        a_t = mean(q(1,min(q.indexmax(),5)));  // do we need log of q0?
        qmu0=a_t;
        //a_t = q0(i);      // do we need log of q0?
      }
      else
        a_t = qmu(i,nt-1);

      if (nt == 1)
      {
        if(fish_flags(i,39))
          var_q = fish_flags(i,39)/100.;
        else
          var_q = 0.01;
      }
      else
        var_q = qvar(i,nt-1) + var_random_walk;

      for (int ii=1;ii<=5;ii++)
      {
        d=q(nt)-a_t;

        dvariable e=0.5*square(d)/(var_E);

        dvariable f= cont1 *exp(-e) + cont/(1.0+e);

        dvariable g= cont1*exp(-e) + cont/square(1.0+e); 

        dvariable u;
        if (nt>1)
          u=qmu(i,nt-1)-a_t;
        else
          u=qmu0-a_t;

        // this is the log of the prior dist for q_i with mean q_{i-1}
        // and variance var_q
        dvariable v=0.5*(square(u)/(var_q)+log(var_q));

        d0f=log(f) - v;

        dvariable d1f = g/f * d/var_E  + 1.0 * u/var_q; 

        d2f = square(g/f)  * -2.*e
            + (cont1*exp(-e) + 2.0*cont/cube(1.0+e))/f * 2.0 * e
            - g/f/var_E -1.0 /var_q ;

        h=d1f/d2f;

        a_t = a_t - d1f/d2f;

      }
      lmax(i,nt)=d0f;
      qmu(i,nt)=a_t;
      qvar(i,nt)= -1.0/d2f;
      qstudent(i,nt)=value(d)/sqrt(value(var_E));
      kf_contrib+=lmax(i,nt)+0.5*log(qvar(i,nt));
      //f-=lmax(i,nt)-0.5*log(qvar(i,nt));
      f+=lmax(i,nt)+0.5*log(qvar(i,nt));

    }  
  }
  cout << "Kalman filter contribution = " << -kf_contrib << endl;
  fd=f;
  f=fd;
  return f;
}

