/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
 #define HOME_VERSION
#include "all.hpp"
#if defined(WIN32)
#  include <windows.h>
#endif
#if defined(USE_ADPVM)
extern "C" {
#include <pvm3.h>
}
#endif //#if defined(USE_ADPVM)
#include <admodel.h>

extern ofstream * clogf;

  
dvar_vector get_x_from_master(void);

int mfget_int_from_master(void);
void mfsend_int_to_slaves(int n);

#include "variable.hpp"
    //int relaxation_sequence_control(void);
  void hess_routine(independent_variables& x,dvar_len_fish_stock_history& fsh,
    int nvar);
  void dep_grad_routine(independent_variables& x,
    dvar_len_fish_stock_history& fsh,int nvar);
  void slave_minimizing_routine(independent_variables& x,
    dvar_len_fish_stock_history& fsh,int nvar);

//  extern int _sd_flag;
void make_stddev_report(void);
void make_stddev_report_noeff(dvar_len_fish_stock_history& fsh);
void make_correlation_report(void);
void start_catchdev_part1(dvar_fish_stock_history & fsh );

dvar_vector get_x_from_master(void);

#if defined(USE_ADPVM)
  void slave_fitting_procedure(dvar_len_fish_stock_history& fsh)
  {
     *clogf << "In slave fitting procedure" << endl;
     {
      // calculate the number of active parameters
      int nvar=fsh.nvcal();
      independent_variables x(1,nvar);
      dvector g(1,nvar);
      // transform the parameters
      fsh.transform_in();
      // now do the minimization
      int svalue=fsh.parest_flags(145);
      //      if (_sd_flag) svalue=_sd_flag;
      switch(svalue)
      {
	case 0:
	{
          *clogf << "starting slave minimizing routine" << endl;
          slave_minimizing_routine(x,fsh,nvar);
          //cerr << "not implemented" << endl;
          //ad_exit(1); 
          //start_catchdev_part1(fsh);
          //fsh.grouped_implicit_catchability_deviations_calc_part1_thread();
     // *********************************************************
     // *********************************************************
     // !!! Dave F Dec 19 02 If you want pvm used in other calculations 
     //  remove this and modify program
     //   ad_exit(0);
     // *********************************************************
     // *********************************************************

          *clogf << "leaving slave minimizing routine" << endl;
	  break;
	}
       /*
	case 1:
	{
	  hess_routine(x,fsh,nvar);
          quit_flag=1;
	  exit(1);
	  break;
	}
	case 2:
	{
          dep_grad_routine(x,fsh,nvar);
          quit_flag=1;
	  exit(1);
	  break;
        }
        case 3:
	{
	  hess_routine(x,fsh,nvar);
          dep_grad_routine(x,fsh,nvar);
          quit_flag=1;
	  exit(1);
	  break;
        }
        case 4:
        {
          if (!sum(column(fsh.fish_flags,55)))
            make_stddev_report();
          else
            make_stddev_report_noeff(fsh);
          quit_flag=1;
	  exit(1);
          break; 
        }
        case 5:
        {
          make_correlation_report();
          quit_flag=1;
	  exit(1);
          break; 
        }
        case -1:
        {
          {
            int hang_flag;
            MY_DOUBLE_TYPE maxg;
            MY_DOUBLE_TYPE ccrit;
            do  
            {
              hang_flag=0;
              maxg=0.0;
              ccrit=0.0;
              slave_minimizing_routine(x,fsh,nvar,quit_flag,hang_flag,
                maxg,ccrit,current_ifn);
            }
            while (hang_flag && (maxg > 10.0*ccrit));
          }
	  hess_routine(x,fsh,nvar);
          quit_flag=1;
	  exit(1);
	  break;
        }
        case -2:
        {
          {
            int hang_flag;
            MY_DOUBLE_TYPE maxg;
            MY_DOUBLE_TYPE ccrit;
            do  
            {
              hang_flag=0;
              maxg=0.0;
              ccrit=0.0;
              slave_minimizing_routine(x,fsh,nvar,quit_flag,hang_flag,
                maxg,ccrit,current_ifn);
            }
            while (hang_flag && (maxg > 10.0*ccrit));
          }
          dep_grad_routine(x,fsh,nvar);
          quit_flag=1;
	  exit(1);
	  break;
        }

        default:
        {
          cerr << "Illegal switch value for fsh.parest_flag(145)" << endl; 
          exit(1);
        } 
      */
      }
      gradient_structure::set_NO_DERIVATIVES();
      dvariable ffpen=0.0;
      fsh.do_everything_calc(ffpen);
      fsh.transform_out();
     }
  }


  void slave_minimizing_routine(independent_variables& x,
    dvar_len_fish_stock_history& fsh,int nvar)
  {
    int ic=0;
    do 
    {
      //mf_pvm->mft.initialize();
      int ii;
      *clogf << "getting int from master" << endl;
      ii=mfget_int_from_master();
      if (ic<5)
        *clogf << "got int from master-- value is " << ii << endl;
      if (ii<=0) break;
      //cout <<  " get_int from master wait = " << 
       //mf_pvm->mft.report() << endl;

      dvector gbest(1,nvar);
      if (ic<2)
        *clogf << "gettinfg x from master" << endl;
      dvar_vector x=get_x_from_master();
      if (ic<2)
        *clogf << "got x from master" << endl;
      ic++;
      MY_DOUBLE_TYPE f=fcomp(fsh,x,nvar,0,gbest,0);
      slave_gradcalc();
    }
    while(1);
  }
#endif //#if defined(USE_ADPVM)
