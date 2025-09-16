/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


#define HOME_VERSION
#include "all.hpp"

#include "variable.hpp"
    //int relaxation_sequence_control(void);
  void hess_routine(independent_variables& x,dvar_len_fish_stock_history& fsh,
    int nvar);
  void dep_grad_routine(independent_variables& x,
    dvar_len_fish_stock_history& fsh,int nvar,int grad_switch);

  void minimizing_routine(independent_variables& x,dvar_len_fish_stock_history& fsh,
    int nvar,int& quit_flag, int& hang_flag, MY_DOUBLE_TYPE& maxg, MY_DOUBLE_TYPE& crit,
    int& current_ifn);
  void make_covariance_report_for_projections(
    const dvar_fish_stock_history& fsh);
  void calculate_variance_by_delta_method_noeff(const char * _s);
  void calculate_variance_by_delta_method(const char * _s);
  void check_hessian_pd(const char * _s);

//extern int _sd_flag;
void write_inverse_hessian(void);
void make_stddev_report(void);
void make_stddev_report_noeff(dvar_len_fish_stock_history& fsh);
void make_correlation_report(dvar_len_fish_stock_history& fsh);
void stitch_parallel_hessian(void);

extern int dertest_flag;

#if defined(_WIN32)
#  if defined(close)
#    undef close
#  endif
#endif

  void fitting_procedure(dvar_len_fish_stock_history& fsh,int& quit_flag)
  {
     {
      // calculate the number of active parameters
      int nvar=fsh.nvcal();
      independent_variables x(1,nvar);

      fsh.construct_indep_vars(nvar);
      
      // get the independent variable scaling vector derived from the positivized Hessian
      if (fsh.parest_flags(385)==1)
      {
        fsh.get_pos_hess_scale(nvar);
      }
      else if (fsh.parest_flags(385)==2)
      {
        fsh.read_psqrt_hess(nvar);
      }
      
      dvector g(1,nvar);
      // transform the parameters
      fsh.transform_in();
      // Get the inital values for "fundamental" parameters
      fsh.get_initial_parameter_values();
      // Get inital x values for active parameters
      fsh.xinit(x);
      if (dertest_flag)
      {
        x(1)=0.9;
        x(2)=1.1;
        x(3)=2.1;
      }
      // now do the minimization
      int svalue=fsh.parest_flags(145);
      int current_ifn=0; 
      //if (_sd_flag) svalue=_sd_flag;
      switch(svalue)
      {
	case 0:
	{
          {
            int hang_flag;
            MY_DOUBLE_TYPE maxg;
            MY_DOUBLE_TYPE ccrit;
            fsh.xinit(x);
            do  
            {
              hang_flag=0;
              maxg=0.0;
              ccrit=0.0;
              test_the_pointer();
              mytestcin();
              minimizing_routine(x,fsh,nvar,quit_flag,hang_flag,
                maxg,ccrit,current_ifn);
            }
            while (hang_flag && (maxg > 10.0*ccrit));
          }
	  break;
	}
	case 1:
	{
	  hess_routine(x,fsh,nvar);
          quit_flag=1;
	  exit(1);
	  break;
	}
	case 2:
	{
          dep_grad_routine(x,fsh,nvar,1);
          quit_flag=1;
	  exit(1);
	  break;
        }
        case 3:
	{
	  hess_routine(x,fsh,nvar);
          dep_grad_routine(x,fsh,nvar,1);
          quit_flag=1;
	  exit(1);
	  break;
        }
        case 4:
        {
          if (!sum(column(fsh.fish_flags,55)))
          {
            calculate_variance_by_delta_method(ad_root);
            //make_stddev_report();
          }
          else
          {
            calculate_variance_by_delta_method_noeff(ad_root);
            //make_stddev_report_noeff(fsh);
          }
          quit_flag=1;
	  exit(1);
          break; 
        }
        case 5:
        {
          make_correlation_report(fsh);
          quit_flag=1;
	  exit(1);
          break; 
        }
        case 6:
        {
          write_inverse_hessian();
          quit_flag=1;
	  exit(1);
          break; 
        }
	case 7:
	{
          dep_grad_routine(x,fsh,nvar,2);
          quit_flag=1;
	  exit(1);
	  break;
        }
        case 8:
        {
          /*
          dvector g(1,nvar);
          gradient_structure::set_YES_DERIVATIVES();
          cout << " Calling fcomp to get BH devs " << endl;
          fcomp(fsh,dvar_vector(x),nvar,0,g,0);
          cout << " Finished fcomp to get BH devs " << endl;
          */
          make_covariance_report_for_projections(fsh);
          quit_flag=1;
	  exit(1);
          break; 
        }
        case 9:
        {
          check_hessian_pd(ad_root);
          break;
        }
        case 10:
        {
          correlation_report(ad_root,fsh);
          break;
        }
         case 11:
        {
          stitch_parallel_hessian();
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
              minimizing_routine(x,fsh,nvar,quit_flag,hang_flag,
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
              minimizing_routine(x,fsh,nvar,quit_flag,hang_flag,
                maxg,ccrit,current_ifn);
            }
            while (hang_flag && (maxg > 10.0*ccrit));
          }
          dep_grad_routine(x,fsh,nvar,1);
          quit_flag=1;
	  exit(1);
	  break;
        }

        default:
        {
          cerr << "Illegal switch value for fsh.parest_flag(145)" << endl; 
          exit(1);
        } 
      }
      gradient_structure::set_NO_DERIVATIVES();
      dvariable ffpen=0.0;
      fsh.do_everything_calc(ffpen);
      fsh.transform_out();
     }
  }

void dvar_fish_stock_history::get_pos_hess_scale(int nvar)
{
  adstring pfname1= ad_root+"_pos_hess_diag";
  ifstream ifs(pfname1);
  int nvar1;
  ifs >> nvar1;
  if (!ifs)
  {
    cerr << "Error reading nvar1 from file " << pfname1 << endl;
    ad_exit(1);
  }
  if (nvar != nvar1)
  {
    cerr << "error reading number of variables from file " << pfname1 << " should be " << nvar
         << " for this model -- the file has " << nvar1 << endl;
    ad_exit(1);
  }
  if (!allocated(pos_hess_scale))
  {
    pos_hess_scale.allocate(1,nvar);
  }
  else
  {
    pos_hess_scale.deallocate();
    pos_hess_scale.allocate(1,nvar);
  }
  ifs >> pos_hess_scale;
  if (!ifs)
  {
    cerr << "Error reading pos_hess_scale from file " << pfname1 << endl;
    ad_exit(1);
  }
}

void dvar_fish_stock_history::read_psqrt_hess(int nvar)
{
  adstring pfname3= ad_root+"_pos_sqrt_hess";
  adstring pfname4= ad_root+"_pos_sqrt_hess_inv";
  uistream p1(pfname3);
  int nvar1;
  p1 >> nvar1;
  if (!p1)
  {
    cerr << "Error reading nvar1 from file " << pfname3 << endl;
    ad_exit(1);
  }
  if (allocated(psqrt_hess))
  {
    psqrt_hess.deallocate();
  }
  psqrt_hess.allocate(1,nvar,1,nvar);
  p1 >> psqrt_hess;
  p1.close();

  uistream p2(pfname4);
  p2 >> nvar1;
  if (!p2)
  {
    cerr << "Error reading nvar1 from file " << pfname4 << endl;
    ad_exit(1);
  }
  if (allocated(psqrt_hess_inv))
  {
    psqrt_hess_inv.deallocate();
  }
  psqrt_hess_inv.allocate(1,nvar,1,nvar);

  p2 >> psqrt_hess_inv;
  p2.close();
}

void dvar_fish_stock_history::construct_indep_vars(int nvar)
{
  if (!allocated(indep_var))
  {
    indep_var.allocate(1,nvar);
    indep_var.initialize();
  }
  else
  {
    indep_var.deallocate();
    indep_var.allocate(1,nvar);
    indep_var.initialize();
  }
  if (!allocated(indep_var_lo))
  {
    indep_var_lo.allocate(1,nvar);
    indep_var_lo.initialize();
  }
  else
  {
    indep_var_lo.deallocate();
    indep_var_lo.allocate(1,nvar);
    indep_var_lo.initialize();
  }
  if (!allocated(indep_var_hi))
  {
    indep_var_hi.allocate(1,nvar);
    indep_var_hi.initialize();
  }
  else
  {
    indep_var_hi.deallocate();
    indep_var_hi.allocate(1,nvar);
    indep_var_hi.initialize();
  }
}

extern adstring full_datafile_path;
extern adstring full_input_parfile_path;
extern adstring full_output_parfile_path;
extern adstring directory_path;
extern adstring ad_root;

void getpaths(adstring& current_path,adstring& root);

void write_inverse_hessian(void)
{
  int nvar;
  int i,j;

  adstring current_path;
  adstring root;

  //getpaths(current_path,root);

  adstring fname= ad_root+".hes";

  uistream ifs(fname);
  if (!ifs)
  {
    cerr << "Error trying to open file " << fname << endl;
  }
  ifs >> nvar;
  if (!ifs)
  {
    cerr << "Error reading nvar from " << fname << endl;
  }
  dmatrix hessinv(1,nvar,1,nvar);
  dvector row(1,nvar);
 cout<<"reading "<<nvar<<" X "<<nvar<<" hessian"<<endl;
  ifs >> hessinv;
 cout<<"hessian read"<<endl;
  if (!ifs)
  {
    cerr << "Error reading hessian from " << fname << endl;
  }
  // check for zero row

  for (i=1;i<=nvar;i++)
  {
    for (int j=1;j<=nvar;j++)
    {
      if (hessinv(i,j)!=0) row(i)=1;
    }
  }
  cout<<"zero check 1  ";
  for (i=1;i<=nvar;i++)
  {
    if (row(i)==0)
    {
      cout << " Error there was a zero row " << endl;
      cout << " in row " << i << endl;
      exit(1);
    }
  }
  cout<<"zero check 2  "<<endl;
  {
    hessinv=inv(hessinv);
  }

  cout<< "writing inverse hessian"<<endl;

  fname= ad_root+".hesinv";
  ofstream ofs(fname);
  if (!ofs)
  {
    cerr << "Error trying to open file " << fname << endl;
  }

  for(i=1;i<=nvar;i++) 
  {
    for(j=1;j<=nvar;j++) 
    {
       ofs <<"  " << setw(9) << setprecision(4) << hessinv(i,j);
    }
    ofs << endl;
  }
  
}


void make_stddev_report(void)
{
  int nvar;
  int nvar1;
  int num_deps;
  int i;

  adstring current_path;
  adstring root;

  //getpaths(current_path,root);

  adstring fname= ad_root+".hes";

  uistream ifs(fname);
  if (!ifs)
  {
    cerr << "Error trying to open file " << fname << endl;
  }
  ifs >> nvar;
  if (!ifs)
  {
    cerr << "Error reading nvar from " << fname << endl;
  }
  dmatrix hessinv(1,nvar,1,nvar);
  dvector row(1,nvar);

  ifs >> hessinv;
  if (!ifs)
  {
    cerr << "Error reading hessian from " << fname << endl;
  }
  // check for zero row
  for (i=1;i<=nvar;i++)
  {
    for (int j=1;j<=nvar;j++)
    {
      if (hessinv(i,j)!=0) row(i)=1;
    }
  }

  for (i=1;i<=nvar;i++)
  {
    if (row(i)==0)
    {
      cout << " Error there was a zero row " << endl;
      cout << " in row " << i << endl;
      exit(1);
    }
  }
  {
    hessinv=inv(hessinv);
  }

  fname=ad_root+".dep";
  uistream ifs1(fname);
  if (!ifs1)
  {
    cerr << "Error trying to open file " << fname << endl;
  }
  ifs1 >> nvar1 >> num_deps;
  if (!ifs1)
  {
    cerr << "Error reading nvar from depgrad.rpt" << endl;
  }
  if (nvar!=nvar1)
  {
    cerr << "number of independent variables is not consistent"
	    " between files" << "hessinv.rpt and depgrad.rpt" << endl;
    exit(1);
  }
  dvector depgrad(1,nvar);
  fname=ad_root+".var";
  ofstream ofs(fname);
  ifstream ifs2("deplabel.tmp");
  line_adstring deplabel;
  MY_DOUBLE_TYPE depvalue=0.0;

  print_identifier_stuff(ofs);
  for (i=1;i<=num_deps;i++)
  {
    ifs1 >> depgrad;
    ifs2 >> depvalue;
    ifs2 >> deplabel;
    if (!ifs1)
    {
      cerr << "Error reading gradient for dependent variable " << i
	   << "  from depgrad.rpt" << endl;
      exit(1);
    }
    MY_DOUBLE_TYPE variance=depgrad*hessinv*depgrad;



    ofs << setw(3) << i; 
    if (variance>=0)
    {
      ofs << "  "  << setw(8) << setprecision(3) << sqrt(variance);
    }
    else
    {
      cout << "variance = " << variance << endl;
      ofs << "  ********";
    }
    ofs << "  "  << setw(9) << setprecision(4) 
         << depvalue << "  "  << deplabel << endl;
  }
}

#undef HOME_VERSION
