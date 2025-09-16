/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#//define __WINDOWS__

#define USE_DD_NOT
#define HOME_VERSION
#include "all.hpp"
//#include <ddfvar.hpp>

#include "variable.hpp"
#ifdef __WINDOWS__
#include <constrea.h>
#endif

#if defined(catch)
#undef catch
#endif
//#define TRUNC_NEWTON

extern int test_adexception_flag;


void inverse_xscale(independent_variables& x,dvector& grad_scale);
dvar_vector xscale(const dvar_vector& x,const dvector& grad_scale);

void estimate_interact(dvar_len_fish_stock_history& fsh,dvar_vector& x);

//extern int num_free_obj;
extern mf_pvm_manager * mf_pvm;

MY_DOUBLE_TYPE davetest(dvar_vector& x);


dvar_vector * const_vb = NULL;
dvar_vector * const_var = NULL;
#include <signal.h>
#include <math.h>
//#include <fstream.h>
extern int ddcheck_flag;

//int global_except_switch=0;

//extern grad_stack_entry * gse;
// #if defined(AD_FPEXCEPTION)
// #include<float.h>
// class ad_float_exception
// {
// public:
//   int errtype;
//   ad_float_exception(int i) : errtype(i){}
// };
// 
// extern "C" void float_except(int k)
// {
//   throw ad_float_exception(1);
// }
// #endif

void print_indep_vars_report(dvector& gbest,dvector& idv,
			     dvector& idv_lo,dvector& idv_hi,int nvar)
{
  adstring fname= "xinit.rpt";
  ifstream ifs2(fname);
  fname= "indepvar.rpt";
  ofstream ofs3(fname);
  adstring_array indep_index(1,nvar);
  adstring_array indep_label(1,nvar);
  for (int indep=1;indep<=nvar;indep++)
  {
    ifs2 >> indep_index(indep);
    ifs2 >> indep_label(indep);
    if (!ifs2)
    {
      cerr << "Error reading indep_var label " << indep
           << "  from xinit.rpt" << endl;
           ad_exit(1);
    }
  }
  ofs3 << " Index  Var_name   Estimate L_bound  U_bound gradient" << endl;
  for (int indep=1;indep<=nvar;indep++)
  {
    ofs3 << indep_index(indep) << " " << indep_label(indep)
         << setw(4) << " "  << setprecision(8)
         << setw(4) << " "   << idv(indep) << " "
         << setw(4) << " "   << idv_lo(indep) << " "
         << setw(4) << " "   << idv_hi(indep) << setw(4) << " "
         << gbest(indep) << endl;

    if (!ifs2)
    {
      cerr << "Error reading indep_var label " << indep
           << "  from xinit.rpt" << endl;
      ad_exit(1);
    }
  }
}

void save_best_estimate(
  int & test_adexception_flag,
  dvector & x,
  dvector & xbest,
  int nvar,
  dvector& gbest,
  dvar_len_fish_stock_history & fsh,
  int ifn)
{
  cerr << "beginning parameter save" << endl;
  //save current best parameter estimates
  test_adexception_flag=0;
  x=xbest;
  gradient_structure::set_NO_DERIVATIVES();
  dvariable f=fcomp(fsh,dvar_vector(x),nvar,0,gbest,0);
  par_ofstream ofs("crash_save");
  ofs << fsh << endl;
  ofstream ofs1("crash_file");
  ofs1 << ifn << endl;
  cerr << "exiting after choleski error and parameter save" << endl;
  ad_exit(1);
}
void save_best_estimate_only(
  int iter,
  dvector & xbest,
  int nvar,
  dvector& gbest,
  dvar_len_fish_stock_history & fsh)
{
  if (fsh.parest_flags(196))
  {
//    ofstream ofs("xbest_values",ios::app);
    ofstream ofs(ad_root+".xbsamples",ios::app);    
    ofs  << setprecision(10) << xbest << endl;
  }
  cout << "beginning parameter save" << endl;
  //save current best parameter estimates
  dvector xtmp(xbest.indexmin(),xbest.indexmax());
  xtmp=xbest;
  gradient_structure::set_NO_DERIVATIVES();
  dvariable f=fcomp(fsh,dvar_vector(xtmp),nvar,0,gbest,0);
  fsh.no_lagrangian_update=1;
  par_ofstream ofs("routine_save");
  ofs << fsh << endl;
  ofs << "# Objective function value" << endl << "   " << -value(f) 
     << endl;
  ofs << "# The number of parameters" << endl << "   " << nvar << endl;
  if (fsh.num_tag_releases)
    ofs << "# Likelihood component for tags " <<fsh.likecomponent[0]<< endl;
  if (!fsh.age_flags(52))
  {
    ofs << "# Maximum magnitude gradient value " << endl  
        << max(sqr(square(gbest))) << endl;
  }
  else
  {
    ofs << "# Maximum magnitude gradient value " << endl  
        << "0.00"  << endl;
  }
  fsh.no_lagrangian_update=0;
  cout << "did routine save" << endl;
  gradient_structure::set_YES_DERIVATIVES();
}


void mfsend_x_to_slaves(const dvar_vector&  x);
/*
{
  // *********  begin variable send block  *************
  for (int i=1; i<=ad_nhost; i++)
  {
    int bufid = adpvm_master_vinitsend( PvmDataDefault );
    adpvm_pack(x);
    adpvm_master_vsend((*ad_stid)(i));
  }
  // *********  end variable send block  *************
}
*/
extern gradient_structure * pgs;
extern int  gss;
extern int stupid_print_switch;


class choleski_exception: public exception
{
  virtual const char* what() const throw()
  {
    return "Choleski exception happened";
  }
};

choleski_exception myex;
extern "C" void adfloat_except(int k)
{
  throw myex;
}


void minimizing_routine(independent_variables& x,dvar_len_fish_stock_history& fsh,
  int nvar,int& quit_flag, int& hang_flag, MY_DOUBLE_TYPE& maxg, MY_DOUBLE_TYPE& ccrit,
  int& current_ifn)
{
    // GGGGGGGGGG
    signal(SIGINT,&adfloat_except);
    //test_adexception_flag=0;
    ofstream ofssg("testgrad");
    dvector xbest(1,nvar);
    xbest=x;
  #ifdef __WINDOWS__
    constream win1;
    win1.window(1,1,40,20);
    win1.clrscr();
  #endif
  //cout <<"newl3, fsh.parest_flags(1): " <<  fsh.parest_flags(1) << endl;
  const_vb = &(fsh.vb_coff);
  const_var = &(fsh.var_coff);
  if (fsh.parest_flags(197) || fsh.parest_flags(1)==0 )
  { 
    if (nvar==0) nvar=1;
    dvariable f=0;
    dvector gbest(1,nvar);
    gbest.fill_seqadd(1.e+50,0.);
    gradient_structure::set_NO_DERIVATIVES();
    f=fcomp(fsh,dvar_vector(x),nvar,1,gbest,0);
    quit_flag=1;
    return;
  } 
  int trunc_newton_flag=fsh.parest_flags(351);
  MY_DOUBLE_TYPE angle_bound=fsh.parest_flags(352)/100.;
  if (fsh.parest_flags(352)==0)
  {
    angle_bound=0.25;
  }
  else if (angle_bound<0.02)
  {
    angle_bound=0.10;
  }
  else if (angle_bound>0.9)
  {
    angle_bound=0.90;
  }

  if (trunc_newton_flag==1)
  {
    // trunc newton
    int nterms=7;
    if (fsh.parest_flags(192)) nterms=fsh.parest_flags(192);
    //fmmt fmc(nvar,nterms);
    fmmt1 fmc(nvar,nterms,angle_bound);
    if (fsh.parest_flags(2)==0)
      fmc.iprint=5;
    else
      fmc.iprint=fsh.parest_flags(2);

    int old_itn=0;
    //  Rationalize meaning of convergence flag -- P.K. 7/17/2004 
    //  Default is now 1.0 instead of 1.e-6
    fmc.crit=pow(10.,double(fsh.parest_flags(50)));

    int nvar_check; 
    dvector grad_scale(1,nvar);
    if (fsh.parest_flags(353)==0)
    {
      gradient_structure::set_YES_DERIVATIVES();
    }
    else
    {
      gradient_structure::set_NO_DERIVATIVES();
    }
    if (fsh.parest_flags(152)>0)
    {
      ifstream ifs("gradient.rpt");
      ifs >> nvar_check;
      if (nvar_check !=nvar)
      {
        cerr << "Number of parameters in gradient scaling file is"
                " not equal to the number" << endl << 
                " of active parameters -- scaling is not enabled" << endl;
        fsh.parest_flags(152)=0;
      }
      else
      {
        ifs >> grad_scale;
        grad_scale= fabs(grad_scale);
      }
    }
    if (fsh.parest_flags(152)>0)
    {
      inverse_xscale(x,grad_scale);
    }
    fmc.scroll_flag=0;
    fmc.imax=40;
    fmc.min_improve=0.;
    fmc.maxfn=fsh.parest_flags(1);
    //cout << " fmc.maxfn = "  << fmc.maxfn  << endl;
    fmc.ifn=current_ifn;
    current_ifn=0;
    int it=fsh.parest_flags(12)+fsh.parest_flags(13)+fsh.parest_flags(14)+
      fsh.parest_flags(15)+fsh.parest_flags(16);
    fsh.parest_flags(143)=0;
    fsh.parest_flags(144)=0;
    MY_DOUBLE_TYPE f;
    MY_DOUBLE_TYPE fbest=1.e+50;;
    dvector g(1,nvar);
    g.initialize();
    dvector gbest(1,nvar);
    greport("very beginning ");
    int do_fmin_flag=0;

    gbest.fill_seqadd(1.e+50,0.);
    {
      if (fsh.age_flags(52)==0)
      {
#if defined(AD_FPEXCEPTION)
        //int oldexception=_control87(0,MCW_EM);
        //_control87(0, 0);
        signal(SIGFPE, &adfloat_except);
#endif
        if (fsh.parest_flags(222)>0)
        {
          fmc.dcheck_flag=fsh.parest_flags(222)-1;
        }
        adtimer xdt;
        while (fmc.ireturn>=0)
        {
          //global_except_switch++;
          int badflag=0;
          xdt.get_elapsed_time_and_reset();
#if defined(AD_FPEXCEPTION)
          try 
#endif
          {
            if (do_fmin_flag==0)
            {
              fmc.fmin(f,x,g);
            }
            int itmp=fmc.itn;
            int itmp1=10;
            if (fsh.parest_flags(196))
            {
              itmp1=fsh.parest_flags(196);
            }
               
            if (fmc.itn>3 && fmc.itn%itmp1==0)
            {
              save_best_estimate_only(itmp,xbest,nvar,gbest,fsh);
            }
          }
#if defined(AD_FPEXCEPTION)
          catch (choleski_exception & myexp)
          {
            save_best_estimate(test_adexception_flag,x,xbest,nvar,
              gbest,fsh,fmc.ifn);
          }
#endif
          cout << "FMIN time = " << xdt.get_elapsed_time_and_reset()/1000. << endl;
          {
            int x;
            //cout << "Enter x" << endl;
            //cin >> x;
          }
#if defined(USE_ADPVM)
          if (mf_pvm->pvm_switch==1)
          {
            mfsend_int_to_slaves(fmc.ireturn);
          }
#endif //#if defined(USE_ADPVM)
          if (fmc.ihang)
          {
            hang_flag=fmc.ihang;
            maxg=max(g);
            ccrit=fmc.crit;
            current_ifn=fmc.ifn;
          }
          if (fmc.ireturn>0)
          {
#if defined(AD_FPEXCEPTION)
            try 
#endif
            {
              /*
              int globalstop;
              ifstream ifsg("globalstop");
              ifsg >> globalstop;
              ifsg.close();
              if (fmc.itn==globalstop)
              {
                ofstream ofsg("globalstop");
                globalstop+=2;
                ofsg << globalstop;
                adfloat_except(1);
              }
              */
              if (fsh.parest_flags(152)>0)
              {
                dvar_vector y=xscale(dvar_vector(x),grad_scale);
                f=fcomp(fsh,y,nvar,0,gbest,0);
                fsh.parest_flags(144)=0;
              }
              else
              {
#if defined(USE_ADPVM)
                if (mf_pvm->pvm_switch==1)
                {
                  dvar_vector vx=dvar_vector(x);
                  mfsend_x_to_slaves(vx);
                  f=fcomp(fsh,vx,nvar,0,gbest,0);
                }
                else
#endif //#if defined(USE_ADPVM)
                {
                  //fsh.do_fishery_projections_flag=1;
                  MY_DOUBLE_TYPE zz=-1;
                  signal(SIGFPE, &adfloat_except);
                  //cout << sqrt(zz) << endl;
                  
                  f=fcomp(fsh,dvar_vector(x),nvar,0,gbest,0);
//                  cout << " LLL f = " << f << endl;
       
                  if(fmc.ireturn ==1 && fsh.parest_flags(1)==1 && sum(fsh.data_fish_flags(2)))   //NMD 20Jan2012
                  {
                    fsh.save_bh_variance = fsh.bh_variance;
                  }                                                                              //NMD 20Jan2012
                  fsh.do_fishery_projections_flag=0;
                }
                fsh.parest_flags(144)=0;
              }
              badflag=1;
              cout << "==============================================" << endl;
	      cout << " f eval " << fmc.ifn << endl;
              char * cptr = (char *) &f;
              for (int iw=0;iw<8;iw++)
              {
                if ( int(*(cptr+iw)) != iw ) badflag=0;
              }
            }
#if defined(AD_FPEXCEPTION)
          
            catch (choleski_exception & myexp)
            {
              save_best_estimate(test_adexception_flag,x,xbest,nvar,
                gbest,fsh,fmc.ifn);
            }
            
            catch (myexception & adfe)
            {
              //save current best parameter estimates
              test_adexception_flag=0;
              x=xbest;
              f=fcomp(fsh,dvar_vector(x),nvar,0,gbest,0);
              par_ofstream ofs("crash_save");
              ofs << fsh << endl;
              cerr << " exiting after choleski error and parameter save" << endl;
              ad_exit(1);
            }
            
            catch (ad_float_exception & adfe)
            {
              badflag=1;
              reset_gradient_stack();

              //save current best parameter estimates
              x=xbest;
              test_adexception_flag=0;
              f=fcomp(fsh,dvar_vector(x),nvar,0,gbest,0);
              par_ofstream ofs("crash_save");
              ofs << fsh << endl;
              
              cerr << " exiting after FPE and parameter save" << endl;
              ad_exit(1);
            }
#endif
            if (badflag)
            {
              f=1.e+95;
            }
#if defined(AD_FPEXCEPTION)
            int in_gradcalc_flag=1;
            try 
#endif
            {
              cout << "calling gradcalc" << endl;
              if (fsh.parest_flags(353)==0)
              {
                gradcalc(nvar,g);
              }
              in_gradcalc_flag=0;
            }
#if defined(AD_FPEXCEPTION)
            catch (choleski_exception & adfe)
            {
              if (in_gradcalc_flag && gradient_structure::get_save_var_flag())
              {
                gradient_structure::restore_arrays();
                gradient_structure::restore_variables();
              }
              save_best_estimate(test_adexception_flag,x,xbest,nvar,
                gbest,fsh,fmc.ifn);
            }
            catch (myexception & adfe)
            {
              //save current best parameter estimates
              test_adexception_flag=0;
              x=xbest;
              f=fcomp(fsh,dvar_vector(x),nvar,0,gbest,0);
              par_ofstream ofs("crash_save");
              ofs << fsh << endl;
              cerr << " exiting after choleski error and parameter save" << endl;
              ad_exit(1);
            }
#endif
            

            if (0)
            {
              random_number_generator rng(435);
              dvector eps(1,nvar);
              eps.fill_randn(rng);
              MY_DOUBLE_TYPE gn=norm(g);
              dvector sg=g/gn;
              
              ofstream ofs("feasible");
              // test for feasible direction
              gradient_structure::set_NO_DERIVATIVES();
              MY_DOUBLE_TYPE fc=f;
              ofs << "current f = " << setscientific() << setprecision(16) << f << " " << fc-f << endl;
              dvector y=x+1.e-6*g;
              MY_DOUBLE_TYPE f1=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
              ofs << "1.e-6 f = " << setscientific() << setprecision(16) << f << " " << fc-f1 << endl;
              /*
              y=x+1.e-9*g;
              MY_DOUBLE_TYPE f2=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
              ofs << "1.e-9 f = " << setscientific() << setprecision(16) << f << " " << fc-f2 << endl;
              y=x+1.e-10*g;
              MY_DOUBLE_TYPE f3=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
              ofs << "1.e-10 f = " << setscientific() << setprecision(16) << f << " " << fc-f3 << endl;
              y=x+1.e-12*g;
              MY_DOUBLE_TYPE f4=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
              ofs << "1.e-12 f = " << setscientific() << setprecision(16) << f << " " << fc-f4 << endl;
              y=x+1.e-13*g;
              MY_DOUBLE_TYPE f5=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
              ofs << "1.e-13 f = " << setscientific() << setprecision(16) << f << " " << fc-f5 << endl;
              */

              for (int i=1;i<=20;i++)
              {
                eps.fill_randn(rng);
                dvector newg=(sg+.05*eps/norm(eps));
                newg/=norm(newg);
                newg*=gn;
  
                y=x+1.e-9*newg;
                MY_DOUBLE_TYPE f8=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
                ofs << "1.e-9 f = " << setscientific() << setprecision(16) << f << " " << fc-f8 << endl;
                y=x+1.e-10*newg;
                MY_DOUBLE_TYPE f9=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
                ofs << "1.e-10 f = " << setscientific() << setprecision(16) << f << " " << fc-f9 << endl;
                y=x+1.e-11*newg;
                MY_DOUBLE_TYPE f10=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
                ofs << "1.e-11 f = " << setscientific() << setprecision(16) << f << " " << fc-f10 << endl;
                y=x+1.e-12*newg;
                MY_DOUBLE_TYPE f11=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
                ofs << "1.e-12 f = " << setscientific() << setprecision(16) << f << " " << fc-f11 << endl;
              }
              ad_exit(1);
            }
       
            cout << "finished gradcalc" << endl;
//            ofssg << g << endl;   //NMD_11May2014 - turned off output
            //cout << "heap 2 = " << heapcheck() << endl;
            //gse->dep_addr=0;
            if (f < fbest)
            {
              fbest=f;
              gbest=g;
              xbest=x;
            }
          }
          if (fsh.parest_flags(333) 
            && 
             (!((fmc.itn)%(fsh.parest_flags(333)))
              || fmc.itn==1) && fmc.itn>old_itn)
          {
            old_itn=fmc.itn;
            //save current best parameter estimates
            int test_adexception_flag_save=test_adexception_flag;
            test_adexception_flag=0;
            gradient_structure::set_NO_DERIVATIVES();
            dvector gtemp(1,nvar);
            MY_DOUBLE_TYPE ftemp=fcomp(fsh,dvar_vector(xbest),nvar,0,gtemp,0);
            int suffix=0;
            if (fmc.itn==1)
              suffix=1;
            else
              suffix=fmc.itn;
            adstring file_name= adstring("intermediate_results_") 
              + str((const int)(suffix));
            par_ofstream ofs(file_name);
            ofs << fsh << endl;
            ofs << "# Objective function value" << endl << "   " << -ftemp 
                << endl;
            ofs << "# The number of parameters" << endl << "   " << nvar << endl;
            adstring command="cp " + file_name + " intermediate_results";
            if (fsh.num_tag_releases)
              ofs << "# Likelihood component for tags " <<fsh.likecomponent[0]<< endl;
            test_adexception_flag=test_adexception_flag_save;
            system(command);
            adstring command1="echo " + str((const int)fmc.itn) +  " > success_itn";
            system(command1);
       
            adstring command2="echo " + str((const int)fmc.ifn) +  " > success_ifn";
            system(command1);
            gradient_structure::set_YES_DERIVATIVES();
          }
        }
        x=xbest;

#if defined(AD_FPEXCEPTION)
        signal(SIGFPE, SIG_DFL);
#endif
      }
     // *********************************************************
     // *********************************************************
     // !!! Dave F Dec 19 02 If you want pvm used in other calculations 
     //  remove this and modify program
        mf_pvm->pvm_save_switch=mf_pvm->pvm_switch;
        mf_pvm->pvm_switch=0;
     // *********************************************************
     // *********************************************************
     // if (!hang_flag)
      {
        //cout << endl << "Calling printing fcomp" <<endl;
        gradient_structure::set_NO_DERIVATIVES();
    
        if (fsh.parest_flags(152)>0)
        {
          dvar_vector y=xscale(dvar_vector(x),grad_scale);
          f=fcomp(fsh,y,nvar,1,gbest,0);
        }
        else
        {
          if (sum(fsh.data_fish_flags(2)))
          {
            fsh.do_fishery_projections_flag=1;
          }
          int num_sims=1;
          if (fsh.do_fishery_projections_flag==1 && fsh.age_flags(20)>0) 
          {
            num_sims=fsh.age_flags(20);
            fsh.proj_output_files[0]=
               new ofstream("other_projected_stuff");
            fsh.proj_output_files[1]=
               new ofstream("projected_randomized_catches");
            fsh.proj_output_files[2]=
               new ofstream("projected_randomized_catch_at_age");
            fsh.proj_output_files[3]=
               new ofstream("projected_numbers_at_age");
            fsh.proj_output_files[4]=
               new ofstream("projected_fmort_at_age_year");    //NMD 8Nov2011
            fsh.proj_output_files[5]=
               new ofstream("projected_numbers_at_age_noeff");    //NMD 10Nov2011
            fsh.proj_output_files[6]=
               new ofstream("Fmults.txt");    //YT_2018-04-26
          }
          cifstream cif("simulated_numbers_at_age");
          int tmp;
          if (!cif && fsh.age_flags(20)>0)  //NMD_Aug6_2018
          {
            cerr << "Error trying to open file simulated_numbers_at_age " << endl;
            ad_exit(1);
          }
          cif >> tmp;
          int error_flag=0;
          if (tmp !=num_sims)
          {
            error_flag=1;
            cerr << "Error -- incorrect number of simulated recruitments"
               " in simulated_numbers_at_age" << endl;
          }
          fsh.simulated_numbers_at_age.allocate
            (1,num_sims,1,fsh.num_regions,2,fsh.nage);
          cif >> fsh.simulated_numbers_at_age;
          fsh.simulated_recruitments.allocate(1,num_sims,1,fsh.num_regions,
            fsh.last_real_year+1,fsh.nyears);
  
          {
            cif >> fsh.simulated_recruitments;
	    
            if (!cif || error_flag)
            {
              fsh.simulated_recruitments.initialize();
            }
         
            fsh.simulated_recruitments=exp(fsh.simulated_recruitments);

            if (!cif)
            {
              cerr << "Error reading from file simulated_numbers_at_age" 
                   << endl;
            }
          }

//NMD 22Feb2012
	  dvector noeff_sw;
	  noeff_sw=column(fsh.fish_flags,55);
          if(sum(noeff_sw))
	  {
            cifstream cif("simulated_numbers_at_age_noeff");
            if (!cif && fsh.age_flags(20)>0)
            {
              cerr << "Error reading from file simulated_numbers_at_age_noeff" 
                 << endl;
              ad_exit(1);  //NMD_aug06_2018
            }
            int tmp;
            cif >> tmp;
            if (tmp !=num_sims)
            {
              cerr << "Error -- incorrect number of simulated recruitments"
                " in simulated_numbers_at_age_noeff" << endl;
            }
            fsh.simulated_numbers_at_age_noeff.allocate
              (1,num_sims,1,fsh.num_regions,2,fsh.nage);
            cif >> fsh.simulated_numbers_at_age_noeff;
	  }
          {
            const imatrix & _simyears=fsh.simyears;
            ADUNCONST(imatrix,simyears)
            if (allocated(simyears))
              simyears.deallocate();

            ifstream ifs("simyears");
            int ns=0; int year1=0; int yearlst;
            ifs >> ns;
            if (ns !=num_sims)
            {
              cerr << "Warning -- number of simulations in file simyears"
               " is different from the number of simulations being done"
               << endl;
              simyears.deallocate();
            }
            else
            {
              ifs >> year1 >> yearlst;
              simyears.allocate(1,ns,year1,yearlst);
              ifs >> simyears;
              if (!ifs)
              {
                cerr << "Did not successfuly read in information from file"
                        " simyears" << endl;
                simyears.deallocate();
              }
            }
          }

//NMD 22Feb2012

          for (int i=1;i<=num_sims;i++)
          {
            cout << "Doing simulation " << i << " of " << num_sims << endl;
            int prswitch=1;
            if (fsh.do_fishery_projections_flag==1 && fsh.age_flags(20)>0)   //NMD 3Mar2012 - allows iterative simulations
//			if (fsh.do_fishery_projections_flag==1&& num_sims>1) 
            {
              fsh.projection_sim_index=i;
              prswitch=1;
            }
            if(fsh.parest_flags(1)==1 && sum(fsh.data_fish_flags(2)))   //NMD 20Jan2012
            {
              fsh.bh_variance = fsh.save_bh_variance;
            }                                                           //NMD 20Jan2012
            stupid_print_switch=1;
            f=fcomp(fsh,dvar_vector(x),nvar,prswitch,gbest,0);
            stupid_print_switch=0;

            if (fsh.parest_flags(246) && fsh.generate_report)
            {
              print_indep_vars_report(gbest,fsh.indep_var,
                fsh.indep_var_lo,fsh.indep_var_hi,nvar);
	    }

//NMD 9Nov2011
            // now call fcomp with fishing mortality turned off
            // in selected fisheries
            // right now turn off in all fisheries
            fsh.q_flag=column(fsh.fish_flags,55);
            d3_array tmp_obs_tot_catch;
            dmatrix tmp_log_effort_by_fishery; //NMD 9Aug2022
            if (sum(fsh.q_flag) && fsh.age_flags(20)>0)  //NMD 18Nov2011
            {
              ivector ff55=column(fsh.fish_flags,55);
              if (sum(ff55))
              {
                // save the real obs_tot_catch data
                tmp_obs_tot_catch.allocate(1,fsh.num_regions,1,fsh.num_fish_periods,
                    1,fsh.num_fish_incidents);
                tmp_obs_tot_catch=fsh.obs_tot_catch;

                // zero out the catches from fisheries that are to be turned off;
                fsh.zero_out_catches();

                // save the real log_effort_by_fishery data
                tmp_log_effort_by_fishery.allocate(1,fsh.num_fisheries,
                  1,fsh.num_fish_times);   //NMD_9Aug2022
                tmp_log_effort_by_fishery=fsh.log_effort_by_fishery;

                // zero out the effort from fisheries that are to be turned off;
                if (fsh.age_flags(92)==2 && sum(fsh.data_fish_flags(2)))
                {
                  fsh.zero_out_log_effort();
                }
      
                // so this will have to be redone since the number of missing
                // catch situations may have changed
                //fsh.set_missing_totcatch_flags();
              }
              if (fsh.parest_flags(152)>0)
              {
                dvar_vector y=xscale(dvar_vector(x),grad_scale);
                //f=fcomp(fsh,y,nvar,1,gbest,0);
                int save1=fsh.af170q0ex;
                int save2=fsh.age_flags(170);
                int save3=fsh.af170q0;
                fsh.af170q0ex=1;
                fsh.af170q0=1;
                fsh.age_flags(170)=1;
                f=fcomp(fsh,y,nvar,1,gbest,0,&fsh.q_flag);
                fsh.af170q0ex=save1;
                fsh.age_flags(170)=save2;
                fsh.af170q0=save3;
              }
              else
              {
               // f=fcomp(fsh,dvar_vector(x),nvar,1,gbest,0);
                int save1=fsh.af170q0ex;
                int save2=fsh.age_flags(170);
                int save3=fsh.af170q0;
                fsh.af170q0ex=1;
                fsh.af170q0=1;
                fsh.age_flags(170)=1;
                fsh.allocate_optional_stuff();
                f=fcomp(fsh,dvar_vector(x),nvar,1,gbest,0,&fsh.q_flag);
                fsh.af170q0ex=save1;
                fsh.age_flags(170)=save2;
                fsh.af170q0=save3;
              }
              if (sum(ff55))
              {
                // get the real obs_tot_catch back
                fsh.obs_tot_catch=tmp_obs_tot_catch;
                // get the real log_effort_by_fishery back
                fsh.log_effort_by_fishery=tmp_log_effort_by_fishery; //NMD_9Aug2022
      
                //and recalculate the missing catch flags etc. 
                //
                //fsh.set_missing_totcatch_flags();
              }
              cout << "===== newl3, returned running fzero" <<endl;
              //Reset q_flag to zero
			  fsh.q_flag = 0;
            }
//NMD 9Nov2011
          }
          if (fsh.proj_output_files[0])
          {
#if defined(close)
#  undef close
            fsh.proj_output_files[0]->close();
#  define close _close
#else
            fsh.proj_output_files[0]->close();
#endif
            delete fsh.proj_output_files[0];
            fsh.proj_output_files[0]=0;
          }
          fsh.do_fishery_projections_flag=0;
        }
    
    
        {
          ofstream ofs("gradient.rpt");
          ofstream ofs1("sorted_gradient.rpt");
          ofs << nvar << endl;
          if (fsh.parest_flags(152)>0)
          {
            ofs << elem_prod(gbest,.01+grad_scale) << endl;
          }
          else
          {
            dmatrix M(1,2,1,gbest.indexmax());
            M(1).fill_seqadd(1,1);
            M(2)=fabs(gbest);
            dmatrix N=sort(trans(M),2);
            ofs << gbest << endl << endl;
            for (int i=1;i<=nvar;i++)
            {
              ofs1 << setw(10) << setw(12) << N(i) << endl;
              ofs << setw(10) << i << " " << setw(12) << gbest(i) << endl;
            }
          }
        }
        if (fmc.quit_flag=='Q')
        {
          quit_flag=1;
          //cout << "===== setting quit_flag = ";
          //cout << quit_flag << endl;
        }
        else if (fmc.quit_flag=='N')
        {
          quit_flag=2;
          //cout << "===== setting quit_flag = ";
          //cout <<"newl3.cpp " << quit_flag << endl;
        }
        else
        {
          quit_flag=0;
          //cout << "===== setting quit_flag = ";
          //cout <<"newl3.cpp " << quit_flag << endl;
        }
        if (fsh.age_flags(52))
        {
          quit_flag=1;
        }
        cout << "===== newl3.cpp,  setting quit_flag = "<< quit_flag << endl;

        // now call fcomp with fishing mortality turned off
        // in selected fisheries
        // right now turn off in all fisheries
        fsh.q_flag=column(fsh.fish_flags,55);
        d3_array tmp_obs_tot_catch;
        dmatrix tmp_log_effort_by_fishery; //NMD 9Aug2022
//        if (sum(fsh.q_flag))
        if (sum(fsh.q_flag) && fsh.age_flags(20)<1)  //NMD 10Nov2011
        {
          ivector ff55=column(fsh.fish_flags,55);
          if (sum(ff55))
          {
            // save the real obs_tot_catch data
            tmp_obs_tot_catch.allocate(1,fsh.num_regions,1,fsh.num_fish_periods,
               1,fsh.num_fish_incidents);
            tmp_obs_tot_catch=fsh.obs_tot_catch;
      
            // zero out the catches from fisheries that are to be turned off;
            fsh.zero_out_catches();

            // save the real log_effort_by_fishery data
            tmp_log_effort_by_fishery.allocate(1,fsh.num_fisheries,
              1,fsh.num_fish_times);   //NMD_9Aug2022
            tmp_log_effort_by_fishery=fsh.log_effort_by_fishery;

            // zero out the effort from fisheries that are to be turned off;
            if (fsh.age_flags(92)==2 && sum(fsh.data_fish_flags(2)))
            {
              fsh.zero_out_log_effort();
            }      
            // so this will have to be redone since the number of missing
            // catch situations may have changed
            //fsh.set_missing_totcatch_flags();
          }
          if (fsh.parest_flags(152)>0)
          {
            dvar_vector y=xscale(dvar_vector(x),grad_scale);
            //f=fcomp(fsh,y,nvar,1,gbest,0);
            int save1=fsh.af170q0ex;
            int save2=fsh.age_flags(170);
            int save3=fsh.af170q0;
            fsh.af170q0ex=1;
            fsh.af170q0=1;
            fsh.age_flags(170)=1;
            if (sum(fsh.data_fish_flags(2)))
            {
              fsh.do_fishery_projections_flag=1;
            }
            f=fcomp(fsh,y,nvar,1,gbest,0,&fsh.q_flag);
            fsh.do_fishery_projections_flag=0;
            fsh.af170q0ex=save1;
            fsh.age_flags(170)=save2;
            fsh.af170q0=save3;
          }
          else
          {
           // f=fcomp(fsh,dvar_vector(x),nvar,1,gbest,0);
            int save1=fsh.af170q0ex;
            int save2=fsh.age_flags(170);
            int save3=fsh.af170q0;
            fsh.af170q0ex=1;
            fsh.af170q0=1;
            fsh.age_flags(170)=1;
            fsh.allocate_optional_stuff();

            if (sum(fsh.data_fish_flags(2)))
            {
              fsh.do_fishery_projections_flag=1;
            }
            f=fcomp(fsh,dvar_vector(x),nvar,1,gbest,0,&fsh.q_flag);
            fsh.do_fishery_projections_flag=0;
            fsh.af170q0ex=save1;
            fsh.age_flags(170)=save2;
            fsh.af170q0=save3;
          }
          if (sum(ff55))
          {
            // get the real obs_tot_catch back
            fsh.obs_tot_catch=tmp_obs_tot_catch;
            // get the real log_effort_by_fishery back
            fsh.log_effort_by_fishery=tmp_log_effort_by_fishery; //NMD_9Aug2022
      
            //and recalculate the missing catch flags etc. 
            //
            //fsh.set_missing_totcatch_flags();
          }
        }
        cout << "===== newl3, returned printing fcomp" <<endl;
      }
      if (mf_pvm->pvm_save_switch)
      {
        mf_pvm->pvm_switch=mf_pvm->pvm_save_switch;
        mf_pvm->pvm_save_switch=0;
      }
    }
  }
  else if (trunc_newton_flag==2)
  {
    // trunc newton
    int nterms=7;
    if (fsh.parest_flags(192)) nterms=fsh.parest_flags(192);
    //fmmt fmc(nvar,nterms);
    //fmmt1 fmc(nvar,nterms,angle_bound);
#if (defined(USE_DD))
    dfmmt1 fmc(nvar,nterms,angle_bound);
#else
    fmmt1 fmc(nvar,nterms,angle_bound);
#endif
    if (fsh.parest_flags(2)==0)
      fmc.iprint=5;
    else
      fmc.iprint=fsh.parest_flags(2);

    int old_itn=0;
    //  Rationalize meaning of convergence flag -- P.K. 7/17/2004 
    //  Default is now 1.0 instead of 1.e-6
    fmc.crit=pow(10.,double(fsh.parest_flags(50)));

    int nvar_check; 
    dvector grad_scale(1,nvar);
    if (fsh.parest_flags(353)==0)
    {
      gradient_structure::set_YES_DERIVATIVES();
    }
    else
    {
      gradient_structure::set_NO_DERIVATIVES();
    }
    if (fsh.parest_flags(152)>0)
    {
      ifstream ifs("gradient.rpt");
      ifs >> nvar_check;
      if (nvar_check !=nvar)
      {
        cerr << "Number of parameters in gradient scaling file is"
                " not equal to the number" << endl << 
                " of active parameters -- scaling is not enabled" << endl;
        fsh.parest_flags(152)=0;
      }
      else
      {
        ifs >> grad_scale;
        grad_scale= fabs(grad_scale);
      }
    }
    if (fsh.parest_flags(152)>0)
    {
      inverse_xscale(x,grad_scale);
    }
    fmc.scroll_flag=0;
    fmc.imax=40;
    fmc.min_improve=0.;
    fmc.maxfn=fsh.parest_flags(1);
    //cout << " fmc.maxfn = "  << fmc.maxfn  << endl;
    fmc.ifn=current_ifn;
    current_ifn=0;
    int it=fsh.parest_flags(12)+fsh.parest_flags(13)+fsh.parest_flags(14)+
      fsh.parest_flags(15)+fsh.parest_flags(16);
    fsh.parest_flags(143)=0;
    fsh.parest_flags(144)=0;
    MY_DOUBLE_TYPE f;
    MY_DOUBLE_TYPE fbest=1.e+50;;
    dvector g(1,nvar);
    g.initialize();
    dvector gbest(1,nvar);
    greport("very beginning ");
    int do_fmin_flag=0;

    gbest.fill_seqadd(1.e+50,0.);
    {
      if (fsh.age_flags(52)==0)
      {
#if defined(AD_FPEXCEPTION)
        //int oldexception=_control87(0,MCW_EM);
        //_control87(0, 0);
        signal(SIGFPE, &adfloat_except);
#endif
        if (fsh.parest_flags(222)>0)
        {
          fmc.dcheck_flag=fsh.parest_flags(222)-1;
        }
        adtimer xdt;
        while (fmc.ireturn>=0)
        {
          //global_except_switch++;
          int badflag=0;
          xdt.get_elapsed_time_and_reset();
#if defined(AD_FPEXCEPTION)
          try 
#endif
          {
            if (do_fmin_flag==0)
            {
              fmc.fmin(f,x,g);
            }
            int itmp=fmc.itn;
            int itmp1=10;
            if (fsh.parest_flags(196))
            {
              itmp1=fsh.parest_flags(196);
            }
               
            if (fmc.itn>3 && fmc.itn%itmp1==0)
            {
              save_best_estimate_only(itmp,xbest,nvar,gbest,fsh);
            }
          }
#if defined(AD_FPEXCEPTION)
          catch (choleski_exception & myexp)
          {
            save_best_estimate(test_adexception_flag,x,xbest,nvar,
              gbest,fsh,fmc.ifn);
          }
#endif
          cout << "FMIN time = " << xdt.get_elapsed_time_and_reset()/1000. << endl;
          {
            int x;
            //cout << "Enter x" << endl;
            //cin >> x;
          }
#if defined(USE_ADPVM)
          if (mf_pvm->pvm_switch==1)
          {
            mfsend_int_to_slaves(fmc.ireturn);
          }
#endif //#if defined(USE_ADPVM)
          if (fmc.ihang)
          {
            hang_flag=fmc.ihang;
            maxg=max(g);
            ccrit=fmc.crit;
            current_ifn=fmc.ifn;
          }
          if (fmc.ireturn>0)
          {
#if defined(AD_FPEXCEPTION)
            try 
#endif
            {
              /*
              int globalstop;
              ifstream ifsg("globalstop");
              ifsg >> globalstop;
              ifsg.close();
              if (fmc.itn==globalstop)
              {
                ofstream ofsg("globalstop");
                globalstop+=2;
                ofsg << globalstop;
                adfloat_except(1);
              }
              */
              if (fsh.parest_flags(152)>0)
              {
                dvar_vector y=xscale(dvar_vector(x),grad_scale);
                f=fcomp(fsh,y,nvar,0,gbest,0);
                fsh.parest_flags(144)=0;
              }
              else
              {
#if defined(USE_ADPVM)
                if (mf_pvm->pvm_switch==1)
                {
                  dvar_vector vx=dvar_vector(x);
                  mfsend_x_to_slaves(vx);
                  f=fcomp(fsh,vx,nvar,0,gbest,0);
                }
                else
#endif //#if defined(USE_ADPVM)
                {
                  //fsh.do_fishery_projections_flag=1;
                  MY_DOUBLE_TYPE zz=-1;
                  signal(SIGFPE, &adfloat_except);
                  //cout << sqrt(zz) << endl;
                  
                  f=fcomp(fsh,dvar_vector(x),nvar,0,gbest,0);
//                  cout << " LLL f = " << f << endl;
       
                  if(fmc.ireturn ==1 && fsh.parest_flags(1)==1 && sum(fsh.data_fish_flags(2)))   //NMD 20Jan2012
                  {
                    fsh.save_bh_variance = fsh.bh_variance;
                  }                                                                              //NMD 20Jan2012
                  fsh.do_fishery_projections_flag=0;
                }
                fsh.parest_flags(144)=0;
              }
              badflag=1;
              cout << "==============================================" << endl;
	      cout << " f eval " << fmc.ifn << endl;
              char * cptr = (char *) &f;
              for (int iw=0;iw<8;iw++)
              {
                if ( int(*(cptr+iw)) != iw ) badflag=0;
              }
            }
#if defined(AD_FPEXCEPTION)
          
            catch (choleski_exception & myexp)
            {
              save_best_estimate(test_adexception_flag,x,xbest,nvar,
                gbest,fsh,fmc.ifn);
            }
            
            catch (myexception & adfe)
            {
              //save current best parameter estimates
              test_adexception_flag=0;
              x=xbest;
              f=fcomp(fsh,dvar_vector(x),nvar,0,gbest,0);
              par_ofstream ofs("crash_save");
              ofs << fsh << endl;
              cerr << " exiting after choleski error and parameter save" << endl;
              ad_exit(1);
            }
            
            catch (ad_float_exception & adfe)
            {
              badflag=1;
              reset_gradient_stack();

              //save current best parameter estimates
              x=xbest;
              test_adexception_flag=0;
              f=fcomp(fsh,dvar_vector(x),nvar,0,gbest,0);
              par_ofstream ofs("crash_save");
              ofs << fsh << endl;
              
              cerr << " exiting after FPE and parameter save" << endl;
              ad_exit(1);

              /*
              //_fpreset();
              signal(SIGFPE, &adfloat_except);
              //_control87(0, 0);
              delete pgs;
              pgs = new gradient_structure(gss);
              //_control87(0,MCW_EM);
              */
            }
#endif
            if (badflag)
            {
              f=1.e+95;
            }
#if defined(AD_FPEXCEPTION)
            int in_gradcalc_flag=1;
            try 
#endif
            {
              cout << "calling gradcalc" << endl;
              if (fsh.parest_flags(353)==0)
              {
                gradcalc(nvar,g);
              }
              in_gradcalc_flag=0;
            }
#if defined(AD_FPEXCEPTION)
            catch (choleski_exception & adfe)
            {
              if (in_gradcalc_flag && gradient_structure::get_save_var_flag())
              {
                gradient_structure::restore_arrays();
                gradient_structure::restore_variables();
              }
              save_best_estimate(test_adexception_flag,x,xbest,nvar,
                gbest,fsh,fmc.ifn);
            }
            catch (myexception & adfe)
            {
              //save current best parameter estimates
              test_adexception_flag=0;
              x=xbest;
              f=fcomp(fsh,dvar_vector(x),nvar,0,gbest,0);
              par_ofstream ofs("crash_save");
              ofs << fsh << endl;
              cerr << " exiting after choleski error and parameter save" << endl;
              ad_exit(1);
            }
#endif
            

            if (0)
            {
              random_number_generator rng(435);
              dvector eps(1,nvar);
              eps.fill_randn(rng);
              MY_DOUBLE_TYPE gn=norm(g);
              dvector sg=g/gn;
              
              ofstream ofs("feasible");
              // test for feasible direction
              gradient_structure::set_NO_DERIVATIVES();
              MY_DOUBLE_TYPE fc=f;
              ofs << "current f = " << setscientific() << setprecision(16) << f << " " << fc-f << endl;
              dvector y=x+1.e-6*g;
              MY_DOUBLE_TYPE f1=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
              ofs << "1.e-6 f = " << setscientific() << setprecision(16) << f << " " << fc-f1 << endl;
              /*
              y=x+1.e-9*g;
              MY_DOUBLE_TYPE f2=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
              ofs << "1.e-9 f = " << setscientific() << setprecision(16) << f << " " << fc-f2 << endl;
              y=x+1.e-10*g;
              MY_DOUBLE_TYPE f3=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
              ofs << "1.e-10 f = " << setscientific() << setprecision(16) << f << " " << fc-f3 << endl;
              y=x+1.e-12*g;
              MY_DOUBLE_TYPE f4=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
              ofs << "1.e-12 f = " << setscientific() << setprecision(16) << f << " " << fc-f4 << endl;
              y=x+1.e-13*g;
              MY_DOUBLE_TYPE f5=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
              ofs << "1.e-13 f = " << setscientific() << setprecision(16) << f << " " << fc-f5 << endl;
              */

              for (int i=1;i<=20;i++)
              {
                eps.fill_randn(rng);
                dvector newg=(sg+.05*eps/norm(eps));
                newg/=norm(newg);
                newg*=gn;
  
                y=x+1.e-9*newg;
                MY_DOUBLE_TYPE f8=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
                ofs << "1.e-9 f = " << setscientific() << setprecision(16) << f << " " << fc-f8 << endl;
                y=x+1.e-10*newg;
                MY_DOUBLE_TYPE f9=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
                ofs << "1.e-10 f = " << setscientific() << setprecision(16) << f << " " << fc-f9 << endl;
                y=x+1.e-11*newg;
                MY_DOUBLE_TYPE f10=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
                ofs << "1.e-11 f = " << setscientific() << setprecision(16) << f << " " << fc-f10 << endl;
                y=x+1.e-12*newg;
                MY_DOUBLE_TYPE f11=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
                ofs << "1.e-12 f = " << setscientific() << setprecision(16) << f << " " << fc-f11 << endl;
              }
              ad_exit(1);
            }
       
            cout << "finished gradcalc" << endl;
            if (f < fbest)
            {
              fbest=f;
              gbest=g;
              xbest=x;
            }
          }
          if (fsh.parest_flags(333) 
            && 
             (!((fmc.itn)%(fsh.parest_flags(333)))
              || fmc.itn==1) && fmc.itn>old_itn)
          {
            old_itn=fmc.itn;
            //save current best parameter estimates
            int test_adexception_flag_save=test_adexception_flag;
            test_adexception_flag=0;
            gradient_structure::set_NO_DERIVATIVES();
            dvector gtemp(1,nvar);
            MY_DOUBLE_TYPE ftemp=fcomp(fsh,dvar_vector(xbest),nvar,0,gtemp,0);
            int suffix=0;
            if (fmc.itn==1)
              suffix=1;
            else
              suffix=fmc.itn;
            adstring file_name= adstring("intermediate_results_") 
              + str((const int)(suffix));
            par_ofstream ofs(file_name);
            ofs << fsh << endl;
            ofs << "# Objective function value" << endl << "   " << -ftemp 
                << endl;
            ofs << "# The number of parameters" << endl << "   " << nvar << endl;
            adstring command="cp " + file_name + " intermediate_results";
            if (fsh.num_tag_releases)
              ofs << "# Likelihood component for tags " <<fsh.likecomponent[0]<< endl;
            test_adexception_flag=test_adexception_flag_save;
            system(command);
            adstring command1="echo " + str((const int)fmc.itn) +  " > success_itn";
            system(command1);
       
            adstring command2="echo " + str((const int)fmc.ifn) +  " > success_ifn";
            system(command1);
            gradient_structure::set_YES_DERIVATIVES();
          }
        }
        x=xbest;

#if defined(AD_FPEXCEPTION)
        signal(SIGFPE, SIG_DFL);
#endif
      }
     // *********************************************************
     // *********************************************************
     // !!! Dave F Dec 19 02 If you want pvm used in other calculations 
     //  remove this and modify program
        mf_pvm->pvm_save_switch=mf_pvm->pvm_switch;
        mf_pvm->pvm_switch=0;
     // *********************************************************
     // *********************************************************
     // if (!hang_flag)
      {
        //cout << endl << "Calling printing fcomp" <<endl;
        gradient_structure::set_NO_DERIVATIVES();
    
        if (fsh.parest_flags(152)>0)
        {
          dvar_vector y=xscale(dvar_vector(x),grad_scale);
          f=fcomp(fsh,y,nvar,1,gbest,0);
        }
        else
        {
          if (sum(fsh.data_fish_flags(2)))
          {
            fsh.do_fishery_projections_flag=1;
          }
          int num_sims=1;
          if (fsh.do_fishery_projections_flag==1 && fsh.age_flags(20)>0) 
          {
            num_sims=fsh.age_flags(20);
            fsh.proj_output_files[0]=
               new ofstream("other_projected_stuff");
            fsh.proj_output_files[1]=
               new ofstream("projected_randomized_catches");
            fsh.proj_output_files[2]=
               new ofstream("projected_randomized_catch_at_age");
            fsh.proj_output_files[3]=
               new ofstream("projected_numbers_at_age");
            fsh.proj_output_files[4]=
               new ofstream("projected_fmort_at_age_year");    //NMD 8Nov2011
            fsh.proj_output_files[5]=
               new ofstream("projected_numbers_at_age_noeff");    //NMD 10Nov2011
            fsh.proj_output_files[6]=
               new ofstream("Fmults.txt");    //YT_2018-04-26
          }
          cifstream cif("simulated_numbers_at_age");
          if (!cif && fsh.age_flags(20)>0)  //NMD_Aug6_2018
          {
            cerr << "Error trying to open file simulated_numbers_at_age " << endl;
            ad_exit(1);
          }
          int tmp;
          cif >> tmp;
          int error_flag=0;
          if (tmp !=num_sims)
          {
            error_flag=1;
            cerr << "Error -- incorrect number of simulated recruitments"
               " in simulated_numbers_at_age" << endl;
          }
          fsh.simulated_numbers_at_age.allocate
            (1,num_sims,1,fsh.num_regions,2,fsh.nage);
          cif >> fsh.simulated_numbers_at_age;
          fsh.simulated_recruitments.allocate(1,num_sims,1,fsh.num_regions,
            fsh.last_real_year+1,fsh.nyears);
          {
      
          cif >> fsh.simulated_recruitments;
          if (!cif || error_flag)
          {
            fsh.simulated_recruitments.initialize();
          }
         
          fsh.simulated_recruitments=exp(fsh.simulated_recruitments);

          if (!cif)
          {
            cerr << "Error reading from file simulated_numbers_at_age" 
                 << endl;
          }
          }

//NMD 22Feb2012
	  dvector noeff_sw;
	  noeff_sw=column(fsh.fish_flags,55);
          if(sum(noeff_sw))
	  {
            cifstream cif("simulated_numbers_at_age_noeff");
            if (!cif && fsh.age_flags(20)>0)
            {
              cerr << "Error reading from file simulated_numbers_at_age_noeff" 
                 << endl;
              ad_exit(1);  //NMD_aug06_2018
            }
            int tmp;
            cif >> tmp;
            if (tmp !=num_sims)
            {
              cerr << "Error -- incorrect number of simulated recruitments"
                " in simulated_numbers_at_age_noeff" << endl;
            }
            fsh.simulated_numbers_at_age_noeff.allocate
              (1,num_sims,1,fsh.num_regions,2,fsh.nage);
            cif >> fsh.simulated_numbers_at_age_noeff;
	  }
          {
            const imatrix & _simyears=fsh.simyears;
            ADUNCONST(imatrix,simyears)
            if (allocated(simyears))
              simyears.deallocate();

            ifstream ifs("simyears");
            int ns=0; int year1=0; int yearlst;
            ifs >> ns;
            if (ns !=num_sims)
            {
              cerr << "Warning -- number of simulations in file simyears"
               " is different from the number of simulations being done"
               << endl;
              simyears.deallocate();
            }
            else
            {
              ifs >> year1 >> yearlst;
              simyears.allocate(1,ns,year1,yearlst);
              ifs >> simyears;
              if (!ifs)
              {
                cerr << "Did not successfuly read in information from file"
                        " simyears" << endl;
                simyears.deallocate();
              }
            }
          }

//NMD 22Feb2012

          for (int i=1;i<=num_sims;i++)
          {
            cout << "Doing simulation " << i << " of " << num_sims << endl;
            int prswitch=1;
            if (fsh.do_fishery_projections_flag==1 && fsh.age_flags(20)>0)   //NMD 3Mar2012 - allows iterative simulations
//			if (fsh.do_fishery_projections_flag==1&& num_sims>1) 
            {
              fsh.projection_sim_index=i;
              prswitch=1;
            }
            if(fsh.parest_flags(1)==1 && sum(fsh.data_fish_flags(2)))   //NMD 20Jan2012
            {
              fsh.bh_variance = fsh.save_bh_variance;
            }                                                           //NMD 20Jan2012
            stupid_print_switch=1;
            f=fcomp(fsh,dvar_vector(x),nvar,prswitch,gbest,0);
            stupid_print_switch=0;

            if (fsh.parest_flags(246) && fsh.generate_report)
            {
              print_indep_vars_report(gbest,fsh.indep_var,
                fsh.indep_var_lo,fsh.indep_var_hi,nvar);
	    }

//NMD 9Nov2011
            // now call fcomp with fishing mortality turned off
            // in selected fisheries
            // right now turn off in all fisheries
            fsh.q_flag=column(fsh.fish_flags,55);
            d3_array tmp_obs_tot_catch;
            dmatrix tmp_log_effort_by_fishery; //NMD 9Aug2022
            if (sum(fsh.q_flag) && fsh.age_flags(20)>0)  //NMD 18Nov2011
            {
              ivector ff55=column(fsh.fish_flags,55);
              if (sum(ff55))
              {
                // save the real obs_tot_catch data
                tmp_obs_tot_catch.allocate(1,fsh.num_regions,1,fsh.num_fish_periods,
                    1,fsh.num_fish_incidents);
                tmp_obs_tot_catch=fsh.obs_tot_catch;
      
                // zero out the catches from fisheries that are to be turned off;
                fsh.zero_out_catches();

                // save the real log_effort_by_fishery data
                tmp_log_effort_by_fishery.allocate(1,fsh.num_fisheries,
                  1,fsh.num_fish_times);   //NMD_9Aug2022
                tmp_log_effort_by_fishery=fsh.log_effort_by_fishery;

                // zero out the effort from fisheries that are to be turned off;
                if (fsh.age_flags(92)==2 && sum(fsh.data_fish_flags(2)))
                {
                  fsh.zero_out_log_effort();
                }
                // so this will have to be redone since the number of missing
                // catch situations may have changed
                //fsh.set_missing_totcatch_flags();
              }
              if (fsh.parest_flags(152)>0)
              {
                dvar_vector y=xscale(dvar_vector(x),grad_scale);
                //f=fcomp(fsh,y,nvar,1,gbest,0);
                int save1=fsh.af170q0ex;
                int save2=fsh.age_flags(170);
                int save3=fsh.af170q0;
                fsh.af170q0ex=1;
                fsh.af170q0=1;
                fsh.age_flags(170)=1;
                f=fcomp(fsh,y,nvar,1,gbest,0,&fsh.q_flag);
                fsh.af170q0ex=save1;
                fsh.age_flags(170)=save2;
                fsh.af170q0=save3;
              }
              else
              {
               // f=fcomp(fsh,dvar_vector(x),nvar,1,gbest,0);
                int save1=fsh.af170q0ex;
                int save2=fsh.age_flags(170);
                int save3=fsh.af170q0;
                fsh.af170q0ex=1;
                fsh.af170q0=1;
                fsh.age_flags(170)=1;
                fsh.allocate_optional_stuff();
                f=fcomp(fsh,dvar_vector(x),nvar,1,gbest,0,&fsh.q_flag);
                fsh.af170q0ex=save1;
                fsh.age_flags(170)=save2;
                fsh.af170q0=save3;
              }
              if (sum(ff55))
              {
                // get the real obs_tot_catch back
                fsh.obs_tot_catch=tmp_obs_tot_catch;
      
                // get the real log_effort_by_fishery back
                fsh.log_effort_by_fishery=tmp_log_effort_by_fishery; //NMD_9Aug2022

                //and recalculate the missing catch flags etc. 
                //
                //fsh.set_missing_totcatch_flags();
              }
              cout << "===== newl3, returned running fzero" <<endl;
              //Reset q_flag to zero
			  fsh.q_flag = 0;
            }
//NMD 9Nov2011
          }
          if (fsh.proj_output_files[0])
          {
#if defined(close)
#  undef close
            fsh.proj_output_files[0]->close();
#  define close _close
#else
            fsh.proj_output_files[0]->close();
#endif
            delete fsh.proj_output_files[0];
            fsh.proj_output_files[0]=0;
          }
          fsh.do_fishery_projections_flag=0;
        }
    
    
        {
          ofstream ofs("gradient.rpt");
          ofstream ofs1("sorted_gradient.rpt");
          ofs << nvar << endl;
          if (fsh.parest_flags(152)>0)
          {
            ofs << elem_prod(gbest,.01+grad_scale) << endl;
          }
          else
          {
            dmatrix M(1,2,1,gbest.indexmax());
            M(1).fill_seqadd(1,1);
            M(2)=fabs(gbest);
            dmatrix N=sort(trans(M),2);
            ofs << gbest << endl << endl;
            for (int i=1;i<=nvar;i++)
            {
              ofs1 << setw(10) << setw(12) << N(i) << endl;
              ofs << setw(10) << i << " " << setw(12) << gbest(i) << endl;
            }
          }
        }
        if (fmc.quit_flag=='Q')
        {
          quit_flag=1;
          //cout << "===== setting quit_flag = ";
          //cout << quit_flag << endl;
        }
        else if (fmc.quit_flag=='N')
        {
          quit_flag=2;
          //cout << "===== setting quit_flag = ";
          //cout <<"newl3.cpp " << quit_flag << endl;
        }
        else
        {
          quit_flag=0;
          //cout << "===== setting quit_flag = ";
          //cout <<"newl3.cpp " << quit_flag << endl;
        }
        if (fsh.age_flags(52))
        {
          quit_flag=1;
        }
        cout << "===== newl3.cpp,  setting quit_flag = "<< quit_flag << endl;

        // now call fcomp with fishing mortality turned off
        // in selected fisheries
        // right now turn off in all fisheries
        fsh.q_flag=column(fsh.fish_flags,55);
        d3_array tmp_obs_tot_catch;
        dmatrix tmp_log_effort_by_fishery; //NMD 9Aug2022
//        if (sum(fsh.q_flag))
        if (sum(fsh.q_flag) && fsh.age_flags(20)<1)  //NMD 10Nov2011
        {
          ivector ff55=column(fsh.fish_flags,55);
          if (sum(ff55))
          {
            // save the real obs_tot_catch data
            tmp_obs_tot_catch.allocate(1,fsh.num_regions,1,fsh.num_fish_periods,
               1,fsh.num_fish_incidents);
            tmp_obs_tot_catch=fsh.obs_tot_catch;
      
            // zero out the catches from fisheries that are to be turned off;
            fsh.zero_out_catches();

            // save the real log_effort_by_fishery data
            tmp_log_effort_by_fishery.allocate(1,fsh.num_fisheries,
              1,fsh.num_fish_times);   //NMD_9Aug2022
            tmp_log_effort_by_fishery=fsh.log_effort_by_fishery;

            // zero out the effort from fisheries that are to be turned off;
            if (fsh.age_flags(92)==2 && sum(fsh.data_fish_flags(2)))
            {
              fsh.zero_out_log_effort();
            }      
            // so this will have to be redone since the number of missing
            // catch situations may have changed
            //fsh.set_missing_totcatch_flags();
          }
          if (fsh.parest_flags(152)>0)
          {
            dvar_vector y=xscale(dvar_vector(x),grad_scale);
            //f=fcomp(fsh,y,nvar,1,gbest,0);
            int save1=fsh.af170q0ex;
            int save2=fsh.age_flags(170);
            int save3=fsh.af170q0;
            fsh.af170q0ex=1;
            fsh.af170q0=1;
            fsh.age_flags(170)=1;
            if (sum(fsh.data_fish_flags(2)))
            {
              fsh.do_fishery_projections_flag=1;
            }
            f=fcomp(fsh,y,nvar,1,gbest,0,&fsh.q_flag);
            fsh.do_fishery_projections_flag=0;
            fsh.af170q0ex=save1;
            fsh.age_flags(170)=save2;
            fsh.af170q0=save3;
          }
          else
          {
           // f=fcomp(fsh,dvar_vector(x),nvar,1,gbest,0);
            int save1=fsh.af170q0ex;
            int save2=fsh.age_flags(170);
            int save3=fsh.af170q0;
            fsh.af170q0ex=1;
            fsh.af170q0=1;
            fsh.age_flags(170)=1;
            fsh.allocate_optional_stuff();

            if (sum(fsh.data_fish_flags(2)))
            {
              fsh.do_fishery_projections_flag=1;
            }
            f=fcomp(fsh,dvar_vector(x),nvar,1,gbest,0,&fsh.q_flag);
            fsh.do_fishery_projections_flag=0;
            fsh.af170q0ex=save1;
            fsh.age_flags(170)=save2;
            fsh.af170q0=save3;
          }
          if (sum(ff55))
          {
            // get the real obs_tot_catch back
            fsh.obs_tot_catch=tmp_obs_tot_catch;

            // get the real log_effort_by_fishery back
            fsh.log_effort_by_fishery=tmp_log_effort_by_fishery; //NMD_9Aug2022

            //and recalculate the missing catch flags etc. 
            //
            //fsh.set_missing_totcatch_flags();
          }
        }
        cout << "===== newl3, returned printing fcomp" <<endl;
      }
      if (mf_pvm->pvm_save_switch)
      {
        mf_pvm->pvm_switch=mf_pvm->pvm_save_switch;
        mf_pvm->pvm_save_switch=0;
      }
    }
  }
  else
  {
    // quasi newton
    fmm fmc(nvar,fsh.parest_flags(198));
    if (ddcheck_flag>-1)  
    {
      fmc.dcheck_flag=ddcheck_flag;
    }
    if (fsh.parest_flags(2)==0)
      fmc.iprint=5;
    else
      fmc.iprint=fsh.parest_flags(2);

    int old_itn=0;
    //  Rationalize meaning of convergence flag -- P.K. 7/17/2004 
    //  Default is now 1.0 instead of 1.e-6
    //fmc.crit=1.e-6;
    //if (fsh.parest_flags(50)!=0) 
    //{
    //  fmc.crit=pow(10.,-double(fsh.parest_flags(50)));
    //}  
    fmc.crit=pow(10.,double(fsh.parest_flags(50)));

    int nvar_check; 
    dvector grad_scale(1,nvar);
    if (fsh.parest_flags(353)==0)
    {
      gradient_structure::set_YES_DERIVATIVES();
    }
    else
    {
      gradient_structure::set_NO_DERIVATIVES();
    }

    if (fsh.parest_flags(152)>0)
    {
      ifstream ifs("gradient.rpt");
      ifs >> nvar_check;
      if (nvar_check !=nvar)
      {
        cerr << "Number of parameters in gradient scaling file is"
                " not equal to the number" << endl << 
                " of active parameters -- scaling is not enabled" << endl;
        fsh.parest_flags(152)=0;
      }
      else
      {
        ifs >> grad_scale;
        grad_scale= fabs(grad_scale);
      }
    }
    if (fsh.parest_flags(152)>0)
    {
      inverse_xscale(x,grad_scale);
    }
    fmc.scroll_flag=0;
    fmc.imax=40;
    fmc.min_improve=0.;
    fmc.maxfn=fsh.parest_flags(1);
    //cout << " fmc.maxfn = "  << fmc.maxfn  << endl;
    fmc.ifn=current_ifn;
    current_ifn=0;
    int it=fsh.parest_flags(12)+fsh.parest_flags(13)+fsh.parest_flags(14)+
      fsh.parest_flags(15)+fsh.parest_flags(16);
    fsh.parest_flags(143)=0;
    fsh.parest_flags(144)=0;
    MY_DOUBLE_TYPE f;
    MY_DOUBLE_TYPE fbest=1.e+50;;
    dvector g(1,nvar);
    g.initialize();
    dvector gbest(1,nvar);
    greport("very beginning ");
    int do_fmin_flag=0;

    gbest.fill_seqadd(1.e+50,0.);
    {
      if (fsh.age_flags(52)==0)
      {
#if defined(AD_FPEXCEPTION)
        //int oldexception=_control87(0,MCW_EM);
        //_control87(0, 0);
        signal(SIGFPE, &adfloat_except);
#endif
        if (fsh.parest_flags(222)>0)
        {
          fmc.dcheck_flag=fsh.parest_flags(222)-1;
        }
        adtimer xdt;
        while (fmc.ireturn>=0)
        {
          //global_except_switch++;
          int badflag=0;
          xdt.get_elapsed_time_and_reset();
#if defined(AD_FPEXCEPTION)
          try 
#endif
          {
            if (do_fmin_flag==0)
            {
              fmc.fmin(f,x,g);
            }
            int itmp=fmc.itn;
            int itmp1=10;
            if (fsh.parest_flags(196))
            {
              itmp1=fsh.parest_flags(196);
            }
               
            if (fmc.itn>3 && fmc.itn%itmp1==0)
            {
              save_best_estimate_only(itmp,xbest,nvar,gbest,fsh);
            }
          }
#if defined(AD_FPEXCEPTION)
          catch (choleski_exception & myexp)
          {
            save_best_estimate(test_adexception_flag,x,xbest,nvar,
              gbest,fsh,fmc.ifn);
          }
#endif
          cout << "FMIN time = " << xdt.get_elapsed_time_and_reset()/1000. << endl;
          {
            int x;
            //cout << "Enter x" << endl;
            //cin >> x;
          }
#if defined(USE_ADPVM)
          if (mf_pvm->pvm_switch==1)
          {
            mfsend_int_to_slaves(fmc.ireturn);
          }
#endif //#if defined(USE_ADPVM)
          if (fmc.ihang)
          {
            hang_flag=fmc.ihang;
            maxg=max(g);
            ccrit=fmc.crit;
            current_ifn=fmc.ifn;
          }
          if (fmc.ireturn>0)
          {
#if defined(AD_FPEXCEPTION)
            try 
#endif
            {
              /*
              int globalstop;
              ifstream ifsg("globalstop");
              ifsg >> globalstop;
              ifsg.close();
              if (fmc.itn==globalstop)
              {
                ofstream ofsg("globalstop");
                globalstop+=2;
                ofsg << globalstop;
                adfloat_except(1);
              }
              */
              if (fsh.parest_flags(152)>0)
              {
                dvar_vector y=xscale(dvar_vector(x),grad_scale);
                f=fcomp(fsh,y,nvar,0,gbest,0);
                fsh.parest_flags(144)=0;
              }
              else
              {
#if defined(USE_ADPVM)
                if (mf_pvm->pvm_switch==1)
                {
                  dvar_vector vx=dvar_vector(x);
                  mfsend_x_to_slaves(vx);
                  f=fcomp(fsh,vx,nvar,0,gbest,0);
                }
                else
#endif //#if defined(USE_ADPVM)
                {
                  //fsh.do_fishery_projections_flag=1;
                  MY_DOUBLE_TYPE zz=-1;
                  signal(SIGFPE, &adfloat_except);
                  //cout << sqrt(zz) << endl;
                  if (fsh.parest_flags(390)>0  && fsh.parest_flags(391)==0)
                  {
                    ofstream ofs("x_vector." 
                      + str(fsh.parest_flags(390)),ios::app);
                    ofs << fmc.itn << " " << fmc.ifn << " " << x.indexmax() 
                        << endl << setscientific() << setprecision(15) << x
                        << endl;
                   }
                  
                  f=fcomp(fsh,dvar_vector(x),nvar,0,gbest,0);
//                  cout << " LLL f = " << f << endl;
       
                  if(fmc.ireturn ==1 && fsh.parest_flags(1)==1 && sum(fsh.data_fish_flags(2)))   //NMD 20Jan2012
                  {
                    fsh.save_bh_variance = fsh.bh_variance;
                  }                                                                              //NMD 20Jan2012
                  fsh.do_fishery_projections_flag=0;
                }
                fsh.parest_flags(144)=0;
              }
              badflag=1;
              cout << "==============================================" << endl;
	      cout << " f eval " << fmc.ifn << endl;
              char * cptr = (char *) &f;
              for (int iw=0;iw<8;iw++)
              {
                if ( int(*(cptr+iw)) != iw ) badflag=0;
              }
            }
#if defined(AD_FPEXCEPTION)
          
            catch (choleski_exception & myexp)
            {
              save_best_estimate(test_adexception_flag,x,xbest,nvar,
                gbest,fsh,fmc.ifn);
            }
            
            catch (myexception & adfe)
            {
              //save current best parameter estimates
              test_adexception_flag=0;
              x=xbest;
              f=fcomp(fsh,dvar_vector(x),nvar,0,gbest,0);
              par_ofstream ofs("crash_save");
              ofs << fsh << endl;
              cerr << " exiting after choleski error and parameter save" << endl;
              ad_exit(1);
            }
            
            catch (ad_float_exception & adfe)
            {
              badflag=1;
              reset_gradient_stack();

              //save current best parameter estimates
              x=xbest;
              test_adexception_flag=0;
              f=fcomp(fsh,dvar_vector(x),nvar,0,gbest,0);
              par_ofstream ofs("crash_save");
              ofs << fsh << endl;
              
              cerr << " exiting after FPE and parameter save" << endl;
              ad_exit(1);

              /*
              //_fpreset();
              signal(SIGFPE, &adfloat_except);
              //_control87(0, 0);
              delete pgs;
              pgs = new gradient_structure(gss);
              //_control87(0,MCW_EM);
              */
            }
#endif
            if (badflag)
            {
              f=1.e+95;
            }
#if defined(AD_FPEXCEPTION)
            int in_gradcalc_flag=1;
            try 
#endif
            {
              cout << "calling gradcalc" << endl;
              if (fsh.parest_flags(353)==0)
              {
                gradcalc(nvar,g);
              }
              in_gradcalc_flag=0;
            }
#if defined(AD_FPEXCEPTION)
            catch (choleski_exception & adfe)
            {
              if (in_gradcalc_flag && gradient_structure::get_save_var_flag())
              {
                gradient_structure::restore_arrays();
                gradient_structure::restore_variables();
              }
              save_best_estimate(test_adexception_flag,x,xbest,nvar,
                gbest,fsh,fmc.ifn);
            }
            catch (myexception & adfe)
            {
              //save current best parameter estimates
              test_adexception_flag=0;
              x=xbest;
              f=fcomp(fsh,dvar_vector(x),nvar,0,gbest,0);
              par_ofstream ofs("crash_save");
              ofs << fsh << endl;
              cerr << " exiting after choleski error and parameter save" << endl;
              ad_exit(1);
            }
#endif
            

            if (0)
            {
              random_number_generator rng(435);
              dvector eps(1,nvar);
              eps.fill_randn(rng);
              MY_DOUBLE_TYPE gn=norm(g);
              dvector sg=g/gn;
              
              ofstream ofs("feasible");
              // test for feasible direction
              gradient_structure::set_NO_DERIVATIVES();
              MY_DOUBLE_TYPE fc=f;
              ofs << "current f = " << setscientific() << setprecision(16) << f << " " << fc-f << endl;
              dvector y=x+1.e-6*g;
              MY_DOUBLE_TYPE f1=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
              ofs << "1.e-6 f = " << setscientific() << setprecision(16) << f << " " << fc-f1 << endl;
              /*
              y=x+1.e-9*g;
              MY_DOUBLE_TYPE f2=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
              ofs << "1.e-9 f = " << setscientific() << setprecision(16) << f << " " << fc-f2 << endl;
              y=x+1.e-10*g;
              MY_DOUBLE_TYPE f3=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
              ofs << "1.e-10 f = " << setscientific() << setprecision(16) << f << " " << fc-f3 << endl;
              y=x+1.e-12*g;
              MY_DOUBLE_TYPE f4=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
              ofs << "1.e-12 f = " << setscientific() << setprecision(16) << f << " " << fc-f4 << endl;
              y=x+1.e-13*g;
              MY_DOUBLE_TYPE f5=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
              ofs << "1.e-13 f = " << setscientific() << setprecision(16) << f << " " << fc-f5 << endl;
              */

              for (int i=1;i<=20;i++)
              {
                eps.fill_randn(rng);
                dvector newg=(sg+.05*eps/norm(eps));
                newg/=norm(newg);
                newg*=gn;
  
                y=x+1.e-9*newg;
                MY_DOUBLE_TYPE f8=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
                ofs << "1.e-9 f = " << setscientific() << setprecision(16) << f << " " << fc-f8 << endl;
                y=x+1.e-10*newg;
                MY_DOUBLE_TYPE f9=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
                ofs << "1.e-10 f = " << setscientific() << setprecision(16) << f << " " << fc-f9 << endl;
                y=x+1.e-11*newg;
                MY_DOUBLE_TYPE f10=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
                ofs << "1.e-11 f = " << setscientific() << setprecision(16) << f << " " << fc-f10 << endl;
                y=x+1.e-12*newg;
                MY_DOUBLE_TYPE f11=fcomp(fsh,dvar_vector(y),nvar,0,gbest,0);
                ofs << "1.e-12 f = " << setscientific() << setprecision(16) << f << " " << fc-f11 << endl;
              }
              ad_exit(1);
            }
       
            cout << "finished gradcalc" << endl;
//            ofssg << g << endl;   //NMD_11May2014 - turned off output
            //cout << "heap 2 = " << heapcheck() << endl;
            //gse->dep_addr=0;
            if (f < fbest)
            {
              fbest=f;
              gbest=g;
              xbest=x;
            }
          }
          if (fsh.parest_flags(333) 
            && 
             (!((fmc.itn)%(fsh.parest_flags(333)))
              || fmc.itn==1) && fmc.itn>old_itn)
          {
            old_itn=fmc.itn;
            //save current best parameter estimates
            int test_adexception_flag_save=test_adexception_flag;
            test_adexception_flag=0;
            gradient_structure::set_NO_DERIVATIVES();
            dvector gtemp(1,nvar);
            MY_DOUBLE_TYPE ftemp=fcomp(fsh,dvar_vector(xbest),nvar,0,gtemp,0);
            int suffix=0;
            if (fmc.itn==1)
              suffix=1;
            else
              suffix=fmc.itn;
            adstring file_name= adstring("intermediate_results_") 
              + str((const int)(suffix));
            par_ofstream ofs(file_name);
            ofs << fsh << endl;
            ofs << "# Objective function value" << endl << "   " << -ftemp 
                << endl;
            ofs << "# The number of parameters" << endl << "   " << nvar << endl;
            adstring command="cp " + file_name + " intermediate_results";
            if (fsh.num_tag_releases)
              ofs << "# Likelihood component for tags " <<fsh.likecomponent[0]<< endl;
            test_adexception_flag=test_adexception_flag_save;
            system(command);
            adstring command1="echo " + str((const int)fmc.itn) +  " > success_itn";
            system(command1);
       
            adstring command2="echo " + str((const int)fmc.ifn) +  " > success_ifn";
            system(command1);
            gradient_structure::set_YES_DERIVATIVES();
          }
        }
        x=xbest;

#if defined(AD_FPEXCEPTION)
        signal(SIGFPE, SIG_DFL);
#endif
      }
     // *********************************************************
     // *********************************************************
     // !!! Dave F Dec 19 02 If you want pvm used in other calculations 
     //  remove this and modify program
        mf_pvm->pvm_save_switch=mf_pvm->pvm_switch;
        mf_pvm->pvm_switch=0;
     // *********************************************************
     // *********************************************************
     // if (!hang_flag)
      {
        //cout << endl << "Calling printing fcomp" <<endl;
        gradient_structure::set_NO_DERIVATIVES();
    
        if (fsh.parest_flags(152)>0)
        {
          dvar_vector y=xscale(dvar_vector(x),grad_scale);
          f=fcomp(fsh,y,nvar,1,gbest,0);
        }
        else
        {
          if (sum(fsh.data_fish_flags(2)))
          {
            fsh.do_fishery_projections_flag=1;
          }
          int num_sims=1;
          if (fsh.do_fishery_projections_flag==1 && fsh.age_flags(20)>0) 
          {
            num_sims=fsh.age_flags(20);
            fsh.proj_output_files[0]=
               new ofstream("other_projected_stuff");
            fsh.proj_output_files[1]=
               new ofstream("projected_randomized_catches");
            fsh.proj_output_files[2]=
               new ofstream("projected_randomized_catch_at_age");
            fsh.proj_output_files[3]=
               new ofstream("projected_numbers_at_age");
            fsh.proj_output_files[4]=
               new ofstream("projected_fmort_at_age_year");    //NMD 8Nov2011
            fsh.proj_output_files[5]=
               new ofstream("projected_numbers_at_age_noeff");    //NMD 10Nov2011
            fsh.proj_output_files[6]=
               new ofstream("Fmults.txt");    //YT_2018-04-26
          }
          cifstream cif("simulated_numbers_at_age");
          if (!cif && fsh.age_flags(20)>0)  //NMD_Aug6_2018
          {
            cerr << "Error trying to open file simulated_numbers_at_age " << endl;
            ad_exit(1);
          }
          int tmp;
          cif >> tmp;
          int error_flag=0;
          if (tmp !=num_sims)
          {
            error_flag=1;
            cerr << "Error -- incorrect number of simulated recruitments"
               " in simulated_numbers_at_age" << endl;
          }
          fsh.simulated_numbers_at_age.allocate
            (1,num_sims,1,fsh.num_regions,2,fsh.nage);
          cif >> fsh.simulated_numbers_at_age;
          fsh.simulated_recruitments.allocate(1,num_sims,1,fsh.num_regions,
            fsh.last_real_year+1,fsh.nyears);
          {
      
          cif >> fsh.simulated_recruitments;
          if (!cif || error_flag)
          {
            fsh.simulated_recruitments.initialize();
          }
         
          fsh.simulated_recruitments=exp(fsh.simulated_recruitments);

          if (!cif)
          {
            cerr << "Error reading from file simulated_numbers_at_age" 
                 << endl;
          }
          }

//NMD 22Feb2012
	  dvector noeff_sw;
	  noeff_sw=column(fsh.fish_flags,55);
          if(sum(noeff_sw))
	  {
            cifstream cif("simulated_numbers_at_age_noeff");
            if (!cif && fsh.age_flags(20)>0)
            {
              cerr << "Error reading from file simulated_numbers_at_age_noeff" 
                 << endl;
              ad_exit(1);    //NMD_aug06_2018
            }
            int tmp;
            cif >> tmp;
            if (tmp !=num_sims)
            {
              cerr << "Error -- incorrect number of simulated recruitments"
                " in simulated_numbers_at_age_noeff" << endl;
            }
            fsh.simulated_numbers_at_age_noeff.allocate
              (1,num_sims,1,fsh.num_regions,2,fsh.nage);
            cif >> fsh.simulated_numbers_at_age_noeff;
	  }
          {
            const imatrix & _simyears=fsh.simyears;
            ADUNCONST(imatrix,simyears)
            if (allocated(simyears))
              simyears.deallocate();

            ifstream ifs("simyears");
            int ns=0; int year1=0; int yearlst;
            ifs >> ns;
            if (ns !=num_sims)
            {
              cerr << "Warning -- number of simulations in file simyears"
               " is different from the number of simulations being done"
               << endl;
              simyears.deallocate();
            }
            else
            {
              ifs >> year1 >> yearlst;
              simyears.allocate(1,ns,year1,yearlst);
              ifs >> simyears;
              if (!ifs)
              {
                cerr << "Did not successfuly read in information from file"
                        " simyears" << endl;
                simyears.deallocate();
              }
            }
          }

//NMD 22Feb2012

          for (int i=1;i<=num_sims;i++)
          {
            cout << "Doing simulation " << i << " of " << num_sims << endl;
            int prswitch=1;
            if (fsh.do_fishery_projections_flag==1 && fsh.age_flags(20)>0)   //NMD 3Mar2012 - allows iterative simulations
//			if (fsh.do_fishery_projections_flag==1&& num_sims>1) 
            {
              fsh.projection_sim_index=i;
              prswitch=1;
            }
            if(fsh.parest_flags(1)==1 && sum(fsh.data_fish_flags(2)))   //NMD 20Jan2012
            {
              fsh.bh_variance = fsh.save_bh_variance;
            }                                                           //NMD 20Jan2012
            stupid_print_switch=1;
//            fsh.censored_gamma_report_flag=1;
            if (fsh.parest_flags(111) >= 5 && fsh.parest_flags(111) <=8)
            {      //NMD_30Sep2019
              fsh.censored_gamma_report_flag=1;
            }
            fsh.censor_report.initialize();
            fsh.print_implicit_effort_flag=1;
            f=fcomp(fsh,dvar_vector(x),nvar,prswitch,gbest,0);
            fsh.print_implicit_effort_flag=0;
            if (fsh.censored_gamma_report_flag)  //NMD_30Sep2019
            {
              ofstream ofs("censored_gamma_report");
              ofs << "Proportion of bins censored " 
                  << fsh.censor_report(1)/sum(fsh.censor_report)
                  << " Cutoff value for censoring " << fsh.censor_report(3) 
                  << " smoothed width for transition " << fsh.censor_report(4) 
                  << endl;
              fsh.censored_gamma_report_flag=0;
              stupid_print_switch=0;
            }

            if (fsh.parest_flags(246) && fsh.generate_report)
            {
              print_indep_vars_report(gbest,fsh.indep_var,
                fsh.indep_var_lo,fsh.indep_var_hi,nvar);
	    }

            //NMD 9Nov2011
            // now call fcomp with fishing mortality turned off
            // in selected fisheries
            // right now turn off in all fisheries
            fsh.q_flag=column(fsh.fish_flags,55);
            d3_array tmp_obs_tot_catch;
            dmatrix tmp_log_effort_by_fishery; //NMD 9Aug2022
            if (sum(fsh.q_flag) && fsh.age_flags(20)>0)  //NMD 18Nov2011
            {
              ivector ff55=column(fsh.fish_flags,55);
              if (sum(ff55))
              {
                // save the real obs_tot_catch data
                tmp_obs_tot_catch.allocate(1,fsh.num_regions,1,fsh.num_fish_periods,
                    1,fsh.num_fish_incidents);
                tmp_obs_tot_catch=fsh.obs_tot_catch;
      
                // zero out the catches from fisheries that are to be turned off;
                fsh.zero_out_catches();

                // save the real log_effort_by_fishery data
                tmp_log_effort_by_fishery.allocate(1,fsh.num_fisheries,
                  1,fsh.num_fish_times);   //NMD_9Aug2022
                tmp_log_effort_by_fishery=fsh.log_effort_by_fishery;

                // zero out the effort from fisheries that are to be turned off;
                if (fsh.age_flags(92)==2 && sum(fsh.data_fish_flags(2)))
                {
                  fsh.zero_out_log_effort();
                }      

                // so this will have to be redone since the number of missing
                // catch situations may have changed
                //fsh.set_missing_totcatch_flags();
              }
              if (fsh.parest_flags(152)>0)
              {
                dvar_vector y=xscale(dvar_vector(x),grad_scale);
                //f=fcomp(fsh,y,nvar,1,gbest,0);
                int save1=fsh.af170q0ex;
                int save2=fsh.age_flags(170);
                int save3=fsh.af170q0;
                fsh.af170q0ex=1;
                fsh.af170q0=1;
                fsh.age_flags(170)=1;
                f=fcomp(fsh,y,nvar,1,gbest,0,&fsh.q_flag);
                fsh.af170q0ex=save1;
                fsh.age_flags(170)=save2;
                fsh.af170q0=save3;
              }
              else
              {
               // f=fcomp(fsh,dvar_vector(x),nvar,1,gbest,0);
                int save1=fsh.af170q0ex;
                int save2=fsh.age_flags(170);
                int save3=fsh.af170q0;
                fsh.af170q0ex=1;
                fsh.af170q0=1;
                fsh.age_flags(170)=1;
                fsh.allocate_optional_stuff();
                f=fcomp(fsh,dvar_vector(x),nvar,1,gbest,0,&fsh.q_flag);
                fsh.af170q0ex=save1;
                fsh.age_flags(170)=save2;
                fsh.af170q0=save3;
              }
              if (sum(ff55))
              {
                // get the real obs_tot_catch back
                fsh.obs_tot_catch=tmp_obs_tot_catch;
                // get the real log_effort_by_fishery back
                fsh.log_effort_by_fishery=tmp_log_effort_by_fishery; //NMD_9Aug2022
      
                //and recalculate the missing catch flags etc. 
                //
                //fsh.set_missing_totcatch_flags();
              }
              cout << "===== newl3, returned running fzero" <<endl;
              //Reset q_flag to zero
			  fsh.q_flag = 0;
            }
//NMD 9Nov2011
          }
          if (fsh.proj_output_files[0])
          {
#if defined(close)
#  undef close
            fsh.proj_output_files[0]->close();
#  define close _close
#else
            fsh.proj_output_files[0]->close();
#endif
            delete fsh.proj_output_files[0];
            fsh.proj_output_files[0]=0;
          }
          fsh.do_fishery_projections_flag=0;
        }
    
    
        {
          ofstream ofs("gradient.rpt");
          ofstream ofs1("sorted_gradient.rpt");
          ofs << nvar << endl;
          if (fsh.parest_flags(152)>0)
          {
            ofs << elem_prod(gbest,.01+grad_scale) << endl;
          }
          else
          {
            dmatrix M(1,2,1,gbest.indexmax());
            M(1).fill_seqadd(1,1);
            M(2)=fabs(gbest);
            dmatrix N=sort(trans(M),2);
            ofs << gbest << endl << endl;
            for (int i=1;i<=nvar;i++)
            {
              ofs1 << setw(10) << setw(12) << N(i) << endl;
              ofs << setw(10) << i << " " << setw(12) << gbest(i) << endl;
            }
          }
        }
        if (fmc.quit_flag=='Q')
        {
          quit_flag=1;
          //cout << "===== setting quit_flag = ";
          //cout << quit_flag << endl;
        }
        else if (fmc.quit_flag=='N')
        {
          quit_flag=2;
          //cout << "===== setting quit_flag = ";
          //cout <<"newl3.cpp " << quit_flag << endl;
        }
        else
        {
          quit_flag=0;
          //cout << "===== setting quit_flag = ";
          //cout <<"newl3.cpp " << quit_flag << endl;
        }
        if (fsh.age_flags(52))
        {
          quit_flag=1;
        }
        cout << "===== newl3.cpp,  setting quit_flag = "<< quit_flag << endl;

        // now call fcomp with fishing mortality turned off
        // in selected fisheries
        // right now turn off in all fisheries
        fsh.q_flag=column(fsh.fish_flags,55);
        d3_array tmp_obs_tot_catch;
        dmatrix tmp_log_effort_by_fishery; //NMD 9Aug2022
//        if (sum(fsh.q_flag))
        if (sum(fsh.q_flag) && fsh.age_flags(20)<1)  //NMD 10Nov2011
        {
          ivector ff55=column(fsh.fish_flags,55);
          if (sum(ff55))
          {
            // save the real obs_tot_catch data
            tmp_obs_tot_catch.allocate(1,fsh.num_regions,1,fsh.num_fish_periods,
               1,fsh.num_fish_incidents);
            tmp_obs_tot_catch=fsh.obs_tot_catch;
      
            // zero out the catches from fisheries that are to be turned off;
            fsh.zero_out_catches();

            // save the real log_effort_by_fishery data
            tmp_log_effort_by_fishery.allocate(1,fsh.num_fisheries,
              1,fsh.num_fish_times);   //NMD_9Aug2022
            tmp_log_effort_by_fishery=fsh.log_effort_by_fishery;

            // zero out the effort from fisheries that are to be turned off;
            if (fsh.age_flags(92)==2 && sum(fsh.data_fish_flags(2)))
            {
              fsh.zero_out_log_effort();
            }      
            // so this will have to be redone since the number of missing
            // catch situations may have changed
            //fsh.set_missing_totcatch_flags();
          }
          if (fsh.parest_flags(152)>0)
          {
            dvar_vector y=xscale(dvar_vector(x),grad_scale);
            //f=fcomp(fsh,y,nvar,1,gbest,0);
            int save1=fsh.af170q0ex;
            int save2=fsh.age_flags(170);
            int save3=fsh.af170q0;
            fsh.af170q0ex=1;
            fsh.af170q0=1;
            fsh.age_flags(170)=1;
            if (sum(fsh.data_fish_flags(2)))
            {
              fsh.do_fishery_projections_flag=1;
            }
            f=fcomp(fsh,y,nvar,1,gbest,0,&fsh.q_flag);
            fsh.do_fishery_projections_flag=0;
            fsh.af170q0ex=save1;
            fsh.age_flags(170)=save2;
            fsh.af170q0=save3;
          }
          else
          {
           // f=fcomp(fsh,dvar_vector(x),nvar,1,gbest,0);
            int save1=fsh.af170q0ex;
            int save2=fsh.age_flags(170);
            int save3=fsh.af170q0;
            fsh.af170q0ex=1;
            fsh.af170q0=1;
            fsh.age_flags(170)=1;
            fsh.allocate_optional_stuff();

            if (sum(fsh.data_fish_flags(2)))
            {
              fsh.do_fishery_projections_flag=1;
            }
            f=fcomp(fsh,dvar_vector(x),nvar,1,gbest,0,&fsh.q_flag);
            fsh.do_fishery_projections_flag=0;
            fsh.af170q0ex=save1;
            fsh.age_flags(170)=save2;
            fsh.af170q0=save3;
          }
          if (sum(ff55))
          {
            // get the real obs_tot_catch back
            fsh.obs_tot_catch=tmp_obs_tot_catch;
      
            // get the real log_effort_by_fishery back
            fsh.log_effort_by_fishery=tmp_log_effort_by_fishery; //NMD_9Aug2022

            //and recalculate the missing catch flags etc. 
            //
            //fsh.set_missing_totcatch_flags();
          }
        }
        cout << "===== newl3, returned printing fcomp" <<endl;
      }
      if (mf_pvm->pvm_save_switch)
      {
        mf_pvm->pvm_switch=mf_pvm->pvm_save_switch;
        mf_pvm->pvm_save_switch=0;
      }
    }
  }
}

void inverse_xscale(independent_variables& x,dvector& grad_scale)
{
  x=elem_prod(x,.1+grad_scale);
}

dvar_vector xscale(const dvar_vector& x,const dvector& grad_scale)
{
  dvar_vector y=elem_div(x,.1+grad_scale);
  return y;
}

#undef HOME_VERSION
