/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#define HOME_VERSION
#include "all.hpp"
#include <admodel.h>
extern dvector mean_weights_kludge;

extern adstring * _current_path;
extern adstring delayed_infile;
extern adstring delayed_outfile;
extern adstring full_datafile_path;
extern adstring full_input_parfile_path;
extern adstring full_output_parfile_path;
extern adstring directory_path;
  
extern int stupid_print_switch;
void send_f_to_master(const dvariable& f);
dvar_vector mfget_f_from_slaves(void);

#include "variable.hpp"
int threaded_catch_devs=0;
extern "C" void adfloat_except(int k);

int dertest_flag=0;

extern mf_pvm_manager * mf_pvm;
    //int relaxation_sequence_control(void);
    dvariable objective_function(dvar_fish_stock_history& fsh);
    dvariable objective_function(dvar_len_fish_stock_history& fsh,
      dvariable& mean_length_constraints,int print_switch,ofstream* of_pen);
  //double fcomp(dvar_len_fish_stock_history& fsh,dvar_vector& x, int nvar,
  //  int print_switch,dvector& gbest,int gradient_switch);
    void write_report(ostream& ofs,dvar_fish_stock_history& fsh);
  dvariable lhmodel_fit(dvar_len_fish_stock_history& fsh);
  void fish_mort_rep(dvar_len_fish_stock_history& fsh,
    ofstream& pof);
 int dep_gradients_calc2(dvar_len_fish_stock_history& fsh);
 int dep_gradients_calc2_for_projections(dvar_len_fish_stock_history& fsh, int noeff_sw);  //NMD 22Feb2012
// int dep_gradients_calc2_for_projections(dvar_len_fish_stock_history& fsh);
 int dep_gradients_calc2_split(dvar_len_fish_stock_history& fsh);
 int dep_gradients_calc_noeff(dvar_len_fish_stock_history& fsh);

 void ests_write(ofstream& of,dvar_len_fish_stock_history& fsh);
 void ests_write1(ofstream& of1,dvar_len_fish_stock_history& fsh,void * pq_flag);
dvector cbiocalc(dvar_len_fish_stock_history& fsh);
dvector cbiocalc_wt(dvar_len_fish_stock_history& fsh);
dvariable implicit_eff_dev_penalty(dvar_fish_stock_history& fsh,int print_switch);

extern int TURNOFF;

  void simreport(dvar_len_fish_stock_history& fsh);

  void greport(const char * s)
  {
    //cout << s << "  " <<  endl;
    //cout << s << "  " << gradient_structure::NUM_GRADSTACK_BYTES_WRITTEN() 
    //  << endl; 
  }

extern dvariable * pdfpen;
 // MY_DOUBLE_TYPE fcomp(dvar_len_fish_stock_history& fsh,dvar_vector& x, int nvar, int print_switch,dvector& gbest,int gradient_switch);
 
void check_slave_crash();
 // extern grad_stack_entry * gse;
void mfsend_x_to_slaves(const dvector&  x);
void mfsend_x_to_slaves(const dmatrix&  x);

  //double get_fmsy_pen_wt(int i)
  //{
  //  if (i==0)
  //    return(MY_DOUBLE_TYPE)( 1000);
  //  else
  //    return(MY_DOUBLE_TYPE)( i);
  //}

MY_DOUBLE_TYPE get_flag(int flg, MY_DOUBLE_TYPE deflt)
  {
    if (flg==0)
      return(MY_DOUBLE_TYPE)( deflt);
    else
      return(MY_DOUBLE_TYPE)( flg);
  }

  MY_DOUBLE_TYPE fcomp(const dvar_len_fish_stock_history& _fsh,
    const dvar_vector& _x, int nvar,int print_switch,const dvector& _gbest,
    int gradient_switch,ivector * pq_flag)
  {
    ADUNCONST(dvar_len_fish_stock_history,fsh)

    //fsh.do_fragment_stuff();

    dvar_vector& x=(dvar_vector&) _x;
    if (dertest_flag)
    {
      dvar_vector o(1,1);
      o(1)=x(1);
      dvar_vector p(1,1);
      p(1)=x(2);
      dvariable tau;
      tau=1.0+exp(x(3));

      //dvar_vector censored_gamma(const prevariable & tau,dvar_vector& o,
        //  const dvar_vector& mu,MY_DOUBLE_TYPE e)
      dvar_vector v=fsh.censored_gamma(tau,o,p,0.99);
      cout << "remove VVVV" << endl;
      return value(v(1));
    }

//    cout << "KKK called fcomp" << endl;
    //adfloat_except(1);
    adtimer mytimer;
    //gse->dep_addr=0;
  if (pq_flag)
    cout <<"===== newl5.cpp, pq_flag: " << pq_flag << endl;
  //cout << "calling heapcheck A" << endl;
  //cout << "heapcheck = " << heapcheck() << endl;
    dvector& gbest=(dvector&) _gbest;
    //dvar_len_fish_stock_history& fsh=(dvar_len_fish_stock_history&) _fsh;
    //cout << "Entered fcomp" << endl;
    char tmp_file_names[19][13]={"tmpout.1","tmpout.2","tmpout.3","tmpout.4",
      "tmpout.5","tmpout.6","tmpout.7","tmpout.8","tmpout.9","tmpout.10",
      "tmpout.11","tmpout.12","tmpout.13","tmpout.14","tmpout.15",
      "tmpout.16","tmpout.17","tmpout.18","tmpout.19"};
    dvariable f=0.0;
    dvariable ppf_tmp=0.0;   //NMD_21nov2023
    //cout << "entered fcomp" << endl;
    greport("beginning fcomp do_every"); 
    dvariable mean_length_constraints=0.0;
    fsh.generate_report=print_switch;
//    f+=fsh.reset(x);
    ppf_tmp=fsh.reset(x);
    f+=ppf_tmp;    
    if (fsh.ppstf && !pq_flag && print_switch)  //NMD_21nov2023
    {
      fsh.ppstf->reset_pen=value(ppf_tmp);
      ppf_tmp=0.0;
    }

    int send_flag=1;
    if (send_flag)
    {
      if (ad_comm::pthread_manager && fsh.threaded_tag_flag==1)
      {
        int ng1=ad_comm::pthread_manager->num_in_group(1);
        for (int i=1;i<=ng1;i++)
          fsh.send_variables_to_tag_slave2(0,i);
      }
    }

//    cout << " #####   Dbug f: after reset: " << setprecision(12) << f <<  endl;  //NMD
    if (fsh.parest_flags(197))
    {
      par_ofstream ofs(full_output_parfile_path);
            cout<<"PARFILE NAME AT BBB: "<<full_output_parfile_path<<endl;
      ofs << fsh;
      exit(0);
    }
    //cout << "constraint penalty = " << f << endl;
    dvariable ffpen=0.0;
    //  cout << " before do_every" << endl;
    ffpen=0.0;

    fsh.set_global_variance();
    fsh.Nsave=fsh.N;
    fsh.N=-100.;
   /*
    fsh.mean_lengths_by_year_for_projections();
    fsh.mean_weights_by_year_for_projections();
    */
   
    fsh.mean_lengths_by_year_calc();
    fsh.mean_weights_by_year_calc();
    
    ppf_tmp=ffpen;  //NMD_21nov2023
    
    fsh.do_everything_calc(ffpen,pq_flag);
    
    ppf_tmp=ffpen-ppf_tmp;  //NMD_21nov2023
    if (fsh.ppstf && !pq_flag && print_switch)  //NMD_21nov2023
    {
      fsh.ppstf->do_every_pen=ppf_tmp;
      ppf_tmp=0.0;
    }


    //ffpen+=100*square(fsh.pmsd->vb_coff(3)-10.);
    ppf_tmp=ffpen;
    if (fsh.parest_flags(373)>0 && fsh.parest_flags(374)>0)
    {
      ffpen+=fsh.kludged_initial_survival_equilibrium_penalty(pq_flag);
    }
    ppf_tmp=ffpen-ppf_tmp;  //NMD_21nov2023
    if (fsh.ppstf && !pq_flag && print_switch)  //NMD_21nov2023
    {
      fsh.ppstf->kludge_surv_pen=ppf_tmp;
      ppf_tmp=0.0;
    }

#if defined(USE_ADPVM)
    if (mf_pvm->pvm_switch ==0 || mf_pvm->pvm_switch ==1)
    {
      f+=ffpen;
    }
#else
      f+=ffpen;
#endif
//    cout << " #####   Dbug f: after ffpen: " << setprecision(12) << f <<  endl;  //NMD
      //cout << " after do_every" << endl;

    greport("after do_every"); 

  //cout << "calling heapcheck B" << endl;
  //cout << "heapcheck = " << heapcheck() << endl;
    if (fsh.age_flags(92))
    {
      // replace kalman filter with robust fit DF feb07 05
      int ffs66=sum(column(fsh.fish_flags,10));
      if (fsh.age_flags(104))
      {
        if (fsh.age_flags(106)==0)
        {
          //f-=fsh.robust_kalman_filter_for_catchability();
#if defined(USE_ADPVM)
          if (mf_pvm)
          {
            if (mf_pvm->pvm_switch==1)
            {
              mfsend_x_to_slaves(value(fsh.implicit_catchability));
            }
            else if (mf_pvm->pvm_switch==0)
            {
              dvariable tmp=fsh.grouped_implicit_catchability_deviations_calc();
              cout << "Catchability deviations likelihood " << tmp << endl;
              f+=tmp;
            }
          }
          else
#endif
          {
            dvariable tmp=fsh.grouped_implicit_catchability_deviations_calc();
            cout << "Catchability deviations likelihood " << tmp << endl;
            f+=tmp;
          }
        }
        else
        {
 //#if defined(__BORLANDC__)
          //fsh.grouped_implicit_catchability_deviations_calc_part1_thread();
           
   //          if (threaded_catch_devs==0)
   //          {
   //            fsh.create_catchability_deviation_thread();
   //          }
   //          pthreads_master_send_signal_to_slave();
   //          threaded_catch_devs=1;
//#else
          cerr << "NOT IMPLEMENTED" << endl;
          ad_exit(1);
//#endif
        }
      }
      if (fsh.age_flags(104)==0 ||  ffs66 !=0)
      {
        //f+=fsh.fit_implicit_catchability();
         cout << "BBBBBBBBBBBBBB" << endl;
        //f+=implicit_eff_dev_penalty(fsh,print_switch);
      }
    }
//    cout << " #####   Dbug f: after eff_devs: " << setprecision(12) << f <<  endl;  //NMD
    //if (TURNOFF==0)

    {
      if (fsh.age_flags(145))
      {
        if (mf_pvm->pvm_switch == 0 || mf_pvm->pvm_switch == 1 )
        {
          dvar_matrix m;
          if(fsh.age_flags(163) == 0)
          {
            if (fsh.pmsd==0)
            {
              dvariable tmp=
                stock_recruit_bh_steep(fsh,(ofstream *)(0),pq_flag);
              cout << "bh_steep contribution = " << tmp << endl;
              if (fsh.ppstf)
              {
                fsh.ppstf->bh_steep_contribution=value(tmp);
              }
              f+=tmp;
            }
            else
            {
              if (!sum(column(fsh.pmsd->species_flags,2)))
              {
                for (int i=1;i<=fsh.pmsd->num_species;i++)
                {
                  fsh.pmsd->current_species=i;
                  dvariable tmp=
                    stock_recruit_bh_steep(fsh,(ofstream *)(0),pq_flag);
                  cout << "species " << i 
                       << " bh_steep contribution = " << tmp << endl;
                     cout  << " fsh.p_MSY "  << fsh.p_MSY << endl;
                  if (fsh.ppstf)
                  {
                    fsh.ppstf->bh_steep_contribution(i)=value(tmp);
                  }
                  f+=tmp;
                } 
              }
              else
              {
                ivector sf2=column(fsh.pmsd->species_flags,2);
                if (fsh.pmsd->num_species !=2 || sum(sf2)!=1 )
                {
                  cerr << "only works at present for 2 species model" << endl;
                  ad_exit(1);
                }
                for (int i=1;i<=fsh.pmsd->num_species;i++)
                {
                  if (sf2(i))
                  {
                    fsh.pmsd->current_species=i;
                    dvariable tmp=
                      stock_recruit_bh_steep(fsh,(ofstream *)(0),pq_flag);
                    cout << "species " << i 
                       << " bh_steep contribution = " << tmp << endl;
                    cout  << " fsh.p_MSY "  << fsh.p_MSY << endl;
                    if (fsh.ppstf)
                    {
                      fsh.ppstf->bh_steep_contribution(i)=value(tmp);
                    }
                    f+=tmp;
                  }
                }
              }
            }
          }
          else
          {
            if (fsh.pmsd==0)
            {
              f+=stock_recruit_bh(fsh,(ofstream *)(0),m,pq_flag);
            }
            else
            {
              cerr <<  " Routine stock_recruit_bh() not yet modified for"
                  " multi-species application - exiting "  << endl;
              ad_exit(1);
            }
          }
        }
      }
    }
//    cout << " #####   Dbug f: after bh_steep: " << setprecision(12) << f <<  endl;  //NMD
    //fsh.catch_equations_calc_for_equilibrium
     // (fsh.sv/*, int nlag*/);
    if (fsh.age_flags(170))
    {
      ivector tmp=column(fsh.fish_flags,55);
      if (sum(tmp))
      {
        fsh.af170q0=1;
        //cout << "Is this a problem???" <<endl;
        //char ch;
        //cin >> ch;
        //fsh.do_everything_calc(ffpen,&tmp);
        // try resetting this flag at the bottom of the routine
        //fsh.af170q0=0;
      }
    }

    if (fsh.age_flags(150))
    {
// John H. 26/10/2001
      fsh.calculate_the_mean_weight();
//
      fsh.calculate_the_biomass();
      fsh.calculate_the_catch_biomass();
      f+=fsh.biomass_dynamics_pt();
      ofstream ofs("testyld");
      //fsh.yield_analysis_pt(&ofs);
      //fsh.yield_analysis_pt((ofstream *)(0));
    }

  //cout << "calling heapcheck C" << endl;
  //cout << "heapcheck = " << heapcheck() << endl;
    //return value(f);
    // !!!!!!!!!  TEST
    //f+=(*pdfpen);
    //f+=pow(norm2(fsh.len_dist),.5L);
    if(fsh.projection_sim_index==0)
    {
      
      if (print_switch==1) 
      {  
        if (!pq_flag) fsh.nb_tag_report_flag=1;
        ofstream * of_pen=0;
        of_pen= new ofstream("contribs");
        //stupid_print_switch=1;
        f+=objective_function(fsh,mean_length_constraints,print_switch,
          of_pen);
        //stupid_print_switch=0;
        if (mf_pvm->pvm_switch == 0 )
          cout << "Total func  " << setfixed() << setprecision(12) 
	       << setw(22) << f << endl;  //NMD
        delete of_pen;
        fsh.nb_tag_report_flag=0;
      }
      else
      {  
        ofstream * of_pen=0;
        f+=objective_function(fsh,mean_length_constraints,print_switch,
          of_pen);
        if (mf_pvm->pvm_switch == 0 )
          cout << "Total func  " << setfixed() << setprecision(12) 
               << setw(22) << f << endl;
      }
    }
    //cout << "after objective fun" << endl;
      MY_DOUBLE_TYPE tt=mytimer.get_elapsed_time();
    if (threaded_catch_devs==1)
    {
                cerr << " not implemented "<< endl;
          ad_exit(1);
      //threaded_catch_devs=0;
      //f+=fsh.grouped_implicit_catchability_deviations_calc_part2_thread_with_signal();

    }
#if defined(USE_ADPVM)
    if (mf_pvm->pvm_switch ==1)
    {
      cout << "!!!!!Need to check this out " << endl;
      if (fsh.age_flags(104))
        f+=fsh.grouped_implicit_catchability_deviations_calc_part2_pvm();
    }
#endif

    if (fsh.parest_flags(246) && !fsh.af170q0 && print_switch)
      fsh.get_indep_vars_for_report(x);

    switch (gradient_switch)
    {
    case 0:
      break;
    case 1:
      if (fsh.age_flags(145))
      {  // Kludge -- need to call this to get yield analysis
        ofstream  of1("plot-" + full_output_parfile_path + ".rep");
        //ofstream  of1("plot.rep");
        cout << "===== Printing to " << "plot-" + full_output_parfile_path + ".rep" << endl;
        ests_write1(of1,fsh,pq_flag);
      }
      // write stuff for simulation analysis
      //if (fsh.parest_flags(191))
      //{
        if (pq_flag)
        {
          cout << " Calling dep_gradients_calc_noeff" << endl;
          int num_dep_vars=dep_gradients_calc_noeff(fsh);
          cout << " Finished dep_gradients_calc_noeff" << endl;
          return(MY_DOUBLE_TYPE)( num_dep_vars) + .00001;
        }
        else
        {
          cout << " Calling dep_gradients_calc2" << endl;
          int num_dep_vars=0;
          if (fsh.parest_flags(229)==0)
          {
            num_dep_vars=dep_gradients_calc2(fsh);
          }
          else
          {
            num_dep_vars=dep_gradients_calc2_split(fsh);
          }
          cout << " Finished dep_gradients_calc2 there were " 
               << num_dep_vars << " dependent variables" << endl;
          return(MY_DOUBLE_TYPE)( num_dep_vars) + .00001;
        }
      break;    
    case 2:
      cout << " Calling dep_gradients_calc2_for_projections" << endl;
        {
          int num_dep_vars=0;
          if (fsh.parest_flags(229)==0)
          {
//NMD 22Feb2012
            if(pq_flag)
			{
//              cout << "      F_zero call to dep_gradients_calc2_for_projections" << endl;
              int noeff_sw = 1;
              num_dep_vars=dep_gradients_calc2_for_projections(fsh, noeff_sw);
			}
			else
			{
              int noeff_sw = 0;
//              cout << "      F_norm call to dep_gradients_calc2_for_projections" << endl;
              num_dep_vars=dep_gradients_calc2_for_projections(fsh, noeff_sw);
			}
//NMD 22Feb2012
          }
          else
          {
            cerr << " Illega value for pf(229)" << endl;
            ad_exit(1);
          }
          cout << " Finished dep_gradients_calc2" << endl;
          return(MY_DOUBLE_TYPE)( num_dep_vars) + .00001;
        }
      break;
    default:
        cerr << "Illegal value for gradient_switch value is " 
             << gradient_switch << endl;
    }
    
    if (print_switch==1)
    {

      int num_files=0;
      if (delayed_outfile.size())
      {
        num_files=2;
      }
      else
      {
        if (fsh.parest_flags(30))
          num_files=2;
        else
          num_files=1;
      }
      if (!pq_flag)
      {
        for (int ifiles=1;ifiles<=num_files;ifiles++)
        {
          adstring tmpout;
          if (fsh.parest_flags(30) && ifiles>1)
            tmpout=full_output_parfile_path + str(fsh.parest_flags(20));
          else 
            tmpout=full_output_parfile_path;

          if(fsh.projection_sim_index<=1)
          {
            cout << "===== Printing to "<<tmpout << endl;
            {
              par_ofstream ofs(tmpout);
              ofs << fsh;
//              ofs << "# Objective function value" << endl << "   " << -value(f) 
              ofs << "# Objective function value" << endl << "   " << value(f) 
                  << endl;    //NMD_8Mar2023
              ofs << "# The number of parameters" << endl << "   " << nvar << endl;
              if (fsh.num_tag_releases)
                ofs << "# Likelihood component for tags " <<fsh.likecomponent[0]<< endl;
              if (!fsh.age_flags(52))
              {
                ofs << "# Penalty for mean length constraints" << endl
    	              << "   "  << value(mean_length_constraints) << endl;
                ofs << "# Maximum magnitude gradient value " << endl  
                    << max(sqr(square(gbest))) << endl;
              }
              else
              {
                ofs << "# Penalty for mean length constraints" << endl
    	              << "   "  << "0.00" << endl;
                ofs << "# Maximum magnitude gradient value " << endl  
                    << "0.00"  << endl;
              }
  
              //---target data---  PK aug01-05
              if (fsh.age_flags(165))
              {
                MY_DOUBLE_TYPE targetratio=double(fsh.age_flags(165))/
                                get_flag(fsh.age_flags(164),100.);
  
                MY_DOUBLE_TYPE tmppen= get_flag(fsh.age_flags(166), 1000.);
  
                if (fsh.age_flags(167)==0)
                {
                  ofs << "# Target FFmsy= " <<  targetratio
                      << "  realized= " <<  (fsh.p_Fmmsy) 
                      << "  pen_wt= " <<  tmppen 
                      << "  rezid_pen= " <<  (fsh.p_Fmmsy_pen) <<endl;
                }
                else if(fsh.age_flags(167)==1)
                {
                  ofs << "# Target BBmsy= " <<  targetratio 
                      << "  realized= " <<  (fsh.p_BBmsy)
                      << "  pen_wt= " <<  tmppen 
                      << "  rezid_pen= " <<   (fsh.p_BBmsy_pen) <<endl;
                  ofs << "# Corresponding Fmsy/F= " <<(fsh.p_Fmmsy) <<endl;
                }
                else if(fsh.age_flags(167)==2)
                {
                  ofs << "# Target sBsBmsy= " <<  targetratio 
                      << "  realized= " <<  (fsh.p_sBsBmsy)
                      << "  pen_wt= " <<  tmppen 
                      << "  rezid_pen= " <<   (fsh.p_sBsBmsy_pen) <<endl;
                  ofs << "# Corresponding Fmsy/F= " <<(fsh.p_Fmmsy) <<endl;
                }
                else
                {
                  cout << "target: " <<fsh.age_flags(167)<< " not implemented yet"<< endl;
                  ad_exit(1);
                }
              }
              //---target data---
              //--- steepness penalty
                //if (fsh.age_flags(153))
                //{
                //  ofs << "# Target steepness= " <<  "9999"
                //      << "  realized= " <<  fsh.steepness
                //      << "  pen_wt= " <<  "9999"
                //      << "  rezid_pen= " <<  (fsh.steep_penalty) <<endl;
                //}
              //---steepness penalty
  
              //---- PK  7-02-07  target fishery impac
              if (fsh.age_flags(175)) {
                dvector ret=fsh.totalbiomass_depletion_for_plot_file();
                ofs << "# Target depletion= "<< ret(1)
                    << "  realized= " <<ret(2)
                    << "  pen_wt= " <<ret(3)
                    << "  resid_pen= "<<ret(4) <<endl;
              } 
  
              if (!fsh.pmsd || (fsh.nage ==fsh.pmsd->nage(2)) )
              {
                dvector ssum(1,fsh.nage);
                MY_DOUBLE_TYPE szsum=0.0;
                ssum.initialize();
      
                for (int ir=1;ir<=fsh.num_regions;ir++)
                {
                  szsum+=size_count(fsh.fish_mort(ir));
                  for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++) 
                  {
                    for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) 
                    {                                 
                      ssum+=exp(value(fsh.fish_mort(ir,ip,fi)));
                    }
                  }
                }
                szsum/=fsh.nage;
                ssum/=szsum;
                MY_DOUBLE_TYPE mult=double(fsh.num_fish_data_recs)/double(fsh.nyears);
                if (fsh.age_flags(37)>0)
                {
                  ofs << "# The target average fish mort is " 
                    << setprecision(4) << fsh.age_flags(37)/(1000.*mult) << endl;
                }
                ofs << "# Average fish mort per fishing incident is " 
                    << setprecision(6) << sum(ssum)/fsh.nage << endl;
                ofs << "# Average fish mort per year is " 
                  << setprecision(6) << sum(ssum)/fsh.nage*mult << endl;
                ofs << "# Average fish mort per year by age class is " << endl 
                  << "#" << setprecision(6) << ssum*mult << endl;
              }
              if (!ofs)
              {
                cerr << "===== Error printing to outfile " << endl;
              }
              else
              {
                cerr <<"===== "<< tmpout<<" finished" << endl;
              }
            }
          }
        }
      }
      delayed_outfile="";
     
      adstring fn;
      adstring fn1;
      adstring fn2;
      if (!pq_flag)
      {
        // simulate new data for length and weight frequencies
        if (fsh.parest_flags(241))
        {
          ofstream * of_pen=0;
          dvariable ffpen=0.0;
          fsh.projected_simulated_data_flags[1]=1;
          fsh.do_everything_calc(ffpen,pq_flag);
          fsh.simulate_length_and_weight_frequencies();
          fsh.generate_simulated_cpue();
        }
//        fsh.simulate_length_and_weight_frequencies();  //NMD16March2017
        fsh.projected_simulated_data_flags[1]=0;
        if (length(directory_path)>0)
        {
          fn= directory_path + adstring("\\ests.rep");
          fn1= directory_path + adstring("\\plot-" + full_output_parfile_path + ".rep");
          //fn1= directory_path + adstring("\\plot.rep");
        }
        else
        {
          fn= adstring("ests.rep");
          fn1= adstring("plot-" + full_output_parfile_path + ".rep");
          //fn1= adstring("plot.rep");
        }
        if (fsh.parest_flags(190))
        {
//          ofstream  of(fn);
          ofstream  of1(fn1);
//          ofstream ofs1("length.fit");
//          if(fsh.projection_sim_index<=1 && fsh.parest_flags(188))
          if(fsh.projection_sim_index<=1)
          {
            if (fsh.parest_flags(188))
            {
              ofstream  of(fn);
              cout << "===== Printing to " << fn << endl;
              ests_write(of,fsh);
              cout << "=====  " << fn <<" finished"<< endl;
            }
            if (fsh.parest_flags(189))
            {
              ofstream ofs1("length.fit");
              fsh.print_pred_frequencies(ofs1);
              if (fsh.nwint && !fsh.parest_flags(181))
              {
                ofstream ofs2("weight.fit");
                fsh.print_pred_wght_frequencies(ofs2);
              }
            }
          }
          cout << "===== Printing to " << fn1 << endl;
          ests_write1(of1,fsh,pq_flag);
          cout << "=====  " << fn1 <<" finished"<< endl;
        }
        // write stuff for simulation analysis
        if (fsh.parest_flags(191))
          simreport(fsh);
      }
      else
      {
        if (length(directory_path)>0)
        {
          fn= directory_path + adstring("\\estsN0.rep");
          fn1= directory_path + adstring("\\plot-" + full_output_parfile_path + ".rep");
          fn1= directory_path + adstring("\\plot.rep");
        }
        else
        {
          fn= adstring("estsN0.rep");
          fn1= adstring("plot-" + full_output_parfile_path + ".rep");
          //fn1= adstring("plot.rep");
        }
        dvar_matrix sel(1,fsh.num_fisheries,1,fsh.nage); // JH 11-Nov-03
        for (int i=1;i<=fsh.num_fisheries;i++)
        {
          int rr=fsh.realization_region(i,1);
          int rp=fsh.realization_period(i,1);
          int ri=fsh.realization_incident(i,1);
          for (int j=1;j<=fsh.nage;j++)
          {
            sel(i,j)=exp(value(fsh.fish_mort(rr,rp,ri,j)))/
                   max(exp(value(fsh.fish_mort(rr,rp,ri))));
          }
        }
        if(fsh.parest_flags(190)) F0biomass_calcs(fn1,fsh,sel);  // JH 11-Nov-03 - appends no-fishing biomass report to plot.rep
        if(fsh.parest_flags(190) && fsh.parest_flags(186)){
          fn2= directory_path + adstring("plotq0-" + full_output_parfile_path + ".rep");
          ofstream ofs((char*)(fn2));
          cout << "===== Printing to "<< fn2 << endl;
          ests_write1(ofs,fsh,pq_flag);
          cout << "=====  " << fn2 <<" finished"<< endl;
        }
      }
    }
    if (fsh.age_flags(170))
    {
      ivector tmp=column(fsh.fish_flags,55);
      if (sum(tmp))
      {
        fsh.af170q0=0;
      }
    }
    
    dvariable f1;
#if defined(USE_ADPVM)
    if (mf_pvm->pvm_switch ==1)
    {
      check_slave_crash();
   
      dvar_vector ff=mfget_f_from_slaves(); 
      f+=sum(ff);
      cout << "Total func  " << setfixed() << setprecision(2) 
           << setw(15) << f << endl;
      
    }
    if (mf_pvm->pvm_switch ==2)
    {
      send_f_to_master(f);
    }
#endif //#if defined(USE_ADPVM)
    f1=f;
    cout << "finished fcomp" << endl;
    return value(f1);
  }

#undef HOME_VERSION
