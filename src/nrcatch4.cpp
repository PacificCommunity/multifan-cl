/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#define USE_DD_NOT
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

extern int global_rmax;
extern int dumflag;

extern ofstream *  pofs_num_fish;

static void df_ddcatchimp(void)
{
  verify_identifier_string("E10");
  dvar_vector_position vpos=restore_dvar_vector_position();
  verify_identifier_string("E11");

  //d3_array_position dFdqpos=restore_d3_array_position();
  //d3_array dFdq= restore_d3_array_value(dFdqpos);
  //verify_identifier_string("E12");

  d3_array_position dsdqpos=restore_d3_array_position();
  d3_array dsdq= restore_d3_array_value(dsdqpos);
  verify_identifier_string("E13");

  dmatrix_position dMdqpos=restore_dmatrix_position();
  dmatrix dMdq= restore_dmatrix_value(dMdqpos);
  verify_identifier_string("E14");

  dmatrix_position dNdqpos=restore_dmatrix_position();
  dmatrix dNdq= restore_dmatrix_value(dNdqpos);
  verify_identifier_string("E15");

  dvar_matrix_position selpos=restore_dvar_matrix_position();
  verify_identifier_string("E16");

  dvar_vector_position Npos=restore_dvar_vector_position();
  verify_identifier_string("E17");

  dvar_vector_position Mpos=restore_dvar_vector_position();
  verify_identifier_string("E18");

  dmatrix_position dwdqpos=restore_dmatrix_position();
  dmatrix dwdq= restore_dmatrix_value(dwdqpos);
  verify_identifier_string("E19");

  dvar_vector_position wpos=restore_dvar_vector_position();
  verify_identifier_string("E20");
  dvector dfu=restore_dvar_vector_derivatives(vpos);

  //dmatrix dfF=dFdq*dfu;
  dmatrix dfs=dsdq*dfu;
  dvector dfN=dNdq*dfu;
  dvector dfM=dMdq*dfu;
  dvector dfw=dwdq*dfu;
  dfw/=1000.;    // weight is divided by 1000.

  if (pofs_num_fish)
  {
    (*pofs_num_fish) << " norm2(dfu)  norm2(dfs)   norm2(dfN)  norm2(dfM) " 
      <<  norm2(dfu) << "  "  << norm2(dfs) << "  " <<  norm2(dfN) 
      << "  " << norm2(dfM) << endl; 
  }
  if (norm2(dfs)> 1.e+30 || norm2(dfN) > 1.e+30 || norm2(dfM) > 1.e+30)
  {
    cout << "B1 posible overflow" << endl;
  } 
  dfs.save_dmatrix_derivatives(selpos);
  dfN.save_dvector_derivatives(Npos);
  dfM.save_dvector_derivatives(Mpos);
  dfw.save_dvector_derivatives(wpos);
}

MY_DOUBLE_TYPE dvar_fish_stock_history::get_disagg_catch
  (int ir, int ip, int fi, int wtf, dvar_vector enf)
{
  /*
  cout << "ir: " << ir << endl;
  cout << "ip: " << ip << endl;
  cout << "fi: " << fi << endl;
  // Need to identify the species associated with the ir, ip ,fi
  int is=pmsd->region_species_pointer(ir);
  int pi=parent(ir,ip,fi);
  int n=pmsd->fisn(ir,ip,fi);   //number of species for this incident
  for (int i=1;i<=n;i++)
  {
    int rr=pmsd->reg_in_catch(ir,ip,fi,i);
    cout << "Mirror i: " << i << " ir: " << rr << " obs_tot_catch:  "
         << obs_tot_catch(rr,ip,fi) << endl;
  }
  */
  int is=pmsd->region_species_pointer(ir);
  int pi=parent(ir,ip,fi);
  int n=pmsd->fisn(ir,ip,fi);   //number of species for this incident
  MY_DOUBLE_TYPE disagg_obs_tot_catch=0.0;

  dvector sp_ratio(1,n);
  dmatrix vuln_spp(1,n);
  sp_ratio.initialize();
  vuln_spp.initialize();
  MY_DOUBLE_TYPE tmp=0.0;
  // Vulnerable sex-specific numbers/weights at age
  for (int i=1;i<=n;i++)
  {
    int rr=pmsd->reg_in_catch(ir,ip,fi,i);
    int ng=get_nage_region(rr);
    vuln_spp(i).allocate(1,ng);
    vuln_spp(i).initialize();
    dvector sel=value(mfexp(incident_sel(rr,ip,fi)));
    for (int j=1;j<=ng;j++)
    {
      if (wtf==0)
      {
        vuln_spp(i,j)=sel(j)*value(enf(j));
      }
      else if (wtf==1)
      {
        vuln_spp(i,j)=sel(j)*value(enf(j))*
                        value(mean_weight(ir,ip,fi,j));
      }
    }
  }
  MY_DOUBLE_TYPE tot_vuln=0.0;
  for (int i=1;i<=n;i++)
  {
    int rr=pmsd->reg_in_catch(ir,ip,fi,i);
    int ng=get_nage_region(rr);
    // summate over age
    for (int j=1;j<=ng;j++)
    {
      tot_vuln+=vuln_spp(i,j);
    }
  }
  for (int i=1;i<=n;i++)
  {
    int rr=pmsd->reg_in_catch(ir,ip,fi,i);
    int ng=get_nage_region(rr);
    MY_DOUBLE_TYPE tmp=0.0;
    // summate over age
    for (int j=1;j<=ng;j++)
    {
      tmp+=vuln_spp(i,j);
    }
    sp_ratio(i)=tmp/tot_vuln;
  }

// apply sex ratio to observed aggregated total catch
  disagg_obs_tot_catch=sp_ratio(is)*obs_tot_catch(ir,ip,fi);

  return disagg_obs_tot_catch;
}


int dvar_fish_stock_history::get_num_for_nr(int ir,int ip)
{
  int fi;
  ivector ff92=column(fish_flags,92);
  ivector ff55=column(fish_flags,55);
  int nmiss=missing_catch_for_period_flag(ir,ip); 
  int leave_out=0;
  int nfi=num_fish_incidents(ir,ip);
  for (fi=1;fi<=nfi;fi++)
  {
    int i=parent(ir,ip,fi);
    int nc=1;
    if (allocated(fml_columns))
    {
      nc=fml_columns(i).indexmax();
    }
    if (ff92(i) || ( nc>0  && (missing_catch_for_incident_flag(ir,ip,fi)==1)))
       // ||  ( af170q0 && ff55(i)) )) )
    {
      leave_out++;
    }
  }
  int nfi1=nfi-leave_out;
  return nfi1;
}

int dvar_fish_stock_history::get_non_zero_effort(int ir,int ip)
{
  int fi;
  ivector ff92=column(fish_flags,92);
  ivector ff55=column(fish_flags,55);
  int nonzero_eff=0;
  int nfi=num_fish_incidents(ir,ip);
  for (fi=1;fi<=nfi;fi++)
  {
    int i=parent(ir,ip,fi);
    if (ff92(i)==0)   // not a survey sample
    {
      int ft=fish_times(ir,ip,fi);
      if (log_effort_by_fishery(i,ft) > -19.5) nonzero_eff++;
    }
  }  
//  cout << nonzero_eff << endl;
  return nonzero_eff;
}

dvar_vector dvar_fish_stock_history::ddnrf_interface(int ir,int ip,
  int nfi1,ivector& wtflag)
{
  int fi;

  MY_DOUBLE_TYPE beta=0.6;
  dvar_vector enf;
  if (!af170q0)
  {
    enf=mfexp(num_fish(ir,ip));
  }
  else
  {
    enf=mfexp(num_fish_q0(ir,ip));
  }

  if (sum(value(enf))< 1.e-10)
  {
    cout << "ddnrf_interface() - population crashed in Newton-Raphson" << endl;
    ad_exit(1);
  }

  // these are dependent variables that we don't need ?
  dvar_matrix sel(1,nfi1,1,nage);
  ivector local_wtflag(1,nfi1);
  dvar_vector other_mort(1,nage);
  other_mort=e_nat_mort_by_period(ir,ip);

  dvector otc(1,nfi1);
  int ii=0;
  int nfi=num_fish_incidents(ir,ip);
  ivector ff55=column(fish_flags,55);
  ivector ff92=column(fish_flags,92);
  ivector ff29=column(fish_flags,29);
  for (fi=1;fi<=nfi;fi++)
  {
    int i=parent(ir,ip,fi);

    if (ff92(i)==0)   // not a survey sample
    {
      if (missing_catch_for_incident_flag(ir,ip,fi)==0)
        // && ( !af170q0 ||!ff55(i)) )
      {
        // normal case with catch
        sel(++ii)=mfexp(incident_sel(ir,ip,fi));
//        otc(ii)=obs_tot_catch(ir,ip,fi);
        local_wtflag(ii)=wtflag(fi);
        if (!pmsd || pmsd->fisn(ir,ip,fi)==1) //NMD 11feb2025
        {
          otc(ii)=obs_tot_catch(ir,ip,fi);
        }
        else if (pmsd->fisn(ir,ip,fi)>1)
        {
          int wtf=local_wtflag(ii);
          otc(ii)=get_disagg_catch(ir,ip,fi,wtf,enf);
        }
      }
      else
      {
        // if (! (af170q0 && ff55(i) ))
        {
          if (parest_flags(377)==0 && parest_flags(378)==0)
          {
//            cout << " HERE!" << endl;
            if (eff_proj_fshry(i)>0)
            {
              dvariable log_fmlevel;
              dvariable fpen=0.0;
              int rmnth=really_true_month(ir,ip);
              for (int fnyri=1;fnyri<=fshtms_finyr(i);fnyri++)
              {
                int pmnth=ptr_fml_Mbs_mnth_finyr(i,fnyri);
                if (rmnth==pmnth)
                {
                  MY_DOUBLE_TYPE leff;
                  leff=log_effort_by_fishery(i,fish_times(ir,ip,fi));
                  log_fmlevel=q_level_finyr(i,fnyri)+leff;
                  break;
                }
              }
              log_fmlevel=0.2-posfun(0.2-log_fmlevel,.01,fpen);
              fm_level(ir,ip,fi)=log_fmlevel;
              //cout << log_fmlevel << " " << exp(log_fmlevel)  << endl;
              dvar_vector logfm=log_fmlevel+incident_sel(ir,ip,fi);
              //cout << "logfm" << endl << logfm << endl << " exp(logfm)" 
              //     << endl   << exp(logfm) << endl;
              other_mort+= exp(log_fmlevel+incident_sel(ir,ip,fi));
              for (int ii=other_mort.indexmin();ii<=other_mort.indexmax();ii++)
              {
                 if (fabs(value(other_mort(ii)))>100.0)
                 {
                   cout << "ddnrf_interface() - other_mort too big" << endl;
                   ad_exit(1);
                 }
              }
            }
          }
          else
          {
          // catch is missing and not set to zero in q0 case so use 
          // regression estimates to produce  a predicted fish mort //
          //which is added to other mortality
            int nc=fml_columns(i).indexmax();
            if (nc>0)
            {
              dvector dv=fml_designvec(ir,ip,fi);
              if (!allocated(dv))
              {
                cerr << "Projection incident has no matching incident in "
                      << "terminal year; fishery: " << i << endl;
                cerr << "Terminating execution" << endl;
                ad_exit(1);
              }
              ivector num_ifmlrp=implicit_fml_bounds(3);
              dvar_vector ests=implicit_fm_level_regression_pars(i)
                (1,num_ifmlrp(i));
              MY_DOUBLE_TYPE leff=log_effort_by_fishery(i,fish_times(ir,ip,fi));
      
              int gi=i;
              if (ff29(i))
              {
                gi=ff29(i);
              }
    
              dvariable tmp2;
              if (!parest_flags(396))
              {
                tmp2 = dv * inv(fml_R(gi)) * ests;
              }
              else
              {
                tmp2 = dv * ests;
              }
              dvariable log_fmlevel= tmp2+leff;
              dvariable fpen=0.0;
              log_fmlevel=0.2-posfun(0.2-log_fmlevel,.01,fpen);
              fm_level(ir,ip,fi)=log_fmlevel;
              //cout << log_fmlevel << " " << exp(log_fmlevel)  << endl;
              dvar_vector logfm=log_fmlevel+incident_sel(ir,ip,fi);
              //cout << "logfm" << endl << logfm << endl << " exp(logfm)" 
              //     << endl   << exp(logfm) << endl;
              other_mort+= exp(log_fmlevel+incident_sel(ir,ip,fi));
              for (int ii=other_mort.indexmin();ii<=other_mort.indexmax();ii++)
              {
                 if (fabs(value(other_mort(ii)))>100.0)
                 {
                   cout << "ddnrf_interface() - other_mort too big" << endl;
                   ad_exit(1);
                 }
              }
            }
          }
        }
        /*
        else
        {
          // catch is set to zero so set fm_level so that exp(fm_level)
          // is close to zero
          fm_level(ir,ip,fi)=-50.0;
        }
        */
      }
    }
    else if (ip>num_real_fish_periods(ir) && do_fishery_projections_flag==1)
    {
// catch is negligible so set fm_level so that exp(fm_level) NMD_12mar2025
// is close to zero
      fm_level(ir,ip,fi)=-15.0;
    }
  }

  dvar_vector w;
  
  w=mean_weight(ir,ip,1);  
  if (sum(column(fish_flags,11) ))
  {
    cerr << "We have fishery dependent mean lengths at age which are"
      " not supported in catch conditioned NR  " << endl;
    ad_exit(1);
  }
  ddnrf ddnrtester(nfi1,nage,value(sel),value(other_mort),
    value(enf),otc,value(w),local_wtflag,beta,
    this,ir,ip,Zmax_fish);

  if (nfi1>0)  //NMD_12mar2025 - only if non-zero catches present
  {
    ddnrtester.testnr();
  }
  if (0)
  {
     ofstream ofs("logfile1");
 
     ofs << "ddnrtester.S" << endl;
     ofs << ddnrtester.get_S() << endl;
     ofs << "ddnrtester.q" << endl;
     ofs << ddnrtester.get_q() << endl;
  }
 

  dvector cv(1,nfi1);
  int ii1=0;
  for (fi=1;fi<=nfi;fi++)
  {
    int i=parent(ir,ip,fi);
    if (ff92(i)==0)   // not a survey sample
    {
      if (missing_catch_for_incident_flag(ir,ip,fi)==0 )
         // && ( !af170q0 ||!ff55(i)) )
      {
        ii1++;
  #if (defined(USE_DD) &&((defined(__MSVC32__) &&  __MSVC32__ >=8) || defined(linux) || defined(__ADMING__)) ||defined(NO_MY_DOUBLE_TYPE) )
        cv(ii1)=to_double(ddnrtester.get_q()(ii1));
  #else
        cv(ii1)=ddnrtester.get_q()(ii1);
  #endif
      }
    }
  }
  dvar_vector v=nograd_assign(cv);


  if (ddnrtester.icount>8)
  {
    cout << "D icount = " << ddnrtester.icount << endl;
  }
  //d3_array & dFdq=make_d3_array(ddnrtester.get_dFdq());
#if defined(USE_DD) || defined(NO_MY_DOUBLE_TYPE)
  const d3_array & dsdq=make_d3_array(ddnrtester.get_dsdq());
  const dmatrix & dNdq=make_dmatrix(ddnrtester.get_dNdq());
  const dmatrix & dMdq=make_dmatrix(ddnrtester.get_dMdq());
  const dmatrix & dwdq=make_dmatrix(ddnrtester.get_dwdq());
#else
  const d3_array & dsdq=ddnrtester.get_dsdq();
  const dmatrix & dNdq=ddnrtester.get_dNdq();
  const dmatrix & dMdq=ddnrtester.get_dMdq();
  const dmatrix & dwdq=ddnrtester.get_dwdq();
#endif
//    save_identifier_string("E20");
  const char * str1;
  str1="E20";
  char* strx1=const_cast <char*> (str1);
  save_identifier_string(strx1);
  w.save_dvar_vector_position();

//    save_identifier_string("E19");
  const char * str2;
  str2="E19";
  char* strx2=const_cast <char*> (str2);
  save_identifier_string(strx2);
  dwdq.save_dmatrix_value();
  dwdq.save_dmatrix_position();
  
//    save_identifier_string("E18");
  const char * str3;
  str3="E18";
  char* strx3=const_cast <char*> (str3);
  save_identifier_string(strx3);
  //e_nat_mort_by_period(ir,ip).save_dvar_vector_position();
  other_mort.save_dvar_vector_position();
  
//    save_identifier_string("E17");
  const char * str4;
  str4="E17";
  char* strx4=const_cast <char*> (str4);
  save_identifier_string(strx4);
  enf.save_dvar_vector_position();
  
//    save_identifier_string("E16");
  const char * str5;
  str5="E16";
  char* strx5=const_cast <char*> (str5);
  save_identifier_string(strx5);
  sel.save_dvar_matrix_position();
  
//    save_identifier_string("E15");
  const char * str6;
  str6="E15";
  char* strx6=const_cast <char*> (str6);
  save_identifier_string(strx6);
  dNdq.save_dmatrix_value();
  dNdq.save_dmatrix_position();

//    save_identifier_string("E14");
  const char * str7;
  str7="E14";
  char* strx7=const_cast <char*> (str7);
  save_identifier_string(strx7);
  dMdq.save_dmatrix_value();
  dMdq.save_dmatrix_position();

//    save_identifier_string("E13");
  const char * str8;
  str8="E13";
  char* strx8=const_cast <char*> (str8);
  save_identifier_string(strx8);
  dsdq.save_d3_array_value();
  dsdq.save_d3_array_position();

  //save_identifier_string("E12");
  //dFdq.save_d3_array_value();
  //dFdq.save_d3_array_position();

//    save_identifier_string("E11");
  const char * str10;
  str10="E11";
  char* strx10=const_cast <char*> (str10);
  save_identifier_string(strx10);
  v.save_dvar_vector_position();

//    save_identifier_string("E10");
  const char * str11;
  str11="E10";
  char* strx11=const_cast <char*> (str11);
  save_identifier_string(strx11);
  gradient_structure::GRAD_STACK1->set_gradient_stack(df_ddcatchimp);
  return v;
} 
