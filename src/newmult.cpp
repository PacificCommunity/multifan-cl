/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
//ofstream check_file("check.tmp");
#ifndef __GNUC__
//#include <mf_menu.h>
#endif

#ifndef CLOGF_TRACE
  #define CLOGF_TRACE
#endif

extern int _file_version_no;
extern int sex_flag;
extern int  _ini_version_no;

#ifdef __MSC__
  MY_DOUBLE_TYPE mfexp(const MY_DOUBLE_TYPE& );
  dvariable mfexp(const prevariable& );
  dvar_vector mfexp(const dvar_vector& );
  dvector mfexp(const dvector& );
#else
  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
  dvariable mfexp(prevariable& );
  dvar_vector mfexp(const dvar_vector& );
  //dvector mfexp(dvector& );
#endif

        d3_array * ptarr1=0;

  extern int _NUMSV;
  const MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;

  void set_default_age_flags(ivector& af)
  {
    af(177)=1;
  }

  void check_number(cifstream& cif,int n,dvar_fish_stock_history & fsh)
  {
    ivector & pf=fsh.parest_flags;
    if (_file_version_no>1054 && pf(197)==0  
      && fsh.historical_parest_flags(195))
    {
      int m=0;
      cif >> m;
      if (n !=m)
      {
         cerr << "error in file read" << endl;
         cerr << " looking for " << n << "  found " << m << endl;
         ad_exit(1);
      }
    }
  }
  void check_number(cifstream& cif,int n,ivector pf)
  {
    if (_file_version_no>1054 && pf(197)==0  && pf(195))
    {
      int m=0;
      cif >> m;
      if (n !=m)
      {
         cerr << "error in file read" << endl;
         cerr << " looking for " << n << "  found " << m << endl;
         ad_exit(1);
      }
    }
  }

  int operator == (const ivector& v1,const ivector& v2)
  {
    int mmin=v1.indexmin();
    int mmax=v1.indexmax();
    for (int i=mmin;i<=mmax;i++)
    {
      if (v1(i) != v2(i))
        return 0;
    }
    return 1;
  }



  dvar_fish_stock_history::dvar_fish_stock_history(int ntg,int nregions, int ng,ivector& nfp,
    imatrix& nfi,int nfsh,int nyrs,ivector& fl ,ivector& par_fl, ivector& nft,
    ivector& _regmin,ivector& _regmax,ivector& _dataswitch,imatrix& _Dflags,
    movement_info& _mo,int _direction_flag,int _mfactor,
    imatrix& _ses_reg_flags,pmulti_species_data & _pmsd,
    ivector& _nage_by_region,ivector& _nage_by_fishery)
  :  
    ff92sum(0),
    initial_orthp_estimate_flag(0),
    Xchk(0,10,1,nregions),
    iloop(0),
    loop_flag(0),
    Zmax_fish(0.0),
    pfml_group(0),
    pgroup_manager_1(0),
    print_implicit_effort_flag(0),
    implicit_flag(0),
    censored_gamma_report_flag(0),
    censor_report(1,4),
    nb_tag_report_flag(0),
    num_regions(nregions),
    first_unfixed_year(1),
    first_unfixed_fish_time(1,nfsh),
    nage_by_region(_nage_by_region),
    nage_by_fishery(_nage_by_fishery),
    pcsa(0),
    pccsa(0),
    simulation_seeds(1,10),
    average_effort(1,nfsh),
    projected_simulated_data_flags(1,10),
    len_sample_size(1,nregions,1,nfp,1,nfi),
    wght_sample_size(1,nregions,1,nfp,1,nfi),
    max_group_ptr(100),
    group_ptr(1,max_group_ptr),
    inv_group_ptr(1,max_group_ptr),
    bl_sel_scaling(10.),
    get_annual_recruitment_flag(0),
    na(0),
    sum_ff71_flag(0),
    sum_ff48_flag(0),
    sum_ff79_flag(0),
    blbsel(1,nfsh),
    set_global_vars_flag(0),
    pmsd(_pmsd),
    ff263flag(0),
    tag_fish_switch(0),
    num_blocks(1,nfsh),
    q_flag(1,nfsh),
    annual_phi(0),
    ppstf(0),
    parest_flags(fl) ,
    thread_f(0) ,
    imppen(0),
    ss2_flag(0),
    have_projection_periods_flag(0) ,
    direction_flag(_direction_flag) ,
    age_flags(1,200) ,
    ses_reg_recr_flags(_ses_reg_flags),
    mo(_mo) ,
    num_ses_reg_active(0) ,
    sv(1,_NUMSV) ,
    year_flags(1,10,1,nyrs),
    num_seasons(_mfactor),
    season_flags(1,10,1,_mfactor),
    num_tag_releases(ntg) ,
    pop_delta(1,nregions),
    epop_delta(1,nregions),
    mean_weight_yr(1,nregions,1,nyrs,1,_nage_by_region),
    mean_weight_yr_alternative(1,nregions,1,nyrs,1,_nage_by_region),
    mean_weight_yr_proj(1,12,1,ng),
    do_fishery_projections_flag(0) ,
    totalcatch_by_numbers(1,nregions) ,
    obstotalcatch_by_numbers(1,nregions) ,
    numtotalcatch_by_numbers(1,nregions) ,
    totalcatch_by_weight(1,nregions) ,
    obstotalcatch_by_weight(1,nregions) ,
    numtotalcatch_by_weight(1,nregions) ,
    old_log_lambda(0),
    effective_len_size(1,nfsh),
    missing_catch_counter(1,nfsh),
    projection_year(1,nfsh),
    projection_month(1,nfsh),
    selmean(1,nfsh),
    ageselmean(1,nfsh),
    fishery_regions(1,nfsh),
    effective_weight_size(1,nfsh),
    gml(1,ng+1) ,
    min_tag_year(0),
    biomass(1,nyrs),
    adult_biomass(1,nyrs),
    biomass_by_region(1,nregions,1,nyrs),
    rel_biomass_by_region(1,nregions,1,nyrs),
    //region_rec_diff_colmeans(1,nregions),
    catch_biomass(1,nyrs),
    catch_biomass_by_fishery(1,nyrs,1,nfsh),
    effort_weight_by_fishery(1,nfsh,1,nft),
    fish_flags(1,nfsh,1,100) ,
    old_fish_flags(1,nfsh,1,100) ,
    num_fish_incidents(nfi) ,
    fishing_period(1,sum(nfi)) ,
    recruitment_period(1,nregions,1,nfp) ,
    regmin(_regmin),
    regmax(_regmax),
    D(1,nregions,1,nregions),
    Dad(1,_mo.num_periods(),1,ng,1,nregions,1,nregions),
    Dad2(1,_mo.num_periods(),1,ng,1,nregions,1,nregions),
    Dflags(_Dflags),
    diff_coffs(1,_mo.num_periods(),1,max(1,sum(_Dflags))),
    xdiff_coffs(1,_mo.num_periods(),1,max(1,sum(_Dflags))),
    y1diff_coffs(1,_mo.num_periods(),1,max(1,sum(_Dflags))/2),
    y2diff_coffs(1,_mo.num_periods(),1,max(1,sum(_Dflags))/2),
    zdiff_coffs(1,_mo.num_periods(),1,max(1,sum(_Dflags))),
    diff_coffs2(1,_mo.num_periods(),1,max(1,sum(_Dflags))),
    diff_coffs3(1,_mo.num_periods(),1,max(1,sum(_Dflags))),
    diff_coffs_prior(1,_mo.num_periods(),1,max(1,sum(_Dflags))),
    diff_coffs2_prior(1,_mo.num_periods(),1,max(1,sum(_Dflags))),
    diff_coffs3_prior(1,_mo.num_periods(),1,max(1,sum(_Dflags))),
    N(1,nregions,1,nyrs,1,_nage_by_region),
    Nsave(1,nregions,1,nyrs,1,_nage_by_region),
    exp_N(1,nregions,1,nyrs,1,_nage_by_region),
    fishing_incident(1,sum(nfi)) ,
    fishing_region(1,nfsh) ,
    num_fish_times(1,nfsh) ,
    num_real_fish_times(1,nfsh) ,
    realization_period(1,nfsh,1,nft) ,
    fishery_realization_index(1,nregions,1,nfp,1,nfi) ,
    missing_effort_by_region_flag(1,nregions,1,nfp,1,nfi) ,
    num_missing_effort_by_region(1,nregions,1,nfp) ,
    num_present_effort_by_region(1,nregions,1,nfp) ,
    sseason(1,nregions,1,nfp,1,nfi) ,
    bblock(1,nregions,1,nfp,1,nfi) ,
    lbsel(1,nfsh,1,_nage_by_fishery),
    realization_incident(1,nfsh,1,nft) ,
    true_effort_by_fishery(1,nfsh,1,nft),
    really_true_effort(1,nregions,1,nfp,1,nfi),
    log_true_effort_by_fishery(1,nfsh,1,nft),
    missing_effort_flag(1,nfsh),
    missing_effort_by_realization_flag(1,nfsh,1,nft),
    effort_by_fishery(1,nfsh,1,nft),
    log_effort_by_fishery(1,nfsh,1,nft),
    normalized_log_effort_by_fishery(1,nfsh,1,nft),
    log_pred_effort_by_fishery(1,nfsh,1,nft),    //NMD_20may2021
    realization_region(1,nfsh,1,nft) ,
    header_record_index(1,nregions,1,nfp,1,nfi) ,
    sel_seasons(1,nregions,1,nfp,1,nfi) ,
    year(1,nregions,1,nfp) ,
    movement_period(1,nregions,1,nfp) ,
    move_index(1,nregions,1,nfp) ,
    region_pars(1,100,1,nregions) ,
    current_biomass_by_year(1,nyrs,1,nregions),
    region_flags(1,10,1,nregions) ,
    region_rec_diffs(1,nyrs,1,nregions),
    rec_delta(1,nyrs),
    region_rec_diff_coffs(1,nyrs,1,nregions),
    region_rec_diff_sums(1,nyrs),
    fraction(1,nregions,1,nfp) ,
    rec_times(1,12) ,
    rec_covars(1,nyrs) ,
    tmprecr(1,nyrs) ,
    avail_coff(1,ng) ,
    orth_recr(1,nregions,0,nyrs-1) ,
    initpop(1,nregions,2,ng) ,
    tmpinitpop(2,ng) ,
    num_fish(1,nregions,1,nfp,1,_nage_by_region) ,
    num_fish0(1,nregions,1,nfp,1,_nage_by_region) ,
    catch(1,nregions,1,nfp,1,nfi,1,_nage_by_region),
    pred_totcatch(1,nregions,1,nfp,1,nfi),
    mean_weight(1,nregions,1,nfp,1,nfi,1,_nage_by_region),
    exp_catch(1,nregions,1,nfp,1,nfi,1,_nage_by_region),
    tot_catch(1,nregions,1,nfp,1,nfi) ,
    obs_catch(1,nregions,1,nfp,1,nfi,1,_nage_by_region),
    prop(1,nregions,1,nfp,1,nfi,1,_nage_by_region),
    catch_numbers(1,nyrs),
    incident_sel(1,nregions,1,nfp,1,nfi,1,_nage_by_region),
    mean_incident_sel(1,nregions,1,nfp,1,nfi),
    between_times(1,nfsh,2,nft) ,
    imp_back(1,nfsh,1,nft) ,
    imp_forward(1,nfsh,1,nft) ,
    fishery_sel(1,nfsh,1,_nage_by_fishery),
    selcoff(1,nfsh,1,_nage_by_fishery),
    ageselcoff(1,nfsh,1,_nage_by_fishery),
    survival(1,nregions,1,nfp,1,_nage_by_region) ,
    actual_init(1,ng) ,
    actual_recruit(1,nyrs) ,
    fish_mort(1,nregions,1,nfp,1,nfi,1,_nage_by_region) ,
    fish_mort_calcs(1,nregions,1,nfp,1,nfi,1,_nage_by_region) ,
    tot_mort(1,nregions,1,nfp,1,_nage_by_region) ,
    e_nat_mort_by_period(1,nregions,1,nfp,1,_nage_by_region) ,
    nat_mort(1,nyrs,1,ng) ,
    catchability(1,nregions,1,nfp,1,nfi) ,
    FF(1,nregions,1,nfp,1,nfi) ,
    effort_normalization_factor(1,nregions,1,nfp,1,nfi),
    effort(1,nregions,1,nfp,1,nfi),
    effort_weight(1,nregions,1,nfp,1,nfi),
    fm_level(1,nregions,1,nfp,1,nfi),
    parent(1,nregions,1,nfp,1,nfi),
    fish_times(1,nregions,1,nfp,1,nfi),
    month(1,nregions,1,nfp),
    true_year(1,nregions,1,nfp),
    true_month(1,nregions,1,nfp),
    really_true_year(1,nregions,1,nfp),
    really_true_month(1,nregions,1,nfp),
    true_week(1,nregions,1,nfp) ,
    week(1,nregions,1,nfp) ,
    q0(1,nfsh) ,
    q0_miss(1,nfsh) ,
    q1(1,nfsh) ,
    effort_dev_coffs(1,nfsh,1,nft) ,
    zero_effdev_flag(1,nfsh,1,nft) ,
    zero_catch_for_period_flag(1,nregions,1,nfp), 
    zero_catch_for_incident_flag(1,nregions,1,nfp,1,nfi), 
    missing_catch_by_fishery_flag(1,nfsh), 
    missing_catch_for_realization(1,nfsh,1,nft), 
    missing_catch_for_period_flag(1,nregions,1,nfp), 
    missing_catch_for_incident_flag(1,nregions,1,nfp,1,nfi), 
    catch_dev_coffs(1,nfsh,2,nft) ,
    implicit_catchability(1,nfsh,1,nft),
    implicit_catchability_ccond(1,nfsh,1,nft),
    new_fm_level(1,nfsh,1,nft),
    fish_pars(1,50,1,nfsh) ,
    implicit_fm_level_regression_pars(1,nfsh,1,25),
    seasonal_catchability_pars(1,nfsh,1,12) ,
    seasonal_catchability_pars_index(1,12,1,2),
    seasonal_catchability_pars_mix(1, 12,1,2),
    age_pars(1,10,1,ng) ,
    biomass_index(1,nfsh,1,nft),
    effort_devs(1,nregions,1,nfp,1,nfi) ,
    sel_devs(1,nregions,1,nfp,1,nfi,1,_nage_by_region) ,
	 corr_wy(1,nfsh) ,
    corr_by(1,nfsh) ,
    corr_wc(1,nfsh) ,
    corr_eff(1,nfsh) ,
    first_survey_time(1,nfsh),
    obs_tot_catch(1,nregions,1,nfp,1,nfi) ,
    obs_region_tot_catch(1,nregions,1,nfp) ,
    obs_prop(1,nregions,1,nfp,1,nfi,1,_nage_by_region),
    pmature(1,ng) ,
    num_fish_periods(nfp),
    num_real_fish_periods(1,nregions),
    total_num_obs(0),
    zeroed_fisheries_flag(0),
    af170q0(0),
    af170q0ex(0),
    recmean(0),
    initmean(0),
    pfmin1(0),
    pamin1(0.8),
    pamin2(0.9),
    projection_sim_index(0),
    sel_dev_coffs(1,nfsh,1,nft,1,ng),
    lagrange_c(1,nregions,1,nfp,1,nfi),
    lagrange_lambda(1,nregions,1,nfp,1,nfi),
    other_lagrange_lambda(1,100),
    other_lagrange_mu(1,100),
    other_lagrange_c(1,100),
    no_lagrangian_update(0),
    num_new_weights(0),
    newfpen(0.0),
    historical_version(0),
    current_version(0)

  {
    tag_fish_mort_calcs.allocate(1,1);
    tag_fish_mort_calcs(1).
      allocate(1,nregions,1,nfp,1,nfi,1,_nage_by_region);
    first_survey_time.initialize();
    kludged_selmean_square=0.0;
    Xchk.initialize();
    fml_designvec.allocate(1,nregions,1,nfp,1,nfi);
    fml_Mbs_finyr.allocate(1,nfsh);  //NMD_11Jul2022
    fshtms_finyr.allocate(1,nfsh);  //NMD_11Jul2022
    eff_proj_fshry.allocate(1,nfsh);  //NMD_11Jul2022
    eff_proj_fshry.initialize();  //NMD_15Jul2022
    catch_proj_fshry.allocate(1,nfsh);
    catch_proj_fshry.initialize();
    q_level_finyr.allocate(1,nfsh);  //NMD_7may2024

    first_unfixed_fish_time=1;
    fish_flags.initialize();
    diff_coffs.initialize();
    diff_coffs2.initialize();
    diff_coffs3.initialize();
    diff_coffs_prior.initialize();
    diff_coffs2_prior.initialize();
    diff_coffs3_prior.initialize();
    new_fm_level.initialize();
    N.initialize();
    testpar=0.0;
    effort_normalization_factor.initialize();
    effort.initialize();
    average_effort.initialize();
    really_true_effort.initialize();
    simulation_seeds.initialize();
    projected_simulated_data_flags.initialize();
    len_sample_size.initialize();
    wght_sample_size.initialize();
    lagrange_lambda.initialize();
    other_lagrange_lambda.initialize();
    pop_delta.initialize();
    epop_delta.initialize();
    int mmin=sel_seasons.indexmin();
    int mmax=sel_seasons.indexmax();
    for (int i=1;i<=nregions;i++)
    {
      for (int j=1;j<=nfp(i);j++)
      {
        num_fish(i,j)=1.e+200;
        fish_mort(i,j)=1.e+200;
      }
    }
    for (int i=mmin;i<=mmax;i++)
    {
      sel_seasons(i)=1;
    }

    if (pmsd)
    {
      //D(1,nregions,1,nregions),
      //Dad(1,_mo.num_periods(),1,ng,1,nregions,1,nregions),
      D.deallocate();
      Dad.deallocate();
      Dad2.deallocate();
      int nrr=pmsd->num_real_regions;
      int ns=pmsd->num_species;
      species_pars.allocate(1,20,1,ns);
      int itmp=nrr*ns;
      D.allocate(1,itmp,1,nrr);
      Dad.allocate(1,_mo.num_periods(),1,ng,1,itmp,1,nrr);
      Dad2.allocate(1,_mo.num_periods(),1,ng,1,itmp,1,nrr);
    }
    else
    {
      species_pars.allocate(1,20,1,1);
    }
    orthogonal_diffusion_matrix=
      get_orthogonal_diffusion_matrix(max(1,sum(_Dflags)));

    new_orthogonal_diffusion_matrix=get_n_season_diffusion_matrix();


    if (num_tag_releases)
    {
      tagmort.allocate(1,num_tag_releases);
      tag_fish_rep.allocate(1,ntg+1,1,nfsh);
      tag_fish_rep_group_flags.allocate(1,ntg+1,1,nfsh);
      tag_fish_rep_active_flags.allocate(1,ntg+1,1,nfsh);
      tag_fish_rep_target.allocate(1,ntg+1,1,nfsh);
      tag_fish_rep_penalty.allocate(1,ntg+1,1,nfsh);
    }
    int i;
    for (i=0;i<10;i++)
    {
      proj_output_files[i]=0;
    }
    for (i=1;i<=nregions;i++)
    {
      for (int j=1;j<=nfp(i);j++)
      {
        catch(i,j)=12345678;
      }
    }
    q_flag.initialize();

    effort_weight_by_fishery.initialize();
    effort_weight.initialize();
    //region_rec_diff_colmeans.initialize();
    bh_numbers.initialize();
    bh_reproductive_biomass.initialize();
    mean_weight.initialize();
    //cout << " we were here" << endl;
    num_regions=nregions;
    num_real_regions=nregions;
    nage=ng;
    nyears=nyrs;
    sel_dev_coffs.initialize();
    corr_wy.initialize();
    num_fish_periods=nfp;
    num_real_fish_periods=num_fish_periods;
    num_fish_incidents=nfi;
    num_fisheries=nfsh;
	ivector stuff;
    num_fish_times = nft;
    num_real_fish_times = nft;
    int nry=(nyears-1)/num_seasons+1;
    orth_recr_tot.allocate(0,nry-1);
    orth_recr_season.allocate(1,num_seasons-1,0,nry-1);   
    orth_recr_all.allocate(1,num_seasons*num_regions,0,nry-1);   
    yearly_recr_all.allocate(1,num_seasons*num_regions,0,nry);   
    orth_recr_region.allocate(1,num_regions-1,0,nry-1);  
    orth_recr_ses_reg.allocate(1,num_seasons-1,1,num_regions-1,0,nry-1);

    orth_recr_tot.initialize();
    orth_recr_season.initialize();
    orth_recr_all.initialize();
    yearly_recr_all.initialize();
    orth_recr_region.initialize();
    orth_recr_ses_reg.initialize();

    sseason.initialize();
    bblock.initialize();

  }

  void error_msg2(const char *s)
  {
    cerr << "Error reading " << s << endl;
  }


  MY_DOUBLE_TYPE age_at_length_calc(MY_DOUBLE_TYPE v,dvector& vb_coff,int nage)
  {
    MY_DOUBLE_TYPE rho=exp(-vb_coff(3));
    MY_DOUBLE_TYPE tmp= (v-vb_coff(1))/(vb_coff(2)-vb_coff(1));
    MY_DOUBLE_TYPE tmp1= 1.-tmp*(1-pow(rho,nage-1));
    if (tmp1<=0.0) tmp1=1.e-20;
    MY_DOUBLE_TYPE age= 1.-log(tmp1)/vb_coff(3);
    if (age<1) age=1;
    if (age>nage) age=nage;
    return age;
  }

  MY_DOUBLE_TYPE daves_kludge(MY_DOUBLE_TYPE x)
  {
    MY_DOUBLE_TYPE cx=x;
    int i=static_cast<int>(cx);
#if !defined(NO_MY_DOUBLE_TYPE)
    if (cx-i <= 0.5L)
#else
    if (cx-i <= 0.5)
#endif
    {
      MY_DOUBLE_TYPE tmp=x-i;
      MY_DOUBLE_TYPE tmp2=tmp*tmp;
      MY_DOUBLE_TYPE tmp3=tmp*tmp*tmp;
      return (24*tmp3-64*tmp3*tmp+48*tmp3*tmp2);
    }
    else
    {
      MY_DOUBLE_TYPE tmp=1-(x-i);
      MY_DOUBLE_TYPE tmp2=tmp*tmp;
      MY_DOUBLE_TYPE tmp3=tmp*tmp*tmp;
      return (1.-24*tmp3+64*tmp3*tmp-48*tmp3*tmp2);
	 }
  }

  fishery_freq_record_array::fishery_freq_record_array(int min,
    int max,int _nlint,MY_DOUBLE_TYPE _shlen,MY_DOUBLE_TYPE _filen,int _nwint,
    MY_DOUBLE_TYPE _wshlen,MY_DOUBLE_TYPE _wfilen)
  {
    nlfreq=max-min+1;
    index_min=min;
    index_max=max;
    nlint=_nlint;
    shlen=_shlen;
    filen=_filen;
    nwint=_nwint;
    wshlen=_wshlen;
    wfilen=_wfilen;
    int sz=size();
    ptr=new fishery_freq_record [sz];
    if (sex_flag)
    {
      ptr1=new fishery_freq_record [sz];
      if (ptr==NULL)
      {
        cerr << "Error allocating memory in fishery_freq_record_array"<<endl;
	  exit(21);
      }
      ptr1-=indexmin();
    }
    else
      ptr1=NULL;  

    if (ptr==NULL)
    {
      cerr << "Error allocating memory in fishery_freq_record_array"<<endl;
	exit(21);
    }
    ptr-=indexmin();
    allocate();
  }

  cifstream& operator >> (cifstream& cif, dvar_fish_stock_history& fsh)
  {
    dvariable xtmp=0.0;
    fsh.old_fish_flags=fsh.fish_flags;

    fsh.sum_ff71_flag=sum(column(fsh.fish_flags,71));
    fsh.sum_ff48_flag=sum(column(fsh.fish_flags,48));
    fsh.sum_ff79_flag=sum(column(fsh.fish_flags,79));
    ivector  ff74=column(fsh.fish_flags,74);
    
    check_number(cif,872,fsh);
    if (!fsh.parest_flags(197) && _file_version_no > 1016 
      && fsh.num_tag_releases)
    {
      // ********************************
      // begin new stuff
      // ********************************
      if (fsh.pmsd)
      {
        cif >> fsh.true_tag_flags;
        if (!cif) error_msg2("true_tag_flags");
        fsh.tag_flags.deallocate();
//        fsh.tag_flags.allocate(1,fsh.num_tag_releases,1,10);
        fsh.tag_flags.allocate(1,fsh.num_tag_releases); //imatrix, unall.vectors
        ivector iind(1,fsh.pmsd->num_species);
        iind.initialize();
        for (int it=1;it<=fsh.old_num_tag_releases;it++)
        {
          for (int is=1;is<=fsh.pmsd->num_species;is++)
          {
            if (fsh.pmsd->tag_species_flag(it,is))
            {
              int ioff=fsh.offset(is)+iind(is)+1;
              fsh.tag_flags(ioff)=fsh.true_tag_flags(it);
              iind(is)++;
            }
          }
        }
      }
      else
      {
        cif >> fsh.tag_flags;
        if (!cif) error_msg2("tag_flags");
        fsh.true_tag_flags=fsh.tag_flags;
      }
      // ********************************
      // end new stuff
      // ********************************
    } 
    else 
    {
      if (fsh.num_tag_releases)
      {
        if (_ini_version_no>6)
        {
          cif >> fsh.tag_flags;
          if (!cif) error_msg2("tag_flags");
          fsh.true_tag_flags=fsh.tag_flags;
        }
        else
        {
          fsh.tag_flags.initialize();
        }
      }
//         fsh.tag_flags.initialize();
    }
    if (!fsh.parest_flags(197) && _file_version_no > 1064 && 
      fsh.num_tag_releases >0)
    {
      cif >> fsh.tagmort;   //RRRR
    }
    else
    {
      if (_ini_version_no>3 && fsh.num_tag_releases >0) //NMD_9jun2025)
      {
        cif >> fsh.tagmort;   //RRRR
      }
      else
      {
        if (allocated(fsh.tagmort))
        {
          fsh.tagmort.initialize();   //RRRR
        }
      }
    }
//    check_number(cif,1872,fsh.parest_flags);
    if (!fsh.parest_flags(197) && _file_version_no > 1041 
      && fsh.num_tag_releases)
    {
      check_number(cif,1872,fsh);  //NMD24jan2018
      cif >> fsh.tag_fish_rep;
      cif >> fsh.tag_fish_rep_group_flags;
      cif >> fsh.tag_fish_rep_active_flags;
      cif >> fsh.tag_fish_rep_target; 
      cif >> fsh.tag_fish_rep_penalty; 
    }
    else
    {
      if (_ini_version_no==0)
      {
        fsh.tag_fish_rep.initialize();
        fsh.tag_fish_rep_group_flags.initialize();
        fsh.tag_fish_rep_active_flags.initialize();
        fsh.tag_fish_rep_target.initialize(); 
        fsh.tag_fish_rep_penalty.initialize(); 
      }
      else if(fsh.num_tag_releases >0) //NMD_9jun2025
      {
        cif >> fsh.tag_fish_rep;
        if (!cif) error_msg2("fsh.tag_fish_rep");
        cif >> fsh.tag_fish_rep_group_flags;
        if (!cif) error_msg2("fsh.tag_fish_rep_group_flags");
        cif >> fsh.tag_fish_rep_active_flags;
        if (!cif) error_msg2("fsh.tag_fish_rep_active_flags");
        cif >> fsh.tag_fish_rep_target;
        if (!cif) error_msg2("fsh.tag_fish_rep_target");
        cif >> fsh.tag_fish_rep_penalty;
        if (!cif) error_msg2("fsh.tag_fish_rep_penalty");
      }
    }
    check_number(cif,2872,fsh);
      
    test_the_pointer();
    {
      int mmin=fsh.fish_flags.rowmin();
      int mmax=fsh.fish_flags.rowmax();
      for (int i=mmin;i<=mmax;i++)
      {
        if(fsh.fish_flags(i,3)==0)
          fsh.fish_flags(i,3)=fsh.nage_by_fishery(i)-1;
        /*
        if(fsh.fish_flags(i,21)>0)
        {
          //fsh.fish_flags(i,3)=fsh.nage-fsh.fish_flags(i,21)+1;
        }
        if(fsh.fish_flags(i,21)==0)
        {
          //fsh.fish_flags(i,3)=fsh.nage-1;
        }
        */
      } 
    }

    if (!fsh.parest_flags(197) )
    {
      if ( _file_version_no > 1011)
      {
        cif >> fsh.region_flags;
        if (!cif) error_msg2("region_flags");
      }
      else
      {
        fsh.region_flags.initialize();
      }
    }
    check_number(cif,873,fsh);
    if (fsh.parest_flags(197) )
    {
      if (_ini_version_no>1)
      {
        cif >> fsh.region_flags;
        if (!cif) error_msg2("region_flags");
      }
      else
      {
        fsh.region_flags.initialize();
      }
    }
    if (fsh.pmsd)
    {
      if (!fsh.parest_flags(197) )
      {
        if ( _file_version_no > 1043)
        {
          cif >> fsh.pmsd->species_flags;
          if (!cif) error_msg2("species_flags");
        }
        else
        {
          fsh.pmsd->species_flags.initialize();
        }
      }
      if (fsh.parest_flags(197))
      {
        if (_ini_version_no>1)
        {
          cif >> fsh.pmsd->species_flags;
          if (!cif) error_msg2("species_flags");
        }
        else
        {
          fsh.pmsd->species_flags.initialize();
        }
      }  
    }
    check_number(cif,1874,fsh);

    if (fsh.parest_flags(197) || _file_version_no > 1017)
    {
      cif >> fsh.pmature;
      for (int i=fsh.pmature.indexmin();i<=fsh.pmature.indexmax();i++)
      {
        if (fsh.pmature(i)<=0.0 || fsh.pmature(i)>1.0)
        {
          cout << "Bounds error reading pmature(" << i << ") in par file "
               << endl << " value is " << fsh.pmature(i)<< endl;
          ad_exit(1);
        }
      }
      if (!cif) error_msg2("pmature");
    }
    else
    {
      fsh.pmature.initialize();
    }

    check_number(cif,874,fsh);


    if (fsh.pmsd)
    {
      if (fsh.parest_flags(197) || _file_version_no > 1041)
      {
        cif >> fsh.pmsd->pmature;
        if (!cif) error_msg2("pmsd->pmature");
        for (int is=2;is<=fsh.pmsd->num_species;is++)
        {
          for (int i=fsh.pmsd->pmature(is).indexmin();i<=fsh.pmsd->pmature(is).indexmax();i++)
          {
            if (fsh.pmsd->pmature(is,i)<=0.0 || fsh.pmsd->pmature(is,i)>1.0)
            {
              cout << "Bounds error reading multi-species pmature(" << i 
               << ") in par file "
               << endl << " value is " << fsh.pmsd->pmature(is,i)<< endl;
              ad_exit(1);
            }
          }
        }
      }
    }
    check_number(cif,875,fsh);

    if (!fsh.parest_flags(197))
    {
      cif >> fsh.totpop_coff;
      if (fsh.pmsd)
      {
        cif >> fsh.pmsd->totpop_coff;
      }
      check_number(cif,876,fsh);
      if (!cif) error_msg2("totpop");
      if (_file_version_no > 1032)
      {
        cif >> fsh.implicit_totpop_coff;
        if (!cif) error_msg2("implicit_totpop");
      }
      else
      {
        fsh.implicit_totpop_coff=15;
      }
      if (_file_version_no > 1024)
      {
        cif >> fsh.rec_init_diff;
        if (fsh.pmsd)
        {
          cif >> fsh.pmsd->rec_init_diff;
        }
        if (!cif) error_msg2("rec_init_diff");
      }
      else
      {
        fsh.rec_init_diff=0;
        if (fsh.pmsd)
        {
          fsh.pmsd->rec_init_diff.initialize();
        }
        if (!cif) error_msg2("rec_init_diff");
      }
#ifndef __GNUC__
//    CLOGF_TRACE(fsh.totpop)
#endif
      check_number(cif,877,fsh);
      if (_file_version_no > 1013)
      {
        cif >> fsh.rec_times;
        if (!cif) error_msg2("rec_times");
      }
      else
      {
        fsh.rec_times.initialize();
      }

      if (_file_version_no > 1050) //NMD_19May2016
      {
        cif >> fsh.recr;
      }
      else
      {
        MY_DOUBLE_TYPE tmprcr;
        cif >> tmprcr;
        cif >> fsh.recr;
      }
//      cif >> fsh.recr; 
      if (!cif) error_msg2("recr");


      if (fsh.pmsd)
      {
        if (_file_version_no > 1050) //NMD_19May2016
        {
          cif >> fsh.pmsd->recr;
        }
        else
        {
          int ns=fsh.pmsd->num_species;
          for (int is=2;is<=ns;is++)
          {

            MY_DOUBLE_TYPE tmprcr;
            cif >> tmprcr;
            cif >> fsh.pmsd->recr(is);
          }
        }
//        cif >> fsh.pmsd->recr;    //NMD_19May2016
      }
    }
    else
    {
      fsh.totpop_coff=25.0;   //NMD_2jun2020
      if (_ini_version_no > 5) //NMD_9jun2025
      {
        double tmp_totpop;
        cif >> tmp_totpop;
        if (tmp_totpop > 0) fsh.totpop_coff=tmp_totpop;
      }
      fsh.rec_times.initialize();
      fsh.recr=1;
      if (fsh.pmsd)
      {
//        fsh.pmsd->totpop_coff=15.0;
        fsh.pmsd->totpop_coff=25.0;    //NMD_2jun2020
        if (_ini_version_no > 5) //NMD_9jun2025
        {
          MY_DOUBLE_TYPE tmp_totpop;
          int ns=fsh.pmsd->num_species;
          for (int is=2;is<=ns;is++)
          {
            cif >> tmp_totpop;
            if (tmp_totpop > 0) fsh.pmsd->totpop_coff(is)=tmp_totpop;
          }
        }
        fsh.pmsd->recr=1.0;
      }
    } 
    check_number(cif,878,fsh);
    if (!fsh.parest_flags(197) && _file_version_no > 1046)
    {
      cif >> fsh.lagrange_lambda;
    }
    else
    {
      fsh.lagrange_lambda.initialize();
    }
    if (!fsh.parest_flags(197) && _file_version_no > 1052)
    {
      cif >> fsh.other_lagrange_lambda;
    }
    else
    {
      fsh.other_lagrange_lambda.initialize();
    }
    if (!fsh.parest_flags(197) && _file_version_no > 1015)
    {
      cif >> fsh.rep_dev_coffs;
      if (!cif) error_msg2("rep_dev_coffs");
    }
    else
    {
      fsh.rep_dev_coffs.initialize();
    }
    if (!fsh.parest_flags(197))
    {
      if (_file_version_no > 999)
      {
        cif >> fsh.avail_coff;
        if (!cif) error_msg2("avail_coff");
      }
      else
      {
        fsh.avail_coff.initialize();
      }
    }
    else
    {
      fsh.avail_coff.initialize();
      fsh.orth_recr.initialize();
    }
    check_number(cif,879,fsh);
    if (_file_version_no > 1036)
    {
      if (fsh.historical_parest_flags(155)>0)
      {
        cif >> fsh.yearly_recr_all;   
        check_number(cif,6875,fsh);
        cif >> fsh.orth_recr_all;   
        check_number(cif,5875,fsh);
        if (!cif) error_msg2("year_recr_all");
      }
      else
      {
        fsh.yearly_recr_all.initialize();   
        fsh.orth_recr_all.initialize();   
        int mmin=fsh.orth_recr_all.indexmin();   
        int mmax=fsh.orth_recr_all.indexmax();   
        for (int i=mmin;i<=mmax;i++)
        {
          fsh.orth_recr_all(1,0)=200.0;   
        }
      }
    }
    else
    {
      fsh.yearly_recr_all.initialize();   
    }
    if (_file_version_no > 1054) 
    {
      if (fsh.pmsd && fsh.pmsd->num_species>1)
      {
        for (int is=2;is<=fsh.pmsd->num_species;is++)
        {
          if (fsh.pmsd->historical_parest_flags(is,155)>0)
          {
            cif >> fsh.pmsd->yearly_recr_all(is);
            cif >> fsh.pmsd->orth_recr_all(is);
          }
        }
      }
      if (!cif) error_msg2("orth_recr_all");
    }

/* //NMD15dec2017
    if (_file_version_no > 1034)
    {
      if (fsh.parest_flags(156)>0)
      {
        cif >> fsh.orth_recr_all;   
        if (!cif) error_msg2("orth_recr_all");
      }
    }
    else 
    {
      fsh.orth_recr_all.initialize();   
    }
*/
    ivector size_stuff=fsh.calculate_blocked_sel_sizes();
    check_number(cif,4875,fsh);
    if (fsh.parest_flags(197)==0)
    {
      cif >> fsh.initpop;
      check_number(cif,3875,fsh);
      if (!cif)
      {
        cerr << "Error reading PAR file"<<endl;
        ad_exit(1);
      }
      if (!cif) error_msg2("initpop");
      ivector ff71=column(fsh.fish_flags,71);
      ivector hff71;
      ivector hff74;
      int tmp71=sum(column(fsh.historical_fish_flags,71));
      int tmp74=sum(column(fsh.historical_fish_flags,74));
      if (_file_version_no < 1048 || (!tmp71 && !tmp74))  //NMD_25Mar2015
      {
        hff71=column(fsh.fish_flags,71);
        hff74=column(fsh.fish_flags,74);
      }
      else
      {
        hff71=column(fsh.historical_fish_flags,71);
        hff74=column(fsh.historical_fish_flags,74);
      }  //NMD_28jan2015
  
  
      int changed_size_flag=0;
      int i;
      for (i=1;i<=fsh.num_fisheries;i++)
      {
        if (ff74(i) <hff74(i) || ff71(i) <hff71(i))
        {
          cerr << "can not decrease the size of ff71 or ff74" << endl;
          ad_exit(1);
        }
        else if (ff74(i) >hff74(i) || ff71(i) >hff71(i))
        {
          changed_size_flag=1;
        }
      }
      fsh.bs_selmean.allocate(1,fsh.num_fisheries,
        1,ff74,1,ff71+1);
      fsh.bs_selcoff.allocate(1,fsh.num_fisheries,
        1,ff74,1,ff71+1,1,size_stuff);
      if (changed_size_flag==0)
      {
        //if (_file_version_no > 1047 )
        {
          read_bs_selcoff(cif,fsh.bs_selcoff,fsh.nage);
        }      
        //else
       // {
         // fsh.bs_selcoff.initialize();
        //}
      }
      else
      {
        d4_array tbs(1,fsh.num_fisheries,1,hff74,1,hff71+1,1,size_stuff);
        if (_file_version_no > 1047 )
        {
          read_bs_selcoff(cif,tbs,fsh.nage);
        }      
        else
        {
          tbs.initialize();
        }
        if (!cif)
        {
          cerr << "Error reading PAR file"<<endl;
          ad_exit(1);
        }
        
        for (i=1;i<=fsh.num_fisheries;i++)
        {
          if (ff71(i)>hff71(i))
          {
            int rem1=(ff71(i)+1)%(hff71(i)+1);
            if (rem1)
            {
              cerr << "must increase number of time blocks by a multiple"
                   << endl;
            }
          }
          if (ff74(i)>hff74(i))
          {
            int rem2=(ff74(i))%(hff74(i));
            if (rem2)
            {
              cerr << "must increase number of seasons by a multiple"
                   << endl;
            }
          }
          for (int ii=1;ii<=ff74(i);ii++)
          {
            for (int jj=1;jj<=ff71(i)+1;jj++)
            {
              int div1=ff74(i)/hff74(i);
              int div2=(1+ff71(i))/(1+hff71(i));
              int off1=(ii-1)/div1+1;
              int off2=(jj-1)/div2+1;
              fsh.bs_selcoff(i,ii,jj)=tbs(i,off1,off2);
            }
          }
        }
      }
    }
    else
    {
      fsh.bs_selcoff.allocate(1,fsh.num_fisheries,
        1,1,1,1,1,size_stuff);
      fsh.bs_selcoff.initialize();
    }

    if (!cif)
    {
      cerr << "Error reading PAR file"<<endl;
      ad_exit(1);
    }
    if (!fsh.parest_flags(197))
    {
      if (_file_version_no < 1030)
      {
        fsh.bs_selcoff=log(fsh.bs_selcoff);
      }
      if (_file_version_no > 1040)
      {
        cif >> fsh.ageselcoff;
      }
      else
      {
        fsh.ageselcoff.initialize();
        fsh.bs_selcoff.initialize();
      }
    }
    cif >> fsh.nat_mort_coff;
    // check to see if number of seasons or number of time blocks
    // has changed
    if (!cif)
    {
      cerr << "Error reading PAR file"<<endl;
      ad_exit(1);
    }
    if (!cif) error_msg2("nat_mort_coff");
	//NMD 02Nov2011
    if (fsh.pmsd)
    {
      cif >> fsh.pmsd->nat_mort_coff;
      if (!cif) error_msg2("pmsd->nat_mort_coff");
	}
//NMD 02Nov2011
    if (!fsh.parest_flags(197))
    {
      cif >> fsh.q0;
      if (!cif) error_msg2("q0");
      cif >> fsh.q1;
      if (!cif) error_msg2("q1");
    }
    else
    {
      fsh.initpop=1;
      fsh.fishery_sel=1.0;
      //fsh.nat_mort_coff=0.2;
      fsh.q0=.1;
      fsh.q1.initialize();
    }

    if (_file_version_no > 1037)
    {
      if (!fsh.parest_flags(197))
      {
        cif >> fsh.q0_miss;
      }
      else
      {
        fsh.q0_miss=-5.0;
      }
    }
    else
    {
      fsh.q0_miss=-5.0;
    }
    check_number(cif,3876,fsh);
    if (fsh.missing_catch_flag)
    {
      if (_file_version_no > 1038 && fsh.parest_flags(197)==0)
      {
        cif >> fsh.fm_level_devs;
      }
      else
      {
        fsh.fm_level_devs.initialize();
      }
    }
    if (!cif)
    {
      cerr << "Error reading PAR file"<<endl;
      ad_exit(1);
    }

    if (allocated(fsh.move_map)) fsh.move_map.deallocate();
    fsh.move_map.allocate(1,fsh.mo.num_periods());
    if (_file_version_no > 1025)
    {
      cif >> fsh.move_map;
      if (!cif) error_msg2("move_map");
      if (_file_version_no < 1057)
      {
        //fsh.parest_flags(356)=-1;
        fsh.parest_flags(357)=1;
//        fish_nonmv_flag=1;   //NMD_14Sep2018
      }
    }
    else
    {
      fsh.move_map=1;
    }
        
    if (!cif)
    {
      cerr << "Error reading PAR file"<<endl;
      ad_exit(1);
    }
    if (_file_version_no > 1009)
    {
      cif >> fsh.diff_coffs;
      if (!cif) error_msg2("diff_coffs");
    }
    else
    {
      fsh.diff_coffs.initialize();
    }
    if (!fsh.parest_flags(197) && _file_version_no > 1058)
    {
      cif >> fsh.xdiff_coffs;
      if (!cif) error_msg2("xdiff_coffs");
    }
    else
    {
      fsh.xdiff_coffs.initialize();
    }
    if (!fsh.parest_flags(197) && _file_version_no > 1061)
    {
      cif >> fsh.y1diff_coffs;
      if (!cif) error_msg2("y1diff_coffs");
      cif >> fsh.y2diff_coffs;
      if (!cif) error_msg2("y2diff_coffs");
    }
    else
    {
      fsh.y1diff_coffs.initialize();
      fsh.y2diff_coffs.initialize();
    }
    if (!fsh.parest_flags(197) && _file_version_no > 1062)
    {
      cif >> fsh.zdiff_coffs;
      if (!cif) error_msg2("xdiff_coffs");
    }
    else
    {
      fsh.zdiff_coffs.initialize();
    }
    if (!fsh.parest_flags(197) && _file_version_no > 1027)
    {
      cif >> fsh.Dad;
      if (!cif) error_msg2("Dad");
    }

    if (!fsh.parest_flags(197))
    {
      if (_file_version_no > 1012)
      {
        cif >> fsh.diff_coffs2;
        if (!cif) error_msg2("diff_coffs2");
        cif >> fsh.diff_coffs3;
        if (!cif) error_msg2("diff_coffs3");
      }
      else
      {
        fsh.diff_coffs2.initialize();
        fsh.diff_coffs3.initialize();
      }
    }
    else
    {
      //fsh.diff_coffs.initialize();
      fsh.diff_coffs2.initialize();
      fsh.diff_coffs3.initialize();
    }
    if (!fsh.parest_flags(197))
    {
      if (_file_version_no > 1031)
      {
        cif >> fsh.diff_coffs_prior;
        cif >> fsh.diff_coffs2_prior;
        cif >> fsh.diff_coffs3_prior;
      }
      else
      {
        fsh.diff_coffs_prior.initialize();
        fsh.diff_coffs2_prior.initialize();
        fsh.diff_coffs3_prior.initialize();
      }
    }
    else
    {
      fsh.diff_coffs_prior=fsh.diff_coffs;
      fsh.diff_coffs2_prior.initialize();
      fsh.diff_coffs3_prior.initialize();
    }  
    if (_file_version_no > 1025)
    {
      if (_file_version_no < 1057)
      {
        fsh.modify_Dad();
        fsh.modify_move_map();
        fsh.modify_all_diff_coffs();
      }
    }
    if (!fsh.parest_flags(197))
    {
      if (_file_version_no > 1010)
      {
        cif >> fsh.region_rec_diff_coffs;
        if (!cif) error_msg2("region_rec_diff_coffs");
      }
      else
      {
        fsh.region_rec_diff_coffs.initialize();
      }
    }
    else
    {
      fsh.region_rec_diff_coffs.initialize();
    }


    if (!fsh.parest_flags(197))
    {
      if (!fsh.parest_flags(120))
      {
        cif >> fsh.effort_dev_coffs;
      }
      if (!cif) error_msg2("effort_dev_coffs");
      cif >> fsh.corr_wy;
    }
    else
    {
      fsh.effort_dev_coffs.initialize();
      fsh.corr_wy.initialize();
    }


    if (!fsh.parest_flags(197) && _file_version_no > 1006)
    {
      if ( _file_version_no > 1049)
      {
        cif >> fsh.fish_pars;
      }
      else if ( _file_version_no > 1026)
      {
        cif >> fsh.fish_pars.sub(1,20);
        fsh.fish_pars.sub(21,50).initialize();
      }
      else
      {
        cif >> fsh.fish_pars.sub(1,10);
        fsh.fish_pars.sub(11,50).initialize();
      }
      if (!cif) error_msg2("fish_pars");
      if (!sum(value(fsh.fish_pars(3))))
      {
        fsh.fish_pars(3)=1;
      }
    }
    else
    {
      fsh.fish_pars.initialize();
    }

    if (!fsh.parest_flags(197) && _file_version_no > 1060)
    {
      cif >> fsh.implicit_fm_level_regression_pars;
    }
    else
    {
      fsh.implicit_fm_level_regression_pars.initialize();
    }

    if (!fsh.parest_flags(197) && _file_version_no > 1051)
    {
      cif >> fsh.species_pars;
    }
    else
    {
      fsh.species_pars.initialize();
    }

    if (!fsh.parest_flags(197) && _file_version_no > 1022)
    {
      cif >> fsh.seasonal_catchability_pars;
      if (!cif) error_msg2("fsh.seasonal_catchability_pars");
    }
    else
    {
      fsh.seasonal_catchability_pars.initialize();
    }



    if (_file_version_no > 1008 )
    {
      fsh.age_pars.initialize();
      if (_file_version_no > 1023 )
      {
        cif >> fsh.age_pars;
//NMD 04Nov2011
        if (fsh.pmsd)
        {
          if (_file_version_no > 1041)
          {
            for (int is=2;is<=fsh.pmsd->num_species;is++)
            {
              cif >> fsh.pmsd->age_pars(is);
              if (!cif) error_msg2("pmsd->age_pars");
            }
          }
        }
//NMD 04Nov2011
      }
      else
      {
        for (int i=1;i<=5;i++) cif >> fsh.age_pars(i);
      }
      if (!cif) error_msg2("age_pars");
    }
    else
    {
      fsh.age_pars.initialize();
    }

    if (fsh.parest_flags(121)>0)
    {
      double sum_aprs5=sum(value(fsh.age_pars(5)));
//      if (fsh.historical_parest_flags(121)==0)
      if (fsh.historical_parest_flags(121)==0 && sum_aprs5==0.0)  //NMD_28apr2023
      {
        int num=fsh.parest_flags(121);
        cout << "ZZZZZ  " << fsh.nat_mort_coff << endl;
        //ad_exit(1);
        fsh.age_pars(5)(1,num)=log(fsh.nat_mort_coff);
      }
    }
    if (fsh.pmsd)   //NMD_7jul2020
    {
      for (int is=2;is<=fsh.pmsd->num_species;is++)
      {
        if (fsh.pmsd->parest_flags(is,121)>0)
        {
          if (fsh.pmsd->historical_parest_flags(is,121)==0)
          {
            int num=fsh.pmsd->parest_flags(is,121);
            cout << "ZZZZZ  " << fsh.pmsd->nat_mort_coff(is) << endl;
            //ad_exit(1);
            fsh.pmsd->age_pars(is,5)(1,num)=log(fsh.pmsd->nat_mort_coff(is));
          }
        }
      }
    }   //NMD_7jul2020

     test_the_pointer();
    if (!fsh.parest_flags(197) && _file_version_no > 1007)
    {
      if (!fsh.parest_flags(197) && _file_version_no > 1048)
      {
        cif >> fsh.region_pars;
      }
      else
      {
        fsh.region_pars.initialize();
        cif >> fsh.region_pars.sub(1,10);
      }
      if (!cif) error_msg2("region_pars");
    }
    else
    {
      fsh.region_pars.initialize();
      cif >> fsh.region_pars(1);
    }

    if (!fsh.parest_flags(197) && _file_version_no > 1007
       && _file_version_no < 1012)
    {
      cif >> fsh.region_flags;
      if (!cif) error_msg2("region_flags");
    }

    if (!fsh.parest_flags(197))
    {
      if (!fsh.parest_flags(120))
      {
        cif >> fsh.catch_dev_coffs;
        if (!cif) error_msg2("effort_dev_coffs");
      }
      else
      {
        int i;
        for (i=1;i<=fsh.num_fisheries;i++)
        {          
          for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
          {
            int rr=fsh.realization_region(i,nt);
            cif >>  fsh.catchability(rr,fsh.realization_period(i,nt),
                       fsh.realization_incident(i,nt));
          }
        }
        for (i=1;i<=fsh.num_fisheries;i++)
        {
          for (int nt=1;nt<fsh.num_fish_times(i);nt++)
          {
            int rr1=fsh.realization_region(i,nt+1);
            int rr=fsh.realization_region(i,nt);
            fsh.catch_dev_coffs(i,nt+1)=
              fsh.catchability(rr1,fsh.realization_period(i,nt+1),
                     fsh.realization_incident(i,nt+1))
             -fsh.catchability(rr,fsh.realization_period(i,nt),
                     fsh.realization_incident(i,nt));
          }
        }
      }
    }
    else
    {
      fsh.catchability.initialize();
      fsh.catch_dev_coffs.initialize();
    }
    
    ivector tmp(1,fsh.num_fisheries);
    int i;
    for (i=1;i<=fsh.num_fisheries;i++)
    {
      tmp(i)=max(1,fsh.fish_flags(i,3));
    }

    if (!(!fsh.delta2))
    {
      fsh.delta2.deallocate();
    }
    fsh.delta2.allocate(1,fsh.num_fisheries,2,fsh.num_fish_times,1,tmp);
    fsh.delta2.initialize();

    xtmp=0.0;
    if (!fsh.parest_flags(197))
    {
      if (!fsh.parest_flags(120))
      {
        ptarr1= new d3_array(1,fsh.num_fisheries,2,fsh.num_fish_times,1,fsh.nage);
        d3_array& tarr=*ptarr1;
        cif >> tarr;
    xtmp=0.0;
        if (!cif) error_msg2("tarr");
        // allocate the selectivity deviations -- delta2
        for (i=1;i<=fsh.num_fisheries;i++)
        {
          for (int nt=2;nt<=fsh.num_fish_times(i);nt++)
          {
            for (int j=1;j<=fsh.fish_flags(i,3);j++)
            {
              fsh.delta2(i,nt,j) = tarr(i,nt,j);
            }
          }
        }
        xtmp=0.0;
      }
      else
      {
        for (i=1;i<=fsh.num_fisheries;i++)
        {
          for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
          {
            int rr=fsh.realization_region(i,nt);
            cif >> fsh.incident_sel(rr,fsh.realization_region(i,nt),
                  fsh.realization_period(i,nt),fsh.realization_incident(i,nt));
          }
        }
        int ir;
        for (ir=1;ir<=fsh.num_regions;ir++)
        {
          fsh.incident_sel(ir)=log(fsh.incident_sel(ir));
        }
        for (i=1;i<=fsh.num_fisheries;i++)
        {
          for (int nt=2;nt<=fsh.num_fish_times(i);nt++)
          {
            for (int j=1;j<=fsh.fish_flags(i,3);j++)
            {
              int rr=fsh.realization_region(i,nt);
              int rr1=fsh.realization_region(i,nt-1);
              fsh.delta2(i,nt,j) =
                fsh.incident_sel(rr,fsh.realization_period(i,nt),
                fsh.realization_incident(i,nt),j)
                -fsh.incident_sel(rr1,fsh.realization_period(i,nt-1),
                   fsh.realization_incident(i,nt-1),j);
            }
          }
        }
        for (ir=1;ir<=fsh.num_regions;ir++)
        {
          fsh.incident_sel(ir)=exp(fsh.incident_sel(ir));
        }
      }
    }
    else
    {
      for (int ir=1;ir<=fsh.num_regions;ir++)
      {
        int mmin=fsh.incident_sel(ir).indexmin();
        int mmax=fsh.incident_sel(ir).indexmax();
        for (int j=mmin;j<=mmax;j++)
        {
          if (allocated(fsh.incident_sel(ir,j)))
            fsh.incident_sel(ir,j)=1.0;
        }
      }
    }
    delete ptarr1;
    ptarr1=0;
    xtmp=0.0;
    if (!fsh.parest_flags(197))
    {
      if (_file_version_no > 1000)
      {
        for (int i=1;i<=fsh.num_fisheries;i++)
        {
          for (int j=1;j<=fsh.num_fish_times(i);j++)
          { 
            cif >> fsh.sel_dev_coffs(i,j);
          }
        }
      }
      if (!cif) fsh.sel_dev_coffs.initialize();
    }
    else
    {
      fsh.sel_dev_coffs.initialize();
    }

  /*
    ivector ff71=column(fsh.fish_flags,71);
    int sumff71=sum(ff71);
    ivector hff71=column(fsh.historical_fish_flags,71);
    int sumhff71=sum(hff71);
    if (_file_version_no > 1047 )
    {
      ivector size_stuff=fsh.calculate_blocked_sel_sizes();
      bs_selcoff.allocate(1,fsh.num_fisheries,1,ff74,1,ff71+1,1,size_stuff);
    }      
   */
   if (_file_version_no < 1048)
   {
     ivector ff57=column(fsh.fish_flags,57);
     ivector ff3=column(fsh.fish_flags,3);
     
     // move old parameters to new place for current model
     for (int i=1;i<=fsh.num_fisheries;i++)
     {
       if (ff57(i)==1)   // logistic
       {
         MY_DOUBLE_TYPE fp10=0.0;
         MY_DOUBLE_TYPE fp9=0.0;
         MY_DOUBLE_TYPE gp10=0.0;
         MY_DOUBLE_TYPE gp9=0.0;
         if (ff3(i)>0)
         {
           MY_DOUBLE_TYPE x=(1.0+fsh.nage)/2.0;
           MY_DOUBLE_TYPE y=(1.0+ff3(i))/2.0;
           MY_DOUBLE_TYPE x1=1.0-x;
           MY_DOUBLE_TYPE y1=1.0-y;
           fp9=value(fsh.fish_pars(9,i));
           fp10=1.0+value(fsh.fish_pars(10,i));
           gp10=x1/y1*fp10;
           gp9=(x-y)/y1+x1/y1*fp9;
         }
         else
         {
           gp10=value(1.0+fsh.fish_pars(10,i));
           gp9=value(fsh.fish_pars(9,i));
         }
         
         fsh.bs_selcoff(i,1,1,1)=gp9*fsh.bl_sel_scaling;
#if !defined(NO_MY_DOUBLE_TYPE)
         fsh.bs_selcoff(i,1,1,2)=(gp10-1.0L)*fsh.bl_sel_scaling;
#else
         fsh.bs_selcoff(i,1,1,2)=(gp10-1.0)*fsh.bl_sel_scaling;
#endif
       }
      
       if (ff57(i)==2)   // doufsh.ble normal
       {
         fsh.bs_selcoff(i,1,1,1)=fsh.fish_pars(9,i)*fsh.bl_sel_scaling;
         fsh.bs_selcoff(i,1,1,2)=fsh.fish_pars(10,i)*fsh.bl_sel_scaling;
         fsh.bs_selcoff(i,1,1,3)=fsh.fish_pars(11,i)*fsh.bl_sel_scaling;
       }
     }
   }

    xtmp=0.0;
    if (_file_version_no > 1033 && fsh.parest_flags(197)==0 )
    {
       
       cif >> fsh.year_flags;
       cif >> fsh.season_flags;
    }
    else
    {
      fsh.year_flags.initialize();
      fsh.season_flags.initialize();
    }
     test_the_pointer();
    return cif;
  }


  void dvar_fish_stock_history::natural_mortality_init()
  {
    nat_mort=nat_mort_coff;
//NMD 03Nov2011
    if (pmsd)
    {
      for (int is=2;is<=pmsd->num_species;is++)
      {
        for (int iy=1;iy<=nyears;iy++)
        {
          int tmp = pmsd->nage(is);
          for (int j=1;j<=tmp;j++)
          {
            MY_DOUBLE_TYPE tmp2 = value(pmsd->nat_mort_coff(is));
            pmsd->nat_mort(is,iy,j) = tmp2;
	  }
        }
      }
    }
//NMD 03Nov2011
    if (age_flags(48)>0)
    {
      for (int iy=1;iy<=nyears;iy++)  // Loop over fishing periods
      {
        MY_DOUBLE_TYPE mm=log(age_flags(49)/10.0);
        for (int j=nage-age_flags(48)+1;j<=nage;j++)
        {
          nat_mort(iy,j)+= mm;
        }
      }
    }
    if (pmsd)
    {
      for (int is=2;is<=pmsd->num_species;is++)
      {
        int ng=pmsd->nage(is);
        if (pmsd->age_flags(is,48)>0)
        {
          for (int iy=1;iy<=nyears;iy++)  // Loop over fishing periods
          {
            MY_DOUBLE_TYPE mm=log(pmsd->age_flags(is,49)/10.0);
            for (int j=ng-pmsd->age_flags(is,48)+1;j<=ng;j++)
            {
              pmsd->nat_mort(is,iy,j)+= mm;
            }
          }
        }
      }
    }
  }

 void dvar_fish_stock_history::get_initial_parameter_values()
 // Initialize all fundamental model parameters
 {
   //fishery_selcoff_init();
   natural_mortality_init();
 }

  void dvar_fish_stock_history::transform_in()
  {
    //fishery_sel=log(1.e-10+fishery_sel);
    {
      int mmin=effort.indexmin();
      int mmax=effort.indexmax();
      for (int i=mmin;i<=mmax;i++)
      {
        int jmin=effort(i).indexmin();
        int jmax=effort(i).indexmax();
        for (int j=jmin;j<=jmax;j++)
        {
          if (allocated(effort(i,j)))
            effort(i,j)=log(1.e-10+effort(i,j));
        }
      }
    }
        
    recr=log(1.e-30+recr);
    if (parest_flags(400)>0)
    {
      recr(last_real_year-parest_flags(400)+1,last_real_year).initialize();
    }
    q0=log(1.e-20+q0);
    initpop=log(1.e-30+initpop);
    nat_mort_coff=log(1.e-20+nat_mort_coff);
//NMD 02Nov2011
    if (pmsd){   // XXXYYY
      pmsd->recr=log(1.e-30+pmsd->recr);
      if (parest_flags(400)>0)  //NMD_19May2016
      {
        int ns=pmsd->num_species;
        for (int is=2;is<=ns;is++)
        {
          pmsd->recr(is)(last_real_year-parest_flags(400)+1,last_real_year).initialize();
        }
      }  //NMD_19May2016
      pmsd->nat_mort_coff=log(1.e-20+pmsd->nat_mort_coff);
	}
//NMD 02Nov2011
	// normalize the relative recruitment
  }

  void dvar_fish_stock_history::transform_out()
  {
    fishery_sel=mfexp(fishery_sel);
    {
      int mmin=effort.indexmin();
      int mmax=effort.indexmax();
      for (int i=mmin;i<=mmax;i++)
      {
        if (allocated(effort(i)))
        {
          int mmin=effort(i).indexmin();
          int mmax=effort(i).indexmax();
          for (int j=mmin;j<=mmax;j++)
          {
            if (allocated(effort(i,j)))
            {
              effort(i,j)=exp(effort(i,j));
            }
          }
        } 
      }
    }
    recr=mfexp(recr);
    q0=exp(q0);
    initpop=mfexp(initpop);
    nat_mort_coff=exp(nat_mort_coff);
//NMD 02Nov2011
    if (pmsd){
      pmsd->nat_mort_coff=exp(pmsd->nat_mort_coff);
	}
//NMD 02Nov2011

  }


  int operator > (fishery_freq_record& r1,fishery_freq_record& r2)
  {
    fishery_header_record& h1=r1;
    fishery_header_record& h2=r2;
    return (h1>h2);
  }

  void fishery_freq_record_array::sort(void)
  {
    fishery_freq_record_pointer_array indx(1,nlfreq);
    int i;
    for (i=indexmin();i<=indexmax();i++)
    {
      indx(i)=&(*this)(i);
    }
    fishery_freq_record_pointer tmp;
    for (i=indexmin();i<=indexmax();i++)
    {
      for (int j=indexmin();j<=indexmax()-i;j++)
      {
        if (*indx(j)>*indx(j+1)) // this users the > for fishery-freq-records
        {
          tmp=indx(j);
          indx(j)=indx(j+1);
          indx(j+1)=tmp;
        }
      }
    }
    // Having sorted the pointers we still need to sort the
    // records pointed to
    fishery_freq_record_array tmparr(indexmin(),indexmax(),
      nlint,shlen,filen,nwint,wshlen,wfilen);
    for (i=indexmin();i<=indexmax();i++)
    {
      tmparr(i)=*indx(i);
    }
    for (i=indexmin();i<=indexmax();i++)
    {
      (*this)(i)=tmparr(i);
    }
  }


    freq_record& freq_record::operator = (freq_record& fr1)
    {
      year = fr1.year;
      month= fr1.month;
      week=  fr1.week;
      freq=  fr1.freq;
      return (*this);
    }

    fishery_freq_record& fishery_freq_record::operator = (fishery_freq_record& fr1)
    {
      year=      fr1.year;
      month=   fr1.month;
      week=    fr1.week;
      if (allocated(fr1.species))
      {
        species.allocate(fr1.species.indexmin(),fr1.species.indexmax());
        species=fr1.species;
      }

      if (allocated(fr1.species1))
      {
        species1.allocate(fr1.species.indexmin(),fr1.species.indexmax());
        species1=fr1.species1;
      }

      if (allocated(fr1.species2))
      {
        species2.allocate(fr1.species.indexmin(),fr1.species.indexmax());
        species2=fr1.species2;
      }

      if (allocated(fr1.species3))
      {
        species3.allocate(fr1.species.indexmin(),fr1.species.indexmax());
        species3=fr1.species3;
      }

      if (allocated(fr1.freq))
      {
        if (!allocated(freq))
        {
          freq.allocate(fr1.freq.indexmin(),fr1.freq.indexmax());
        }
        freq=fr1.freq;
      }
      else
      {
        if (allocated(freq))
          freq.deallocate();
      }

      // *********************************************************
      // *********************************************************
      // *********************************************************
      if (allocated(fr1.wfreq))
      {
        if (!allocated(wfreq))
        {
          wfreq.allocate(fr1.wfreq.indexmin(),fr1.wfreq.indexmax());
        }
        wfreq=fr1.wfreq;
      }
      else
      {
        if (allocated(wfreq))
          wfreq.deallocate();
      }
      // *********************************************************
      // *********************************************************
      // *********************************************************

      if (allocated(fr1.age_freq))
        age_freq=fr1.age_freq;


      fishery= fr1.fishery;
      obs_tot_catch = fr1.obs_tot_catch;
      effort = fr1.effort;
      effort_weight = fr1.effort_weight;
      return (*this);
    }

  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& x)
  {
    if (x<50.0) return exp(x);
    return exp(50.0);
  }


  void fishery_freq_record_array::allocate()
  {
    for (int i=indexmin();i<=indexmax();i++)
    {
      ((*this)(i)).allocate(nlint,nwint);
    }
  }

  fishery_freq_record& fishery_freq_record_array::operator [] (int i)
  {
    #ifdef SAFE_ARRAYS
    check_index(index_min,index_max,i," fishery_freq_record& operator [] (int i)");
    #endif
    return ptr[i];
  }

  fishery_freq_record& fishery_freq_record_array::operator () (int i)
  {
    #ifdef SAFE_ARRAYS
      check_index(index_min,index_max,i," freq_record& operator [] (int i)");
    #endif
    return ptr[i];
  }

  // !! Dave F april 22 02   we can do this as they are read in.
  void fishery_freq_record::allocate(int nlint,int nwint)
  {
    //freq.allocate(1,nlint);
    //freq.initialize();
    //wfreq.allocate(1,nwint);
    //wfreq.initialize();
  }

  fishery_freq_record_array::~fishery_freq_record_array()
  {
    if (ptr)
    {
      ptr+=indexmin();
      delete [] ptr;
      ptr=NULL;
    }
  }

  dmatrix colsub(const dmatrix& _M,int lb,int ub)
  {
    ADUNCONST(dmatrix,M)
    int mmin=M.indexmin();
    int mmax=M.indexmax();
    dmatrix tmp(mmin,mmax);

    for (int i=mmin;i<=mmax;i++)
    {
      tmp(i)=M(i)(lb,ub);
    }
    return tmp;
  }

  dvar_matrix colsub(const dvar_matrix& _M,int lb,int ub)
  {
    ADUNCONST(dvar_matrix,M)
    int mmin=M.indexmin();
    int mmax=M.indexmax();
    dvar_matrix tmp(mmin,mmax);

    for (int i=mmin;i<=mmax;i++)
    {
      M(i)(lb,ub);
      tmp(i)=M(i)(lb,ub);
    }
    return tmp;
  }

void dvar_fish_stock_history::fish_flags_sanity_check(void)
{
  ivector ff24=column(fish_flags,24);
  int ff263=check(column(fish_flags,26),3);
  ivector ff75=column(fish_flags,75);
  if (sum(ff24)>0)
  {
    if (sum(ff75)>0)
    {
      grouped_flag_sanity_check(ff24,ff75);
    }
  }
  if (age_flags(177)==0 && parest_flags(155)==0)
  {
    cerr << " age_flags(177)==0 and  age_flags(155)==0 means"
            " mean+devs recruitment parameterisation is mis-specified" << endl;
    ad_exit(1);
  }
  ivector ff92=column(fish_flags,92);
  ivector ff94=column(fish_flags,94);
  if (sum(ff94)>0 && sum(ff92)==0)
  {
    cerr << " override fish_flags(94) active but fish_flags(92)==0 means"
            " no unequal CPUE penalty weights have been assigned" << endl;
    ad_exit(1);
  }
}

void dvar_fish_stock_history::tag_flags_sanity_check(void)
{
  ivector tf1=column(tag_flags,1);
  ivector tf2=column(tag_flags,2);
  int mmax=tf1.indexmax();
  for (int it=1;it<=mmax;it++)
  {
    if (tf1(it)==0)
    {
      cerr << " ERROR: mixing period = 0 for tag release event:  " << it
      << " Since version 2.2.7.5 this is no longer permitted, exiting" << endl;
//      ad_exit(1);
    }
    else
    {
      if (tf2(it)==0)
      {
        cerr << " WARNING: reporting rates for tag release event:  " << it
        <<  " applied during the mixing period. This is not advisable!" << endl;
        cerr << " Consider setting tag_flags(it,2) = 1 to disable  "
        <<  " application during the mixing period" << endl;
      }
    }
  }
  cout << "HERE" << endl;
}

int check_non_zero(const ivector& v)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    if (v(i)) return 1;
  }
  return 0;
}
void dvar_fish_stock_history::data_fish_flags_sanity_check(int index,int group_index)
{
  ivector flags=data_fish_flags(index);
  ivector gp=column(fish_flags,group_index);
  if (sum(gp)>0)
  {
    if (sum(flags)>0)
    {
      grouped_flag_sanity_check(flags,gp,index,group_index);
    }
  }
}

void dvar_fish_stock_history::fish_flags_sanity_check(int index,int group_index)
{
  ivector flags=column(fish_flags,index);
  ivector gp=column(fish_flags,group_index);
  if (sum(gp)>0)
  {
    if (sum(flags)>0)
    {
      grouped_flag_sanity_check(flags,gp,index,group_index);
    }
  }
}

void dvar_fish_stock_history::fish_flags_sanity_check(int index,int group_index,int override_index)
{
  ivector flags=column(fish_flags,index);
  ivector gp=column(fish_flags,group_index);
  ivector ovrd=column(fish_flags,override_index);
  if (sum(gp)>0)
  {
    if (sum(flags)>0)
    {
      grouped_flag_sanity_check(flags,gp,ovrd,index,group_index,override_index);
    }
  }
}


void grouped_flag_sanity_check(const ivector flags,const ivector& gp,int index,int group_index)
{
  int mmin=gp.indexmin();
  int mmax=gp.indexmax();
  if (mmin !=flags.indexmin() || mmax !=flags.indexmax() )
  {
    cerr << "Incompatible vector shapes in"
      " grouped_flag_sanity_check for grouping flags" <<  group_index << " and fish flags " << index << endl;
    ad_exit(1);
  }
  if (!check_non_zero(gp))
    return;
  imatrix M(1,2,mmin,mmax);
  M(1)=gp;
  M(2)=flags;
  M=trans(sort(trans(M),1));
  ivector gp1=M(1);
  ivector fl1=M(2);
  for (int i=mmin+1;i<=mmax;i++)
  {
    // check for holes
    if (gp1(i)>gp1(i-1)+1)
    {
      cerr << "grouping flags hole for grouping flags" <<  group_index << " and fish flags " << index << endl;
      ad_exit(1);
    }
  }
  
  for (int i=mmin+1;i<=mmax;i++)
  {
    if (gp1(i)==gp1(i-1) && fl1(i) !=fl1(i-1))
    {
      cerr << "grouping flags error for grouping flags" <<  group_index << " and fish flags " << index << endl;
      ad_exit(1);
    }
  }
    
}

void grouped_flag_sanity_check(const ivector flags,const ivector& gp,const ivector& ovrd,
       int index,int group_index,int override_index)
{
  int mmin=gp.indexmin();
  int mmax=gp.indexmax();
  if (mmin !=flags.indexmin() || mmax !=flags.indexmax() )
  {
    cerr << "Incompatible vector shapes in"
      " grouped_flag_sanity_check for grouping flags" <<  group_index << " and fish flags " << index << endl;
    ad_exit(1);
  }
  if (!check_non_zero(gp))
    return;
  imatrix M(1,3,mmin,mmax);
  M(1)=gp;
  M(2)=flags;
  M(3)=ovrd;
  M=trans(sort(trans(M),1));
  ivector gp1=M(1);
  ivector fl1=M(2);
  ivector ovr1=M(3);
  for (int i=mmin+1;i<=mmax;i++)
  {
    // check for holes
    if (gp1(i)>gp1(i-1)+1)
    {
      cerr << "grouping flags hole for grouping flags" <<  group_index << " and fish flags " << index << endl;
      ad_exit(1);
    }
  }
  
  for (int i=mmin+1;i<=mmax;i++)
  {
    if (gp1(i)==gp1(i-1) && ovr1(i)==ovr1(i-1) && ovr1(i)!=0)
    {
      //within group over-ride of equal flags is active - do nothing
    }
    else
    {
      if (gp1(i)==gp1(i-1) && fl1(i) !=fl1(i-1))
      {
        cerr << "grouping flags error for grouping flags" <<  group_index << " and fish flags " << index << endl;
        ad_exit(1);
      }
    }
  }
    
}


void grouped_flag_sanity_check(const ivector gp,const ivector& flags)
{
  int mmin=gp.indexmin();
  int mmax=gp.indexmax();
  if (mmin !=flags.indexmin() || mmax !=flags.indexmax() )
  {
    cerr << "Incompatible vector shapes in"
      " grouped_flag_sanity_check " << endl;
    ad_exit(1);
  }
  if (!check_non_zero(gp))
    return;
  imatrix M(1,2,mmin,mmax);
  M(1)=gp;
  M(2)=flags;
  M=trans(sort(trans(M),1));
  ivector gp1=M(1);
  ivector fl1=M(2);
  for (int i=mmin+1;i<=mmax;i++)
  {
    // check for holes
    if (gp1(i)>gp1(i-1)+1)
    {
      cerr << "grouping flags hole " << endl;
      ad_exit(1);
    }
  }
  
  for (int i=mmin+1;i<=mmax;i++)
  {
    if (gp1(i)==gp1(i-1) && fl1(i) !=fl1(i-1))
    {
      cerr << "grouped flags error " << endl;
      ad_exit(1);
    }
  }
    
}

void dvar_fish_stock_history::orthpoly_recr_flags_sanity_check(void)
{
  MY_DOUBLE_TYPE sum_reg;
  MY_DOUBLE_TYPE diff_sum_reg;
  dvector reg_pars;
  if (!pmsd)      // XXXXYYY
  {
    reg_pars=value(region_pars(1));
    sum_reg=sum(reg_pars);
    if (!parest_flags(216) && sum_reg>0)
    {
      cerr << "Error: parest_flags(216)=0 and region_pars(1)>0 "
      << "Not estimating orth-poly recruitment regional "
      << "distribution, but applying a prior penalty function  "
      << "upon a specified regional recruitment distribution."
      << endl;
      cerr << "Exiting... "
      << endl;
      ad_exit(1);
    }
    if (parest_flags(216)>0 && sum_reg>0)
    {
      diff_sum_reg=fabs(1.0-sum_reg);
      if (diff_sum_reg>=0.1)
      {
        cerr << "Error: sum of region_pars(1) not equal to 1 " 
             << "while estimating orth-poly recruitment regional recrs"
             << " 1.0 - sum(region_pars(1)) =  " << diff_sum_reg
             << endl;
        cerr  << "    Exiting... " << endl;
        ad_exit(1);
      }
    }
  }
  else
  {
    for (int is=1;is<=pmsd->num_species;is++)
    {
      int mmin=pmsd->region_bounds(is,1); 
      int mmax=pmsd->region_bounds(is,2); 
      reg_pars=value(region_pars(1)(mmin,mmax));
      sum_reg=sum(reg_pars);
      if (!parest_flags(216) && sum_reg>0)
      {
        cerr << "Error: species/sex =  " << is
        << "; parest_flags(216)=0 and region_pars(1)>0 "
        << "Not estimating orth-poly recruitment regional "
        << "distribution, but applying a prior penalty function  "
        << "upon a specified regional recruitment distribution."
        << endl;
        cerr << "Exiting... "
        << endl;
        ad_exit(1);
      }
      if (parest_flags(216)>0 && sum_reg>0)
      {
        diff_sum_reg=fabs(1.0-sum_reg);
        if (diff_sum_reg>0.1)
        {
          cerr << "Error: species/sex =  " << is
               << "Error: sum of region_pars(1) not equal to 1 " 
               << "while estimating orth-poly recruitment regional recrs"
               << " 1.0 - sum(region_pars(1)) =  " << diff_sum_reg
               << endl;
          cerr  << "    Exiting... " << endl;
          ad_exit(1);
        }
      }
    }
  }
}


void dvar_fish_stock_history::read_fish_flags(cifstream *pinfile)
{
  cifstream & cif = *pinfile;
  if (!parest_flags(197))
  {
    if (!cif) error_msg2("age_flags");

    if (_file_version_no > 1018)
    {
      fish_flags.initialize();
      if (_file_version_no > 1028)
        cif >> fish_flags;
      else if (_file_version_no > 1024)
      {
        for (int i=1;i<=num_fisheries;i++)
          cif >> fish_flags(i)(1,60);
      }
      else
      {
        for (int i=1;i<=num_fisheries;i++)
          cif >> fish_flags(i)(1,50);
      }
    }
    else if (_file_version_no > 1003)
    {
      fish_flags.initialize();
      int mmin=fish_flags.rowmin();
      int mmax=fish_flags.rowmax();
      for (int i=mmin;i<=mmax;i++)
      {
        for (int j=1;j<=40;j++)
        {
          cif >> fish_flags(i,j);
        }
      }
    }
    else 
    {
      fish_flags.initialize();
      int mmin=fish_flags.rowmin();
      int mmax=fish_flags.rowmax();
      for (int i=mmin;i<=mmax;i++)
      {
        for (int j=1;j<=20;j++)
        {
          cif >> fish_flags(i,j);
        }
      }
    }
    if (!cif) error_msg2("fish_flags");
  }
  else
  {
    age_flags.initialize();
    set_default_age_flags(age_flags);
    fish_flags.initialize();
  }
  for (int i=1;i<=num_fisheries;i++)
  {
    if (fish_flags(i,74)==0)
    {
      fish_flags(i,74)=1;
    }
  } 
  if (!age_flags(193) || !pmsd)  //kludge for multi-sex model with length-sel
  {
    fish_flags_sanity_check();
  }
  else
  {
  //  XXXXXXXXXXXXXXXXXXXX  This needs fixing   //NMD_23jun2017
    cerr << " Warning: fish_flags_sanity_check not operating " << endl;
    cerr << " in this case - multi-sex model with length-selectivity " << endl;
  }
}

ivector dvar_fish_stock_history::calculate_blocked_sel_sizes(void)
{
  ivector itmp(1,num_fisheries);
  ivector ff57=column(fish_flags,57);
  ivector ff75=column(fish_flags,75);
  ivector ff61=column(fish_flags,61);
  ivector ff3=column(fish_flags,3);
  for (int i=1;i<=num_fisheries;i++)
  {
    switch (ff57(i))
    {
    case 0:
      if (ff3(i)==0)
      {
        pmsd_error();
        itmp(i)=nage-1-ff75(i);
      }
      else
      {
        itmp(i)=ff3(i)-ff75(i);
      }
      break;
    case 1:
      itmp(i)=2;
      break;
    case 2:
      itmp(i)=3;
      break;
    case 3:
      int num_nodes=0;
      if (ff61(i)>=0)
      {
        num_nodes=ff61(i);
      }
      else
      {
        num_nodes=-ff61(i);
      }
      itmp(i)=num_nodes;
      break;
    }
  }
  return itmp;
}

// New code for checking age_flags for recruitment penalties
// implementation, insert to newmult.cpp:1716 
void dvar_fish_stock_history::age_flags_sanity_check(void)
{
  if (age_flags(129) && parest_flags(149))
  {
    cerr << "Incompatible age_flags(129) and parest_flags(149)"
         << endl;
    ad_exit(1);
  }
  if (age_flags(129) && age_flags(111))
  {
    cerr << "Incompatible age_flags(129) and age_flags(111)"
         << endl;
    ad_exit(1);
  }
  if (age_flags(129) && age_flags(72))
  {
    cerr << "Incompatible age_flags(129) and age_flags(72)"
         << endl;
    ad_exit(1);
  }
  if (age_flags(57) && age_flags(95))
  {
    if (age_flags(57) >= age_flags(95))
    {
      cerr << "Incompatible age_flags(57) and age_flags(95)"
              " age_flags(95) must be >  age_flags(57)"
           << endl;
      ad_exit(1);
    }
  }
}

void dvar_fish_stock_history::modify_Dad(void)
{
  int mmin=Dad.indexmin();
  int mmax=Dad.indexmax();
  int nr1;
  int nr2;
  if (pmsd)   //NMD_aug28-2018
  {
    nr1=Dad(mmin,1).indexmax();
    int nr1_min=Dad(mmin,1).indexmin();    
    nr2=Dad(mmin,1,nr1_min).indexmax();
  }
  else
  {
    nr1=num_regions;
    nr2=num_regions;
  }
//  d3_array tmp(1,nage,1,num_regions,1,num_regions);
  d3_array tmp(1,nage,1,nr1,1,nr2);   //NMD_aug28-2018  
  tmp=value(Dad(mmax));
  for (int i=mmax-1;i>=mmin;i--)
  {
    Dad(i+1)=Dad(i);
  }
  Dad(mmin)=tmp;
}

void dvar_fish_stock_history::modify_move_map(void)
{
  int mmin=move_map.indexmin();
  int mmax=move_map.indexmax();
  int tmp=move_map(mmax);
  for (int i=mmax-1;i>=mmin;i--)
  {
    move_map(i+1)=move_map(i);
  }
  move_map(mmin)=tmp;
}

void dvar_fish_stock_history::modify_diff_coffs(dvar_matrix & m)
{
  int mmin=m.indexmin();
  int mmax=m.indexmax();
  int mmin1=m(mmax).indexmin();
  int mmax1=m(mmax).indexmax();
  dvector tmp(mmin1,mmax1);
  tmp=value(m(mmax));
  for (int i=mmax-1;i>=mmin;i--)
  {
    m(i+1)=m(i);
  }
  m(mmin)=tmp;
}


void dvar_fish_stock_history::modify_all_diff_coffs(void)
{
  modify_diff_coffs(diff_coffs);
  modify_diff_coffs(diff_coffs2);
  modify_diff_coffs(diff_coffs3);
  modify_diff_coffs(diff_coffs_prior);
  modify_diff_coffs(diff_coffs2_prior);
  modify_diff_coffs(diff_coffs3_prior);
}

void dvar_fish_stock_history::get_fishery_group_pointers(void)
{
  int mmin=fish_flags.indexmin();
  int mmax=fish_flags(mmin).indexmax();
  if(!allocated(fishery_group_ptr))
  {
    fishery_group_ptr.allocate(1,mmax);
  }

  
  if ((sum(column(fish_flags,10)) || sum(column(fish_flags,93)) ||
       sum(column(fish_flags,81))) && sum(column(fish_flags,29)))
  {
    make_fishery_group_pointer(29);
  }
  
  if (sum(column(fish_flags,99)))
  {
    make_fishery_group_pointer(99);
  }
  if (sum(column(fish_flags,24)))
  {
    make_fishery_group_pointer(24);
  }
  if (sum(column(fish_flags,68)))
  {
    make_fishery_group_pointer(68);
  }
  if (sum(column(fish_flags,43)))
  {
    make_fishery_group_pointer(44);
  }
}


void dvar_fish_stock_history::make_fishery_group_pointer(int k)
{
  ivector group=column(fish_flags,k);
  int mmin=group.indexmin();
  int mmax=group.indexmax();
  int numgroups=0;
  for (int i=mmin;i<=mmax;i++)
  {
    if (group(i)<0)
    {
      cerr << " group(" << i << ") <= 0 fishery pointer is implausible" 
           << endl;
      ad_exit(1);
    }
  }
  if (sum(group))
  {
    ivector key(mmin,mmax);
    sort(group,key);
    check_group_for_holes(group);
    numgroups=max(group)-min(group)+1;
    if(!allocated(fishery_group_ptr(k)))
    {
      fishery_group_ptr(k).allocate(1,numgroups);
    }
    ivector num_in_group(1,numgroups);
    num_in_group.initialize();
    int ir=0;
    int ic=1;
    for (int i=1;i<mmax;i++)
    {
      if (group(key(i+1))>group(key(i)))
      {
        ir++;
        num_in_group(ir)=ic;
        //group_members(ir).allocate(1,ic);
        fishery_group_ptr(k,ir).allocate(1,ic);
        ic=1;
      }
      else
      {
        ic++;
      }
    }
    num_in_group(numgroups)=ic;
    if(!allocated(fishery_group_ptr(k,numgroups)))
    {
      fishery_group_ptr(k,numgroups).allocate(1,ic);
    }

    int ii=0;
    for (int i=1;i<=numgroups;i++)
    {
      for (int j=1;j<=num_in_group(i);j++)
      {
       fishery_group_ptr(k,i,j)=key(++ii);
      }
    }
  }
}

void dvar_fish_stock_history::do_all_sanity_checks(void)
{
  fish_flags_sanity_check();
  fish_flags_sanity_check(73,29);
//  fish_flags_sanity_check(92,99);
  fish_flags_sanity_check(92,99,94);  //NMD_11April2023
  fish_flags_sanity_check(66,99);
  fish_flags_sanity_check(48,24);
  data_fish_flags_sanity_check(1,99);
  if (parest_flags(155) && !parest_flags(197)) orthpoly_recr_flags_sanity_check(); //NMD_12dec2023
  if (parest_flags(392) && !parest_flags(361))
  {
    cout << "WARNING: parest_flags(392)=1 disables estimation of ALL "
         << " parameters and will crash" << endl;
    cout << " the model run. Recommend enabling parest_flags(361) "
         << " for a single proxy null parameter to be included " << endl;
  }
  if ((parest_flags(377)>0 && !parest_flags(378)) ||
      (parest_flags(378)>0 && !parest_flags(377)))
  {
    cout << "ERROR: parest_flags(377) and parest_flags(378)"
         << " must both be set == 0 or both be set > 0" << endl;
    cout << " Review the requirements for the fml_effort_rltnshp  "
         << " -- the model will now stop " << endl;
    ad_exit(1);
  }
  if (num_tag_releases>0) tag_flags_sanity_check();  
}
