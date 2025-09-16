/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

int check(const ivector& v,int n)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  int flag=0;
  for (int i=mmin;i<=mmax;i++)
  {
    if (v(i)==n)
    {
      flag=1;
      break;
    }
  }
  return flag;
}
  
void block_error()
{
  cerr << "this selectivity option is not modified for blocking" << endl;
  ad_exit(1);
}


extern dvar_len_fish_stock_history * pcfsh;

  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
  dvariable mfexp(dvariable& v1);
  //dvar_vector mfexp(dvar_vector& );
  dvector mfexp(const dvector& );
  extern dvar_vector * psv;
  MY_DOUBLE_TYPE spline_fit(dvar_vector& splcoffs,dvector& oldsel,dvector& xx,
    dvector& x1)
  {
    int mmin=splcoffs.indexmin();
    int mmax=splcoffs.indexmax();
    dvariable mn=log(mean(exp(splcoffs)));
    vcubic_spline_function csf(x1,splcoffs-mn);
    dvariable f=norm2(oldsel-csf(xx));
    f+=1000.*square(mn);
    int i;
    for (i=mmin;i<=mmax;i++)
    {
      if (splcoffs(i)<-19.99)
      {
        f+=1000.* square(splcoffs(i)+19.99);
      }
    }
    return value(f);
  }

  dvector get_new_selcoffs(dvector& oldsel,dvector& tmpsel,dvector& xx,
    dvector& x1)
  {
    int nvar=tmpsel.indexmax();
    independent_variables x(1,nvar);
    x=tmpsel;
    fmm fmc(nvar);
    MY_DOUBLE_TYPE f;
    dvector g(1,nvar);
    g.initialize();
    dvector xbest(1,nvar);
    MY_DOUBLE_TYPE fbest=1.e+50;;
    dvector gbest(1,nvar);
    gbest.fill_seqadd(1.e+50,0.);
    MY_DOUBLE_TYPE maxg;
    fmc.iprint=5;
    fmc.crit=1.e-6;
    gradient_structure::set_YES_DERIVATIVES();
    fmc.scroll_flag=0;
    fmc.imax=40;
    fmc.min_improve=0.;
    fmc.maxfn=100;
    fmc.ifn=0;
    {
      while (fmc.ireturn>=0)
      {
        int badflag=0;
        fmc.fmin(f,x,g);
        {
          maxg=max(g);
        }
        if (fmc.ireturn>0)
        {
          dvar_vector vx=dvar_vector(x);
          f=spline_fit(vx,oldsel,xx,x1);
          gradcalc(nvar,g);
        }
      }
      if (f < fbest)
      {
        fbest=f;
        gbest=g;
        xbest=x;
      }
    }
    return xbest;
  }

  /*
  dvector dvar_fish_stock_history::get_new_orthpolys(void)
  {
    int nvar=recinpop_orth_size_count();
    independent_variables x(1,nvar);
    x.initialize();
    fmm fmc(nvar);
    MY_DOUBLE_TYPE f;
    dvector g(1,nvar);
    g.initialize();
    dvector xbest(1,nvar);
    MY_DOUBLE_TYPE fbest=1.e+50;;
    dvector gbest(1,nvar);
    gbest.fill_seqadd(1.e+50,0.);
    MY_DOUBLE_TYPE maxg;
    fmc.iprint=5;
    fmc.crit=1.e-6;
    gradient_structure::set_YES_DERIVATIVES();
    fmc.scroll_flag=0;
    fmc.imax=40;
    fmc.min_improve=0.;
    fmc.maxfn=100;
    fmc.ifn=0;
    {
      while (fmc.ireturn>=0)
      {
        int badflag=0;
        fmc.fmin(f,x,g);
        {
          maxg=max(g);
        }
        if (fmc.ireturn>0)
        {
          dvar_vector vx=dvar_vector(x);
          f=orth_poly_fit(vx);
          gradcalc(nvar,g);
        }
      }
      if (f < fbest)
      {
        fbest=f;
        gbest=g;
        xbest=x;
      }
    }
    return xbest;
  }
  */


void dvar_fish_stock_history::add_cs_selectivity_coffs(void)
{
  ivector ff74=column(fish_flags,74);
  ivector ff71=column(fish_flags,71);
  int i,j;
  ivector nd=column(fish_flags,3);
  for (i=1;i<=num_fisheries;i++)
  {
    if (nd(i)==0) nd(i)=nage;
  }
  for (i=1;i<=num_fisheries;i++)
  {
    if (fish_flags(i,57)==3)
    {
      if (fish_flags(i,62)>0)
      {
        // old number of nodes
        int sd=fish_flags(i,61);
        //new number of nodes;
        int sd1=fish_flags(i,61)+fish_flags(i,62);
        for (int is=1;is<=ff74(i);is++)
        {
          for (int ib=1;ib<=ff71(i)+1;ib++)
          {
            // old node parameters
            if (bs_selcoff(i,is,ib).indexmax() !=fish_flags(i,61))
            {
              cerr << "This can't happen" << endl;
              ad_exit(1);
            }
            dvector tsel=value(bs_selcoff(i,is,ib));
            tsel-=log(mean(exp(tsel)));
            dvector x(1,sd);
            dvector x1(1,sd1);
            // old x values
            x.fill_seqadd(0,1.0/(sd-1));
            // new x values
            x1.fill_seqadd(0,1.0/(sd1-1));
            cubic_spline_function csf(x,tsel);
            // nodes for current selectivity values
            dvector xx(1,nd(i));
            xx.fill_seqadd(0,1.0/(nd(i)-1));
            dvector oldsel=csf(xx);     
            // initial values for new selpars
            dvector tmpsel=csf(x1);     
            // find new selpars giving best old values for selectivities
            tmpsel=get_new_selcoffs(oldsel,tmpsel,xx,x1);
            bs_selcoff(i,is,ib).deallocate();
            bs_selcoff(i,is,ib).allocate(1,sd1);
            bs_selcoff(i,is,ib)=tmpsel;
            fish_flags(i,61)+=fish_flags(i,62);
            fish_flags(i,62)=0.0;
          }
        }
      }
    }
  }
}

void dvar_fish_stock_history::add_equilibrium_selectivity_coffs(void)
{
  pmsd_error();
  int sd=parest_flags(374);
  int sd1=parest_flags(374)+parest_flags(376);
  parest_flags(374)+=parest_flags(376);
  parest_flags(376)=0;
  ivector mps = kludged_get_equilibrium_movements_per_season();

  for (int is=1;is<=age_flags(57);is++)
  {
    for (int im=1;im<=mps(is);im++)
    {
      for (int ir=1;ir<=num_regions;ir++)
      {
        //nt is=1; int im=1;  int ir=1;
        dvector tsel=value(kludged_equilib_coffs(is,im,ir)(1,sd));
        //tsel-=log(mean(exp(tsel)));
        // note that I think we don't need to add the regional level
        // to fit for the new coffs since both would have the same 
        // thing added to them.
        dvector x(1,sd);
        dvector x1(1,sd1);
        // old x values
        x.fill_seqadd(0,1.0/(sd-1));
        // new x values
        x1.fill_seqadd(0,1.0/(sd1-1));
        cubic_spline_function csf(x,tsel);
        // nodes for current selectivity values
        dvector xx(1,nage);
        xx.fill_seqadd(0,1.0/(nage-1));
        //dvector oldsel=csf(xx);     
        // initial values for new selpars
        dvector tmpsel=csf(x1);     
        // find new selpars giving best old values for selectivities
        //tmpsel=get_new_selcoffs(oldsel,tmpsel,xx,x1);
        kludged_equilib_coffs(is,im,ir)(1,sd1)=tmpsel;
      }
    }
  }
}

void dvar_fish_stock_history::fishing_mortality_calc(d3_array& len_sample_size,
  d3_array& wght_sample_size)
{
  dvar_vector tmp_effort;
  dvar_vector tmp_bio=fish_pars(8);
  ivector ff52=column(fish_flags,52);
  int sflag=sum(ff52);
  if (age_flags(156)) tmp_effort=exp(fish_pars(7));
  if (age_flags(157)) 
  {
    pcfsh->set_global_variance();
    set_global_vars_flag=1;
    pcfsh->mean_lengths_by_year_calc();
    pcfsh->mean_weights_by_year_calc();
    pcfsh->calculate_the_biomass_by_region();
  }
  int nya=1; // num_years_for_average
  //   VVVVVV
  if (age_flags(125))
  {
    nya=age_flags(95); // num_years_for_average
    if (nya==0)
    {
      cerr << "Error age_flags(95) must be >0 " << endl;
      ad_exit(1);
    }
  }


  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
    {
      int yr=year(ir,ip);
      dvar_matrix cd_mult=region_pars.sub(2,nya);
      if (af170q0==0 && age_flags(125) && yr>nya) break;
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        dvar_vector* pfm=0;
        dvar_vector * pcat=0;
        //const prevariable* pcat=0;
        if (af170q0==0)
        {
          pcat=&(catchability(ir,ip));
          pfm=&(fish_mort(ir,ip,fi));
        }
        else
        {
          pcat=&(catchability_q0(ir,ip));
          pfm=&(fish_mort_q0(ir,ip,fi));
        }

        int i=parent(ir,ip,fi);
        int rr=realization_region(i,1);
        int rp=realization_period(i,1);
        int ri=realization_incident(i,1);
        int gfp=global_fishing_periods(ir,ip);
        dvar_vector& fs=incident_sel(ir,ip,fi);
        MY_DOUBLE_TYPE eff=effort(ir,ip,fi);
        if (!age_flags(156)) 
        { 
          *pfm=fs+eff+(*pcat)(fi);
//          if (age_flags(125) && yr>1)
          if (af170q0==0 && age_flags(125) && yr>1)  //NMD_21feb2019
          {
            *pfm+=tmp_bio(i)*cd_mult(yr,ir);
          }
        }
        else
        {
          if (sflag)
#if !defined(NO_MY_DOUBLE_TYPE)
            *pfm=fs+(tmp_effort(i)-1.0L)*grouped_effort(ff52(i),gfp)+eff
#else
            *pfm=fs+(tmp_effort(i)-1.0)*grouped_effort(ff52(i),gfp)+eff
#endif
              +(*pcat)(fi);
          else
            *pfm=fs+tmp_effort(i)*eff+ (*pcat)(fi);
        }
      }
    }
  }
  for (int i=1;i<=num_fisheries;i++)
  {
    dvar_vector* pfm=0;
    if (fish_flags(i,19)>0)
    {
      int num_ages=min(fish_flags(i,19),nage);
      for (int j=1;j<=num_fish_times(i);j++)
      {
        int rr=realization_region(i,j);
        int rp=realization_period(i,j);
        int ri=realization_incident(i,j);
         // VVVVVVVVV
        if (af170q0==0 && age_flags(125) && year(rr,rp)>nya) break;
        if (af170q0==0)
        {
          pfm=&(fish_mort(rr,rp,ri));
        }
        else
        {
          pfm=&(fish_mort_q0(rr,rp,ri));
        }
        if (len_sample_size(rr,rp,ri)>0 || wght_sample_size(rr,rp,ri)>0)
        {
          for (int jj=1;jj<=num_ages;jj++)
          {
            (*pfm)(jj)+=sel_dev_coffs(i,j,jj);
          }
        }
      }
    }
  }

  if (age_flags(34)>0)
  {
    for (int i=1;i<=num_fisheries;i++)
    {
      dvar_vector* pfm=0;
      if (fish_flags(i,4)>0)
      {
        for (int j=1;j<=num_fish_times(i);j++)
        {
          int rr=realization_region(i,j);
          int rp=realization_period(i,j);
         // VVVVVVVVV
          if (af170q0==0 && age_flags(125) && year(rr,rp)>nya) break;
          int ri=realization_incident(i,j);
          if (af170q0==0)
          {
            pfm=&(fish_mort(rr,rp,ri));
          }
          else
          {
            pfm=&(fish_mort_q0(rr,rp,ri));
          }
          (*pfm)+=effort_devs(rr,rp,ri);
        }
      }
    }
  }
}
void dvar_fish_stock_history::put_values_in_fish_mort(void)
{
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
    {
      int yr=year(ir,ip);
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        fish_mort(ir,ip,fi)=1.e+100;
      }
    }
  }
}

void dvar_fish_stock_history::fishing_mortality_calc(d3_array& len_sample_size,
  int break_year)
{
  dvar_vector tmp_effort;
  dvar_vector tmp_bio;
  ivector ff52=column(fish_flags,52);
  int sflag=sum(ff52);
  if (age_flags(156)) tmp_effort=exp(fish_pars(7));
  if (age_flags(157)) 
  {
    pcfsh->set_global_variance();
    set_global_vars_flag=1;
    tmp_bio=fish_pars(8);
    pcfsh->mean_lengths_by_year_calc();
    pcfsh->mean_weights_by_year_calc();
    pcfsh->calculate_the_biomass_by_region();
  }

  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
    {
      int yr=year(ir,ip);
      if (age_flags(125) && yr>1) break;
      if (break_year<yr) break;
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        dvar_vector& fm=fish_mort(ir,ip,fi);
        int i=parent(ir,ip,fi);
        int rr=realization_region(i,1);
        int rp=realization_period(i,1);
        int ri=realization_incident(i,1);
        int gfp=global_fishing_periods(ir,ip);
        dvar_vector& fs=incident_sel(rr,rp,ri);
        MY_DOUBLE_TYPE eff=effort(ir,ip,fi);
        const prevariable& cat=catchability(ir,ip,fi);
        if (!age_flags(156)) 
        { 
          fm=fs+eff+cat;
        }
        else
        {
          if (sflag)
#if !defined(NO_MY_DOUBLE_TYPE)
            fm=fs+(tmp_effort(i)-1.0L)*grouped_effort(ff52(i),gfp)+eff+cat;
#else
            fm=fs+(tmp_effort(i)-1.0)*grouped_effort(ff52(i),gfp)+eff+cat;
#endif
          else
            fm=fs+tmp_effort(i)*eff+cat;
        }
      }
    }
  }
  for (int i=1;i<=num_fisheries;i++)
  {
    if (fish_flags(i,19)>0)
    {
      int num_ages=min(fish_flags(i,19),nage);
      for (int j=1;j<=num_fish_times(i);j++)
      {
        int rr=realization_region(i,j);
        int rp=realization_period(i,j);
        int ri=realization_incident(i,j);
        if (break_year<year(rr,rp)) break;
        if (len_sample_size(rr,rp,ri)>0)
        {
          for (int jj=1;jj<=num_ages;jj++)
          {
            fish_mort(rr,rp,ri,jj)+=sel_dev_coffs(i,j,jj);
          }
        }
      }
    }
  }

  if (age_flags(34)>0)
  {
    for (int i=1;i<=num_fisheries;i++)
    {
      if (fish_flags(i,4)>0)
      {
        for (int j=1;j<=num_fish_times(i);j++)
        {
          int rr=realization_region(i,j);
          int rp=realization_period(i,j);
          int ri=realization_incident(i,j);
          if (break_year<year(rr,rp)) break;
          fish_mort(rr,rp,ri)+=effort_devs(rr,rp,ri);
        }
      }
    }
  }
}

  void dvar_fish_stock_history::catchability_calc(void)
  {
    if(!sum(column(fish_flags,29)))
      catchability_devs_calc();
    else
      grouped_catchability_calc();

    int s47=sum(column(fish_flags,47));
    
    int s27=sum(column(fish_flags,27));

    if (s47 && s27)
    {
      cerr << "You can't have seasonal_catchability and explicit seasonal catchability"
	      " at the same time " << endl;
      ad_exit(1);
    }

    if (s27)
      seasonal_catchability_calc();

    if (s47)
      explicit_seasonal_catchability_calc();

  }

  void dvar_fish_stock_history::catchability_devs_calc(void)
  {
    dvar3_array* pc=0;
    if (af170q0)
    {
      pc=&catchability_q0;
    }
    else
    {
      pc=&catchability;
    }
    for (int i=1;i<=num_fisheries;i++)
    {
      int rr=realization_region(i,1);
      (*pc)(rr,realization_period(i,1),realization_incident(i,1))=q0(i);
      for (int nt=2;nt<=num_fish_times(i);nt++)
      {
        int rr=realization_region(i,nt);
        int rr1=realization_region(i,nt-1);
        int rp=realization_period(i,nt);
        int rp1=realization_period(i,nt-1);
        (*pc)(rr,rp,realization_incident(i,nt))
          =(*pc)(rr1,rp1,realization_incident(i,nt-1))
            + catch_dev_coffs(i,nt);
      }
    }
  }

  void dvar_fish_stock_history::seasonal_catchability_calc(void)
  {
    const MY_DOUBLE_TYPE tpi=2.e0*3.14159e0;
    dvar3_array * pc=0;
    if (af170q0==0)
    {
      pc=&(catchability);
    }
    else
    {
      pc=&(catchability_q0);
    }
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
      {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
        {                                         // incidents for this period
          int i=parent(ir,ip,fi);
          if (fish_flags(i,27))
          {
            (*pc)(ir,ip,fi)+=fish_pars(1,i)
              *sin(tpi*(true_month(ir,ip)/12.-fish_pars(2,i)));
          } 
        }
      }
    }
  }

  void dvar_fish_stock_history::explicit_seasonal_catchability_calc(void)
  {
    dvar3_array * pc=0;
    if (af170q0==0)
    {
      pc=&(catchability);
    }
    else
    {
      pc=&(catchability_q0);
    }
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
      {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
        {                                         // incidents for this period
          int i=parent(ir,ip,fi);
	  int tm=true_month(ir,ip);
	  ivector spair=seasonal_catchability_pars_index(tm);
	  dvector smix=seasonal_catchability_pars_mix(tm);
          if (fish_flags(i,47))
          {
            (*pc)(ir,ip,fi)+=seasonal_catchability_pars(i,spair(1))*smix(1)
              +seasonal_catchability_pars(i,spair(2))*smix(2);
	  }  
        }
      }
    }
  }

  dvariable dvar_fish_stock_history::normalize_seasonal_catchability(void)
  {
    dvariable tmp=0.0;
    for (int fi=1;fi<=num_fisheries;fi++) // Loop over fishing
    {                                         // incidents for this period
      int ii=fish_flags(fi,47);
      if (ii>0)
      {
        dvariable ssum=mean(seasonal_catchability_pars(fi)(1,ii));
        seasonal_catchability_pars(fi)(1,ii)-=ssum;
        tmp+=100.0*square(ssum);
      }
    }  
    return tmp;
  }

  void dvar_fish_stock_history::natural_mortality_calc(void)
  {
//    cout << " Inside natural_mortality_calc " << endl;
    dvariable sv1 = (*psv)(1);
    int nsame=0;
    if (age_flags(81)>1) nsame=age_flags(81)-1;
    dvar_vector nmdevs=age_pars(2)(1,nage-nsame)-
      mean(age_pars(2)(1,nage-nsame));
    for (int iy=1;iy<=nyears;iy++)  // Loop over fishing periods
    {
      for (int j=1;j<=nage;j++)         // Loop over age classes
      {
        nat_mort(iy,j)=nat_mort_coff;
        //nat_mort(iy,j)=-0.5;
      }
      // age-dependent changes in nat mort;
//      if (age_flags(73))
//      {
        nat_mort(iy)(1,nage-nsame)+=nmdevs;
        if (nsame>0)
        {
          nat_mort(iy)(nage-nsame+1,nage)+=nmdevs(nage-nsame);
        }
//      }
      if (age_flags(48)>0)
      {
        MY_DOUBLE_TYPE mm=log(age_flags(49)/10.0);
        for (int j=nage-age_flags(48)+1;j<=nage;j++)  
        {
          nat_mort(iy,j)+= mm;
        }
      }
      if (age_flags(47)>0)
      {
        //for (int j=nage-1;j<=nage;j++)   
        for (int j=1;j<=age_flags(47);j++)  
        {
          nat_mort(iy,j)+=sv1;
        }
      }
    }
    if(pmsd)    //NMD 4Nov2011
    {
      for (int is=2;is<=pmsd->num_species;is++)
      {
        int nsame=0;
        if (pmsd->age_flags(is,81)>1) nsame=pmsd->age_flags(is,81)-1;
        dvar_vector nmdevs=pmsd->age_pars(is,2)(1,pmsd->nage(is)-nsame)-
        mean(pmsd->age_pars(is,2)(1,pmsd->nage(is)-nsame));

        for (int iy=1;iy<=nyears;iy++)  // Loop over fishing periods
        {
          for (int j=1;j<=pmsd->nage(is);j++)         // Loop over age classes
          {
            pmsd->nat_mort(is,iy,j)=pmsd->nat_mort_coff(is);
          }
          // age-dependent changes in nat mort;
//          if (age_flags(73))
//          {
            pmsd->nat_mort(is,iy)(1,pmsd->nage(is)-nsame)+=nmdevs;
            if (nsame>0)
            {
              pmsd->nat_mort(is,iy)(pmsd->nage(is)-nsame+1,pmsd->nage(is))+=
                   nmdevs(pmsd->nage(is)-nsame);
            }
//          }
          if (age_flags(48)>0)
          {
            MY_DOUBLE_TYPE mm=log(age_flags(49)/10.0);
            for (int j=pmsd->nage(is)-age_flags(48)+1;j<=pmsd->nage(is);j++)  
            {
              pmsd->nat_mort(is,iy,j)+= mm;
            }
          }
          if (age_flags(47)>0)
          {
            //for (int j=nage-1;j<=nage;j++)   
            for (int j=1;j<=age_flags(47);j++)  
            {
              pmsd->nat_mort(is,iy,j)+=sv1;
            }
          }
        }
//        cout << "Check of multi-species nat_mort. species: " << is << endl;
//        cout << pmsd->nat_mort(is,nyears) << endl;
      }
    }           //NMD 4Nov2011
    if (age_flags(47)>0)
    {
      //cout << " nat_mort =" << endl << nat_mort(1) << endl;
    }
  }

  void dvar_fish_stock_history::natural_mortality_calc2(void)
  {
    int nsame=0;
    if (age_flags(81)>1) nsame=age_flags(81)-1;
 
    for (int iy=1;iy<=nyears;iy++)  // Loop over fishing periods
    {
      nat_mort(iy)(1,nage-nsame)=age_pars(5)(1,nage-nsame);
      if (nsame>0)
      {
        nat_mort(iy)(nage-nsame+1,nage)+=age_pars(5)(nage-nsame);
      }
    }
    if (pmsd && pmsd->num_species>1)
    {
      for (int is=2;is<=pmsd->num_species;is++)
      {
        int ng=pmsd->nage(is);
        int nsame=0;
        if (pmsd->age_flags(is,81)>1) nsame=pmsd->age_flags(is,81)-1;
        for (int iy=1;iy<=nyears;iy++)  // Loop over fishing periods
        {
          pmsd->nat_mort(is,iy)(1,ng-nsame)=
            pmsd->age_pars(is,5)(1,ng-nsame);
          if (nsame>0)
          {
            pmsd->nat_mort(is,iy)(ng-nsame+1,ng)
               +=pmsd->age_pars(is,5)(ng-nsame);
          }
        }
      }
    }
  }


void dvar_fish_stock_history::proportion_at_age_calc(void)
{
  for (int ir=1;ir<=num_regions;ir++)
  {
    int tmp_nfp;
    if (do_fishery_projections_flag==0)
    {
      tmp_nfp=num_real_fish_periods(ir);  
    }
    else
    {
      tmp_nfp=num_fish_periods(ir);  
    }
    for (int ip=1;ip<=tmp_nfp;ip++)  // Loop over fishing periods
    {
      dvariable totcatch;
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)  // Loop over fishing
      {                                    // incidents for this period
        totcatch=0.;
        
        int mmin=catch(ir,ip,fi).indexmin();
        int mmax=catch(ir,ip,fi).indexmax();
        exp_catch(ir,ip,fi)=mfexp(catch(ir,ip,fi))+1.e-4;
        totcatch=sum(exp_catch(ir,ip,fi));
        if (totcatch<=0.0)
        {
          cout << "catch <=0" << endl;
        }
           
        prop(ir,ip,fi)=mfexp(catch(ir,ip,fi))/totcatch;
        tot_catch(ir,ip,fi)=totcatch;
      }
    }
  }
}

extern dvar_vector * const_vb;
extern dvar_vector * const_var;

void dvar_fish_stock_history::incident_selectivity_calc(void)
{
  //double min_select=.025;
  //double min_select=.0001;
  MY_DOUBLE_TYPE rnage=sqrt(double(nage));
  MY_DOUBLE_TYPE min_select=.000001;
  MY_DOUBLE_TYPE one_plus=1.+min_select;
  MY_DOUBLE_TYPE alpha=min_select/exp(1.);
  dvar_vector& sel_corr = corr_wy;
  dvar_vector root_sel_corr=sqr(1.0-square(sel_corr));
  int num_age=0;
  for (int i=1;i<=num_fisheries;i++)
  {
    if (fish_flags(i,3)>1)
    {
      int rr=realization_region(i,1);
      int rp=realization_period(i,1);
      int ri=realization_incident(i,1);
      // Is selectivity fixed to be 1 for last age classes
      if (fish_flags(i,12)==0) 
      {
        num_age=nage;
      }
      else
      {
        num_age=nage-2;
        for (int j=max(1,nage-1);j<=nage;j++)
        {
          incident_sel(rr,rp,ri,j) = 1.;
        }
      }
      int j;
      for (j=1;j<=num_age;j++)
      {
        incident_sel(rr,rp,ri,j) = exp(selcoff(i,min(fish_flags(i,3),j)));
        if (incident_sel(rr,rp,ri,j)< min_select)
        {
          incident_sel(rr,rp,ri,j)=alpha*exp(incident_sel(rr,rp,ri,j)/
            min_select);
        }
      }

      if (!parest_flags(162))
      {
        incident_sel(rr,rp,ri)/=mean(incident_sel(rr,rp,ri));
      }
      else
      {
        incident_sel(rr,rp,ri)/=(norm(incident_sel(rr,rp,ri))*rnage);
      }

      for (j=1;j<=nage;j++)
      {
        if (incident_sel(rr,rp,ri,j)< min_select)
        {
          incident_sel(rr,rp,ri,j)=alpha*exp(incident_sel(rr,rp,ri,j)/
            min_select);
        }
      }

      for (int nt=2;nt<=num_fish_times(i);nt++)
      {
        int rr=realization_region(i,nt);
        int rr1=realization_region(i,nt-1);
        int rp=realization_period(i,nt);
        int ri=realization_incident(i,nt);
        // fishing selectivity in the next realization of the fishery is the
        // same as this realization plus...
        incident_sel(rr,rp,ri) = incident_sel(rr1,realization_period(i,nt-1),
                                       realization_incident(i,nt-1));
        // ... some random noise which is correlated across age classes
        for (j=1;j<=fish_flags(i,3);j++)
        {
          incident_sel(rr,rp,ri,j)+=delta2(i,nt,j);
          if (incident_sel(rr,rp,ri,j)< min_select)
          {
            incident_sel(rr,rp,ri,j)=alpha*exp(incident_sel(rr,rp,ri,j)/
              min_select);
          }
        }

        for (j=fish_flags(i,3)+1;j<=num_age;j++)
        {
          incident_sel(rr,rp,ri,j) = incident_sel(rr,rp,ri,j-1);
        }
        if (fish_flags(i,12)!=0) 
        {
          for (int j=max(1,nage-1);j<=nage;j++)
          {
            incident_sel(rr,rp,ri,j) = 1.;
          }  
        }
        if (!parest_flags(162))
        {
          incident_sel(rr,rp,ri)/=mean(incident_sel(rr,rp,ri));
        }
        else
        {
          incident_sel(rr,rp,ri)/=(norm(incident_sel(rr,rp,ri))*rnage);
        }

        for (j=1;j<=nage;j++)
        {
          if (incident_sel(rr,rp,ri,j)< min_select)
          {
            incident_sel(rr,rp,ri,j)=alpha*exp(incident_sel(rr,rp,ri,j)/
              min_select);
          }
        }
      }
    }
  }

  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)  // Loop over fishing
      {
        incident_sel(ir,ip,fi)=log(incident_sel(ir,ip,fi));
      }
    }
  }
}


void dvar_fish_stock_history::calculate_nat_mort_by_period(void)
{
  if (!pmsd)
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int ip=1;ip<=num_fish_periods(ir);ip++)//Loop over fishing periods
      {
        e_nat_mort_by_period(ir,ip)=
            mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
      }
    }
  }
  else
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int ip=1;ip<=num_fish_periods(ir);ip++)//Loop over fishing periods
      {
        e_nat_mort_by_period(ir,ip)=
            mfexp(get_nat_mort_region(ir)(year(ir,ip))+fraction(ir,ip));
      }
    }
  }
}

void dvar_fish_stock_history::get_first_unfixed_year(void)
{
  int break_flag=0;
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)
    {
      int ty=really_true_year(ir,ip)+year1-1;
      if (parest_flags(243)==ty)
      {
        first_unfixed_year=year(ir,ip);
        break_flag=1;
        break;
      }
    }
    if (break_flag) break;
  }
  for (int i=1;i<=num_fisheries;i++)
  {
    int break_flag=0;
    for (int j=1;j<=num_fish_times(i);j++)
    {
      int rr=realization_region(i,j);
      int rp=realization_period(i,j);
      int ty=really_true_year(rr,rp)+year1-1;
      if (parest_flags(243)==ty)
      {
        first_unfixed_fish_time(i)=j;
        break_flag=1;
      }
      if (break_flag) break;
    }
  }
}
