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


void dvar_len_fish_stock_history::calculate_the_catch_biomass(void)
{
  int fmin=1;
  int fmax=num_fisheries;
  if (pmsd)
  {
    int cs=pmsd->current_species;
    fmin=(cs-1)*pmsd->num_real_fisheries+1;
    fmax=(cs)*pmsd->num_real_fisheries;
  }
  catch_biomass.initialize();
  catch_biomass_by_fishery.initialize();
  for (int i=fmin;i<=fmax;i++)
  {
    int nft=0;
    if (do_fishery_projections_flag==0)
      nft=num_real_fish_times(i);
    else
      nft=num_fish_times(i);

    for (int nt=1;nt<=nft;nt++)
    {
      int rr=realization_region(i,nt);
      int rp=realization_period(i,nt);
      int ri=realization_incident(i,nt);
      int yr=year(rr,rp);
      if (!age_flags(112))
      {
        catch_biomass(yr)+=exp(catch(rr,rp,ri))*
          mean_weight(rr,rp,ri)/1000.;
        catch_biomass_by_fishery(yr,i)+=exp(catch(rr,rp,ri))*
          mean_weight(rr,rp,ri)/1000.;
      }
      else
      {
        catch_biomass(yr)+=sum(exp(catch(rr,rp,ri)))/1000.;
        catch_biomass_by_fishery(yr,i)+=sum(exp(catch(rr,rp,ri)))/1000.;
      }
    }
  }
}
void dvar_len_fish_stock_history::calculate_the_catch_numbers(void) // JH 27/03/02
{
  int fmin=1;
  int fmax=num_fisheries;
  if (pmsd)
  {
    int cs=pmsd->current_species;
    fmin=(cs-1)*pmsd->num_real_fisheries+1;
    fmax=(cs)*pmsd->num_real_fisheries;
  }

  catch_numbers.initialize();
  for (int i=fmin;i<=fmax;i++)
  {
    for (int nt=1;nt<=num_fish_times(i);nt++)
    {
      int rr=realization_region(i,nt);
      int rp=realization_period(i,nt);
      int ri=realization_incident(i,nt);
      int yr=year(rr,rp);
      catch_numbers(yr)+=sum(exp(catch(rr,rp,ri)));
    }
  }
}

void dvar_len_fish_stock_history::calculate_the_mean_weight(void)
{
  for (int i=1;i<=num_fisheries;i++)
  {
    for (int nt=1;nt<=num_fish_times(i);nt++)
    {
      int ir=realization_region(i,nt);
      int ip=realization_period(i,nt);
      int fi=realization_incident(i,nt);
      dvariable sv27=get_sv_region(ir,27);
      dvariable sv28=get_sv_region(ir,28);
      if (value(sv28)==3.0)
      {
        pmsd_error();
        mean_weight(ir,ip,fi)=len_wt_coff*
          (pow(mean_length(ir,ip,fi),3)+
           3.0*elem_prod(mean_length(ir,ip,fi),vars(ir,ip,fi)));
      }
      else
      {
        for (int j=1;j<=nage;j++)
        {  
          mean_weight(ir,ip,fi,j)=normal_length_to_weight(0.5,-3.5,
            mean_length(ir,ip,fi,j),sqrt(vars(ir,ip,fi,j)),
            value(sv27),value(sv28));
        }
      }
    }
  }
}
void get_msy_pt(dvariable& Bmsy,dvariable & Msy,dvariable& k,
  const dvariable& m, dvariable& r)
{
#if !defined(NO_MY_DOUBLE_TYPE)
  Bmsy=k/pow(m,1.0/(m-1.0L));
#else
  Bmsy=k/pow(m,1.0/(m-1.0));
#endif
  Msy=r*Bmsy*(1.0-1.0/m);
  cout << "MSY = " << Msy << endl;
}


void dvar_len_fish_stock_history::calculate_the_biomass(void)
{
  pmsd_error();
  biomass.initialize();
  adult_biomass.initialize();
  MY_DOUBLE_TYPE lwc=len_wt_coff;
  for (int ir=1;ir<=num_regions;ir++)
  {
    if (!age_flags(112))
    {
      for (int i=1;i<=nyears;i++)  
      {
        biomass(i)+=exp(N(ir,i))*mean_weight_yr(ir,i)/1000.;
        adult_biomass(i)+=elem_prod(exp(N(ir,i)),pmature)*mean_weight_yr(ir,i)/1000.;
      }
    }
    else
    {
      for (int i=1;i<=nyears;i++)  
      {
        biomass(i)+=sum(exp(N(ir,i)))/1000.;
      }
    }
  }
}


void dvar_len_fish_stock_history::calculate_the_biomass_by_region(void)
{
  int ir;
  biomass_by_region.initialize();
  MY_DOUBLE_TYPE lwc=len_wt_coff;
  for (ir=1;ir<=num_regions;ir++)
  {
    if (!age_flags(112))
    {
      for (int i=1;i<=nyears;i++)  
      {
        biomass_by_region(ir,i)+=exp(N(ir,i))*mean_weight_yr(ir,i);
      }
    }
    else
    {
      for (int i=1;i<=nyears;i++)  
      {
        biomass_by_region(ir,i)+=sum(exp(N(ir,i)));
      }
    }
  }
  for (ir=1;ir<=num_regions;ir++)
  {
    rel_biomass_by_region(ir)=
      log(biomass_by_region(ir)/mean(biomass_by_region(ir)));
  }
}

dvariable dvar_len_fish_stock_history::biomass_dynamics_pt(void)
{
  dvariable f;
  dvariable r=sv(22)+0.2;
  dvariable eta=sv(23)+1.0;
// John H. 26/10/2001
  dvariable m1=sv(24)+1.0;
//  dvariable m1=sv(24)+2.0;
  dvariable Bmsy;
  dvariable Msy;
  
  int pen_wt=age_flags(150);
  f=0.0;
  const MY_DOUBLE_TYPE cutoff=0.8;
  dvar_vector pred_biomass(1,nyears+1);
  pred_biomass(1)=biomass(1);
  dvariable k=biomass(1)*eta;
#if !defined(NO_MY_DOUBLE_TYPE)
  get_msy_pt(Bmsy,Msy,k,m1+1.0L,r);
#else
  get_msy_pt(Bmsy,Msy,k,m1+1.0,r);
#endif
// John H. 26/10/2001
  if (!age_flags(153) && m1==1)
//  if (!age_flags(153) && m1==2)
  {
    for (int i=1;i<=nyears;i++)  
    {
      dvariable F=catch_biomass(i)/biomass(i);
      if (F>cutoff)
      {
        F=.99-posfun(.99-F,cutoff,f);
      }
      pred_biomass(i+1)=pred_biomass(i)*(1.0+r)/
        (1.0+r*pred_biomass(i)/(pred_biomass(1)*eta)+F);
    }
  }
  else
  {
    for (int i=1;i<=nyears;i++)  
    {
      dvariable F=catch_biomass(i)/biomass(i);
      if (F>cutoff)
      {
        F=.99-posfun(.99-F,cutoff,f);
      }
      dvariable tmp= r*pow(pred_biomass(i)/(pred_biomass(1)*eta),m1);
      pred_biomass(i+1)=pred_biomass(i)*(1.0+r)/(1.0+tmp+F);
    }
  }
  f+=pen_wt*norm2(log(elem_div(0.1+biomass,0.1+pred_biomass(1,nyears))));
  f+=square(log(m1));
  return f;
}
  
void dvar_len_fish_stock_history::yield_analysis_pt(ofstream * pof)
{
  if (pmsd)
  {
    cerr << "dvar_len_fish_stock_history::yield_analysis_pt(ofstream * pof)"
         << " not modified for multi-species" << endl;
    //ad_exit(1);
    return;
  }
  dvariable r=sv(22)+0.2;
  dvariable eta=sv(23)+1.0;
  dvariable k=biomass(1)*eta;
  dvariable m=sv(24)+2.0;
  get_fishing_mortality_by_age_by_year();
  dvar_matrix& Fay=F_by_age_by_year;
  dvar_vector n1(1,nage);
  dvector tmp1(0,500);
  dvar_vector tmp2(0,500);
  int i;
  int ii=1;
  tmp1(0)=0.0;
  tmp2(0)=0.0;
  int navg=age_flags(148);
  int tmult=age_flags(57);
  if (!navg)navg=tmult;

  dvar_vector F(1,nage);
  F.initialize();
  for (i=1;i<=navg;i++)
    F+=Fay(nyears-i+1);
  F/=double(navg);

  dvar_vector nm=get_nat_mort_species();
  for (i=1;i<=500;i++)
  {
    MY_DOUBLE_TYPE lambda=i/10.0;
    dvar_vector Z=lambda*F+exp(nm);
    n1(1)=1.0;
    for (int j=2;j<=nage-1;j++)  
    {
      n1(j)=n1(j-1)*exp(-Z(j-1));
    }
    n1(nage)=n1(nage-1)*exp(-Z(nage-1))/(1.0-exp(-Z(nage)));
    
    dvariable b1;
    b1=n1*mean_weight_yr(1,1)/1000.;
   
    dvar_vector C1=elem_prod(elem_div(lambda*F,Z),
       elem_prod(1.0-exp(-Z),n1));

    dvariable c1=C1*mean_weight_yr(1,1)/1000.;

#if !defined(NO_MY_DOUBLE_TYPE)
    //dvariable n=pow(k/b1,m-1.0L)*(1.0-c1/(r*b1));
#else
    //dvariable n=pow(k/b1,m-1.0)*(1.0-c1/(r*b1));
#endif
    dvariable n;
    if (c1<r*b1)
#if !defined(NO_MY_DOUBLE_TYPE)
      n=k/b1*pow(1.0-c1/(r*b1),1.0/(m-1.0L));
#else
      n=k/b1*pow(1.0-c1/(r*b1),1.0/(m-1.0));
#endif
    else
      n=0.0;

    dvariable Cn=n*c1;

    if (value(Cn)>0.0)
    {
      tmp1(ii)=lambda;
      tmp2(ii++)=Cn;
    } 
    else 
    {
      tmp1(ii)=lambda;
      tmp2(ii++)=0.0;
      break;
    }
  } 
  
  dvariable Bmsy;
  dvariable Msy;
  get_msy_pt(Bmsy,Msy,k,m,r);
  predicted_yield_pt.deallocate();
  predicted_yield_pt.allocate(0,ii-1);
  predicted_yield_pt=tmp2(0,ii-1);
  predicted_yield_pt_x.deallocate();
  predicted_yield_pt_x.allocate(0,ii-1);
  predicted_yield_pt_x=tmp1(0,ii-1);
  if (pof) {
    (*pof) << "# Pella_t yield analysis report " 
           << "  Bmsy = " << Bmsy << "  "  << "  Msy = " << Msy << endl;
    (*pof) << setw(10) << tmp1(0,ii-1) << endl;
    (*pof) << setw(10) << tmp2(0,ii-1) << endl;
  }
} 



void dvar_len_fish_stock_history::get_fishing_mortality_by_age_by_year(void)
{
  //dvariable r=sv(22)+0.2;
  int cs=get_current_species();
  int ng=get_current_nage();
  dvariable fpen=0.0;
  if (!allocated(F_by_age_by_year))
    F_by_age_by_year.allocate(1,nyears,1,nage);
  ivector itmp=get_region_bounds();
  int rmin=itmp(1);
  int rmax=itmp(2);
  if (pmsd)
  {
    if (!allocated(pmsd->F_by_age_by_year))
    {
      ivector nn(2,pmsd->num_species);
      nn(2,pmsd->num_species)=pmsd->nage(2,pmsd->num_species);
      pmsd->F_by_age_by_year.allocate(2,pmsd->num_species,1,nyears,1,nn);
    }
  }
  
  int iy,ir;
  // ----  old stuff ----------    
  dvar3_array nwts(rmin,rmax,1,nyears,1,ng);
  
  //int ir,iy;
  for (ir=rmin;ir<=rmax;ir++)
  {
    if (af170q0 != 0)   //NMD_17may2022
    {
      nwts(ir)=exp(N_q0(ir));
    }
    else
    {
      nwts(ir)=exp(N(ir));
    }
//    nwts(ir)=exp(N(ir));
  }

  if (pmsd==0)
  {
    for (iy=1;iy<=nyears;iy++)
    {
      for (int j=1;j<=ng;j++)
      {
        dvariable ssum=0.0;
        for (ir=rmin;ir<=rmax;ir++)
          ssum+=nwts(ir,iy,j);
        for (ir=rmin;ir<=rmax;ir++)
          nwts(ir,iy,j)/=(1.e-10+ssum);
      }
    }
  }
  else
  {
    for (iy=1;iy<=nyears;iy++)
    {
      for (int j=1;j<=ng;j++)
      {
        dvariable ssum=0.0;
        for (ir=rmin;ir<=rmax;ir++)
          ssum+=nwts(ir,iy,j);
        for (ir=rmin;ir<=rmax;ir++)
          nwts(ir,iy,j)/=ssum;
      }
    }
  }

// definition of fishery-specific F multipliers for yield analysis
  if (cs==1)
  {
    F_by_age_by_year.initialize();
  }
  else
  {
    pmsd->F_by_age_by_year.initialize();
  }
  if (colsum(fish_flags,70) == 0)
  {
    for (ir=rmin;ir<=rmax;ir++)
    {
      for (int ip=1;ip<=num_fish_periods(ir);ip++) 
      {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {                                 
          int yr=year(ir,ip);
          if (cs==1)
          {
            if (af170q0 != 0)   //NMD_17may2022
            {
              F_by_age_by_year(yr)+=
                elem_prod(exp(fish_mort_q0(ir,ip,fi)),nwts(ir,yr));
            }
            else
            {
              F_by_age_by_year(yr)+=
                elem_prod(exp(fish_mort(ir,ip,fi)),nwts(ir,yr));
            }
//            F_by_age_by_year(yr)+=
//               elem_prod(exp(fish_mort(ir,ip,fi)),nwts(ir,yr));
          }
          else
          {
            if (af170q0 != 0)   //NMD_17may2022
            {
              pmsd->F_by_age_by_year(cs,yr)+=
                 elem_prod(exp(fish_mort_q0(ir,ip,fi)),nwts(ir,yr));
            }
            else
            {
              pmsd->F_by_age_by_year(cs,yr)+=
                 elem_prod(exp(fish_mort(ir,ip,fi)),nwts(ir,yr));
 
            }
//            pmsd->F_by_age_by_year(cs,yr)+=
//               elem_prod(exp(fish_mort(ir,ip,fi)),nwts(ir,yr));
          }
        }
      }
    }
  }
  else
  {
    dvar_vector fmult(1,num_fisheries);
    for (int iff=1;iff<=num_fisheries;iff++)
    {
      fmult(iff)=fish_flags(iff,70)/100.;
    }
    cout << "fmult = " << fmult;
    for (ir=1;ir<=num_regions;ir++)
    {
      for (int ip=1;ip<=num_fish_periods(ir);ip++) 
      {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {                                 
          int yr=year(ir,ip);
          if (cs==1)
          {
            if (af170q0 != 0)   //NMD_17may2022
            {
              F_by_age_by_year(yr)+=elem_prod(
                 exp(fish_mort_q0(ir,ip,fi))*fmult(parent(ir,ip,fi)),nwts(ir,yr));
            }
            else
            {
              F_by_age_by_year(yr)+=elem_prod(
                 exp(fish_mort(ir,ip,fi))*fmult(parent(ir,ip,fi)),nwts(ir,yr));
            }
//            F_by_age_by_year(yr)+=elem_prod(
//               exp(fish_mort(ir,ip,fi))*fmult(parent(ir,ip,fi)),nwts(ir,yr));
          }
          else
          {
            if (af170q0 != 0)   //NMD_17may2022
            {
              pmsd->F_by_age_by_year(cs,yr)+=elem_prod(
                 exp(fish_mort_q0(ir,ip,fi))*fmult(parent(ir,ip,fi)),nwts(ir,yr));
            }
            else
            {
              pmsd->F_by_age_by_year(cs,yr)+=elem_prod(
                 exp(fish_mort(ir,ip,fi))*fmult(parent(ir,ip,fi)),nwts(ir,yr));
            }    
//            pmsd->F_by_age_by_year(cs,yr)+=elem_prod(
//               exp(fish_mort(ir,ip,fi))*fmult(parent(ir,ip,fi)),nwts(ir,yr));
          }
        }
      }
    }
  }
}

void dvar_len_fish_stock_history::get_fishing_mortality_by_age_by_year(int is)
{
  if (is>1 && pmsd==0)
  {
    cerr << "Can not have species value > 1 when pmsd=0" << endl;
    ad_exit(1);
  }
  int rmin,rmax;
  if (pmsd==0)
  {
    rmin=1;
    rmax=num_regions;
  }
  else
  {
    rmin=pmsd->region_bounds(is,1);
    rmax=pmsd->region_bounds(is,2);
    if (!allocated(pmsd->F_by_age_by_year))
    {
      ivector nn(2,pmsd->num_species);
      nn=pmsd->nage(2,pmsd->num_species);
      pmsd->F_by_age_by_year.allocate(2,pmsd->num_species,1,nyears,1,nn);
    }
  }
  //dvariable r=sv(22)+0.2;
  dvariable fpen=0.0;
  if (!allocated(F_by_age_by_year))
    F_by_age_by_year.allocate(1,nyears,1,nage);

  //--- new stuff: PK and JH, Nov. 2002 
  
  int iy,ir;
  // ----  old stuff ----------    
  dvar3_array nwts(1,num_regions,1,nyears,1,nage);
  
  //int ir,iy;
  for (ir=rmin;ir<=rmax;ir++)
  {
    if (af170q0 != 0)   //NMD_17may2022
    {
      nwts(ir)=exp(N_q0(ir));
    }
    else
    {
      nwts(ir)=exp(N(ir));
    }
//    nwts(ir)=exp(N(ir));
  }

  for (iy=1;iy<=nyears;iy++)
  {
    for (int j=1;j<=nage;j++)
    {
      dvariable ssum=0.0;
      for (ir=rmin;ir<=rmax;ir++)
        ssum+=nwts(ir,iy,j);
      for (ir=rmin;ir<=rmax;ir++)
        nwts(ir,iy,j)/=ssum;
    }
  }

// definition of fishery-specific F multipliers for yield analysis
  if (is==1)
  {
    F_by_age_by_year.initialize();
  }
  else
  {
    pmsd->F_by_age_by_year(is).initialize();
  }
  
  if (colsum(fish_flags,70) == 0)
  {
    for (ir=rmin;ir<=rmax;ir++)
    {
      for (int ip=1;ip<=num_fish_periods(ir);ip++) 
      {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {                                 
          int yr=year(ir,ip);
          if (is==1)
          {
            if (af170q0 != 0)   //NMD_17may2022
            {
              F_by_age_by_year(yr)+=
                elem_prod(exp(fish_mort_q0(ir,ip,fi)),nwts(ir,yr));
            }
            else
            {
              F_by_age_by_year(yr)+=
                elem_prod(exp(fish_mort(ir,ip,fi)),nwts(ir,yr));
            }
//            F_by_age_by_year(yr)+=
//               elem_prod(exp(fish_mort(ir,ip,fi)),nwts(ir,yr));
          }
          else
          {
            if (af170q0 != 0)   //NMD_17may2022
            {
              pmsd->F_by_age_by_year(is,yr)+=
                 elem_prod(exp(fish_mort_q0(ir,ip,fi)),nwts(ir,yr));
            }
            else
            {
              pmsd->F_by_age_by_year(is,yr)+=
                 elem_prod(exp(fish_mort(ir,ip,fi)),nwts(ir,yr));
 
            }
//            pmsd->F_by_age_by_year(is,yr)+=
//               elem_prod(exp(fish_mort(ir,ip,fi)),nwts(ir,yr));
          }
        }
      }
    }
  }
  else
  {
    dvar_vector fmult(1,num_fisheries);
    for (int iff=1;iff<=num_fisheries;iff++)
    {
      fmult(iff)=fish_flags(iff,70)/100.;
    }
    cout << "fmult = " << fmult;
    for (ir=rmin;ir<=rmax;ir++)
    {
      for (int ip=1;ip<=num_fish_periods(ir);ip++) 
      {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {                                 
          int yr=year(ir,ip);
          if (is==1)
          {
            if (af170q0 != 0)   //NMD_17may2022
            {
              F_by_age_by_year(yr)+=elem_prod(
                 exp(fish_mort_q0(ir,ip,fi))*fmult(parent(ir,ip,fi)),nwts(ir,yr));
            }
            else
            {
              F_by_age_by_year(yr)+=elem_prod(
                 exp(fish_mort(ir,ip,fi))*fmult(parent(ir,ip,fi)),nwts(ir,yr));
            }
//            F_by_age_by_year(yr)+=elem_prod(
//               exp(fish_mort(ir,ip,fi))*fmult(parent(ir,ip,fi)),nwts(ir,yr));
          }
          else
          {
            if (af170q0 != 0)   //NMD_17may2022
            {
              pmsd->F_by_age_by_year(is,yr)+=elem_prod(
                 exp(fish_mort_q0(ir,ip,fi))*fmult(parent(ir,ip,fi)),nwts(ir,yr));
            }
            else
            {
              pmsd->F_by_age_by_year(is,yr)+=elem_prod(
                 exp(fish_mort(ir,ip,fi))*fmult(parent(ir,ip,fi)),nwts(ir,yr));
            }    
//            pmsd->F_by_age_by_year(is,yr)+=elem_prod(
//               exp(fish_mort(ir,ip,fi))*fmult(parent(ir,ip,fi)),nwts(ir,yr));
          }
        }
      }
    }
  }
}


void dvar_len_fish_stock_history::
  get_fishing_mortality_by_age_by_year_by_region(void)
{
  F_by_age_by_year_by_region.deallocate();
  F_by_age_by_year_by_region.allocate(1,num_regions,1,nyears,1,nage);
  F_by_age_by_year_by_region.initialize();
  int ir;

  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++) 
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {                                 
        if (af170q0 != 0)   //NMD_17may2022
        {
          F_by_age_by_year_by_region(ir,year(ir,ip))+=
                   exp(value(fish_mort_q0(ir,ip,fi)));
        }
        else
        {
          F_by_age_by_year_by_region(ir,year(ir,ip))+=
                   exp(value(fish_mort(ir,ip,fi)));
        }
//        F_by_age_by_year_by_region(ir,year(ir,ip))+=
//                 exp(value(fish_mort(ir,ip,fi)));
      }
    }
  }
}

void dvar_len_fish_stock_history::calculate_the_biomass_by_region(int ir,
  int i)
{
  if (!age_flags(112))
  {
    biomass_by_region(ir,i)=exp(N(ir,i))*mean_weight_yr(ir,i);
  }
  else
  {
    biomass_by_region(ir,i)=sum(exp(N(ir,i)));
  }
  biomass_by_region(ir,i)=log(biomass_by_region(ir,i));
}

