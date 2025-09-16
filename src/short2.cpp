/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
#if !defined(ZFEV)
extern dvar_len_fish_stock_history * pcfsh;
#endif
dvar_vector *  psv=NULL;
  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
int do_maturity_from_length_to_age_flag=1;

  long int diff(long int& x, const long int& y)
  {
     long int tmp=y-x;
     x=y;
     return tmp;
  }
MY_DOUBLE_TYPE& value(MY_DOUBLE_TYPE& x)
{
  return x;
}
void setup_gml(int nage,dvar_vector& lcoffs,dvar_vector& gml,const prevariable& sv20);
void setup_gmlxx(int nage,dvar_vector& vb_coffs,dvar_vector& gml,
  dvar_vector& sv);
extern mf_pvm_manager * mf_pvm;

void dvar_fish_stock_history::get_exp_N(void)
{
  //fsh.exp_N.allocate(1,num_regions,1,fsh.nyears,1,nage),
  dvar3_array * pN=0;
  if (af170q0==0)
  {
    pN=&N;
  }
  else
  {
    pN=&N_q0;
  }
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int i=1;i<=nyears;i++)
    {
      exp_N(ir,i)=mfexp((*pN)(ir,i));
    }
  }
}

void dvar_len_fish_stock_history::do_everything_calc(d3_array& len_sample_size,
  dvar_vector& vb_coff,dvar_vector& var_coff,dvar_vector& sv,
  dvariable& ffpen,ivector *pq_flag)
  {

// NMD_8mar2023
    if (parest_flags(157) && !pq_flag)
    {
      dvar3_array * pN=0;
      pN=&N;
      if (parest_flags(155)==0)
      {
        recinpop_standard(pq_flag,pN,4);
      }
      else if (parest_flags(155)<0)
      {
        orth_pf155_less_than_zero_why(pN);
      }
    }
// NMD_8mar2023

//    ff263flag=check(column(fish_flags,26),3);  //NMD_13jan2023
    if (ff263flag)
    {
#  if !defined(ZFEV)
      if (!parest_flags(143))
#  endif
      {
        if (!parest_flags(175))
          mean_length_calc(0);
        else
          mean_length_calcx();
      }
    }
    switch (age_flags(188))
    {
    case 0:
      break;
    case 1:
      pmature=maturity_length_to_age_simple_spline(1);
      //cout << pmature << endl;
      if (pmsd)   //NMD_19Sep2018
      {
        int num_species=pmsd->num_species;
        for (int is=2;is<=num_species;is++)
        {
          pmsd->pmature(is)=maturity_length_to_age_simple_spline(is);
        }
      }   //NMD_19Sep2018
      break;
    case 2:
      pmature=maturity_length_to_age_weighted_spline(1);
      //cout << pmature << endl;
      if (pmsd)   //NMD_19Sep2018
      {
        int num_species=pmsd->num_species;
        for (int is=2;is<=num_species;is++)
        {
          pmsd->pmature(is)=maturity_length_to_age_weighted_spline(is);
        }
      }   //NMD_19Sep2018
      break;
    default:
      cerr << "illegal value for age_flags(188)" << endl;
      ad_exit(1);
    }

    greport("before incident_selectivity_calc"); 
     new_incident_selectivity_calc(len_sample_size,vb_coff,var_coff);
    greport("after incident_selectivity_calc"); 

    int i=1;

    if (!age_flags(92))    //NMD_19jun2024
    {
    // DF moved this outside Feb 05 05
      catchability_calc();
     // set catchability to 0 in selected fisheries
      set_catchabilities_to_zero(pq_flag);
    }

     if (!age_flags(92)) 
     {
       if (age_flags(34)>0)
       {
         effort_devs_calc();
       }
       //if (age_flags(125)==0)   // cobb douglas
       {
         if (age_flags(157)==0 || age_flags(94)==3)
           fishing_mortality_calc(len_sample_size,wght_sample_size);
         else if (age_flags(94)>0)
           fishing_mortality_calc(len_sample_size,age_flags(95));
       }
     }
    /*
     if (age_flags(188))
     {
       int mult=100;
       dvar_vector alc=age_at_length_calc(nlint,shlen,filen,vb_coff,
         nage,parest_flags,mult);
       //cpmature_at_length=1.0;

       dvar_vector maturity=maturity_length_to_age(alc,cpmature_at_length,
         nage,nlint,mult);
       cout << maturity << endl;
     }
    */

     if (have_projection_periods_flag)
     {
       mean_length_calc(1);
       mean_weights_calc(1);
     }

    greport("after fishing_mortality_calc"); 
     switch (age_flags(109))
     {
     case 0:
       natural_mortality_calc();
       break;
     case 1:
       natural_mortality_calc2();
       break;
     case 2:
       natural_mortality_splines();
       break;
     case 3:
       natural_mortality_lorenzen();
       break;
     case 4:
       natural_mortality_double_normal();
       break;
     default:
       cerr << "Illegal value for age_flags(109)" << endl;
       break;
     }

     calculate_nat_mort_by_period();

     if (num_tag_releases)
     {
       rep_rate_devs_calc();
     }
     /* can't have this here as tot_mort has not been calculated*/
    /*
     if (!age_flags(92) && !age_flags(157))
     {
       do_fish_mort_intermediate_calcs();
     }
    */
     age_flags(181)=0;
     if (!age_flags(54))
     {
       greport("before catch_equations_calc");
       //if (age_flags(158))
       //{
       //  catch_equations_calc_movement(sv);
       //}
       //else
       {
         if (!age_flags(92))
         {
           imppen=0.0;
           if (!age_flags(157))
           {
             if (ad_comm::pthread_manager && threaded_tag_flag==1)
             {
               catch_equations_calc1(sv,pq_flag,imppen);
               int mmin=ad_comm::pthread_manager->gmin(1);
               int mmax=ad_comm::pthread_manager->gmax(1);

               for (int i=mmin;i<=mmax;i++)
               {
                 send_variables_to_tag_slave1(0,i);
               }
               catch_equations_calc2(sv,pq_flag,imppen);
             }
             else
             {
               catch_equations_calc(sv,pq_flag,imppen);
             }
             if (imppen>0.0)
             {
               cout << " AA Excess fishing mortality penalty = " 
                    << imppen << endl;
             }
               cout << " BB Excess fishing mortality penalty = " 
                    << imppen << endl;     //NMD 27Sep2012
             ffpen+=imppen;
           }
           else
           {
             xcatch_equations_calc(sv);
           }
         }
         else
         {
           imppen=0.0;
           age_flags(181)=0;
           switch (age_flags(92))
           {
             case 1:
               // Don't allow this any more
               cerr << "Illegal value for age-flags(92)" << endl;
               ad_exit(1);  
               catch_equations_calc_implicit(sv,imppen);
               ffpen+=imppen;
               do_fish_mort_intermediate_calcs();
               break;
             case 2:
             case 3:
             case 4:
               new_catch_equations_calc_implicit_experiment_loop
                 (sv,pq_flag,imppen);
               ffpen+=imppen;
               break;
             default:
                cerr << "Illegal value for age-flags(92)" << endl;
                ad_exit(1);
            }
           // need this because it is not calculated above
           //do_fish_mort_intermediate_calcs();
           if (imppen>0.0)
           {
             if (age_flags(181)==2)
             {
               cout << "Fish not conserved --";
             }
             age_flags(181)=1;
           }
         }
       }
       if (parest_flags(241) && simulation_seeds(5)  && !pq_flag
	   && projection_sim_index > 0 && projected_simulated_data_flags(1)) //NMD_1Dec2017
       {
         generate_simulated_tags();
         if (parest_flags(242))
         {
           dvar_vector alc=age_at_length_calc(nlint,shlen,filen,vb_coff,
             nage,parest_flags,1);

           generate_real_simulated_tagsx();
         }
       }
       if (!pq_flag)
       {
         if (!ad_comm::pthread_manager || threaded_tag_flag==0)
         {

           if (num_tag_releases && age_flags(122)==0)
           {
             if (!age_flags(96))
             {
               MY_DOUBLE_TYPE svffpen;     //NMD 27Sep2012
               svffpen = value(ffpen);  //NMD 27Sep2012
               ffpen+=tag_catch_equations_calc(sv);
             }
             else
             {
               ffpen+=xtag_catch_equations_calc_pooled(sv);
    #if defined(USE_ADPVM)
               if (mf_pvm->pvm_switch==2)
               {
                 send_dv3_to_master(epooled_tagnum_fish_recr);   
               }
               if (mf_pvm->pvm_switch==1)
               {
                 // add pooled tags from slaves);   
                 mfget_dv3_sum_from_slaves(epooled_tagnum_fish_recr);   
               }
    #endif
             
               if (mf_pvm->pvm_switch==0 || mf_pvm->pvm_switch==1)
               {
                 xpooled_tag_catch_equations_calc(sv);
                 if (parest_flags(241) && simulation_seeds(5)
                     && projection_sim_index > 0)
                 {
                   //sim_xpooled_tag_catch_equations_calc(sv);
                 }
               }
             }
           }
         }
       }
       greport("after catch_equations_calc");
     }
    greport("before proportion_at_age_calc"); 
     proportion_at_age_calc();
    greport("after proportion_at_age_calc"); 
    get_exp_N();
  }
  void dvar_len_fish_stock_history::do_everything_calc(dvariable& ffpen,
    ivector* xpq_flag)
  {
    long int tmp=gradient_structure::totalbytes();
    long int tmp1;
    if (parest_flags(175))
      setup_gml(nage,age_pars(4),gml,sv(20));
    if (parest_flags(174))
      setup_gmlxx(nage,vb_coff,gml,sv);

    psv=&sv;
    greport("before fish_stock_history_do_every"); 
    dvar_vector rel_area(1,num_regions);
    if (num_tag_releases>0)
    {
      if (!ad_comm::pthread_manager || threaded_tag_flag==0)
      {
        var_convert_tag_lengths_to_age();
        if (!age_flags(96))
          observed_tags_by_age_from_length();
        else
          observed_tags_by_age_from_length_pooled();
      }
    }
      
    // proportions by region for recruitment
    get_pop_delta();

    variance_calc();

    if (af170q0==0 && xpq_flag==0)
    {
      dvar_len_fish_stock_history::do_everything_calc(len_sample_size,
         vb_coff,var_coff,sv,ffpen,0);
    }
    else if (af170q0==1 || xpq_flag )
    {
  
      ivector q_flag=column(fish_flags,55);
      ivector * pq_flag=&q_flag;
      //if (age_flags(170) && pq_flag)
      {
        if (sum(*pq_flag)==0)
        {
          cerr << "Incompatible flags " << endl;
          ad_exit(1);
        }
        dvar_len_fish_stock_history::do_everything_calc(len_sample_size,
           vb_coff,var_coff,sv,ffpen,pq_flag);
      }
    }

    if (ad_comm::pthread_manager && threaded_tag_flag==1)
    {
      if (ad_comm::pthread_manager->ngroups>2)
      {
        int mmin=ad_comm::pthread_manager->gmin(3);
        int mmax=ad_comm::pthread_manager->gmin(3);
        for (int i=mmin;i<=mmax;i++)
        {
          send_variables_to_tag_slave3(0,i);
        }
      }
    }


     tmp1=diff(tmp,gradient_structure::totalbytes());
#if !defined(ZFEV)
     if (!parest_flags(143))
#endif
     {
       if (!parest_flags(175))
         ;//mean_length_calc();
       else
         mean_length_calcx();
        // don't need this for just catch at age data?

       if (parest_flags(240))
       {
         all_age_length_calcs();
       }

       if (age_flags(40)==0)
       {
         if (!age_flags(61))
         {
           fast_pred_frequency_calc();
         }
         else
         {
           fast_pred_frequency_calc_len_based();
           length_basedsel_calc(lengthbsel,*this,vb_coff,var_coff,
             nlint,fmid);
           for (int i=1;i<=num_fisheries;i++)
           {
             MY_DOUBLE_TYPE  tt=value(mean(lengthbsel(i)));
             cout <<"short2.cpp " << tt << endl;
           }
         }
       }
       if (parest_flags(311))
       {
         // compute tail compressed predcited length frequencies
         tail_compress_predicted_length_frequencies();
       }
       if (age_flags(61))
       {
         big_fish_catch();
       }
     }
#if !defined(ZFEV)
     else
     {
        // do this once at the beginning of the run
        if (parest_flags(144))
        {
        pcfsh->mean_length_calc();
        pcfsh->variance_calc();
      }
      if (age_flags(40)==0)
      {
        cfast_pred_frequency_calc(*pcfsh);
      }
    }
 #endif
    if (nwint && !parest_flags(181) ) 
    {
      //mean_weights_calc(0);    // mode flag =0 because not projection
      fast_weight_pred_frequency_calc();
    }

    if (parest_flags(301))
    {
      // compute tail compressed predicted wieght frequencies
      tail_compress_predicted_weight_frequencies();
    }
  /*  
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
      {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
        {                                         // incidents for this period
          ffpen+=log(norm2(wtprob(ir,ip,fi)));		 
	}
      }
    }      
   */  
   
  }
void dvar_len_fish_stock_history::predicted_frequency_calc(void)
{
  tprob.initialize();
  int ir;
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
     {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
        {                                         // incidents for this period
          for (int j=1;j<=nage;j++)         // Loop over age classes
          {
             for (int i=1; i<=nlint; i++) /* L1200  */
             {
                dvariable temp=(fmid(i)-mean_length(ir,ip,fi,j))/sdevs(ir,ip,fi,j);
                dvariable fdiff2=temp*temp/2.e0;
                dvariable prob = exp(-fdiff2) / sdevs(ir,ip,fi,j);
                tprob(ir,ip,fi,i) += prop(ir,ip,fi,j) * prob;
             }
          }
        }
     }
  }
  for (ir=1;ir<=num_regions;ir++)
  {
     for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
     {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
        {                                         // incidents for this period
          tprob(ir,ip,fi) = tprob(ir,ip,fi) / sum(tprob(ir,ip,fi));
        }
     }
  }
}
void dvar_len_fish_stock_history::variance_calc(void)
{
  dvariable rho=exp(-vb_coff(3));
  dvariable temp2=1.-pow(rho,nage-1);
  dvariable xn1inv=1.0/temp2;
  if (!pmsd)
  {
    dvar_matrix u(1,nage,1,12);
    //dvariable sd=var_coff(1);
    MY_DOUBLE_TYPE v;
    for (int tmonth=1;tmonth<=12;tmonth++)
    {
      if (parest_flags(226)==0)   //NMD_5Oct2018
      {
        for (int j=1;j<=nage;j++)         // Loop over age classes
        {
          if (parest_flags(34) == 0)
          {
            dvariable tmp=(1.-pow(rho,j-1.e0+(tmonth-1)/12.))*xn1inv;
            u(j,tmonth)=(-1.0+2.0*tmp);
          }
          else if (parest_flags(34) == 1)
          {
            u(j,tmonth)= ( 2.*(j-1.e0+(tmonth-1)/12.)-(nage-1.0))/(nage-1.e0);
          }
          else if (parest_flags(34) == 2)
          {
            cout << "WARNING: incorrect scaled mean-length calculated!!!" << endl;
            v= j-1.e0+(tmonth-1)/12.;
            u(j,tmonth)= -1.e0+2.e0*(1.0-pow(vb_coff(3),v))/
                   (1.0-pow(vb_coff(3),nage-1.0));
          }
          else
          {
            cerr << "ERROR: incorrect option for parest_flags(34)" << endl;
            ad_exit(1);
          }
        }
      }
      else
      {
        dvariable T;
        if (parest_flags(226)==1)
          T=exp(vb_coff(4));   // change sign to make simpler
        else
          T=-exp(vb_coff(4));
        dvariable c1=pow(vb_coff(1),1.0/T);
        dvariable cN=pow(vb_coff(2),1.0/T);
        dvariable diff=cN-c1;
        dvariable rho=exp(-vb_coff(3));
        dvariable temp2=1.-pow(rho,nage-1);
        dvariable xn1inv=1.0/temp2;
        for (int j=1;j<=nage;j++)
        {
          v= j-1.e0+(tmonth-1)/12.;
          dvariable tmp=(1.-pow(rho,v))*xn1inv;
          dvariable tt=c1+diff*tmp;
          u(j,tmonth)=pow(tt,T);
        }
      }   //NMD_5Oct2018
    }
    dvariable tt;
    //dvariable v1=var_coff(1)*(vb_coff(2)-vb_coff(1));
    dvariable v1=var_coff(1);
    dvariable v2=var_coff(2);
    //int set_global_vars_flag=0;
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
      {
        if (parest_flags(226)==0)    //NMD_5Oct2018
        {
          for (int j=1;j<=nage;j++)         // Loop over age classes
          {
            if (parest_flags(35)==0)  //log-linear relationship
            {
              tt=v1*exp(v2*u(j,month(ir,ip)));
            }
            else if (parest_flags(35)==1)   //normal linear relationship
            {
              tt=v1+(v2*u(j,month(ir,ip)));
            }
            else
            {
              cerr << "ERROR: incorrect option for parest_flags(35)" << endl;
              ad_exit(1);
            }
            for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
            {                                         // incidents for this period
              sdevs(ir,ip,fi,j)=tt;
              //vars(ir,ip,fi,j)=square(tt);
            }
          }
        }
        else
        {
          dvar_vector scaled_ml=setm11(column(u,month(ir,ip)));   //scale the ml between -1 and 1
          for (int j=1;j<=nage;j++)         // Loop over age classes
          {
            tt=v1*exp(v2*scaled_ml(j));
            for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
            {                                         // incidents for this period
              sdevs(ir,ip,fi,j)=tt;
              //vars(ir,ip,fi,j)=square(tt);
            }
          }
        }    //NMD_5Oct2018
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
        {                                         // incidents for this period
          vars(ir,ip,fi)=square(sdevs(ir,ip,fi));
          if (!set_global_vars_flag)
          {
            set_global_vars_flag=1;
            global_vars=vars(ir,ip,fi);
          }
        }
      }
    }
    char ch;
  }
  else
  {
    int ns=pmsd->num_species;
//    dvar3_array u(1,ns,1,nage,1,12);
    dvar3_array u(1,ns);
    //dvariable sd=var_coff(1);
    MY_DOUBLE_TYPE v;
    
    for (int is=1;is<=ns;is++)
    {
      int ng=get_nage_species(is);       //NMD_10Oct2018
      u(is).allocate(1,ng,1,12);          // dimension the partially allocated array
      u(is).initialize();
      dvar_vector vbc=get_vb_coff(is);
//      dvariable vb3=get_vb_coff(is)(3);
      for (int tmonth=1;tmonth<=12;tmonth++)
      {
        if (parest_flags(226)==0)       //NMD_10Oct2018
        {
          for (int j=1;j<=ng;j++)         // Loop over age classes
          {
//NMD_26oct2023
            if (parest_flags(34) == 0)
            {
              dvariable tmp=(1.-pow(rho,j-1.e0+(tmonth-1)/12.))*xn1inv;
              u(is,j,tmonth)=(-1.0+2.0*tmp);
            }
            else if (parest_flags(34) == 1)
            {
              u(is,j,tmonth)= ( 2.*(j-1.e0+(tmonth-1)/12.)-(ng-1.0))/(ng-1.e0);
            }
            else if (parest_flags(34) == 2)
            {
              cout << "WARNING: incorrect scaled mean-length calculated!!!" << endl;
              v= j-1.e0+(tmonth-1)/12.;
              u(is,j,tmonth)= -1.e0+2.e0*(1.0-pow(vb_coff(3),v))/
                     (1.0-pow(vb_coff(3),ng-1.0));
            }
            else
            {
              cerr << "ERROR: incorrect option for parest_flags(34)" << endl;
              ad_exit(1);
            }
//NMD_26oct2023
          }
        }             //NMD_10Oct2018
        else
        {
          dvariable T;
          if (parest_flags(226)==1)
            T=exp(vbc(4));   // change sign to make simpler
          else
            T=-exp(vbc(4));
          dvariable c1=pow(vbc(1),1.0/T);
          dvariable cN=pow(vbc(2),1.0/T);
          dvariable diff=cN-c1;
          dvariable rho=exp(-vbc(3));
          dvariable temp2=1.-pow(rho,ng-1);
          dvariable xn1inv=1.0/temp2;
          for (int j=1;j<=ng;j++)       //NMD_10Oct2018
          {
            v= j-1.e0+(tmonth-1)/12.;
            dvariable tmp=(1.-pow(rho,v))*xn1inv;
            dvariable tt=c1+diff*tmp;
            u(is,j,tmonth)=pow(tt,T);
          }
        }             //NMD_10Oct2018
      }
    }
    dvariable tt;
    //int set_global_vars_flag=0;
    for (int ir=1;ir<=num_regions;ir++)
    {
      int is=pmsd->region_species_pointer(ir);
      int ng=get_nage_species(is);       //NMD_10Oct2018
      dvariable v1=get_var_coff_species(is)(1);
      dvariable v2=get_var_coff_species(is)(2);
      for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
      {
        if (parest_flags(226)==0)             //NMD_10Oct2018
        {
          for (int j=1;j<=ng;j++)         // Loop over age classes       //NMD_10Oct2018
          {
            tt=v1*exp(v2*u(is,j,month(ir,ip)));
            for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
            {                           // incidents for this period
              sdevs(ir,ip,fi,j)=tt;
              //vars(ir,ip,fi,j)=square(tt);
            }
          }
        }             //NMD_10Oct2018
        else
        {
          dvar_vector scaled_ml=setm11(column(u(is),month(ir,ip)));   //scale the ml between -1 and 1
          for (int j=1;j<=ng;j++)         // Loop over age classes       //NMD_10Oct2018
          {
            tt=v1*exp(v2*scaled_ml(j));
            for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
            {                                         // incidents for this period
              sdevs(ir,ip,fi,j)=tt;
              //vars(ir,ip,fi,j)=square(tt);
            }
          }
        }
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
        {                                         // incidents for this period
          vars(ir,ip,fi)=square(sdevs(ir,ip,fi));
          if (!set_global_vars_flag)
          {
            set_global_vars_flag=1;
//            global_vars=vars(ir,ip,fi);
            pmsd->global_vars(is)=vars(ir,ip,fi);   //NMD_2Oct2018
          }
        }
      }
    }
    char ch;
  }
}


void length_sel_calc(dvar_len_fish_stock_history& fsh)
{
  dvar_matrix esel=exp(fsh.selcoff);
  for (int i=1;i<=fsh.num_fisheries;i++)  // Loop over fisheries
  {
	 int ff3=fsh.fish_flags(i,3);
	 for (int ii=1;ii<=fsh.nlint;ii++)
	 {
		dvariable al=age_at_length_calc(fsh.fmid(i),fsh.vb_coff,fsh.nage);
		int il=static_cast<int>(value(al));
		dvariable b=daves_kludge(al-il);
		fsh.length_sel(i,ii)=((1-b)*esel(i,min(ff3,ii))+b*esel(i,min(ff3,il+1)));
	 }
  }
}

  void dvar_len_fish_stock_history::mean_length_calc(void)
  {
    int ir;
     MY_DOUBLE_TYPE tpi=2.e0*3.14159e0;
     dvariable tmppow;
    dvariable rho=exp(-vb_coff(3));
    dvariable temp1,temp2,tmp;
    dvar_vector rel_strength(1,nyears);
    rel_strength.initialize();
    int numyrs= actual_recruit.indexmax();
    if (parest_flags(157)>0)
    {
      if (age_flags(92) || af170q0)
      {
        cerr << "Incompatible flags af92 and pf157 or af170q0" << endl;
        ad_exit(1);
      }
      for (int i=1;i<=nyears;i++)
      {
        for (ir=1;ir<=num_regions;ir++)
        {
          rel_strength(i)+=exp(N(ir,i,1));
        }
      }
      rel_strength=rel_strength-mean(rel_strength);
      rel_strength/=sqrt((.001+norm2(rel_strength))/rel_strength.size());
    }
    temp2=1.-pow(rho,nage-1.0);
    for (ir=1;ir<=num_regions;ir++)
    {
      for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
      {
        MY_DOUBLE_TYPE xx=-1+(month(ir,ip)-1.)/12.+(week(ir,ip)-1.)/48.;
        for (int j=1;j<=nage;j++)         // Loop over age classes
        {
          tmppow=j+xx;
          if (parest_flags(21))
          {
            tmppow+=sv(1)/tpi*sin(tpi*(true_month(ir,ip)/12.-sv(2)));
          }
          if (parest_flags(157)>0 || parest_flags(163)>0)
          {
            int cohort=year(ir,ip)-j+1;
            tmp=0.0;
                if (parest_flags(157)>0)
            {
              if (cohort>0)
              {
                tmp=sv(3)*rel_strength(cohort);
                if (parest_flags(158)>0 && cohort >1)
                {
                  tmp+=sv(4)*rel_strength(cohort-1);
                }
                if (parest_flags(159)>0 && cohort < numyrs)
                {
                  tmp+=sv(5)*rel_strength(cohort+1);
                }
              }
            }
            if (parest_flags(163)>0 )
            {
              if (cohort>0) tmp+=growth_dev(cohort);
            }
            tmp=1./(1.+exp(-tmp));
            tmp=1.9*(tmp-.5);
            tmppow+=tmp;
          }
          if (parest_flags(160)>0 && j==nage)
          {
            tmppow += sv(6);
          }
          if (parest_flags(168)) 
          {
            if (!parest_flags(169)) 
              tmppow=pow(tmppow/nage,sv(16))*nage;
            else if (!parest_flags(170))
            {
              dvariable tmp=tmppow/nage;
              tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5))*nage;
            }
            else if (!parest_flags(171))
            {
              dvariable tmp=tmppow/nage;
              tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5)
                 +sv(18)*square(tmp-0.5))*nage;
            }
            else if (!parest_flags(172))
            {
              dvariable tmp=tmppow/nage;
              tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5)
                 +sv(18)*square(tmp-0.5)
                 +sv(19)*cube(tmp-0.5))*nage;
            }
            else
            {
              dvariable tmp=tmppow/nage;
              tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5)
                 +sv(18)*square(tmp-0.5)
                 +sv(20)*square(square(tmp-0.5))
                 +sv(19)*cube(tmp-0.5))*nage;
            }
          }
          if (parest_flags(174))
          {
            dvariable temp4=sv(19)/(sv(18)*2.506628275)*
                         exp(-0.5*square((j-sv(17))/sv(18)));
            rho=exp(-(vb_coff(3)-temp4));
            temp2=1.-pow(rho,nage-1.0);
          }
          
          temp1=1.-pow(rho,tmppow);

          dvariable tt;
          length_calc(tt,vb_coff(1),vb_coff(2),temp1,temp2);
          
          for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
          {
            mean_length(ir,ip,fi,j)=tt;
          }
        }
      }
    }
    if (parest_flags(173))
    {
      pmsd_error();
      int num=parest_flags(173);
      for (ir=1;ir<=num_regions;ir++) {
        for (int ip=1;ip<=num_fish_periods(ir);ip++) {
          for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) { 
            for (int j=1;j<num;j++) {
              mean_length(ir,ip,fi,j+1)+= age_pars(3,j);
	    }
          }
        }
      }
    }
    for (ir=1;ir<=num_regions;ir++)
     {
      for (int ip=1;ip<=num_fish_periods(ir);ip++) 
      {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
        {                                         // incidents for this period
          if (fish_flags(parent(ir,ip,fi),11) > 0)
          {
            mean_length(ir,ip,fi,1)+=vb_bias(parent(ir,ip,fi))*
              (12.-month(ir,ip))/12.;
          }
        }
      }
    }
  }

void length_calc(dvariable& tt,const prevariable& vb1,
  const prevariable& vb2,const prevariable& t1,const prevariable& t2)
{
  
  MY_DOUBLE_TYPE diff=value(vb2)-value(vb1);
  MY_DOUBLE_TYPE quot=value(t1)/value(t2);
  value(tt)=value(vb1)+diff*quot;
  gradient_structure::GRAD_STACK1->set_gradient_stack(default_evaluation4ind,
    address(tt),address(vb1),1-quot,address(vb2),
    quot,address(t1),diff/value(t2),address(t2),-diff*quot/value(t2));
  
}
MY_DOUBLE_TYPE davetest(dvar_vector& x)
{
  dvariable tt;
  dvariable u=0.0;
  for (int i=1;i<=100;i++)
  {
    for (int j=1;j<=10;j++)
    {
      length_calc(tt,x(1+i),x(2+i),x(3+i),x(4+i));
      u+=square(tt-100.0);
    }
  }
  return value(u);
}

  dvector dvar_len_fish_stock_history::vb_length_calc(int ir,int ip,int fi)
  {
     MY_DOUBLE_TYPE tpi=2.e0*3.14159e0;
    MY_DOUBLE_TYPE tmppow;
    MY_DOUBLE_TYPE rho=exp(-value(vb_coff(3)));
    MY_DOUBLE_TYPE temp1,temp2,tmp;

    dvector vb_lengths(1,nage);
    temp2=1.-pow(value(rho),nage-1.0);
    MY_DOUBLE_TYPE xx=-1+(month(ir,ip)-1.)/12.+(week(ir,ip)-1.)/48.;
    for (int j=1;j<=nage;j++)         // Loop over age classes
    {
      tmppow=j+xx;
      if (parest_flags(21))
      {
        tmppow+=value(sv(1))/tpi*sin(tpi*(true_month(ir,ip)/12.-value(sv(2))));
      }
      temp1=1.-pow(value(rho),tmppow);
      if (parest_flags(168)) tmppow=pow(tmppow/nage,value(sv(16)))*nage;
      vb_lengths(j)=value(vb_coff(1))+(value(vb_coff(2))-value(vb_coff(1)))*
        temp1/temp2;
    }
    if (fish_flags(parent(ir,ip,fi),11) > 0)
    {
      vb_lengths(1)+=value(vb_bias(parent(ir,ip,fi)))*(12.-month(ir,ip))/12.;
    }
    return vb_lengths;
  }

  void dvar_len_fish_stock_history::mean_length_calcx(void)
  {
    dvar_vector& lcoffs=age_pars(4);
    dvar_vector mls(1,nage);
    mls(1)=lcoffs(1);
    for (int j=1;j<nage;j++)         // Loop over age classes
      mls(j+1)=mls(j)+lcoffs(j+1);

    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
      {
        MY_DOUBLE_TYPE xx=(month(ir,ip)-1.)/12.+(week(ir,ip)-1.)/48.;
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {
          for (int j=1;j<nage;j++)         // Loop over age classes
          {
            mean_length(ir,ip,fi,j)=mls(j)+xx*lcoffs(j+1);
          }
          mean_length(ir,ip,fi,nage)=mls(nage)+xx*sv(20);
        }
      }
    }
  }

  void setup_gml(int nage,dvar_vector& lcoffs,dvar_vector& gml,const prevariable& sv20)
  {  
    gml(1)=lcoffs(1);
    for (int j=1;j<nage;j++)         // Loop over age classes
      gml(j+1)=gml(j)+lcoffs(j+1);
      gml(nage+1)=gml(nage)+sv20;
  }

  void setup_gmlxx(int nage,dvar_vector& vb_coff,dvar_vector& gml,
    dvar_vector& sv)
  {
    dvariable tmppow;
    dvariable rho=exp(-vb_coff(3));
    dvariable temp2=1.-pow(rho,nage-1.0);
    for (int j=1;j<=nage+1;j++)         // Loop over age classes
    {
      tmppow=j;
      dvariable temp4=sv(19)/(sv(18)*2.506628275)*
        exp(-0.5*square((j-sv(17))/sv(18)));
        rho=exp(-(vb_coff(3)-temp4));
        temp2=1.-pow(rho,nage-1.0);
      dvariable temp1=1.-pow(rho,tmppow);
      dvariable tt;
      length_calc(tt,vb_coff(1),vb_coff(2),temp1,temp2);
      gml(j)=tt;
    }
  }

  void dvar_fish_stock_history::set_catchabilities_to_zero(ivector * pq_flag)
  {
    zeroed_fisheries_flag=0;
    dvar3_array* pc=0;
    if (af170q0)
    {
      pc=&catchability_q0;
    }
    else
    {
      pc=&catchability;
    }
    if (pq_flag)
    {
      for (int fi=1;fi<=num_fisheries;fi++)
      {
        if ((*pq_flag)(fi)) 
        {
          zeroed_fisheries_flag=1;
          for (int i=1;i<=num_fish_times(fi);i++)
          {
            int rr=realization_region(fi,i);
            int rp=realization_period(fi,i);
            int ri=realization_incident(fi,i);
            (*pc)(rr,rp,ri)=-50;
          }
        }
      }
    }
  }    
 
#undef HOME_VERSION
