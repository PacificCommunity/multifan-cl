/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

dvector simulate_multinomial(const dvector & truep,int n,
  const random_number_generator& rng)
 {
   int mmin=truep.indexmin();
   int mmax=truep.indexmax();
   ivector musamp(1,n);
   musamp.fill_multinomial(rng,truep/(1.e-10+sum(truep)));
   dvector hist(mmin,mmax);
   hist.initialize();
   for (int ii=1;ii<=n;ii++)
   {
     hist(musamp(ii))+=1;
   }
   return hist;
 }


void dvar_len_fish_stock_history::simulate_length_and_weight_frequencies(void)
{
  int iseed1=8001;
  if (simulation_seeds(1))
    iseed1=simulation_seeds(1);
  random_number_generator rng1(iseed1+2*projection_sim_index);
  ofstream ofs("test_lw_sim",ios::app);
  ofstream ofs1("effort_sim",ios::app);

  std::unique_ptr<std::ofstream> pofs;
  pofs=std::unique_ptr<std::ofstream>(new std::ofstream("effort_sim_true", 
       std::ios::app));
  
  ofstream ofs2("catch_sim",ios::app);
  ofs << "# Simulated length frequencies" << endl;
  ofs << "# projection " << projection_sim_index << endl;
  ofs << "# seed " << iseed1+2*projection_sim_index << endl;
  ofstream ofsalt("test_lw_sim_alt",ios::app); // FS Jan 19
  ofsalt << "# Simulated length frequencies" << endl;
  ofsalt << "# projection seed year month week fishery sum(lensamp) data" << endl;
  ivector start_period(1,num_regions);
  start_period=1;
  for (int ir=1;ir<=num_regions;ir++)
  {
    if (parest_flags(242)==0)
    {
      start_period(ir)=num_real_fish_periods(ir)+1;
    }
  }

  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=start_period(ir);ip<=num_fish_periods(ir);ip++)
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)
      {
        int ty=really_true_year(ir,ip)+year1-1;
        if (len_sample_size(ir,ip,fi)>0)
        {
          dvector lensamp=simulate_multinomial(value(tprob(ir,ip,fi)), 
            static_cast<int>(len_sample_size(ir,ip,fi)),rng1);
          ofs << ty << "  " << really_true_month(ir,ip) <<  "  " << true_week(ir,ip) << " " <<  parent(ir,ip,fi) << endl;
          ofs << len_sample_size(ir,ip,fi) << endl;
          ofsalt << projection_sim_index << " " << iseed1+2*projection_sim_index << " " << ty << " " <<
            really_true_month(ir,ip) <<" " << true_week(ir,ip) << " " <<  parent(ir,ip,fi) << " " <<
            sum(lensamp); // FS Jan 19
          if (sum(lensamp)>0.001)     //NMD_13Sep2018
          {
            ofs << lensamp << endl;
            ofsalt << lensamp << endl;
          }         //NMD_13Sep2018
          else // If sum of samples is too low then just print empty vector
          {
            lensamp = 0;
            ofs << lensamp << endl;
            ofsalt << lensamp << endl;
          }
        }
      }
    }
  }
  int iseed2=8003;
  if (simulation_seeds(2))
    iseed2=simulation_seeds(2);
  random_number_generator rng2(iseed2+2*projection_sim_index);
  ofs << "# Simulated weight frequencies" << endl;
  ofs << "# projection " << projection_sim_index << endl;
  ofs << "# seed " << iseed2+2*projection_sim_index << endl;
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=start_period(ir);ip<=num_fish_periods(ir);ip++)
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)
      {
        int ty=really_true_year(ir,ip)+year1-1;
        if (wght_sample_size(ir,ip,fi)>0)
        {
          dvector wtsamp=simulate_multinomial(value(wtprob(ir,ip,fi)),
            static_cast<int>(wght_sample_size(ir,ip,fi)),rng2);
          if (sum(wtsamp)>0.001)
          {
            ofs << ty << "  " << really_true_month(ir,ip) <<  "  "
                << true_week(ir,ip) << " " <<  parent(ir,ip,fi) << endl;
            ofs << wght_sample_size(ir,ip,fi) << endl;
            ofs << wtsamp << endl;
          }
        }
      }
    }
  }

  ofs1 << "# projection " << projection_sim_index << endl;
  ofs2 << "# projection " << projection_sim_index << endl;
  if (parest_flags(244))
  {
    *pofs << "# projection " << projection_sim_index
                              << endl;
  }

  int iseed3=14003;
  int iseed4=16003;
  if (simulation_seeds(3))
    iseed3=simulation_seeds(3);
  if (simulation_seeds(4))
    iseed4=simulation_seeds(4);
  int the_seed=iseed3+2*projection_sim_index;
  int the_seed4=iseed4+2*projection_sim_index;
  ofs1 << "# seed " << the_seed << endl;
  ofs2 << "# seed " << the_seed4 << endl;
  if (parest_flags(244))
  {
    *pofs << "# seed " << the_seed  << endl;
  }

  random_number_generator rng3(the_seed);
  random_number_generator rng4(the_seed4);
  for (int i=1;i<=num_fisheries;i++)
  {
    int label_switch=0;
    for (int nt=1;nt<=num_fish_times(i);nt++)
    {
      int ir=realization_region(i,nt);
      int ip=realization_period(i,nt);
      int fi=realization_incident(i,nt);  
      if (ip>=start_period(ir))
      {
        if (label_switch==0)
        {
          ofs1 << setw(5) << i;
          ofs2 << "# fishery " << i << endl;
          if (parest_flags(244))
          {
            *pofs << setw(5) << i ;
          }
          label_switch=1;
        }
        MY_DOUBLE_TYPE eps=0.0;
        MY_DOUBLE_TYPE eps1=0.0;
        if (age_flags(186)>0)
        {
          MY_DOUBLE_TYPE sd=age_flags(186)/100.;
          MY_DOUBLE_TYPE bias=0.5*sd*sd;
          eps=sd*randn(rng3)-bias;
        }
        if (age_flags(187)>0)
        {
          MY_DOUBLE_TYPE sd=age_flags(187)/100.;
          MY_DOUBLE_TYPE bias=0.5*sd*sd;
          eps1=sd*randn(rng4)-bias;
        }
        if (fm_level(ir,ip,fi)>1.e-20)
        {
          if (!age_flags(92))
          {
            ofs1 << setscientific() << setprecision(4) << " "
                 << fm_level(ir,ip,fi)/exp(catchability(ir,ip,fi)) 
                   *effort_normalization_factor(ir,ip,fi)
                   *exp(eps);
            if (parest_flags(244))
            {
              *pofs << setscientific() << setprecision(4) << " "
                   << fm_level(ir,ip,fi)/exp(catchability(i,ip,fi)) 
                     *effort_normalization_factor(ir,ip,fi);
            }
          }
          else
          {
            if (eff_proj_fshry(i) || catch_proj_fshry(i))
            {
              ofs1 << setscientific() << setprecision(4) << " "
                   << fm_level(ir,ip,fi)
                     /exp(implicit_catchability_ccond(i,nt)) 
                     *effort_normalization_factor(ir,ip,fi)
                     *exp(eps);
              if (parest_flags(244))
              {
                *pofs << setscientific() << setprecision(4) << " "
                     << fm_level(ir,ip,fi)
                       /exp(implicit_catchability_ccond(i,nt)) 
                       *effort_normalization_factor(ir,ip,fi);
              }
            }

          }
        }
        else
        {
          ofs1 << setscientific() << setprecision(4) << " "
               <<  exp(effort(ir,ip,fi)+effort_devs(ir,ip,fi))*
                   effort_normalization_factor(ir,ip,fi)*exp(eps);
          if (parest_flags(244))
          {
            *pofs << setscientific() << setprecision(4) << " "
                  <<  exp(effort(ir,ip,fi)+effort_devs(ir,ip,fi))*
                      effort_normalization_factor(ir,ip,fi);
          }


        }
        if (!data_fish_flags(1,parent(ir,ip,fi)))
        {
          ofs2 << setscientific() << setprecision(4) << " "
               << tot_catch(ir,ip,fi) << " "
               << sum(exp_catch(ir,ip,fi))*exp(eps1) 
               << endl;
        }
        else
        {
          dvariable ecatch=0.0;
          dvariable sv27=get_sv_region(ir,27);
          dvariable sv28=get_sv_region(ir,28);
          if (value(sv28)==3.0)
          {
            ecatch=len_wt_coff/1000.* (exp_catch(ir,ip,fi)*
              (pow(mean_length(ir,ip,fi),3)+
	       3.0*elem_prod(mean_length(ir,ip,fi),vars(ir,ip,fi))));
          }
          else
          {
            for (int j=1;j<=nage;j++)
            {  
              ecatch+=exp_catch(ir,ip,fi,j)*
                normal_length_to_weight(0.5,-3.5,
                mean_length(ir,ip,fi,j),sqrt(vars(ir,ip,fi,j)),
                value(sv27),value(sv28));
            }
            ecatch/=1000.;
          }
          if (ecatch<=0.0) 
          {
            cout << "negative exp_catch for ir " << ir << " ip " << ip 
                 << "  fi " << fi << endl;
          }
          ofs2 << setscientific() << setprecision(4) << " "
               << ecatch << " " << ecatch*exp(eps1) << endl;
        }
      }
    }
    if (label_switch==1)
    {
      ofs1 << endl;
      if (parest_flags(244))
      {
        *pofs << endl;
      }
    }
  }
}

void dvar_len_fish_stock_history::generate_simulated_cpue(void)
{
  ivector ff92(1,num_fisheries);
  ff92=column(fish_flags,92);

  ivector start_time(1,num_fisheries);
  start_time=1;
  for (int fi=1;fi<=num_fisheries;fi++)
  {
    if (parest_flags(242)==0)
    {
      start_time(fi)=num_real_fish_times(fi)+1;
    }
  }

  ofstream ofs("cpue_sim",ios::app);
  std::unique_ptr<std::ofstream> pofs;
  pofs=std::unique_ptr<std::ofstream>(new std::ofstream("cpue_sim_true", 
       std::ios::app));

  dvector grpd_mls(1,num_fisheries);
  dvar_vector grpd_mlv(1,num_fisheries);
  grpd_mls.initialize();
  grpd_mlv.initialize();
  ivector nft;
  ivector nrft;
  nft=num_fish_times;
  nrft=num_real_fish_times;

  if (allocated(fishery_group_ptr) &&  allocated(fishery_group_ptr(99)))
  {
// - derive mean over estimation periods only
    int num_groups=fishery_group_ptr(99).indexmax();
    for (int ig=1;ig<=num_groups;ig++)
    {
      if (ff92(fishery_group_ptr(99,ig,1)))
      {
        int nobs=0;
        int nobs1=0;
        ivector fgp99=fishery_group_ptr(99,ig);
        for (int i=1;i<=fgp99.indexmax();i++)
        {
          int fi=fgp99(i);
          nobs+= nrft(fi)-first_survey_time(fi)+1;
          nobs1+=survey_index(fi).indexmax()-survey_index(fi).indexmin()+1;
        }
        if (nobs != nobs1)
        {
          cerr << " size error for grouped survey data vs cpue" << endl;
          ad_exit(1);
        }
        if (data_fish_flags(1,fgp99(1))==0)
        {
          dvar_vector grouped_vulnerable_numbers(1,nobs);
          dvector grouped_survey_index(1,nobs);
          grouped_survey_index.initialize();
          grouped_vulnerable_numbers.initialize();
          int offset=0;
          for (int ii=1;ii<=fgp99.indexmax();ii++)
          {
            int fi=fgp99(ii);
            for (int i=first_survey_time(fi);i<=nrft(fi);i++)
            {
              int i1=i+offset-first_survey_time(fi)+1;
              grouped_vulnerable_numbers(i1)=sum(vul_at_age(fi,i));  
              grouped_survey_index(i1)=survey_index(fi,i);
            }
            offset+=nrft(fi)-first_survey_time(fi)+1;
          }
          dvector log_groupedsurv=log(1.e-10+grouped_survey_index);
          dvar_vector log_groupedvuln=log(1.e-10+grouped_vulnerable_numbers);
          double mls=mean(log_groupedsurv);
          dvariable mlv=mean(log_groupedvuln);
          for (int ii=1;ii<=fgp99.indexmax();ii++)
          {
            int fi=fgp99(ii);
            grpd_mls(fi)=mls;
            grpd_mlv(fi)=mlv;
          }
        }
        else
        {
          dvar_vector grouped_vulnerable_biomass(1,nobs);
          dvector grouped_survey_index(1,nobs);
          grouped_survey_index.initialize();
          grouped_vulnerable_biomass.initialize();
          int offset=0;
          for (int ii=1;ii<=fgp99.indexmax();ii++)
          {
            int fi=fgp99(ii);
            for (int i=first_survey_time(fi);i<=nrft(fi);i++)
            {
              int i1=i+offset-first_survey_time(fi)+1;
              int rr=realization_region(fi,i); 
              int rp=realization_period(fi,i); 
              int ri=realization_incident(fi,i); 
              grouped_vulnerable_biomass(i1)=mean_weight(rr,rp,ri)*
                vul_at_age(fi,i);
              grouped_survey_index(i1)=survey_index(fi,i);
            }
            offset+=nrft(fi)-first_survey_time(fi)+1;
          }
          dvector log_groupedsurv=log(1.e-10+grouped_survey_index);
          dvar_vector log_groupedvuln=log(1.e-10+grouped_vulnerable_biomass);
          double mls=mean(log_groupedsurv);
          dvariable mlv=mean(log_groupedvuln);
          for (int ii=1;ii<=fgp99.indexmax();ii++)
          {
            int fi=fgp99(ii);
            grpd_mls(fi)=mls;
            grpd_mlv(fi)=mlv;
          }
        }
      }
    }
  }

  int iseed5=19003;
  if (simulation_seeds(6))
    iseed5=simulation_seeds(6);
  int the_seed5=iseed5+2*projection_sim_index;
  random_number_generator rng4(the_seed5);
  ofs << "# Simulated CPUE indices" << endl;
  ofs << "# projection " << projection_sim_index << endl;
  ofs << "# seed " << the_seed5 << endl;

  if (parest_flags(244))
  {
    *pofs << "# Simulated CPUE indices - no error" << endl;
    *pofs << "# projection " << projection_sim_index << endl;
    *pofs << "# seed " << the_seed5 << endl;
  }

  MY_DOUBLE_TYPE sd=0.0;
  MY_DOUBLE_TYPE bias=0.0;
  if (age_flags(26)>0)
  {
    sd=age_flags(26)/100.;
    bias=0.5*sd*sd;
  }
  MY_DOUBLE_TYPE eps1=0.0;

  for (int fi=1;fi<=num_fisheries;fi++)
  {
    if (ff92(fi))
    {
      dvector cpue_obs(first_survey_time(fi),nft(fi));
      dvector cpue_pred(first_survey_time(fi),nft(fi));

      if (data_fish_flags(1,fi)==0)
      {
// - derive mean over estimation periods only
        dvar_vector vulnerable_numbers(first_survey_time(fi),nrft(fi));
        vulnerable_numbers.initialize();
        for (int i=first_survey_time(fi);i<=nrft(fi);i++)
        {
          vulnerable_numbers(i)=sum(vul_at_age(fi,i));  
        }
        dvector logsurv=log(1.e-10+survey_index(fi));
        dvar_vector logvuln=log(1.e-10+vulnerable_numbers);
        double mls;
        dvariable mlv;
        if (allocated(fishery_group_ptr) &&  allocated(fishery_group_ptr(99)))
        {
          mls=grpd_mls(fi);
          mlv=grpd_mlv(fi);
        }
        else
        {
          mls=mean(logsurv);
          mlv=mean(logvuln);
        }
// scale the full time periods by the means
        dvar_vector vulnerable_numbers_projn(first_survey_time(fi),nft(fi));
        vulnerable_numbers_projn.initialize();
        for (int i=first_survey_time(fi);i<=nft(fi);i++)
        {
          vulnerable_numbers_projn(i)=sum(vul_at_age(fi,i));  
        }
        dvar_vector logvuln_projn=log(1.e-10+vulnerable_numbers_projn);
        cpue_pred=value(logvuln_projn-mlv+mls);
      }
      else
      {
// - derive mean over estimation periods only
        dvar_vector  vulnerable_biomass(first_survey_time(fi),nrft(fi));
        vulnerable_biomass.initialize();
        for (int i=first_survey_time(fi);i<=nrft(fi);i++)
        {
          int rr=realization_region(fi,i); 
          int rp=realization_period(fi,i); 
          int ri=realization_incident(fi,i); 
          vulnerable_biomass(i)=mean_weight(rr,rp,ri)
            *vul_at_age(fi,i);
        }
        dvector logsurv=log(1.e-10+survey_index(fi));
        dvar_vector logvuln=log(1.e-10+vulnerable_biomass);
        double mls=mean(logsurv);
        dvariable mlv=mean(logvuln);
        if (allocated(fishery_group_ptr) &&  allocated(fishery_group_ptr(99)))
        {
          mls=grpd_mls(fi);
          mlv=grpd_mlv(fi);
        }
        else
        {
          mls=mean(logsurv);
          mlv=mean(logvuln);
        }
// scale the full time periods by the means
        dvar_vector vulnerable_biomass_projn(first_survey_time(fi),nft(fi));
        vulnerable_biomass_projn.initialize();
        for (int i=first_survey_time(fi);i<=nft(fi);i++)
        {
          int rr=realization_region(fi,i); 
          int rp=realization_period(fi,i); 
          int ri=realization_incident(fi,i); 
          vulnerable_biomass_projn(i)=mean_weight(rr,rp,ri)
            *vul_at_age(fi,i);
        }
        dvar_vector logvuln_projn=log(1.e-10+vulnerable_biomass_projn);
        cpue_pred=value(logvuln_projn-mlv+mls);
      }
// express predictions in normal space
      dvector ecpue_pred(1,nft(fi));
      ecpue_pred=exp(cpue_pred);
// add pseudo-observation error
      dvector ecpue_pred_err(1,nft(fi));
      eps1=0.0;
      if (age_flags(26)>0){
        for (int i=start_time(fi);i<=nft(fi);i++)
        {
          eps1=sd*randn(rng4)-bias;
          ecpue_pred_err(i)=ecpue_pred(i)*exp(eps1);
        }
      }

      if (!af170q0)
      {
        ofs << " Fishery: " << fi << " " << setfixed() << setprecision(4)
	    << setw(7) << ecpue_pred_err << endl;
        if (parest_flags(244))
        {
          *pofs << " Fishery: " << fi << " " << setfixed() << setprecision(4)
                << setw(7) << ecpue_pred << endl;
        }
      }

    }
  }

}
