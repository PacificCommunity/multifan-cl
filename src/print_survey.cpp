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
#include "newmprot.hpp"

extern mf_pvm_manager * mf_pvm;
void dvar_fish_stock_history::print_survey_fishery_fits(ivector& ff92)
{
  // check if catchability/variance for surver fishery catch cpue vs 
  // survey cpue  is grouped and print grouped fit reports
  //
  dvector grpd_mls(1,num_fisheries);
  dvar_vector grpd_mlv(1,num_fisheries);
  grpd_mls.initialize();
  grpd_mlv.initialize();
  ivector local_fishery_pointer;
  ivector nft;    //NMD_28jun2022
  if (sum(data_fish_flags(2)))
  {
    nft=num_real_fish_times;
  }
  else
  {
    nft=num_fish_times;
  }

  if (allocated(fishery_group_ptr) &&  allocated(fishery_group_ptr(99)))
  {
    std::unique_ptr<std::ofstream> ofs1;
    ofs1=std::unique_ptr<std::ofstream>(new std::ofstream("grouped_cpue_obs",
      std::ios::app));
    std::unique_ptr<std::ofstream> ofs2;
    ofs2=std::unique_ptr<std::ofstream>(new std::ofstream("grouped_cpue_pred",
      std::ios::app));

    int num_groups=fishery_group_ptr(99).indexmax();
    ivector srtd_grp;
    ivector srtd_grp_sub;
    for (int ig=1;ig<=num_groups;ig++)
    {
      if (ff92(fishery_group_ptr(99,ig,1)))
      {
        int nobs=0;
        int nobs1=0;
        ivector fgp99=fishery_group_ptr(99,ig);
        if (!pmsd)
        {
          for (int i=1;i<=fgp99.indexmax();i++)
          {
            int fi=fgp99(i);
            nobs+= nft(fi)-first_survey_time(fi)+1;
            nobs1+=survey_index(fi).indexmax()-survey_index(fi).indexmin()+1;
          }
        }
        else
        {
          srtd_grp.allocate(1,fgp99.indexmax());
          srtd_grp.initialize();
          srtd_grp=sort(fgp99,1);
          int spp1_max;
          for (int ii=1;ii<=srtd_grp.indexmax();ii++)
          {
            if (srtd_grp(ii)>pmsd->num_real_fisheries)
            {
              break;
            }
            spp1_max=ii;
          }

          srtd_grp_sub.allocate(1,spp1_max);
          srtd_grp_sub=srtd_grp(1,spp1_max);
          for (int i=1;i<=srtd_grp_sub.indexmax();i++)
          {
            int fi=srtd_grp_sub(i);
            nobs+= nft(fi)-first_survey_time(fi)+1;
            nobs1+=survey_index(fi).indexmax()-survey_index(fi).indexmin()+1;
          }
        }
        if (nobs != nobs1)
        {
          cerr << " size error for grouped survey data vs cpue" << endl;
          ad_exit(1);
        }
        //
        dvector cpue_obs(1,nobs);
        dvector cpue_pred(1,nobs);
        cpue_obs.initialize();
        cpue_pred.initialize();
        if (data_fish_flags(1,fgp99(1))==0)
        {
          dvar_vector grouped_vulnerable_numbers(1,nobs);
          dvector grouped_survey_index(1,nobs);

          grouped_survey_index.initialize();
          grouped_vulnerable_numbers.initialize();
          int offset=0;
          if (!pmsd)
          {
            for (int ii=1;ii<=fgp99.indexmax();ii++)
            {
              int fi=fgp99(ii);

              for (int i=first_survey_time(fi);i<=nft(fi);i++)
              {
                int i1=i+offset-first_survey_time(fi)+1;
                grouped_vulnerable_numbers(i1)=sum(vul_at_age(fi,i));  
                grouped_survey_index(i1)=survey_index(fi,i);
              }
              offset+=nft(fi)-first_survey_time(fi)+1;
            }
	  }
          else
          {
//          Aggregated predictions over kludged fisheries in group	    
            for (int ii=1;ii<=srtd_grp_sub.indexmax();ii++)
            {
              int fi=srtd_grp_sub(ii);
              int irr=fishery_regions(fi);
              int isp=pmsd->region_species_pointer(irr);
              local_fishery_pointer.allocate(1,pmsd->num_species);
              int irr_sp=0;
              if (isp==1)
              {
                local_fishery_pointer(isp)=fi;
              }
              offset=0;
              for (int i=2;i<=pmsd->num_species;i++) 
              {
                offset=(i-1)*pmsd->num_real_fisheries;
                int irr_sp=fishery_regions(fi+offset);
                local_fishery_pointer(i)=fi+offset;
              }

// Aggregate over the kludged fisheries for species in fi
              int ng=vul_at_age(fi,1).indexmax();
              dvar_matrix tmp_spp(1,nft(fi),1,ng);
              tmp_spp.initialize();
              for (int i=first_survey_time(fi);i<=nft(fi);i++)
              {
                for (int isp=1;isp<=pmsd->num_species;isp++) 
                {
                  int fisp=local_fishery_pointer(isp);
                  tmp_spp(i)+=vul_at_age(fisp,i);
                }
              }
// Replace into spp 1 fishery in group for each member
              for (int i=first_survey_time(fi);i<=nft(fi);i++)
              {
                vul_at_age(fi,i)=tmp_spp(i);
              }
              tmp_spp.deallocate();
            }
// Concatenation all fisheries in the group using srtd_grp_sub order
// as applied for obs
            offset=0;
            for (int ii=1;ii<=srtd_grp_sub.indexmax();ii++)
            {
              int fi=srtd_grp_sub(ii);
              for (int i=first_survey_time(fi);i<=nft(fi);i++)
              {
                int i1=i+offset-first_survey_time(fi)+1;
                grouped_vulnerable_numbers(i1)=sum(vul_at_age(fi,i));  
                grouped_survey_index(i1)=survey_index(fi,i);
              }
              offset+=nft(fi)-first_survey_time(fi)+1;
            }
          }
          dvector log_groupedsurv=log(1.e-10+grouped_survey_index);
          dvar_vector log_groupedvuln=log(1.e-10+grouped_vulnerable_numbers);
          double mls=mean(log_groupedsurv);
          dvariable mlv=mean(log_groupedvuln);
          cpue_obs=log_groupedsurv-mls;
          cpue_pred=value(log_groupedvuln-mlv);
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
          if (!pmsd)
          {
            for (int ii=1;ii<=fgp99.indexmax();ii++)
            {
              int fi=fgp99(ii);
              for (int i=first_survey_time(fi);i<=nft(fi);i++)
              {
                int i1=i+offset-first_survey_time(fi)+1;
                int rr=realization_region(fi,i); 
                int rp=realization_period(fi,i); 
                int ri=realization_incident(fi,i); 
                grouped_vulnerable_biomass(i1)=mean_weight(rr,rp,ri)*
                  vul_at_age(fi,i);
                grouped_survey_index(i1)=survey_index(fi,i);
              }
              offset+=nft(fi)-first_survey_time(fi)+1;
            }
          }
          else
          {
//          Aggregated predictions over kludged fisheries in group	    
            for (int ii=1;ii<=srtd_grp_sub.indexmax();ii++)
            {
              int fi=srtd_grp_sub(ii);
              int irr=fishery_regions(fi);
              int isp=pmsd->region_species_pointer(irr);
              local_fishery_pointer.allocate(1,pmsd->num_species);
              int irr_sp=0;
              if (isp==1)
              {
                local_fishery_pointer(isp)=fi;
              }
              offset=0;
              for (int i=2;i<=pmsd->num_species;i++) 
              {
                offset=(i-1)*pmsd->num_real_fisheries;
                int irr_sp=fishery_regions(fi+offset);
                local_fishery_pointer(i)=fi+offset;
              }
// Aggregate over the kludged fisheries for species in fi
              int ng=vul_at_age(fi,1).indexmax();
              dvar_matrix tmp_spp(1,nft(fi),1,ng);
              tmp_spp.initialize();
              for (int i=first_survey_time(fi);i<=nft(fi);i++)
              {
                for (int isp=1;isp<=pmsd->num_species;isp++) 
                {
                  int fisp=local_fishery_pointer(isp);
                  tmp_spp(i)+=vul_at_age(fisp,i);
                }
              }
// Replace into spp 1 fishery in group for each member
              for (int i=first_survey_time(fi);i<=nft(fi);i++)
              {
                vul_at_age(fi,i)=tmp_spp(i);
              }
              tmp_spp.deallocate();
            }
// Concatenation of all fisheries in the group using srtd_grp_sub order
            offset=0;
            for (int ii=1;ii<=srtd_grp_sub.indexmax();ii++)
            {
              int fi=srtd_grp_sub(ii);
              for (int i=first_survey_time(fi);i<=nft(fi);i++)
              {
                int rr=realization_region(fi,i); 
                int rp=realization_period(fi,i); 
                int ri=realization_incident(fi,i); 
                int i1=i+offset-first_survey_time(fi)+1;
                grouped_vulnerable_biomass(i1)=mean_weight(rr,rp,ri)
                                               *vul_at_age(fi,i);
                grouped_survey_index(i1)=survey_index(fi,i);
              }
              offset+=nft(fi)-first_survey_time(fi)+1;
            }
          }
          dvector log_groupedsurv=log(1.e-10+grouped_survey_index);
          dvar_vector log_groupedvuln=log(1.e-10+grouped_vulnerable_biomass);
          double mls=mean(log_groupedsurv);
          dvariable mlv=mean(log_groupedvuln);
          cpue_obs=log_groupedsurv-mls;
          cpue_pred=value(log_groupedvuln-mlv);
          for (int ii=1;ii<=fgp99.indexmax();ii++)
          {
            int fi=fgp99(ii);
            grpd_mls(fi)=mls;
            grpd_mlv(fi)=mlv;
          }
        }
        if (!af170q0)
        {
          *ofs1 << " Group: " << ig << " " << setfixed() << setprecision(4)
               << setw(7) << cpue_obs << endl;    //NMD_30jul2021
          *ofs2 << " Group: " << ig << " " << setfixed() << setprecision(4)
               << setw(7) << cpue_pred << endl;   //NMD_30jul2021
        }
      }
    }
  }

  std::unique_ptr<std::ofstream> ofs1;
  ofs1=std::unique_ptr<std::ofstream>(new std::ofstream("cpue_obs",
    std::ios::app));
  std::unique_ptr<std::ofstream> ofs2;
  ofs2=std::unique_ptr<std::ofstream>(new std::ofstream("cpue_pred",
    std::ios::app));
  std::unique_ptr<std::ofstream> ofs3;
  ofs3=std::unique_ptr<std::ofstream>(new std::ofstream("cpue_obs_mls",
    std::ios::app));
  std::unique_ptr<std::ofstream> ofs4;
  ofs4=std::unique_ptr<std::ofstream>(new std::ofstream("cpue_pred_mlv",
    std::ios::app));

  for (int fi=1;fi<=num_fisheries;fi++)
  {
    int fspp1=1;
    if (ff92(fi) && pmsd)
    {
      if (fi<=pmsd->num_real_fisheries)
      {
        int irr=fishery_regions(fi);
        int isp=pmsd->region_species_pointer(irr);
        local_fishery_pointer.allocate(1,pmsd->num_species);
        if (isp==1)
        {
          local_fishery_pointer(isp)=fi;
        }
        for (int i=2;i<=pmsd->num_species;i++) 
        {
          int offset=(i-1)*pmsd->num_real_fisheries;
          local_fishery_pointer(i)=fi+offset;
        }
        fspp1=1;
      }
      else
      {
        fspp1=0;
      }
    }
    int nobs=survey_index(fi).indexmax()-survey_index(fi).indexmin()+1;
    if (ff92(fi) && fspp1==1)
    {
      dvector cpue_obs(first_survey_time(fi),nft(fi));
      dvector cpue_pred(first_survey_time(fi),nft(fi));
      double tmps=0.0;
      double tmpv=0.0;
      if (data_fish_flags(1,fi)==0)
      {  
        dvar_vector vulnerable_numbers(first_survey_time(fi),nft(fi));
        vulnerable_numbers.initialize();
        if (!pmsd)
        {
          for (int i=first_survey_time(fi);i<=nft(fi);i++)
          {
            vulnerable_numbers(i)=sum(vul_at_age(fi,i));  
          }
        }
        else
        {
          for (int i=first_survey_time(fi);i<=nft(fi);i++)
          {
            if (allocated(fishery_group_ptr) &&
                allocated(fishery_group_ptr(99)))
            {
              vulnerable_numbers(i)=sum(vul_at_age(fi,i));  
            }
            else
            {
              for (int isp=1;isp<=pmsd->num_species;isp++) 
              {
                int fisp=local_fishery_pointer(isp);
                vulnerable_numbers(i)+=sum(vul_at_age(fisp,i));
              }
            }
          }
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
        cpue_obs=logsurv-mls;
        cpue_pred=value(logvuln-mlv);
        survey_cpue_obs(fi)=cpue_obs;
        survey_cpue_pred(fi)=cpue_pred;
	tmps=mls;
	tmpv=value(mlv);
      }
      else
      {
        dvar_vector  vulnerable_biomass(first_survey_time(fi),nft(fi));
        vulnerable_biomass.initialize();
        if (!pmsd)
        {
          for (int i=first_survey_time(fi);i<=nft(fi);i++)
          {
            int rr=realization_region(fi,i); 
            int rp=realization_period(fi,i); 
            int ri=realization_incident(fi,i); 
            vulnerable_biomass(i)=mean_weight(rr,rp,ri)
              *vul_at_age(fi,i);
          }
        }
        else
        {
          for (int i=first_survey_time(fi);i<=nft(fi);i++)
          {
            int rr=realization_region(fi,i); 
            int rp=realization_period(fi,i); 
            int ri=realization_incident(fi,i); 
            if (allocated(fishery_group_ptr) &&
                allocated(fishery_group_ptr(99)))
            {
              vulnerable_biomass(i)=mean_weight(rr,rp,ri)
                *vul_at_age(fi,i);
            }
            else
            {
              for (int isp=1;isp<=pmsd->num_species;isp++) 
              {
                int fisp=local_fishery_pointer(isp);
                vulnerable_biomass(i)+=mean_weight(rr,rp,ri)
                  *vul_at_age(fisp,i);
              }
            }
          }
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
        cpue_obs=logsurv-mls;
        cpue_pred=value(logvuln-mlv);
        survey_cpue_obs(fi)=cpue_obs;	
        survey_cpue_pred(fi)=cpue_pred;	
	tmps=mls;
	tmpv=value(mlv);
      }
      if (!af170q0)
      {
        *ofs1 << " Fishery: " << fi << " " << setfixed() << setprecision(4)
              << setw(7) << cpue_obs << endl;    //NMD_30jul2021
        *ofs2 << " Fishery: " << fi << " " << setfixed() << setprecision(4)
              << setw(7) << cpue_pred << endl;   //NMD_30jul2021
        *ofs3 << " Fishery: " << fi << " " << setfixed() << setprecision(4)
              << setw(7) << tmps << " "    //NMD_05apr2051
              << setw(7) << cpue_obs << endl;    //NMD_30jul2021
        *ofs4 << " Fishery: " << fi << " " << setfixed() << setprecision(4)
              << setw(7) << tmpv << " "    //NMD_05apr2051
              << setw(7) << cpue_pred << endl;   //NMD_30jul2021
      }
    }
  }
}

#undef HOME_VERSION
