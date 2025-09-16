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
dvariable dvar_fish_stock_history::fit_survey_fishery_indices
(ivector& ff92, int print_switch)
{
  ivector ff66=column(fish_flags,66);
  ivector ff94=column(fish_flags,94);
  ivector nft;    //NMD_28jun2022
  if (sum(data_fish_flags(2)))
  {
    nft=num_real_fish_times;
  }
  else
  {
    nft=num_fish_times;
  }
  if (print_switch && !af170q0)
  {
    print_survey_fishery_fits(ff92);
  }
  MY_DOUBLE_TYPE pen=0.22;
  int use_self_scaling=0;   // use pen wght instead of self scaling likelihood
  dvariable ptmp=0.0;
  dvariable xy=0.0;
  ivector local_reg_pointer;
  ivector local_fishery_pointer;

  // check if catchability/variance for surver fishery catch cpue vs sruvey cpue
  // is grouped
  if (allocated(fishery_group_ptr) &&  allocated(fishery_group_ptr(99)))
  {
    int num_groups=fishery_group_ptr(99).indexmax();
    for (int ig=1;ig<=num_groups;ig++)
    {
      if (ff92(fishery_group_ptr(99,ig,1)))
      {
        int nobs=0;
        int nobs1=0;
        ivector fgp99=fishery_group_ptr(99,ig);
        ivector fshry_obs_ptr;
        ivector srtd_grp;
        ivector srtd_grp_sub;
        dvector grouped_survey_index;
        if (!pmsd)
        {
          fshry_obs_ptr.allocate(1,fgp99.indexmax());
          fshry_obs_ptr.initialize();
          for (int i=1;i<=fgp99.indexmax();i++)
          {
            int fi=fgp99(i);
            nobs+= nft(fi)-first_survey_time(fi)+1;
            nobs1+=survey_index(fi).indexmax()-survey_index(fi).indexmin()+1;
            fshry_obs_ptr(i)=nobs;
          }
          if (nobs != nobs1)
          {
            cerr << " size error for grouped survey data vs cpue" << endl;
            ad_exit(1);
          }
          grouped_survey_index.allocate(1,nobs);
          grouped_survey_index.initialize();
          int offset=0;
          for (int ii=1;ii<=fgp99.indexmax();ii++)
          {
            int fi=fgp99(ii);
            for (int i=first_survey_time(fi);i<=nft(fi);i++)
            {
              int i1=i+offset-first_survey_time(fi)+1;
              grouped_survey_index(i1)=survey_index(fi,i);
            }
            offset+=nft(fi)-first_survey_time(fi)+1;
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
	    
          fshry_obs_ptr.allocate(1,srtd_grp_sub.indexmax());
          fshry_obs_ptr.initialize();
          for (int i=1;i<=srtd_grp_sub.indexmax();i++)
          {
            int fi=srtd_grp_sub(i);
            nobs+= nft(fi)-first_survey_time(fi)+1;
            nobs1+=survey_index(fi).indexmax()-survey_index(fi).indexmin()+1;
            fshry_obs_ptr(i)=nobs;
          }
          if (nobs != nobs1)
          {
            cerr << " size error for grouped survey data vs cpue" << endl;
            ad_exit(1);
          }
          grouped_survey_index.allocate(1,nobs);
          grouped_survey_index.initialize();
          int offset=0;
          for (int ii=1;ii<=srtd_grp_sub.indexmax();ii++)
          {
            int fi=srtd_grp_sub(ii);
            for (int i=first_survey_time(fi);i<=nft(fi);i++)
            {
              int i1=i+offset-first_survey_time(fi)+1;
              grouped_survey_index(i1)=survey_index(fi,i);
            }
            offset+=nft(fi)-first_survey_time(fi)+1;
          }
        }

        if (data_fish_flags(1,fgp99(1))==0)
        {  
          dvar_vector grouped_vulnerable_numbers(1,nobs);
          grouped_vulnerable_numbers.initialize();
          int offset=0;
          dvector lambda(1,nobs);
          lambda.initialize();
          if (!pmsd)
          {
            for (int ii=1;ii<=fgp99.indexmax();ii++)
            {
              int fi=fgp99(ii);
              for (int i=first_survey_time(fi);i<=nft(fi);i++)
              {
                int i1=i+offset-first_survey_time(fi)+1;
                grouped_vulnerable_numbers(i1)=sum(vul_at_age(fi,i));  
                lambda(i1)=effort_weight_by_fishery(fi,i);
              }
              offset+=nft(fi)-first_survey_time(fi)+1;
            }
            for (int ii=1;ii<=fgp99.indexmax();ii++)
            {
              int fi=fgp99(ii);
              int strt;
              int fin;
              if (ii==1)
              {
                strt=1;
                fin=fshry_obs_ptr(ii);
              }
              else
              {
                strt=fshry_obs_ptr(ii-1)+1;
                fin=fshry_obs_ptr(ii);
              }
              double lambda_mn=mean(lambda(strt,fin));
              for (int i=strt;i<=fin;i++)
              {
                lambda(i)=lambda(i)/lambda_mn;
              }
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
              local_reg_pointer.allocate(1,pmsd->num_species);
              local_fishery_pointer.allocate(1,pmsd->num_species);
              int irr_sp=0;
              if (isp==1)
              {
                local_reg_pointer(isp)=irr;
                local_fishery_pointer(isp)=fi;
              }
              offset=0;
              for (int i=2;i<=pmsd->num_species;i++) 
              {
                offset=(i-1)*pmsd->num_real_fisheries;
                int irr_sp=fishery_regions(fi+offset);
                local_reg_pointer(i)=irr_sp;
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


// Concatenation code for fisheries in the group using srtd_grp_sub order
// as applied for obs
            offset=0;
            for (int ii=1;ii<=srtd_grp_sub.indexmax();ii++)
            {
              int fi=srtd_grp_sub(ii);
              for (int i=first_survey_time(fi);i<=nft(fi);i++)
              {
                int i1=i+offset-first_survey_time(fi)+1;
                grouped_vulnerable_numbers(i1)=sum(vul_at_age(fi,i));  
                lambda(i1)=effort_weight_by_fishery(fi,i);
              }
              offset+=nft(fi)-first_survey_time(fi)+1;
            }
// Normalise lambdas within each fishery of the group
            for (int ii=1;ii<=srtd_grp_sub.indexmax();ii++)
            {
              int fi=srtd_grp_sub(ii);
              int strt;
              int fin;
              if (ii==1)
              {
                strt=1;
                fin=fshry_obs_ptr(ii);
              }
              else
              {
                strt=fshry_obs_ptr(ii-1)+1;
                fin=fshry_obs_ptr(ii);
              }
              double lambda_mn=mean(lambda(strt,fin));
              for (int i=strt;i<=fin;i++)
              {
                lambda(i)=lambda(i)/lambda_mn;
              }
            }
          }
          dvector log_groupedsurv=log(1.e-10+grouped_survey_index);
          dvar_vector log_groupedvuln=log(1.e-10+grouped_vulnerable_numbers);
          double mls=mean(log_groupedsurv);
          dvariable mlv=mean(log_groupedvuln);
          if (!ff66(fgp99(1)))
          {
            int ff94sum=sum(ff94(fishery_group_ptr(99,ig)));
            if (!ff94sum)
            {
              dvariable diff2=norm2(log_groupedsurv-log_groupedvuln-(mls-mlv));
              int fi=fgp99(1);  // first fishery in this group
              if (ff92(fi) != 0)
                pen=1.0/(2.0*square(ff92(fi)/100.));
              ptmp=pen*diff2;
//////////    new code for constant term
              dvar_vector sigma(1,nobs);
              sigma.initialize();
              int fshgp_max;
              if (!pmsd)
              {
                fshgp_max=fgp99.indexmax();
              }
              else
              {
                fshgp_max=srtd_grp_sub.indexmax();
              }
//              for (int ii=1;ii<=fgp99.indexmax();ii++)
              for (int ii=1;ii<=fshgp_max;ii++)
              {
                int fi;
                if (!pmsd)
                {
                  fi=fgp99(ii);
                }
                else
                {
                  fi=srtd_grp_sub(ii);
                }
                int strt;
                int fin;
                if (ii==1)
                {
                  strt=1;
                  fin=fshry_obs_ptr(ii);
                }
                else
                {
                  strt=fshry_obs_ptr(ii-1)+1;
                  fin=fshry_obs_ptr(ii);
                }
                for (int i=strt;i<=fin;i++)
                {
                  sigma(i)=ff92(fi)/100.;
                }
              }
              dvariable vhat=0.0;
//              for (int ii=1;ii<=fgp99.indexmax();ii++)
              for (int ii=1;ii<=fshgp_max;ii++)
              {
                int fi;
                if (!pmsd)
                {
                  fi=fgp99(ii);
                }
                else
                {
                  fi=srtd_grp_sub(ii);
                }
                int strt;
                int fin;
                if (ii==1)
                {
                  strt=1;
                  fin=fshry_obs_ptr(ii);
                }
                else
                {
                  strt=fshry_obs_ptr(ii-1)+1;
                  fin=fshry_obs_ptr(ii);
                }
                for (int i=strt;i<=fin;i++)
                {
                  vhat+=log(square(sigma(i)));
                }
              }
              vhat=0.5*vhat;
              ptmp=ptmp+vhat;
//////////    new code for constant term
            }
            else
            {
              int fshgp_max;
              if (!pmsd)
              {
                fshgp_max=fgp99.indexmax();
              }
              else
              {
                fshgp_max=srtd_grp_sub.indexmax();
              }
              dvar_matrix log_grpdvuln_sub(1,fshgp_max);
              dmatrix log_grpdsurv_sub(1,fshgp_max);
              offset=0;
//              for (int ii=1;ii<=fgp99.indexmax();ii++)
              for (int ii=1;ii<=fshgp_max;ii++)
              {
                int fi;
                if (!pmsd)
                {
                  fi=fgp99(ii);
                }
                else
                {
                  fi=srtd_grp_sub(ii);
                }
                int strt;
                if (ii==1)
                {
                  strt=1;
                }
                else
                {
                  strt=fshry_obs_ptr(ii-1)+1;
                }
                int nsub=fshry_obs_ptr(ii)-strt+1;
                log_grpdvuln_sub(ii).allocate(1,nsub);
                log_grpdsurv_sub(ii).allocate(1,nsub);
                int icnt=0;
                for (int i=first_survey_time(fi);i<=nft(fi);i++)
                {
                  int i1=i+offset-first_survey_time(fi)+1;
                  icnt++;
                  log_grpdvuln_sub(ii,icnt)=log_groupedvuln(i1);
                  log_grpdsurv_sub(ii,icnt)=log_groupedsurv(i1);
                }
                offset+=nft(fi)-first_survey_time(fi)+1;
              }
              dvariable diff2=0.0;
//              for (int ii=1;ii<=fgp99.indexmax();ii++)
              for (int ii=1;ii<=fshgp_max;ii++)
              {
                int fi;
                if (!pmsd)
                {
                  fi=fgp99(ii);
                }
                else
                {
                  fi=srtd_grp_sub(ii);
                }
                if (ff92(fi) != 0)
                  pen=1.0/(2.0*square(ff92(fi)/100.));
                dvector otmp=log_grpdsurv_sub(ii);
                dvar_vector vtmp=log_grpdvuln_sub(ii);
                diff2+=pen*norm2(otmp-vtmp-(mls-mlv));
              }
              ptmp=diff2;
//////////    new code for constant term
              dvar_vector sigma(1,nobs);
              sigma.initialize();
//              for (int ii=1;ii<=fgp99.indexmax();ii++)
              for (int ii=1;ii<=fshgp_max;ii++)
              {
                int fi;
                if (!pmsd)
                {
                  fi=fgp99(ii);
                }
                else
                {
                  fi=srtd_grp_sub(ii);
                }
                int strt;
                int fin;
                if (ii==1)
                {
                  strt=1;
                  fin=fshry_obs_ptr(ii);
                }
                else
                {
                  strt=fshry_obs_ptr(ii-1)+1;
                  fin=fshry_obs_ptr(ii);
                }
                for (int i=strt;i<=fin;i++)
                {
                  sigma(i)=ff92(fi)/100.;
                }
              }
              dvariable vhat=0.0;
//              for (int ii=1;ii<=fgp99.indexmax();ii++)
              for (int ii=1;ii<=fshgp_max;ii++)
              {
                int fi;
                if (!pmsd)
                {
                  fi=fgp99(ii);
                }
                else
                {
                  fi=srtd_grp_sub(ii);
                }
                int strt;
                int fin;
                if (ii==1)
                {
                  strt=1;
                  fin=fshry_obs_ptr(ii);
                }
                else
                {
                  strt=fshry_obs_ptr(ii-1)+1;
                  fin=fshry_obs_ptr(ii);
                }
                for (int i=strt;i<=fin;i++)
                {
                  vhat+=log(square(sigma(i)));
                }
              }
              vhat=0.5*vhat;
              ptmp=ptmp+vhat;
//////////    new code for constant term
            }
            xy+=ptmp;
            if (ppstf && !af170q0)
            {
              ppstf->survey_index_like_by_group(ig)=value(ptmp);
            }
          }
          else
          {
            int ff94sum=sum(ff94(fishery_group_ptr(99,ig)));
            dvar_vector sigma(1,nobs);
            sigma.initialize();
            if (ff94sum)
            {
              int fshgp_max;
              if (!pmsd)
              {
                fshgp_max=fgp99.indexmax();
              }
              else
              {
                fshgp_max=srtd_grp_sub.indexmax();
              }
//              for (int ii=1;ii<=fgp99.indexmax();ii++)
              for (int ii=1;ii<=fshgp_max;ii++)
              {
                int fi;
                if (!pmsd)
                {
                  fi=fgp99(ii);
                }
                else
                {
                  fi=srtd_grp_sub(ii);
                }
                int strt;
                int fin;
                if (ii==1)
                {
                  strt=1;
                  fin=fshry_obs_ptr(ii);
                }
                else
                {
                  strt=fshry_obs_ptr(ii-1)+1;
                  fin=fshry_obs_ptr(ii);
                }
                for (int i=strt;i<=fin;i++)
                {
                  if (ff92(fishery_group_ptr(99,ig))(ii)!=0)
                  {
                    sigma(i)=ff92(fishery_group_ptr(99,ig))(ii)/100.;
                  }
                  else
                  {
                    sigma(i)=pen;
                  }
                }
              }
            }
            else
            {
              int fi=fgp99(1);  // first fishery in this group
              if (ff92(fi) != 0)
                pen=ff92(fi)/100.;
            }
            dvariable vhat=0.0;
            dvariable vhat_est=0.0;
            dvariable tmp=0.0;
            for (int i=1;i<=nobs;i++)
            {
              if (!ff94sum)
              {
                if (parest_flags(371)==0)
                {
                  vhat+=log(lambda(i) * square(pen));
                }
                else
                {
                  vhat+=log(lambda(i) + square(pen));
                }
              }
              else
              {
                if (parest_flags(371)==0)
                {
                  vhat+=log(lambda(i) * square(sigma(i)));
                }
                else
                {
                  vhat+=log(lambda(i) + square(sigma(i)));
                }
              }
            }

            dvariable vhat_wtd=0.0;
            dvariable norm2chk=0.0;
            for (int i=1;i<=nobs;i++)
            {
              if (!ff94sum)
              {
                if (parest_flags(371)==0)
                {
                  vhat_wtd+=square(log_groupedsurv(i)-log_groupedvuln(i)-
                    (mls-mlv))/(lambda(i)*square(pen));
                }
                else
                {
                  vhat_wtd+=square(log_groupedsurv(i)-log_groupedvuln(i)-
                    (mls-mlv))/(lambda(i)+square(pen));
                }
              }
              else
              {
                if (parest_flags(371)==0)
                {
                  vhat_wtd+=square(log_groupedsurv(i)-log_groupedvuln(i)-
                    (mls-mlv))/(lambda(i)*square(sigma(i)));
                  norm2chk+=square(log_groupedsurv(i)-log_groupedvuln(i)-
                    (mls-mlv));
                }
                else
                {
                  vhat_wtd+=square(log_groupedsurv(i)-log_groupedvuln(i)-
                    (mls-mlv))/(lambda(i)+square(sigma(i)));
                }
              }
            }
            ptmp=(0.5*vhat)+(0.5*vhat_wtd);
            xy+=ptmp;  
            if (ppstf && !af170q0)
            {
              ppstf->survey_index_like_by_group(ig)=value(ptmp);
            }
            if (print_switch && !af170q0)
            {
              for (int i=1;i<=nobs;i++)
              {
//             - estimated sigma from fitted likelihood
                vhat_est+=square(log_groupedsurv(i)-log_groupedvuln(i)-
                  (mls-mlv))/(nobs*lambda(i));
              }
              cout << " Estimated CPUE sigma for grouped fisheries "
                   << fishery_group_ptr(99,ig)
                   << " is : " << vhat_est << endl;
            }
          }
        }    
        else
        {
          dvar_vector  grouped_vulnerable_biomass(1,nobs);
          dvector lambda(1,nobs);
          int offset=0;
          if (!pmsd)
          {
//             dvar_vector  grouped_vulnerable_biomass(1,nobs);
            grouped_vulnerable_biomass.initialize();
//            dvector lambda(1,nobs);
            lambda.initialize();
//            int offset=0;
            for (int ii=1;ii<=fgp99.indexmax();ii++)
            {
              int fi=fgp99(ii);
              for (int i=first_survey_time(fi);i<=nft(fi);i++)
              {
                int rr=realization_region(fi,i); 
                int rp=realization_period(fi,i); 
                int ri=realization_incident(fi,i); 
                int i1=i+offset-first_survey_time(fi)+1;
                grouped_vulnerable_biomass(i1)=mean_weight(rr,rp,ri)
                  *vul_at_age(fi,i);
                lambda(i1)=effort_weight_by_fishery(fi,i);  //NMD_nov29_2021
              }
              offset+=nft(fi)-first_survey_time(fi)+1;
            }
            for (int ii=1;ii<=fgp99.indexmax();ii++)
            {
              int fi=fgp99(ii);
              int strt;
              int fin;
              if (ii==1)
              {
                strt=1;
                fin=fshry_obs_ptr(ii);
              }
              else
              {
                strt=fshry_obs_ptr(ii-1)+1;
                fin=fshry_obs_ptr(ii);
              }
              double lambda_mn=mean(lambda(strt,fin));
              for (int i=strt;i<=fin;i++)
              {
                lambda(i)=lambda(i)/lambda_mn;
              }
            }
          }
          else //pmsd
          {
//          Aggregated predictions over kludged fisheries in group	    
            for (int ii=1;ii<=srtd_grp_sub.indexmax();ii++)
            {
              int fi=srtd_grp_sub(ii);
              int irr=fishery_regions(fi);
              int isp=pmsd->region_species_pointer(irr);
              local_reg_pointer.allocate(1,pmsd->num_species);
              local_fishery_pointer.allocate(1,pmsd->num_species);
              int irr_sp=0;
              if (isp==1)
              {
                local_reg_pointer(isp)=irr;
                local_fishery_pointer(isp)=fi;
              }
              offset=0;
              for (int i=2;i<=pmsd->num_species;i++) 
              {
                offset=(i-1)*pmsd->num_real_fisheries;
                int irr_sp=fishery_regions(fi+offset);
                local_reg_pointer(i)=irr_sp;
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


// Concatenation code for fisheries in the group using srtd_grp_sub order
// as applied for obs
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
                lambda(i1)=effort_weight_by_fishery(fi,i);
              }
              offset+=nft(fi)-first_survey_time(fi)+1;
            }
// Normalise lambdas within each fishery of the group
            for (int ii=1;ii<=srtd_grp_sub.indexmax();ii++)
            {
              int fi=srtd_grp_sub(ii);
              int strt;
              int fin;
              if (ii==1)
              {
                strt=1;
                fin=fshry_obs_ptr(ii);
              }
              else
              {
                strt=fshry_obs_ptr(ii-1)+1;
                fin=fshry_obs_ptr(ii);
              }
              double lambda_mn=mean(lambda(strt,fin));
              for (int i=strt;i<=fin;i++)
              {
                lambda(i)=lambda(i)/lambda_mn;
              }
            }
          }  //end of pmsd

          dvector log_groupedsurv=log(1.e-10+grouped_survey_index);
          dvar_vector log_groupedvuln=
            log(1.e-10+grouped_vulnerable_biomass);
          double mls=mean(log_groupedsurv);
          dvariable mlv=mean(log_groupedvuln);
          if (!ff66(fgp99(1)))
          {
            int ff94sum=sum(ff94(fishery_group_ptr(99,ig)));
            if (!ff94sum)
            {
              dvariable diff2=norm2(log_groupedsurv-log_groupedvuln-(mls-mlv));
              int fi=fgp99(1);  // first fishery in this group
              if (ff92(fi) !=  0)
                pen=1.0/(2.0*square(ff92(fi)/100.0));
              ptmp=pen*diff2;

//////////    new code for constant term
              dvar_vector sigma(1,nobs);
              sigma.initialize();
              int fshgp_max;
              if (!pmsd)
              {
                fshgp_max=fgp99.indexmax();
              }
              else
              {
                fshgp_max=srtd_grp_sub.indexmax();
              }
//              for (int ii=1;ii<=fgp99.indexmax();ii++)
              for (int ii=1;ii<=fshgp_max;ii++)
              {
                int fi;
                if (!pmsd)
                {
                  fi=fgp99(ii);
                }
                else
                {
                  fi=srtd_grp_sub(ii);
                }
                int strt;
                int fin;
                if (ii==1)
                {
                  strt=1;
                  fin=fshry_obs_ptr(ii);
                }
                else
                {
                  strt=fshry_obs_ptr(ii-1)+1;
                  fin=fshry_obs_ptr(ii);
                }
                for (int i=strt;i<=fin;i++)
                {
                  sigma(i)=ff92(fi)/100.0;		  
                }
              }
              dvariable vhat=0.0;
//              for (int ii=1;ii<=fgp99.indexmax();ii++)
              for (int ii=1;ii<=fshgp_max;ii++)
              {
                int fi;
                if (!pmsd)
                {
                  fi=fgp99(ii);
                }
                else
                {
                  fi=srtd_grp_sub(ii);
                }
//                int fi=fgp99(ii);
                int strt;
                int fin;
                if (ii==1)
                {
                  strt=1;
                  fin=fshry_obs_ptr(ii);
                }
                else
                {
                  strt=fshry_obs_ptr(ii-1)+1;
                  fin=fshry_obs_ptr(ii);
                }
                for (int i=strt;i<=fin;i++)
                {
                  vhat+=log(square(sigma(i)));
                }
              }
              vhat=0.5*vhat;
              ptmp=ptmp+vhat;
//////////    new code for constant term
            }
            else
            {
              int fshgp_max;
              if (!pmsd)
              {
                fshgp_max=fgp99.indexmax();
              }
              else
              {
                fshgp_max=srtd_grp_sub.indexmax();
              }
              dvar_matrix log_grpdvuln_sub(1,fshgp_max);
              dmatrix log_grpdsurv_sub(1,fshgp_max);
              offset=0;
//              for (int ii=1;ii<=fgp99.indexmax();ii++)
              for (int ii=1;ii<=fshgp_max;ii++)
              {
                int fi;
                if (!pmsd)
                {
                  fi=fgp99(ii);
                }
                else
                {
                  fi=srtd_grp_sub(ii);
                }
                int strt;
                if (ii==1)
                {
                  strt=1;
                }
                else
                {
                  strt=fshry_obs_ptr(ii-1)+1;
                }
                int nsub=fshry_obs_ptr(ii)-strt+1;
                log_grpdvuln_sub(ii).allocate(1,nsub);
                log_grpdsurv_sub(ii).allocate(1,nsub);
                int icnt=0;
                for (int i=first_survey_time(fi);i<=nft(fi);i++)
                {
                  int i1=i+offset-first_survey_time(fi)+1;
                  icnt++;
                  log_grpdvuln_sub(ii,icnt)=log_groupedvuln(i1);
                  log_grpdsurv_sub(ii,icnt)=log_groupedsurv(i1);
                }
                offset+=nft(fi)-first_survey_time(fi)+1;
              }
              dvariable diff2=0.0;
//              for (int ii=1;ii<=fgp99.indexmax();ii++)
              for (int ii=1;ii<=fshgp_max;ii++)
              {
                int fi;
                if (!pmsd)
                {
                  fi=fgp99(ii);
                }
                else
                {
                  fi=srtd_grp_sub(ii);
                }
                if (ff92(fi) != 0)
                  pen=1.0/(2.0*square(ff92(fi)/100.));
                dvector otmp=log_grpdsurv_sub(ii);
                dvar_vector vtmp=log_grpdvuln_sub(ii);
                diff2+=pen*norm2(otmp-vtmp-(mls-mlv));
              }
              ptmp=diff2;

//////////    new code for constant term
              dvar_vector sigma(1,nobs);
              sigma.initialize();
//              for (int ii=1;ii<=fgp99.indexmax();ii++)
              for (int ii=1;ii<=fshgp_max;ii++)
              {
                int fi;
                if (!pmsd)
                {
                  fi=fgp99(ii);
                }
                else
                {
                  fi=srtd_grp_sub(ii);
                }
                int strt;
                int fin;
                if (ii==1)
                {
                  strt=1;
                  fin=fshry_obs_ptr(ii);
                }
                else
                {
                  strt=fshry_obs_ptr(ii-1)+1;
                  fin=fshry_obs_ptr(ii);
                }
                for (int i=strt;i<=fin;i++)
                {
                  sigma(i)=ff92(fi)/100.;
                }
              }
              dvariable vhat=0.0;
//              for (int ii=1;ii<=fgp99.indexmax();ii++)
              for (int ii=1;ii<=fshgp_max;ii++)
              {
                int fi;
                if (!pmsd)
                {
                  fi=fgp99(ii);
                }
                else
                {
                  fi=srtd_grp_sub(ii);
                }
                int strt;
                int fin;
                if (ii==1)
                {
                  strt=1;
                  fin=fshry_obs_ptr(ii);
                }
                else
                {
                  strt=fshry_obs_ptr(ii-1)+1;
                  fin=fshry_obs_ptr(ii);
                }
                for (int i=strt;i<=fin;i++)
                {
                  vhat+=log(square(sigma(i)));
                }
              }
              vhat=0.5*vhat;
              ptmp=ptmp+vhat;
//////////    new code for constant term
            }
            xy+=ptmp;
            if (ppstf && !af170q0)
            {
              ppstf->survey_index_like_by_group(ig)=value(ptmp);
            }
          }
          else
          {
            int ff94sum=sum(ff94(fishery_group_ptr(99,ig)));
            dvar_vector sigma(1,nobs);
            sigma.initialize();
            if (ff94sum)
            {
              int fshgp_max;
              if (!pmsd)
              {
                fshgp_max=fgp99.indexmax();
              }
              else
              {
                fshgp_max=srtd_grp_sub.indexmax();
              }
//              for (int ii=1;ii<=fgp99.indexmax();ii++)
              for (int ii=1;ii<=fshgp_max;ii++)
              {
                int fi;
                if (!pmsd)
                {
                  fi=fgp99(ii);
                }
                else
                {
                  fi=srtd_grp_sub(ii);
                }
                int strt;
                int fin;
                if (ii==1)
                {
                  strt=1;
                  fin=fshry_obs_ptr(ii);
                }
                else
                {
                  strt=fshry_obs_ptr(ii-1)+1;
                  fin=fshry_obs_ptr(ii);
                }
                for (int i=strt;i<=fin;i++)
                {
                  if (ff92(fishery_group_ptr(99,ig))(ii)!=0)
                  {
                    sigma(i)=ff92(fishery_group_ptr(99,ig))(ii)/100.;
                  }
                  else
                  {
                    sigma(i)=pen;
                  }
                }
              }
            }
            else
            {
              int fi=fgp99(1);  // first fishery in this group
              if (ff92(fi) != 0)
                pen=ff92(fi)/100.;
            }
            dvariable vhat=0.0;
            dvariable vhat_est=0.0;
            dvariable tmp=0.0;
            for (int i=1;i<=nobs;i++)
            {
              if (!ff94sum)
              {
                if (parest_flags(371)==0)
                {
                  vhat+=log(lambda(i) * square(pen));
                }
                else
                {
                  vhat+=log(lambda(i) + square(pen));
                }
              }
              else
              {
                if (parest_flags(371)==0)
                {
                  vhat+=log(lambda(i) * square(sigma(i)));
                }
                else
                {
                  vhat+=log(lambda(i) + square(sigma(i)));
                }
              }
            }
	    
            dvariable vhat_wtd=0.0;
            for (int i=1;i<=nobs;i++)
            {
              if (!ff94sum)
              {
                if (parest_flags(371)==0)
                {
                  vhat_wtd+=square(log_groupedsurv(i)-log_groupedvuln(i)-
                    (mls-mlv))/(lambda(i)*square(pen));
                }
                else
                {
                  vhat_wtd+=square(log_groupedsurv(i)-log_groupedvuln(i)-
                    (mls-mlv))/(lambda(i)+square(pen));
                }
              }
              else
              {
                if (parest_flags(371)==0)
                {
                  vhat_wtd+=square(log_groupedsurv(i)-log_groupedvuln(i)-
                    (mls-mlv))/(lambda(i)*square(sigma(i)));
                }
                else
                {
                  vhat_wtd+=square(log_groupedsurv(i)-log_groupedvuln(i)-
                    (mls-mlv))/(lambda(i)+square(sigma(i)));
                }
              }
            }
            ptmp=(0.5*vhat)+(0.5*vhat_wtd);
            xy+=ptmp;  
            if (ppstf && !af170q0)
            {
              ppstf->survey_index_like_by_group(ig)=value(ptmp);
            }
            if (print_switch && !af170q0)
            {
              for (int i=1;i<=nobs;i++)
              {
//             - estimated sigma from fitted likelihood
                vhat_est+=square(log_groupedsurv(i)-log_groupedvuln(i)-
                  (mls-mlv))/(nobs*lambda(i));
              }
	      cout << " Estimated CPUE sigma for grouped fisheries "
                   << fishery_group_ptr(99,ig)
                   << " is : " << vhat_est << endl;
            }
          }
        }    
      }
    }
  }
  else
  {
    for (int fi=1;fi<=num_fisheries;fi++)
    {
      int fspp1=1;
      if (ff92(fi) && pmsd)
      {
        if (fi<=pmsd->num_real_fisheries)
        {
          int irr=fishery_regions(fi);
          int isp=pmsd->region_species_pointer(irr);
          local_reg_pointer.allocate(1,pmsd->num_species);
          local_fishery_pointer.allocate(1,pmsd->num_species);
          int irr_sp=0;
          if (isp==1)
          {
            local_reg_pointer(isp)=irr;
            local_fishery_pointer(isp)=fi;
          }
          for (int i=2;i<=pmsd->num_species;i++) 
          {
            int offset=(i-1)*pmsd->num_real_fisheries;
            int irr_sp=fishery_regions(fi+offset);
            local_reg_pointer(i)=irr_sp;
            local_fishery_pointer(i)=fi+offset;
          }
          fspp1=1;
        }
        else
        {
          fspp1=0;
        }
      }
      if (ff92(fi) && fspp1==1)
      {
        int nobs=survey_index(fi).indexmax()-survey_index(fi).indexmin()+1;
        dvector cpue_obs(first_survey_time(fi),nft(fi));
        dvector cpue_pred(first_survey_time(fi),nft(fi));
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
              for (int isp=1;isp<=pmsd->num_species;isp++) 
              {
                int fisp=local_fishery_pointer(isp);
                vulnerable_numbers(i)+=sum(vul_at_age(fisp,i));
              }
            }
          }
          dvector logsurv=log(1.e-10+survey_index(fi));
          dvar_vector logvuln=log(1.e-10+vulnerable_numbers);
          double mls=mean(logsurv);
          dvariable mlv=mean(logvuln);
          cpue_obs=logsurv-mls;
          cpue_pred=value(logvuln-mlv);
          if (!fish_flags(fi,66))
          {    
            if (ff92(fi) != 0)
              pen=1.0/(2.*square(ff92(fi)/100.0));
            dvariable vhat_wtd=norm2(logsurv-logvuln-(mls-mlv));
            ptmp=pen*vhat_wtd;   //NMD_22sep2023

//////////    new code for constant term
            dvar_vector sigma(1,nobs);
            sigma.initialize();
            for (int i=1;i<=nobs;i++)
            {
              sigma(i)=ff92(fi)/100.0;
            }
            dvariable vhat=0.0;
            for (int i=1;i<=nobs;i++)
            {
              vhat+=log(square(sigma(i)));
            }
            vhat=0.5*vhat;
            ptmp=ptmp+vhat;
//////////    new code for constant term
            xy+=ptmp;
            if (ppstf && !af170q0)
            {
              ppstf->survey_index_like_by_fishery(fi)=value(ptmp);
            }
          }
          else
          {  
            dvariable sigma;
            sigma=pen;
            if (ff92(fi) != 0)
              sigma=ff92(fi)/100.0;
            dvector lambda(first_survey_time(fi),nft(fi));
            lambda=
              effort_weight_by_fishery(fi)(first_survey_time(fi),nft(fi));
            double lambda_mn=mean(lambda);
            lambda/=lambda_mn;
            dvariable vhat=0.0;
            dvariable vhat_est=0.0;
            for (int i=first_survey_time(fi);i<=nft(fi);i++)
            {
              if (parest_flags(371)==0)
              {
                vhat+=log(lambda(i) * square(sigma));
              }
              else
              {
                vhat+=log(lambda(i) + square(sigma));
              }
            }
  	  
            dvariable vhat_wtd=0.0;
            for (int i=first_survey_time(fi);i<=nft(fi);i++)
            {
              if (parest_flags(371)==0)
              {
                vhat_wtd+=square(logsurv(i)-logvuln(i)-
                  (mls-mlv))/(lambda(i)*square(sigma));
              }
              else
              {
                vhat_wtd+=square(logsurv(i)-logvuln(i)-
                  (mls-mlv))/(lambda(i)+square(sigma));
              }
            }
            ptmp=(0.5*vhat)+(0.5*vhat_wtd);
            xy+=ptmp;
            if (ppstf && !af170q0)
            {
              ppstf->survey_index_like_by_fishery(fi)=value(ptmp);
            }
            if (print_switch && !af170q0)
            {
              for (int i=1;i<=nobs;i++)
              {
//             - estimated sigma from fitted likelihood
                vhat_est+=square(logsurv(i)-logvuln(i)-
                  (mls-mlv))/(nobs*lambda(i));
              }
	      cout << " Estimated CPUE sigma for fishery "
                   << fi << " is : " << vhat_est << endl;
            }
          }
        }
        else
        {
          dvar_vector  vulnerable_biomass(first_survey_time(fi),nft(fi));
          vulnerable_biomass.initialize();
//          for (int i=first_survey_time(fi);i<=nft(fi);i++)
//          {
//            int rr=realization_region(fi,i); 
//            int rp=realization_period(fi,i); 
//            int ri=realization_incident(fi,i); 
//            vulnerable_biomass(i)=mean_weight(rr,rp,ri)
//              *vul_at_age(fi,i);
//          }
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
              for (int isp=1;isp<=pmsd->num_species;isp++) 
              {
                int fisp=local_fishery_pointer(isp);
                vulnerable_biomass(i)+=mean_weight(rr,rp,ri)
                  *vul_at_age(fisp,i);
              }
            }
          }

          dvector logsurv=log(1.e-10+survey_index(fi));
          dvar_vector logvuln=log(1.e-10+vulnerable_biomass);
          double mls=mean(logsurv);
          dvariable mlv=mean(logvuln);
          cpue_obs=logsurv-mls;
          cpue_pred=value(logvuln-mlv);
          if (!fish_flags(fi,66))
          {    
            if (ff92(fi) != 0)
              pen=1.0/(2.*square(ff92(fi)/100.0));
            dvariable vhat_wtd=norm2(logsurv-logvuln-(mls-mlv));
            ptmp=pen*vhat_wtd;

//////////    new code for constant term
            dvar_vector sigma(1,nobs);
            sigma.initialize();
            for (int i=1;i<=nobs;i++)
            {
              sigma(i)=ff92(fi)/100.0;
            }
            dvariable vhat=0.0;
            for (int i=1;i<=nobs;i++)
            {
              vhat+=log(square(sigma(i)));
            }
            vhat=0.5*vhat;
            ptmp=ptmp+vhat;
//////////    new code for constant term

            xy+=ptmp;
            if (ppstf && !af170q0)
            {
              ppstf->survey_index_like_by_fishery(fi)=value(ptmp);
            }
          }
          else
          {  
            dvariable sigma;
            sigma=pen;
            if (ff92(fi) != 0)
              sigma=ff92(fi)/100.0;
            dvector lambda(first_survey_time(fi),nft(fi));
            lambda=
              effort_weight_by_fishery(fi)(first_survey_time(fi),nft(fi));
            double lambda_mn=mean(lambda);
            lambda/=lambda_mn;
            dvariable vhat=0.0;
            dvariable vhat_est=0.0;
            dvariable tmp=0.0;
            for (int i=first_survey_time(fi);i<=nft(fi);i++)
            {
              if (parest_flags(371)==0)
              {
                vhat+=log(lambda(i) * square(sigma));
              }
              else
              {
                vhat+=log(lambda(i) + square(sigma));
              }
            }
  	  
            dvariable vhat_wtd=0.0;
            for (int i=first_survey_time(fi);i<=nft(fi);i++)
            {
              if (parest_flags(371)==0)
              {
                vhat_wtd+=square(logsurv(i)-logvuln(i)-
                  (mls-mlv))/(lambda(i)*square(sigma));
              }
              else
              {
                vhat_wtd+=square(logsurv(i)-logvuln(i)-
                  (mls-mlv))/(lambda(i)+square(sigma));
              }
            }
  	  
            ptmp=(0.5*vhat)+(0.5*vhat_wtd);
            xy+=ptmp;
            if (ppstf && !af170q0)
            {
              ppstf->survey_index_like_by_fishery(fi)=value(ptmp);
            }
            if (print_switch && !af170q0)
            {
              for (int i=1;i<=nobs;i++)
              {
//             - estimated sigma from fitted likelihood
                vhat_est+=square(logsurv(i)-logvuln(i)-
                  (mls-mlv))/(nobs*lambda(i));
              }
	      cout << " Estimated CPUE sigma for fishery "
                   << fi << " is : " << vhat_est << endl;
            }
          }
        }
      }
    }
  }
  return xy;
}
  

dvariable rel_bio_fit(dvar_len_fish_stock_history& fsh)
{
  dvariable x=0.0;

  dvar_vector ssum(1,fsh.num_fisheries);
  ssum.initialize();
  int i;
  for (i=1;i<=fsh.num_fisheries;i++)
  {
    int ii=0;
    int j;
    for (j=1;j<=fsh.num_fish_times(i);j++)
    {
      int rr=fsh.realization_region(i,j);
      int rp=fsh.realization_period(i,j);
      int ri=fsh.realization_incident(i,j);
              
      if (fsh.biomass_index(i,j)>0)
      {
        ii++;
        ssum(i)+=fsh.RB(rr,rp,ri);
      }
    }
    ssum(i)/=ii;
    for (j=1;j<=fsh.num_fish_times(i);j++)
    {
      int rr=fsh.realization_region(i,j);
      int rp=fsh.realization_period(i,j);
      int ri=fsh.realization_incident(i,j);
              
      if (fsh.biomass_index(i,j)>0)
      {
        fsh.RB(rr,rp,ri)/=ssum(i);
      }
    }
  }
  for (i=1;i<=fsh.num_fisheries;i++)
  {
    for (int j=1;j<=fsh.num_fish_times(i);j++)
    {
      int rr=fsh.realization_region(i,j);
      int rp=fsh.realization_period(i,j);
      int ri=fsh.realization_incident(i,j);
              
      if (fsh.biomass_index(i,j)>0)
      {
        MY_DOUBLE_TYPE wt=fsh.age_flags(61);
        const prevariable& x1=fsh.RB(rr,rp,ri);
        MY_DOUBLE_TYPE& x2=fsh.biomass_index(i,j);
        x+=wt*square(log(x1/x2));
        //x+=wt*square(log(fsh.RB(rr,rp,ri)/fsh.biomass_index(i,j)));
      }
    }
  }
  return x;
}

//dvariable dvar_len_fish_stock_history::fit_totals(int acf)
dvariable dvar_len_fish_stock_history::fit_totals(int print_switch,int acf)
{
  dvariable xy=0.0;

  if (!age_flags(92))
  {
    if (parest_flags(41)==0)
    {
      if (!sum(data_fish_flags(1)))
        xy+=total_catch_fit(*this,0,acf);
      else
      {
        //if (!pmsd)
          xy+=total_catch_or_weight_fit(*this,print_switch,acf);
        //else
        //  xy+=total_catch_or_weight_fit_ms(*this,0,acf);
      }
    }
    else
    {
      xy+=total_weight_fit(*this,0,acf);
    }
  }


  return xy;
}

dvariable objective_function(dvar_len_fish_stock_history& fsh,
  dvariable& mean_length_constraints,int print_switch,ofstream* of_pen)
{
  //cout << "A0" << flush;
  adtimer adt2;
  char ch;
  dvar_vector eff_dev_pen(1,fsh.num_fisheries);
  dvariable xy=0.;
  dvariable tmp1=0.;
  dvariable tmp2=0.;
  dvariable tmp_slave=0.0;
  dvariable tmp_slave_3=0.0;
  dvariable tmp_slave_2=0.0;

  //cout << " A" << flush;
  // do all the penalty function stuff
    // biomass depletion target 
      if (fsh.age_flags(170)>0)
      {
        if (fsh.age_flags(175)>0)
        {
          xy+=fsh.calculate_the_totalbiomass_depletion();
        }
      }
  
  if ( mf_pvm->pvm_switch ==0 || mf_pvm->pvm_switch ==1 ) 
  {
    // these are penalties that have no effect on the
    // solution. they are there to prevent signularities in the Hessian 
    tmp1=no_penalties(fsh,mean_length_constraints,print_switch,of_pen);
    cout << "After no_penalties pen = " << setprecision(12) <<  tmp1 << endl; //NMD
    xy+=tmp1;

    if (fsh.parest_flags(219)==0)
    {
      tmp1=good_penalties(fsh,mean_length_constraints,print_switch,of_pen);
      cout << "After good_penalties pen = " << setprecision(12) <<   tmp1 << endl;//NMD
      xy+=tmp1;
    }

    if (!fsh.age_flags(52) && fsh.parest_flags(165)==0)
    {
      if (fsh.parest_flags(220)==0)
      {
        tmp1=call_penalties(fsh,mean_length_constraints,print_switch,of_pen);
        cout << "After call_penalties pen = " << setprecision(12) <<  tmp1 << endl;//NMD
        xy+=tmp1;
      }
    }
  
    tmp1=0.0;
    // is this length or age data?
    int tm2=static_cast<int>(adt2.get_elapsed_time_and_reset());

    if (fsh.age_flags(40)==0)
    {
      if (fsh.nlint)
        xy+=catch_at_length_fit(fsh,print_switch);
  
      if (fsh.nwint && !fsh.parest_flags(181))
        xy+=catch_at_weight_fit(fsh,print_switch);
    }
    else
    {
      //cout << " D" << flush;
      xy+=catch_at_age_fit(fsh,print_switch);
    }

    if (fsh.age_nage)
    {
      switch(fsh.age_flags(124))
      {
      case 0:
        xy+=catch_at_age_fit(fsh,print_switch);
        break;
      case 1:
        xy+=catch_at_age_fit_dirichlet_multinomial(fsh,print_switch);
        break;
      default:
        cerr << "Illegal value for flags age_flags(124) " << fsh.age_flags(124)
             << endl;
        ad_exit(1);
      }
    }

    if (!ad_comm::pthread_manager || ad_comm::pthread_manager->ngroups<3)
    {
      xy+=fsh.fit_totals(print_switch,0);
    }
  }
   adtimer adt1;

  if (fsh.parest_flags(240))
  {
    xy+=fsh.fit_age_length_data();
  }
 
  if (fsh.num_tag_releases && fsh.age_flags(122) ==0)
  {
    dvariable tmp=0.0;
    // for now using switch 121 for survival anysis type analysis
    if (fsh.age_flags(121))
    {
      tmp=fsh.fit_tag_returns_survival_analysis();
      fsh.likecomponent[0]=value(tmp);
      xy+=tmp;
      cout << " Tag data = " << setprecision(12) <<  fsh.likecomponent[0] << endl;//NMD
    }
    else
    {
      int get_slave_flag=1;
      int tm2=static_cast<int>(adt2.get_elapsed_time_and_reset());
//      cout << "spent " << tm2/1000.
//           << " seconds at the end  of objective function" << endl;
      if (!ad_comm::pthread_manager || fsh.threaded_tag_flag==0 ||
          get_slave_flag==0)
      {
        if (!fsh.parest_flags(101))
        {
          if (!fsh.age_flags(100))
          {
            if (!fsh.pmsd || fsh.pmsd->combined_tags_flag==0)
            {
              //char s1[31]="start fit tag return123456789";
              char s1[31]="ABCDEFGHIJKLMNOP";
              //set_grad_profiler_30(s1);

              if (!fsh.parest_flags(249))
                tmp=fsh.fit_tag_returns();
              else
                tmp=fsh.fit_tag_returns_like_ss3();

              //char s2[31]="end fit tag return123456789";
              //char s2[31]="12345678901234567";
              //set_grad_profiler_30(s2);
            }
            else
            {
              if (!fsh.parest_flags(249))
                tmp=fsh.fit_tag_returns2();
              else
              {
                 cerr << "not implemented" << endl;
                 ad_exit(1);
              }
            }
          }
          else
            tmp=fsh.fit_tag_returns_mix();
        }
        else
        {
          tmp=fsh.fit_tag_returns_sqrt();
        }
      }
      else
      {
        {
          // get tag stuff from slave this should be done as late as possible
          adt1.get_elapsed_time_and_reset();
          int mmin=ad_comm::pthread_manager->gmin(1);
          int mmax=ad_comm::pthread_manager->gmax(1);
          if (get_slave_flag)
          {
            for (int sn=mmin;sn<=mmax;sn++)
            {
              ad_comm::pthread_manager->read_lock_buffer(sn);
              
              tmp_slave+=ad_comm::pthread_manager->get_dvariable(sn);
              ad_comm::pthread_manager->read_unlock_buffer(sn);
            }
            int tm1=static_cast<int>(adt1.get_elapsed_time_and_reset());
            cout << "waited " << tm1/1000. << " seconds for group 1" << endl;
            if (ad_comm::pthread_manager->ngroups>1)
            {
              int mmin=ad_comm::pthread_manager->gmin(2);
              int mmax=ad_comm::pthread_manager->gmax(2);
              for (int sn=mmin;sn<=mmax;sn++)
              {
                ad_comm::pthread_manager->read_lock_buffer(sn);
                tmp_slave_2+=ad_comm::pthread_manager->get_dvariable(sn);
                cout << "From group 2 tmp_slave = " << tmp_slave 
                     << " for slave number " << sn << endl;
                ad_comm::pthread_manager->read_unlock_buffer(sn);
              }
              int tm1=static_cast<int>(adt1.get_elapsed_time_and_reset());
              cout << "waited " << tm1/1000. << " seconds for group 2" << endl;
            }
            if (ad_comm::pthread_manager->ngroups>2)
            {
              int mmin=ad_comm::pthread_manager->gmin(3);
              int mmax=ad_comm::pthread_manager->gmax(3);
              for (int sn=mmin;sn<=mmax;sn++)
              {
                ad_comm::pthread_manager->read_lock_buffer(sn);
                tmp_slave_3+=ad_comm::pthread_manager->get_dvariable(sn);
                cout << "From group 3 tmp_slave_3 = " << tmp_slave_3 << endl;
                ad_comm::pthread_manager->read_unlock_buffer(sn);
              }
              int tm1=static_cast<int>(adt1.get_elapsed_time_and_reset());
              cout << "waited " << tm1/1000. << " seconds for group 3" << endl;
            }
          }
        }
      }
  
      if (!fsh.parest_flags(177))
      {
        fsh.likecomponent[0]=value(tmp);
        xy+=tmp;
      }
      else
      {
        MY_DOUBLE_TYPE wt=fsh.parest_flags(177)/1000.;
        fsh.likecomponent[0]=wt*value(tmp);
        xy+=wt*tmp;
      }
      cout << " Tag data = " << setprecision(12) <<  fsh.likecomponent[0] << endl;//NMD
    }
  }
 

  if ( mf_pvm->pvm_switch ==0 || mf_pvm->pvm_switch ==1 ) 
  {
    if (fsh.age_flags(61))
    {
      xy+=rel_bio_fit(fsh);
    }
  }

  ivector ff92=column(fsh.fish_flags,92);
//  int ff92sum=sum(ff92);
  if (fsh.ff92sum)
  {
    dvariable tmp=0.0;
    tmp=fsh.fit_survey_fishery_indices(ff92,print_switch);
    xy+=tmp;
      cout << " Survey index term = " << setprecision(12) <<  tmp << endl;//NMD
  }

  /*  this is now done in alldevpen */
  /*
  if (fsh.age_flags(92)>0)
  {
    xy+=fsh.fit_implicit_q();
  }
  */
  //  cout << "ZYY " << tmp_slave   << endl;
  //  cout << "DDY " << tmp_slave_3 << endl;

  xy+=tmp_slave+tmp_slave_2+tmp_slave_3;
//  xy+=fsh.newfpen;
//  fsh.newfpen=0.0;
  return xy;
}

#if defined(__MSVC32__)
MY_DOUBLE_TYPE mmiinn(const MY_DOUBLE_TYPE& x,MY_DOUBLE_TYPE y)
#else
MY_DOUBLE_TYPE mmiinn(const MY_DOUBLE_TYPE& x,MY_DOUBLE_TYPE y)
#endif
{
  if (x<y) return x;
  return y;
}

#undef HOME_VERSION
