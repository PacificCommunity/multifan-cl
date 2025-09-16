/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#define HOME_VERSION
#include "all.hpp"
//#include "here.h"
#include "variable.hpp"
#include "newmprot.hpp"

extern dvar_vector * psv;

dvariable calculate_overall_exploitation(dvar_len_fish_stock_history& fsh);
dvariable avail_coff_penalty(dvar_fish_stock_history& fsh);
dvariable selectivity_form_penalty(dvar_len_fish_stock_history& fsh,
                                   int print_switch);
dvariable recrpen(dvar_len_fish_stock_history& fsh,int is,int print_switch,
		  MY_DOUBLE_TYPE recr_weight); //NMD14May2012
dvariable recruitment_autocorrelation(dvar_len_fish_stock_history& fsh,int is,
  int print_switch, MY_DOUBLE_TYPE recr_weight);

dvariable recruitment_autocorrelation_moment_estimator
  (dvar_len_fish_stock_history& fsh,int is,int print_switch,
   MY_DOUBLE_TYPE recr_weight);

dvariable bh_recruitment_autocorrelation_moment_estimator
  (dvar_len_fish_stock_history& fsh,int is,int print_switch,
   MY_DOUBLE_TYPE wght);

dvariable adhoc(dvar_len_fish_stock_history& fsh)
{
  char ch;
  dvariable tmp =  fsh.sv(10)*exp(fsh.sv(13));
  dvariable tmp1;
  if (tmp<.2) tmp1=square(log(0.2/(1.e-10+tmp))); 
  return tmp1;

}

dvariable call_penalties(dvar_len_fish_stock_history& fsh,
  dvariable& mean_length_constraints,int print_switch,ofstream* pof_pen)
{
  dvariable xy=0.;
  dvariable ppf_tmp=0.0;   //NMD_21nov2023

  if (fsh.parest_flags(388)) 
  {
    int pen=fsh.parest_flags(388); 
    xy+=pen*fsh.kludged_selmean_square;
  }

  if (fsh.parest_flags(395))  
  {
    xy+=norm2(fsh.diff_coffs-.01);
  }

  if (fsh.parest_flags(393))  //NMD_27may2021
  {
    int sd=fsh.parest_flags(374);
    if (fsh.parest_flags(379) && sd>4)
    {
      int pen=10;
      if (fsh.parest_flags(379))
      {
        pen=fsh.parest_flags(379);
      }
      ivector  mps=fsh.get_equilibrium_movements_per_season();
      for (int is=1;is<=fsh.age_flags(57);is++)
      {
        for (int im=1;im<=mps(is);im++)
        {
          for (int ir=1;ir<=fsh.num_regions;ir++)
          {
            dvar_vector tsel=
              fsh.kludged_equilib_coffs(is,im,ir)(sd-4,sd);
            xy+=pen*norm2(tsel);
          }
        }
      }
    }
  }

  if (fsh.parest_flags(377))  //NMD_27may2021
  {
    ppf_tmp=xy;
    ivector ff73=column(fsh.fish_flags,73);
    ivector ff27=column(fsh.fish_flags,27);
    for (int fi=1;fi<=fsh.num_fisheries;fi++)  
    {
      int icol=1;
      if (ff27(fi)>0)
      {
        icol+=2;
      }
      int lb=icol+2;
      if (ff73(fi)>0)
      {
        icol+=ff73(fi);
      }
      int ub=icol;

      if (ub>=lb)
      {
        int pen=10;
        if (fsh.parest_flags(383))
        {
          pen=fsh.parest_flags(383);
        }
        xy+=pen*norm2(fsh.implicit_fm_level_regression_pars(fi)(lb,ub));
      }
    }
    ppf_tmp=value(xy)-ppf_tmp;  //NMD_21nov2023
    if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_21nov2023
    {
      fsh.ppstf->impl_fm_level_regr_pars_pen=ppf_tmp;
      ppf_tmp=0.0;
    }
  }  //NMD27may2021

  if (fsh.age_flags(125))
  {
    int nya=1; // num_years_for_average
    //   VVVVVV
    nya=fsh.age_flags(95); // num_years_for_average
    if (nya==0)
    {
      cerr << "Error age_flags(95) must be >0 " << endl;
      ad_exit(1);
    }
    dvar_matrix tmp=fsh.region_pars.sub(2,nya);
    dvar_matrix tmp1=fsh.current_biomass_by_year.sub(2,nya);
    MY_DOUBLE_TYPE pen=1000.;
    if (fsh.age_flags(127))
      pen=fsh.age_flags(127);
    xy+=pen*norm2(tmp-tmp1);
    cout << "BBB " << xy << endl;
    #if defined(__MINGW32__) || defined (_WIN32)
      Sleep(1);
    #else
      sleep(1);
    #endif
//    sleep(1);

  }

  if (fsh.age_flags(123))
  {
    dvariable t=fsh.age_flags(123)*norm2(2.0+fsh.q0);
    cout << " q0 penalty = " << t << endl;
    xy+=t;
  }
  if (fsh.parest_flags(339)>0)
  {
    xy+=fsh.parest_flags(339)/1000.0*norm2(exp(fsh.fish_pars(14)));
  }
  if (fsh.parest_flags(340)>0)
  {
    xy+=fsh.parest_flags(340)*norm2(fsh.fish_pars(16));
  }
  if (fsh.parest_flags(341)>0)
  {
    xy+=fsh.parest_flags(341)*norm2(exp(fsh.fish_pars(18)));
  }

  if (fsh.parest_flags(74)>0 || sum(column(fsh.fish_flags,72)))
  {
    xy+=bs_sel_pen(fsh);
  }

  if (fsh.age_flags(170)>0)
  {
    if (fsh.age_flags(175)==0)
    {
      //cout << "Not callingcalculate_the_totalbiomass_depletion()"
      //   << endl << "because af(175) (target biomass) =0" << endl;
    }
    else
    {
      xy+=fsh.calculate_the_totalbiomass_depletion();
    }
  }
  if (fsh.parest_flags(346)==1)
  {
    xy+=fsh.calculate_the_relative_biomass_depletion();
  }
  else if (fsh.parest_flags(346)==2)
  {
    xy+=fsh.calculate_the_average_biomass();
  }
 /*
  if (fsh.age_flags(103>0)
    xy+=fsh.age_flags(103)*norm2(fsh.initpop-mean(fsh.initpop));
  */

  // ad hoc penalty
  // this penalty has now been put with the code for the
  // likelihood contribution for the tag returns
  //if (fsh.parest_flags(110))
  //{
  //  for (int i=30;i<=104;i++)
  //    xy+=100*fsh.parest_flags(110)
  //      *square(fsh.rep_rate(1,i,2)-fsh.rep_rate(1,i,3));
  //}
  if (fsh.parest_flags(169)>0)
  {
    dvariable pen=square(fsh.sv(17));
    cout <<"callpen.cpp " << pen << endl;
    xy+=square(pen);
  }
  // ******************************************************** 
  // ******************************************************** 

  if (fsh.pmsd)
  {
    fsh.pmsd->current_species=1;
  }
  dvar_vector ml=vget_generic_mean_lengths(fsh);
  dvariable xtmp=0.0;
  dvariable ytmp=0.0;
  MY_DOUBLE_TYPE pwght=fsh.parest_flags(182)/10.0;
  if (pwght==0) pwght=.01;
  if (fsh.parest_flags(173))
  {
    int num=fsh.parest_flags(173);
    for (int j=1;j<=num;j++) 
    {
      xtmp+= pwght*square(fsh.age_pars(3,j));
      dvariable diff=ml(j+1)-ml(j);
      if (diff<0.0)
      {
        ytmp+=1e+3*square(diff);
      }
    }
  }
  xy+=xtmp;
  //    cout << "#### Dbug callpen: age_pars3 after xtmp= " << setprecision(15) << xy << endl; //NMD
  xy+=ytmp;
  //    cout << "#### Dbug callpen: age_pars3 after ytmp= " << setprecision(15) << xy << endl; //NMD
  ppf_tmp=value(xtmp+ytmp);  //NMD_21nov2023
  if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_21nov2023
  {
    fsh.ppstf->vonb_devs_pen=ppf_tmp;
    ppf_tmp=0.0;
  }

  if (fsh.pmsd)
  {
    int ns=fsh.pmsd->num_species;
    for (int is=2;is<=ns;is++)
    {
      fsh.pmsd->current_species=is;
      dvar_vector ml=vget_generic_mean_lengths(fsh);
      dvariable xtmp=0.0;
      dvariable ytmp=0.0;
      MY_DOUBLE_TYPE pwght=fsh.parest_flags(182)/10.0;
      if (pwght==0) pwght=.01;
      if (fsh.parest_flags(173))
      {
        int num=fsh.parest_flags(173);
        dvar_vector ap3=fsh.get_age_pars_species(is,3);
        for (int j=1;j<=num;j++) 
        {
	  //          xtmp+= pwght*square(fsh.age_pars(3,j));
          xtmp+= pwght*square(ap3(j));        //NMD 17Aug2012
          dvariable diff=ml(j+1)-ml(j);
          if (diff<0.0)
          {
            ytmp+=1e+3*square(diff);
          }
        }
      }
      xy+=xtmp;
      xy+=ytmp;
    }
  }
  //    cout << "#### Dbug callpen: age_pars3 after xtmp, ytmp spp2 = " << setprecision(12) << xy << endl; //NMD

  if (fsh.parest_flags(175)>0)
  {
    int wt=fsh.parest_flags(176);
    if (wt)
    {
      dvar_vector tmp=fsh.gml/(fsh.gml(fsh.nage+1)-fsh.gml(1));
      xy+=wt*norm2(first_difference(first_difference(first_difference(tmp))));
    }
  }

  if (fsh.age_flags(100)>0)
  {
    dvar_vector& q=fsh.fish_pars(5);
    xy+=10.*norm2(q);
  }
      
  if (fsh.parest_flags(170)>0)
  {
    dvariable pen=square(fsh.sv(18));
    cout <<"callpen.cpp " << pen << endl;
    xy+=square(pen);
  }
 
  if (fsh.parest_flags(171)>0)
    xy+=square(fsh.sv(19));
 
  if (fsh.parest_flags(172)>0)
    xy+=square(fsh.sv(20));

  /*  NMD_12dec2023 - implementation of sv(15) age-specific diffusion obsolete
  if (fsh.parest_flags(173)>0)
  {
    xy+=.1*square(fsh.sv(15));
    if (fsh.pmsd)
    {
      for (int is=2;is<=fsh.pmsd->num_species;is++)
      {
        xy+=.1*square(fsh.pmsd->sv(is,15));
      }
    }
  }
  */
  //    cout << "#### Dbug callpen: age_pars3 after sv15= " << setprecision(12) << xy << endl; //NMD
  // Bayesian prior on albacore diffusion
  // *****************************************************
  // *****************************************************
  if (fsh.age_flags(87))
  {
    xy+=fsh.age_flags(87)*
      square(log(10.*fsh.sv(10)/double(fsh.age_flags(86))));
  }
  char ch;
  dvariable tmp1=0.;
  ivector cutflags=column(fsh.fish_flags,31);  
  if (sum(cutflags))
  {
    fsh.pmsd_error();
    dvar_vector cs=cutoff_sel(fsh,fsh.vb_coff,fsh.var_coff);
    tmp1=5000.*norm2(first_difference(log(1.e-10+cs)));
    xy+=tmp1;
    tmp1=0;
  }
  int i;
  for (i=1;i<=fsh.num_fisheries;i++)  
  {
    if (fsh.fish_flags(i,37)&&fsh.fish_flags(i,33))
    {
      xy+=100.*norm2(fsh.rep_dev_coffs(i));
      for (int nt=2;nt<=fsh.num_real_fish_times(i);nt++)
      {
        int rr=fsh.realization_region(i,nt);
        int rp=fsh.realization_period(i,nt);
        int ri=fsh.realization_incident(i,nt);
        if(fsh.rep_rate(rr,rp,ri)>0.98)
          xy+=250.*square(fsh.rep_rate(rr,rp,ri)-0.98);
      }
    }
  }
  if (fsh.age_flags(72))
  {
    dvar_vector rec(1,fsh.nyears);
    rec.initialize();
    for (int ir=1;ir<=fsh.num_regions;ir++)
    {
      for (int i=1;i<=fsh.nyears;i++)
      {
        rec(i)+=exp(fsh.N(ir,i,1));
      }
    }
    rec=log(rec);
    rec-=mean(rec);
    rec/=norm(rec);
    //cout << " rec = " << rec << endl;
    //cout << " env = " << fsh.rec_covars << endl;
    //exit(1);
    dvariable rpen;
    if (fsh.age_flags(72)>0)
    {
      dvariable rpen= fsh.age_flags(72)/(-100.0) * (1.-square(rec*fsh.rec_covars));
    }
    else
    {
      rpen= 0.5*fsh.nyears*log(1.001-square(rec*fsh.rec_covars));
    }
    //cout << " corr and penalty for environmental effect on recruitment = " 
      //   << value(rec)*value(fsh.rec_covars) << "  " << rpen;
    if (print_switch)
    {
      cout << " penalty for environmental effect on recruitment = " << rpen;
    }
    xy+=rpen;
  }
  ppf_tmp=value(xy);  //NMD_21nov2023
  for (i=1;i<=fsh.num_fisheries;i++)
  {
    if (fsh.fish_flags(i,27))
    {
      xy+=.1*square(fsh.fish_pars(1,i));
    }
  }
  ppf_tmp=value(xy)-ppf_tmp;  //NMD_21nov2023
  if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_21nov2023
  {
    fsh.ppstf->seasonal_q_pen=ppf_tmp;
    ppf_tmp=0.0;
  }

  //  cout << "#### Dbug callpen: ff27 after fish_pars1= " << setprecision(12) << xy << endl; //NMD
  ppf_tmp=value(xy);  //NMD_21nov2023
  if (fsh.age_flags(198)==0)
  {

  // NMD added the following
    int num_tag_fish_groups = fsh.fish_flags(1,34);
    for (i=2 ; i<= fsh.num_fisheries ; i++) 
    {
      int num2=fsh.fish_flags(i,34);
      num_tag_fish_groups=max(num_tag_fish_groups,num2);
    }
    ivector gflag(1,num_tag_fish_groups);
    gflag.initialize();
// NMD added the above

    for (i=1;i<=fsh.num_fisheries;i++)
    {
      if (fsh.fish_flags(i,35))
      {
        // NMD added the following
        if (gflag(fsh.fish_flags(i,34))==0)
        {
          gflag(fsh.fish_flags(i,34))=1;
          MY_DOUBLE_TYPE bf=value(xy);                  //NMD
          // NMD added the above
          xy+= fsh.fish_flags(i,35)*
            square(fsh.fish_pars(3,i)-fsh.fish_flags(i,36)/100.);
          MY_DOUBLE_TYPE af=value(xy);                 //NMD
        }
      }
    }
  }
  else 
  {
    int num_tag_fish_groups=max(fsh.tag_fish_rep_group_flags);
    if (num_tag_fish_groups)
    {
      ivector gflag(1,num_tag_fish_groups);
      gflag.initialize();
      for (int i=1;i<=fsh.num_tag_releases+1;i++)
      {
        for (int j=1;j<=fsh.num_fisheries;j++)
        {
          if (gflag(fsh.tag_fish_rep_group_flags(i,j))==0)
          {
            if (fsh.tag_fish_rep_penalty(i,j)>0)
            {
              gflag(fsh.tag_fish_rep_group_flags(i,j))=1;
              MY_DOUBLE_TYPE bf=value(xy);
              xy+= fsh.tag_fish_rep_penalty(i,j)*
                square(fsh.tag_fish_rep(i,j)
                  -fsh.tag_fish_rep_target(i,j)/100.);
              MY_DOUBLE_TYPE af=value(xy);
            }
          }
        }
      }
    }
    else    // just in cvase there is no grouping
    {
      for (int i=1;i<=fsh.num_tag_releases+1;i++)
      {
        for (int j=1;j<=fsh.num_fisheries;j++)
        {
          if (fsh.tag_fish_rep_penalty(i,j)>0)
          {
            xy+= fsh.tag_fish_rep_penalty(i,j)*
              square(fsh.tag_fish_rep(i,j)
                -fsh.tag_fish_rep_target(i,j)/100.);
          }
        }
      }
    }
  }
  ppf_tmp=value(xy)-ppf_tmp;  //NMD_21nov2023
  if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_21nov2023
  {
    fsh.ppstf->tag_rep_rate_pen=ppf_tmp;
    ppf_tmp=0.0;
  }

  //  cout << "#### Dbug callpen: after tagRRpen= " << setprecision(12) << xy << endl; //NMD
  if (fsh.age_flags(74))
  {
    xy+=exploitation_penalty(fsh);
  }

  ppf_tmp=value(xy);  //NMD_21nov2023
  MY_DOUBLE_TYPE pen=0.1;
  if (fsh.age_flags(110)> 0) pen=fsh.age_flags(110)/10.;
  if (fsh.age_flags(70))
  {
    // put these in no penalties
    //xy+=norm2(fsh.region_rec_diff_sums);

    //xy+=10.0*norm2(fsh.region_rec_diff_colmeans);

    xy+=pen*norm2(fsh.region_rec_diffs);
  }
  ppf_tmp=value(xy)-ppf_tmp;  //NMD_21nov2023
  if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_21nov2023
  {
    fsh.ppstf->region_recr_dev_pen=ppf_tmp;
    ppf_tmp=0.0;
  }

  if (sum(fsh.region_flags(2))> 0)    // JH 14/6/03 This is to penalise 
  {                                   // region_rec_diffs in particular regions
    MY_DOUBLE_TYPE penr=0.0;
    for (int ir=1;ir<=fsh.num_regions;ir++)
    {
      if (fsh.region_flags(2,ir)> 0) 
      {
        int iymax=fsh.nyears;
        if(fsh.parest_flags(183)) iymax=fsh.parest_flags(183);
        penr=double(fsh.region_flags(2,ir))/10.;
        for (int iy=1;iy<=iymax;iy++)
        {
          xy+=penr*square(fsh.region_rec_diffs(iy,ir));
        }
      }
    }
  }

  if (fsh.age_flags(67))
  {
    xy+=10.*norm2(first_difference(first_difference(fsh.age_pars(1))));
    xy+=.1*norm2(fsh.age_pars(1));
  }

  if (fsh.age_flags(73))
  {
    /*
    old code
    xy+=25.*norm2(first_difference(first_difference(fsh.age_pars(2))));
    xy+=5.*norm2(first_difference(fsh.age_pars(2)));
    dvariable mn=mean(fsh.age_pars(2));
    xy+=10.*norm2(fsh.age_pars(2)-mn);
    xy+=10.0*square(mn);
    */
    // newcode
    int wt1=25;
    int wt2=5;
    int wt3=10.;
    int wt4=10.0;
    if (fsh.age_flags(77)) wt1= fsh.age_flags(77);
    if (fsh.age_flags(78)) wt2= fsh.age_flags(78);
    if (fsh.age_flags(79)) wt3= fsh.age_flags(79);
    if (fsh.age_flags(80)) wt4= fsh.age_flags(80);
    /*
    xy+=wt1*norm2(first_difference(first_difference(fsh.age_pars(2))));
    xy+=wt2*norm2(first_difference(fsh.age_pars(2)));
    dvariable mn=mean(fsh.age_pars(2));
    xy+=wt3*norm2(fsh.age_pars(2)-mn);
    xy+=wt4*square(mn);
    */
    // Dave Fournier Nov 7 01
    // This should be the new code to handle the option of having the
    //  nat mort deviations the same for termnal age classes
    //int num=nage-fsh.age_pars(81)+1; JohnH Nov 19 01
    // JohnH Dec 17 01 - stops array bound being exceeded when af(81)=0
    int af81=1;
    if (fsh.age_flags(81)>0) af81=fsh.age_flags(81);
    int num=fsh.nage-af81+1;
 //   int num=fsh.nage-fsh.age_flags(81)+1;
    if (num<=1)
    {
      xy+=wt1*norm2(first_difference(first_difference(fsh.age_pars(2))));
      xy+=wt2*norm2(first_difference(fsh.age_pars(2)));
      dvariable mn=mean(fsh.age_pars(2));
      xy+=wt3*norm2(fsh.age_pars(2)-mn);
      xy+=wt4*square(mn);
    }
    else
    {
      const dvar_vector& ap2=fsh.age_pars(2)(1,num);
      xy+=wt1*norm2(first_difference(first_difference(ap2)));
      xy+=wt2*norm2(first_difference(ap2));
      dvariable mn=mean(ap2);
      xy+=wt3*norm2(ap2-mn);
      xy+=wt4*square(mn);
     }
  }

  if (fsh.age_flags(64))
  {
	 xy+=fsh.sv(14)*fsh.sv(14);
  }

  if (fsh.age_flags(66))
  {
	 xy+=fsh.sv(15)*fsh.sv(15);
  }

  if ( fsh.age_flags(63)>0 )
  {
    xy+=fsh.sv(13)*fsh.sv(13);
  }

  if ( fsh.age_flags(89)>0 )
  {
    if (fsh.age_flags(28)==0)
    {
      xy+=5.*norm2(fsh.diff_coffs2);
    }
    else if (fsh.age_flags(28)>0)
    {
      xy+=fsh.age_flags(28)/10.*norm2(fsh.diff_coffs2);
    }
    else if (fsh.age_flags(28)<0)
    {
      xy-=fsh.age_flags(28)/10.0*norm2(fsh.diff_coffs2-fsh.diff_coffs2_prior);
    }
  }

  if ( fsh.age_flags(91)>0 )
  {
    if (fsh.age_flags(29)==0)
    {
      xy+=5.*norm2(fsh.diff_coffs3);
    }
    else if (fsh.age_flags(29)>0)
    {
      xy+=fsh.age_flags(29)/10.*norm2(fsh.diff_coffs3);
    }
    else if (fsh.age_flags(29)<0)
    {
      xy-=fsh.age_flags(29)/10.0*norm2(fsh.diff_coffs3-fsh.diff_coffs3_prior);
    }
  }

  ppf_tmp=value(xy);  //NMD_21nov2023
  if (!fsh.parest_flags(155))
  {
    if (sum(fsh.region_flags(1)))
    {
      if (!fsh.pmsd)      // XXXXYYY
      {
        xy+=100.*square(sum(fsh.region_pars(1))-1.);
      }
      else
      {
        for (int is=1;is<=fsh.pmsd->num_species;is++)
        {
          int mmin=fsh.pmsd->region_bounds(is,1); 
          int mmax=fsh.pmsd->region_bounds(is,2); 
          xy+=100.*square(sum(fsh.region_pars(1)(mmin,mmax))-1.);
        }
      }
    }
  }
  else
  {
    dvariable tmp=0.0;
    if (!fsh.pmsd)
    {
      double sum_reg=value(sum(fsh.region_pars(1)));
      if (fsh.parest_flags(216)>0 && sum_reg>0.0 && fsh.num_regions>1)
      {
        MY_DOUBLE_TYPE pen=0.1;
        if (fsh.age_flags(110)> 0) pen=fsh.age_flags(110)/10.;
        dvar_matrix reg_recr(1,fsh.num_regions,1,fsh.nyears);
        dvar_vector mnreg_recr(1,fsh.num_regions);
        dvar_vector preg_recr(1,fsh.num_regions);
        dvariable tot_recr=0.0;
        for (int ir=1;ir<=fsh.num_regions;ir++)
        {
          for (int iy=1;iy<=fsh.nyears;iy++)
          {
            reg_recr(ir,iy)=exp(fsh.N(ir,iy,1));
          }
        }
        for (int ir=1;ir<=fsh.num_regions;ir++)
        {
          mnreg_recr(ir)=mean(reg_recr(ir));
        }
        tot_recr=sum(mnreg_recr);
        for (int ir=1;ir<=fsh.num_regions;ir++)
        {
          preg_recr(ir)=mnreg_recr(ir)/tot_recr;
        }
// now calculate the normal prior penalty
        dvar_vector  tmpdif=preg_recr-fsh.region_pars(1);
        tmp=pen*norm2(tmpdif);
        cout << " preg_recr 1 = " << preg_recr(1) << "  preg_recr 1 = " << preg_recr(2)  << endl;
        cout << " tmp = " << tmp << "   pen = " << pen << endl;
      }
    }
    else
    {
      int mmin=0;
      int mmax=0;
      for (int is=1;is<=fsh.pmsd->num_species;is++)
      {
        mmin=fsh.pmsd->region_bounds(is,1); 
        mmax=fsh.pmsd->region_bounds(is,2); 
        int nr=mmax-mmin+1;
	dvar_vector reg_pars(1,nr);
	reg_pars=fsh.region_pars(1)(mmin,mmax);
        double sum_reg=value(sum(reg_pars));
        dvar_matrix reg_recr(1,nr,1,fsh.nyears);
        dvar_vector mnreg_recr(1,nr);
        dvar_vector preg_recr(1,nr);
        dvariable tot_recr=0.0;
        if (fsh.parest_flags(216)>0 && sum_reg>0.0)
        {
          MY_DOUBLE_TYPE pen=0.1;
          if (fsh.age_flags(110)>0) pen=fsh.age_flags(110)/10.;
          dvar_vector tmpdif(1,fsh.nyears);
	  tmpdif.initialize();
          for (int ir=mmin;ir<=mmax;ir++)
          {
            for (int iy=1;iy<=fsh.nyears;iy++)
            {
              reg_recr(ir,iy)=exp(fsh.N(ir,iy,1));
            }
          }
          for (int ir=mmin;ir<=mmax;ir++)
          {
            mnreg_recr(ir)=mean(reg_recr(ir));
          }
          tot_recr=sum(mnreg_recr);
          for (int ir=mmin;ir<=mmax;ir++)
          {
            preg_recr(ir)=mnreg_recr(ir)/tot_recr;
          }
          tmpdif=preg_recr-reg_pars;
          tmp+=pen*norm2(tmpdif);
        }
      }
    }
    xy+=tmp;
  }
  ppf_tmp=value(xy)-ppf_tmp;  //NMD_21nov2023
  if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_21nov2023
  {
    fsh.ppstf->mean_region_recr_pen=ppf_tmp;
    ppf_tmp=0.0;
  }

  for (i=1;i<=fsh.num_fisheries;i++)
  {
    if (fsh.fish_flags(i,19)>0)
    {
      MY_DOUBLE_TYPE pen=1;
      if (fsh.fish_flags(i,20)) pen=fsh.fish_flags(i,20)/10.;
      int num_ages=min(fsh.fish_flags(i,19),fsh.nage);
      for (int j=1;j<=fsh.num_real_fish_times(i);j++)
      {
        int rp=fsh.realization_period(i,j);
        int rr=fsh.realization_region(i,j);
        int ri=fsh.realization_incident(i,j);
        if (fsh.len_sample_size(rr,rp,ri)>0 || 
            fsh.wght_sample_size(rr,rp,ri)>0)
        {
          for (int jj=1;jj<=num_ages;jj++)
          {
            xy+=pen*square(fsh.sel_dev_coffs(i,j,jj));
          }
        }
      }
    }
  }
  //  cout << "#### Dbug callpen: after fish_flags19= " << setprecision(12) << xy << endl; //NMD
  // cout << " N1" << flush;
  if (print_switch)
  {
	 //cout << "after sel devs penalty xy = "<< xy << endl;
  }

  // cout << " N2" << flush;
  if (fsh.parest_flags(163)>0)
  {
	 MY_DOUBLE_TYPE wt=1.0;
	 if (fsh.parest_flags(164)>0)
	 {
		wt=double(fsh.parest_flags(164));
	 }
	 xy+=wt*norm2(fsh.growth_dev);
  }
  // cout << " N3" << flush;
  if (fsh.parest_flags(166)>0)
  {
    MY_DOUBLE_TYPE wt=(MY_DOUBLE_TYPE)( fsh.parest_flags(166))/10.0;
    {
      xy+=wt*square(fsh.var_coff(1));
    }
  }
  if (fsh.parest_flags(167)>0)
  {
    MY_DOUBLE_TYPE wt=(MY_DOUBLE_TYPE)( fsh.parest_flags(167))/10.0;
    {
      xy+=wt*square(fsh.var_coff(2));
    }
  }
  // ************************************************************
  // ************************************************************
  // ************************************************************


 //cout << "pen = " << xy << endl;
 //cout << "pen = " << xy << endl;
  
  MY_DOUBLE_TYPE risp_wt=1000;
  if (fsh.parest_flags(75))
  {
    if (fsh.parest_flags(75)>0)
      risp_wt= 100*fsh.parest_flags(75);
    else
      risp_wt=0;
  }
  dvariable risp=risp_wt*recr_initpop_sum_penalty(fsh,print_switch);
  xy+=risp;
  ppf_tmp=value(risp);  //NMD_21nov2023
  if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_21nov2023
  {
    fsh.ppstf->init_recr_pen=ppf_tmp;
    ppf_tmp=0.0;
  }

  if (print_switch)
  {
    //cout << "recr initpop sum penalty penalty = "<< risp << endl;
  }
  {
    if (fsh.nage>=3)
    {
      dvariable tmp3=0.0;
      {
        MY_DOUBLE_TYPE wt=0.5;
        MY_DOUBLE_TYPE pen_wght=.1;
        if (fsh.age_flags(76)>0) 
        {
          pen_wght=fsh.age_flags(76);
        }
        for (int ir=1;ir<=fsh.num_regions;ir++)
        {
          tmp3+=pen_wght*norm2(first_difference(first_difference(fsh.N(ir,1))));
        }
      }
      xy+=tmp3;
      cout << "fishery initial population curvature penalty = "<< tmp3 << endl;
      ppf_tmp=value(tmp3);  //NMD_21nov2023
    }
    if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_21nov2023
    {
      fsh.ppstf->init_agecomp_curve_pen=ppf_tmp;
      ppf_tmp=0.0;
    }
  }
  if (fsh.age_flags(41) || fsh.pmsd)
  {
    int ns=1;
    if (fsh.pmsd)
      ns=fsh.pmsd->num_species;
    MY_DOUBLE_TYPE recr_weight=1.e-6;
    if (fsh.parest_flags(149))
    {
      recr_weight=fsh.parest_flags(149)/10.;
    }
    if (!fsh.pmsd)
    {
      int i=1;
      xy+=recrpen(fsh,i,print_switch,recr_weight);
    }
    else
    {
      if (!sum(column(fsh.pmsd->species_flags,2)))
      {
        for (int i=1;i<=ns;i++)
        {
          xy+=recrpen(fsh,i,print_switch,recr_weight);
        }
      }
      else
      {
        ivector sf2=column(fsh.pmsd->species_flags,2);
        if (ns !=2 || sum(sf2)!=1 )
        {
          cerr << "only works at present for 2 species model" << endl;
          ad_exit(1);
        }
        for (int i=1;i<=ns;i++)
        {
          if (sf2(i))
          {
            xy+=recrpen(fsh,i,print_switch,recr_weight);
          }
        }
      }
    }

    if (fsh.age_flags(129)==1)
    {
      if (fsh.age_flags(135)==1)
      {
        cerr << " ERROR: age_flags 135 and 129 are incompatible" << endl;
        ad_exit(1);
      }
      for (int is=1;is<=ns;is++)
      {
        MY_DOUBLE_TYPE wght=1.0;
        xy+=recruitment_autocorrelation_moment_estimator(fsh,is,print_switch,wght);
      }
    }
    xy+=recr_weight*norm2(fsh.initpop);
  }
  else
  {
    fsh.pmsd_error();
    dvar_vector ttmprecr=fsh.tmprecr(1,fsh.last_real_year);
    dvariable tmp4=.01*norm2(first_difference(ttmprecr));
    if (print_switch)
    {
      cout << tmp4 << " for first diff tmprecr"<< endl;
    }
    MY_DOUBLE_TYPE tmprecr_weight=1.e-6;
    if (fsh.parest_flags(149))
    {
      tmprecr_weight=fsh.parest_flags(149);
    }
    dvariable tmp_stuff;
    tmp_stuff=tmprecr_weight*norm2(fsh.tmprecr-mean(ttmprecr));
    if (print_switch)
    {
     cout << "norm2(recr-mean) penalty"  << tmp_stuff << endl;
    }
    tmp4+=tmp_stuff;

    tmp_stuff=tmprecr_weight*norm2(fsh.tmpinitpop);
    tmp4+=tmp_stuff;

    dvector trend(ttmprecr.indexmin(),ttmprecr.indexmax());
    trend.fill_seqadd(-1.,2./(trend.size()-1.));

    if (fsh.parest_flags(153)>=0)
    {
      MY_DOUBLE_TYPE trend_wt=.01;
      if (fsh.parest_flags(153)!=0)
      {
        trend_wt=fsh.parest_flags(153)/10.;
      }
      dvariable tmp6=trend_wt*square(ttmprecr*trend);
  
      if (print_switch)
      {
        cout << "fishery recruitment trend penalty = "<< setprecision(12) << tmp6 << endl;//NMD
      }
      xy+=tmp6;
      xy+=tmp4;
      if (print_switch)
      {
        cout << " after recruitment trend penalty xy = " << xy << endl;
      }
    }
    dvariable r2pen=.01*norm2(first_difference(first_difference(ttmprecr)));
    if (print_switch)
    {
      cout << "fishery recruitment curv penalty = "<< r2pen << endl;
    }
    xy+=r2pen;
    if (fsh.parest_flags(148))
    {
      dvariable last_wt=fsh.parest_flags(148)/10.0;
      dvariable rlast_pen=last_wt*square(fsh.tmprecr(fsh.last_real_year)-
                          fsh.tmprecr(fsh.last_real_year-1));
      xy+=rlast_pen;
      cout << "Last recruit. penalty = " << rlast_pen << endl; 
    }
  }

  if (fsh.age_flags(92)==0)
  {
    xy+=eff_dev_penalty(fsh,print_switch);
  }

  if (fsh.age_flags(92)==0)   // DF  apr16 07
  {
    ivector grouping=column(fsh.fish_flags,29);
    const ivector& xcatflags=(column(fsh.fish_flags,15));  
    ivector& catflags=(ivector&) xcatflags;  
    
    if (!sum(grouping))
    {
      dvar_vector tt(1,fsh.num_fisheries);
      if (!sum(catflags))
      {
        for (i=1;i<=fsh.num_fisheries;i++)
        {
          //tmp1=50.*norm2(elem_prod(sqr(fsh.between_times),
          //  fsh.catch_dev_coffs));
          tt(i)=50.*norm2(elem_prod(sqr(fsh.between_times(i)),
            fsh.catch_dev_coffs(i)));
        }
      }
      else
      {
        for (i=1;i<=fsh.num_fisheries;i++)
        {
          if (catflags(i)==0) catflags(i)=50;
          tt(i)=catflags(i)*
            norm2(elem_prod(sqr(fsh.between_times(i)),
              fsh.catch_dev_coffs(i)));
        }
      }
      tmp1=sum(tt);
      if (fsh.ppstf)
      {
        fsh.ppstf->catchability_dev_penalty_by_fishery=value(tt);
      }
    }
    else
    {
      dvar_vector tt(1,fsh.ngroups);
      tt.initialize();
      for (i=1;i<=fsh.ngroups;i++)
      {
        int i1=fsh.gfish_index(i,1);
        if (catflags(i1)==0)catflags(i1)=50;
        tt(i)+=catflags(i1)*
          norm2(elem_prod(sqr(fsh.grouped_between_times(i)),
          fsh.grouped_catch_dev_coffs(i)));
      }
      if (fsh.ppstf)
      {
        fsh.ppstf->catchability_dev_penalty_by_group=value(tt);
      }
      tmp1=sum(tt);
    }
    xy+=tmp1;
    tmp1=0.0;
  }

  xy+=first_length_bias_penalty(fsh,print_switch);

  ppf_tmp=value(xy);  //NMD_27dec2023
  if (!fsh.parest_flags(143))
  {
    if (fsh.nfmbound>0)
    {
      mean_length_constraints=fsh.fmeanpen();
      xy+=mean_length_constraints;
    }
  }
  ppf_tmp=value(xy)-ppf_tmp;  //NMD_27dec2023
  if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_27dec2023
  {
    fsh.ppstf->mean_len_constr_pen=ppf_tmp;
    ppf_tmp=0.0;
  }


  ppf_tmp=value(xy);  //NMD_22nov2023
  xy+=incident_sel_curvature_penalty(fsh,print_switch);
  ppf_tmp=value(xy)-ppf_tmp;  //NMD_22nov2023
  if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_22nov2023
  {
    fsh.ppstf->incid_sel_curve_pen=ppf_tmp;
    ppf_tmp=0.0;
  }
  // now only for ageselmean
  ppf_tmp=value(xy);  //NMD_22nov2023
  xy+=selmean_penalty(fsh,print_switch);
  ppf_tmp=value(xy)-ppf_tmp;  //NMD_22nov2023
  if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_22nov2023
  {
    fsh.ppstf->mean_sel_pen=ppf_tmp;
    ppf_tmp=0.0;
  }
  
  ppf_tmp=value(xy);  //NMD_22nov2023
  xy+=bs_selmean_penalty(fsh,print_switch);
  ppf_tmp=value(xy)-ppf_tmp;  //NMD_22nov2023
  if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_22nov2023
  {
    fsh.ppstf->time_block_sel_mean_pen=ppf_tmp;
    ppf_tmp=0.0;
  }

  ppf_tmp=value(xy);  //NMD_22nov2023
  xy+=selectivity_deviations_penalty(fsh,print_switch);
  ppf_tmp=value(xy)-ppf_tmp;  //NMD_22nov2023
  if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_22nov2023
  {
    fsh.ppstf->sel_devs_pen=ppf_tmp;
    ppf_tmp=0.0;
  }

  tmp1=0.0;

  ivector ff57=column(fsh.fish_flags,57);
  for (i=1;i<=fsh.num_fisheries;i++)
  {
    // don't do this if we have cubic splines
    if (ff57(i) !=3)
    {
// NMD_17Mar2015
      ppf_tmp=value(xy);  //NMD_22nov2023
      int jmin=fsh.bs_selcoff(i).indexmin();
      int jmax=fsh.bs_selcoff(i).indexmax();
      for (int j=jmin;j<=jmax;j++)
      {
        int kmin=fsh.bs_selcoff(i,j).indexmin();
        int kmax=fsh.bs_selcoff(i,j).indexmax();
        for (int k=kmin;k<=kmax;k++)
        {
          tmp1=.0001*norm2(fsh.bs_selcoff(i,j,k));
          if (print_switch)
          {
            cout << "fishery selectivity norm penalty = "<< tmp1 << endl;
          }
          xy+=tmp1;
          int sz=fsh.bs_selcoff(i,j,k).indexmax()-
	    fsh.bs_selcoff(i,j,k).indexmin()+1;
          tmp1=0.0;
          if (sz>2)
          {
            tmp1+=.1*norm2(first_difference(first_difference(fsh.bs_selcoff(i,j,k))));
          }
          if (print_switch)
          {
           //cout << "fishery selectivity curvature penalty = "<< tmp1 << endl;
          }
          xy+=tmp1;
          ppf_tmp=value(xy)-ppf_tmp;  //NMD_22nov2023
          if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_22nov2023
          {
            fsh.ppstf->nonspline_sel_curve_pen(i)=value(ppf_tmp);
            ppf_tmp=0.0;
          }
        }
      }
      
// NMD_17Mar2015
    }
    else
    {
      ppf_tmp=value(xy);  //NMD_22nov2023
      // add small penalty to keep selcoffs from getting stuck at -20.
      int jmin=fsh.bs_selcoff(i).indexmin();
      int jmax=fsh.bs_selcoff(i).indexmax();
      
      MY_DOUBLE_TYPE spen=0.0;  //NMD_27jul2021
      if (fsh.parest_flags(359)>0)
      {
        spen=fsh.parest_flags(359)/10000.0;  //NMD_27jul2021
        for (int j=jmin;j<=jmax;j++)
        {
          int kmin=fsh.bs_selcoff(i,j).indexmin();
          int kmax=fsh.bs_selcoff(i,j).indexmax();
          for (int k=kmin;k<=kmax;k++)
          {
            int lmin=fsh.bs_selcoff(i,j,k).indexmin();
            int lmax=fsh.bs_selcoff(i,j,k).indexmax();
            for (int l=lmin;l<=lmax;l++)
            {
              if (fsh.bs_selcoff(i,j,k,l)<=-15.0)
              {
                dvariable tmp=spen*square(fsh.bs_selcoff(i,j,k,l)+15.0);
                xy+=tmp;
              }
            }
          }
        }
      }
      ppf_tmp=value(xy)-ppf_tmp;  //NMD_22nov2023
      if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_22nov2023
      {
        fsh.ppstf->spline_sel_bound_pen(i)=value(ppf_tmp);
        ppf_tmp=0.0;
      }
    }
  }

  dvariable penx=selectivity_form_penalty(fsh,print_switch);
  if (print_switch)
  {
    if(pof_pen)
    {
      (*pof_pen) << " selectivity_form_penalty  " << penx << endl;
    }
  }
  xy+=penx;
  ppf_tmp=value(penx);  //NMD_22nov2023
  if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_22nov2023
  {
    fsh.ppstf->sel_form_pen=ppf_tmp;
    ppf_tmp=0.0;
  }

  pen=0.0;
  dvariable fmpen=0.;
  for (i=1;i<=fsh.num_fisheries;i++)
  {
    if (fsh.fish_flags(i,14)>0)
    {
      MY_DOUBLE_TYPE max_fm=log(fsh.fish_flags(i,14)/10.);
      MY_DOUBLE_TYPE max_fm1=max_fm+.01;
      for (int j=1;j<=fsh.num_real_fish_times(i);j++)
      {
        int rr=fsh.realization_region(i,j);
        int rp=fsh.realization_period(i,j);
        int ri=fsh.realization_incident(i,j);
        dvar_vector& fm=fsh.fish_mort(rr,rp,ri);
        for (int jj=1;jj<=fsh.nage;jj++)
        {
          if (fm(jj)>=max_fm)
          {
            fmpen+=1000.*square(fm(jj)- max_fm);
            if (fm(jj)>=max_fm1)
            {
              fmpen+=10000.*square(fm(jj)- max_fm1);
            }
          }
        }
      }
    }
  }
  xy+=fmpen;
  if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_11dec2023
  {
    fsh.ppstf->fmort_max_pen=fmpen;
    ppf_tmp=0.0;
  }
  
  fmpen=0.0;  //NMD_30nov2023
  for (i=1;i<=fsh.num_fisheries;i++)
  {
    if (fsh.fish_flags(i,17)>0)
    {
      int lag=fsh.fish_flags(i,17);
      MY_DOUBLE_TYPE multiplier=log(fsh.fish_flags(i,18)/10.);
      for (int j=1;j<=fsh.num_real_fish_times(i);j++)
      {
        int rr=fsh.realization_region(i,j);
        int rp=fsh.realization_period(i,j);
        int ri=fsh.realization_incident(i,j);
        dvar_vector& fm=fsh.fish_mort(rr,rp,ri);
        if (fm(fsh.nage-lag)<fm(fsh.nage)+multiplier)
        {
          fmpen+=1000.*square(fm(fsh.nage)-fm(fsh.nage-lag)+multiplier);
        }
      }
    }
  }
  xy+=fmpen;
  if (fsh.ppstf && !fsh.af170q0 && print_switch)  //NMD_11dec2023
  {
    fsh.ppstf->fmort_agelag_max_pen=fmpen;
    ppf_tmp=0.0;
  }

  // set target average fishing mortality
  if (fsh.age_flags(37)>0)
  {
  /*  dvariable ssum=0.;
  //  dvariable szsum=0.;
  //  for (int ir=1;ir<=fsh.num_regions;ir++)
  //  {
  //    ssum+=sum(exp(fsh.fish_mort(ir)));
  //    szsum+=size_count(fsh.fish_mort(ir));
  //  }
  //  ssum/=szsum;
  */
  //*************************************
    dvariable av_exploitation;
    dvariable totnum;
    dvariable totcatch;
    totcatch.initialize();
    totnum.initialize();
    int age1=1;
    if (fsh.age_flags(42)>0) age1=fsh.age_flags(42);
    int ir;
    for (ir=1;ir<=fsh.num_regions;ir++)
    {
      for (int ip=1;ip<=fsh.num_real_fish_periods(ir);ip++) 
      {
        if (age1>1)
        {
          for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++) 
          {                                 
            //for (int ia=age1;ia<=fsh.nage;ia++)
            // totcatch is the total catch summed over year and age (from age1)
            totcatch+=sum(exp(fsh.catch(ir,ip,fi)(age1,fsh.nage)));
          }
        }
        else
        {
          if (allocated(fsh.catch(ir,ip)))
            // cout<<"callpen.cpp-913 "<<fsh.catch(ir,ip)<<endl;
            totcatch+=sum(exp(fsh.catch(ir,ip)));
        }
      }
    }
    for (ir=1;ir<=fsh.num_regions;ir++)
    {
      if (age1>1)
      {
        for (int iy=1;iy<=fsh.nyears;iy++)
        {
          totnum+=sum(exp(fsh.N(ir,iy)(age1,fsh.nage)));
        }
      }
      else
      {
        for (int iy=1;iy<=fsh.nyears;iy++)
        {
          totnum+=sum(exp(fsh.N(ir,iy)));
        }
      }
    }


    av_exploitation=totcatch/totnum;
  //*************************************
    cout << "Av. F = " << setprecision(6) << av_exploitation << endl;
    MY_DOUBLE_TYPE wt=1000;
    //double wt=10;
    if (fsh.age_flags(39))
    {
 //     wt=100.*fsh.age_flags(39);
      wt=1.*fsh.age_flags(39);
    }
    MY_DOUBLE_TYPE mult=double(fsh.num_fish_data_recs)/double(fsh.nyears);
  //  dvariable tmp2=wt*square(log(ssum/(fsh.age_flags(37)/(1000.*mult))));
    dvariable tmp2=wt*square(log(av_exploitation/(fsh.age_flags(37)/1000.)));
    xy+=tmp2;
    cout << "The fish mort penalty is " << tmp2 << endl;
    if (print_switch)
    {
 //cout << "The fish mort penalty is " << tmp2 << endl;
 //cout << "average fish mort is " << ssum << endl;
 //cout << " after fish mort penalty xy = " << xy << endl;
    }
  }

  if (fsh.age_flags(43)>0)
  {
    dvariable ssum=0.;
    for (int ir=1;ir<=fsh.num_regions;ir++)
    {
      for (int ip=1;ip<=fsh.num_fish_periods(ir);ip++) 
      {
        if (fsh.year(ir,ip) >= (fsh.nyears-fsh.age_flags(43)+1) )
        {
          for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
	  {
            ssum+=sum(exp(fsh.fish_mort(ir,ip,fi)));
          }
        }
      }
    }
    ssum/=double(fsh.age_flags(43));
    MY_DOUBLE_TYPE wt;
    if (fsh.age_flags(45) > 0)
    {
      wt=fsh.age_flags(45);
    }
    else
    {
      wt=1000.;
    }
    MY_DOUBLE_TYPE level=fsh.age_flags(44)/100.;
    dvariable last_pen=wt*square(ssum-level);
    xy+=last_pen;
    if (print_switch)
    {
 cout << "Target level average mortality is " << level << endl;
 cout << "The penalty for terminal years fishing mortality is "
	 << last_pen << endl;
 cout << " after terminal years penalty xy = " << xy << endl;
    }
  }


/*
  if (fsh.age_flags(55)>0)
  {
    dvariable tmp2=10.*shmodel_fit(fsh,print_switch);
    if (fsh.age_flags(56)>0)
    {
      tmp2+=4.*square(log(fsh.sv(9)-1.));
    }
    xy+=tmp2;
  }
*/

  
  if (fsh.age_flags(47)>0)
  {
    xy+=10.*square((*psv)(1));
  }

  if (fsh.age_flags(54)>1)
  {
    xy+=avail_coff_penalty(fsh);
  }

  if (fsh.parest_flags(160))
  {
    xy+=.1*fsh.sv(6)*fsh.sv(6);
  }
  if (fsh.age_flags(58)>0)
  {
    xy+=0.1*fsh.sv(10)*fsh.sv(10);
  }
  if (fsh.age_flags(59)>0)
  {
    xy+=0.1*fsh.sv(11)*fsh.sv(11);
    xy+=0.1*fsh.sv(12)*fsh.sv(12);
  }

  if (fsh.age_flags(82)>0 )
  {
    MY_DOUBLE_TYPE tnm=fsh.age_flags(82)/100.0;

    int mmin=1;
    int mmax=fsh.nage;
    if(fsh.age_flags(83)>0) mmin=fsh.age_flags(83);
    if(fsh.age_flags(85)>0) mmax=fsh.age_flags(85);
    dvariable avgmort=mean(exp(fsh.nat_mort(1)(mmin,mmax)));
    dvariable mpen=fsh.age_flags(84)*square(log((1.e-10+avgmort)/tnm));
    cout << "avg nat mort " << avgmort << " target avgnm " << tnm
         << " penalty " << mpen << endl;
    xy+=mpen;
  }
  if (fsh.pmsd)
  {
    for (int is=2;is<=fsh.pmsd->num_species;is++)
    {
      ivector & af=fsh.pmsd->age_flags(is);
      if (af(82))
      {
        MY_DOUBLE_TYPE tnm=af(82)/100.0;

        int mmin=1;
        int mmax=fsh.nage;
        if(af(83)>0) mmin=af(83);
        if(af(85)>0) mmax=af(85);
        dvariable avgmort=mean(exp(fsh.get_nat_mort_species(is)(1)(mmin,mmax)));
        dvariable mpen=af(84)*square(log((1.e-10+avgmort)/tnm));
        cout << "avg nat mort " << avgmort << " target avgnm " << tnm
             << " penalty " << mpen << endl;
        xy+=mpen;
      }
    }
  }
  if (fsh.age_flags(97)>0)
  {
    dvariable bratio=ratio_first_last_biocalc(fsh);
    dvariable tmp=bratio/(fsh.age_flags(97)/100.0);
    dvariable wt=1000.0;
    if (fsh.age_flags(98)>0) wt=float(fsh.age_flags(98));
    tmp=wt*square(log(tmp));
    xy+=tmp;
    cout << "Target biomass ratio" 
      << setfixed() << setprecision(4) << setw(8) 
      << (fsh.age_flags(97)/100.0)
      << " Actual biomass ratio " 
      << setfixed() << setprecision(4) << setw(8) 
      << bratio
      << " Penalty " 
      << setfixed() << setprecision(2) << setw(8) 
      << tmp
      << endl;
  }
  
  if (fsh.age_flags(130)>0)
  {
    fsh.pmsd_error();
    int a1=1;
    if(fsh.age_flags(131)>0) a1=fsh.age_flags(131);
    int a2=fsh.nage;
    if(fsh.age_flags(132)>0) a2=fsh.age_flags(132);
    dvariable mkratio=mean(exp(fsh.nat_mort(1)(a1,a2)))/fsh.vb_coff(3);
    dvariable tmp=mkratio/(fsh.age_flags(130)/100.0);
    dvariable wt=fsh.age_flags(133);
    if (wt < 0.0) wt=0.0;
    tmp=wt*square(log(tmp+1.e-10));
    xy+=tmp;
    cout << "Target M/K ratio" 
      << setfixed() << setprecision(4) << setw(8) 
      << (fsh.age_flags(130)/100.0)
      << " Actual M/K ratio " 
      << setfixed() << setprecision(4) << setw(8) 
      << mkratio
      << " Penalty " 
      << setfixed() << setprecision(2) << setw(8) 
      << tmp
      << endl;
  }
  
  if (fsh.age_flags(107)>0)
  {	 
    MY_DOUBLE_TYPE targetexp=fsh.age_flags(108)/100.0;
    dvariable exprate=calculate_overall_exploitation(fsh);
    dvariable pen=double(fsh.age_flags(107))*square(log(exprate/targetexp));
    cout <<"callpen.cpp "   << setfixed() << setprecision(5) << setw(10) 
	     << "exprate = " << exprate << " target = " << targetexp
         << " penalty = " << pen << endl;	    
    xy+=pen;
  }  

  xy+=fsh.normalize_seasonal_catchability();

  xy+=fsh.recinpop_orth_penalties();

  return xy;
}

dvariable dvar_fish_stock_history::calculate_the_totalbiomass_depletion
  (void)
{
  int i,ir;
  int startyr=nyears-age_flags(173)+1;
  int endyr=nyears-age_flags(174)+1;
  startyr=min(startyr,nyears);
  endyr=min(endyr,nyears);
  dvariable bio=0.0;
  dvariable bio_q0=0.0;
  for (i=startyr;i<=endyr;i++)
  {
    for (ir=1;ir<=num_regions;ir++)
    {
      if (age_flags(172)==0)
      {
        bio+=(exp(N(ir,i))*mean_weight_yr(ir,i))/1000.;
        bio_q0+=(exp(N_q0(ir,i))*mean_weight_yr(ir,i))/1000.;
      }
      else
      {
        bio+=(pmature*(elem_prod(exp(N(ir,i)),
           mean_weight_yr(ir,i))))/1000.;
        bio_q0+=(pmature*(elem_prod(exp(N_q0(ir,i)),
           mean_weight_yr(ir,i))))/1000.;
      }
    }
  } 
  dvariable depletion=bio/bio_q0;
  MY_DOUBLE_TYPE target=age_flags(175)/1000.;
  MY_DOUBLE_TYPE penwt=(MY_DOUBLE_TYPE)age_flags(176);
  dvariable pen=penwt*square(log(depletion/target));
  cout << "Target depletion "
    << setfixed() << setprecision(4) << setw(8)
    << target
    << "  Actual depletion "
    << setfixed() << setprecision(4) << setw(8)
    << depletion
    << "  Penalty "
    << setfixed() << setprecision(4) << setw(8)
    << pen
    << endl;
  return pen;
}

// PK 6-28-05
dvector dvar_fish_stock_history::totalbiomass_depletion_for_plot_file
  (void)
{
  int i,ir;
  int startyr=nyears-age_flags(173)+1;
  int endyr=nyears-age_flags(174)+1;
  startyr=min(startyr,nyears);
  endyr=min(endyr,nyears);
  dvariable bio=0.0;
  dvariable bio_q0=0.0;
  for (i=startyr;i<=endyr;i++)
  {
    for (ir=1;ir<=num_regions;ir++)
    {
      if (age_flags(172)==0)
      {
        bio+=(exp(N(ir,i))*mean_weight_yr(ir,i))/1000.;
        bio_q0+=(exp(N_q0(ir,i))*mean_weight_yr(ir,i))/1000.;
      }
      else
      {
        bio+=(pmature*(elem_prod(exp(N(ir,i)),
           mean_weight_yr(ir,i))))/1000.;
        bio_q0+=(pmature*(elem_prod(exp(N_q0(ir,i)),
           mean_weight_yr(ir,i))))/1000.;
      }
    }
  }
  dvariable depletion=bio/bio_q0;
  MY_DOUBLE_TYPE target=age_flags(175)/1000.;
  MY_DOUBLE_TYPE penwt=(MY_DOUBLE_TYPE)age_flags(176);
  dvariable pen=penwt*square(log(depletion/target));
  dvector ret(1,4);
  ret(1)=target;
  ret(2)=value(depletion);
  ret(3)=penwt;
  ret(4)=value(pen);
  return ret;
}

void  dvar_fish_stock_history::pmsd_error(void)
{
  if (pmsd)
  {
    cerr << "this code not yet modified for mult-species" << endl;
    ad_exit(1);
  }
}


void  dvar_fish_stock_history::indepvars_error(void)
{
  if (parest_flags(246))
  {
    cerr << "this code not yet developed for the "
         << "independent variables report"  << endl;
    ad_exit(1);
  }
}

dvariable dvar_fish_stock_history::calculate_the_relative_biomass_depletion
  (void)
{
  int i,ir;
  int startyr=nyears-age_flags(173)+1;
  int endyr=nyears-age_flags(174)+1;
  startyr=min(startyr,nyears);
  endyr=min(endyr,nyears);
  dvariable init_bio=0.0;
  dvariable term_bio=0.0;
  int initl=1;
  int initu=5;
  int terml=nyears-6;
  int termu=nyears;
  if (age_flags(173)) terml=startyr;
  if (age_flags(174)) termu=endyr;
  for (i=initl;i<=initu;i++)
  {
    for (ir=1;ir<=num_regions;ir++)
    {
      if (age_flags(172)==0)
      {
        init_bio+=(exp(N(ir,i))*mean_weight_yr(ir,i))/1000.;
      }
      else
      {
        init_bio+=(pmature*(elem_prod(exp(N(ir,i)),
           mean_weight_yr(ir,i))))/1000.;
      }
    }
  } 
  init_bio/=(initu-initl+1);
  for (i=terml;i<=termu;i++)
  {
    for (ir=1;ir<=num_regions;ir++)
    {
      if (age_flags(172)==0)
      {
        term_bio+=(exp(N(ir,i))*mean_weight_yr(ir,i))/1000.;
      }
      else
      {
        term_bio+=(pmature*(elem_prod(exp(N(ir,i)),
           mean_weight_yr(ir,i))))/1000.;
      }
    }
  } 
  term_bio/=(termu-terml+1);

  dvariable depletion=term_bio/init_bio;
  
  cout << depletion << endl;
  if (parest_flags(347)==0)
  {
    ofstream ofs("relative_depletion");
    ofs << depletion << endl;
    ad_exit(1);
  }
  
  MY_DOUBLE_TYPE penwt=(MY_DOUBLE_TYPE)parest_flags(348);  // aug lag uses 0.5* penwt
  MY_DOUBLE_TYPE target=parest_flags(347)/1000.;
  //other_lagrange_c(1)=log(depletion/target);
  //other_lagrange_mu(1)=2.0*(MY_DOUBLE_TYPE) penwt;
  //dvariable tmp=0.5*other_lagrange_mu(1)*square(other_lagrange_c(1));
  //if (parest_flags(349)>0)  // use aug lagrange otherwise just quad pen
  //{
  //  tmp-=other_lagrange_lambda(1)*other_lagrange_c(1);
  //}

  dvariable pen=penwt*square(log(depletion/target));
    
  cout << "Target depletion "
    << setfixed() << setprecision(4) << setw(8)
    << target
    << "  Actual depletion "
    << setfixed() << setprecision(4) << setw(8)
    << depletion
    << "  Penalty "
    << setfixed() << setprecision(4) << setw(8)
    << pen
    << endl;
  return pen;
}

dvariable dvar_fish_stock_history::calculate_the_average_biomass(void)
{
  int i,ir;
  int startyr;
  int endyr;
  dvariable avg_bio=0.0;
  startyr=1; 
  endyr=nyears;
  if (age_flags(173)) startyr=nyears-age_flags(173)+1;
  if (age_flags(174)) endyr=nyears-age_flags(174)+1;
  startyr=min(startyr,nyears);
  endyr=min(endyr,nyears);

  for (i=startyr;i<=endyr;i++)
  {
    for (ir=1;ir<=num_regions;ir++)
    {
      if (age_flags(172)==0)
      {
        avg_bio+=(exp(N(ir,i))*mean_weight_yr(ir,i))/1000.;
      }
      else
      {
        avg_bio+=(pmature*(elem_prod(exp(N(ir,i)),
           mean_weight_yr(ir,i))))/1000.;
      }
    }
  } 
  avg_bio/=(endyr-startyr+1);

  
  cout << avg_bio << endl;
  if (parest_flags(347)==0)
  {
    ofstream ofs("avg_bio");
    ofs << avg_bio << endl;
    ad_exit(1);
  }
  
  MY_DOUBLE_TYPE penwt=(MY_DOUBLE_TYPE)parest_flags(348);  // aug lag uses 0.5* penwt
  MY_DOUBLE_TYPE target=parest_flags(347);
  other_lagrange_c(1)=log(avg_bio/target);
  other_lagrange_mu(1)=2.0*(MY_DOUBLE_TYPE) penwt;
  dvariable tmp=0.5*other_lagrange_mu(1)*square(other_lagrange_c(1));
  if (parest_flags(349)>0)  // use aug lagrange otherwise just quad pen
  {
    /*
    if (fabs(other_lagrange_lambda(1))<1.e-20 &&
      historical_parest_flags(349)==0)
    {
      other_lagrange_lambda(1)=
        value(other_lagrange_mu(1)*other_lagrange_c(1));
    }
    */
    tmp-=other_lagrange_lambda(1)*other_lagrange_c(1);
  }
  dvariable pen=penwt*square(log(avg_bio/target));

    
  cout << "Target average biomass "
    << setfixed() << setprecision(2) << setw(7)
    << target
    << " Actual average biomass "
    << setfixed() << setprecision(2) << setw(7)
    << avg_bio
    << " Other_lagrange_c(1) "
    << setscientific() << setprecision(3) << setw(8)
    <<  other_lagrange_c(1)
    << " Other lagrange_lambda(1)*other_lagrange_c(1) "
    << setscientific() << setprecision(3) << setw(8)
    <<  other_lagrange_lambda(1)*other_lagrange_c(1)
    << " Penalty "
    << setfixed() << setprecision(2) << setw(7)
    << pen
    << " Penalty weight "
    << setfixed() << setprecision(0) << setw(7)
    << penwt
    << endl;
  return tmp;
}

#undef HOME_VERSION

