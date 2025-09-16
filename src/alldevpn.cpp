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
//#include "newmprot.hpp"
static dvector operator - (const dvector& v)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  dvector tmp(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    tmp(i)=-v(i);
  }
  return tmp;
}
 
dvariable eff_dev_penalty(dvar_fish_stock_history& fsh,int print_switch) 
{
  dvar_vector tmppen(1,fsh.num_fisheries);
  tmppen.initialize();
  for (int fi=1;fi<=fsh.num_fisheries;fi++)
  {
    if (fsh.fish_flags(fi,83)==0)   //NMD_16jun2020
    {
      if (fsh.missing_effort_flag(fi)==0)
      {
        if (fsh.fish_flags(fi,65)==0)
        {
          int nrft=fsh.num_real_fish_times(fi);
          if (fsh.fish_flags(fi,13)==0)
          {
            if (fsh.fish_flags(fi,66)==0)
            {
              MY_DOUBLE_TYPE eff_wt=10.;
              dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
              dvariable eff_dev_pen=-sum(log(1.e-10+exp(-eff_wt*tmp_square)+
                .05*exp(-eff_wt/5.*tmp_square)));
              tmppen(fi)+=eff_dev_pen;
            }
            else
            {
              dvector eff_wt=10.*fsh.effort_weight_by_fishery(fi)(1,nrft);
              dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
              dvariable eff_dev_pen=
                -sum(log(1.e-10+exp(-elem_prod(eff_wt,tmp_square))+
                .05*exp(-elem_prod(eff_wt/5.,tmp_square))));
              tmppen(fi)+=eff_dev_pen;
            }
          }
          else if (fsh.fish_flags(fi,13)>0)
          {
            if (fsh.fish_flags(fi,66)==0)
            {
              MY_DOUBLE_TYPE eff_wt=fsh.fish_flags(fi,13);
              dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
              dvariable eff_dev_pen=-sum(log(1.e-10+exp(-eff_wt*tmp_square)+
                .05*exp(-eff_wt/5.*tmp_square)));
              tmppen(fi)+=eff_dev_pen;
            }
            else
            {
              dvector eff_wt=fsh.fish_flags(fi,13)*
                fsh.effort_weight_by_fishery(fi)(1,nrft);
              dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
              dvariable eff_dev_pen=
                -sum(log(1.e-10+exp(-elem_prod(eff_wt,tmp_square))+
                .05*exp(-elem_prod(eff_wt/5.,tmp_square))));
              tmppen(fi)+=eff_dev_pen;
            }
          }
          else if (fsh.fish_flags(fi,13)<0)
          {
            dvector eff_wt;
            if (fsh.fish_flags(fi,66)==0)
              eff_wt=-fsh.fish_flags(fi,13)*
                sqrt(.01+fsh.effort_by_fishery(fi)(1,nrft));
            else
              eff_wt=-fsh.fish_flags(fi,13)*
                elem_prod(fsh.effort_weight_by_fishery(fi)(1,nrft),
                sqrt(.01+fsh.effort_by_fishery(fi)(1,nrft)));
    
    
            dvariable ssum=mean(fsh.effort_dev_coffs(fi)(1,nrft));
            dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
            dvariable eff_dev_pen=-sum(log(1.e-10+exp(elem_prod(-eff_wt,tmp_square))+
              .05*exp(elem_prod(-eff_wt/5,tmp_square))));
            tmppen(fi)+=eff_dev_pen;
          }
        }
        else
        {
          int nrft=fsh.num_real_fish_times(fi);
          if (fsh.fish_flags(fi,13)==0)
          {
            MY_DOUBLE_TYPE eff_wt=10.;
            dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
            dvariable eff_dev_pen=-sum(log(1.e-10+exp(-eff_wt*tmp_square)+
              .01/(1.0+eff_wt/3.*tmp_square)));
            tmppen(fi)+=eff_dev_pen;
          }
          else if (fsh.fish_flags(fi,13)>0)
          {
            MY_DOUBLE_TYPE eff_wt=fsh.fish_flags(fi,13);
            dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
            dvariable eff_dev_pen=-sum(log(1.e-10+exp(-eff_wt*tmp_square)+
              .01/(1.0+eff_wt/3.*tmp_square)));
            tmppen(fi)+=eff_dev_pen;
          }
          else if (fsh.fish_flags(fi,13)<0)
          {
            dvar_vector eff_wt=-fsh.fish_flags(fi,13)*
              sqrt(.01+fsh.effort_by_fishery(fi)(1,nrft));
            dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
            dvariable eff_dev_pen=
               -sum(log(1.e-10+exp(elem_prod(-eff_wt,tmp_square))+
              .01/(1.0+elem_prod(eff_wt/5,tmp_square))));
            tmppen(fi)+=eff_dev_pen;
          }
        }
      }
      else   // have missing effort so need to remove those from the
      {      // penalties on effort_devs  DF Oct 27 2008
        int nrft=fsh.num_real_fish_times(fi);
        ivector merf=fsh.missing_effort_by_realization_flag(fi)(1,nrft);
        if (fsh.fish_flags(fi,65)==0)
        {
          if (fsh.fish_flags(fi,13)==0)
          {
            if (fsh.fish_flags(fi,66)==0)
            {
              MY_DOUBLE_TYPE eff_wt=10.;
              dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
              dvariable eff_dev_pen=-merf*(log(1.e-10+exp(-eff_wt*tmp_square)+
                .05*exp(-eff_wt/5.*tmp_square)));
              tmppen(fi)+=eff_dev_pen;
            }
            else
            {
              dvector eff_wt=10.*fsh.effort_weight_by_fishery(fi)(1,nrft);
              dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
              dvariable eff_dev_pen=
                -merf*(log(1.e-10+exp(-elem_prod(eff_wt,tmp_square))+
                .05*exp(-elem_prod(eff_wt/5.,tmp_square))));
              tmppen(fi)+=eff_dev_pen;
            }
          }
          else if (fsh.fish_flags(fi,13)>0)
          {
            if (fsh.fish_flags(fi,66)==0)
            {
              MY_DOUBLE_TYPE eff_wt=fsh.fish_flags(fi,13);
              dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
              dvariable eff_dev_pen=-merf*(log(1.e-10+exp(-eff_wt*tmp_square)+
                .05*exp(-eff_wt/5.*tmp_square)));
              tmppen(fi)+=eff_dev_pen;
            }
            else
            {
              dvector eff_wt=fsh.fish_flags(fi,13)*
                fsh.effort_weight_by_fishery(fi)(1,nrft);
              dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
              dvariable eff_dev_pen=
                -merf*(log(1.e-10+exp(-elem_prod(eff_wt,tmp_square))+
                .05*exp(-elem_prod(eff_wt/5.,tmp_square))));
              tmppen(fi)+=eff_dev_pen;
            }
          }
          else if (fsh.fish_flags(fi,13)<0)
          {
            dvector eff_wt;
            if (fsh.fish_flags(fi,66)==0)
              eff_wt=-fsh.fish_flags(fi,13)*
                sqrt(.01+fsh.effort_by_fishery(fi)(1,nrft));
            else
              eff_wt=-fsh.fish_flags(fi,13)*
                elem_prod(fsh.effort_weight_by_fishery(fi)(1,nrft),
                sqrt(.01+fsh.effort_by_fishery(fi)(1,nrft)));
    
    
            dvariable ssum=mean(fsh.effort_dev_coffs(fi)(1,nrft));
            dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
            dvariable eff_dev_pen=-merf*(log(1.e-10+exp(elem_prod(-eff_wt,tmp_square))+
              .05*exp(elem_prod(-eff_wt/5,tmp_square))));
            tmppen(fi)+=eff_dev_pen;
          }
        }
        else
        {
          if (fsh.fish_flags(fi,13)==0)
          {
            MY_DOUBLE_TYPE eff_wt=10.;
            dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
            dvariable eff_dev_pen=-merf*(log(1.e-10+exp(-eff_wt*tmp_square)+
              .01/(1.0+eff_wt/3.*tmp_square)));
            tmppen(fi)+=eff_dev_pen;
          }
          else if (fsh.fish_flags(fi,13)>0)
          {
            MY_DOUBLE_TYPE eff_wt=fsh.fish_flags(fi,13);
            dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
            dvariable eff_dev_pen=-merf*(log(1.e-10+exp(-eff_wt*tmp_square)+
              .01/(1.0+eff_wt/3.*tmp_square)));
            tmppen(fi)+=eff_dev_pen;
          }
          else if (fsh.fish_flags(fi,13)<0)
          {
            dvar_vector eff_wt=-fsh.fish_flags(fi,13)*
              sqrt(.01+fsh.effort_by_fishery(fi)(1,nrft));
            dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
            dvariable eff_dev_pen=
               -merf*(log(1.e-10+exp(elem_prod(-eff_wt,tmp_square))+
              .01/(1.0+elem_prod(eff_wt/5,tmp_square))));
            tmppen(fi)+=eff_dev_pen;
          }
        }
      }
      if (fsh.ppstf)
      {
        fsh.ppstf->effort_dev_penalty_by_fishery(fi)=value(tmppen(fi)); 
      }
    
    } //NMD ff83

  }
  dvariable sumpen=sum(tmppen);
  cout << "effort dev penalty " << sumpen << endl;
  return sumpen;
}
static void error_message(void)
{
  cerr << "this option not implemented for implicit catchability" << endl;
  ad_exit(1);
}

dvariable implicit_eff_dev_penalty(dvar_fish_stock_history& fsh,int print_switch) 
{
  dvariable tmppen=0.;
  const MY_DOUBLE_TYPE tpi=2.0*3.14159;
  
  for (int fi=1;fi<=fsh.num_fisheries;fi++)
  { 
    // check whether this fishery is done in the times series changes 
    // in catchability 
    if (fsh.age_flags(104)==0 || fsh.fish_flags(fi,10)==0)
    {
      // ff(65) selects for fati-tailed cauchy dist if ==1
      // doesn't work with catch-conditioned option
      if (fsh.fish_flags(fi,65)==0)
      {
        int nrft=fsh.num_real_fish_times(fi);
        if (fsh.fish_flags(fi,13)==0)
        {
          // XXXXXXXXXXXXXXXXXXXXXXxx
          //error_message();
          MY_DOUBLE_TYPE eff_wt=10.;
          dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
          dvariable eff_dev_pen=-sum(log(1.e-10+exp(-eff_wt*tmp_square)+
            .05*exp(-eff_wt/5.*tmp_square)));
          tmppen+=eff_dev_pen;
        }
        else if (fsh.fish_flags(fi,13)>0)
        {
          error_message();
          MY_DOUBLE_TYPE eff_wt=fsh.fish_flags(fi,13);
          dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
          dvariable eff_dev_pen=-sum(log(1.e-10+exp(-eff_wt*tmp_square)+
            .05*exp(-eff_wt/5.*tmp_square)));
          tmppen+=eff_dev_pen;
        }
        else if (fsh.fish_flags(fi,13)<0)
        {
          dvar_vector eff_wt=-fsh.fish_flags(fi,13)*
            sqrt(.01+fsh.effort_by_fishery(fi)(1,nrft));
          //dvar_vector iq=fsh.implicit_catchability(fi);
          //iq-=mean(iq);
          //iq-=fsh.q0(fi);
          //fsh.effort_dev_coffs(fi)(1,nrft)=iq;
          
          if (fsh.fish_flags(fi,27)==1)  // seasonal catchability
          {
            int rr=fsh.realization_region(fi,1);
            dvar_vector  sce=
              fsh.get_seasonal_catchability_effect(fi);
            for (int it=1;it<=fsh.num_real_fish_times(fi);it++)
            {
              int rp=fsh.realization_period(fi,it);
              fsh.effort_dev_coffs(fi,it)=
                fsh.implicit_catchability(fi,it)-fsh.q0(fi)
                -sce(fsh.true_month(rr,rp));
            }
            
            int j=fsh.fish_flags(fi,29);  // the group the fishery fi belongs to
            MY_DOUBLE_TYPE cfp1=value(fsh.fish_pars(1,fi));
            MY_DOUBLE_TYPE cfp2=value(fsh.fish_pars(2,fi));

            dvector seasonal= 
              cfp1*sin(tpi*(fsh.grouped_true_month(j)/12.-cfp2));

            fsh.grouped_catchability_coffs(j)=value(fsh.q0(fi))+seasonal;
           
          }
          else
         
          {
            int j=fsh.fish_flags(fi,29);  // the group the fishery fi belongs to
            fsh.effort_dev_coffs(fi)(1,nrft)=
              fsh.implicit_catchability(fi)-fsh.q0(fi);
            //fsh.grouped_catchability_coffs(j)=value(fsh.q0(fi));
          } 
  
          dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
          dvariable eff_dev_pen=-sum(log(1.e-10+exp(elem_prod(-eff_wt,tmp_square))+
            .05*exp(elem_prod(-eff_wt/5,tmp_square))));
          tmppen+=eff_dev_pen;
        }
      }
      else
      {
        error_message();
        int nrft=fsh.num_real_fish_times(fi);
        if (fsh.fish_flags(fi,13)==0)
        {
          MY_DOUBLE_TYPE eff_wt=10.;
          dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
          dvariable eff_dev_pen=-sum(log(1.e-10+exp(-eff_wt*tmp_square)+
            .01/(1.0+eff_wt/3.*tmp_square)));
          tmppen+=eff_dev_pen;
        }
        else if (fsh.fish_flags(fi,13)>0)
        {
          MY_DOUBLE_TYPE eff_wt=fsh.fish_flags(fi,13);
          dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
          dvariable eff_dev_pen=-sum(log(1.e-10+exp(-eff_wt*tmp_square)+
            .01/(1.0+eff_wt/3.*tmp_square)));
          tmppen+=eff_dev_pen;
        }
        else if (fsh.fish_flags(fi,13)<0)
        {
          dvar_vector eff_wt=-fsh.fish_flags(fi,13)*
            sqrt(.01+fsh.effort_by_fishery(fi)(1,nrft));
          dvar_vector tmp_square=square(fsh.effort_dev_coffs(fi)(1,nrft));
          dvariable eff_dev_pen=
             -sum(log(1.e-10+exp(elem_prod(-eff_wt,tmp_square))+
            .01/(1.0+elem_prod(eff_wt/5,tmp_square))));
          tmppen+=eff_dev_pen;
        }
      }
    }
  }
  cout << "implict effdev pen = " << tmppen << endl;
  return tmppen;
}

dvariable selectivity_deviations_penalty(dvar_fish_stock_history& fsh,
  int print_switch)
{
  dvariable tmp2=0.;
  for (int i=1;i<=fsh.num_fisheries;i++)
  {
    if (fsh.fish_flags(i,3)>0)
    {
      for (int nt=2;nt<=fsh.num_fish_times(i);nt++)
      {
	int rp=fsh.realization_period(i,nt);
	int ri=fsh.realization_incident(i,nt);
        if (nt<=2)
        {
	  for (int j=1;j<=fsh.fish_flags(i,3);j++)
	  {
	    tmp2+=(1000.* fsh.between_times(i,nt)) 
                        * square(fsh.delta2(i,nt,j));
          }
        }
        else
        {
	  for (int j=1;j<=fsh.fish_flags(i,3);j++)
	  {
	    tmp2+=(1000.* fsh.between_times(i,nt)) 
              * square(fsh.delta2(i,nt,j)- fsh.delta2(i,nt-1,j));
          }
	}
      }
    }
  }
  return tmp2;
}
 

dvariable selectivity_form_penalty(dvar_len_fish_stock_history& fsh,
  int print_switch)
{
  dvariable pen=0.0; 
  dvariable pen1=0.0;
  ivector ff71=column(fsh.fish_flags,71);
  ivector ff74=column(fsh.fish_flags,74);
  for (int i=1;i<=fsh.num_fisheries;i++)
  {
    if (fsh.fish_flags(i,16)==1)
    {
      MY_DOUBLE_TYPE penwt=1.e+6;
      if (fsh.fish_flags(i,56)!=0) penwt=fsh.fish_flags(i,56);
      // for splines lendth dep sel put the penalty on the
      // spline values
      if (fsh.fish_flags(i,26)==3)
      {
        for (int j=2;j<=fsh.nlint;j++)
        {
          dvar_vector & spl=fsh.bstempsel(i,1,1);
          if (spl(j)<spl(j-1))
          {
            dvariable diff=spl(j-1)-spl(j);
            pen+=penwt*diff*diff*diff;
          }
          if (fsh.nwint>0)
          {
            dvar_vector & wspl=fsh.bswtempsel(i,1,1);
            if (wspl(j)<wspl(j-1))
            {
              dvariable diff=wspl(j-1)-wspl(j);
              pen+=penwt*diff*diff*diff;
            }
          }
        }
      }
      else
      {
        if ((ff71(i)>0) || (ff74(i)>1))
        {
          for (int is=1;is<=ff74(i);is++)
          {
            if (fsh.selseas_ff16(i,is))
            {
              for (int ib=1;ib<=fsh.num_blocks(i);ib++)
              {
                if (fsh.selblks_ff16(i,ib))
                {
                  dvar_vector & spl=fsh.bstempsel(i,is,ib);
                  int mmin=spl.indexmin();
                  int mmax=spl.indexmax();
                  for (int j=mmin+1;j<=mmax;j++)
                  {
                    if (spl(j)<spl(j-1))
                    {
                      dvariable diff=spl(j-1)-spl(j);
                      pen+=penwt*diff*diff*diff;
                    }
                  }
                }
              }
            }
          }
        }
        else
        {
          int rr=fsh.realization_region(i,1);
          int rp=fsh.realization_period(i,1);
          int ri=fsh.realization_incident(i,1);
          dvar_vector& insel=fsh.incident_sel(rr,rp,ri);
          // penalty to make selectivity a non decreasing function of age
          for (int j=2;j<=fsh.nage;j++)
          {
            if (insel(j)<insel(j-1))
            {
              dvariable diff=insel(j-1)-insel(j);
              pen+=penwt*diff*diff*diff;
            }
          }
        }
      }
    }
    if (fsh.fish_flags(i,16)==2)
    {
      int rr=fsh.realization_region(i,1);
      int rp=fsh.realization_period(i,1);
      int ri=fsh.realization_incident(i,1);
   //   dvar_vector& insel1=fsh.incident_sel(rr,rp,ri);
      const dvar_vector& insel1=mfexp(fsh.fish_mort(rr,rp,ri))/
             max(mfexp(fsh.fish_mort(rr,rp,ri)));
   //   insel1=exp(insel1);
   //   insel1/=max(insel1);
      // penalty to make selectivity tend towards zero for older age classes
      int jj=fsh.fish_flags(i,3);
      MY_DOUBLE_TYPE xpen=1.e+6;
      // this penalty is too large for these forms of selectivity parameterization
      if (fsh.fish_flags(i,57)>0)
      {
        xpen=100;
//      }
        for (int j=jj;j<=fsh.nage;j++)
        {
          pen1+=xpen*insel1(j);
        }
      }
    }
  }
  cout << "dome selectivity pen = " << pen1 << endl;
  cout << "non decreasing  selectivity pen = " << pen << endl;
  pen+=pen1;
  return pen;
}
dvariable selmean_penalty(dvar_fish_stock_history& fsh,
  int print_switch)
{
  dvariable tmp=0.0;
  for (int fi=1;fi<=fsh.num_fisheries;fi++)
  {
   /*  use bs_selmean now
    if (fsh.fish_flags(fi,57)==3)
    {
      tmp+=1000*square(fsh.selmean(fi));
    }
  */
    if (fsh.fish_flags(fi,63))
    {
      tmp+=1000*square(fsh.ageselmean(fi));
    }
  }  
  return tmp;
}

dvariable bs_selmean_penalty(dvar_fish_stock_history& fsh,
  int print_switch)
{
  ivector nd=column(fsh.fish_flags,3);
  ivector ff57=column(fsh.fish_flags,57);
  ivector ff74=column(fsh.fish_flags,74);
  ivector ff71=column(fsh.fish_flags,71);
  dvariable tmp=0.0;
  for (int i=1;i<=fsh.num_fisheries;i++)
  {
    for (int is=1;is<=ff74(i);is++)
    {
      for (int ib=1;ib<=ff71(i)+1;ib++)
      {
        if (ff57(i)==3)
        {
          tmp+=1000*square(fsh.bs_selmean(i,is,ib));
        }
      }
    }
  }
  cout << "selmean pen = " << tmp << endl;
  return tmp;
}

 

dvariable incident_sel_curvature_penalty(dvar_fish_stock_history& fsh,
  int print_switch)
{
  dvariable tmp1=0.;
  dvariable tmp2=0.;
  for (int fishery=1;fishery<=fsh.num_fisheries;fishery++)
  {
    int ir=fsh.realization_region(fishery,1);
    int ip=fsh.realization_period(fishery,1);
    int fi=fsh.realization_incident(fishery,1);

    if (fsh.fish_flags(fishery,3)>0)
    {
      const dvar_vector& fis=mfexp(fsh.incident_sel(ir,ip,fi));
      MY_DOUBLE_TYPE penwt1=0.01;  // !!! old penalty wt for second diff
      MY_DOUBLE_TYPE penwt2=0.01;  // !!! penalty wt for third difference
      if (fsh.fish_flags(fishery,41)>0)
        penwt1=fsh.fish_flags(fishery,41)/10.0;
      if (fsh.fish_flags(fishery,42)>0)
        penwt2=fsh.fish_flags(fishery,42)/10.0;
           
      dvar_vector vec=first_difference(
        first_difference(log(.10+fis)));

      if (fsh.fish_flags(fishery,41)>=0)
        tmp1+=penwt1*norm2(vec);

      if (fsh.fish_flags(fishery,42)>=0)
      {
        dvar_vector vec2=first_difference(vec);
        tmp2+=penwt2*norm2(vec2);
      }
    }
  }
  return tmp1+tmp2;
}


dvariable recr_initpop_sum_penalty(dvar_fish_stock_history& fsh,
  int print_switch)
{
  dvariable tmp1;    // XXXYYY
  if (fsh.age_flags(41))
  {
    tmp1=square(fsh.recmean)+square(fsh.initmean);
    if (fsh.pmsd)
    {
      tmp1+=norm2(fsh.pmsd->recmean);
    }
  }
  else
  {
    tmp1=square(fsh.recmean)+square(fsh.initmean);
    if (fsh.pmsd)
    {
      tmp1+=norm2(fsh.pmsd->recmean);
    }
  }
  return tmp1;
}

dvariable first_length_bias_penalty(dvar_len_fish_stock_history& fsh,
  int print_switch)
{
  // small penalty on first length bias parameters
  dvariable tmp1=.001*norm2(fsh.vb_bias);
  // penalty on catchability_deviations
  for (int i=1;i<=fsh.num_fisheries;i++)
  {
    int rr=fsh.fishery_regions(i);   //NMD_2Oct2018
    dvar_vector vvar=fsh.get_global_vars_region(rr);   
//    if (value(fsh.vb_bias(i)) > 1.75*sqrt(value(fsh.global_vars(1))))
    if (value(fsh.vb_bias(i)) > 1.75*sqrt(value(vvar(1))))  //
    {
#if !defined(NO_MY_DOUBLE_TYPE)
//      tmp1+=100.0*square(fsh.vb_bias(i)-1.75*pow(1.e-10+fsh.global_vars(1),.5L));
#else
//      tmp1+=100.0*square(fsh.vb_bias(i)-1.75*pow(1.e-10+fsh.global_vars(1),.5));
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
      tmp1+=100.0*square(fsh.vb_bias(i)-1.75*pow(1.e-10+vvar(1),.5L));   //NMD_2Oct2018
#else
      tmp1+=100.0*square(fsh.vb_bias(i)-1.75*pow(1.e-10+vvar(1),.5));   //NMD_2Oct2018
#endif
    }
  }
  return tmp1;
}

#undef HOME_VERSION
