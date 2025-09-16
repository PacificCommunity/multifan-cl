/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"

  
void block_error();
extern dvar_vector * const_vb;
extern dvar_vector * const_var;

void dvar_len_fish_stock_history::incident_selectivity_calc
  (d3_array& len_sample_size,dvar_vector& vb_coff,dvar_vector& var_coff)
{
  ivector nd=column(fish_flags,3);
  ivector ff26=column(fish_flags,26);
  ivector ff57=column(fish_flags,57);
  int i,j;
  for (i=1;i<=num_fisheries;i++)
  {
    if (nd(i)==0) nd(i)=nage;
  }
  ivector tmplength(1,num_fisheries);
  for (i=1;i<=num_fisheries;i++)
  {
    if (ff26(i)==3)
      tmplength(i)=nlint;
    else
      tmplength(i)=nd(i);
  }
  
  dvar_matrix tempsel(1,num_fisheries,1,tmplength);
  tempsel.initialize();
  dvariable rho=0.0;
  dvar_matrix tlength(1,num_fisheries,1,nd);
  for (i=1;i<=num_fisheries;i++)
  {
    int ir=fishery_regions(i);
    if (!pmsd  || pmsd->region_species_pointer(ir)==1)
      rho=exp(-vb_coff(3));
    else
      rho=exp(-pmsd->vb_coff(pmsd->region_species_pointer(ir),3));

    tlength(i,1)=0.0;
    dvariable div=1.0/(1.0-pow(rho,nd(i)-1));
    for (j=2;j<=nd(i);j++)
    {
      tlength(i,j)=(1.0-pow(rho,j-1))*div;
    }
  }
  //dvar_matrix lbsel(1,num_fisheries,1,nage);
  if (age_flags(185)==0)
  {
    lbsel_calc(lbsel,*this,vb_coff,var_coff);
  }
  else
  {
    lbsel_calc(lbsel,*this,*const_vb,*const_var);
  }
  //cout << "lbsel(10)" << endl;
  //cout <<"rshort1.cpp " << lbsel(10) << endl;
  //double min_select=.025;
  //double min_select=.0001;
  MY_DOUBLE_TYPE min_select=.000001;
  MY_DOUBLE_TYPE rnage=sqrt(double(nage));
  MY_DOUBLE_TYPE one_plus=1.+min_select;
  MY_DOUBLE_TYPE alpha=min_select/exp(1.);
  dvar_vector& sel_corr = corr_wy;
  dvar_vector root_sel_corr=sqr(1.0-square(sel_corr));
  int num_age=0;
  dvar3_array btempsel;
  ivector ff71=column(fish_flags,71);
  if (sum_ff71_flag)
  {
    int nb=sel_block_index.indexmax();
    ivector  xtmpl(1,nb);
    for (int i=1;i<=nb;i++)
    {
      xtmpl(i)=tmplength(sel_block_fisheries(i));
    }
    btempsel.allocate(1,nb,1,sel_block_index,1,xtmpl);
  }
  for (i=1;i<=num_fisheries;i++)
  {
    for (int ic=0;ic<=ff71(i);ic++)
    {
      if (fish_flags(i,57) !=3 )
      {
        if (ic>0)
        {
          if (fish_flags(i,57) !=1 )
          block_error();
        }
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
          if (ff26(i)==0)
          {
            for (int j=1;j<=num_age;j++)
            {
              incident_sel(rr,rp,ri,j) = exp(selcoff(i,min(fish_flags(i,3),j)));
              if (incident_sel(rr,rp,ri,j)< min_select)
              {
                incident_sel(rr,rp,ri,j)=alpha*exp(incident_sel(rr,rp,ri,j)/
                  min_select);
              }
            }
            //cout <<"rshort1.cpp " << incident_sel(rr,rp,ri) << endl;
          }
          else
          {
            if (sum_ff71_flag==0 || fish_flags(i,71)==0)
            {
              for (int j=1;j<=num_age;j++)
              {
                incident_sel(rr,rp,ri,j) = lbsel(i,j);
                if (incident_sel(rr,rp,ri,j)< min_select)
                {
                  incident_sel(rr,rp,ri,j)=alpha*exp(incident_sel(rr,rp,ri,j)/
                    min_select);
                }
              }
            }
            else
            {
              for (int j=1;j<=num_age;j++)
              {
                incident_sel(rr,rp,ri,j) = lbsel(i,j);
                if (incident_sel(rr,rp,ri,j)< min_select)
                {
                  incident_sel(rr,rp,ri,j)=alpha*exp(incident_sel(rr,rp,ri,j)/
                    min_select);
                }
              }
            }
          }
          if (!parest_flags(162))
          {
            if (fish_flags(i,30)==0)
            {
              mean_incident_sel(rr,rp,ri)=mean(incident_sel(rr,rp,ri));
              incident_sel(rr,rp,ri)/=mean_incident_sel(rr,rp,ri);
            }
            else
            {
              incident_sel(rr,rp,ri)/=incident_sel(rr,rp,ri,nage);
            }
          }
          else
          {
#if !defined(NO_MY_DOUBLE_TYPE)
            dvariable tmp=pow(sum(pow(incident_sel(rr,rp,ri),20.L)),.05L);
#else
            dvariable tmp=pow(sum(pow(incident_sel(rr,rp,ri),20.)),.05);
#endif
            incident_sel(rr,rp,ri)/=tmp;
          }
          for (int j=1;j<=nage;j++)
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
            if (sum_ff71_flag==0 || fish_flags(i,71)==0)
            {
              incident_sel(rr,rp,ri) = 
                incident_sel(rr1,realization_period(i,nt-1),
                realization_incident(i,nt-1));
            }
            else
            {
              int bfi=block_ptr(i);
              int bl=break_block(bfi,nt);
              if (bl==0)
              {
                incident_sel(rr,rp,ri) = 
                  incident_sel(rr1,realization_period(i,nt-1),
                  realization_incident(i,nt-1));
              }
              else
              {
                incident_sel(rr,rp,ri) = blbsel(i,bl);
              }
            }
           
            if (len_sample_size(rr,rp,ri)>0)
            {
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
                if (fish_flags(i,30)==0)
                {
                  incident_sel(rr,rp,ri)/=mean(incident_sel(rr,rp,ri));
                }
                else
                {
                  incident_sel(rr,rp,ri)/=incident_sel(rr,rp,ri,nage);
                }
              }
              else
              {
#if !defined(NO_MY_DOUBLE_TYPE)
                dvariable tmp=pow(sum(pow(incident_sel(rr,rp,ri),20.L)),.05L);
#else
                dvariable tmp=pow(sum(pow(incident_sel(rr,rp,ri),20.)),.05);
#endif
                incident_sel(rr,rp,ri)/=tmp;
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
      }
      else
      {
        int sd=fish_flags(i,61);
        if (sd>=0)
        {
          dvector x(1,sd);
          dvector xx(1,nd(i));
          MY_DOUBLE_TYPE longlen=fmid(nlint)+0.5*filen;
    
          if (ff26(i)==3)
          {
            x.fill_seqadd(shlen,(longlen-shlen)/(sd-1));
          }
          else
          {
            x.fill_seqadd(0,1.0/(sd-1));
            xx.fill_seqadd(0,1.0/(nd(i)-1));
          }
          dvar_vector tsel(1,sd);
          tsel=selcoff(i)(1,sd);
          if (ic>0)
          {
            tsel+=sel_dev_coffs(i,ic)(1,sd);
          }
          selmean(i)=log(mean(exp(tsel)));
          tsel-=selmean(i);
          vcubic_spline_function csf(x,tsel);
          switch (ff26(i))
          {
          case 0:
           if (ic>0)
            {
              block_error();
            }
            tempsel(i)=csf(xx);     
            break;
          case 1:
          case 2:
            if (ic==0)
            {
              tempsel(i)=csf(tlength(i));
            }
            else
            {
              int ib=block_ptr(i);
              btempsel(ib,ic)=csf(tlength(i));
            }
            break;
          case 3:
            if (ic>0)
            {
              block_error();
            }
            tempsel(i)=csf(fmid);
            splinesel(i)=mfexp(tempsel(i));
      
            if (nwint>0)
            {
              wsplinesel(i)=mfexp(csf(wmid));
            }
            break;
          default:
            cerr << "illegal value for fish_flags(" << i <<",61)"
                 << endl;
            ad_exit(1);
          }
        }
        else
        {
          switch(sd)
          {
          case -1:
            {
              const prevariable & fp9=fish_pars(9,i);
              dvariable fp10=1.0+fish_pars(10,i);
              dvariable nineteen=19.0;
              {
                MY_DOUBLE_TYPE jmid=0.5*(1+nlint);
                MY_DOUBLE_TYPE jmid1=1.0-jmid;
                
                splinesel(i)= 1.0/(1+pow(nineteen,((fmid-jmid)/jmid1-fp9)/fp10));
              }
              {
                MY_DOUBLE_TYPE jmid=0.5*(1+nwint);
                MY_DOUBLE_TYPE jmid1=1.0-jmid;
                wsplinesel(i)= 1.0/(1+pow(nineteen,((wmid-jmid)/jmid1-fp9)/fp10));
              }
            }
            break;
          case -2:
            {
              const prevariable & fp9=fish_pars(9,i);
              dvariable fp10=1.0+fish_pars(10,i);
              const dvariable & efp11=exp(fish_pars(11,i));
              MY_DOUBLE_TYPE jmid=0.5*(1+nlint);
              MY_DOUBLE_TYPE jmid1=1.0-jmid;
              int j;
              if (!allocated(splinesel(i)))
                splinesel(i).allocate(1,nlint);
                
              for (j=1;j<=nlint;j++)
              {
                MY_DOUBLE_TYPE fj=(j-jmid)/jmid1;
                if (fj<=value(fp9))
                {
                  splinesel(i,j)= pow(2,-square(fj/(fp10*efp11)));
                }
                else
                {
                  splinesel(i,j)= pow(2,-square(fj*efp11/fp10));
                }
              }
              MY_DOUBLE_TYPE wjmid=0.5*(1+nwint);
              MY_DOUBLE_TYPE wjmid1=1.0-jmid;
              if (!allocated(wsplinesel(i)))
                wsplinesel(i).allocate(1,nwint);
              for (j=1;j<=nwint;j++)
              {
                MY_DOUBLE_TYPE fj=(j-wjmid)/wjmid1;
                if (fj<=value(fp9))
                {
                  wsplinesel(i,j)= pow(2,-square(fj/(fp10*efp11)));
                }
                else
                {
                  wsplinesel(i,j)= pow(2,-square(fj*efp11/fp10));
                }
              }
            }    
            break;
          default:
            cerr << "illegal value for fish_flags(" << i <<",61)"
                 << endl;
            ad_exit(1);
          }
        }
      } 
    }
  }

  if (ff263flag)
  {
    all_length_dist_calcs_long();
    if (nwint>0)
    {
      all_weight_dist_calcs_long();
    }
  }


  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)  // Loop over fishing
      {
        int pi=parent(ir,ip,fi);
        if (fish_flags(pi,57) !=3)
        {
          for (int j=1;j<=nage;j++)
          {
            if (value(incident_sel(ir,ip,fi,j))<=0.0)
            {
              cout << ir << " " << ip << " " << fi << " " 
                   << pi << endl;
              cout <<"rshort1.cpp " << incident_sel(ir,ip,fi) << endl;
            }
          }
          incident_sel(ir,ip,fi)=log(incident_sel(ir,ip,fi));
        }
        else // for now only stationary selectivities for this option
        {
            if (fish_flags(pi,26) !=3)
            {
              if (ff71(pi)==0)
              {
                incident_sel(ir,ip,fi)(1,nd(pi))=tempsel(pi);
                if (nd(pi)<nage)
                  incident_sel(ir,ip,fi)(nd(pi),nage)=tempsel(pi,nd(pi));
              }
              else
              {
                int is=fishery_realization_index(ir,ip,fi);
                int bfi=block_ptr(pi);
                int bl=break_block(bfi,is);
                if (bl==0)
                {
                  incident_sel(ir,ip,fi)(1,nd(pi))=tempsel(pi);
                  if (nd(pi)<nage)
                    incident_sel(ir,ip,fi)(nd(pi),nage)=tempsel(pi,nd(pi));
                }
                else
                {
                  incident_sel(ir,ip,fi)(1,nd(pi))=btempsel(bfi,bl);
                  if (nd(pi)<nage)
                    incident_sel(ir,ip,fi)(nd(pi),nage)=
                    btempsel(bfi,bl,nd(pi));
                }
              }
            }
            else
            {
              int mn=month(ir,ip);
              int wk=week(ir,ip);
              int pi=parent(ir,ip,fi);
              incident_sel(ir,ip,fi)=log(1.e-20+incident_sel(ir,ip,fi));
            }
        }
      }
    }
  }
  ivector ff63=column(fish_flags,63);

  int ff63flag=sum(ff63);
  // include age-based part of selectivity
  if (ff63flag)
  {
    dvar_matrix agetempsel(1,num_fisheries,1,nd);
    for (i=1;i<=num_fisheries;i++)
    {
      // so ff63 determine the number of nodes of the cubic splines
      // for age effect in length-based selectivity. If 0 there is no effect
      int sd=ff63(i);
      if (sd>0)
      {
        dvector x(1,sd);
        dvector xx(1,nd(i));
        x.fill_seqadd(0,1.0/(sd-1));
        xx.fill_seqadd(0,1.0/(nd(i)-1));
        dvar_vector tagesel=ageselcoff(i)(1,sd);
        ageselmean(i)=log(mean(exp(tagesel)));
        tagesel-=ageselmean(i);
        vcubic_spline_function csf(x,tagesel);
        agetempsel(i)=csf(xx);     
      } 
    }
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
      {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)  // Loop over fishing
        {
          int pi=parent(ir,ip,fi);
          if (ff63(pi))
          {
            incident_sel(ir,ip,fi)(1,nd(pi))+=agetempsel(pi);
            if (nd(pi)<nage)
              incident_sel(ir,ip,fi)(nd(pi),nage)+=agetempsel(pi,nd(pi));
          }
        }
      }
    }
  }
}
