/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


#include "all.hpp"
#ifdef __MSC__
  MY_DOUBLE_TYPE mfexp(const MY_DOUBLE_TYPE& );
  dvariable mfexp(const prevariable& );
  dvar_vector mfexp(const dvar_vector& );
  dvector mfexp(const dvector& );
#else
  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
  dvariable mfexp(prevariable& );
  //dvar_vector mfexp(dvar_vector& );
  //dvector mfexp(dvector& );
#endif


void dvar_len_fish_stock_history::
  do_newton_raphson_missing_effort(int ir,int ip,dvariable& ffpen)
{
  int nfi=num_fish_incidents(ir,ip);
  dvar_vector Z(1,nage);
  Z.initialize();
  int ii=0;
  int nfi1=num_missing_effort_by_region(ir,ip);
  const dvar_vector& enf=mfexp(num_fish(ir,ip));
  dvariable tot_num_fish=sum(enf);
  dvar_matrix sel(1,nfi1,1,nage);
  dvar_matrix logsel(1,nfi1,1,nage);
  dvar_matrix fm(1,nfi1);
  dvariable ortc=0;
  dvar_vector otc(1,nfi1);
  dvariable tnf=0.0;
  for (int fi=1;fi<=nfi;fi++)
  {
    if (missing_effort_by_region_flag(ir,ip,fi)==0)
    {
      // have effort data so just add fishing mortality to total mortality
      Z+=fish_mort(ir,ip,fi);
    }
    else
    {
      int i=parent(ir,ip,fi);
      int rr=realization_region(i,1);
      int rp=realization_period(i,1);
      int ri=realization_incident(i,1);
      logsel(ii)=incident_sel(rr,rp,ri);
      sel(ii)=mfexp(logsel(ii));
      fm(ii)=fish_mort(ir,ip,fi);

      if (data_fish_flags(1,parent(ir,ip,fi))) 
      {
        otc(ii)=obs_tot_catch(ir,ip,fi);
      }
      else // convert from weigth to numbers
      {
        otc(ii)=obs_tot_catch(ir,ip,fi)/get_mean_weight(ir,ip,fi);
      }
      tnf+=enf(ii);
      ii++;
    }
  }
}    
    
  //    
  //    MY_DOUBLE_TYPE sb=0.1;
  //    if (age_flags(159))
  //    {
  //      sb=age_flags(159)/10.;
  //    }
  //    MY_DOUBLE_TYPE sb1=0.1;
  //    int nmiss=missing_catch_for_period_flag(ir,ip); 
  //    MY_DOUBLE_TYPE nfi=num_fish_incidents(ir,ip);
  //    int nfi1=nfi-nmiss;
  //    int fi;
  //    int ii=1;
  //    if (nfi1==0) 
  //    {
  //      for (fi=1;fi<=nfi;fi++)
  //      {
  //        FF(ir,ip,fi)=catchability(ir,ip,fi)+effort(ir,ip,fi);
  //      }
  //      else
  //      {
  //        int i=parent(ir,ip,fi);
  //        int rr=realization_region(i,1);
  //        int rp=realization_period(i,1);
  //        int ri=realization_incident(i,1);
  //        logsel(ii)=incident_sel(rr,rp,ri);
  //        sel(ii)=mfexp(logsel(ii));
  //        fm(ii)=fish_mort(ir,ip,fi);
  //  
  //        if (data_fish_flags(1,parent(ir,ip,fi))) 
  //        {
  //          otc(ii)=obs_tot_catch(ir,ip,fi);
  //        }
  //        else // convert from weigth to numbers
  //        {
  //          otc(ii)=obs_tot_catch(ir,ip,fi)/get_mean_weight(ir,ip,fi);
  //        }
  //          tnf+=enf(ii);
  //          ii++;
  //        }
  //        ortc=sum(otc);
  //      }
  //  
  //  
  //  
  //      }
  //    }
  //    else
  //    {
  //      const dvar_vector& enf=mfexp(num_fish(ir,ip));
  //      dvariable tot_num_fish=sum(enf);
  //      dvar_matrix sel(1,nfi1,1,nage);
  //      dvar_matrix logsel(1,nfi1,1,nage);
  //      dvar_matrix fm(1,nfi1);
  //      dvariable ortc=0;
  //      dvar_vector otc(1,nfi1);
  //      int ii=1;
  //      dvariable tnf=0.0;
  //      for (fi=1;fi<=nfi;fi++)
  //      {
  //        if (missing_catch_for_incident_flag(ir,ip,fi)==0)
  //        {
  //          int i=parent(ir,ip,fi);
  //          int rr=realization_region(i,1);
  //          int rp=realization_period(i,1);
  //          int ri=realization_incident(i,1);
  //          logsel(ii)=incident_sel(rr,rp,ri);
  //          sel(ii)=mfexp(logsel(ii));
  //          fm(ii)=fish_mort(ir,ip,fi);
  //  
  //          if (data_fish_flags(1,parent(ir,ip,fi))) 
  //          {
  //            otc(ii)=obs_tot_catch(ir,ip,fi);
  //          }
  //          else // convert from weigth to numbers
  //          {
  //            otc(ii)=obs_tot_catch(ir,ip,fi)/get_mean_weight(ir,ip,fi);
  //          }
  //          tnf+=enf(ii);
  //          ii++;
  //        }
  //        ortc=sum(otc);
  //      }
  //    
  //      dvar_matrix M(1,nfi1,1,nfi1);
  //      M.initialize();
  //      dvar_matrix C(1,nfi1,1,nage);
  //      if (nfi>1)
  //      {
  //        //cout <<"test2.cpp " << nfi << endl;
  //      }
  //      dvariable surv_rate;
  //      if (tnf>1.e-8)
  //        surv_rate=(1.0-obs_region_tot_catch(ir,ip)/tnf);
  //      else
  //        surv_rate=(1.0-obs_region_tot_catch(ir,ip)/1.e-8);
  //      //cout << "surv_rate " << surv_rate << endl;
  //      dvar_vector kc(1,nfi1);
  //      if (surv_rate<=0.2)
  //      {
  //        dvariable ks=sb+posfun(surv_rate-sb,sb1,ffpen);
  //    
  //        dvariable kr= (1.0-ks)*tnf/ortc;
  //    
  //        kc=kr*otc+1.e-10;
  //    
  //        //dvariable kludge=tnf*posfun(surv_rate,0.2,ffpen)/
  //          //obs_region_tot_catch(ir,ip);
  //        //kc=kludge*obs_tot_catch(ir,ip)+1.e-10;
  //      }
  //      else
  //      {
  //        kc=otc+1.e-10;
  //      }
  //      dvar_vector qq=elem_div(kc,sel*enf);
  //      qq/=4.0;
  //      dvar_vector TC(1,nfi1);
  //      dvar_vector diff(1,nfi1);
  //      //for (int itt=1;itt<=2;itt++)
  //      int itt=0;
  //      int badflag=0;
  //      do
  //      {
  //        int yr=year(ir,ip);  
  //        dvar_vector Z=qq*sel+mfexp(nat_mort(yr)+fraction(ir,ip));
  //        dvar_vector S=exp(-Z);
  //        dvar_vector S1N=elem_prod(1.0-S,enf);
  //        for (int fi=1;fi<=nfi1;fi++)
  //        {
  //          dvar_vector t1=elem_div(qq(fi)*sel(fi),Z);
  //          C(fi)=elem_prod(t1,S1N);
  //          TC(fi)=sum(C(fi));
  //          dvariable Dq=sum(C(fi))/qq(fi);
  //          dvar_vector DZ=qq(fi)*elem_prod(elem_div(sel(fi),Z),enf)-
  //            elem_prod(C(fi),1.0+1.0/Z);
  //          M(fi,fi)=Dq;
  //          for (int fj=1;fj<=nfi1;fj++)
  //          {
  //            M(fi,fj)+=DZ*sel(fj);
  //          }
  //        }
  //        diff=TC-kc;
  //        //cout <<"test2.cpp " <<  elem_div(diff,TC)  << endl;
  //        int pflag=1;
  //        if (!badflag)
  //        {
  //          // this is newton raphson for q
  //          qq-=inv(M)*diff;
  //          pflag=check_pos(qq,nfi1);
  //          if (!pflag)
  //          {
  //    	badflag=1;
  //            qq+=inv(M)*diff;
  //            qq=elem_div(qq,mfexp(elem_div(inv(M)*diff,qq)));
  //          }
  //          else
  //          {
  //            itt++;
  //          }
  //        }
  //        else
  //        {
  //          // this i newton raphson for log(q)
  //          qq=elem_div(qq,1.e-20+mfexp(elem_div(inv(M)*diff,qq)));
  //          //cout << " " << qq(1);
  //          itt++;
  //          pflag=check_pos(qq,nfi);
  //        }
  //        if ( (itt>2 && pflag && !badflag) || (itt>2 && pflag) )break;
  //      }
  //      while(1);
  //      MY_DOUBLE_TYPE nd=norm(elem_div(value(diff),1.e-4+value(kc)));
  //      if (nd>.1)
  //        cout << "norm(diff) = " << nd << endl;
  //      //cout <<"test2.cpp " << endl;
  //      ii=1;
  //      for (fi=1;fi<=nfi;fi++)
  //      {
  //        if (missing_catch_for_incident_flag(ir,ip,fi)==0)
  //        {
  //          FF(ir,ip,fi)=log(qq(ii));
  //          ii++;
  //        }
  //        else
  //        {
  //          FF(ir,ip,fi)=catchability(ir,ip,fi)+effort(ir,ip,fi);
  //        }
  //      }
  //    }
  //    dvar_matrix totallogsel(1,nfi,1,nage);
  //    for (fi=1;fi<=nfi;fi++)
  //    {
  //      int i=parent(ir,ip,fi);
  //      int rr=realization_region(i,1);
  //      int rp=realization_period(i,1);
  //      int ri=realization_incident(i,1);
  //      totallogsel(fi)=incident_sel(rr,rp,ri);
  //    }
  //    dvar_matrix& totalfm=fish_mort(ir,ip);
  //    for (fi=1;fi<=nfi;fi++)
  //    {
  //      totalfm(fi)=totallogsel(fi)+FF(ir,ip,fi);
  //    }
  //    dvar_vector& tm=tot_mort(ir,ip);
  //    tm.initialize();
  //    for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
  //    {                                         // incidents for this period
  //      tm+=mfexp(totalfm(fi));
  //    }
  //    tm+=mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
  //    survival(ir,ip)=mfexp(-tm);
  //    //cout <<"test2.cpp " << survival(ir,ip) << endl;
  //  }
  //  */
  //  
  //  
  //  
  //  dvariable dvar_len_fish_stock_history::get_mean_weight(int ir,int ip,int fi)
  //  {
  //    dvariable mw;
  //    dvariable sv27=get_sv_region(ir,27);  //NMD 12Dec2011
  //    dvariable sv28=get_sv_region(ir,28);  //NMD 12Dec2011
  //    if (value(sv28)==3.0)  //NMD 12Dec2011
  //    {
  //      mw=len_wt_coff/1000.* 
  //         (prop(ir,ip,fi)*
  //         (pow(mean_length(ir,ip,fi),3)+
  //         3.0*elem_prod(mean_length(ir,ip,fi),vars(ir,ip,fi))));
  //    }
  //    else
  //    {
  //      mw=0.0;
  //      for (int j=1;j<=nage;j++)
  //      {  
  //        mw+=prop(ir,ip,fi,j)*
  //          normal_length_to_weight(0.5,-3.5,
  //          mean_length(ir,ip,fi,j),
  //            sqrt(vars(ir,ip,fi,j)),value(sv27),  //NMD 12Dec2011
  //          value(sv28));  //NMD 12Dec2011
  //      }
  //      mw/=1000.; 
  //    }
  //    return mw;
  //  }
  //  
  //  dvar_vector dvar_len_fish_stock_history::get_kludged_catch_with_totcatch(int ir,
  //    int ip,dvariable& ffpen,dvar_vector& enf)
  //  {
  //    int fi;
  //    dvariable surv_rate;
  //  
  //    int nfi=num_fish_incidents(ir,ip);
  //  
  //    dvar_vector otc(1,nfi);
  //    
  //    for (fi=1;fi<=nfi;fi++)
  //    {
  //      if (data_fish_flags(1,parent(ir,ip,fi))) 
  //      {
  //        otc(fi)=obs_tot_catch(ir,ip,fi)/get_mean_weight(ir,ip,fi);
  //      }
  //      else
  //      {
  //        otc(fi)=obs_tot_catch(ir,ip,fi);
  //      }
  //    }
  //    dvariable ortc=sum(otc);
  //  
  //  
  //    dvariable tnf=sum(enf);
  //    if (tnf>1.e-8)
  //      surv_rate=(1.0-ortc/tnf);
  //    else
  //      surv_rate=(1.0-ortc/1.e-8);
  //    
  //    MY_DOUBLE_TYPE sb=0.1;
  //    if (age_flags(159))
  //    {
  //      sb=age_flags(159)/10.;
  //    }
  //    MY_DOUBLE_TYPE sb1=0.1;
  //    dvar_vector kc(1,nfi);
  //    if (surv_rate<=sb)
  //    {
  //      dvariable ks=sb+posfun(surv_rate-sb,sb1,ffpen);
  //  
  //      dvariable kr= (1.0-ks)*tnf/ortc;
  //  
  //      kc=kr*otc+1.e-10;
  //    }
  //    else
  //    {
  //      kc=otc+1.e-10;
  //    }
  //    return kc;
  //  }
