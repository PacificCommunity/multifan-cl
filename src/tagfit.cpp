/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
//#include <f:\linad99\fvar.hpp>
#include "all.hpp"
#include "df12fun.h"
//extern "C" {
//#include <pvm3.h>
//}
#include <admodel.h>
//#include <adpvm.h>
#if !defined(linux)
#include <windows.h>
#endif

#define  __declspec(dllexport) 

#ifdef __MSC__
  MY_DOUBLE_TYPE mfexp(const MY_DOUBLE_TYPE& );
  dvariable mfexp(const prevariable& );
  dvar_vector mfexp(const dvar_vector& );
  dvector mfexp(const dvector& );
#else
  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
  dvariable mfexp(BOR_CONST prevariable& );
  //dvar_vector mfexp(dvar_vector& );
  //dvector mfexp(dvector& );
#endif

extern mf_pvm_manager * mf_pvm;

ofstream ofstags("tagcont");
//ofstream ofstags1("tagcont1");
dvariable dvar_fish_stock_history::tag_rep_rate
  (int it,int ir,int rp,int fi)
{
   int pp=parent(ir,rp,fi);
   return tag_fish_rep(it,pp);
}

dvar_vector censored_lognormal(const prevariable & tau,dvar_vector& o,
  const dvar_vector& mu,MY_DOUBLE_TYPE e);

void  set_value(dvar5_array & v,MY_DOUBLE_TYPE x)
{
  int mmin1=v.indexmin();
  int mmax1=v.indexmax();
  for (int i1=mmin1;i1<=mmax1;i1++)
  {
    int mmin2=v(i1).indexmin();
    int mmax2=v(i1).indexmax();
    for (int i2=mmin2;i2<=mmax2;i2++)
    {
      int mmin3=v(i1,i2).indexmin();
      int mmax3=v(i1,i2).indexmax();
      for (int i3=mmin3;i3<=mmax3;i3++)
      {
        v(i1,i2,i3)=x;
      }
    }
  }
}
    
dvariable dvar_len_fish_stock_history::fit_tag_returns(void)
{
  dvariable f=0.0;
  MY_DOUBLE_TYPE empvar=0.0; //NMD_12May2020
  int empvar_nelem=0;
  dvariable gp_pen=0.0;
  int fi;
  dvar_matrix obsgroupedcatch;
  dvar_matrix groupedcatch;

  ivector group_flags32=column(fish_flags,32);
  int gsum32=sum(group_flags32);
  int gmax32=Max(group_flags32);
  ivector gp_fish32(1,gmax32);
  ivector iflag(1,gmax32);
  //if (allocated(grouped_fishery_projection_flag))
   // grouped_fishery_projection_flag.deallocate();
  //grouped_fishery_projection_flag.allocate(1,num_regions,
   // 1,num_fish_periods,1,ngroups);
  if (gmax32)
  {
    gp_fish32.initialize();
    for (fi=1;fi<=num_fisheries;fi++)
      gp_fish32(group_flags32(fi))=fi;
    if (allocated(obsgroupedcatch))
      obsgroupedcatch.deallocate();
    obsgroupedcatch.allocate(1,gmax32,1,nage);
    if (allocated(groupedcatch))
      groupedcatch.deallocate();
    groupedcatch.allocate(1,gmax32,1,nage);
  }
  else
  {
    if (allocated(obsgroupedcatch))
      obsgroupedcatch.deallocate();
    obsgroupedcatch.allocate(1,num_fisheries,1,nage);
    if (allocated(groupedcatch))
      groupedcatch.deallocate();
    groupedcatch.allocate(1,num_fisheries,1,nage);
  }
  dvariable totobstags=sum(obstagcatch);
  dvariable totpredtags=0.0;

  int mmin,mmax;
  if (mf_pvm->pvm_switch)
  {
    mmin=min_tag_group;
    mmax=max_tag_group;
  }
  else
  {
    mmin=1;
    mmax=num_tag_releases;
  }
  for (int it=mmin;it<=mmax;it++)
  {
    int flatflag=check_flat_tag_time_periods(it);
    if (flatflag)
    {
      cout << "non flatfor tag release " << it << endl;
    }
  }
  if (!gsum32)   // no fishery grouping
  {
    if (ppstf && allocated(ppstf->tag_like))
      ppstf->tag_like.deallocate();

    const int max_index=1000;

    d3_array tmp_tag_like;
    if (ppstf)
    {
      tmp_tag_like.allocate(mmin,mmax,1,num_fisheries,1,max_index);
      tmp_tag_like.initialize();
    }
    imatrix uicount(mmin,mmax,1,num_fisheries);
    uicount.initialize();
    for (int it=mmin;it<=mmax;it++)
    {
      int ng=nage_by_tag_release(it);
      for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
      {
        int ub;
        if (!age_flags(96))
          ub=num_fish_periods(ir);
        else
          ub=terminal_tag_period(it,ir);
        // DF  july 14 05
        if (ub>num_real_fish_periods(ir))
        {
          ub=num_real_fish_periods(ir);
        }
        for (int ip=initial_tag_period(it,ir);ip<=ub;ip++)
        {
          for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
          {
            if (fishery_projection_flag(ir,ip,fi)==0)
            {
              dvar_vector& tc=tagcatch(it,ir,ip,fi);
              int jmin=tc.indexmin();
              for (int j=jmin;j<=ng;j++)
              {
                if (tagcatch(it,ir,ip,fi,j)<0)
                {
                  cerr << tagcatch(it,ir,ip,fi,j) << endl;
                }
                if (tagcatch(it,ir,ip,fi,j)>1000)
                {
                  cerr << tagcatch(it,ir,ip,fi,j) << endl;
                }
              }
              int pp1=parent(ir,ip,fi);
              dvar_vector rtc;
              if (age_flags(198))
                rtc=tag_rep_rate(it,ir,ip,fi)*tagcatch(it,ir,ip,fi);
              else
                rtc=rep_rate(ir,ip,fi)*tagcatch(it,ir,ip,fi);
          
              totpredtags+=sum(rtc);
              int pf111=parest_flags(111);
              switch (pf111)
              {
                case 0:
                  f+=sum(elem_div(square(obstagcatch(it,ir,ip,fi)-rtc),
                     .01+rtc));
                  break;
                case 1:
                  f-=sum(log(exp(elem_div(
                    square(obstagcatch(it,ir,ip,fi)-rtc),
                    -.01-rtc ))+0.01));
                  break;
                case 2:
                {
                  dvar_vector lrtc=log(1.e-11+rtc);
                    f+=sum(rtc) - obstagcatch(it,ir,ip,fi)*lrtc;
               // for (int i=1;i<=nage;i++)
                for (int j=jmin;j<=ng;j++)
#if !defined(NO_MY_DOUBLE_TYPE)
                      f+=gammln(obstagcatch(it,ir,ip,fi,j)+1.0L);
#else
                      f+=gammln(obstagcatch(it,ir,ip,fi,j)+1.0);
#endif
                }
                  break;
                case 3:
                  {
                    dvariable a=fish_pars(4,pp1)+50.0001;
                    rtc+=1.e-11;
                    int ns=rtc.indexmax()-rtc.indexmin()+1; 
                    dvar_vector ap=log(a+rtc);
                    f+=sum(a*ap)+obstagcatch(it,ir,ip,fi)*ap;
                    f-=ns*a*log(a)+obstagcatch(it,ir,ip,fi)*log(rtc);
                    f-=sum(gammln(a+obstagcatch(it,ir,ip,fi)));
#if !defined(NO_MY_DOUBLE_TYPE)
                    f+=sum(gammln(obstagcatch(it,ir,ip,fi)+1.0L));
#else
                    f+=sum(gammln(obstagcatch(it,ir,ip,fi)+1.0));
#endif
                    f+=ns*gammln(a);
                  }
                  break;
                case 4:
                  {
                    rtc+=1.e-11;
                    dvar_vector a;
                    if (parest_flags(305)==0)
                    {
                      a=(fish_pars(4,pp1)+50.0001)*rtc;
                    }
                    else
                    {
                      dvariable tau=1.0+exp(fish_pars(4,pp1));
#if !defined(NO_MY_DOUBLE_TYPE)
                      a=rtc/(1.e-20+(tau-1.0L));
#else
                      a=rtc/(1.e-20+(tau-1.0));
#endif
                    }
                    // use wiki terms
                    /*
                    {
                       dvariable tau=1.0+exp(fish_pars(4,pp1));
#if !defined(NO_MY_DOUBLE_TYPE)
                       dvar_vector r = rtc/(tau-1.0L);
#else
                       dvar_vector r = rtc/(tau-1.0);
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
                       dvariable p=(tau-1.0L)/tau;
#else
                       dvariable p=(tau-1.0)/tau;
#endif
                       dvar_vector& O=obstagcatch(it,ir,ip,fi);
                       dvariable tmp=sum(gammln(O+r)-gammln(O+1)-gammln(r)
#if !defined(NO_MY_DOUBLE_TYPE)
                        - (r+O)*log(tau) + O*log(tau-1.0L));
#else
                        - (r+O)*log(tau) + O*log(tau-1.0));
#endif
                    }
                    */
                    int ns=rtc.indexmax()-rtc.indexmin()+1; 
                    dvar_vector ap=log(a+rtc);
                    dvariable tmpf;
                    tmpf=a*ap+obstagcatch(it,ir,ip,fi)*ap;
                    tmpf-=a*log(a)+obstagcatch(it,ir,ip,fi)*log(rtc);
                    tmpf-=sum(gammln(a+obstagcatch(it,ir,ip,fi)));
#if !defined(NO_MY_DOUBLE_TYPE)
                    tmpf+=sum(gammln(obstagcatch(it,ir,ip,fi)+1.0L));
#else
                    tmpf+=sum(gammln(obstagcatch(it,ir,ip,fi)+1.0));
#endif
                    tmpf+=sum(gammln(a));

                    uicount(it,pp1)++;
                    if (uicount(it,pp1)>max_index)
                    {
                      cerr << "need to increase max_index" << endl;
                      ad_exit(1);
                    }
                    tmp_tag_like(it,pp1,uicount(it,pp1))=value(tmpf);

                    f+=tmpf;
                    int nelem=rtc.indexmax()-rtc.indexmin()+1;
                    empvar_nelem+=nelem;
                    dvector crtc=value(rtc);
                    dvector cobstagcatch=value(obstagcatch(it,ir,ip,fi));
                    empvar+=sum(elem_div(square(cobstagcatch-crtc),crtc));   //NMD_12May2020
                  }
    
                  break;
                case 5:
                case 6:
                case 7:
                case 8:
                  {
                    rtc+=1.e-11;
                    dvar_vector obs1=obstagcatch(it,ir,ip,fi);
                    dvariable tau=exp(fish_pars(30,fi));
                    MY_DOUBLE_TYPE censor=1.1;
                    if (parest_flags(325))
                    {
                      censor=parest_flags(325)/100.0;
                    }
                    dvar_vector vtmp=censored_gamma_or_lognormal(tau,
                      obs1,rtc,censor,parest_flags(111));

                    uicount(it,pp1)++;
                    if (uicount(it,pp1)>max_index)
                    {
                      cerr << "need to increase max_index" << endl;
                      ad_exit(1);
                    }
                    dvariable tmpf=-sum(log(1.e-10+vtmp));
                    tmp_tag_like(it,pp1,uicount(it,pp1))=value(tmpf);

                    //f_by_tag(it,icount(it))-=sum(log(1.e-10+vtmp));
                    f+=tmpf;
                    int nelem=rtc.indexmax()-rtc.indexmin()+1;
                    empvar_nelem+=nelem;
                    dvector crtc=value(rtc);
                    dvector cobstagcatch=value(obstagcatch(it,ir,ip,fi));
                    empvar+=sum(elem_div(square(cobstagcatch-crtc),crtc));   //NMD_12May2020
                 }
                  break;
                default:
                  cerr << "Illegal value for parest_flags(111) value = "
                       <<  parest_flags(111)  << endl;
              
              }
            }
          }
        }
      }
    }
    if (ppstf)
    {
      if (allocated(ppstf->ungrouped_tag_like))
        ppstf->ungrouped_tag_like.deallocate();

      ppstf->ungrouped_tag_like.allocate(mmin,mmax,1,num_fisheries);
      for (int it=mmin;it<=mmax;it++)
      {
        for (int fi=1;fi<=num_fisheries;fi++)
        {
          ppstf->ungrouped_tag_like(it,fi)=
            tmp_tag_like(it,fi)(1,uicount(it,fi));
        }
      }
    }
  }
  else   // have grouping
  {
    //ofstream ofsz1("var1");
    //ofstream ofsz2("n1");
    //ofstream ofsz3("tau");
    dvar_matrix f_by_tag(mmin,mmax,1,10000);
    ivector icount(mmin,mmax);
    imatrix gicount(mmin,mmax,1,gmax32);
    icount.initialize();
    gicount.initialize();
    f_by_tag.initialize();
    int rmin=1;
    int rmax=num_regions;
    if (ppstf && allocated(ppstf->grouped_tag_like))
      ppstf->grouped_tag_like.deallocate();

    const int max_index=1000;

    d3_array tmp_grouped_tag_like;
    if (ppstf)
    {
      tmp_grouped_tag_like.allocate(mmin,mmax,1,gmax32,1,max_index);
      tmp_grouped_tag_like.initialize();
    }
    for (int it=mmin;it<=mmax;it++)
    {
      if (pmsd)
      {
        int cs=pmsd->tag_species_pointer(it);
        rmin=pmsd->region_bounds(cs,1);
        rmax=pmsd->region_bounds(cs,2);
      }

      for (int ir=rmin;ir<=rmax;ir++)
      {
        int ub;
        if (!age_flags(96))
          ub=num_real_fish_periods(ir);
        else
          ub=terminal_tag_period(it,ir);
        // DF  july 14 05
        if (ub>num_real_fish_periods(ir))
        {
          ub=num_real_fish_periods(ir);
        }
        for (int ip=initial_tag_period(it,ir);ip<=ub;ip++)
        {
          if (num_fish_incidents(ir,ip)>0)
          { 
            // deal with the number of age classes for this tag group
            // at this time
            if (allocated(groupedcatch))
              groupedcatch.deallocate();
            if (allocated(obsgroupedcatch))
              obsgroupedcatch.deallocate();
            int jmin1=tagcatch(it,ir,ip,1).indexmin();
            obsgroupedcatch.allocate(1,gmax32,jmin1,nage);
            groupedcatch.allocate(1,gmax32,jmin1,nage);

            groupedcatch.initialize();
            obsgroupedcatch.initialize();
            groupedcatch.initialize();
            iflag.initialize();
            ivector& pi=parent(ir,ip);
            dvector xtmp(1,9);
            for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
            {
              int pp1=pi(fi);
              int pp=group_flags32(pp1);

              dvar_vector rtc;
              if (age_flags(198))
              {
                rtc=tag_rep_rate(it,ir,ip,fi)*tagcatch(it,ir,ip,fi);
              }
              else
              {
                rtc=rep_rate(ir,ip,fi)*tagcatch(it,ir,ip,fi);
              }

              int jmin=rtc.indexmin();
              totpredtags+=sum(rtc);
              iflag(pp)=1;
              //obsgroupedcatch(pp)(jmin,nage)+=obstagcatch(it,ir,ip,fi); 
              obsgroupedcatch(pp)+=obstagcatch(it,ir,ip,fi); 
              //groupedcatch(pp)(jmin,nage)+=rtc;
              groupedcatch(pp)+=rtc;
            }
    
            int ig;
            for (ig=1;ig<=gmax32;ig++)
            {
              if(iflag(ig))
              {
                if (grouped_fishery_projection_flag(ir,ip,ig)==0)
                {
                  int pf111=parest_flags(111);
                  switch (pf111)
                  {
                    case 0:
                      f+=sum(elem_div(
                        square(obsgroupedcatch(ig)-groupedcatch(ig)),
                       .01+groupedcatch(ig) ));
                      break;
                    case 1:
                      f-=sum(log(exp(elem_div(
                        square(obsgroupedcatch(ig)-groupedcatch(ig)),
                        -.01-groupedcatch(ig) ))+0.01));
                      break;
                    case 2:
                      {
                        dvar_vector lrtc=log(1.e-11+groupedcatch(ig));
                          f+=sum(groupedcatch(ig)) - obsgroupedcatch(ig)*lrtc;
                      for (int i=1;i<=nage;i++)
#if !defined(NO_MY_DOUBLE_TYPE)
                            f+=gammln(obsgroupedcatch(ig,i)+1.0L);
#else
                            f+=gammln(obsgroupedcatch(ig,i)+1.0);
#endif
                    }
                      break;
                    case 3:
                      {
                        dvar_vector rtc=1.e-11+groupedcatch(ig);
                        dvariable a=fish_pars(4,gp_fish32(ig))+50.0001;
                        int ns=rtc.indexmax()-rtc.indexmin()+1; 
                        dvar_vector ap=log(a+rtc);
                        f+=sum(a*ap)+obsgroupedcatch(ig)*ap;
                        f-=ns*a*log(a)+obsgroupedcatch(ig)*log(rtc);
                        f-=sum(gammln(a+obsgroupedcatch(ig)));
#if !defined(NO_MY_DOUBLE_TYPE)
                        f+=sum(gammln(obsgroupedcatch(ig)+1.0L));
#else
                        f+=sum(gammln(obsgroupedcatch(ig)+1.0));
#endif
                        f+=ns*gammln(a);
                      }
      
                      break;
                    case 4:
                      {
                        dvar_vector obsgroupedcatch1=obsgroupedcatch(ig);
                        dvar_vector rtc;
                        if (parest_flags(372)>0)
                        {
                          rtc=1.e-11+parest_flags(372)/1000.0 + groupedcatch(ig);
                        }
                        else
                        {
                          rtc=1.e-11+groupedcatch(ig);
                        }
                        dvar_vector a;
                        if (parest_flags(305)==0)
                           a=(fish_pars(4,gp_fish32(ig)) +50.0001)*rtc;
                        else
                        {
                          dvariable tau=1.0+exp(fish_pars(4,gp_fish32(ig)));
#if !defined(NO_MY_DOUBLE_TYPE)
                          a=rtc/(1.e-20+(tau-1.0L));
#else
                          a=rtc/(1.e-20+(tau-1.0));
#endif
                        }
                        int ns=rtc.indexmax()-rtc.indexmin()+1; 
                        dvar_vector ap=log(a+rtc);
                      
                        icount(it)++;
                        gicount(it,ig)++;
                        if (gicount(it,ig)>max_index)
                        {
                          cerr << "need to increase max_index" << endl;
                          ad_exit(1);
                        }
                        // use wiki terms
                        // this is LL not neg LL
                        /*
                        {
                           dvariable tau=1.0+exp(fish_pars(4,gp_fish32(ig)));
#if !defined(NO_MY_DOUBLE_TYPE)
                           dvar_vector r = rtc/(tau-1.0L);
#else
                           dvar_vector r = rtc/(tau-1.0);
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
                           dvariable p=(tau-1.0L)/tau;
#else
                           dvariable p=(tau-1.0)/tau;
#endif
                           dvar_vector O=obsgroupedcatch(ig);
                           dvariable tmp=sum(gammln(O+r)-gammln(O+1)-gammln(r)
#if !defined(NO_MY_DOUBLE_TYPE)
                            - (r+O)*log(tau) + O*log(tau-1.0L));
#else
                            - (r+O)*log(tau) + O*log(tau-1.0));
#endif
                        }
                        */

                        f_by_tag(it,icount(it))+=a*ap+obsgroupedcatch1*ap;
                        f_by_tag(it,icount(it))-=a*log(a)+obsgroupedcatch1*log(rtc);

                        f_by_tag(it,icount(it))-=sum(gammln(a+obsgroupedcatch1));
#if !defined(NO_MY_DOUBLE_TYPE)
                        f_by_tag(it,icount(it))+=sum(gammln(obsgroupedcatch1+1.0L));
#else
                        f_by_tag(it,icount(it))+=sum(gammln(obsgroupedcatch1+1.0));
#endif
                        f_by_tag(it,icount(it))+=sum(gammln(a));
                        int nelem=rtc.indexmax()-rtc.indexmin()+1;
                        empvar_nelem+=nelem;
                        dvector crtc=value(rtc);
                        dvector cobstagcatch=value(obsgroupedcatch1);
                        /*
                        dvector tmp=elem_div(square(cobstagcatch-crtc),crtc);   //NMD_12May2020
                        ofsz1<< tmp << endl;
                        ofsz2<< crtc << endl;
                        dvariable tau=1.0+exp(fish_pars(4,gp_fish32(ig)));
                        for (int i=1;i<=nelem;i++)
                        {
                          ofsz3 << tau << " " << endl;
                        }
                        ofsz3 << endl;
                        */

                        empvar+=sum(elem_div(square(cobstagcatch-crtc),crtc));   //NMD_12May2020
                        if (ppstf)
                        {
                          tmp_grouped_tag_like(it,ig,gicount(it,ig))=
                            value(f_by_tag(it,icount(it)));
                        }
                      }
                      break;
                    case 5:
                    case 6:
                    case 7:
                    case 8:
                      {
                        dvar_vector obsgroupedcatch1=obsgroupedcatch(ig);
                        dvar_vector rtc=1.e-11+groupedcatch(ig);
                        dvariable tau=exp(fish_pars(30,gp_fish32(ig)));
                        MY_DOUBLE_TYPE censor=1.1;
                        if (parest_flags(325))
                        {
                          censor=parest_flags(325)/100.0;
                        }
//                        dvar_vector vtmp=censored_gamma_or_lognormal(tau,
//                          obsgroupedcatch1,1.e-2+groupedcatch(ig),censor,pf111);
                        dvar_vector vtmp=censored_gamma_or_lognormal(tau,
                          obsgroupedcatch1,rtc,censor,pf111); //NMD_12May2020

                        icount(it)++;
                        gicount(it,ig)++;
                        if (gicount(it,ig)>max_index)
                        {
                          cerr << "need to increase max_index" << endl;
                          ad_exit(1);
                        }

                        f_by_tag(it,icount(it))-=sum(log(1.e-10+vtmp));

                        int nelem=rtc.indexmax()-rtc.indexmin()+1;
                        empvar_nelem+=nelem;
                        dvector crtc=value(rtc);
                        dvector cobstagcatch=value(obsgroupedcatch1);
                        empvar+=sum(elem_div(square(cobstagcatch-crtc),crtc));   //NMD_12May2020
                        if (ppstf)   //NMD_30Sep2019
                        {
                          tmp_grouped_tag_like(it,ig,gicount(it,ig))=
                            value(f_by_tag(it,icount(it)));
                        }

                      }
                      break;
                    default:
                      cerr << "Illegal value for parest_flags(111) value = "
                         <<  parest_flags(111)  << endl;
                  }
                } 
              }
            }
          }
        }
      }
      //            ofstags << "it: " << it << setw(12) << sum(f_by_tag(it)(1,icount(it))) 
      //      	      << f_by_tag(it)(1,icount(it)) << endl;    ///NMD 22Aug2012
    }
    f+=sum(f_by_tag);

    if (ppstf)
    {
      if (allocated(ppstf->grouped_tag_like))
      {
        //ppstf->grouped_tag_like.deallocate();
      }

      if (!allocated(ppstf->grouped_tag_like))
        ppstf->grouped_tag_like.allocate(mmin,mmax,1,gmax32);

      for (int it=mmin;it<=mmax;it++)
      {
        for (int ig=1;ig<=gmax32;ig++)
        {
          ppstf->grouped_tag_like(it,ig)=
            tmp_grouped_tag_like(it,ig)(1,gicount(it,ig));
        }
      }
    }

    cout << " tagfit before pooling " << f << endl;
    //    ofstags << setw(12) << sum(f_by_tag(17)(1,icount(17))) 
    //            << f_by_tag(17)(1,icount(17)) << endl;
    //ofstags1 << setw(12) << f_by_tag(17)(1,icount(18)) << endl;

  }

  cout <<  "  Empirical variance calculation: " << empvar/empvar_nelem << endl; //NMD_12May2020
  
  const dvar_vector& a=fish_pars(4)+50.0001;
  //cout << "a = " << a << endl;
  dvar_vector& q=fish_pars(5);
  //cout << "q = " << q << endl;
  if (sum(column(fish_flags,34)))
    gp_pen=grouped_tag_reporting_rate_penalty();
  if (age_flags(105)>0)
  {
    dvariable tot_pen=age_flags(105)/100.*square(totobstags-totpredtags);
    cout << "Total observed tags = " << totobstags << endl;
    cout << "Total predicted tags = " << totpredtags << endl;
    cout << "Total tags penalty =   " << tot_pen << endl;
    f+=tot_pen;
  }

  if (age_flags(96))
  {
    if (mf_pvm->pvm_switch == 0 || mf_pvm->pvm_switch == 1)
    {
      dvariable tmp=fit_pooled_tag_returns();
      cout << " tagfit pooling " << tmp << endl;
      f+=tmp;
    }
  }
  f+=gp_pen;
  return f;
}


dvariable dvar_len_fish_stock_history::fit_tag_returns_mix(void)
{
  dvariable f=0.0;
  dvariable gp_pen=0.0;
  int fi;
  dvar_matrix obsgroupedcatch;
  dvar_matrix groupedcatch;

  ivector group_flags32=column(fish_flags,32);
  int gsum32=sum(group_flags32);
  int gmax32=Max(group_flags32);
  ivector gp_fish32(1,gmax32);
  ivector iflag(1,gmax32);
  if (gmax32)
  {
    gp_fish32.initialize();
    for (fi=1;fi<=num_fisheries;fi++)
      gp_fish32(group_flags32(fi))=fi;
    obsgroupedcatch.allocate(1,gmax32,1,nage);
    groupedcatch.allocate(1,gmax32,1,nage);
  }
  else
  {
    obsgroupedcatch.allocate(1,num_fisheries,1,nage);
    groupedcatch.allocate(1,num_fisheries,1,nage);
  }

  int mmin,mmax;
  if (mf_pvm->pvm_switch)
  {
    mmin=min_tag_group;
    mmax=max_tag_group;
  }
  else
  {
    mmin=1;
    mmax=num_tag_releases;
  }
  if (!gsum32)   // no fishery grouping
  {
    int rmin=1;
    int rmax=num_regions;
    for (int it=mmin;it<=mmax;it++)
    {
      if (pmsd)
      {
        int cs=pmsd->tag_species_pointer(it);
        rmin=pmsd->region_bounds(cs,1);
        rmax=pmsd->region_bounds(cs,2);
      }

      for (int ir=rmin;ir<=rmax;ir++)
      {
        int ub;
        if (!age_flags(96))
          ub=num_fish_periods(ir);
        else
          ub=terminal_tag_period(it,ir);
        for (int ip=initial_tag_period(it,ir);ip<=ub;ip++)
        {
          for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
          {
            for (int j=1;j<=nage;j++)
            {
              if (tagcatch(it,ir,ip,fi,j)<0)
              {
                cerr << tagcatch(it,ir,ip,fi,j) << endl;
              }
            }
            int pp1=parent(ir,ip,fi);
            dvar_vector rtc;
            if (age_flags(198))
              rtc=tag_rep_rate(it,ir,ip,fi)*tagcatch(it,ir,ip,fi);
            else
              rtc=rep_rate(ir,ip,fi)*tagcatch(it,ir,ip,fi);
            switch (parest_flags(111))
            {
              case 0:
                f+=sum(elem_div(square(obstagcatch(it,ir,ip,fi)-rtc),
                   .01+rtc));
                break;
              case 1:
                f-=sum(log(exp(elem_div(
                  square(obstagcatch(it,ir,ip,fi)-rtc),
                  -.01-rtc ))+0.01));
                break;
              case 2:
              {
                dvar_vector lrtc=log(1.e-11+rtc);
                  f+=sum(rtc) - obstagcatch(it,ir,ip,fi)*lrtc;
              for (int i=1;i<=nage;i++)
#if !defined(NO_MY_DOUBLE_TYPE)
                    f+=gammln(obstagcatch(it,ir,ip,fi,i)+1.0L);
#else
                    f+=gammln(obstagcatch(it,ir,ip,fi,i)+1.0);
#endif
                }
                break;
              case 3:
                {
                  dvariable a=fish_pars(4,pp1)+50.0001;
                  rtc+=1.e-11;
                  int ns=rtc.indexmax()-rtc.indexmin()+1; 
                  dvar_vector ap=log(a+rtc);
                  f+=sum(a*ap)+obstagcatch(it,ir,ip,fi)*ap;
                  f-=ns*a*log(a)+obstagcatch(it,ir,ip,fi)*log(rtc);
                  f-=sum(gammln(a+obstagcatch(it,ir,ip,fi)));
#if !defined(NO_MY_DOUBLE_TYPE)
                  f+=sum(gammln(obstagcatch(it,ir,ip,fi)+1.0L));
#else
                  f+=sum(gammln(obstagcatch(it,ir,ip,fi)+1.0));
#endif
                  f+=ns*gammln(a);
                }
  
                break;
              case 4:
                {
                  rtc+=1.e-11;
                  dvar_vector a=(fish_pars(4,fi)+50.0001)*rtc;
                  //dvar_vector a=(fish_pars(4,fi)+50.0001)*sqrt(rtc);
                  int pp1=parent(ir,ip,fi);
                  dvar_vector zerop=
                    fish_pars(5,pp1)/(1.0+fish_pars(6,pp1)*rtc);
                  dvar_vector zerop1=1-zerop;
                  rtc=elem_div(rtc,zerop1);
                  int ns=rtc.indexmax()-rtc.indexmin()+1; 
                  dvar_vector ap=log(a+rtc);
               
                  dvar_vector f1=elem_prod(a,ap)+elem_prod(obstagcatch(it,ir,ip,fi),ap);
                  f1-=elem_prod(a,log(a))+elem_prod(obstagcatch(it,ir,ip,fi),log(rtc));
                  f1-=gammln(a+obstagcatch(it,ir,ip,fi));
#if !defined(NO_MY_DOUBLE_TYPE)
                  f1+=gammln(obstagcatch(it,ir,ip,fi)+1.0L);
#else
                  f1+=gammln(obstagcatch(it,ir,ip,fi)+1.0);
#endif
                  f1+=gammln(a);
                  dvar_vector lzerop1=log(1.0-zerop);
              for (int j=1;j<=nage;j++)
                  {
                    if (obstagcatch(it,ir,ip,fi,j)>0)
                      f-=lzerop1(j)-f1(j);
                    else
                      f-=log(zerop(j)+zerop1(j)*exp(-f1(j))); 
                  }
                } 
                break;
            default:
              cerr << "Illegal value for parest_flags(111) value = "
                   <<  parest_flags(111)  << endl;
            }
          }
        }
      }
    }
  }
  else   // have grouping
  {
    int rmin=1;
    int rmax=num_regions;
    for (int it=mmin;it<=mmax;it++)
    {
      if (pmsd)
      {
        int cs=pmsd->tag_species_pointer(it);
        rmin=pmsd->region_bounds(cs,1);
        rmax=pmsd->region_bounds(cs,2);
      }

      for (int ir=rmin;ir<=rmax;ir++)
      {
        int ub;
        if (!age_flags(96))
          ub=num_fish_periods(ir);
        else
          ub=terminal_tag_period(it,ir);
        for (int ip=initial_tag_period(it,ir);ip<=ub;ip++)
        {
          obsgroupedcatch.initialize();
          groupedcatch.initialize();
          iflag.initialize();
          ivector& pi=parent(ir,ip);
          for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
          {
            int pp1=pi(fi);
            int pp=group_flags32(pp1);
            dvar_vector rtc;
            if (age_flags(198))
              rtc=tag_rep_rate(it,ir,ip,fi)*tagcatch(it,ir,ip,fi);
            else
              rtc=rep_rate(ir,ip,fi)*tagcatch(it,ir,ip,fi);
            iflag(pp)=1;
            int jmin=obstagcatch(it,ir,ip,fi).indexmin(); 
            obsgroupedcatch(pp)(jmin,nage)+=obstagcatch(it,ir,ip,fi); 
            int jmin1=rtc.indexmin(); 
            groupedcatch(pp)(jmin1,nage)+=rtc;
          }
  
          for (fi=1;fi<=gmax32;fi++)
          {
            if(iflag(fi))
            {
              switch (parest_flags(111))
              {
                case 0:
                  f+=sum(elem_div(
                    square(obsgroupedcatch(fi)-groupedcatch(fi)),
                   .01+groupedcatch(fi) ));
                  break;
                case 1:
                  f-=sum(log(exp(elem_div(
                    square(obsgroupedcatch(fi)-groupedcatch(fi)),
                    -.01-groupedcatch(fi) ))+0.01));
                  break;
                case 2:
                  {
                    dvar_vector lrtc=log(1.e-11+groupedcatch(fi));
                      f+=sum(groupedcatch(fi)) - obsgroupedcatch(fi)*lrtc;
                  for (int i=1;i<=nage;i++)
#if !defined(NO_MY_DOUBLE_TYPE)
                        f+=gammln(obsgroupedcatch(fi,i)+1.0L);
#else
                        f+=gammln(obsgroupedcatch(fi,i)+1.0);
#endif
                }
                  break;
                case 3:
                  {
                    dvar_vector rtc=1.e-11+groupedcatch(fi);
                    dvariable a=fish_pars(4,gp_fish32(fi))+50.0001;
                    int ns=rtc.indexmax()-rtc.indexmin()+1; 
                    dvar_vector ap=log(a+rtc);
                    f+=sum(a*ap)+obsgroupedcatch(fi)*ap;
                    f-=ns*a*log(a)+obsgroupedcatch(fi)*log(rtc);
                    f-=sum(gammln(a+obsgroupedcatch(fi)));
#if !defined(NO_MY_DOUBLE_TYPE)
                    f+=sum(gammln(obsgroupedcatch(fi)+1.0L));
#else
                    f+=sum(gammln(obsgroupedcatch(fi)+1.0));
#endif
                    f+=ns*gammln(a);
                  }
  
                  break;
                case 4:
                  {
                    dvar_vector rtc=groupedcatch(fi);
                    rtc+=1.e-11;
                    dvar_vector a=(fish_pars(4,gp_fish32(fi))+50.0001)*rtc;
                    //dvar_vector a=(fish_pars(4,gp_fish32(fi))+50.0001)*sqrt(rtc);
                    dvar_vector zerop= fish_pars(5,gp_fish32(fi)) /
                      (1.0+fish_pars(6,gp_fish32(fi))*rtc);
                    dvar_vector zerop1=1-zerop;
                    rtc=elem_div(rtc,zerop1);
                    int ns=rtc.indexmax()-rtc.indexmin()+1; 
                    dvar_vector ap=log(a+rtc);
               
                    dvar_vector f1=elem_prod(a,ap)
                       +elem_prod(obsgroupedcatch(fi),ap);
                    f1-=elem_prod(a,log(a))+elem_prod(obsgroupedcatch(fi),log(rtc));
                    f1-=gammln(a+obsgroupedcatch(fi));
#if !defined(NO_MY_DOUBLE_TYPE)
                    f1+=gammln(obsgroupedcatch(fi)+1.0L);
#else
                    f1+=gammln(obsgroupedcatch(fi)+1.0);
#endif
                    f1+=gammln(a);
                    dvar_vector lzerop1=log(1.0-zerop);
                for (int j=1;j<=nage;j++)
                    {
                      if (obsgroupedcatch(fi,j)>0)
                        f-=lzerop1(j)-f1(j);
                      else
                        f-=log(zerop(j)+zerop1(j)*exp(-f1(j))); 
                    }
                  } 
  
                  break;
                default:
                  cerr << "Illegal value for parest_flags(111) value = "
                     <<  parest_flags(111)  << endl;
              } 
            }
          }
        }
      }
    }
  }
  const dvar_vector& a=fish_pars(4)+50.0001;
  dvar_vector& q=fish_pars(5);
  //cout << "q = " << q << endl;
  if (sum(column(fish_flags,34)))
    gp_pen=grouped_tag_reporting_rate_penalty();

  if (age_flags(96))
  {
    if (mf_pvm->pvm_switch == 0 || mf_pvm->pvm_switch == 1)
    {
      f+fit_pooled_tag_returns_mix();
    }
  }
  f+=gp_pen;
  return f;
}

dvar_vector dvar_len_fish_stock_history::censored_gamma_or_lognormal
  (const prevariable & tau,dvar_vector& o,const dvar_vector& mu,MY_DOUBLE_TYPE e,
   int pf111)
{
  switch (pf111)
  {
  case 5:
    {
      return censored_gamma(tau,o,mu,e);
    }
    break;
  case 6:
    {
      return censored_gamma_mix(tau,o,mu,e);
    }
    break;
  case 7:
    {
      return censored_lognormal(tau,o,mu,e);
    }
    break;
  case 8:
    {
      return censored_lognormal_mix(tau,o,mu,e);
    }
    break;
  default:
   cerr << "Illegal value for pf111  value = " << pf111 << endl;
   ad_exit(1);
  }
}

dvariable  my_cumd_gamma(MY_DOUBLE_TYPE x,const prevariable& mu,const prevariable& tau)
{
  //return cumd_gamma(x/tau,mu/tau);
  return xcumd_gamma(x/tau,mu/tau);
}


 MY_DOUBLE_TYPE mysmoother(MY_DOUBLE_TYPE x,MY_DOUBLE_TYPE e,MY_DOUBLE_TYPE eps)
 {
   MY_DOUBLE_TYPE y=x-e;
   if (y<=-eps)
     return 0.0;
   else if (y>=eps)
     return 1.0;
   else if (y>0.0)
   {
     MY_DOUBLE_TYPE alpha=1.0/eps;
     MY_DOUBLE_TYPE beta=-0.5*alpha*alpha;
     return 0.5+alpha*y+beta*y*y;
   }
   else 
   {
     MY_DOUBLE_TYPE alpha=1.0/eps;
     MY_DOUBLE_TYPE beta=-0.5*alpha*alpha;
     MY_DOUBLE_TYPE z=-y;
     return 0.5-alpha*z-beta*z*z;
   }
}

 dvariable mysmoother(const prevariable& x,MY_DOUBLE_TYPE e,MY_DOUBLE_TYPE eps)
 {
   dvariable y=x-e;
   if (y<=-eps)
     return 0.0;
   else if (y>=eps)
     return 1.0;
   else if (y>0.0)
   {
     MY_DOUBLE_TYPE alpha=1.0/eps;
     MY_DOUBLE_TYPE beta=-0.5*square(alpha);
     return 0.5+alpha*y+beta*square(y);
   }
   else 
   {
     MY_DOUBLE_TYPE alpha=1.0/eps;
     MY_DOUBLE_TYPE beta=-0.5*square(alpha);
     dvariable z=-y;
     return 0.5-alpha*z-beta*square(z);
   }
}

dvar_vector dvar_len_fish_stock_history::censored_gamma_mix(
  const prevariable & tau,dvar_vector& o,const dvar_vector& mu,MY_DOUBLE_TYPE e)
{
  // this is a kludged mixture just to test the ideas
  const MY_DOUBLE_TYPE eps=0.1;
  int mmin=o.indexmin();
  int mmax=o.indexmax();
  dvar_vector tmp(mmin,mmax);
  tmp.initialize();
  MY_DOUBLE_TYPE pcon=0.05;
  if(parest_flags(326))
  {
    pcon=parest_flags(326)/1000.0;
  }
  MY_DOUBLE_TYPE tmult=9.0;
  for (int i=mmin;i<=mmax;i++)
  {
    if (o(i)<=e-eps)
    {
      tmp(i)=(1.0-pcon)*my_cumd_gamma(e,mu(i),tau)
        +pcon*my_cumd_gamma(e,mu(i),tmult*tau);
    }
    else if (o(i)>=e+eps)
    {
      tmp(i)=(1.0-pcon)*my_gamma_density(o(i),mu(i),tau)
        +pcon*my_gamma_density(o(i),mu(i),tmult*tau);
    }
    else
    {
      dvariable p=mysmoother(o(i),e,eps);
      //tmp(i)=my_cumd_gamma(e,mu(i),tau);
      tmp(i)=(1.0-p)*( (1.0-pcon)*my_cumd_gamma(e,mu(i),tau)
        +pcon*my_cumd_gamma(e,mu(i),tmult*tau))
        +p*((1.0-pcon)*my_gamma_density(o(i),mu(i),tau)
        +pcon*my_gamma_density(o(i),mu(i),tmult*tau));
    }
  }
  if (censored_gamma_report_flag)  //NMD_24Apr2020
  {
    //censor_report.initialize();
    for (int i=mmin;i<=mmax;i++)
    {
      if (o(i)<=e-eps)
      {
        censor_report(1)+=1;
      }
      else if (o(i)>=e+eps)
      {
        censor_report(2)+=1;
      }
      else
      {
        MY_DOUBLE_TYPE p=mysmoother(value(o(i)),e,eps);
        //tmp(i)=my_cumd_gamma(e,mu(i),tau);
        tmp(i)=(1.0-p)*( (1.0-pcon)*my_cumd_gamma(e,mu(i),tau)
          +pcon*my_cumd_gamma(e,mu(i),tmult*tau))
          +p*((1.0-pcon)*my_gamma_density(o(i),mu(i),tau)
          +pcon*my_gamma_density(o(i),mu(i),tmult*tau));
        censor_report(1)+=1.0-p;
        censor_report(2)+=p;
        censor_report(3)=e;
        censor_report(4)=eps;
      }
    }
    //ofstream ofs("censored_gamma_report");
    //ofs << "Proportion of bins censored " << censor_report(1)/sum(censor_report)
    //  << " Cutoff value for censoring " << e << " smoothed width for transition "
    //    << eps << endl;
  }  //NMD_24Apr2020
  return tmp;
}

dvar_vector dvar_len_fish_stock_history::censored_gamma(const prevariable & tau,dvar_vector& o,
  const dvar_vector& mu,MY_DOUBLE_TYPE e)
{
  const MY_DOUBLE_TYPE eps=0.1;
  int mmin=o.indexmin();
  int mmax=o.indexmax();
  dvar_vector tmp(mmin,mmax);
  tmp.initialize();
  for (int i=mmin;i<=mmax;i++)
  {
    if (o(i)<=e-eps)
    {
      tmp(i)=my_cumd_gamma(e,mu(i),tau);
    }
    else if (o(i)>=e+eps)
    {
      tmp(i)=my_gamma_density(o(i),mu(i),tau);
    }
    else
    {
      dvariable p=mysmoother(o(i),e,eps);
      //tmp(i)=my_cumd_gamma(e,mu(i),tau);
      tmp(i)=(1.0-p)*my_cumd_gamma(e,mu(i),tau)
        + p*my_gamma_density(o(i),mu(i),tau);
    }
  }
  if (censored_gamma_report_flag)
  {
    //censor_report.initialize();
    for (int i=mmin;i<=mmax;i++)
    {
      if (o(i)<=e-eps)
      {
        censor_report(1)+=1;
      }
      else if (o(i)>=e+eps)
      {
        censor_report(2)+=1;
      }
      else
      {
        MY_DOUBLE_TYPE p=mysmoother(value(o(i)),e,eps);
        //tmp(i)=my_cumd_gamma(e,mu(i),tau);
        tmp(i)=(1.0-p)*my_cumd_gamma(e,mu(i),tau)
          + p*my_gamma_density(o(i),mu(i),tau);
        censor_report(1)+=1.0-p;
        censor_report(2)+=p;
        censor_report(3)=e;
        censor_report(4)=eps;
      }
    }
    //ofstream ofs("censored_gamma_report");
    //ofs << "Proportion of bins censored " << censor_report(1)/sum(censor_report)
    //  << " Cutoff value for censoring " << e << " smoothed width for transition "
    //    << eps << endl;
  }
  return tmp;
}

dvariable xgamma_density(const prevariable& _x,const prevariable& _r,
  const  prevariable& _mu)
{
  ADUNCONST(prevariable,x)
  ADUNCONST(prevariable,r)
  ADUNCONST(prevariable,mu)
  //return exp(r*log(mu) + (r-1)*log(x)-mu*x-gammln(r));
  return exp(r*log(mu) + (r-1)*log(x)-mu*x-xgammln(r));
}


dvariable my_gamma_density(const prevariable&  x,const prevariable & mu,
  const prevariable & tau)
{
  return xgamma_density(x,mu/tau,1.0/tau);
  //return gamma_density(x,mu/tau,1.0/tau);
}

dvariable xgammln(const prevariable& _x)
{
   ADUNCONST(dvariable, x)
   dvariable y=0.0;
   init_df1_two_variable xx(x);
   init_df1_two_variable yy(y);
   df1_two_variable zz=xgammln(xx);
   dvariable tmp;
   tmp=zz;
   return tmp;
}

MY_DOUBLE_TYPE my_lognormal_density(MY_DOUBLE_TYPE x,MY_DOUBLE_TYPE mu,MY_DOUBLE_TYPE sigma)
{
  const MY_DOUBLE_TYPE c=sqrt(2.0*3.1415926535);
  return 1.0/
    (x*sigma*c)*exp(-0.5*square((log(x)-(log(mu)-0.5*square(sigma))/sigma)));
  //return 1.0/(x*sigma*c)*exp(-0.5*square((log(x)-log(mu))/sigma));
}

dvariable my_lognormal_density(const prevariable& x,const prevariable&  mu,
  const prevariable& sigma)
{
  const MY_DOUBLE_TYPE c=sqrt(2.0*3.1415926535);
  return 1.0/
    (x*sigma*c)*exp(-0.5*square( (log(x)-(log(mu)-0.5*square(sigma))/sigma)));
  //return 1.0/(x*sigma*c)*exp(-0.5*square( (log(x)-log(mu)) /sigma));
}

dvariable erf(const prevariable& _x);

dvariable my_cumd_lognormal(MY_DOUBLE_TYPE x,const prevariable& mu,
  const prevariable& sigma)
{
  const MY_DOUBLE_TYPE rt=sqrt(2.0);
  return 0.5*(1.0+erf((log(x)-(log(mu)-0.5*square(sigma))/(rt*sigma))));
  //return 0.5*(1.0+erf((log(x)-log(mu))/(rt*sigma)));
}

dvar_vector dvar_len_fish_stock_history::censored_lognormal_mix(
  const prevariable & sigma,dvar_vector& o,const dvar_vector& mu,MY_DOUBLE_TYPE e)
{
  // this is a kludged mixture just to test the ideas
  // !! VVVTEST
  dvariable emu=fish_pars(29,1);
  const MY_DOUBLE_TYPE eps=0.1;
  int mmin=o.indexmin();
  int mmax=o.indexmax();
  dvar_vector tmp(mmin,mmax);
  tmp.initialize();
  MY_DOUBLE_TYPE pcon=0.05;
  if(parest_flags(326))
  {
    pcon=parest_flags(326)/1000.0;
  }
  MY_DOUBLE_TYPE smult=3.0;
  // !! VVVTEST
  for (int i=mmin;i<=mmax;i++)
  {
    // !! VVVTEST
    //dvariable sigmai=sigma*pow(mu(i),emu);
    if (o(i)<=e-eps)
    {
      //tmp(i)=(1.0-pcon)*my_cumd_gamma(e,mu(i),tau)
      //  +pcon*my_cumd_gamma(e,mu(i),tmult*tau);
      tmp(i)=(1.0-pcon)*my_cumd_lognormal(e,mu(i),sigma)
        +pcon*my_cumd_lognormal(e,mu(i),smult*sigma);
    }
    else if (o(i)>=e+eps)
    {
      //tmp(i)=(1.0-pcon)*my_gamma_density(o(i),mu(i),tau)
      //  +pcon*my_gamma_density(o(i),mu(i),tmult*tau);
      tmp(i)=(1.0-pcon)*my_lognormal_density(o(i),mu(i),sigma)
        +pcon*my_lognormal_density(o(i),mu(i),smult*sigma);
    }
    else
    {
      dvariable p=mysmoother(o(i),e,eps);
      
      tmp(i)=(1.0-p)*( (1.0-pcon)*my_cumd_lognormal(e,mu(i),sigma)
        +pcon*my_cumd_lognormal(e,mu(i),smult*sigma))
        +p*((1.0-pcon)*my_lognormal_density(o(i),mu(i),sigma)
        +pcon*my_lognormal_density(o(i),mu(i),smult*sigma));
    }
  }
  return tmp;
}


dvar_vector censored_lognormal(const prevariable & sigma,dvar_vector& o,
  const dvar_vector& mu,MY_DOUBLE_TYPE e)
{
  int mmin=o.indexmin();
  int mmax=o.indexmax();
  dvar_vector tmp(mmin,mmax);
  const MY_DOUBLE_TYPE eps=0.1;
  for (int i=mmin;i<=mmax;i++)
  {
    if (o(i)<=e-eps)
    {
      tmp(i)=my_cumd_lognormal(e,mu(i),sigma);
    }
    else if (o(i)>=e+eps)
    {
      tmp(i)=my_lognormal_density(o(i),mu(i),sigma);
    }
    else
    {
      dvariable p=mysmoother(o(i),e,eps);
      tmp(i)=(1.0-p)*my_cumd_lognormal(e,mu(i),sigma)
        +p*my_lognormal_density(o(i),mu(i),sigma);
    }
  }
  return tmp;
}

dvariable erf(const prevariable& _x)
{
   ADUNCONST(dvariable, x)
   dvariable x2=x*x;
   dvariable u=0.5;
   init_df1_two_variable vx2(x2);
   init_df1_two_variable vu(u);
   df1_two_variable vy;
   vy = cumd_gamma(vx2,vu);
   dvariable z;
   z = vy;
   return z;
}
