/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

dvariable dvar_len_fish_stock_history::fit_tag_returns2(void)
{
  dvariable f=0.0;
  int mmin,mmax;
#if defined(USE_ADPVM)
  if (mf_pvm->pvm_switch)
  {
    mmin=min_tag_group;
    mmax=max_tag_group;
  } 
  else
#endif
  {
    mmin=1;
    mmax=num_tag_releases;
  }
  dvar_matrix f_by_tag(mmin,mmax,1,10000);
  ivector icount(mmin,mmax);
  icount.initialize();
  f_by_tag.initialize();
  int rmin=1;
  int rmax=num_regions;
  ivector group_flags32=column(fish_flags,32);
  int gsum32=sum(group_flags32);
  int gmax32=Max(group_flags32);
  ivector gp_fish32(1,gmax32);

  if (gmax32)
  {
    gp_fish32.initialize();
    for (int fi=1;fi<=num_fisheries;fi++)
      gp_fish32(group_flags32(fi))=fi;
  }

  for (int itt=1;itt<=num_tag_releases;itt++)
  {
    int ngt=pmsd->tag_species_index(itt).indexmax();
   
    dvar4_array xxx(1,num_regions,1,
       num_real_fish_periods,1,gmax32,1,nage);
    dvar4_array yyy(1,num_regions,1,
       num_real_fish_periods,1,gmax32,1,nage);
    i3_array iflag2(1,num_regions,1,num_fish_periods,1,gmax32);
    xxx.initialize();
    yyy.initialize();
    iflag2.initialize();

        dvariable sumtmp1; //NMD
        dvariable sumtmp2; //NMD
	sumtmp1=0;
	sumtmp2=0;
    for (int ic=1;ic<=ngt;ic++)
    {
      int it=pmsd->tag_species_index(itt,ic);
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
        if (ub>num_real_fish_periods(ir))
        {
          ub=num_real_fish_periods(ir);
        }
        for (int ip=initial_tag_period(it,ir);ip<=ub;ip++)
        {
          if (num_fish_incidents(ir,ip)>0)
          { 
            ivector& pi=parent(ir,ip);
            dvector xtmp(1,9);
            for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
            {
              int pp1=pi(fi);
              int pp=group_flags32(pp1);

              dvar_vector rtc;
              if (age_flags(198))
                rtc=tag_rep_rate(it,ir,ip,fi)*tagcatch(it,ir,ip,fi);
              else
                rtc=rep_rate(ir,ip,fi)*tagcatch(it,ir,ip,fi);

              int jmin=rtc.indexmin();
              //totpredtags+=sum(rtc);
//              if (itt==pmsd->tag_species_index(itt,1))  //NMD5Dec2013
              if (it==pmsd->tag_species_index(itt,1))
              {
                int mmin=obstagcatch(it,ir,ip,fi).indexmin();
                int mmax=obstagcatch(it,ir,ip,fi).indexmax();
                yyy(ir,ip,pp)(mmin,mmax)+=obstagcatch(it,ir,ip,fi); 
              }

              int it1=pmsd->tag_species_index(itt,1);
              int it2=pmsd->tag_species_index(itt,ic);
              int sp1=pmsd->tag_species_pointer(it1);
              int sp2=pmsd->tag_species_pointer(it2);

              int ir1=(sp1-sp2)*pmsd->num_real_regions+ir;
              int mmin=rtc.indexmin();
              int mmax=rtc.indexmax();
              xxx(ir1,ip,pp)(mmin,mmax)+=rtc;
              iflag2(ir1,ip,pp)=1;
              sumtmp1+=sum(obstagcatch(it,ir,ip,fi)); //NMD
              sumtmp2+=sum(xxx(ir1,ip,pp)); //NMD
            }
          }
        }
      }
//      if (itt==2)  //NMD
//      {
      //        cout << " Debug check of tags" << endl;
      //        cout << " itt:  " << itt << "  ic: " << ic << " obstagcatch " << sumtmp1 << "  rtc  " << sumtmp2 << endl;
//      }

    }

   
    if (itt==pmsd->tag_species_index(itt,1))
    {
      int sp=pmsd->tag_species_pointer(itt);
      
      ivector rb=get_region_bounds(sp);
      for (int ir=rb(1);ir<=rb(2);ir++)
      {
        int ub;
        if (!age_flags(96))
          ub=num_real_fish_periods(ir);
        else
          ub=terminal_tag_period(itt,ir);
        if (ub>num_real_fish_periods(ir))
        {
          ub=num_real_fish_periods(ir);
        }
        for (int ip=initial_tag_period(itt,ir);ip<=ub;ip++)
        {
          if (num_fish_incidents(ir,ip)>0)
          { 
            int ig;
            for (ig=1;ig<=gmax32;ig++)
            {
              if (iflag2(ir,ip,ig))
              {
                if (grouped_fishery_projection_flag(ir,ip,ig)==0)
                {
                  switch (parest_flags(111))
                  {
                    case 0:
                     // f+=sum(elem_div(
                      //  square(obsgroupedcatch(ig)-groupedcatch(ig)),
                       //.01+groupedcatch(ig) ));
                      break;
                    case 1:
                      //f-=sum(log(exp(elem_div(
                        //square(obsgroupedcatch(ig)-groupedcatch(ig)),
                        //-.01-groupedcatch(ig) ))+0.01));
                      break;
                    case 2:
                    {
                      //dvar_vector lrtc=log(1.e-11+groupedcatch(ig));
                        //f+=sum(groupedcatch(ig)) - obsgroupedcatch(ig)*lrtc;
                      //for (int i=1;i<=nage;i++)
#if !defined(NO_MY_DOUBLE_TYPE)
                        //  f+=gammln(obsgroupedcatch(ig,i)+1.0L);
#else
                        //  f+=gammln(obsgroupedcatch(ig,i)+1.0);
#endif
                    }
                      break;
                    case 3:
                      {
                       // dvar_vector rtc=1.e-11+groupedcatch(ig);
                       // dvariable a=fish_pars(4,gp_fish32(ig))+50.0001;
                        //int ns=rtc.indexmax()-rtc.indexmin()+1; 
                        //dvar_vector ap=log(a+rtc);
                        //f+=sum(a*ap)+obsgroupedcatch(ig)*ap;
                        //f-=ns*a*log(a)+obsgroupedcatch(ig)*log(rtc);
                        //f-=sum(gammln(a+obsgroupedcatch(ig)));
#if !defined(NO_MY_DOUBLE_TYPE)
                        //f+=sum(gammln(obsgroupedcatch(ig)+1.0L));
#else
                        //f+=sum(gammln(obsgroupedcatch(ig)+1.0));
#endif
                        //f+=ns*gammln(a);
                      }
    
                      break;
                    case 4:
                      {
                        //dvar_vector obsgroupedcatch1=obsgroupedcatch(ig);
                        dvar_vector obsgroupedcatch1=yyy(ir,ip,ig);
                        dvar_vector rtc=1.e-11+xxx(ir,ip,ig);
                        //dvar_vector rtc=1.e-11+groupedcatch(ig);
                        dvar_vector a=(fish_pars(4,gp_fish32(ig)) +50.0001)*rtc;
                        int ns=rtc.indexmax()-rtc.indexmin()+1; 
                        dvar_vector ap=log(a+rtc);
                        icount(itt)++;
                        f_by_tag(itt,icount(itt))+=a*ap+obsgroupedcatch1*ap;
                        f_by_tag(itt,icount(itt))-=a*log(a)+obsgroupedcatch1*log(rtc);
                        f_by_tag(itt,icount(itt))-=sum(gammln(a+obsgroupedcatch1));
#if !defined(NO_MY_DOUBLE_TYPE)
                        f_by_tag(itt,icount(itt))+=sum(gammln(obsgroupedcatch1+1.0L));
#else
                        f_by_tag(itt,icount(itt))+=sum(gammln(obsgroupedcatch1+1.0));
#endif
                        f_by_tag(itt,icount(itt))+=sum(gammln(a));

                        if (itt==1 && icount(itt)==1)
                        {
  //                           cout << "HERE" << endl;
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
    }
  }
  dvariable tmp=sum(f_by_tag);
  f+=tmp;
  return f;
}


