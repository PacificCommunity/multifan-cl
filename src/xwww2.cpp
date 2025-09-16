/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
//# define COUNT_FISH 
extern mf_pvm_manager * mf_pvm;
dvar_matrix fast_diffusion_calcs(int nage,int num_regions,dvar3_array& N,
  dvar3_array& Dad,ivector& rip,int maxage,MY_DOUBLE_TYPE delta)
{ 
  //const MY_DOUBLE_TYPE delta=1.e-6;
  const MY_DOUBLE_TYPE nd=num_regions*delta;
  int j;
  if (maxage) nage=min(nage,maxage);
  int jmin=N(1,rip(1)).indexmin();
  dvar_matrix EN(jmin,nage,1,num_regions);
  dvar_matrix tmp(jmin,nage,1,num_regions);
  tmp.initialize();

  for (j=jmin;j<=nage;j++)
  {
    for (int is=1;is<=num_regions;is++)
    {
      EN(j,is)=mfexp(N(is,rip(is),j));
    }
  }
  for (j=jmin;j<=nage;j++)
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int is=1;is<=num_regions;is++)
      {
        tmp(j,ir)+=Dad(j)(ir,is)*EN(j,is);
      }
    }
    dvariable ss=sum(tmp(j));
    tmp(j)+=delta;
    tmp(j)*=ss/(ss+nd);
  }
  return trans(tmp);
}

dvariable dvar_fish_stock_history::xtag_catch_equations_calc_pooled(dvar_vector& sv)
{
  int break_flag=0;
  dvariable ffpen=0.0;
  epooled_tagnum_fish_recr.initialize();
  pooledtagN=-20;
  int it;
  int mmin,mmax;
  //ofstream xpofs("xww2_movestuff");
  //ivector tmp_mp(1,num_regions);
  /*
  ivector bug_tmp_mp(1,num_regions);
  ivector tmp_yr(1,num_regions);
  ivector tmp_mn(1,num_regions);
  ivector tmp_ip(1,num_regions);
  bug_tmp_mp.initialize();
  //tmp_mp.initialize();
  tmp_yr.initialize();
  tmp_mn.initialize();
  tmp_ip.initialize();
  */
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
  for (int itt=mmin;itt<=mmax;itt++)
  {
    int it=itt;
    ivector tmp_mp(tag_region_bounds(1,it),tag_region_bounds(2,it));
    ivector tmp_ip(tag_region_bounds(1,it),tag_region_bounds(2,it));
    ivector tmp_yr(tag_region_bounds(1,it),tag_region_bounds(2,it));
    ivector tmp_mn(tag_region_bounds(1,it),tag_region_bounds(2,it));
    tmp_mp.initialize();
    tmp_yr.initialize();
    tmp_mn.initialize();
    tmp_ip.initialize();
   int ng=get_nage_tag(it);

    get_initial_tag_population(sv,it);

    int current_year=0;
    ivector rip(tag_region_bounds(1,it),tag_region_bounds(2,it));

    current_year=tag_year(it);
    rip=initial_tag_period(it);

    int finished_flag=0;
    int eflag=0;
    int ir;
    for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
    {
      tagnum_fish(it,ir)=-20.0;
    }
    for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
    {
      tagnum_fish(it,ir,rip(ir))=tagN(it,ir,year(ir,rip(ir)));
    }

    do
    {
      finished_flag=1;
      eflag=0;
      for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
      {
        break_flag=0;
        int& ip=rip(ir);
        int ttp;
        ttp=terminal_tag_period(it,ir);
        do
        {
          if ( ip>=num_fish_periods(ir) || ip>ttp ) 
          {
            eflag=1;
            break;
          }
          finished_flag=0;
          if (ip < ttp)
          {
            if (year(ir,ip+1)==year(ir,ip))
            { 
              if (!num_fish_incidents(ir,ip) || !tag_flags(it,1) 
                 || ip >= initial_tag_period(it,ir)+tag_flags(it,1))
              {
                int mmin=tagnum_fish(it,ir,ip).indexmin();
                int mmax=tagnum_fish(it,ir,ip).indexmax();
                if (parest_flags(360)==0)
                {
                  tagnum_fish(it,ir,ip+1)=
                    tagnum_fish(it,ir,ip)-tot_mort(ir,ip)(mmin,mmax);
                }
                else
                {
                  tagnum_fish(it,ir,ip+1)=
                    tagnum_fish(it,ir,ip)-tot_mort(ir,ip)(mmin,mmax)
                     -tagmort(it);
                }
              }
              else
              {
                do_newton_raphson_for_tags(it,ir,ip,ffpen);
//                tagnum_fish(it,ir,ip+1)=tagnum_fish(it,ir,ip)-nrtm(it,ir,ip);
                if (parest_flags(360)==0)   //NMD_7feb2022
                {
                  tagnum_fish(it,ir,ip+1)=tagnum_fish(it,ir,ip)-nrtm(it,ir,ip);
                }
                else
                {
                  tagnum_fish(it,ir,ip+1)=tagnum_fish(it,ir,ip)-nrtm(it,ir,ip)
                    -tagmort(it);
                }
              }
// Correction - remove this unnecessary step that over-writes the mixing calcs
              /*         //NMD_14feb2022
              int mmin=tagnum_fish(it,ir,ip).indexmin();
              int mmax=tagnum_fish(it,ir,ip).indexmax();
              if (parest_flags(360)==0)
              {
                tagnum_fish(it,ir,ip+1)=
                  tagnum_fish(it,ir,ip)-tot_mort(ir,ip)(mmin,mmax);
              }
              else
              {
                tagnum_fish(it,ir,ip+1)=
                  tagnum_fish(it,ir,ip)-tot_mort(ir,ip)(mmin,mmax)
                     -tagmort(it);
              }
              */
            }
            else
            {
              if (!num_fish_incidents(ir,ip) || !tag_flags(it,1) 
                || ip >= initial_tag_period(it,ir)+tag_flags(it,1))
              {
                dvar_vector  tnf;
                dvar_vector  tnfp;
                tnf=tagnum_fish(it,ir,ip);
                tnfp=tagnum_fish(it,ir,ip+1);
                int jmin=tnf.indexmin();
                if (nage>2 && jmin<nage-1)
                {
                  if (parest_flags(360)==0)
                  {
                    --tnfp(jmin+1,nage-1)
                      =tnf(jmin,nage-2)
                      -tot_mort(ir,ip)(jmin,nage-2);
                  }
                  else
                  {
                    --tnfp(jmin+1,nage-1)
                      =tnf(jmin,nage-2)
                      -tot_mort(ir,ip)(jmin,nage-2)-tagmort(it);
                  }
                }

                if (jmin<nage)
                {
                  if (parest_flags(360)==0)
                  {
                    tnfp(nage) 
                      = log(1.e-20
                        +mfexp(tnf(nage)-tot_mort(ir,ip,nage))
                        + mfexp(tnf(nage-1)-tot_mort(ir,ip,nage-1)));
                  }
                  else
                  {
                    tnfp(nage) 
                      = log(1.e-20
                        +mfexp(tnf(nage)-tot_mort(ir,ip,nage)-tagmort(it))
                        + mfexp(tnf(nage-1)-tot_mort(ir,ip,nage-1)
                         -tagmort(it)));
                  }
                }
                else
                {
                  if (parest_flags(360)==0)
                  {
                    tnfp(nage) 
                      = tnf(nage)-tot_mort(ir,ip,nage);
                  }
                  else
                  {
                    tnfp(nage) 
                      = tnf(nage)-tot_mort(ir,ip,nage)-tagmort(it);
                  }
                }
              }
              else
              {
                do_newton_raphson_for_tags(it,ir,ip,ffpen);

                int minage=nrtm(it,ir,ip).indexmin();
                // cout << " XXXXXXXXXXXXXXXX  test change here " << endl;
//                if (nage>minage+1)
//                  --tagnum_fish(it,ir,ip+1)(minage+1,nage-1)
//                   =log(1.e-5+mfexp(tagnum_fish(it,ir,ip)(minage,nage-2)
//                   -nrtm(it,ir,ip)(minage,nage-2)));
//                tagnum_fish(it,ir,ip+1,nage) 
//                  = log( 1.e-20
//                     +(mfexp(tagnum_fish(it,ir,ip,nage)-nrtm(it,ir,ip,nage))
//                      )
//                  + mfexp(tagnum_fish(it,ir,ip,nage-1)-nrtm(it,ir,ip,nage-1)));
                if (parest_flags(360)==0)   //NMD_7feb2022
                {
                  if (nage>minage+1)
                    --tagnum_fish(it,ir,ip+1)(minage+1,nage-1)
                      =log(1.e-5+mfexp(tagnum_fish(it,ir,ip)(minage,nage-2)
                      -nrtm(it,ir,ip)(minage,nage-2)));
                  tagnum_fish(it,ir,ip+1,nage) 
                    = log( 1.e-20
                       +(mfexp(tagnum_fish(it,ir,ip,nage)-nrtm(it,ir,ip,nage))
                        )
                    + mfexp(tagnum_fish(it,ir,ip,nage-1)
                    -nrtm(it,ir,ip,nage-1)));
                }
                else
                {
                  if (nage>minage+1)
                    --tagnum_fish(it,ir,ip+1)(minage+1,nage-1)
                      =log(1.e-5+mfexp(tagnum_fish(it,ir,ip)(minage,nage-2)
                      -nrtm(it,ir,ip)(minage,nage-2)-tagmort(it)));
                  tagnum_fish(it,ir,ip+1,nage) 
                    = log( 1.e-20
                      +(mfexp(tagnum_fish(it,ir,ip,nage)-nrtm(it,ir,ip,nage)
      	              -tagmort(it)))
                      + mfexp(tagnum_fish(it,ir,ip,nage-1)
                      -nrtm(it,ir,ip,nage-1)-tagmort(it)));
                }
              }
            }
          }
          else
          {
            eflag=1;
            if (!tag_flags(it,1) || ip >= initial_tag_period(it,ir)
              +tag_flags(it,1))
            {
              int jmin;
              jmin=tagnum_fish(it,ir,ip).indexmin();
              
              if (nage>2 && jmin<nage-1)
              {
                if (parest_flags(360)==0)
                {
                  --epooled_tagnum_fish_recr(ir,ip+1)(jmin+1,nage-1)
                    +=mfexp(tagnum_fish(it,ir,ip)(jmin,nage-2)
                    -tot_mort(ir,ip)(jmin,nage-2));
                }
                else
                {
                  --epooled_tagnum_fish_recr(ir,ip+1)(jmin+1,nage-1)
                    +=mfexp(tagnum_fish(it,ir,ip)(jmin,nage-2)
                    -tot_mort(ir,ip)(jmin,nage-2)-tagmort(it));
                }
              }
    
              if (nage>2 && jmin<nage)
              {
                if (parest_flags(360)==0)
                {
                  epooled_tagnum_fish_recr(ir,ip+1,nage) 
                    += mfexp(tagnum_fish(it,ir,ip,nage)-tot_mort(ir,ip,nage))
                    +mfexp(tagnum_fish(it,ir,ip,nage-1)-tot_mort(ir,ip,nage-1));
                }
                else
                {
                  epooled_tagnum_fish_recr(ir,ip+1,nage) 
                    += mfexp(tagnum_fish(it,ir,ip,nage)
                        -tot_mort(ir,ip,nage)-tagmort(it))
                    +mfexp(tagnum_fish(it,ir,ip,nage-1)
                        -tot_mort(ir,ip,nage-1)-tagmort(it));
                }
              }
              else
              {
                if (parest_flags(360)==0)
                {
                  epooled_tagnum_fish_recr(ir,ip+1,nage) 
                    += mfexp(tagnum_fish(it,ir,ip,nage)-tot_mort(ir,ip,nage));
                }
                else
                {
                  epooled_tagnum_fish_recr(ir,ip+1,nage) 
                    += mfexp(tagnum_fish(it,ir,ip,nage)
                       -tot_mort(ir,ip,nage)-tagmort(it));
                }
              }
            }
            else
            {
              do_newton_raphson_for_tags(it,ir,ip,ffpen);

//              if (nage>2)
//                --epooled_tagnum_fish_recr(ir,ip+1)(2,nage-1)
//                 +=mfexp(tagnum_fish(it,ir,ip)(1,nage-2)
//                 -nrtm(it,ir,ip)(1,nage-2));
//   
//              epooled_tagnum_fish_recr(ir,ip+1,nage) 
//                += mfexp(tagnum_fish(it,ir,ip,nage)-nrtm(it,ir,ip,nage))
//                + mfexp(tagnum_fish(it,ir,ip,nage-1)-nrtm(it,ir,ip,nage-1));

	      if (parest_flags(360)==0) //NMD_7feb2022
              {
                if (nage>2)
                  --epooled_tagnum_fish_recr(ir,ip+1)(2,nage-1)
                   +=mfexp(tagnum_fish(it,ir,ip)(1,nage-2)
                   -nrtm(it,ir,ip)(1,nage-2));
   
                epooled_tagnum_fish_recr(ir,ip+1,nage) 
                  += mfexp(tagnum_fish(it,ir,ip,nage)-nrtm(it,ir,ip,nage))
                  + mfexp(tagnum_fish(it,ir,ip,nage-1)-nrtm(it,ir,ip,nage-1));
              }
	      else
              {
                if (nage>2)
                  --epooled_tagnum_fish_recr(ir,ip+1)(2,nage-1)
                   +=mfexp(tagnum_fish(it,ir,ip)(1,nage-2)
                   -nrtm(it,ir,ip)(1,nage-2)-tagmort(it));
   
                epooled_tagnum_fish_recr(ir,ip+1,nage) 
                  += mfexp(tagnum_fish(it,ir,ip,nage)-nrtm(it,ir,ip,nage)
                  -tagmort(it))
                  + mfexp(tagnum_fish(it,ir,ip,nage-1)-nrtm(it,ir,ip,nage-1)
                  -tagmort(it));
              }
            }
          }
          if (move_flags(ir,ip+1))
          {
            tmp_ip(ir)=ip+1;
            tmp_yr(ir)=year(ir,ip+1);
            tmp_mn(ir)=month(ir,ip+1);
            tmp_mp(ir)=move_index(ir,ip+1);
          
            if (ir==num_regions)
            {
              //xpofs << "tag release group " << it << "  ";
              //print_movement_stuff(xpofs,tmp_mp,tmp_yr,tmp_mn,tmp_ip,
              //  num_regions);
            }
            break_flag=1;
          } 
          ip++;
        }
        while (!break_flag); 
      }
      // move the tags  -- need to deal with pooled ones
      {
        if (num_regions>1 && !eflag)
        {
          check_sanity(tmp_mp,tagnum_fish(it),rip,it,*this);
          if (pmsd) pmsd->tag_index=it;
          dvar_matrix tmp=fast_diffusion_calcs(ng,num_regions,
            tagnum_fish(it),Dad(tmp_mp(tmp_mp.indexmin())),rip,0,pmsd);
            //tagnum_fish(it),Dad(tmp_mp(1)),rip);
          if (pmsd) pmsd->tag_index=0;
          //check_sanity(tmp_mp);
          //int tp=tmp_mp(1);
          int tp=tmp_mp(tag_region_bounds(1,it));
//          if (parest_flags(356))
//            tp=bug_tmp_mp(1);
          //dvar_matrix tmp=
          //  fast_diffusion_calcs(nage,num_regions,tagnum_fish(it),
          //  Dad(tp),rip,0,pmsd);
          for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
          {
            tagnum_fish(it,ir,rip(ir))=log(5.e-10+tmp(ir));
          }
        }
      }
      tmp_mp=0;
    } // need to decide when to quit
    while (!finished_flag);
  }

  const MY_DOUBLE_TYPE one_plus=1.e0+1.e-12;

  int id=0;
  for (int itt=mmin;itt<=mmax;itt++)
  {
    int it=itt;
    imatrix * pitp;
    imatrix * pttp;
    pitp=&initial_tag_period;
    pttp=&terminal_tag_period;

    id=0;
    for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
    {
      for (int ip=(*pitp)(it,ir);ip<=(*pttp)(it,ir);ip++)
      {
        if (!tag_flags(it,1) || 
          ip >= initial_tag_period(it,ir)+tag_flags(it,1))
        {
          for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
          {        
            dvar_vector& tc=tagcatch(it,ir,ip,fi);
            int jmin=tc.indexmin();
            if (parest_flags(360)==0)
            {
              tc=exp(fish_mort_calcs(ir,ip,fi)(jmin,nage)
               +tagnum_fish(it,ir,ip));
            }
            else
            {
              tc=exp(tag_fish_mort_calcs(1,ir,ip,fi)(jmin,nage)
               +tagnum_fish(it,ir,ip));
            }
          }
        }
        else
        {
          for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
          {        
            dvar_vector& tc=tagcatch(it,ir,ip,fi);
            int jmin=tc.indexmin();
            nrsurv(it,ir,ip)(jmin,nage);
//            tc=exp(nrfm(it,ir,ip,fi)-log(1.e-10+nrtm(it,ir,ip))+
//              log(one_plus-nrsurv(it,ir,ip)(jmin,nage))+tagnum_fish(it,ir,ip));
            dvar_vector tmp1;
            dvar_vector tmp2;
            if (parest_flags(360)==0)
            {
              tmp1=log(1.e-10+nrtm(it,ir,ip));
              tmp2=log(one_plus-nrsurv(it,ir,ip)(jmin,nage));
            }
            else
            {
              tmp1=log(1.e-10+nrtm(it,ir,ip)) + tagmort(it);
              tmp2=log(one_plus-nrsurv(it,ir,ip)(jmin,nage)) - tagmort(it);
            }		
            tc=exp(nrfm(it,ir,ip,fi) - tmp1 + tmp2 + tagnum_fish(it,ir,ip));
          }
        }
      }
    }
  }
/*
# if defined(COUNT_FISH)
  print_tag_accounting_info();
  exit(1);
#endif
*/
 /*
  ofstream ofs("tagnumfish");
  ofs << tagnum_fish << endl;
  ad_exit(1);
 */
  return ffpen;
}
