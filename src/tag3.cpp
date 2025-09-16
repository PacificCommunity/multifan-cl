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
  dvariable mfexp(BOR_CONST prevariable& );
  //dvar_vector mfexp(dvar_vector& );
  //dvector mfexp(dvector& );
#endif

#ifdef __MSVC32__
  dvariable age_at_length_calcxx(MY_DOUBLE_TYPE& v,dvar_vector& gml,int nslots);
#else
  dvariable age_at_length_calcxx(const dvariable& v,dvar_vector& gml,int nslots);
#endif
extern mf_pvm_manager * mf_pvm;

MY_DOUBLE_TYPE maxZ[5];

ofstream ofsdd("XXX");

dvar_matrix ADJ_get_jac2(int nfi,int jmin,int nage,dvar_vector _q,
  dvar_matrix _sel,const dvar_vector& _N,dvar_vector& _M,MY_DOUBLE_TYPE rmax);

dvar_matrix get_jac2(int nfi,int jmin,int nage,dvar_vector q,
  dvar_matrix sel,const dvar_vector& N,dvar_vector& M,MY_DOUBLE_TYPE rmax)
{
  int i,j,k;
  dvar_matrix J(1,nfi,1,nfi);
  J.initialize();
  dvar_matrix F(1,nfi);
  dvar_vector Z(jmin,nage);
  for (i=1;i<=nfi;i++)
  {
    F(i)=q(i)*sel(i);
  }
  Z=colsum(F)+M;
  dvar_vector S=exp(-Z);
  dvar_vector NZ=elem_div(N,Z);
  for (i=1;i<=nfi;i++)
  {
    for (k=1;k<=nfi;k++)
    {
      for (j=jmin;j<=nage;j++)
      {
        dvariable tmp=0;
        tmp-=F(i,j)/Z(j)*sel(k,j);
        if (i==k)
        {
          tmp+=sel(k,j);
        }
        tmp*=(1.0-S(j));
        if (Z(j)<=rmax)
        {
          tmp+=F(i,j)*S(j)*sel(k,j);  //baranov
        }
        else
        {
          tmp+=F(i,j)*exp(-rmax)*sel(k,j);  //kludged
        }
        J(i,k)+=tmp*NZ(j);
      }
    }
  }
  return J;
}

/*
dvar_matrix get_jacy(int nfi,int jmin,int nage,dvar_vector q,
  dvar_matrix sel,const dvar_vector& N,dvar_vector& M,MY_DOUBLE_TYPE rmax)
{
  int i,j,k;
  dvar_matrix J(1,nfi,1,nfi);
  J.initialize();
  dvar_matrix F(1,nfi);
  dvar_vector Z(jmin,nage);
  for (i=1;i<=nfi;i++)
  {
    F(i)=q(i)*sel(i);
  }
  Z=colsum(F)+M;
  dvar_vector S=exp(-Z);
  dvar_vector NZ=elem_div(N,Z);
  for (i=1;i<=nfi;i++)
  {
    for (k=1;k<=nfi;k++)
    {
      for (j=jmin;j<=nage;j++)
      {
        dvariable tmp=0;
        //tmp-=F(i,j)/Z(j)*sel(k,j);
        if (i==k)
        {
          tmp+=sel(k,j);
        }
        tmp*=(1.0-S(j));
        if (Z(j)<=rmax)
        {
          //tmp+=F(i,j)*S(j)*sel(k,j);  //baranov
        }
        else
        {
          //tmp+=F(i,j)*exp(-rmax)*sel(k,j);  //kludged
        }
        //J(i,k)+=tmp*NZ(j);
      }
    }
  }
  return J;
}

dvar_matrix get_jacx(int nfi,int jmin,int nage,dvar_vector q,
  dvar_matrix sel,const dvar_vector& N,dvar_vector& M,MY_DOUBLE_TYPE rmax)
{
  int i,j,k;
  dvar_matrix J(1,nfi,1,nfi);
  J.initialize();
  dvar_matrix F(1,nfi);
  dvar_vector Z(jmin,nage);
  return J;
}
*/


dvariable dvar_len_fish_stock_history::tag_catch_equations_calc(dvar_vector& sv)
{
  maxZ[1]=0.0;
  dvariable ffpen=0.0;
  MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
  int break_flag=0;
  ivector tmp_mp(1,num_regions);
  int mmin,mmax;
  int ir,ip,it;
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

  if (simulation_seeds(5))
  {
    cout << "XXX" << endl;
    ad_exit(1);
  }
  if (pmsd) rescale_initial_tag_population();  //NMD 11Nov13

  for (it=mmin;it<=mmax;it++)
  {

    int ng=get_nage_tag(it);
    get_initial_tag_population(sv,it);

    int current_year=tag_year(it);
    //ivector rip(1,num_regions);
    ivector rip(tag_region_bounds(1,it),tag_region_bounds(2,it));
    //ivector skip_flag(1,num_regions);
    ivector skip_flag(tag_region_bounds(1,it),tag_region_bounds(2,it));
    rip=initial_tag_period(it);
    skip_flag.initialize();
    ivector tmp_mp(tag_region_bounds(1,it),tag_region_bounds(2,it));
    ivector bug_tmp_mp(tag_region_bounds(1,it),tag_region_bounds(2,it));
    ivector tmp_ip(tag_region_bounds(1,it),tag_region_bounds(2,it));
    ivector tmp_yr(tag_region_bounds(1,it),tag_region_bounds(2,it));
    ivector tmp_mn(tag_region_bounds(1,it),tag_region_bounds(2,it));
    //ivector tmp_mp(1,num_regions);
    //ivector bug_tmp_mp(1,num_regions);
    //ivector tmp_yr(1,num_regions);
    //ivector tmp_mn(1,num_regions);
    //ivector tmp_ip(1,num_regions);
    tmp_mp.initialize();
    tmp_yr.initialize();
    tmp_mn.initialize();
    tmp_ip.initialize();
    bug_tmp_mp.initialize();

    //for (ir=1;ir<=num_regions;ir++)
    for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
    {
      //if (rip(ir)==num_fish_periods(ir))
      {
        if (move_flags(ir,rip(ir)))
        {
          skip_flag(ir)=1;
          tmp_mp(ir)=move_index(ir,rip(ir));
        }
      }
    }

    int finished_flag=0;
    //for (ir=1;ir<=num_regions;ir++)
    for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
    {
      tagnum_fish(it,ir,rip(ir))=tagN(it,ir,year(ir,rip(ir)));
    }
    do
    {
      ofstream xpofs("tag3__movestuff");
      finished_flag=1;
      int myprintflag=0;
      //for (int ir=1;ir<=num_regions;ir++)
      for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
      {
        break_flag=0;
        int& ip=rip(ir);
        do
        {
          if (ip>=num_fish_periods(ir)) break;
          finished_flag=0;
          if (skip_flag(ir))
          {
            skip_flag(ir)=0;
            break;
          }
          int jmin=tagnum_fish(it,ir,ip).indexmin();
          if (year(ir,ip+1)==year(ir,ip))
          {
            if (!num_fish_incidents(ir,ip) || !tag_flags(it,1) || 
              ip >= initial_tag_period(it,ir)+tag_flags(it,1))
            {
//              tagnum_fish(it,ir,ip+1)=tagnum_fish(it,ir,ip)
//                -tot_mort(ir,ip)(jmin,ng);
              if (parest_flags(360)==0)  //NMD_11Feb2022
              {
                tagnum_fish(it,ir,ip+1)=tagnum_fish(it,ir,ip)
                  -tot_mort(ir,ip)(jmin,ng);
              }
              else
              {
                tagnum_fish(it,ir,ip+1)=tagnum_fish(it,ir,ip)
                  -tot_mort(ir,ip)(jmin,ng)-tagmort(it);
              }

            }
            else
            {
              //new_do_newton_raphson_for_tags(it,ir,ip,ffpen);
              do_newton_raphson_for_tags2(it,ir,ip,ffpen);
              // new code to do NR properly DF May 08 08
              //new_do_newton_raphson_for_tags(it,ir,ip,ffpen);
//              tagnum_fish(it,ir,ip+1)=tagnum_fish(it,ir,ip)
//                -nrtm(it,ir,ip)(jmin,ng);
              if (parest_flags(360)==0)   //NMD_11feb2022
              {
                tagnum_fish(it,ir,ip+1)=tagnum_fish(it,ir,ip)
                  -nrtm(it,ir,ip)(jmin,ng);
              }
              else
              {
                tagnum_fish(it,ir,ip+1)=tagnum_fish(it,ir,ip)
                  -nrtm(it,ir,ip)(jmin,ng)-tagmort(it);
              }
            }

              if (myprintflag==1)
              {
		ofsdd << endl << "  it: " << it << "  ir: " << ir << "  ip: " << ip << endl;  //NMD28Sep2012
                ofsdd << endl << tagnum_fish(it,ir,ip) << endl;
                ofsdd << tagnum_fish(it,ir,ip+1) << endl;
                ofsdd << nrtm(it,ir,ip)(jmin,ng) << endl;
                ofsdd << tot_mort(ir,ip,ng-1) << endl;
                ofsdd << ffpen << endl;
                myprintflag=0;
              }


          }
          else
          {
            // age the tags
            //tagnum_fish(it,ir,ip+1,1)=-15.0;

            if (!num_fish_incidents(ir,ip) || !tag_flags(it,1) || 
              ip >= initial_tag_period(it,ir)+tag_flags(it,1))
            {
              if (ng>2 && jmin<ng-1)
              {
//                --tagnum_fish(it,ir,ip+1)(jmin+1,ng-1)=
//                  tagnum_fish(it,ir,ip)(jmin,ng-2)
//                    -tot_mort(ir,ip)(jmin,ng-2);
                if (parest_flags(360)==0)
                {
                  --tagnum_fish(it,ir,ip+1)(jmin+1,ng-1)=
                    tagnum_fish(it,ir,ip)(jmin,ng-2)
                      -tot_mort(ir,ip)(jmin,ng-2);
                }
                else
                {
                  --tagnum_fish(it,ir,ip+1)(jmin+1,ng-1)=
                    tagnum_fish(it,ir,ip)(jmin,ng-2)
                      -tot_mort(ir,ip)(jmin,ng-2)-tagmort(it);
                }
              }
              if (jmin<ng)
              {
//                tagnum_fish(it,ir,ip+1,ng)=
//                  log(1.e-12 + mfexp(tagnum_fish(it,ir,ip,ng-1)
//                  -tot_mort(ir,ip,ng-1))
//                  + mfexp(tagnum_fish(it,ir,ip,ng)-tot_mort(ir,ip,ng)) );
                if (parest_flags(360)==0)  //NMD_11Feb2022
                {
                  tagnum_fish(it,ir,ip+1,ng)=
                    log(1.e-12 + mfexp(tagnum_fish(it,ir,ip,ng-1)
                    -tot_mort(ir,ip,ng-1))
                    + mfexp(tagnum_fish(it,ir,ip,ng)-tot_mort(ir,ip,ng)) );
                }
                else
                {
                  tagnum_fish(it,ir,ip+1,ng)=
                    log(1.e-12 + mfexp(tagnum_fish(it,ir,ip,ng-1)
                    -tot_mort(ir,ip,ng-1)-tagmort(it))
                    + mfexp(tagnum_fish(it,ir,ip,ng)-tot_mort(ir,ip,ng)
                   -tagmort(it)) );
                }
              }
              else
              {
//                tagnum_fish(it,ir,ip+1,ng)=
//                  tagnum_fish(it,ir,ip,ng)-tot_mort(ir,ip,ng);
                if (parest_flags(360)==0)
                {
                  tagnum_fish(it,ir,ip+1,ng)=
                    tagnum_fish(it,ir,ip,ng)-tot_mort(ir,ip,ng);
                }
                else
                {
                  tagnum_fish(it,ir,ip+1,ng)=
                    tagnum_fish(it,ir,ip,ng)-tot_mort(ir,ip,ng)-tagmort(it);
                }
              }
            }
            else
            {
              //new_do_newton_raphson_for_tags(it,ir,ip,ffpen);
              do_newton_raphson_for_tags2(it,ir,ip,ffpen);
              // new code to do NR properly DF May 08 08
              //new_do_newton_raphson_for_tags(it,ir,ip,ffpen);



              int n1=ng-1;
              if (ng>2) 
              {
                if (jmin<n1)
                {
//                  --tagnum_fish(it,ir,ip+1)(jmin+1,ng-1)=
//                  tagnum_fish(it,ir,ip)(jmin,ng-2)-nrtm(it,ir,ip)(jmin,ng-2);
                  if (parest_flags(360)==0)   //NMD_11feb2022
                  {
                    --tagnum_fish(it,ir,ip+1)(jmin+1,ng-1)=
                    tagnum_fish(it,ir,ip)(jmin,ng-2)-nrtm(it,ir,ip)(jmin,ng-2);
                  }
                  else
                  {
                    --tagnum_fish(it,ir,ip+1)(jmin+1,ng-1)=
                    tagnum_fish(it,ir,ip)(jmin,ng-2)-nrtm(it,ir,ip)(jmin,ng-2)
                    -tagmort(it);
                  }
                }
              }
              if (jmin<ng)
              {
//                tagnum_fish(it,ir,ip+1,ng)=
//                  log(1.e-10 + mfexp(tagnum_fish(it,ir,ip,ng-1)
//                  -tot_mort(ir,ip,ng-1))
//                  + mfexp(tagnum_fish(it,ir,ip,ng)-nrtm(it,ir,ip,ng)));
                if (parest_flags(360)==0)   //NMD_11feb2022
                {
                  tagnum_fish(it,ir,ip+1,ng)=
                    log(1.e-10 + mfexp(tagnum_fish(it,ir,ip,ng-1)
                    -nrtm(it,ir,ip,ng-1))
                    + mfexp(tagnum_fish(it,ir,ip,ng)-nrtm(it,ir,ip,ng)));
                }
                else
                {
                  tagnum_fish(it,ir,ip+1,ng)=
                    log(1.e-10 + mfexp(tagnum_fish(it,ir,ip,ng-1)
                    -nrtm(it,ir,ip,ng-1)-tagmort(it))
                    + mfexp(tagnum_fish(it,ir,ip,ng)-nrtm(it,ir,ip,ng)
                   -tagmort(it)));
                }
              }
              else
              {
//                tagnum_fish(it,ir,ip+1,ng)=
//                  tagnum_fish(it,ir,ip,ng)-nrtm(it,ir,ip,ng);
                if (parest_flags(360)==0)   //NMD_11feb2022
                {
                  tagnum_fish(it,ir,ip+1,ng)=
                    tagnum_fish(it,ir,ip,ng)-nrtm(it,ir,ip,ng);
                }
                else
                {
                  tagnum_fish(it,ir,ip+1,ng)=
                    tagnum_fish(it,ir,ip,ng)-nrtm(it,ir,ip,ng)-tagmort(it);
                }
              }
              if (myprintflag==1)
              {
		ofsdd << endl << "  it: " << it << "  ir: " << ir << "  ip: " << ip << endl;  //NMD28Sep2012
                ofsdd << endl << tagnum_fish(it,ir,ip) << endl;
                ofsdd << tagnum_fish(it,ir,ip+1) << endl;
                ofsdd << nrtm(it,ir,ip)(jmin,ng) << endl;
                ofsdd << tot_mort(ir,ip,ng-1) << endl;
                ofsdd << ffpen << endl;
                myprintflag=0;

              }

            }
          }
          if (move_flags(ir,ip+1))
          {
            tmp_ip(ir)=ip+1;
            tmp_yr(ir)=year(ir,ip+1);
            tmp_mn(ir)=month(ir,ip+1);
            tmp_mp(ir)=move_index(ir,ip+1);
            bug_tmp_mp(ir)=move_index(ir,ip);

            if (ir==num_regions)
            {
              xpofs << "tag release group " << it << "  ";
              print_movement_stuff(xpofs,tmp_mp,tmp_yr,tmp_mn,tmp_ip,
                num_regions);
            }

            break_flag=1;
          } 
          ip++;
        }
        while (!break_flag); 
      }
      // move the tags
      if (sum(tmp_mp))
      {
        if (num_regions>1)
        {
          check_sanity(tmp_mp,tagnum_fish(it),rip,it,*this);
          if (pmsd) pmsd->tag_index=it;
          dvar_matrix tmp=fast_diffusion_calcs(ng,num_regions,
            tagnum_fish(it),Dad(tmp_mp(tmp_mp.indexmin())),rip,0,pmsd);
            //tagnum_fish(it),Dad(tmp_mp(1)),rip);
          if (pmsd) pmsd->tag_index=0;
          //for (int ir=1;ir<=num_regions;ir++)
          for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
          {
            tagnum_fish(it,ir,rip(ir))=log(1.e-12+tmp(ir));
          }
        }
      }
      tmp_mp=0;
    } // need to decide when to quit
    while (!finished_flag);
  }

  
  for (it=mmin;it<=mmax;it++)
  {
    int rmin=1;
    int rmax=num_regions;
    int ng=get_nage_tag(it);
 
    for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
    {
      for (int ip=initial_tag_period(it,ir);ip<=num_fish_periods(ir);ip++)  
      {
        if (!tag_flags(it,1) || 
          ip >= initial_tag_period(it,ir)+tag_flags(it,1))
        {
          for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
          {        
            dvar_vector& tc=tagcatch(it,ir,ip,fi);
            int jmin=tc.indexmin();
//            tc=mfexp(fish_mort_calcs(ir,ip,fi)(jmin,ng)
//              +tagnum_fish(it,ir,ip));
            if (parest_flags(360)==0)  //NMD_11Feb2022
            {
              tc=mfexp(fish_mort_calcs(ir,ip,fi)(jmin,ng)
                +tagnum_fish(it,ir,ip));
            }
            else
            {
              tc=mfexp(tag_fish_mort_calcs(1,ir,ip,fi)(jmin,ng)
                +tagnum_fish(it,ir,ip));
            }
          }
        }
        else
        {
          //double ssum=0;
          for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
          {        
//            tagcatch(it,ir,ip,fi)=
//              mfexp(nrfm(it,ir,ip,fi)-log(1.e-10+nrtm(it,ir,ip))+
//              log(one_plus-nrsurv(it,ir,ip))+tagnum_fish(it,ir,ip));
            dvar_vector tmp1;  //NMD_11Feb2022
            dvar_vector tmp2;
            if (parest_flags(360)==0)
            {
              tmp1=log(1.e-10+nrtm(it,ir,ip));
              tmp2=log(one_plus-nrsurv(it,ir,ip));
            }
            else
            {
              tmp1=log(1.e-10+nrtm(it,ir,ip)) + tagmort(it);
              tmp2=log(one_plus-nrsurv(it,ir,ip)) - tagmort(it);
            }		
            tagcatch(it,ir,ip,fi)=
              mfexp(nrfm(it,ir,ip,fi) - tmp1 +
              tmp2 + tagnum_fish(it,ir,ip));
          }
        }
      }
    }
   
   /*
    if (num_tag_releases>50 && it==67)
    {
      ofstream ofs("ntr67");
      for (int ir=7;ir<=12;ir++)
      {
        ofs <<"region " << ir-6 << endl;
        ofs << tagcatch(it,ir) << endl;
      }
      ad_exit(1);
    }
    if (num_tag_releases<50 && it==18)
    {
      ofstream ofs("ntr18");
      for (int ir=1;ir<=6;ir++)
      {
        ofs <<"region " << ir << endl;
        ofs << tagcatch(it,ir) << endl;
      }
      ad_exit(1);
    }
   */
  }
  cout << " ffpen = " << ffpen << endl;
  return ffpen;
}


void check_for_early_returns(const i3_array& tr,
  const ivector& true_tag_year,const ivector& true_tag_month,int ntr,
  ivector& itr)
{
  int bad_flag=0;
  for (int it=1;it<=ntr;it++)
  {
    if (itr(it)>0)
    {
      int mmin=tr(it).indexmin();
      int mmax=tr(it).indexmax();
      int mmb=12*true_tag_year(it)+true_tag_month(it);
      for (int i=mmin;i<=mmax;i++)
      {
        int yr=tr(it,i,3);
        int mn=tr(it,i,4);
        int mm=12*yr+mn;
        if (mm<mmb)
       {
          cout << "Tag returns record " << i << " for tag group " << it
               << " were caught before they were released " << endl;
          bad_flag=1;
        }
      }
    }
  }
  if (bad_flag) 
  {
    cout << " That does it -- I quit!!" << endl;
    exit(1);
  }  
}  

void dvar_fish_stock_history::rescale_initial_tag_population(void)
{
  for (int itt=1;itt<=num_tag_releases;itt++)
  {
    int ngt=pmsd->tag_species_index(itt).indexmax();
    dvar_vector weights(1,ngt);
    if (ngt>1)
    {
      if (itt==pmsd->tag_species_index(itt,1))
      {
        for (int ic=1;ic<=ngt;ic++)
        {
          int it=pmsd->tag_species_index(itt,ic);
          int ir=tag_region(it);
          int ip=initial_tag_period(it,ir);
          dvar_vector en=exp(num_fish(ir,ip));
          weights(ic)=sum(en);
        }
        weights/=sum(weights);
        weights/=ngt;
        for (int ic=1;ic<=ngt;ic++)
        {
          int it=pmsd->tag_species_index(itt,ic);
          initial_tag_release_by_age(it)*=weights(ic);
        }
      }
    }
  }
}  

void dvar_fish_stock_history::get_initial_tag_population(dvar_vector& sv, int it)
{
  tagN(it).initialize();
  int ir;

  int yr1=tag_year(it);
  //for (ir=1;ir<=num_regions;ir++)
  for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
  {
    tagN(it,ir)=-15.0;
  }
  ir=tag_region(it);
  tagN(it,ir,yr1)=log(initial_tag_release_by_age(it)+1.e-12);

  // here is where we "diffuse " the tagged fish to move them between 
  // the regions
  if (pmsd) pmsd->tag_index=it;
  if (num_regions>1) do_the_diffusion(yr1,sv,tagN(it));
  if (pmsd) pmsd->tag_index=0;
}


void dvar_len_fish_stock_history::var_convert_tag_lengths_to_age(void)
{
  initial_tag_release_by_age.initialize();
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
  if (!pmsd)
  {
    for (int it=mmin;it<=mmax;it++)
    {
      dvector& itrl=initial_tag_release_by_length(it);
      for (int il=1;il<=tag_nlint;il++)
      {
        if (itrl(il))
        {        
          MY_DOUBLE_TYPE len=tag_shlen+(il-0.5)*tag_filen;
          dvariable age;
          if (!parest_flags(175) && !parest_flags(174))
          {
            age=age_at_length_calc(len,vb_coff,nage,parest_flags);
          }
          else
            age=age_at_length_calcxx(len,gml,nage+1);
          MY_DOUBLE_TYPE cage=value(age);
          if (cage<=1)
          {
            initial_tag_release_by_age(it,1)+=itrl(il);
          }
          else if (cage>=nage)
          {
            initial_tag_release_by_age(it,nage)+=itrl(il);
          }
          else
          {
            dvariable sf;
            sf=daves_kludge1(age);
            int jj=int(cage);
            dvariable tp= sf*itrl(il);
  
            initial_tag_release_by_age(it,jj)+=
              itrl(il);
            initial_tag_release_by_age(it,jj)-=tp;
  
            initial_tag_release_by_age(it,jj+1)+=tp;
          }
        } 
      }
    }
  }
  else
  {
    for (int it=mmin;it<=mmax;it++)
    {
      int isp=pmsd->tag_species_pointer(it);
      //int offset=4*(isp-2);
        
      dvector& itrl=initial_tag_release_by_length(it);
      for (int il=1;il<=tag_nlint;il++)
      {
        if (itrl(il))
        {        
          MY_DOUBLE_TYPE len=tag_shlen+(il-0.5)*tag_filen;
          dvariable age;
          if (!parest_flags(175) && !parest_flags(174))
          {
            if (isp==1)
            {
              age=age_at_length_calc(len,vb_coff,nage,parest_flags);
            }
            else
            {
              dvar_vector tvb=pmsd->vb_coff(isp);
//              age=age_at_length_calc(len,tvb,nage,parest_flags);
              age=age_at_length_calc(len,tvb,pmsd->nage(isp),parest_flags);   //NMD_27Sep2018
            }
          }
          else
            age=age_at_length_calcxx(len,gml,nage+1);
          MY_DOUBLE_TYPE cage=value(age);
          if (cage<=1)
          {
            initial_tag_release_by_age(it,1)+=itrl(il);
          }
          else if (cage>=nage)
          {
            initial_tag_release_by_age(it,nage)+=itrl(il);
          }
          else
          {
            dvariable sf;
            sf=daves_kludge1(age);
            int jj=int(cage);
            dvariable tp= sf*itrl(il);
  
            initial_tag_release_by_age(it,jj)+=
              itrl(il);
            initial_tag_release_by_age(it,jj)-=tp;
  
            initial_tag_release_by_age(it,jj+1)+=tp;
          }
        } 
      }
    }
  }
}

void dvar_len_fish_stock_history::tot_tags_catch(void)
{
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
    //for (int ir=1;ir<=num_regions;ir++)
    for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
    {
      int ipmax=num_fish_periods(ir);  
      if (age_flags(96)) ipmax=terminal_tag_period(it,ir);
      for (int ip=initial_tag_period(it,ir);ip<=ipmax;ip++)
      {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {
          tot_tag_catch(it,ir,ip,fi)= sum(obstagcatch_by_length(it,ir,ip,fi));
        }
      }
    }
  }
  if (age_flags(96))
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=minttp(ir)+1;ip<=num_fish_periods(ir);ip++)  
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {
        pooledtot_tag_catch(ir,ip,fi)
          = sum(pooledobstagcatch_by_length(ir,ip,fi));
      }
    }
  }
}

void dvar_fish_stock_history::tag_catch_equations_calc_mc(dvar_vector& sv)
{
  MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;

  for (int it=1;it<=num_tag_releases;it++)
  {
    get_initial_tag_population(sv,it);

    int current_recruitment_period=initial_tag_recruitment_period(it,1);
    ivector rip(1,num_regions);
    rip=initial_tag_period(it);
    do
    {
      for (int ir=1;ir<=num_regions;ir++)
      {
        int& ip=rip(ir);
        if (recruitment_period(ir,ip)==current_recruitment_period)
        {
          tagnum_fish(it,ir,ip)=tagN(it,ir,current_recruitment_period);
          do
          {
            if (ip>=num_fish_periods(ir)) break;
            if (recruitment_period(ir,ip+1)==current_recruitment_period)
            {
              tagnum_fish(it,ir,ip+1)=tagnum_fish(it,ir,ip)-tot_mort(ir,ip);
              ip++;
            }
            else
            {
              tagnum_fish(it,ir,ip+1,1)=-15.0;

              --tagnum_fish(it,ir,ip+1)(2,nage)=
                tagnum_fish(it,ir,ip)(1,nage-1)-tot_mort(ir,ip)(1,nage-1);
 
              tagnum_fish(it,ir,ip+1,nage)=
                log(1.e-10 + mfexp(tagnum_fish(it,ir,ip,nage-1)-tot_mort(ir,ip,nage-1))
                  + mfexp(tagnum_fish(it,ir,ip,nage)-tot_mort(ir,ip,nage)) );
  
              if (current_recruitment_period <num_recruitment_periods)
              {
                tagN(it,ir,current_recruitment_period+1)=tagnum_fish(it,ir,ip+1);
              }
              ip++;
              break;
            }
          }
          while (1);
        }
        else   // there were no fisheries in this region for this year
        {
          if (current_recruitment_period <num_recruitment_periods)
          {
            if (!pmsd)
            {
              for (int j=1;j<nage;j++)      // Loop over age classes
              {
                tagN(it,ir,current_recruitment_period+1,j+1)=
                  tagN(it,ir,current_recruitment_period,j)-
                    exp(nat_mort(current_recruitment_period,j));
              }
            }
            else
            {
              const dvar_matrix& nm=get_nat_mort_region(ir);
              for (int j=1;j<nage;j++)      // Loop over age classes
              {
                tagN(it,ir,current_recruitment_period+1,j+1)=
                  tagN(it,ir,current_recruitment_period,j)-
                    exp(nm(current_recruitment_period,j));
              }
            }
          }
        }
      }
      if (current_recruitment_period <num_recruitment_periods)
      {
      // Changed af(57) to af(53) J.H. 27 Nov 01
        if (age_flags(53))
        {
          if ( !((current_recruitment_period)%age_flags(53)) )
          {
            if (num_regions>1)
              do_the_diffusion(current_recruitment_period+1,sv,tagN(it));
          }
        }
        else
        {
          if (num_regions>1)
            do_the_diffusion(current_recruitment_period+1,sv,tagN(it));
        }
      }
      current_recruitment_period++;
    }
    while (current_recruitment_period <=num_recruitment_periods);
  
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int ip=initial_tag_period(it,ir);ip<=num_fish_periods(ir);ip++)  
      {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {        
          tagcatch(it,ir,ip,fi)=exp(fish_mort(ir,ip,fi)-log(1.e-10+tot_mort(ir,ip))+
            log(one_plus-survival(ir,ip))+tagnum_fish(it,ir,ip));
        }
      }
    }
  }
}


void dvar_fish_stock_history::print_tagging_fit_info(ofstream& of)
{
  ivector group_flags=column(fish_flags,32);
  ivector ff37=column(fish_flags,37);
  int i,ir,ff,ip;
  of << "# Tag reporting rates" << endl;
  of << "# Grouping indicator (0 = no grouping, >0 = grouping)" << endl;
  of << group_flags << endl;
  of << "# Time series variation in reporting rates (0 = no, >0 = yes)" << endl;
  of << sum(ff37) << endl;
  
  if (!sum(ff37))
  {
    if (!age_flags(198))
    {
      of << "# Reporting rates by fishery (no time series variation)" << endl;
      for (int ifish=1;ifish<=num_fisheries;ifish++)
      {
        of << fish_pars(3,ifish) << endl;
      }
    }
    else
    {
        of << "# tag fish rep group flags" << endl;       //NMD 29Mar2012
//NMD 29Mar2012
		for (int ifish=1;ifish<=num_fisheries;ifish++)
		{
			for(int irgrp=tag_fish_rep_group_flags.indexmin(); irgrp<=tag_fish_rep_group_flags.indexmax(); irgrp++)
			{
				of << tag_fish_rep_group_flags(irgrp,ifish) << " " ;				
			}
			of << endl;
		}
//        of <<  tag_fish_rep_group_flags << endl;          //NMD 29Mar2012
        of << "# Reporting rates by fishery by tag group"
           " (no time series variation)" << endl;
		//NMD 29Mar2012
      for (int ifish=1;ifish<=num_fisheries;ifish++)
      {
//        of << column(tag_fish_rep,ifish) << endl;
        of << setfixed() << setprecision(7) << setw(12) << column(tag_fish_rep,ifish) << endl;    //NMD7Nov2011
      }
    }
  }
  else
  {
    of << "# Reporting rates by fishery by time period (across)" << endl;
    for (int i=1;i<=num_fisheries;i++)
    {
      for (int nt=itind(i);nt<=num_fish_times(i);nt++)
      {
        int rr=realization_region(i,nt);
        int rp=realization_period(i,nt);
        of <<  setprecision(3) 
           <<  rep_rate(rr,rp,realization_incident(i,nt)) << " ";
      }
      of << endl;
    }
  }
  int tmult=1;
  if (age_flags(57)) tmult=age_flags(57);
  int min_tag_year=nyears; 
  int mtt_reg=1;
  for (ir=1;ir<=num_regions;ir++)
  {
    int cur_tag_year=(true_year(ir,minimum_initial_tag_period(ir))-1)*tmult+1;
    if (min_tag_year > cur_tag_year)
    {
      min_tag_year=cur_tag_year;
      mtt_reg=ir;
    }
  }
  int ip_min=min_tag_year;
  int ip_max=nyears;
  if(!sum(group_flags))
  {
    dmatrix ob(1,num_fisheries,ip_min,ip_max);
    dmatrix pr(1,num_fisheries,ip_min,ip_max);
    ob.initialize();
    pr.initialize();
    for (ir=1;ir<=num_regions;ir++)
    {
      for (int ip=minimum_initial_tag_period(ir);ip<=num_fish_periods(ir);ip++) 
      {
        int iper=(true_year(ir,ip)-1)*tmult+(true_month(ir,ip)-1)/(12/tmult)+1;
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
        {                                         // incidents for this period
          int ff=parent(ir,ip,fi);
          for (int it=1;it<=num_tag_releases;it++)
          {
            if (age_flags(96))
            {    
              if (initial_tag_period(it,ir)<=ip
                && terminal_tag_period(it,ir)>=ip)
              {
                ob(ff,iper)
                  +=sum(value(obstagcatch(it,ir,ip,fi)));
                pr(ff,iper)
                  +=sum(get_rep_rate_correction(it,ir,ip,fi)*value(tagcatch(it,ir,ip,fi)));
              }
            }
            else
            {
              if (initial_tag_period(it,ir)<=ip)
              {
                ob(ff,iper)
                  +=sum(value(obstagcatch(it,ir,ip,fi)));
                pr(ff,iper)
                  +=sum(get_rep_rate_correction(it,ir,ip,fi)*value(tagcatch(it,ir,ip,fi)));
              }
            }
          }
        }
      }
    }
    of << "# No. of time periods associated with tag returns" << endl;
    of << setfixed() << setprecision(0) << setw(5)
       << ip_max-ip_min+1 << endl;

   of << "# Time periods associated with tag returns" << endl;
    for (int iper=ip_min;iper<=ip_max;iper++)
    {
      MY_DOUBLE_TYPE ctime=double((iper-1))/double(tmult)+0.5/double(tmult)+year1;
      of << setfixed() << setprecision(3) << setw(10) << ctime; 
    }
    of << endl;
    of << "# Observed tag returns by time period (across) by fishery (down)" << endl;
    for (ff=1;ff<=num_fisheries;ff++)
    {
      for (int iper=ip_min;iper<=ip_max;iper++)
      {
        of << setfixed() << setprecision(0) << setw(10) << ob(ff,iper); 
      }
      of << endl;
    }
    of << "# Predicted tag returns by time period (across) by fishery (down)" << endl;
    for (ff=1;ff<=num_fisheries;ff++)
    {
      for (int iper=ip_min;iper<=ip_max;iper++)
      {
        of << setfixed() << setprecision(3) << setw(10) << pr(ff,iper); 
      }
      of << endl;
    }
  }
  else
  {
    int maxg=max(group_flags);
    int ming=min(group_flags);
    dmatrix ob(ming,maxg,ip_min,ip_max);
    ivector reg_group(ming,maxg);
    dmatrix pr(ming,maxg,ip_min,ip_max);
    ob.initialize();
    pr.initialize();
    for (int ir=1;ir<=num_regions;ir++)
    {
      reg_group.initialize();
      for (ip=minimum_initial_tag_period(ir);ip<=num_fish_periods(ir);ip++) 
      {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
        {                                                 // incidents for this period
          int ff=parent(ir,ip,fi);
          reg_group(group_flags(ff))=1;
        }
      }
      for (ip=minimum_initial_tag_period(ir);ip<=num_fish_periods(ir);ip++) 
      {
        int iper=(true_year(ir,ip)-1)*tmult+(true_month(ir,ip)-1)/(12/tmult)+1;
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
        {                                                 // incidents for this period
          int ff=parent(ir,ip,fi);
          for (int it=1;it<=num_tag_releases;it++)
          {
//            if (ir<tag_region_bounds(1,it) || ir >tag_region_bounds(2,it)) 
//              break;
            if (ir >= tag_region_bounds(1,it) && ir <= tag_region_bounds(2,it))
            {       //NMD_13Sep2018
              if (age_flags(96))
              {    
                if (initial_tag_period(it,ir)<=ip
                  && terminal_tag_period(it,ir)>=ip)
                {
                  ob(group_flags(ff),iper)
                    +=sum(value(obstagcatch(it,ir,ip,fi)));
                  pr(group_flags(ff),iper)
                    +=sum(get_rep_rate_correction(it,ir,ip,fi)*value(tagcatch(it,ir,ip,fi)));
                }
              }
              else
              {
                if (initial_tag_period(it,ir)<=ip)
                {
                  ob(group_flags(ff),iper)
                    +=sum(value(obstagcatch(it,ir,ip,fi)));
                  pr(group_flags(ff),iper)
                    +=sum(get_rep_rate_correction(it,ir,ip,fi)*value(tagcatch(it,ir,ip,fi)));
                }
              }
            }   //NMD_13Sep2018
          }
        }
      }
    }
    of << "# No. of time periods associated with tag returns" << endl;
    of << setfixed() << setprecision(0) << setw(5)
       << ip_max-ip_min+1 << endl;

   of << "# Time periods associated with grouped tag returns" << endl;
    for (int iper=ip_min;iper<=ip_max;iper++)
    {
      MY_DOUBLE_TYPE ctime=double((iper-1))/double(tmult)+0.5/double(tmult)+year1;
      of << setfixed() << setprecision(3) << setw(10) << ctime; 
    }
    of << endl;
    of << "# Observed tag returns by time period (across) by fishery groupings (down)" << endl;
    int ig;
    for (ig=ming;ig<=maxg;ig++)
    {
      for (int iper=ip_min;iper<=ip_max;iper++)
      {
        of << setfixed() << setprecision(0) << setw(10) << ob(ig,iper); 
      }
      of << endl;
    }
    of << "# Predicted tag returns by time period (across) by fishery groupings (down)" << endl;
    for (ig=ming;ig<=maxg;ig++)
    {
      for (int iper=ip_min;iper<=ip_max;iper++)
      {
        of << setfixed() << setprecision(3) << setw(10) << pr(ig,iper); 
      }
      of << endl;
    }
  }
}

void dvar_fish_stock_history::new_print_tagging_fit_info(ofstream& of)
{
  ivector group_flags=column(fish_flags,32);
  ivector ff37=column(fish_flags,37);
  int i,ir,ff,ip;
  of << "# Tag reporting rates" << endl;
  of << "# Grouping indicator (0 = no grouping, >0 = grouping)" << endl;
  of << group_flags << endl;
  of << "# Time series variation in reporting rates (0 = no, >0 = yes)" << endl;
  of << sum(ff37) << endl;
  
  if (!sum(ff37))
  {
    of << "# tag fish rep group flags" << endl;       //NMD 29Mar2012
    for (int ifish=1;ifish<=num_fisheries;ifish++)
    {
      for(int irgrp=tag_fish_rep_group_flags.indexmin(); irgrp<=tag_fish_rep_group_flags.indexmax(); irgrp++)
      {
        of << tag_fish_rep_group_flags(irgrp,ifish) << " " ;
      }
      of << endl;
    }
    of << "# Reporting rates by fishery by tag group"
           " (no time series variation)" << endl;
    for (int ifish=1;ifish<=num_fisheries;ifish++)
    {
      of << setfixed() << setprecision(7) << setw(12) << column(tag_fish_rep,ifish) << endl;    //NMD7Nov2011
    }
  }
  else
  {
    of << "# Reporting rates by fishery by time period (across)" << endl;
    for (int i=1;i<=num_fisheries;i++)
    {
      for (int nt=itind(i);nt<=num_fish_times(i);nt++)
      {
        int rr=realization_region(i,nt);
        int rp=realization_period(i,nt);
        of <<  setprecision(3) 
           <<  rep_rate(rr,rp,realization_incident(i,nt)) << " ";
      }
      of << endl;
    }
  }

  int tmult=1;
  if (age_flags(57)) tmult=age_flags(57);
  int min_tag_year=nyears; 
  int mtt_reg=1;
  for (ir=1;ir<=num_regions;ir++)
  {
    int cur_tag_year=(true_year(ir,minimum_initial_tag_period(ir))-1)*tmult+1;
    if (min_tag_year > cur_tag_year)
    {
      min_tag_year=cur_tag_year;
      mtt_reg=ir;
    }
  }
  int ip_min=min_tag_year;
  int ip_max=nyears;
  if(!sum(group_flags))
  {
    dmatrix ob(1,num_fisheries,ip_min,ip_max);
    dmatrix pr(1,num_fisheries,ip_min,ip_max);
    ob.initialize();
    pr.initialize();
    for (ir=1;ir<=num_regions;ir++)
    {
      for (int ip=minimum_initial_tag_period(ir);ip<=num_fish_periods(ir);ip++) 
      {
        int iper=(true_year(ir,ip)-1)*tmult+(true_month(ir,ip)-1)/(12/tmult)+1;
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
        {                                         // incidents for this period
          int ff=parent(ir,ip,fi);
          for (int it=1;it<=num_tag_releases;it++)
          {
            if (age_flags(96))
            {    
              if (initial_tag_period(it,ir)<=ip
                && terminal_tag_period(it,ir)>=ip)
              {
                ob(ff,iper)
                  +=sum(value(obstagcatch(it,ir,ip,fi)));
                pr(ff,iper)
                  +=sum(value(tag_rep_rate(it,ir,ip,fi))*value(tagcatch(it,ir,ip,fi)));
              }
            }
            else
            {
              if (initial_tag_period(it,ir)<=ip)
              {
                ob(ff,iper)
                  +=sum(value(obstagcatch(it,ir,ip,fi)));
                pr(ff,iper)
                  +=sum(value(tag_rep_rate(it,ir,ip,fi))*value(tagcatch(it,ir,ip,fi)));
              }
            }
          }
        }
      }
    }
    of << "# No. of time periods associated with tag returns" << endl;
    of << setfixed() << setprecision(0) << setw(5)
       << ip_max-ip_min+1 << endl;

   of << "# Time periods associated with tag returns" << endl;
    for (int iper=ip_min;iper<=ip_max;iper++)
    {
      MY_DOUBLE_TYPE ctime=double((iper-1))/double(tmult)+0.5/double(tmult)+year1;
      of << setfixed() << setprecision(3) << setw(10) << ctime; 
    }
    of << endl;
    of << "# Observed tag returns by time period (across) by fishery (down)" << endl;
    for (ff=1;ff<=num_fisheries;ff++)
    {
      for (int iper=ip_min;iper<=ip_max;iper++)
      {
        of << setfixed() << setprecision(0) << setw(10) << ob(ff,iper); 
      }
      of << endl;
    }
    of << "# Predicted tag returns by time period (across) by fishery (down)" << endl;
    for (ff=1;ff<=num_fisheries;ff++)
    {
      for (int iper=ip_min;iper<=ip_max;iper++)
      {
        of << setfixed() << setprecision(3) << setw(10) << pr(ff,iper); 
      }
      of << endl;
    }
  }
  else
  {
    int maxg=max(group_flags);
    int ming=min(group_flags);
    dmatrix ob(ming,maxg,ip_min,ip_max);
    ivector reg_group(ming,maxg);
    dmatrix pr(ming,maxg,ip_min,ip_max);
    ob.initialize();
    pr.initialize();
    for (int ir=1;ir<=num_regions;ir++)
    {
      reg_group.initialize();
      for (ip=minimum_initial_tag_period(ir);ip<=num_fish_periods(ir);ip++) 
      {
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
        {                                                 // incidents for this period
          int ff=parent(ir,ip,fi);
          reg_group(group_flags(ff))=1;
        }
      }
      for (ip=minimum_initial_tag_period(ir);ip<=num_fish_periods(ir);ip++) 
      {
        int iper=(true_year(ir,ip)-1)*tmult+(true_month(ir,ip)-1)/(12/tmult)+1;
        for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
        {                                                 // incidents for this period
          int ff=parent(ir,ip,fi);
          for (int it=1;it<=num_tag_releases;it++)
          {
            if (age_flags(96))
            {    
              if (initial_tag_period(it,ir)<=ip
                && terminal_tag_period(it,ir)>=ip)
              {
                ob(group_flags(ff),iper)
                  +=sum(value(obstagcatch(it,ir,ip,fi)));
                pr(group_flags(ff),iper)
                  +=sum(value(tag_rep_rate(it,ir,ip,fi))*value(tagcatch(it,ir,ip,fi)));
              }
            }
            else
            {
              if (initial_tag_period(it,ir)<=ip)
              {
                ob(group_flags(ff),iper)
                  +=sum(value(obstagcatch(it,ir,ip,fi)));
                pr(group_flags(ff),iper)
                  +=sum(value(tag_rep_rate(it,ir,ip,fi))*value(tagcatch(it,ir,ip,fi)));
              }
            }
          }
        }
      }
    }
    of << "# No. of time periods associated with tag returns" << endl;
    of << setfixed() << setprecision(0) << setw(5)
       << ip_max-ip_min+1 << endl;

    of << "# Time periods associated with grouped tag returns" << endl;
    for (int iper=ip_min;iper<=ip_max;iper++)
    {
      MY_DOUBLE_TYPE ctime=double((iper-1))/double(tmult)+0.5/double(tmult)+year1;
      of << setfixed() << setprecision(3) << setw(10) << ctime; 
    }
    of << endl;
    of << "# Observed tag returns by time period (across) by fishery groupings (down)" << endl;
    int ig;
    for (ig=ming;ig<=maxg;ig++)
    {
      for (int iper=ip_min;iper<=ip_max;iper++)
      {
        of << setfixed() << setprecision(0) << setw(10) << ob(ig,iper); 
      }
      of << endl;
    }
    of << "# Predicted tag returns by time period (across) by fishery groupings (down)" << endl;
    for (ig=ming;ig<=maxg;ig++)
    {
      for (int iper=ip_min;iper<=ip_max;iper++)
      {
        of << setfixed() << setprecision(3) << setw(10) << pr(ig,iper); 
      }
      of << endl;
    }
  }
}

/*  comment this out whuile we try the ungrouped catch equation dynamics
from initial_tag_period to terminal_tag_period
*/

void dvar_fish_stock_history::get_initial_tag_fishery_realization_index(void)
{
  for (int i=1;i<=num_fisheries;i++)
  {
    itind(i)=num_fish_times(i);
    for (int nt=1;nt<=num_fish_times(i);nt++)
    {
      int rr=realization_region(i,nt);
      int rp=realization_period(i,nt);
      if(rp>=minimum_initial_tag_period(rr))
      {
        itind(i)=nt;
        break;
      }
    }
  }
}


//ofstream ofsff("nr_tags");   //NMD1Oct2012
/*
void dvar_fish_stock_history::
  do_newton_raphson_for_tags(int it,int ir,int ip,dvariable& _ffpen)
{
  //get the number of non zero tag catches
  int nfi=num_fish_incidents(ir,ip);
  int fi;
  int nzfi=0;
  for (fi=1;fi<=nfi;fi++)
  {
    if (tot_tag_catch(it,ir,ip,fi)>1.e-10) nzfi++; 
  }
  ivector isub(1,nzfi);
  int ii=0;
  for (fi=1;fi<=nfi;fi++)
  {
    if (tot_tag_catch(it,ir,ip,fi)>1.e-10) 
    {
      ii++; 
      isub(ii)=fi;
    }
  }
  dvariable ffpen=0.0;
  int jmin=tagnum_fish(it,ir,ip).indexmin();
  const dvar_vector& enf=mfexp(tagnum_fish(it,ir,ip)(isub));
  dvariable tot_num_fish=sum(enf);
  dvar_matrix sel(1,nzfi,jmin,nage);
  dvar_matrix logsel(1,nzfi,jmin,nage);
  int fi;
  for (ii=1;ii<=nzfi;ii++)
  {
    fi=isub(ii);
    int i=parent(ir,ip,fi);
    int rr=realization_region(i,1);
    int rp=realization_period(i,1);
    int ri=realization_incident(i,1);
    logsel(ii)=incident_sel(rr,rp,ri)(jmin,nage);
    sel(ii)=mfexp(logsel(fi));
  }

  dvar_matrix& fm=nrfm(it,ir,ip)(isub);
  dvar_matrix M(1,nfi,1,nzfi);
  M.initialize();
  dvar_matrix C(1,nfi,jmin,nage)(isub);
  dvariable tnf=sum(enf);
  //ofsff << tnf << " " << it << " " << ir << " " << ip << endl;
  if (nfi>1)
  {
    //cout <<"tag3.cpp " << nfi << endl;
  }
  dvar_vector actual_tag_catch=
    elem_div(tot_tag_catch(it,ir,ip)(isub),1.e-6+rep_rate(ir,ip)(isub));
  
  dvariable region_tot_tag_catch=sum(actual_tag_catch);

  const dvar_vector& ctm=nrtm(it,ir,ip)(jmin,nage);
  dvar_vector& tm=(dvar_vector&)ctm;
  if (region_tot_tag_catch>0)
  {
    dvariable surv_rate=(1.0-region_tot_tag_catch/tnf);
  
    dvar_vector kc(1,nfi);


    MY_DOUBLE_TYPE cut=0.2;
    MY_DOUBLE_TYPE fringe=0.02;
    if (age_flags(118)>0)
    {
      cut=age_flags(118)/100.;
    }
    if (age_flags(119)>0)
    {
      fringe=age_flags(119)/100.;
    }
    MY_DOUBLE_TYPE c2=cut+fringe;

    
    if (surv_rate<=c2)
    {
      dvariable tmp=0.0;
      dvariable ks=fringe+posfun(surv_rate-fringe,cut,tmp);

      MY_DOUBLE_TYPE penwt=1.0;
      if (age_flags(117)>0) penwt=age_flags(117)/1.e+8;
      if (age_flags(117)<0) penwt=0.0;
        ffpen+=penwt*tmp;

      if (value(ks)>1.0)
      {
        cout << "fringe " << fringe << " surv rate  " << surv_rate 
             << " cut " << cut << endl;
        cout << "enter number " << endl;
        int num;
        cin >> num;
        cout << "this can't happen ks = " << ks << endl;
        ad_exit(1);
      }
      dvariable kr= (1.0-ks)*tnf/region_tot_tag_catch;
      kc=kr*actual_tag_catch+1.e-5;
    }
    else
    {
      kc=actual_tag_catch+1.e-5;
    }
    dvar_vector qq=elem_div(kc,sel*enf);
    dvar_vector TC(1,nfi);
    int itt=0;
    int badflag=0;
    MY_DOUBLE_TYPE normd=0.0;
    do
    {
      int yr=year(ir,ip);  
      dvar_vector Z=qq*sel+mfexp(nat_mort(yr)(jmin,nage)+fraction(ir,ip));
      dvar_vector S=exp(-Z);
      dvar_vector S1N=elem_prod(1.0-S,enf);
      for (int fi=1;fi<=nfi;fi++)
      {
        dvar_vector t1=elem_div(qq(fi)*sel(fi),Z);
        C(fi)=elem_prod(t1,S1N);
        TC(fi)=sum(C(fi));
        dvariable Dq=sum(C(fi))/qq(fi);
        dvar_vector DZ=qq(fi)*elem_prod(elem_div(sel(fi),Z),enf)-
          elem_prod(C(fi),1.0+1.0/Z);
        M(fi,fi)=Dq;
        for (int fj=1;fj<=nfi;fj++)
        {
          M(fi,fj)+=DZ*sel(fj);
        }
      }
      dvar_vector diff=TC-kc;
      normd=norm(value(diff));
      if (break_flag==1) break; 
      //cout <<"tag3.cpp " <<  elem_div(diff,TC)  << endl;
      int pflag=1;
      if (!badflag)
      {
        // this is newton raphson for q
        qq-=solve(M,diff);
    //     int mmin=qq.indexmin();
    //     int mmax=qq.indexmax();
    //     dvariable fp1=0.0;
    //     for (int i=mmin;i<=mmax;i++)
    //     {
    //       qq(i)=posfun(qq(i),1.e-10,fp1);
    //     }
    //     ffpen+=fp1;
        //pflag=check_pos(qq,nfi);
        if (!pflag)
        {
          dvariable tmp=0.0;
          throw function_minimizer_exception();
          badflag=1;
          qq+=inv(M)*diff;
          qq=elem_div(qq,mfexp(elem_div(inv(M)*diff,qq)));
        }
        else
        {
          itt++;
        }
      }
      else
      {
        // this i newton raphson for log(q)
        qq=elem_div(qq,1.e-20+mfexp(elem_div(inv(M)*diff,qq)));
        //cout << " " << qq(1);
        itt++;
        pflag=check_pos(qq,nfi);
      }
      if ( (itt>5 && pflag && !badflag) || (itt>5 && pflag) )break;
    }
    while(1);
    dvariable fp1=0.0;
    int mmin=qq.indexmin();
    int mmax=qq.indexmax();
    for (int i=mmin;i<=mmax;i++)
    {
      qq(i)=posfun(qq(i),1.e-10,fp1);
    }
    ffpen+=fp1;
    if (normd>.1)
    {
      cout << "difficult tag NR B " << normd;
      cout << "  it= " << it ;
      cout << "  ir= " << ir ;
      cout << "  ip= " << ip << endl;
      cout << "Predicted tag catch           " << setw(10) << setprecision(2) << setfixed() << TC  << endl; 
      cout << "Expanded observed tag catch kc" << setw(10) << setprecision(2) << setfixed() <<  kc << endl;
      cout << "total number of tags present  " << setw(10) << setprecision(2) << setfixed() <<  tnf << endl;
    }
    //cout <<"tag3.cpp " << endl;
    dvar_vector tmpFF=log(qq);
    for (fi=1;fi<=nfi;fi++)
    {
      fm(fi)(jmin,nage)=logsel(fi)+tmpFF(fi);
    }
    tm.initialize();
    for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
    {                                         // incidents for this period
      tm+=mfexp(fm(fi));
    }
    tm+=mfexp(nat_mort(year(ir,ip))(jmin,nage)+fraction(ir,ip));
  }
  else
  {
    for (fi=1;fi<=nfi;fi++)
    {
      fm(fi)=-10.0;
    }
    tm=mfexp(nat_mort(year(ir,ip))(jmin,nage)+fraction(ir,ip));
  }
  nrsurv(it,ir,ip)=mfexp(-tm);
  switch (it)
  {
  case 72:
  case 21:
  case 112:
    break;
  default:
    ;
    _ffpen+=ffpen;
  }
}
*/
void dvar_fish_stock_history::
  do_newton_raphson_for_tags2(int it,int ir,int ip,dvariable& _ffpen)
{
  dvariable ffpen=0.0;
  int jmin=tagnum_fish(it,ir,ip).indexmin();
  const dvar_vector& enf=mfexp(tagnum_fish(it,ir,ip));
  dvariable tot_num_fish=sum(enf);
  int nfi=num_fish_incidents(ir,ip);
  dvar_matrix sel(1,nfi,jmin,nage);
  dvar_matrix logsel(1,nfi,jmin,nage);
  int fi;
  for (fi=1;fi<=nfi;fi++)
  {
    int i=parent(ir,ip,fi);
    int rr=realization_region(i,1);
    int rp=realization_period(i,1);
    int ri=realization_incident(i,1);
    logsel(fi)=incident_sel(rr,rp,ri)(jmin,nage);
    sel(fi)=mfexp(logsel(fi));
    /*
  if(it == 2 && ir == 3){  //NMD2Oct2012
  //if(it == 51 && ir == 9){  //NMD2Oct2012
    ofsff << endl << " it = " << it << "   ir = " << ir << "  ip: " << ip << "  fi: " << fi << endl;
    ofsff << " sel(fi) = " <<  sel(fi) << endl;
    ofsff << " logsel(fi) = " <<  logsel(fi) << endl;
  }  //NMD2Oct2012
    */
  }

  dvar_matrix& fm=nrfm(it,ir,ip);
  dvar_matrix M(1,nfi,1,nfi);
  M.initialize();
  dvar_matrix C(1,nfi,jmin,nage);
  dvariable tnf=sum(enf);
  //ofsff << tnf << " " << it << " " << ir << " " << ip << endl;
  if (nfi>1)
  {
    //cout <<"tag3.cpp " << nfi << endl;
  }
  dvar_vector actual_tag_catch=
    elem_div(tot_tag_catch(it,ir,ip),1.e-6+rep_rate(ir,ip));
  
  dvariable region_tot_tag_catch=sum(actual_tag_catch);

  const dvar_vector& ctm=nrtm(it,ir,ip)(jmin,nage);
  dvar_vector& tm=(dvar_vector&)ctm;
  int yr=year(ir,ip);  
  dvar_vector MM=mfexp(get_nat_mort_region(ir)(yr)(jmin,nage)+fraction(ir,ip));
  MY_DOUBLE_TYPE normd=0.0;
  int noprogressflag=0;
  if (region_tot_tag_catch>0)
  {
    dvariable surv_rate=(1.0-region_tot_tag_catch/(1.e-8+tnf));
  
    dvar_vector kc(1,nfi);

    kc=actual_tag_catch+1.e-5;
    
    dvar_vector qq=elem_div(kc,sel*(1.e-8+enf));
    dvar_vector TC(1,nfi);
    int itt=0;
    int badflag=0;
    MY_DOUBLE_TYPE rmax=0.7;
    if (age_flags(116)!=0)
    {
      rmax=age_flags(116)/100.;
    }
    dvar_vector Z(jmin,nage);
    do
    {
      Z=qq*sel+MM;
      dvar_vector S(jmin,nage);
      int Zmaxflag=0;
      int ia;
      for (ia=jmin;ia<=nage;ia++)
      {
        if (value(Z(ia))<=rmax)
        {
          S(ia)=exp(-Z(ia));
        }
        else
        {
          Zmaxflag=1;
          dvariable dd=Z(ia)-rmax;
          S(ia)=exp(-rmax)-exp(-rmax)*dd;
        }
      }
      dvar_vector S1N=elem_prod(1.0-S,enf);
      dvar_matrix t1(1,nfi);
      for (int fi=1;fi<=nfi;fi++)
      {
        t1(fi)=elem_div(qq(fi)*sel(fi),Z);
        C(fi)=elem_prod(t1(fi),S1N);
        TC(fi)=sum(C(fi));
      }
      /*
      if (Zmaxflag==0)
      {
        for (int fi=1;fi<=nfi;fi++)
        {
          dvariable Dq=sum(C(fi))/qq(fi);
          dvar_vector DZ=qq(fi)*elem_prod(elem_div(sel(fi),Z),enf)-
            elem_prod(C(fi),1.0+1.0/Z);
          M(fi,fi)=Dq;
          for (int fj=1;fj<=nfi;fj++)
          {
            M(fi,fj)+=DZ*sel(fj);
          }
        }
      }
      else
      {
        for (int fi=1;fi<=nfi;fi++)
        {
          dvariable Dq=sum(C(fi))/qq(fi);
          dvar_vector DZ(jmin,nage);
          for (ia=jmin;ia<=nage;ia++)
          {
            if (value(Z(ia))<=rmax)
            {
              DZ(ia)=qq(fi)*sel(fi,ia)/Z(ia)*enf(ia)-
                C(fi,ia)*(1.0+1.0/Z(ia));
            }
            else
            {
              DZ(ia)=qq(fi)*sel(fi,ia)/Z(ia)*enf(ia)-
                C(fi,ia)/Z(ia)
                +C(fi,ia)/(1.0-S(ia))*exp(-rmax);
            }
          }
          M(fi,fi)=Dq;
          for (int fj=1;fj<=nfi;fj++)
          {
            M(fi,fj)+=DZ*sel(fj);
          }
        }
      }
      */
      int repeat_flag=0;
      do
      {
        dvar_matrix M2=ADJ_get_jac2(nfi,jmin,nage,qq,sel,enf,MM,rmax);
        //dvar_matrix M2=get_jac2(nfi,jmin,nage,qq,sel,enf,MM,rmax);
        M=M2;
      }
      while(repeat_flag);
     /*
      cout <<"tag3.cpp " << M2-M << endl << endl;
      cout <<"tag3.cpp " << trans(M2)-M << endl << endl;
      cout <<"tag3.cpp " << setw(10) << setfixed() << setprecision(3) << norm2(M2-M) << endl;
      */
      dvar_vector diff=TC-kc;
      normd=norm(value(diff));
      //cout <<"tag3.cpp " <<  elem_div(diff,TC)  << endl;
      int pflag=1;
      if (!badflag)
      {
        // this is newton raphson for q
        dvariable fp1=0.0;
        //dvar_matrix invM=inv(M);
        //ofsff << "M" << endl;
        //ofsff << setw(10) << setprecision(4) << M << endl;
        //ofsff << "invM" << endl;
        //ofsff << setscientific() << setw(10)<< setprecision(4) << invM << endl;
        dvector cqold=value(qq);
        dvar_vector hh=solve(M,diff);
        qq-=hh;
        if (itt>3)
        {
          if (norm(value(qq)-cqold)<1.e-20)
          {
            noprogressflag=1;
          }
        }
        
        int mmin=qq.indexmin();
        int mmax=qq.indexmax();
        for (int i=mmin;i<=mmax;i++)
        {
          qq(i)=posfun(qq(i),1.e-10,fp1);
          if (value(qq(i))<=0.0)
          {
            cerr << "This can't happen " << qq(i) << endl;
            ad_exit(1);
          }
        }
        //ffpen+=fp1;
        //pflag=check_pos(qq,nfi);
        if (!pflag)
        {
          dvariable tmp=0.0;
          throw function_minimizer_exception();
          badflag=1;
          qq+=inv(M)*diff;
          qq=elem_div(qq,mfexp(elem_div(inv(M)*diff,qq)));
        }
        else
        {
          itt++;
        }
      }
      else
      {
        // this i newton raphson for log(q)
        qq=elem_div(qq,1.e-20+mfexp(elem_div(inv(M)*diff,qq)));
        //cout << " " << qq(1);
        itt++;
        pflag=check_pos(qq,nfi);
      }
      if (itt>7) break;
      if (noprogressflag==1) 
      {
        noprogressflag=0;
        break;
      }
      //if ( (itt>1 && pflag && !badflag) || (itt>3 && pflag) )break;
    }
    while(1);
    if (normd>0.5)
    {
      cout << "difficult tag NR A " << normd;
      cout << "  it= " << it ;
      cout << "  region = " << ir ;
      cout << "  fishing period = " << ip << endl;
      cout << " There are " << nfi << " fisheries ";
      for (int i=1;i<=nfi;i++)
        cout << parent(ir,ip,i) << " ";
      cout << endl;
    }
    //cout <<"tag3.cpp " << endl;
    int j;
    //dvar_vector lambda(jmin,nage);
    //lambda.initialize();
    for (j=jmin;j<=nage;j++)
    {
      if (value(Z(j))>rmax)
      {
        if (value(Z(j)) > maxZ[1])
        {
          maxZ[1]=value(Z(j));
          maxZ[2]=it;
          maxZ[3]=ir;
          maxZ[4]=ip;
        }
        age_flags(181)=2;
        //fpen+=100.*square(Z(j)-rmax);
        dvariable dd=Z(j)-rmax;
        dvariable ppen=0.0;
        dvariable Zstar=rmax-posfun(0.2-dd,0.1,ppen)+0.2;
        ffpen+=ppen;

	/*
        //if(it == 2 && ir == 3)  //NMD1Oct2012
	if(it == 51 && ir == 9)  //NMD1Oct2012
	  ofsff << endl << " it = " << it << "   ir = " << ir << "   j = " << j << endl; endl;
	  ofsff << " Zstar = " << Zstar << endl;
	  ofsff << " dd = " << dd << endl;
	  ofsff << " rmax = " << rmax << endl;
	  ofsff << " dd = " << dd << endl;
	  ofsff << " ppen = " << ppen << endl;
	  ofsff << " MM(j) = " << MM(j) << endl;
	  ofsff << " Zstar-MM(j) = " << Zstar-MM(j) << endl;
          //NMD1Oct2012
	*/
  
        dvariable tmp=Zstar-MM(j);
        if (value(tmp)<=0.0)
        {
          ffpen+=100.0*square(log(1.0-tmp));
          //cerr << "Error zstar too small" << endl;
          //ad_exit(1);
        }
        //lambda(j)=log(tmp/(Z(j)-MM(j)));
      }
    }
    dvar_vector tmpFF=log(qq);
    for (fi=1;fi<=nfi;fi++)
    {
      fm(fi)(jmin,nage)=logsel(fi)+tmpFF(fi);  // +lambda;
    }
    tm.initialize();
    for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
    {                                         // incidents for this period
      tm+=mfexp(fm(fi));
    }
    tm+=MM;
  }
  else
  {
    for (fi=1;fi<=nfi;fi++)
    {
      fm(fi)=-10.0;
    }
    tm=MM;
  }
  nrsurv(it,ir,ip)=mfexp(-tm);
  switch (it)
  {
  case 72:
  case 21:
  case 112:
    break;
  default:
    ;
    _ffpen+=ffpen;
  }
}


void dvar_fish_stock_history::allocate_some_tag_stuff(void)
{
  tag_return_probability.allocate(1,num_tag_releases);
  
  if (!age_flags(96)) 
  {
    min_tag_age4.allocate(1,num_tag_releases, // 1,num_regions,
      tag_region_bounds(1),tag_region_bounds(2),
      initial_tag_period,tag_num_fish_periods,
      //initial_tag_period,imatrix(1,num_tag_releases,num_fish_periods),
      1,num_tagfish_incidents);
  
  /*
   // for (int it=1;it<=num_tag_releases;it++)
   // {
   //   for (int ir=1;ir<=num_regions;ir++)
   //   {
   //     int min_age=1;
   //     min_tag_age4(it,ir,initial_tag_period(it,ir))=min_age;
   //     for (int ip=initial_tag_period(it,ir)+1;
   //       ip<=num_fish_periods(ir);ip++)
   //     {
   //       if (year(ir,ip)>year(ir,ip-1)) 
   //       {
   //         if (min_age<nage) min_age++;
   //       }
   //       for (int fi=1;fi<=num_tagfish_incidents(it,ir,ip);fi++)
   //       {
   //         min_tag_age4(it,ir,ip,fi)=min_age;
   //       }
   //     }
   //   }
   // }
    */
    for (int it=1;it<=num_tag_releases;it++)
    {
      //for (int ir=1;ir<=num_regions;ir++)
      for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
      {
        int min_age=1;
        for (int ip=initial_tag_period(it,ir);
          ip<=num_fish_periods(ir);ip++)
        {
          for (int fi=1;fi<=num_tagfish_incidents(it,ir,ip);fi++)
          {
            min_tag_age4(it,ir,ip,fi)=
              min(nage_by_tag_release(it),
              min_age+year(ir,ip)-initial_tag_year(it));
          }
        }
      }
    }
    
   /*
    ivector rmin(1,num_tag_releases);
    ivector rmax(1,num_tag_releases);
    if (pmsd)
    {
      for (it=1;it<=num_tag_releases;it++)
      {
        int cs=pmsd->tag_species_pointer(it);
        rmin(it)=pmsd->region_bounds(cs,1);
        rmax(it)=pmsd->region_bounds(cs,2);
      }
    }
    else
    {
      rmin=1;
      rmax=num_regions;
    }
    */
    //int rmin=1;
    //int rmax=num_regions;

    tagcatch.allocate(1,num_tag_releases,
       //  rmin,rmax,
      tag_region_bounds(1),tag_region_bounds(2),
      initial_tag_period,tag_num_fish_periods,
      //imatrix(1,num_tag_releases,num_fish_periods),
      1,num_tagfish_incidents,min_tag_age4,nage_by_tag_release);

    obstagcatch.allocate(1,num_tag_releases,
      //rmin,rmax,
      tag_region_bounds(1),tag_region_bounds(2),
      initial_tag_period,tag_num_fish_periods,
      //imatrix(1,num_tag_releases,num_fish_periods),
      1,num_tagfish_incidents,min_tag_age4,nage_by_tag_release);

    obstagcatch1.allocate(1,num_tag_releases,
      //rmin,rmax,
      tag_region_bounds(1),tag_region_bounds(2),
      initial_tag_period,tag_num_fish_periods,
      //imatrix(1,num_tag_releases,num_fish_periods),
      1,num_tagfish_incidents,min_tag_age4,nage_by_tag_release);
  }
  if (num_tag_releases)
  {
    if (sum(column(tag_flags,1)))
    {
      ivector tf1=column(tag_flags,1);
      imatrix tmp(1,num_tag_releases,tag_region_bounds(1),tag_region_bounds(2));
      tmp=initial_tag_period;
      int it;
      for (it=1;it<=num_tag_releases;it++)
      {
        if (tf1(it)) tmp(it)+=tf1(it)-1;
        for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
        {
          if (tmp(it,ir)>num_fish_periods(ir)) {
            tmp(it,ir)=num_fish_periods(ir);
        }
          }
      }
      i3_array itmp(1,num_tag_releases,tag_region_bounds(1),
        tag_region_bounds(2),initial_tag_period,tmp);
      for (it=1;it<=num_tag_releases;it++)
      {
        for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
        {
          for (int ip=initial_tag_period(it,ir);ip<=tmp(it,ir);ip++)
          {
            itmp(it,ir,ip)=num_fish_incidents(ir,ip);
          }
        }
      }

      min_tag_age2.allocate(1,num_tag_releases,
        tag_region_bounds(1),tag_region_bounds(2),
        initial_tag_period,tmp,1,itmp);

      /*
      for (it=1;it<=num_tag_releases;it++)
      {
        for (int ir=1;ir<=num_regions;ir++)
        {
          int min_age=1;
          min_tag_age2(it,ir,initial_tag_period(it,ir))=min_age;
          for (int ip=initial_tag_period(it,ir)+1;ip<=tmp(it,ir);ip++)
          {
            if (year(ir,ip)>year(ir,ip-1)) 
            {
              if (min_age<nage) min_age++;
            }
            for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)
            {
              min_tag_age2(it,ir,ip,fi)=min_age;
            }
          }
        }
      }
     */
      for (it=1;it<=num_tag_releases;it++)
      {
        for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
        {
          int min_age=1;
          for (int ip=initial_tag_period(it,ir);ip<=tmp(it,ir);ip++)
          {
            for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)
            {
              min_tag_age2(it,ir,ip,fi)=
                min(nage_by_tag_release(it),
                min_age+year(ir,ip)-initial_tag_year(it));
            }
          }
        }
      }
    
      min_tag_age3.allocate(1,num_tag_releases,
        tag_region_bounds(1),tag_region_bounds(2),
        initial_tag_period,tmp);
     /*
      for (it=1;it<=num_tag_releases;it++)
      {
        for (int ir=1;ir<=num_regions;ir++)
        {
          int min_age=1;
          min_tag_age3(it,ir,initial_tag_period(it,ir))=min_age;
          for (int ip=initial_tag_period(it,ir)+1;ip<=tmp(it,ir);ip++)
          {
            if (year(ir,ip)>year(ir,ip-1)) 
            {
              if (min_age<nage) min_age++;
            }
            min_tag_age3(it,ir,ip)=min_age;
          }
        }
      }
     */

      for (it=1;it<=num_tag_releases;it++)
      {
        for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
        {
          int min_age=1;
          for (int ip=initial_tag_period(it,ir);ip<=tmp(it,ir);ip++)
          {
            min_tag_age3(it,ir,ip)=
              min(nage_by_tag_release(it),
              min_age+year(ir,ip)-initial_tag_year(it));
          }
        }
      }

      const index_type& tmp1=index_type(min_tag_age2);
      const index_type& tmp2=tmp1(1)(1)(initial_tag_period(1,1))(itmp(1,1,initial_tag_period(1,1)));
      //const ad_integer& tmp3=tmp2;
      cout <<"tag3.cpp " << min_tag_age2(1)(1)(initial_tag_period(1,1))(itmp(1,1,initial_tag_period(1,1))) << endl;
    
      nrfm.allocate(1,num_tag_releases,
        tag_region_bounds(1),tag_region_bounds(2),
        initial_tag_period,tmp,
        1,itmp,index_type(min_tag_age2),nage);
      nrtm.allocate(1,num_tag_releases,
        tag_region_bounds(1),tag_region_bounds(2),
        initial_tag_period,tmp,
        min_tag_age3,nage);
      nrsurv.allocate(1,num_tag_releases,
        tag_region_bounds(1),tag_region_bounds(2),
        initial_tag_period,tmp,
        min_tag_age3,nage);
    }
  }
}


void dvar_fish_stock_history::print_tag_return_by_time_at_liberty(ofstream& of)
{
  int maxtal=0;
  int ir;
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=minimum_initial_tag_period(ir);ip<=num_fish_periods(ir);ip++) 
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        for (int it=1;it<=num_tag_releases;it++)
        {
//          if (ir<tag_region_bounds(1,it) || ir >tag_region_bounds(2,it)) 
//              break;
          if (ir >= tag_region_bounds(1,it) && ir <= tag_region_bounds(2,it))
          {       //NMD_13Sep2018
            if (age_flags(96))
            {    
              if (initial_tag_period(it,ir)<=ip
                && terminal_tag_period(it,ir)>=ip)
              {
                int itp=initial_tag_period(it,ir);
                int tal=12*(true_year(ir,ip)-true_year(ir,itp))
                  +(true_month(ir,ip)-true_month(ir,itp));
                tal=(tal)/3+1;
                if (tal>maxtal) maxtal=tal;
              }
            }
            else
            {
              if (initial_tag_period(it,ir)<=ip)
              {
                int itp=initial_tag_period(it,ir);
                int tal=12*(true_year(ir,ip)-true_year(ir,itp))
                  +(true_month(ir,ip)-true_month(ir,itp));
                tal=(tal)/3+1;
                if (tal>maxtal) maxtal=tal;
              }
            }
          }   //NMD_13Sep2018
        }
      }
    }
  }
  dvector ob(1,maxtal);
  dvector pr(1,maxtal);
  ob.initialize();
  pr.initialize();
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=minimum_initial_tag_period(ir);ip<=num_fish_periods(ir);ip++) 
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        for (int it=1;it<=num_tag_releases;it++)
        {
//          if (ir<tag_region_bounds(1,it) || ir >tag_region_bounds(2,it)) 
//              break;
          if (ir >= tag_region_bounds(1,it) && ir <= tag_region_bounds(2,it))
          {       //NMD_13Sep2018
            if (age_flags(96))
            {    
              if (initial_tag_period(it,ir)<=ip
                && terminal_tag_period(it,ir)>=ip)
              {
                int itp=initial_tag_period(it,ir);
                int tal=12*(true_year(ir,ip)-true_year(ir,itp))
                  +(true_month(ir,ip)-true_month(ir,itp));
                tal=(tal)/3+1;
                ob(tal)
                  +=sum(value(obstagcatch(it,ir,ip,fi)));
                pr(tal)
                  +=sum(get_rep_rate_correction(it,ir,ip,fi)
                    *value(tagcatch(it,ir,ip,fi)));
                  //+=sum(value(rep_rate(ir,ip,fi))*value(tagcatch(it,ir,ip,fi)));
              }
            }
            else
            {
              if (initial_tag_period(it,ir)<=ip)
              {
                int itp=initial_tag_period(it,ir);
                int tal=12*(true_year(ir,ip)-true_year(ir,itp))
                  +(true_month(ir,ip)-true_month(ir,itp));
                if (tal<0)
                {
                   cout << "tal =  " << tal << " true year(ir,ip) = "
                        << true_year(ir,ip) << endl
                        << " true_year(ir,itp) = "
                        << true_year(ir,itp) 
                        << " true month(ir,ip) = "
                        << true_month(ir,ip) << endl
                        << " true_month(ir,itp) = "
                        << true_month(ir,itp) << endl;
                }   
                tal=(tal)/3+1;
                ob(tal)
                  +=sum(value(obstagcatch(it,ir,ip,fi)));
                pr(tal)
                  +=sum(get_rep_rate_correction(it,ir,ip,fi)
                    *value(tagcatch(it,ir,ip,fi)));
                  //+=sum(value(rep_rate(ir,ip,fi))*value(tagcatch(it,ir,ip,fi)));
              }
            }
          }  //NMD_13Sep2018
        }
      }
    }
  }
  of << "# Maximum time at liberty" << endl;
  of << "  " << maxtal << endl;
  of << "# Observed vs predicted tag returns by time at liberty" << endl;
  for (int i=1;i<=maxtal;i++)
  {
    of << "   " << setfixed() << setprecision(0) << setw(5) << ob(i) << " ";
    of << "   " << setfixed() << setprecision(1) << setw(7) << pr(i) << endl;
  }
}

MY_DOUBLE_TYPE dvar_fish_stock_history::get_rep_rate_correction
  (int it,int ir,int ip,int fi)
{
  // check if we are using tag newton raphson for this period
  if (!num_fish_incidents(ir,ip) || !tag_flags(it,1)
     || ip >= initial_tag_period(it,ir)+tag_flags(it,1))
  {
    if (age_flags(198))
    {
      return value(tag_rep_rate(it,ir,ip,fi));
    }
    else
    {
      return value(rep_rate(ir,ip,fi));
    }
  }
  else  // doing the newton raphson
  {
    if (tag_flags(it,2)==1)
    {
      return 1.0;
    }
    else
    {
      if (age_flags(198))
      {
        return value(tag_rep_rate(it,ir,ip,fi));
      }
      else
      {
        return value(rep_rate(ir,ip,fi));
      }
    }
  }
}

void dvar_fish_stock_history::print_movement_report(ofstream& of)
{
  of << "# Movement analysis" << endl;
  if (!pmsd)
  {
    dvector N(1,num_regions);
    int nmp=mo.num_periods();
    for (int ir=1;ir<=num_regions;ir++)
    {
      N.initialize();
      N(ir)=1000.0;
      of << "# Region " << ir << endl;
      of << setfixed() << setprecision(1) << setw(5) << N << endl;
      int j;
      for (j=1;j<nage;j++)
      {
        for (int mp=1;mp<=nmp;mp++)
        {
          N=value(Dad(move_map(mp),j)*N);
        }
        if (j<nage-1)
          of << setfixed() << setprecision(1) << setw(5) << N << endl;
      }
      for (j=1;j<=nage+2;j++)
      {
        for (int mp=1;mp<=nmp;mp++)
        {
          N=value(Dad(move_map(mp),nage)*N);
        }
        of << setfixed() << setprecision(1) << setw(5) << N << endl;
      }
    }
  }
  else
  {
    
    int nrr=pmsd->num_real_regions;
    int ns=pmsd->num_species;
    for (int isp=1;isp<=ns;isp++)
    {
      dvector N(1,nrr);
      dvector Zeroes(1,nrr);
      Zeroes.initialize();
      int nmp=mo.num_periods();
      int offset=(isp-1)*nrr;
      for (int ir=1;ir<=nrr;ir++)
      {
        N.initialize();
        N(ir)=1000.0;
        of << "# Region " << ir +offset << endl;
        // printing extra zeroes so that MF viewer will be happy
        for (int is1=1;is1<=ns;is1++)
        {
          if (is1==isp)
            of << setfixed() << setprecision(1) << setw(5) << N << " ";
          else
            of << setfixed() << setprecision(1) << setw(5) << Zeroes << " ";
        }
        of << endl;
        int j;
        for (j=1;j<nage;j++)
        {
          for (int mp=1;mp<=nmp;mp++)
          {
            dvar_matrix Dd=
              Dad(move_map(mp),j).sub(1+offset,nrr+offset).shift(1);
            N=value(Dd)*N;
          }
          if (j<nage-1)
          {
            for (int is1=1;is1<=ns;is1++)
            {
              if (is1==isp)
                of << setfixed() << setprecision(1) << setw(5) << N << " ";
              else
                of << setfixed() << setprecision(1) << setw(5) << Zeroes << " ";
            }
            of << endl;
          }
        }
        for (j=1;j<=nage+2;j++)
        {
          for (int mp=1;mp<=nmp;mp++)
          {
            dvar_matrix Dd=
              Dad(move_map(mp),nage).sub(1+offset,nrr+offset).shift(1);
            N=value(Dd)*N;
          }
          for (int is1=1;is1<=ns;is1++)
          {
            if (is1==isp)
              of << setfixed() << setprecision(1) << setw(5) << N << " ";
            else
              of << setfixed() << setprecision(1) << setw(5) << Zeroes << " ";
          }
          of << endl;
        }
      }
    }
  }
}
void dvar_fish_stock_history::do_fish_mort_intermediate_calcs(void)
{
  const MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
  dvar_matrix * ptmp=0;
      
  //if (age_flags(92) !=3)
  switch(age_flags(92))
  {
  case 0:
  case 2:
    {
      for (int ir=1;ir<=num_regions;ir++)
      {
        for (int ip=1;ip<=num_real_fish_periods(ir);ip++)  
        {
          int yr=year(ir,ip);
          if (age_flags(125) &&year(ir,ip)>1) 
            break;
          dvar_matrix* pfm;
          dvar_vector* ptotmort;
          dvar_vector* psurv;
          if (af170q0==0)
          {
            ptotmort=&(tot_mort(ir,ip));
            psurv=&(survival(ir,ip));
            ptmp=&(fish_mort_calcs(ir,ip));
            pfm=&(fish_mort(ir,ip));
          }
          else
          {
            ptotmort=&(tot_mort_q0(ir,ip));
            psurv=&(survival_q0(ir,ip));
            ptmp=&(fish_mort_calcs_q0(ir,ip));
            pfm=&(fish_mort_q0(ir,ip));
          }

          if (age_flags(115)==0)
          {
            const dvar_vector& tmp1=log(1.e-10+ *ptotmort);
            const dvar_vector& tmp2=log(one_plus- *psurv);
            const dvar_vector& tmp3=tmp2-tmp1;
            for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
            {        
              (*ptmp)(fi)=(*pfm)(fi)+tmp3;
            }
          }
          else   // ss2 parameterization
          {
            for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
            {        
              (*ptmp)(fi)=(*pfm)(fi)-0.5*e_nat_mort_by_period(ir,ip);
            }
          }
        }
      } 
    }
    break;
  case 4:
    {
      for (int ir=1;ir<=num_regions;ir++)
      {
        for (int ip=1;ip<=num_real_fish_periods(ir);ip++)  
        {
          if (age_flags(125) &&year(ir,ip)>1) 
            break;
          if (af170q0==0)
          {
            ptmp=&(fish_mort_calcs(ir,ip));
          }
          else
          {
            ptmp=&(fish_mort_calcs_q0(ir,ip));
          }
          dvar_vector em=e_nat_mort_by_period(ir,ip);
          const dvar_vector& zz=tot_mort(ir,ip)-em;
          int mmin=zz.indexmin();
          int mmax=zz.indexmax();
          for (int j=mmin;j<=mmax;j++)
          {
            if (zz(j)<0.0)
                cout<< "HERE 0" << endl;
          }
           
          const dvar_vector& tmp1=log(1.e-10+zz);
          const dvar_vector& tmp2=log(one_plus-exp(-zz));
          const dvar_vector& tmp3=tmp2-tmp1;
          dvar_matrix* pfm;
          if (af170q0)
          {
            pfm=&(fish_mort_q0(ir,ip));
          }
          else
          {
            pfm=&(fish_mort(ir,ip));
          }
          for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
          {        
            (*ptmp)(fi)=(*pfm)(fi)+tmp3-0.5*em;
          }
        }
      } 
    }
    break;
  case 3:
    {
      for (int ir=1;ir<=num_regions;ir++)
      {
        for (int ip=1;ip<=num_real_fish_periods(ir);ip++)  
        {
          if (age_flags(125) &&year(ir,ip)>1) 
            break;
          if (af170q0==0)
          {
            ptmp=&(fish_mort_calcs(ir,ip));
          }
          else
          {
            ptmp=&(fish_mort_calcs_q0(ir,ip));
          }
          dvar_matrix* pfm;
          if (af170q0)
          {
            pfm=&(fish_mort_q0(ir,ip));
          }
          else
          {
            pfm=&(fish_mort(ir,ip));
          }
          for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
          {        
            (*ptmp)(fi)=(*pfm)(fi)-0.5*e_nat_mort_by_period(ir,ip);
          }
        }
      } 
    }
    break;
  default: 
    cerr << "Illegal value of " << age_flags(92) << " for age_flags(92)" 
         << endl;
    ad_exit(1);
  }
}

void dvar_fish_stock_history::do_fish_mort_intermediate_projection_calcs(void)
{
  const MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
  dvar_matrix * ptmp=0;
      
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=num_real_fish_periods(ir)+1;ip<=num_fish_periods(ir);ip++)
    {
      if (af170q0==0)
      {
        ptmp=&(fish_mort_calcs(ir,ip));
      }
      else
      {
        ptmp=&(fish_mort_calcs_q0(ir,ip));
      }
      const dvar_vector& tmp1=log(1.e-10+tot_mort(ir,ip));
      const dvar_vector& tmp2=log(one_plus-survival(ir,ip));
      const dvar_vector& tmp3=tmp2-tmp1;
      dvar_matrix* pfm;
      if (af170q0)
      {
        pfm=&(fish_mort_q0(ir,ip));
      }
      else
      {
        pfm=&(fish_mort(ir,ip));
      }
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {        
        (*ptmp)(fi)=(*pfm)(fi)+tmp3;
      }
    }
  } 
}

void dvar_fish_stock_history::do_fish_mort_intermediate_calcs(int ir,int ip)
{
  const MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
  const dvar_vector& tmp1=log(1.e-10+tot_mort(ir,ip));
  const dvar_vector& tmp2=log(one_plus-survival(ir,ip));
  const dvar_vector& tmp3=tmp2-tmp1;
  for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
  {        
//    fish_mort_calcs_q0(ir,ip,fi)=fish_mort(ir,ip,fi)+tmp3;  //NMD_2Dec2019
    fish_mort_calcs(ir,ip,fi)=fish_mort(ir,ip,fi)+tmp3;
  }
}

void dvar_fish_stock_history::do_tag_loss_fish_mort_intermediate_calcs(int ir,int ip)
{
  const MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
  int num_tag_loss_groups=tag_mortality_group.get_ngroups();
  ivector inv_group_ptr=tag_mortality_group.get_inv_group_ptr();
  //for (int ig=1;ig<=num_tag_loss_groups;ig++) 
  for (int ig=1;ig<=1;ig++) 
  {        
    const dvar_vector& tmp1=log(1.e-10+tot_mort(ir,ip))
      +tagmort(1);
      //+tagmort(inv_group_ptr(ig));
    const dvar_vector& tmp2=log(one_plus-survival(ir,ip))
      -tagmort(1);
      //-tagmort(inv_group_ptr(ig));
    const dvar_vector& tmp3=tmp2-tmp1;
    
    for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
    {        
      tag_fish_mort_calcs(ig,ir,ip,fi)=fish_mort(ir,ip,fi)+tmp3;
    }
  }
}

void dvar_fish_stock_history::do_fish_mort_intermediate_calcs_q0(int ir,int ip)
{
  const MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
  const dvar_vector& tmp1=log(1.e-10+tot_mort_q0(ir,ip));
  const dvar_vector& tmp2=log(one_plus-survival_q0(ir,ip));
  const dvar_vector& tmp3=tmp2-tmp1;
  for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
  {        
//    fish_mort_calcs(ir,ip,fi)=fish_mort_q0(ir,ip,fi)+tmp3;   //NMD_2Dec2019
    fish_mort_calcs_q0(ir,ip,fi)=fish_mort_q0(ir,ip,fi)+tmp3;
  }
}

void dvar_fish_stock_history::read_recruitment_env_data(adstring& root)
{
  adstring tmpstring=root+adstring(".rnv");
  cifstream cifs((char *)(tmpstring));

  if (!cifs)
  {
    cerr << "Error trying to open file " << tmpstring << endl;
    exit(1);
  }
  cifs >> rec_covars;
  if (!cifs)
  {
    cerr << "Error reading recruitment environment data file " << tmpstring << endl;
    exit(1);
  }
  // normalize the environment data
  rec_covars-=mean(rec_covars);
  rec_covars/=(1.e-10+std_dev(rec_covars));
}

void check_region_number(int tag_region,int nu,int it,int yr,int mn,
  int num_regions)
{
  if (tag_region<1 || tag_region> num_regions)
  {
    cerr << "Illegal region for tag release " << it << endl;
    cerr << "tag_region = " << tag_region << endl;
    cerr << "tag_year = " << yr << endl;
    cerr << "tag_month = " << mn << endl;
    ad_exit(1);
  }
}

void dvar_fish_stock_history::new_print_tag_return_by_time_at_liberty(ofstream& of)
{
  int maxtal=0;
  int ir;
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=minimum_initial_tag_period(ir);ip<=num_fish_periods(ir);ip++) 
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        for (int it=1;it<=num_tag_releases;it++)
        {
          if (age_flags(96))
          {    
            if (initial_tag_period(it,ir)<=ip
              && terminal_tag_period(it,ir)>=ip)
            {
              int itp=initial_tag_period(it,ir);
              int tal=12*(true_year(ir,ip)-true_year(ir,itp))
                +(true_month(ir,ip)-true_month(ir,itp));
              tal=(tal)/3+1;
              if (tal>maxtal) maxtal=tal;
            }
          }
          else
          {
            if (initial_tag_period(it,ir)<=ip)
            {
              int itp=initial_tag_period(it,ir);
              int tal=12*(true_year(ir,ip)-true_year(ir,itp))
                +(true_month(ir,ip)-true_month(ir,itp));
              tal=(tal)/3+1;
              if (tal>maxtal) maxtal=tal;
            }
          }
        }
      }
    }
  }
  dvector ob(1,maxtal);
  dvector pr(1,maxtal);
  ob.initialize();
  pr.initialize();
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=minimum_initial_tag_period(ir);ip<=num_fish_periods(ir);ip++) 
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        for (int it=1;it<=num_tag_releases;it++)
        {
          if (age_flags(96))
          {    
            if (initial_tag_period(it,ir)<=ip
              && terminal_tag_period(it,ir)>=ip)
            {
              int itp=initial_tag_period(it,ir);
              int tal=12*(true_year(ir,ip)-true_year(ir,itp))
                +(true_month(ir,ip)-true_month(ir,itp));
              tal=(tal)/3+1;
              ob(tal)
                +=sum(value(obstagcatch(it,ir,ip,fi)));
              pr(tal)
                +=sum(value(tag_rep_rate(it,ir,ip,fi))*value(tagcatch(it,ir,ip,fi)));
            }
          }
          else
          {
            if (initial_tag_period(it,ir)<=ip)
            {
              int itp=initial_tag_period(it,ir);
              int tal=12*(true_year(ir,ip)-true_year(ir,itp))
                +(true_month(ir,ip)-true_month(ir,itp));
              if (tal<0)
              {
                 cout << "tal =  " << tal << " true year(ir,ip) = "
                      << true_year(ir,ip) << endl
                      << " true_year(ir,itp) = "
                      << true_year(ir,itp) 
                      << " true month(ir,ip) = "
                      << true_month(ir,ip) << endl
                      << " true_month(ir,itp) = "
                      << true_month(ir,itp) << endl;
              }   
              tal=(tal)/3+1;
              ob(tal)
                +=sum(value(obstagcatch(it,ir,ip,fi)));
              pr(tal)
                +=sum(value(tag_rep_rate(it,ir,ip,fi))*value(tagcatch(it,ir,ip,fi)));
            }
          }
        }
      }
    }
  }
  of << "# Maximum time at liberty" << endl;
  of << "  " << maxtal << endl;
  of << "# Observed vs predicted tag returns by time at liberty" << endl;
  for (int i=1;i<=maxtal;i++)
  {
    of << "   " << setfixed() << setprecision(0) << setw(5) << ob(i) << " ";
    of << "   " << setfixed() << setprecision(1) << setw(7) << pr(i) << endl;
  }
}
void dvar_fish_stock_history::
  do_newton_raphson_for_tags(int it,int ir,int ip,dvariable& _ffpen)
{
  dvariable ffpen=0.0;
  int jmin=tagnum_fish(it,ir,ip).indexmin();
  const dvar_vector& enf=mfexp(tagnum_fish(it,ir,ip));
  dvariable tot_num_fish=sum(enf);
  int nfi=num_fish_incidents(ir,ip);
  dvar_matrix sel(1,nfi,jmin,nage);
  dvar_matrix logsel(1,nfi,jmin,nage);
  int fi;
  for (fi=1;fi<=nfi;fi++)
  {
    int i=parent(ir,ip,fi);
    int rr=realization_region(i,1);
    int rp=realization_period(i,1);
    int ri=realization_incident(i,1);
    logsel(fi)=incident_sel(rr,rp,ri)(jmin,nage);
    sel(fi)=mfexp(logsel(fi));
  }

  dvar_matrix& fm=nrfm(it,ir,ip);
  dvar_matrix M(1,nfi,1,nfi);
  dvar_matrix C(1,nfi,jmin,nage);
  dvariable tnf=sum(enf);
  //ofsff << tnf << " " << it << " " << ir << " " << ip << endl;
  if (nfi>1)
  {
    //cout <<"tag3.cpp " << nfi << endl;
  }

  dvar_vector actual_tag_catch(tot_tag_catch(it,ir,ip).indexmin(),
    tot_tag_catch(it,ir,ip).indexmax());

  if (!num_fish_incidents(ir,ip) || !tag_flags(it,1) 
     || ip >= initial_tag_period(it,ir)+tag_flags(it,1))
  {
    cout << "need to deal with this case" << endl;
    ad_exit(1);
  }
  else
  {
    if (tag_flags(it,2)==0)
    {
      if (age_flags(198))
      {
        int nfi=num_fish_incidents(ir,ip);
        dvar_vector tmp(1,nfi);
        for (fi=1;fi<=nfi;fi++)
        {
          tmp(fi)=tag_rep_rate(it,ir,ip,fi);
        }
        actual_tag_catch=
          elem_div(tot_tag_catch(it,ir,ip),1.e-6+tmp);
      }
      else
      {
        actual_tag_catch=
          elem_div(tot_tag_catch(it,ir,ip),1.e-6+rep_rate(ir,ip));
      }
    }
    else
    {
      actual_tag_catch=tot_tag_catch(it,ir,ip);
    }
  }
  dvariable region_tot_tag_catch=sum(actual_tag_catch);

  const dvar_vector& ctm=nrtm(it,ir,ip)(jmin,nage);
  dvar_vector& tm=(dvar_vector&)ctm;
  int break_flag=0;
  if (region_tot_tag_catch>0)
  {
    dvariable surv_rate=(1.0-region_tot_tag_catch/tnf);
  
    dvar_vector kc(1,nfi);


    MY_DOUBLE_TYPE cut=0.2;
    MY_DOUBLE_TYPE fringe=0.02;
    if (age_flags(118)>0)
    {
      cut=age_flags(118)/100.;
    }
    if (age_flags(119)>0)
    {
      fringe=age_flags(119)/100.;
    }
    MY_DOUBLE_TYPE c2=cut+fringe;

    
    if (surv_rate<=c2)
    {
      dvariable tmp=0.0;
      dvariable ks=fringe+posfun(surv_rate-fringe,cut,tmp);

      MY_DOUBLE_TYPE penwt=1.0;
      if (age_flags(117)>0) penwt=age_flags(117)/1.e+8;
      if (age_flags(117)<0) penwt=0.0;
        ffpen+=penwt*tmp;

      if (value(ks)>1.0)
      {
        cout << "fringe " << fringe << " surv rate  " << surv_rate 
             << " cut " << cut << endl;
        cout << "enter number " << endl;
        int num;
        cin >> num;
        cout << "this can't happen ks = " << ks << endl;
        ad_exit(1);
      }
      dvariable kr= (1.0-ks)*tnf/region_tot_tag_catch;
      kc=kr*actual_tag_catch+1.e-10;
    }
    else
    {
      kc=actual_tag_catch+1.e-10;
    }
    dvar_vector qq=elem_div(kc,sel*enf);
    dvar_vector TC(1,nfi);
    int itt=0;
    int badflag=0;
    MY_DOUBLE_TYPE normd=0.0;
    do
    {
      // moved M.initialize() to here
      M.initialize();
      int yr=year(ir,ip);  
      dvar_vector Z=qq*sel+mfexp(get_nat_mort_region(ir)(yr)(jmin,nage)
        +fraction(ir,ip));
      dvar_vector S=exp(-Z);
      dvar_vector S1N=elem_prod(1.0-S,enf);
      for (int fi=1;fi<=nfi;fi++)
      {
        dvar_vector t1=elem_div(qq(fi)*sel(fi),Z);
        C(fi)=elem_prod(t1,S1N);
        TC(fi)=sum(C(fi));
        dvariable Dq=sum(C(fi))/qq(fi);
        dvar_vector DZ=qq(fi)*elem_prod(elem_div(sel(fi),Z),enf)-
          elem_prod(C(fi),1.0+1.0/Z);
        M(fi,fi)=Dq;
        for (int fj=1;fj<=nfi;fj++)
        {
          M(fi,fj)+=DZ*sel(fj);
        }
      }
      for (int fi=1;fi<=nfi;fi++)
      {
        if (kc(fi)<1.e-5)
        {
          for (int fj=1;fj<=nfi;fj++)
          {
            M(fi,fj)=0.0;
            M(fj,fi)=0.0;
          }
          M(fi,fi)=1000.0;
        }
      }
      dvar_vector diff=TC-kc;
      normd=norm(value(diff));
      if (break_flag==1) break;
      //cout <<"tag3.cpp " <<  elem_div(diff,TC)  << endl;
      int pflag=1;
      if (!badflag)
      {
        // this is newton raphson for q
        dvar_vector hs=solve(M,diff);
        for (int fi=1;fi<=nfi;fi++)
        {
          if (kc(fi)<1.e-5) hs(fi)=0.0;
        }
        qq-=hs;
       /*
        int mmin=qq.indexmin();
        int mmax=qq.indexmax();
        dvariable fp1=0.0;
        for (int i=mmin;i<=mmax;i++)
        {
          qq(i)=posfun(qq(i),1.e-10,fp1);
        }
        ffpen+=fp1;
       */
        //pflag=check_pos(qq,nfi);
        if (!pflag)
        {
          dvariable tmp=0.0;
          throw function_minimizer_exception();
          badflag=1;
          qq+=inv(M)*diff;
          qq=elem_div(qq,mfexp(elem_div(inv(M)*diff,qq)));
        }
        else
        {
          itt++;
        }
      }
      else
      {
        // this i newton raphson for log(q)
        qq=elem_div(qq,1.e-20+mfexp(elem_div(inv(M)*diff,qq)));
        //cout << " " << qq(1);
        itt++;
        pflag=check_pos(qq,nfi);
      }
      if ( (itt>8 && pflag && !badflag) || (itt>8 && pflag) || normd<1.e-6 )break_flag=1;
    }
    while(1);
    dvariable fp1=0.0;
    int mmin=qq.indexmin();
    int mmax=qq.indexmax();
    for (int i=mmin;i<=mmax;i++)
    {
      qq(i)=posfun(qq(i),1.e-10,fp1);
    }
    ffpen+=fp1;
    if (normd>.05)
    {
      cout << "difficult tag NR B " << normd;
      cout << "  it= " << it ;
      cout << "  ir= " << ir ;
      cout << "  ip= " << ip << endl;
      cout << "Predicted tag catch           " << setw(10) << setprecision(2) << setfixed() << TC  << endl; 
      cout << "Expanded observed tag catch kc" << setw(10) << setprecision(2) << setfixed() <<  kc << endl;
      cout << "total number of tags present  " << setw(10) << setprecision(2) << setfixed() <<  tnf << endl;
    }
    //cout <<"tag3.cpp " << endl;
    dvar_vector tmpFF=log(qq);
    for (fi=1;fi<=nfi;fi++)
    {
      fm(fi)(jmin,nage)=logsel(fi)+tmpFF(fi);
    }
    tm.initialize();
    for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
    {                                         // incidents for this period
      tm+=mfexp(fm(fi));
    }
    tm+=mfexp(get_nat_mort_region(ir)(year(ir,ip))(jmin,nage)
      +fraction(ir,ip));
  }
  else
  {
    for (fi=1;fi<=nfi;fi++)
    {
      fm(fi)=-10.0;
    }
    tm=mfexp(get_nat_mort_region(ir)(year(ir,ip))(jmin,nage)
      +fraction(ir,ip));
  }
  nrsurv(it,ir,ip)=mfexp(-tm);
  switch (it)
  {
  case 72:
  case 21:
  case 112:
    break;
  default:
    ;
    _ffpen+=ffpen;
  }
}

void dfjet_jac2(void);

dvar_matrix ADJ_get_jac2(int nfi,int jmin,int nage,dvar_vector _q,
  dvar_matrix _sel,const dvar_vector& _N,dvar_vector& _M,MY_DOUBLE_TYPE rmax)
{
  int i,j,k;
  dmatrix J(1,nfi,1,nfi);
  J.initialize();
  dmatrix F(1,nfi);
  dvector q=value(_q);
  dvector Z(jmin,nage);
  dmatrix sel=value(_sel);
  dvector N=value(_N);
  dvector M=value(_M);
  for (i=1;i<=nfi;i++)
  {
    F(i)=q(i)*sel(i);
  }
  Z=colsum(F)+M;
  dvector S=exp(-Z);
  dvector NZ=elem_div(N,Z);
  for (i=1;i<=nfi;i++)
  {
    for (k=1;k<=nfi;k++)
    {
      for (j=jmin;j<=nage;j++)
      {
        MY_DOUBLE_TYPE tmp0=0;
        tmp0-=F(i,j)/Z(j)*sel(k,j);
        if (i==k)
        {
          tmp0+=sel(k,j);
        }
        MY_DOUBLE_TYPE tmp=tmp0*(1.0-S(j));
        if (Z(j)<=rmax)
        {
          tmp+=F(i,j)*S(j)*sel(k,j);  //baranov
        }
        else
        {
          tmp+=F(i,j)*exp(-rmax)*sel(k,j);  //kludged
        }
        J(i,k)+=tmp*NZ(j);
      }
    }
  }
  dvar_matrix VJ=nograd_assign(J);
//    save_identifier_string("X1");
  const char * str1;
  str1="X1";
  char* strx1=const_cast <char*> (str1);
  save_identifier_string(strx1);
  _sel.save_dvar_matrix_value();
  _sel.save_dvar_matrix_position();
//    save_identifier_string("Y0");
  const char * str2;
  str2="Y0";
  char* strx2=const_cast <char*> (str2);
  save_identifier_string(strx2);
  _q.save_dvar_vector_value();
  _q.save_dvar_vector_position();
  _N.save_dvar_vector_value();
  _N.save_dvar_vector_position();
//    save_identifier_string("Y1");
  const char * str3;
  str3="Y1";
  char* strx3=const_cast <char*> (str3);
  save_identifier_string(strx3);
  _M.save_dvar_vector_value();
  _M.save_dvar_vector_position();
  VJ.save_dvar_matrix_position();
//    save_identifier_string("Z1");
  const char * str4;
  str4="Z1";
  char* strx4=const_cast <char*> (str4);
  save_identifier_string(strx4);
  save_int_value(nfi);
  save_int_value(jmin);
  save_int_value(nage);
  save_double_value(rmax);
//    save_identifier_string("Z2");
  const char * str5;
  str5="Z2";
  char* strx5=const_cast <char*> (str5);
  save_identifier_string(strx5);
  gradient_structure::GRAD_STACK1->
            set_gradient_stack(dfjet_jac2);
  return VJ;
}
static void dfjet_jac2(void)
{
  int i,j,k;
  verify_identifier_string("Z2");
  MY_DOUBLE_TYPE rmax=restore_double_value();
  int nage=restore_int_value();
  int jmin=restore_int_value();
  int nfi=restore_int_value();
  verify_identifier_string("Z1");
  dvar_matrix_position VJpos=restore_dvar_matrix_position();
  dmatrix dfJ=restore_dvar_matrix_derivatives(VJpos);
  dvar_vector_position Mpos=restore_dvar_vector_position();
  dvector M=restore_dvar_vector_value(Mpos);
  dvector dfM(M.indexmin(),M.indexmax());
  dfM.initialize();

  verify_identifier_string("Y1");
  dvar_vector_position Npos=restore_dvar_vector_position();
  dvector N=restore_dvar_vector_value(Npos);
  dvector dfN(N.indexmin(),N.indexmax());
  dfN.initialize();

  dvar_vector_position qpos=restore_dvar_vector_position();
  dvector q=restore_dvar_vector_value(qpos);
  dvector dfq(q.indexmin(),q.indexmax());
  dfq.initialize();

  verify_identifier_string("Y0");

  dvar_matrix_position selpos=restore_dvar_matrix_position();
  dmatrix sel=restore_dvar_matrix_value(selpos);
  int mmin=sel.indexmin();
  int mmax=sel.indexmax();
  ivector sel_indexmin(mmin,mmax);
  ivector sel_indexmax(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    sel_indexmin(i)=sel(i).indexmin();
    sel_indexmax(i)=sel(i).indexmax();
  }

  dmatrix dfF(mmin,mmax,sel_indexmin,sel_indexmax);
  dfF.initialize();

  dmatrix dfsel(mmin,mmax,sel_indexmin,sel_indexmax);
  dfsel.initialize();

  verify_identifier_string("X1");

  dmatrix J(1,nfi,1,nfi);
  J.initialize();
  dmatrix F(1,nfi);
  dvector Z(jmin,nage);
  dvector dfZ(jmin,nage);
  dfZ.initialize();
  for (i=1;i<=nfi;i++)
  {
    F(i)=q(i)*sel(i);
  }
  Z=colsum(F)+M;
  dvector S=exp(-Z);
  dvector dfS(S.indexmin(),S.indexmax());
  dfS.initialize();
  dvector NZ=elem_div(N,Z);
  dvector dfNZ(jmin,nage);
  dfNZ.initialize();
  for (i=1;i<=nfi;i++)
  {
    for (k=1;k<=nfi;k++)
    {
      for (j=jmin;j<=nage;j++)
      {
        MY_DOUBLE_TYPE tmp0=0;
        MY_DOUBLE_TYPE dftmp0=0;
        MY_DOUBLE_TYPE u=F(i,j)/Z(j);
        MY_DOUBLE_TYPE dfu;
        dfu=0.0;
        tmp0-=u*sel(k,j);
        if (i==k)
        {
          tmp0+=sel(k,j);
        }
        MY_DOUBLE_TYPE tmp=tmp0*(1.0-S(j));
        MY_DOUBLE_TYPE dftmp=0;
        if (Z(j)<=rmax)
        {
          tmp+=F(i,j)*S(j)*sel(k,j);  //baranov
        }
        else
        {
          tmp+=F(i,j)*exp(-rmax)*sel(k,j);  //kludged
        }
        J(i,k)+=tmp*NZ(j);
       
        //J(i,k)+=tmp*NZ(j);
        dftmp+=dfJ(i,k)*NZ(j);
        dfNZ(j)+=dfJ(i,k)*tmp;
        if (Z(j)<=rmax)
        {
          //tmp+=F(i,j)*S(j)*sel(k,j);  //baranov
          dfsel(k,j)+=dftmp*F(i,j)*S(j);
          dfF(i,j)+=dftmp*S(j)*sel(k,j);
          dfS(j)+=dftmp*F(i,j)*sel(k,j);
        }
        else
        {
          //tmp+=F(i,j)*exp(-rmax)*sel(k,j);  //kludged
          dfsel(k,j)+=dftmp*F(i,j)*exp(-rmax);
          dfF(i,j)+=dftmp*exp(-rmax)*sel(k,j);
        }
        //double tmp=tmp0*(1.0-S(j));
        dftmp0+=dftmp*(1.0-S(j));
        dfS(j)-=dftmp*tmp0;
        dftmp=0;
        if (i==k)
        {
          //tmp0+=sel(k,j);
          dfsel(k,j)+=dftmp0;
        }
        //tmp0-=u*sel(k,j);
        dfu-=dftmp0*sel(k,j);
        dfsel(k,j)-=dftmp0*u;
        //double u=F(i,j)/Z(j);
        dfF(i,j)+=dfu/Z(j);
        dfZ(j)-=dfu*F(i,j)/square(Z(j));;
        dfu=0.0;
      }
    }
  }
  //cout << "norm2(J)   "  <<  norm2(J) << endl;
 
  //dvector NZ=elem_div(N,Z);
  dfN+=elem_div(dfNZ,Z);
  dfZ-=elem_prod(dfNZ,elem_div(N,square(Z)));
  dfNZ.initialize();
  //dvector S=exp(-Z);
  dfZ-=elem_prod(dfS,exp(-Z));
  dfS.initialize();
  //Z=colsum(F)+M;
  int rowmin=F.indexmin();
  int rowmax=F.indexmax();
  int cmin=F(rowmin).indexmin();
  int cmax=F(rowmin).indexmax();
  for (int ll=rowmin;ll<=rowmax;ll++)
  {
    for (int kk=cmin;kk<=cmax;kk++)
    {
      //Z(kk)+=F(ll,kk);
      dfF(ll,kk)+=dfZ(kk);
    }
  }
  for (int kk=cmin;kk<=cmax;kk++)
  {
    //Z(kk)=M(kk);
    dfM(kk)+=dfZ(kk);
    dfZ(kk)=0.0;
  }
  for (i=1;i<=nfi;i++)
  {
    int mmin=sel(i).indexmin();
    int mmax=sel(i).indexmax();
    for (int j=mmin;j<=mmax;j++)
    {
      //F(i,j)=q(i)*sel(i,j);
      dfq(i)+=dfF(i,j)*sel(i,j);
      dfsel(i,j)+=dfF(i,j)*q(i);
      dfF(i,j)=0.0;
    }
  }

  dfq.save_dvector_derivatives(qpos);
  dfM.save_dvector_derivatives(Mpos);
  dfN.save_dvector_derivatives(Npos);
  dfsel.save_dmatrix_derivatives(selpos);

  // *********************************************************
  // *********************************************************
}
//void dvar_fish_stock_history::do_fish_mort_intermediate_calcs(void)
void dvar_fish_stock_history::do_tag_loss_fish_mort_intermediate_calcs(void)
{
  const MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
  dvar_matrix * ptmp=0;
  dvar_matrix * ptmp1=0;
      
  //if (age_flags(92) !=3)
  switch(age_flags(92))
  {
  case 0:
  case 2:
    {
      for (int ir=1;ir<=num_regions;ir++)
      {
        for (int ip=1;ip<=num_real_fish_periods(ir);ip++)  
        {
          int yr=year(ir,ip);
          if (age_flags(125) &&year(ir,ip)>1) 
            break;
          dvar_matrix* pfm;
          dvar_vector* ptotmort;
          dvar_vector* psurv;
          if (af170q0==0)
          {
            ptotmort=&(tot_mort(ir,ip));
            psurv=&(survival(ir,ip));
            ptmp1=&(fish_mort_calcs(ir,ip));
            ptmp=&(tag_fish_mort_calcs(1,ir,ip));
            pfm=&(fish_mort(ir,ip));
          }
          else
          {
            ptotmort=&(tot_mort_q0(ir,ip));
            psurv=&(survival_q0(ir,ip));
            ptmp=&(fish_mort_calcs_q0(ir,ip));
            pfm=&(fish_mort_q0(ir,ip));
          }

          if (age_flags(115)==0)
          {
            const dvar_vector& tmp1=log(1.e-10+ *ptotmort)
              +tagmort(1);
            const dvar_vector& tmp2=log(one_plus- *psurv)
              -tagmort(1);
            const dvar_vector& tmp3=tmp2-tmp1;
            for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
            {        
              (*ptmp)(fi)=(*pfm)(fi)+tmp3;
            }
          }
          else   // ss2 parameterization
          {
            for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
            {        
              (*ptmp)(fi)=(*pfm)(fi)-0.5*e_nat_mort_by_period(ir,ip);
            }
          }
        }
      } 
    }
    break;
  case 4:
    {
      for (int ir=1;ir<=num_regions;ir++)
      {
        for (int ip=1;ip<=num_real_fish_periods(ir);ip++)  
        {
          if (age_flags(125) &&year(ir,ip)>1) 
            break;
          if (af170q0==0)
          {
            ptmp=&(fish_mort_calcs(ir,ip));
          }
          else
          {
            ptmp=&(fish_mort_calcs_q0(ir,ip));
          }
          dvar_vector em=e_nat_mort_by_period(ir,ip);
          const dvar_vector& zz=tot_mort(ir,ip)-em;
          int mmin=zz.indexmin();
          int mmax=zz.indexmax();
          for (int j=mmin;j<=mmax;j++)
          {
            if (zz(j)<0.0)
                cout<< "HERE 0" << endl;
          }
           
          const dvar_vector& tmp1=log(1.e-10+zz);
          const dvar_vector& tmp2=log(one_plus-exp(-zz));
          const dvar_vector& tmp3=tmp2-tmp1;
          dvar_matrix* pfm;
          if (af170q0)
          {
            pfm=&(fish_mort_q0(ir,ip));
          }
          else
          {
            pfm=&(fish_mort(ir,ip));
          }
          for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
          {        
            (*ptmp)(fi)=(*pfm)(fi)+tmp3-0.5*em;
          }
        }
      } 
    }
    break;
  case 3:
    {
      for (int ir=1;ir<=num_regions;ir++)
      {
        for (int ip=1;ip<=num_real_fish_periods(ir);ip++)  
        {
          if (age_flags(125) &&year(ir,ip)>1) 
            break;
          if (af170q0==0)
          {
            ptmp=&(fish_mort_calcs(ir,ip));
          }
          else
          {
            ptmp=&(fish_mort_calcs_q0(ir,ip));
          }
          dvar_matrix* pfm;
          if (af170q0)
          {
            pfm=&(fish_mort_q0(ir,ip));
          }
          else
          {
            pfm=&(fish_mort(ir,ip));
          }
          for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
          {        
            (*ptmp)(fi)=(*pfm)(fi)-0.5*e_nat_mort_by_period(ir,ip);
          }
        }
      } 
    }
    break;
  default: 
    cerr << "Illegal value of " << age_flags(92) << " for age_flags(92)" 
         << endl;
    ad_exit(1);
  }
}
