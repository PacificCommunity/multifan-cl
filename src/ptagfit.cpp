/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
//#include <f:\linad99\fvar.hpp>
#include "all.hpp"
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

dvariable dvar_len_fish_stock_history::fit_pooled_tag_returns(void)
{
  //ofstream ofsyy("pooledcont");
  int fi;
  ivector group_flags32=column(fish_flags,32);
  int gsum32=sum(group_flags32);
  int gmax32=Max(group_flags32);
  ivector gp_fish32(1,gmax32);
  ivector iflag(1,gmax32);
  dvar_matrix pooledobsgroupedcatch;
  dvar_matrix pooledgroupedcatch;
  if (gmax32)
  {
    gp_fish32.initialize();
    for (fi=1;fi<=num_fisheries;fi++)
      gp_fish32(group_flags32(fi))=fi;
    pooledobsgroupedcatch.allocate(1,gmax32,1,nage);
    pooledgroupedcatch.allocate(1,gmax32,1,nage);
  }
  else
  {
    pooledobsgroupedcatch.allocate(1,num_fisheries,1,nage);
    pooledgroupedcatch.allocate(1,num_fisheries,1,nage);
  }

  dvariable f=0.0;
  dvariable gp_pen=0.0;
  dvar_vector& rep_rate1=fish_pars(3);
  if (!gsum32) // no grouping
  {
    const int max_index=10000;
    dmatrix tmp_tag_like;
    if (ppstf)
    {
      tmp_tag_like.allocate(1,num_fisheries,1,max_index);
      tmp_tag_like.initialize();
    }
    ivector uicount(1,num_fisheries);
    uicount.initialize();

    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int ip=minttp(ir)+1;ip<=num_real_fish_periods(ir);ip++)
      {
        for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {
          if (fishery_projection_flag(ir,ip,fi)==0)
          {
            for (int j=1;j<=nage;j++)
            {
              if (pooled_tagcatch(ir,ip,fi,j)<0)
              {
                cerr << pooled_tagcatch(ir,ip,fi,j) << endl;
              }
            }
            int pp1=parent(ir,ip,fi);
            const dvar_vector& xrtc=rep_rate(ir,ip,fi)*pooled_tagcatch(ir,ip,fi);
            dvar_vector& rtc=(dvar_vector&) xrtc;
            int pf111=parest_flags(111);
            switch (pf111)
	      {
            case 0:
              f+=sum(elem_div(square(pooledobstagcatch(ir,ip,fi)-rtc),
                 .01+rtc));
              break;
            case 1:
              f-=sum(log(exp(elem_div(
                square(pooledobstagcatch(ir,ip,fi)-rtc),
                -.01-rtc ))+0.01));
              break;
            case 2:
	        {
	          dvar_vector lrtc=log(1.e-10+rtc);
                  f+=sum(rtc) - pooledobstagcatch(ir,ip,fi)*lrtc;
		  for (int i=1;i<=nage;i++)
#if !defined(NO_MY_DOUBLE_TYPE)
                  f+=gammln(pooledobstagcatch(ir,ip,fi,i)+1.0L);
#else
                  f+=gammln(pooledobstagcatch(ir,ip,fi,i)+1.0);
#endif
              }
              break;
            case 3:
              {
                dvariable a=fish_pars(4,pp1)+50.0001;
                rtc+=1.e-10;
                int ns=rtc.indexmax()-rtc.indexmin()+1; 
                dvar_vector ap=log(a+rtc);
                f+=sum(a*ap)+pooledobstagcatch(ir,ip,fi)*ap;
                f-=ns*a*log(a)+pooledobstagcatch(ir,ip,fi)*log(rtc);
                f-=sum(gammln(a+pooledobstagcatch(ir,ip,fi)));
#if !defined(NO_MY_DOUBLE_TYPE)
                f+=sum(gammln(pooledobstagcatch(ir,ip,fi)+1.0L));
#else
                f+=sum(gammln(pooledobstagcatch(ir,ip,fi)+1.0));
#endif
                f+=ns*gammln(a);
              }
              break;
            case 4:
              {
                rtc+=1.e-10;
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
                int ns=rtc.indexmax()-rtc.indexmin()+1; 
                dvar_vector ap=log(a+rtc);
                f+=a*ap+pooledobstagcatch(ir,ip,fi)*ap;
                f-=a*log(a)+pooledobstagcatch(ir,ip,fi)*log(rtc);
                f-=sum(gammln(a+pooledobstagcatch(ir,ip,fi)));
#if !defined(NO_MY_DOUBLE_TYPE)
                f+=sum(gammln(pooledobstagcatch(ir,ip,fi)+1.0L));
#else
                f+=sum(gammln(pooledobstagcatch(ir,ip,fi)+1.0));
#endif
                f+=sum(gammln(a));
              }
              break;

            case 5:
            case 6:
            case 7:
            case 8:
              {
//                dvariable tau=1.0+exp(fish_pars(4,pp1));
                dvariable tau=exp(fish_pars(30,pp1));   //NMD_18Sep_19
                MY_DOUBLE_TYPE censor=1.1;
                dvar_vector pred =rep_rate(ir,ip,fi)*pooled_tagcatch(ir,ip,fi);
                if (parest_flags(325))
                {
                  censor=parest_flags(325)/100.0;
                }
                dvar_vector vtmp=censored_gamma_or_lognormal(tau,
                  pooledobstagcatch(ir,ip,fi), // pooledobsgroupedcatch(fi),
                  1.e-2+pred, censor,pf111);
                f-=sum(log(1.e-10+vtmp));

                dvariable tmpf;       //NMD_30Sep2019
                tmpf=-sum(log(1.e-10+vtmp));
                uicount(pp1)++;
                if (uicount(pp1)>max_index)
                {
                  cerr << "need to increase max_index" << endl;
                  ad_exit(1);
                }
                tmp_tag_like(pp1,uicount(pp1))=value(tmpf);
              }       //NMD_30Sep2019
              break;

	      default:
	        cerr << "Illegal value for parest_flags(111) value = "
	           <<  parest_flags(111)  << endl;
            }
          }
	}
      }
    }
    //NMD_30Sep2019
    if (ppstf)
    {
      if (allocated(ppstf->ungrouped_pooled_tag_like))
        ppstf->ungrouped_pooled_tag_like.deallocate();

      ppstf->ungrouped_pooled_tag_like.allocate(1,num_fisheries);
      for (int fi=1;fi<=num_fisheries;fi++)
      {
        ppstf->ungrouped_pooled_tag_like(fi)=
          tmp_tag_like(fi)(1,uicount(fi));
      }
    }
    //NMD_30Sep2019

  }
  else  // grouping
  {
    int icount=0;
    ivector gicount(1,gmax32);    //NMD_30Sep2019
    gicount.initialize();
    if (ppstf && allocated(ppstf->grouped_pooled_tag_like))
      ppstf->grouped_pooled_tag_like.deallocate();

//    const int max_index=1000;
    const int max_index=2000;  //NMD_19jan2023
    dvar_vector f_by_tag(1,10000);
    f_by_tag.initialize();      
    dmatrix tmp_grouped_tag_like;
    if (ppstf)
    {
      tmp_grouped_tag_like.allocate(1,gmax32,1,max_index);
      tmp_grouped_tag_like.initialize();
    }

    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int ip=minttp(ir)+1;ip<=num_real_fish_periods(ir);ip++)
      {
        pooledobsgroupedcatch.initialize();
        pooledgroupedcatch.initialize();
        iflag.initialize();
        ivector& pi=parent(ir,ip);
        for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {
          int pp1=pi(fi);
          int pp=group_flags32(pp1);
          const dvar_vector& rtc=rep_rate(ir,ip,fi)*pooled_tagcatch(ir,ip,fi);
          iflag(pp)=1;
          pooledobsgroupedcatch(pp)+=pooledobstagcatch(ir,ip,fi);
          pooledgroupedcatch(pp)+=rtc;
        }

        for (fi=1;fi<=gmax32;fi++)
        {
          if (grouped_fishery_projection_flag(ir,ip,fi)==0)
          {
            if(iflag(fi))
            {
              int pf111=parest_flags(111);
              switch (pf111)
              {
              case 0:
                f+=sum(elem_div(
                  square(pooledobsgroupedcatch(fi)-pooledgroupedcatch(fi)),
                 .01+pooledgroupedcatch(fi) ));
                break;
              case 1:
                f-=sum(log(exp(elem_div(
                  square(pooledobsgroupedcatch(fi)-pooledgroupedcatch(fi)),
                  -.01-pooledgroupedcatch(fi) ))+0.01));
                break;
              case 2:
                {
                  dvar_vector lrtc=log(1.e-10+pooledgroupedcatch(fi));
                  f+=sum(pooledgroupedcatch(fi))
                    - pooledobsgroupedcatch(fi)*lrtc;
                  for (int i=1;i<=nage;i++)
#if !defined(NO_MY_DOUBLE_TYPE)
                  f+=gammln(pooledobsgroupedcatch(fi,i)+1.0L);
#else
                  f+=gammln(pooledobsgroupedcatch(fi,i)+1.0);
#endif
                }
                break;
              case 3:
                {
                  dvar_vector rtc=1.e-10+pooledgroupedcatch(fi);
                  dvariable a=fish_pars(4,gp_fish32(fi))+50.0001;
                  int ns=rtc.indexmax()-rtc.indexmin()+1; 
                  dvar_vector ap=log(a+rtc);
                  f+=sum(a*ap)+pooledobsgroupedcatch(fi)*ap;
                  f-=ns*a*log(a)+pooledobsgroupedcatch(fi)*log(rtc);
                  f-=sum(gammln(a+pooledobsgroupedcatch(fi)));
#if !defined(NO_MY_DOUBLE_TYPE)
                  f+=sum(gammln(pooledobsgroupedcatch(fi)+1.0L));
#else
                  f+=sum(gammln(pooledobsgroupedcatch(fi)+1.0));
#endif
                  f+=ns*gammln(a);
                }
                break;
              case 4:
                {
                  // !! Dave dvar_vector rtc=1.e-10+pooledgroupedcatch(fi);
                  dvar_vector rtc=1.e-10+pooledgroupedcatch(fi);
                  dvar_vector a;
                  if (parest_flags(305)==0)
                  {
                    a=(fish_pars(4,gp_fish32(fi))+50.0001)*rtc;
                  }
                  else
                  {
                    dvariable tau=1.0+exp(fish_pars(4,gp_fish32(fi)));
#if !defined(NO_MY_DOUBLE_TYPE)
                    a=rtc/(1.e-20+(tau-1.0L));
#else
                    a=rtc/(1.e-20+(tau-1.0));
#endif
                  }
                  dvariable f1=0.0;
                  ++icount;
       	          //NMD_14nov2023
                  gicount(fi)++;
		  if (gicount(fi)>max_index)
                  {
                    cerr << "need to increase max_index" << endl;
                    ad_exit(1);
                  }		  
                  int ns=rtc.indexmax()-rtc.indexmin()+1; 
                  dvar_vector ap=log(a+rtc);
                  f1+=a*ap+pooledobsgroupedcatch(fi)*ap;
                  f1-=a*log(a)+pooledobsgroupedcatch(fi)*log(rtc);
                  f1-=sum(gammln(a+pooledobsgroupedcatch(fi)));
#if !defined(NO_MY_DOUBLE_TYPE)
                  f1+=sum(gammln(pooledobsgroupedcatch(fi)+1.0L));
#else
                  f1+=sum(gammln(pooledobsgroupedcatch(fi)+1.0));
#endif
                  f1+=sum(gammln(a));
                  //ofsyy << icount << " " << f1 << " " << f << endl;
                  f+=f1;
                  f_by_tag(icount)=f1;
                  if (ppstf)   //NMD_14nov2023
                  {
                    tmp_grouped_tag_like(fi,gicount(fi))=
                      value(f_by_tag(icount));
                  }
                }
                break;
              case 5:
              case 6:
              case 7:
              case 8:
                {
                  // !! Dave dvar_vector rtc=1.e-10+pooledgroupedcatch(fi);
                  dvar_vector rtc=1.e-10+pooledgroupedcatch(fi);
                  dvariable tau=exp(fish_pars(30,gp_fish32(fi)));
                  MY_DOUBLE_TYPE censor=1.1;
                  if (parest_flags(325))
                  {
                    censor=parest_flags(325)/100.0;
                  }
                  icount++;
                  gicount(fi)++;
		  if (gicount(fi)>max_index)
                  {
                    cerr << "need to increase max_index" << endl;
                    ad_exit(1);
                  }		  
                  dvar_vector vtmp=censored_gamma_or_lognormal(tau,
                    pooledobsgroupedcatch(fi),1.e-2+pooledgroupedcatch(fi),
                    censor,pf111);
                  //NMD_30Sep2019
                  f_by_tag(icount)-=sum(log(1.e-10+vtmp));
      	          f+=f_by_tag(icount);

                  if (ppstf)   //NMD_30Sep2019
                  {
                    tmp_grouped_tag_like(fi,gicount(fi))=
                      value(f_by_tag(icount));
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
    //NMD_30Sep2019
    if (ppstf)
    {
      if (allocated(ppstf->grouped_pooled_tag_like))
        ppstf->grouped_pooled_tag_like.deallocate();

      ppstf->grouped_pooled_tag_like.allocate(1,num_fisheries);
      for (int fi=1;fi<=gmax32;fi++)
      {
        ppstf->grouped_pooled_tag_like(fi)=
          tmp_grouped_tag_like(fi)(1,gicount(fi));
      }
    }
    //NMD_30Sep2019

  }
  f+=gp_pen;
  //check_derivative_values("Bpos0");
  return f;
}

dvariable dvar_len_fish_stock_history::fit_pooled_tag_returns_mix(void)
{
  int fi;
  ivector group_flags32=column(fish_flags,32);
  int gsum32=sum(group_flags32);
  int gmax32=Max(group_flags32);
  ivector gp_fish32(1,gmax32);
  ivector iflag(1,gmax32);
  dvar_matrix pooledobsgroupedcatch;
  dvar_matrix pooledgroupedcatch;
  if (gmax32)
  {
    gp_fish32.initialize();
    for (fi=1;fi<=num_fisheries;fi++)
      gp_fish32(group_flags32(fi))=fi;
    pooledobsgroupedcatch.allocate(1,gmax32,1,nage);
    pooledgroupedcatch.allocate(1,gmax32,1,nage);
  }
  else
  {
    pooledobsgroupedcatch.allocate(1,num_fisheries,1,nage);
    pooledgroupedcatch.allocate(1,num_fisheries,1,nage);
  }

  dvariable f=0.0;
  dvariable gp_pen=0.0;
  dvar_vector& rep_rate1=fish_pars(3);
  if (!gsum32) // no grouping
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int ip=minttp(ir)+1;ip<=num_real_fish_periods(ir);ip++)
      {
        for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {
          for (int j=1;j<=nage;j++)
          {
            if (pooled_tagcatch(ir,ip,fi,j)<0)
            {
              cerr << pooled_tagcatch(ir,ip,fi,j) << endl;
            }
            int pp1=parent(ir,ip,fi);
            const dvar_vector& xrtc=rep_rate(ir,ip,fi)*pooled_tagcatch(ir,ip,fi);
            dvar_vector& rtc=(dvar_vector&) xrtc;

            int pf111=parest_flags(111);
            switch (pf111)
	    {
            case 0:
              f+=sum(elem_div(square(pooledobstagcatch(ir,ip,fi)-rtc),
                 .01+rtc));
              break;
            case 1:
              f-=sum(log(exp(elem_div(
                square(pooledobstagcatch(ir,ip,fi)-rtc),
                -.01-rtc ))+0.01));
              break;
            case 2:
	      {
	        dvar_vector lrtc=log(1.e-10+rtc);
                f+=sum(rtc) - pooledobstagcatch(ir,ip,fi)*lrtc;
		for (int i=1;i<=nage;i++)
#if !defined(NO_MY_DOUBLE_TYPE)
                  f+=gammln(pooledobstagcatch(ir,ip,fi,i)+1.0L);
#else
                  f+=gammln(pooledobstagcatch(ir,ip,fi,i)+1.0);
#endif
              }
              break;
            case 3:
              {
                dvariable a=fish_pars(4,pp1)+50.0001;
                rtc+=1.e-10;
                int ns=rtc.indexmax()-rtc.indexmin()+1; 
                dvar_vector ap=log(a+rtc);
                f+=sum(a*ap)+pooledobstagcatch(ir,ip,fi)*ap;
                f-=ns*a*log(a)+pooledobstagcatch(ir,ip,fi)*log(rtc);
                f-=sum(gammln(a+pooledobstagcatch(ir,ip,fi)));
#if !defined(NO_MY_DOUBLE_TYPE)
                f+=sum(gammln(pooledobstagcatch(ir,ip,fi)+1.0L));
#else
                f+=sum(gammln(pooledobstagcatch(ir,ip,fi)+1.0));
#endif
                f+=ns*gammln(a);
              }
              break;
            case 4:
              {
                rtc+=1.e-10;
                int pp1=parent(ir,ip,fi);
                dvar_vector a=(fish_pars(4,pp1)+50.0001)*rtc;
                //dvar_vector a=(fish_pars(4,pp1)+50.0001)*sqrt(rtc);
                dvar_vector zerop=fish_pars(5,pp1)/
                  (1.0+fish_pars(6,pp1)*rtc);
                dvar_vector zerop1=1-zerop;
                rtc=elem_div(rtc,zerop1);
                int ns=rtc.indexmax()-rtc.indexmin()+1; 
                dvar_vector ap=log(a+rtc);
               
                dvar_vector f1=elem_prod(a,ap)
                  +elem_prod(pooledobstagcatch(ir,ip,fi),ap);
                f1-=elem_prod(a,log(a))
                  +elem_prod(pooledobstagcatch(ir,ip,fi),log(rtc));
                f1-=gammln(a+pooledobstagcatch(ir,ip,fi));
#if !defined(NO_MY_DOUBLE_TYPE)
                f1+=gammln(pooledobstagcatch(ir,ip,fi)+1.0L);
#else
                f1+=gammln(pooledobstagcatch(ir,ip,fi)+1.0);
#endif
                f1+=gammln(a);
                dvar_vector lzerop1=log(1.0-zerop);
		for (int j=1;j<=nage;j++)
                {
                  if (pooledobstagcatch(ir,ip,fi,j)>0)
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
  else  // grouping
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int ip=minttp(ir)+1;ip<=num_real_fish_periods(ir);ip++)
      {
        pooledobsgroupedcatch.initialize();
        pooledgroupedcatch.initialize();
        iflag.initialize();
        ivector& pi=parent(ir,ip);
        for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {
          int pp1=pi(fi);
          int pp=group_flags32(pp1);
          const dvar_vector& rtc=rep_rate(ir,ip,fi)*pooled_tagcatch(ir,ip,fi);
          iflag(pp)=1;
          pooledobsgroupedcatch(pp)+=pooledobstagcatch(ir,ip,fi);
          pooledgroupedcatch(pp)+=rtc;
        }

        for (fi=1;fi<=gmax32;fi++)
        {
          if(iflag(fi))
          {
            switch (parest_flags(111))
            {
            case 0:
              f+=sum(elem_div(
                square(pooledobsgroupedcatch(fi)-pooledgroupedcatch(fi)),
               .01+pooledgroupedcatch(fi) ));
              break;
            case 1:
              f-=sum(log(exp(elem_div(
                square(pooledobsgroupedcatch(fi)-pooledgroupedcatch(fi)),
                -.01-pooledgroupedcatch(fi) ))+0.01));
              break;
            case 2:
	      {
                dvar_vector lrtc=log(1.e-10+pooledgroupedcatch(fi));
                f+=sum(pooledgroupedcatch(fi))
	          - pooledobsgroupedcatch(fi)*lrtc;
	        for (int i=1;i<=nage;i++)
#if !defined(NO_MY_DOUBLE_TYPE)
                f+=gammln(pooledobsgroupedcatch(fi,i)+1.0L);
#else
                f+=gammln(pooledobsgroupedcatch(fi,i)+1.0);
#endif
              }
              break;
            case 3:
              {
                dvar_vector rtc=1.e-10+pooledgroupedcatch(fi);
                dvariable a=fish_pars(4,gp_fish32(fi))+50.0001;
                int ns=rtc.indexmax()-rtc.indexmin()+1; 
                dvar_vector ap=log(a+rtc);
                f+=sum(a*ap)+pooledobsgroupedcatch(fi)*ap;
                f-=ns*a*log(a)+pooledobsgroupedcatch(fi)*log(rtc);
                f-=sum(gammln(a+pooledobsgroupedcatch(fi)));
#if !defined(NO_MY_DOUBLE_TYPE)
                f+=sum(gammln(pooledobsgroupedcatch(fi)+1.0L));
#else
                f+=sum(gammln(pooledobsgroupedcatch(fi)+1.0));
#endif
                f+=ns*gammln(a);
              }
              break;
            case 4:
              {
                dvar_vector rtc=1.e-10+pooledgroupedcatch(fi);
                dvar_vector a=(fish_pars(4,gp_fish32(fi))+50.0001)*rtc;
                //dvar_vector a=(fish_pars(4,gp_fish32(fi))+50.0001)*sqrt(rtc);
                dvar_vector zerop=fish_pars(5,gp_fish32(fi))/
                  (1.0+fish_pars(6,gp_fish32(fi))*rtc);
                dvar_vector zerop1=1-zerop;
                rtc=elem_div(rtc,zerop1);
                int ns=rtc.indexmax()-rtc.indexmin()+1; 
                dvar_vector ap=log(a+rtc);
               
                dvar_vector f1=elem_prod(a,ap)
                  +elem_prod(pooledobsgroupedcatch(fi),ap);
                f1-=elem_prod(a,log(a))
                  +elem_prod(pooledobsgroupedcatch(fi),log(rtc));
                f1-=gammln(a+pooledobsgroupedcatch(fi));
#if !defined(NO_MY_DOUBLE_TYPE)
                f1+=gammln(pooledobsgroupedcatch(fi)+1.0L);
#else
                f1+=gammln(pooledobsgroupedcatch(fi)+1.0);
#endif
                f1+=gammln(a);
                dvar_vector lzerop1=log(1.0-zerop);
		for (int j=1;j<=nage;j++)
                {
                  if (pooledobsgroupedcatch(fi,j)>0)
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
  f+=gp_pen;
  //check_derivative_values("Bpos0");
  return f;
}
