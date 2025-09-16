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
/*
dvariable dvar_len_fish_stock_history::fit_tag_returns(void)
{
  //check_derivative_values("Apos1");
  ivector group_flags=column(fish_flags,44);
  int gsum44=sum(group_flags);
  int gmax44=Max(group_flags);
  dvariable f=0.0;
  dvariable gp_pen=0.0;
  int fi;
  dvar_vector& rep_rate1=fish_pars(3);
  ivector gp_fish44(1,gmax44);
  dvar_matrix obsgroupedcatch;
  dvar_matrix groupedcatch;
  ivector iflag(1,gmax44);
  if (gmax44)
  {
    gp_fish44.initialize();
    for (fi=1;fi<=num_fisheries;fi++)
      gp_fish44(group_flags(fi))=fi;
    obsgroupedcatch.allocate(1,gmax44,1,nage);
    groupedcatch.allocate(1,gmax44,1,nage);
  }
  else
  {
    obsgroupedcatch.allocate(1,num_fisheries,1,nage);
    groupedcatch.allocate(1,num_fisheries,1,nage);
  }

  for (int it=1;it<=num_tag_releases;it++)
  {
    cout << "it = " << it << endl;
    for (int ir=1;ir<=num_regions;ir++)
    {
      cout << "ir = " << ir << endl;
      int ub;
      if (!age_flags(96))
        ub=num_fish_periods(ir);
      else
        ub=terminal_tag_period(it,ir);
      for (int ip=initial_tag_period(it,ir);ip<=ub;ip++)
      {
        cout << "ip = " << ip << endl;
        for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {
          cout << "fi = " << fi << endl;
          if (!gsum44)
          {
            // don't group fisheries
            for (int j=1;j<=nage;j++)
            {
              if (tagcatch(it,ir,ip,fi,j)<0)
              {
                cerr << tagcatch(it,ir,ip,fi,j) << endl;
              }
            }
            int pp1=parent(ir,ip,fi);
            if (rep_rate1(pp1) != rep_rate(ir,ip,fi))
            {
            //  cout <<  rep_rate1(pp1) << " " <<  rep_rate(ir,ip,fi) << endl;
            }
            //  dvar_vector& rtc=rep_rate(pp1)*tagcatch(it,ir,ip,fi);
            dvar_vector& rtc=rep_rate(ir,ip,fi)*tagcatch(it,ir,ip,fi);
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
	        dvar_vector lrtc=log(1.e-10+rtc);
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
                dvariable a=fish_pars(4,fi)+50.0001;
                rtc+=1.e-10;
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
                rtc+=1.e-10;
                //dvar_vector a=(fish_pars(4,fi)+50.0001)*sqrt(rtc);
                dvar_vector a=(fish_pars(4,fi)+50.0001)*rtc;
                int ns=rtc.indexmax()-rtc.indexmin()+1; 
                dvar_vector ap=log(a+rtc);
                f+=a*ap+obstagcatch(it,ir,ip,fi)*ap;
                f-=a*log(a)+obstagcatch(it,ir,ip,fi)*log(rtc);
                f-=sum(gammln(a+obstagcatch(it,ir,ip,fi)));
#if !defined(NO_MY_DOUBLE_TYPE)
                f+=sum(gammln(obstagcatch(it,ir,ip,fi)+1.0L));
#else
                f+=sum(gammln(obstagcatch(it,ir,ip,fi)+1.0));
#endif
                f+=sum(gammln(a));
              }

              break;
	    default:
	      cerr << "Illegal value for parest_flags(111) value = "
	           <<  parest_flags(111)  << endl;
	    }
          }
          else
          {
            obsgroupedcatch.initialize();
            groupedcatch.initialize();
            iflag.initialize();
            ivector& pi=parent(ir,ip);
            for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
            {
              int pp1=pi(fi);
              int pp=group_flags(pp1);
              dvar_vector& rtc=rep_rate(ir,ip,fi)*tagcatch(it,ir,ip,fi);
              iflag(pp)=1;
              obsgroupedcatch(pp)+=obstagcatch(it,ir,ip,fi); 
              groupedcatch(pp)+=rtc;
            }

            for (fi=1;fi<=gmax44;fi++)
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
                    dvar_vector lrtc=log(1.e-2+groupedcatch(fi));
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
                    dvar_vector rtc=1.e-2+groupedcatch(fi);
                    dvariable a=fish_pars(4,gp_fish44(fi))+50.0001;
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
                    // !!! Dave dvar_vector rtc=1.e-10+groupedcatch(fi);
                    dvar_vector obsgroupedcatch1=obsgroupedcatch(fi);
                    dvar_vector rtc=1.e-2+groupedcatch(fi);
                    //dvar_vector a=(fish_pars(4,gp_fish44(fi))+50.0001)*sqrt(rtc);
                    dvar_vector a=(fish_pars(4,gp_fish44(fi))+50.0001)*rtc;
                    int ns=rtc.indexmax()-rtc.indexmin()+1; 
                    dvar_vector ap=log(a+rtc);
                    f+=a*ap+obsgroupedcatch1*ap;
                    f-=a*log(a)+obsgroupedcatch1*log(rtc);
                    f-=sum(gammln(a+obsgroupedcatch1));
#if !defined(NO_MY_DOUBLE_TYPE)
                    f+=sum(gammln(obsgroupedcatch1+1.0L));
#else
                    f+=sum(gammln(obsgroupedcatch1+1.0));
#endif
                    f+=sum(gammln(a));
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



  dvar_vector& a=fish_pars(4)+50.0001;
  //cout << "a = " << a << endl;
  dvar_vector& q=fish_pars(5);
  //cout << "q = " << q << endl;
  if (sum(column(fish_flags,34)))
    gp_pen=grouped_tag_reporting_rate_penalty();

  if (age_flags(96))
  {
    f+=fit_pooled_tag_returns();
  }
  f+=gp_pen;
  cout << "Tag data    " << setfixed() << setprecision(2) << setw(15) << f << endl;
  //check_derivative_values("Apos0");
  return f;
}
xxx`
dvariable dvar_len_fish_stock_history::fit_pooled_tag_returns(void)
{
  //check_derivative_values("Bpos1");
  int fi;
  ivector group_flags=column(fish_flags,44);
  int gsum44=sum(group_flags);
  int gmax44=Max(group_flags);
  ivector gp_fish44(1,gmax44);
  ivector iflag(1,gmax44);
  dvar_matrix pooledobsgroupedcatch;
  dvar_matrix pooledgroupedcatch;
  if (gmax44)
  {
    gp_fish44.initialize();
    for (fi=1;fi<=num_fisheries;fi++)
      gp_fish44(group_flags(fi))=fi;
    pooledobsgroupedcatch.allocate(1,gmax44,1,nage);
    pooledgroupedcatch.allocate(1,gmax44,1,nage);
  }
  else
  {
    pooledobsgroupedcatch.allocate(1,num_fisheries,1,nage);
    pooledgroupedcatch.allocate(1,num_fisheries,1,nage);
  }

  dvariable f=0.0;
  dvariable gp_pen=0.0;
  dvar_vector& rep_rate1=fish_pars(3);
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=minttp(ir)+1;ip<=num_fish_periods(ir);ip++)
    {
      for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {
        if (!gsum44)
        {
          // don't group fisheries
          for (int j=1;j<=nage;j++)
          {
            if (pooled_tagcatch(ir,ip,fi,j)<0)
            {
              cerr << pooled_tagcatch(ir,ip,fi,j) << endl;
            }
            int pp1=parent(ir,ip,fi);
            dvar_vector& rtc=rep_rate(ir,ip,fi)*pooled_tagcatch(ir,ip,fi);
            switch (parest_flags(111))
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
	        dvar_vector lrtc=log(1.e-2+rtc);
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
                dvariable a=fish_pars(4,fi)+50.0001;
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
                dvar_vector a=(fish_pars(4,fi)+50.0001)*rtc;
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

	    default:
	      cerr << "Illegal value for parest_flags(111) value = "
	           <<  parest_flags(111)  << endl;
            }
          }
	}
        else
        {
          pooledobsgroupedcatch.initialize();
          pooledgroupedcatch.initialize();
          iflag.initialize();
          ivector& pi=parent(ir,ip);
          for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
          {
            int pp1=pi(fi);
            int pp=group_flags(pp1);
            dvar_vector& rtc=rep_rate(ir,ip,fi)*pooled_tagcatch(ir,ip,fi);
            iflag(pp)=1;
            pooledobsgroupedcatch(pp)+=pooledobstagcatch(ir,ip,fi);
            pooledgroupedcatch(pp)+=rtc;
          }
          for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
          {
            for (int fi1=1;fi1<=fi;fi1++)
            {
              if (group_flags(pi(fi))== group_flags(pi(fi1)))
                gp_pen+=100.*square(log(
                  (1.e-5+rep_rate(ir,ip,fi))/
                  (1.e-5+rep_rate(ir,ip,fi1)) ));
             }
          }

          for (fi=1;fi<=gmax44;fi++)
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
                  dvariable a=fish_pars(4,fi)+50.0001;
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
                  //dvar_vector a=(fish_pars(4,gp_fish44(fi))+50.0001)*sqrt(rtc);
                  dvar_vector a=(fish_pars(4,gp_fish44(fi))+50.0001)*rtc;
                  int ns=rtc.indexmax()-rtc.indexmin()+1; 
                  dvar_vector ap=log(a+rtc);
                  f+=a*ap+pooledobsgroupedcatch(fi)*ap;
                  f-=a*log(a)+pooledobsgroupedcatch(fi)*log(rtc);
                  f-=sum(gammln(a+pooledobsgroupedcatch(fi)));
#if !defined(NO_MY_DOUBLE_TYPE)
                  f+=sum(gammln(pooledobsgroupedcatch(fi)+1.0L));
#else
                  f+=sum(gammln(pooledobsgroupedcatch(fi)+1.0));
#endif
                  f+=sum(gammln(a));
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
  f+=gp_pen;
  //check_derivative_values("Bpos0");
  return f;
}
*/
void dvar_fish_stock_history::rep_rate_devs_calc(void)
{
  rep_rate.initialize();
  for (int i=1;i<=num_fisheries;i++)
  {
    int rr=realization_region(i,itind(i));
    int rp=realization_period(i,itind(i));
    int ri=realization_incident(i,itind(i));
    rep_rate(rr,rp,ri)=fish_pars(3,i);
    for (int nt=itind(i)+1;nt<=num_fish_times(i);nt++)
    {
      int rr=realization_region(i,nt);
      int rr1=realization_region(i,nt-1);
      int rp=realization_period(i,nt);
      int rp1=realization_period(i,nt-1);
      rep_rate(rr,rp,realization_incident(i,nt))
        =rep_rate(rr1,rp1,realization_incident(i,nt-1))
          * mfexp(rep_dev_coffs(i,nt));
    }
  }
}

dvariable dvar_len_fish_stock_history::fit_tag_returns_sqrt(void)
{
  ivector group_flags=column(fish_flags,44);
  int gsum44=sum(group_flags);
  dvariable f=0.0;
  int fi;
  dvar_vector& rep_rate1=fish_pars(3);
  for (int it=1;it<=num_tag_releases;it++)
  {
    for (int ir=1;ir<=num_regions;ir++)
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
          if (!gsum44)
          {
            // don't group fisheries
            for (int j=1;j<=nage;j++)
            {
              if (tagcatch(it,ir,ip,fi,j)<0)
              {
                cerr << tagcatch(it,ir,ip,fi,j) << endl;
              }
            }
            int pp1=parent(ir,ip,fi);
            if (rep_rate1(pp1) != rep_rate(ir,ip,fi))
            {
            //  cout <<  rep_rate1(pp1) << " " <<  rep_rate(ir,ip,fi) << endl;
            }
            //  dvar_vector& rtc=rep_rate(pp1)*tagcatch(it,ir,ip,fi);
              dvar_vector rtc= sqrt(0.375+obstagcatch(it,ir,ip,fi))-
                sqrt(0.375+rep_rate(ir,ip,fi)*tagcatch(it,ir,ip,fi));
              f-=sum(log(exp(square(rtc)/(-0.4))+.01));
            }
          else
          {
            dvar_matrix obsgroupedcatch(1,num_fisheries,1,nage);
            dvar_matrix groupedcatch(1,num_fisheries,1,nage);
            ivector iflag(1,num_fisheries);
            obsgroupedcatch.initialize();
            groupedcatch.initialize();
            iflag.initialize();
            for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
            {
              int pp1=parent(ir,ip,fi);
              //dvar_vector& rtc=rep_rate(pp1)*tagcatch(it,ir,ip,fi);
              if (rep_rate1(pp1) != rep_rate(ir,ip,fi))
              {
                //cout <<  rep_rate1(pp1) << " " <<  rep_rate(ir,ip,fi) << endl;
              }
              int pp=group_flags(pp1);
              const dvar_vector& rtc=rep_rate(ir,ip,fi)*tagcatch(it,ir,ip,fi);
              iflag(pp)=1;
              obsgroupedcatch(pp)+=obstagcatch(it,ir,ip,fi); 
              groupedcatch(pp)+=rtc;
            }
            for (fi=1;fi<=num_fisheries;fi++)
            {
              if(iflag(fi))
              {
                dvar_vector rtc= sqrt(0.375+obsgroupedcatch(fi))-
                  sqrt(0.375+groupedcatch(fi));
                f-=sum(log(exp(square(rtc)/(-0.4))+.01));
              }
            } 
          }
        }
      }
    }
  }
  return f;
}

void dvar_len_fish_stock_history::print_tag_data(int ir,int ip,
  ofstream& of)
{
  ivector group_flags=column(fish_flags,44);
  int gsum44=sum(group_flags);
  dvariable f=0.0;
  int fi;
  dvar_vector& rep_rate1=fish_pars(3);
  for (int it=1;it<=num_tag_releases;it++)
  {
    if (ip<initial_tag_period(it,ir)) break;
    if (!(!(terminal_tag_period)))
    {
      if (ip>terminal_tag_period(it,ir)) break;
    }
    int ub;
    if (!age_flags(96))
      ub=num_fish_periods(ir);
    else
      ub=terminal_tag_period(it,ir);
    {
      for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {
        if (!gsum44)
        {
          // don't group fisheries
          for (int j=1;j<=nage;j++)
          {
            if (tagcatch(it,ir,ip,fi,j)<0)
            {
              cerr << tagcatch(it,ir,ip,fi,j) << endl;
            }
          }
          int pp1=parent(ir,ip,fi);
          const dvector& rtc=value(rep_rate(ir,ip,fi))*value(tagcatch(it,ir,ip,fi));
          of << endl << "     " << sum(rtc) << " " << sum(obstagcatch(it,ir,ip,fi)) << " "
             << sum(elem_div(square(value(obstagcatch(it,ir,ip,fi))-rtc),
               .01+rtc))
             << endl << setfixed() << setprecision(2)
	     << rtc/(1.e-10+sum(rtc))
             << setfixed() << setprecision(2)  << endl
             << obstagcatch(it,ir,ip,fi)/(1.e-20+sum(obstagcatch(it,ir,ip,fi)));
        }
        else
        {
          dvar_matrix obsgroupedcatch(1,num_fisheries,1,nage);
          dvar_matrix groupedcatch(1,num_fisheries,1,nage);
          ivector iflag(1,num_fisheries);
          obsgroupedcatch.initialize();
          groupedcatch.initialize();
          iflag.initialize();
          for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
          {
            int pp1=parent(ir,ip,fi);
            int pp=group_flags(pp1);
            const dvar_vector& rtc=rep_rate(ir,ip,fi)*tagcatch(it,ir,ip,fi);
            iflag(pp)=1;
            obsgroupedcatch(pp)+=obstagcatch(it,ir,ip,fi); 
            groupedcatch(pp)+=rtc;
          }
          of << endl;
          for (fi=1;fi<=num_fisheries;fi++)
          {
            if(iflag(fi))
            {
              int yr=year(ir,ip);
              int mnth=month(ir,ip);
              if (month_factor!=0 && month_factor!=1)
              {         
                date_struc newdate=get_old_time(year(ir,ip),
                  month(ir,ip),1,month_factor,first_time);
                yr=newdate.year;
                mnth=newdate.month;
              }
              of << setw(4) << "    " << it << " " 
                 << setw(4) << yr << " " 
                 << setw(4) << mnth << "   "
                 << setw(10) << sum(groupedcatch(fi)) << " " 
                 << setw(10) << sum(obsgroupedcatch(fi)) << "   "
                  << setw(10) << sum(elem_div(
                    square(value(obsgroupedcatch(fi))-value(groupedcatch(fi))),
                 .01+value(groupedcatch(fi)))) << "  ";
               of  <<  value(obsgroupedcatch(fi))*
                     log(elem_div(value(1.e-10+groupedcatch(fi)),
                         value(1.e-10+obsgroupedcatch(fi))))
                   << endl;
            }
          } 
        }
      }
    }
  }
}

// penalties to keep reporting rates equal for different fisheries
dvariable dvar_len_fish_stock_history::grouped_tag_reporting_rate_penalty(void)
{
  ivector group_flags=column(fish_flags,34);
  dvariable gp_pen=0.0;
  int fi;
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=min_init_tag_period(ir);ip<=num_fish_periods(ir);ip++)
    {
      for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {
        ivector& pi=parent(ir,ip);
        for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
        {
          for (fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
          {
            for (int fi1=1;fi1<=fi;fi1++)
            {
              if (group_flags(pi(fi))== group_flags(pi(fi1)))
	      {
                gp_pen+=100.*square(log(
                  (1.e-5+rep_rate(ir,ip,fi))/
                  (1.e-5+rep_rate(ir,ip,fi1)) ));
              }
            }
          }
        }
      }
    }
  }
  return gp_pen;
}


