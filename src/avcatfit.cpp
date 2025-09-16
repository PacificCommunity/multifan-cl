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

#include "newmprot.hpp"

 ofstream ofsuu("totalcatchpen");

int dvar_fish_stock_history::is_projection_period(int ir,int ip)
{
  if (year(ir,ip) < projection_year(1) )
    return 0;
  if ( year(ir,ip) == projection_year(1) && 
        month(ir,ip) < projection_month(1) )
  {
    return 0;
  }
  else
  {
    return 1;
  }

  if ( (year(ir,ip) < projection_year(1) )
    || ( year(ir,ip) == projection_year(1) && 
        month(ir,ip) < projection_month(1) ) )
  {
    return 0;
  }
  else
  {
     return 1;
  }
}


dvariable total_catch_or_weight_fit(dvar_len_fish_stock_history& fsh,
  int print_switch,int avg_calc_flag)
{
  dvariable xy=0.;
  dvariable tmp1=0.;
  dvariable tmp=0.;
  int ntimes;
  int catflg=sum(column(fsh.fish_flags,45));
  MY_DOUBLE_TYPE total_catch_error_sq=0.0;
  MY_DOUBLE_TYPE quadpen=0.0;
  for (int ir=1;ir<=fsh.num_regions;ir++)
  {
    if (fsh.parest_flags(142)==0)
    {
      ntimes=fsh.num_real_fish_periods(ir);
    }
    else
    {
      ntimes=min(fsh.parest_flags(142),fsh.num_real_fish_periods(ir));
    }
    for (int ip=1;ip<=ntimes;ip++)
    {
      for (int fi=1;fi<=fsh.num_fish_incidents(ir,ip);fi++)
      {
        if (fsh.obs_tot_catch(ir,ip,fi)>-0.5L)
        {
          if (!fsh.data_fish_flags(1,fsh.parent(ir,ip,fi))) 
          {
            if (catflg)
            {
              if (avg_calc_flag)
              {
                fsh.totalcatch_by_numbers(ir)+=fsh.tot_catch(ir,ip,fi);
                fsh.obstotalcatch_by_numbers(ir)+=fsh.obs_tot_catch(ir,ip,fi);
                fsh.numtotalcatch_by_numbers(ir)+=1;
              }
              else
              {
                tmp=.01*fsh.fish_flags(fsh.parent(ir,ip,fi),45)
                      *square(log(1.+fsh.tot_catch(ir,ip,fi))
                   -log(1.+fsh.obs_tot_catch(ir,ip,fi)));
                tmp1+=tmp;
//                ofsuu << tmp << " ";
                if (fsh.ppstf && !sum(fsh.q_flag) && print_switch)  //NMD_6sep2021
                {
                  fsh.ppstf->tot_catch_like_by_realization(ir,ip,fi)
                    =value(tmp);
                }
              }
            }
            else
            {
              if (avg_calc_flag)
              {
                fsh.totalcatch_by_numbers(ir)+=fsh.tot_catch(ir,ip,fi);
                fsh.obstotalcatch_by_numbers(ir)+=fsh.obs_tot_catch(ir,ip,fi);
                fsh.numtotalcatch_by_numbers(ir)+=1;
              }
              else
              {
                if (!fsh.pmsd || fsh.pmsd->fisn(ir,ip,fi)==1)
                {
                  if (fsh.parest_flags(399))
                  {
                    fsh.lagrange_c(ir,ip,fi)=
                      log(1.+fsh.tot_catch(ir,ip,fi))
                      -log(1.+fsh.obs_tot_catch(ir,ip,fi));
                    total_catch_error_sq
                      +=square(value(fsh.lagrange_c(ir,ip,fi)));
                    fsh.lagrange_mu=.02*fsh.age_flags(144);
                    tmp=0.5*fsh.lagrange_mu*square(fsh.lagrange_c(ir,ip,fi));
                    quadpen+=value(tmp);
                    tmp-=fsh.lagrange_lambda(ir,ip,fi)*fsh.lagrange_c(ir,ip,fi);
                  }
                  else
                  {
                    tmp=.01*fsh.age_flags(144)*
                      square(log(1.+fsh.tot_catch(ir,ip,fi))
                       -log(1.+fsh.obs_tot_catch(ir,ip,fi)));
//                    if (fsh.ppstf)
                    if (fsh.ppstf && !sum(fsh.q_flag) && print_switch)  //NMD_13jul2017
                    {
                      fsh.ppstf->tot_catch_like_by_realization(ir,ip,fi)
                        =value(tmp);
                    }
                  }
                  tmp1+=tmp;
//                  ofsuu << tmp << " ";
                }
                else
                { 
                  if (ir==fsh.pmsd->reg_in_catch(ir,ip,fi,1))
                  {
                    dvariable tmptot=fsh.tot_catch(ir,ip,fi);
                    int n=fsh.pmsd->fisn(ir,ip,fi);
                    for (int i=2;i<=n;i++)
                    {
                      int rr=fsh.pmsd->reg_in_catch(ir,ip,fi,i);
                      tmptot+=fsh.tot_catch(rr,ip,fi);
                    }
                    tmp=.01*fsh.age_flags(144)*
                      square(log(1.+tmptot)
                      -log(1.+fsh.obs_tot_catch(ir,ip,fi)));
                    if (fsh.ppstf && !sum(fsh.q_flag) && print_switch)  //NMD_2oct2019
                    {
                      fsh.ppstf->tot_catch_like_by_realization(ir,ip,fi)
                        =value(tmp);
                    }
                    tmp1+=tmp;  //NMD_2oct2019
                  }
                }
              }
            }
          }
          else
          {
            if (!fsh.pmsd || fsh.pmsd->fisn(ir,ip,fi)==1)
            {
              fsh.pred_totcatch(ir,ip,fi)=0.0; 
              dvariable sv27=fsh.get_sv_region(ir,27);  //NMD 9Dec2011
              dvariable sv28=fsh.get_sv_region(ir,28);  //NMD 9Dec2011
              if (value(sv28)==3.0)
              {
                fsh.pmsd_error(); 
                fsh.pred_totcatch(ir,ip,fi)=fsh.len_wt_coff/1000.* 
                   (fsh.exp_catch(ir,ip,fi)*
                   (pow(fsh.mean_length(ir,ip,fi),3)+
                   3.0*elem_prod(fsh.mean_length(ir,ip,fi),fsh.vars(ir,ip,fi))));
              }
              else
              {
                for (int j=1;j<=fsh.nage_by_region(ir);j++)
                {  
                  fsh.pred_totcatch(ir,ip,fi)+=fsh.exp_catch(ir,ip,fi,j)*
                    normal_length_to_weight(0.5,-3.5,
                    fsh.mean_length(ir,ip,fi,j),
                      sqrt(fsh.vars(ir,ip,fi,j)),value(sv27),value(sv28)); 
                }
  
                fsh.pred_totcatch(ir,ip,fi)/=1000.; 
              }
              if (fsh.pred_totcatch(ir,ip,fi)<0.0)
              {
                cout <<"avcatfit.cpp " << fsh.len_wt_coff << endl;
                cout <<"avcatfit.cpp " << fsh.catch(ir,ip,fi) << endl;
                cout <<"avcatfit.cpp " << fsh.mean_length(ir,ip,fi) << endl;
              }
              if (catflg)
              {
                if (avg_calc_flag)
                {
                  fsh.totalcatch_by_weight(ir)+=fsh.pred_totcatch(ir,ip,fi);
                  fsh.obstotalcatch_by_weight(ir)+=fsh.obs_tot_catch(ir,ip,fi);
                  fsh.numtotalcatch_by_weight(ir)+=1;
                }
                else
                {
                  tmp=.01*fsh.fish_flags(fsh.parent(ir,ip,fi),45)
                        *square(log(1.+fsh.pred_totcatch(ir,ip,fi))
                     -log(1.+fsh.obs_tot_catch(ir,ip,fi)));
                  tmp1+=tmp;
//                  ofsuu << tmp << " ";
                  if (fsh.ppstf && !sum(fsh.q_flag) && print_switch)  //NMD_6sep2021
                  {
                    fsh.ppstf->tot_catch_like_by_realization(ir,ip,fi)
                      =value(tmp);
                  }
                }
              }
              else
              {
                if (avg_calc_flag)
                {
                  fsh.totalcatch_by_weight(ir)+=fsh.pred_totcatch(ir,ip,fi);
                  fsh.obstotalcatch_by_weight(ir)+=fsh.obs_tot_catch(ir,ip,fi);
                  fsh.numtotalcatch_by_weight(ir)+=1;
                }
                else
                {
                  tmp=.01*fsh.age_flags(144)*square(log(1.+fsh.pred_totcatch(ir,ip,fi))
                     -log(1.+fsh.obs_tot_catch(ir,ip,fi)));
                  tmp1+=tmp;
//                  ofsuu << tmp << " ";
                  if (fsh.ppstf && !sum(fsh.q_flag) && print_switch)  //NMD_13jul2017
                  {
                    fsh.ppstf->tot_catch_like_by_realization(ir,ip,fi)
                      =value(tmp);
                  }
                }
              }
            }
            else
            {
              if (ir==fsh.pmsd->reg_in_catch(ir,ip,fi,1))
              {
                dvariable tmptot=0.0; 
                int n=fsh.pmsd->fisn(ir,ip,fi);
                for (int i=1;i<=n;i++)
                {
                  int rr=fsh.pmsd->reg_in_catch(ir,ip,fi,i);

                  dvariable sv27=fsh.get_sv_region(rr,27);  //NMD 9Dec2011
                  dvariable sv28=fsh.get_sv_region(rr,28);  //NMD 9Dec2011
                  if (value(sv28)==3.0)
                  {
                    fsh.pmsd_error(); 
                    tmptot+=fsh.len_wt_coff/1000.* 
                      (fsh.exp_catch(rr,ip,fi)*
                      (pow(fsh.mean_length(rr,ip,fi),3)+
                      3.0*elem_prod(fsh.mean_length(rr,ip,fi),fsh.vars(rr,ip,fi))));
                  }
                  else
                  {
                    for (int j=1;j<=fsh.nage;j++)
                    {  
                      tmptot+=fsh.exp_catch(rr,ip,fi,j)*
                        normal_length_to_weight(0.5,-3.5,
                        fsh.mean_length(rr,ip,fi,j),
                        sqrt(fsh.vars(rr,ip,fi,j)),value(sv27),value(sv28)); 
                    }
                  }
                }
                fsh.pred_totcatch(ir,ip,fi)=tmptot/1000.; 
                if (fsh.pred_totcatch(ir,ip,fi)<0.0)
                {
                  cout <<"avcatfit.cpp " << fsh.len_wt_coff << endl;
                  cout <<"avcatfit.cpp " << fsh.catch(ir,ip,fi) << endl;
                  cout <<"avcatfit.cpp " << fsh.mean_length(ir,ip,fi) << endl;
                }
                if (catflg)
                {
                  if (avg_calc_flag)
                  {
                    fsh.totalcatch_by_weight(ir)+=fsh.pred_totcatch(ir,ip,fi);
                    fsh.obstotalcatch_by_weight(ir)+=fsh.obs_tot_catch(ir,ip,fi);
                    fsh.numtotalcatch_by_weight(ir)+=1;
                  }
                  else
                  {
                    tmp=.01*fsh.fish_flags(fsh.parent(ir,ip,fi),45)
                          *square(log(1.+fsh.pred_totcatch(ir,ip,fi))
                       -log(1.+fsh.obs_tot_catch(ir,ip,fi)));
                    tmp1+=tmp;
//                    ofsuu << tmp << " ";
                    if (fsh.ppstf && !sum(fsh.q_flag) && print_switch)  //NMD_2oct2019
                    {
                      fsh.ppstf->tot_catch_like_by_realization(ir,ip,fi)
                        =value(tmp);
                    }
                  }
                }
                else
                {
                  if (avg_calc_flag)
                  {
                    fsh.totalcatch_by_weight(ir)+=fsh.pred_totcatch(ir,ip,fi);
                    fsh.obstotalcatch_by_weight(ir)+=fsh.obs_tot_catch(ir,ip,fi);
                    fsh.numtotalcatch_by_weight(ir)+=1;
                  }
                  else
                  {
                    tmp=.01*fsh.age_flags(144)*square(log(1.+fsh.pred_totcatch(ir,ip,fi))
                       -log(1.+fsh.obs_tot_catch(ir,ip,fi)));
                    tmp1+=tmp;
//                    ofsuu << tmp << " ";
                    if (fsh.ppstf && !sum(fsh.q_flag) && print_switch)  //NMD_2oct2019
                    {
                      fsh.ppstf->tot_catch_like_by_realization(ir,ip,fi)
                        =value(tmp);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  xy+=tmp1;
  /*
  if (print_switch)
  {
    cout << " after total catch_or_weight fit = " << xy << endl;
    cout << "exp(q0)= " << setprecision(10) <<exp(fsh.q0) << "  "
       << setprecision(10) << exp(fsh.q0(1)) << endl;
    cout << "total catch contribution = "<< tmp1 << endl;
  }
  */
//  ofsuu << "Total catch main " << setfixed() << setprecision(9) << setw(12) << tmp1 << endl << endl;
  cout << "Total catch main " << setfixed() << setprecision(9) << setw(12) << tmp1 << 
   "   Total_catch_quadpen " << setfixed() << setprecision(9) << setw(12) <<  quadpen << 
   "   Total_catch_error_sq " << setfixed() << setprecision(9) << setw(12) <<  total_catch_error_sq << endl;//NMD
  return xy;
}
#undef HOME_VERSION
