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

dvariable total_catch_or_weight_fit_ms(dvar_len_fish_stock_history& fsh,
  int print_switch,int avg_calc_flag)
{
  dvariable xy=0.;
  dvariable tmp1=0.;
  dvariable tmp=0.;
  int ntimes;
  int catflg=sum(column(fsh.fish_flags,45));
  int nr=fsh.pmsd->num_real_regions;
  for (int ir=1;ir<=nr;ir++)
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
              fsh.pmsd_error();
              if (avg_calc_flag)
              {
                fsh.totalcatch_by_numbers(ir)+=fsh.tot_catch(ir,ip,fi);
                fsh.obstotalcatch_by_numbers(ir)+=fsh.obs_tot_catch(ir,ip,fi);
                fsh.numtotalcatch_by_numbers(ir)+=1;
              }
              else
              {
                dvariable tot=0.0;
                for (int is=1;is<=fsh.pmsd->num_species;is++)
                {
                  int irr=ir+fsh.pmsd->rowoffset(is);
                  if (fsh.pmsd->fisc(irr,ip,fi,is)>0)
                  {
                    tot+=fsh.tot_catch(irr,ip,fi);
                  }
                }
                tmp=.01*fsh.fish_flags(fsh.parent(ir,ip,fi),45)
                   *square(log(1.+tot)
                   -log(1.+fsh.obs_tot_catch(ir,ip,fi)));
                tmp1+=tmp;
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
                dvariable tot=0.0;
                for (int is=1;is<=fsh.pmsd->num_species;is++)
                {
                  int irr=ir+fsh.pmsd->rowoffset(is);
                  if (fsh.pmsd->fisc(irr,ip,fi,is)>0)
                  {
                    tot+=fsh.tot_catch(irr,ip,fi);
                  }
                }
                tmp=.01*fsh.age_flags(144)*square(log(1.+tot)
                   -log(1.+fsh.obs_tot_catch(ir,ip,fi)));

                tmp1+=tmp;
              }
            }
          }
          else
          {
            dvariable tot=0.0;
            for (int is=1;is<=fsh.pmsd->num_species;is++)
            {
              int irr=ir+fsh.pmsd->rowoffset(is);
              fsh.pred_totcatch(irr,ip,fi)=0.0; 
              dvariable sv27=fsh.get_sv_region(irr,27);  //NMD 9Dec2011
              dvariable sv28=fsh.get_sv_region(irr,28);  //NMD 9Dec2011
              dvariable lwc=fsh.get_len_wt_coff_region(irr);
              if (value(sv28)==3.0)
              {
                fsh.pred_totcatch(irr,ip,fi)=fsh.len_wt_coff/1000.* 
                   (fsh.exp_catch(irr,ip,fi)*
                   (pow(fsh.mean_length(irr,ip,fi),3)+
                   3.0*elem_prod(fsh.mean_length(irr,ip,fi),fsh.vars(irr,ip,fi))));
              }
              else
              {
                for (int j=1;j<=fsh.nage;j++)
                {  
                  fsh.pred_totcatch(irr,ip,fi)+=fsh.exp_catch(irr,ip,fi,j)*
                    normal_length_to_weight(0.5,-3.5,
                    fsh.mean_length(irr,ip,fi,j),
                      sqrt(fsh.vars(irr,ip,fi,j)),value(sv27),value(sv28)); 
                }
  
                fsh.pred_totcatch(irr,ip,fi)/=1000.; 
              }
              if (fsh.pred_totcatch(irr,ip,fi)<0.0)
              {
                cout <<"avcatfit.cpp " << fsh.len_wt_coff << endl;
                cout <<"avcatfit.cpp " << fsh.catch(irr,ip,fi) << endl;
                cout <<"avcatfit.cpp " << fsh.mean_length(irr,ip,fi) << endl;
              }
              if (catflg)
              {
                fsh.pmsd_error();
                if (avg_calc_flag)
                {
                  fsh.totalcatch_by_weight(irr)+=fsh.pred_totcatch(irr,ip,fi);
                  fsh.obstotalcatch_by_weight(irr)+=fsh.obs_tot_catch(irr,ip,fi);
                  fsh.numtotalcatch_by_weight(irr)+=1;
                }
                else
                {
                  if (fsh.pmsd->fisc(irr,ip,fi,is)>0)
                  {
                    tot+=fsh.pred_totcatch(irr,ip,fi);
                  }
                  if (is==fsh.pmsd->num_species)
                  {
                    tmp=.01*fsh.fish_flags(fsh.parent(irr,ip,fi),45)
                        *square(log(1.+tot)
                     -log(1.+fsh.obs_tot_catch(ir,ip,fi)));
                    tmp1+=tmp;
                  }
                }
              }
              else
              {
                if (avg_calc_flag)
                {
                  fsh.totalcatch_by_weight(irr)+=fsh.pred_totcatch(irr,ip,fi);
                  fsh.obstotalcatch_by_weight(irr)+=fsh.obs_tot_catch(irr,ip,fi);
                  fsh.numtotalcatch_by_weight(irr)+=1;
                }
                else
                {
                  if (fsh.pmsd->fisc(irr,ip,fi,is)>0)
                  {
                    tot+=fsh.pred_totcatch(irr,ip,fi);
                  }
                  if (is==fsh.pmsd->num_species)
                  {
                    tmp=.01*fsh.age_flags(144)*square(log(1.+tot)
                     -log(1.+fsh.obs_tot_catch(ir,ip,fi)));
                    tmp1+=tmp;
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
  cout << "tmp1=" << tmp1 << endl;
  if (print_switch)
  {
    cout << " after total catch_or_weight fit = " << xy << endl;
    cout << "exp(q0)= " << setprecision(10) <<exp(fsh.q0) << "  "
       << setprecision(10) << exp(fsh.q0(1)) << endl;
    cout << "total catch contribution = "<< tmp1 << endl;
  }
  cout << "Total catch " << setfixed() << setprecision(2) << setw(15) << tmp1 << endl;
  return xy;
}
