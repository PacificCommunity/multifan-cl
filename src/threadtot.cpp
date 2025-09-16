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

#include "threaded_tag3.h"
#include "adthread.h"

/*
class xxmspd
{
public:
    //i4_array fisc; 
    i3_array fisn; 
    //i4_array sp_in_catch;
    i4_array reg_in_catch;
    //dvar_matrix rec_delta;
    //dvar_vector recmean;
    //dvar_vector totpop_coff;
};


class  xxtotcatwt
{
public:
  // constant
  imatrix      fish_flags; 
  ivector      num_real_fish_periods;
  int          num_regions; 
  ivector      parest_flags; 
  d3_array     obs_tot_catch;
  imatrix      num_fish_incidents; 
  ivector      age_flags; 
  dvector obstotalcatch_by_weight;
  ivector numtotalcatch_by_weight;
  int nage;
  dvector obstotalcatch_by_numbers;
  ivector numtotalcatch_by_numbers;
  i3_array     parent;  
  imatrix      data_fish_flags;     
  MY_DOUBLE_TYPE len_wt_coff;

  // variables
  dvar3_array      tot_catch;
  dvar_vector totalcatch_by_weight;
  dvar3_array     pred_totcatch; 
  dvar4_array     catch;   
  dvar4_array mean_length;
  dvar4_array vars;
  dvar_vector totalcatch_by_numbers;
  dvar4_array     exp_catch;    
  dvar_vector      q0;  
  xxmspd * pmsd;
  void pmsd_error();
  dvariable get_sv_region(int is,int svind);
  dvar_vector sv;
};
*/


dvariable total_catch_or_weight_fit(xxtotcatwt& fsh,
  int print_switch,int avg_calc_flag)
{
  dvariable xy=0.;
  dvariable tmp1=0.;
  dvariable tmp=0.;
  int ntimes;
  int catflg=sum(column(fsh.fish_flags,45));
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
                //if (!fsh.pmsd || fsh.pmsd->fisn(ir,ip,fi)==1)
                {
                  tmp=.01*fsh.age_flags(144)*
                    square(log(1.+fsh.tot_catch(ir,ip,fi))
                    -log(1.+fsh.obs_tot_catch(ir,ip,fi)));
                  tmp1+=tmp;
                }
               /*
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
                  }
                }
               */
              }
            }
          }
          else
          {
            //if (!fsh.pmsd || fsh.pmsd->fisn(ir,ip,fi)==1)
            {
              fsh.pred_totcatch(ir,ip,fi)=0.0; 
              //dvariable sv27=fsh.get_sv_region(ir,27);  //NMD 9Dec2011
              //dvariable sv28=fsh.get_sv_region(ir,28);  //NMD 9Dec2011
              dvariable sv27=fsh.sv(27);  
              dvariable sv28=fsh.sv(28); 
              if (value(sv28)==3.0)
              {
                //fsh.pmsd_error(); 
                fsh.pred_totcatch(ir,ip,fi)=fsh.len_wt_coff/1000.* 
                   (fsh.exp_catch(ir,ip,fi)*
                   (pow(fsh.mean_length(ir,ip,fi),3)+
                   3.0*elem_prod(fsh.mean_length(ir,ip,fi),fsh.vars(ir,ip,fi))));
              }
              else
              {
                for (int j=1;j<=fsh.nage;j++)
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
                }
              }
            }
           /*
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
                    //fsh.pmsd_error(); 
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
                  }
                }
              }
            }
           */
          }
        }
      }
    }
  }
  xy+=tmp1;
  if (print_switch)
  {
    cout << " after total catch_or_weight fit = " << xy << endl;
    cout << "exp(q0)= " << setprecision(10) <<exp(fsh.q0) << "  "
       << setprecision(10) << exp(fsh.q0(1)) << endl;
    cout << "total catch contribution = "<< tmp1 << endl;
  }
  cout << "Total catch in threadtot" << setfixed() << setprecision(9) << setw(19) << tmp1 << "  "  << xy   << endl;//NMD
  return xy;
}

void admb_catch_wt(void *ptr)
{
  new_thread_data * tptr = (new_thread_data *) ptr;
  int sno=tptr->thread_no;
  ad_comm::pthread_manager->set_slave_number(sno);

  //get information from master
  pthread_mutex_lock(&ad_comm::pthread_manager->copy_mutex);
  xxtotcatwt cwt;
  pthread_mutex_unlock(&ad_comm::pthread_manager->copy_mutex);
  cwt.catch_wt_loop_3();
}


void xxtotcatwt::catch_wt_loop_3(void)
{
  int tn=0;
  get_constant_data_3(tn);
  gradient_structure::set_MAX_NVAR_OFFSET(10000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(20000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(120000000);
  gradient_structure::set_MAX_NVAR_OFFSET(30000);
  gradient_structure gs(100000000);
  adtimer adt1;

  allocate();
  do
  {
    int mmin=ad_comm::pthread_manager->gmin(1);
    int mmax=ad_comm::pthread_manager->gmax(1);
    int iflag;

    int tn=0;
    iflag=get_variables_from_master3(tn);

    if (iflag==1) break;
    adt1.get_elapsed_time_and_reset();
    dvar_vector sv;
    dvariable ffpen=0.0;
    adt1.get_elapsed_time_and_reset();
    dvariable tmp=total_catch_or_weight_fit(*this,0,0);
    MY_DOUBLE_TYPE tm=adt1.get_elapsed_time_and_reset();
    cout << "total_catch_or_weight_fit in loop 3 "  << tmp << endl;
    int nt=0;
    // send results of fit_pooled_tags to master thread
    ad_comm::pthread_manager->write_lock_buffer(nt);
    ad_comm::pthread_manager->send_dvariable(tmp,nt); 
    ad_comm::pthread_manager->write_unlock_buffer(nt);
    slave_gradcalc();
  }
  while(1);
}
void xxtotcatwt::allocate(void)
{
  pred_totcatch.allocate(1,num_regions,1,num_fish_periods,1,num_fish_incidents);
}

void xxtotcatwt::get_constant_data_3(int tn)
{
  ad_comm::pthread_manager->cread_lock_buffer(tn);
  num_regions=ad_comm::pthread_manager->get_int(tn); 
  num_fish_incidents=ad_comm::pthread_manager->get_imatrix(tn); 
  nage=ad_comm::pthread_manager->get_int(tn); 
  age_flags=ad_comm::pthread_manager->get_ivector(tn); 
  fish_flags=ad_comm::pthread_manager->get_imatrix(tn);  
  num_real_fish_periods=ad_comm::pthread_manager->get_ivector(tn); 
  parest_flags=ad_comm::pthread_manager->get_ivector(tn);
  obs_tot_catch=ad_comm::pthread_manager->get_d3_array(tn);
  parent=ad_comm::pthread_manager->get_i3_array(tn);
  data_fish_flags=ad_comm::pthread_manager->get_imatrix(tn);
  len_wt_coff=ad_comm::pthread_manager->get_double(tn);
  avg_calc_flag=ad_comm::pthread_manager->get_int(tn);
  num_fish_periods=ad_comm::pthread_manager->get_ivector(tn);
  ad_comm::pthread_manager->cread_unlock_buffer(tn);


}


void dvar_len_fish_stock_history::send_constant_data_3(int sno)
{
  ad_comm::pthread_manager->cwrite_lock_buffer(sno);
  ad_comm::pthread_manager->send_int(num_regions,sno); 
  ad_comm::pthread_manager->send_imatrix(num_fish_incidents,sno); 
  ad_comm::pthread_manager->send_int(nage,sno); 
  ad_comm::pthread_manager->send_ivector(age_flags,sno); 
  ad_comm::pthread_manager->send_imatrix(fish_flags,sno);  
  ad_comm::pthread_manager->send_ivector(num_real_fish_periods,sno); 
  ad_comm::pthread_manager->send_ivector(parest_flags,sno);
  ad_comm::pthread_manager->send_d3_array(obs_tot_catch,sno);
  ad_comm::pthread_manager->send_i3_array(parent,sno);
  ad_comm::pthread_manager->send_imatrix(data_fish_flags,sno);
  ad_comm::pthread_manager->send_double(len_wt_coff,sno);
  ad_comm::pthread_manager->send_int(0,sno);
  ad_comm::pthread_manager->send_ivector(num_fish_periods,sno);
  ad_comm::pthread_manager->cwrite_unlock_buffer(sno);
}

int dvar_len_fish_stock_history::send_variables_to_tag_slave3(int flag,
  int slave_number)
{
  ad_comm::pthread_manager->write_lock_buffer(slave_number);
  ad_comm::pthread_manager->send_int(flag,slave_number);
  if (pmsd)
  {
  //  ad_comm::pthread_manager->send_dvar_matrix(pmsd->vb_coff,slave_number); 
  }
  ad_comm::pthread_manager->send_dvar4_array(mean_length,slave_number);
  ad_comm::pthread_manager->send_dvar4_array(vars,slave_number);
  ad_comm::pthread_manager->send_dvar4_array(exp_catch,slave_number);
  ad_comm::pthread_manager->send_dvar_vector(sv,slave_number);
  ad_comm::pthread_manager->send_dvar3_array(tot_catch,slave_number);
  long int offset=ad_comm::pthread_manager->get_offset(slave_number);
  cout << "sent " << offset << " bytes in send slave 3" << endl;
  ad_comm::pthread_manager->write_unlock_buffer(slave_number);
  return 0;
}

int xxtotcatwt::get_variables_from_master3(int tn)
{
  ad_comm::pthread_manager->read_lock_buffer(tn);
  int flag=ad_comm::pthread_manager->get_int(tn);
  mean_length=ad_comm::pthread_manager->get_dvar4_array(tn);
  vars=ad_comm::pthread_manager->get_dvar4_array(tn);
  exp_catch=ad_comm::pthread_manager->get_dvar4_array(tn);
  sv=ad_comm::pthread_manager->get_dvar_vector(tn);
  tot_catch=ad_comm::pthread_manager->get_dvar3_array(tn);
  ad_comm::pthread_manager->read_unlock_buffer(tn);
  return 0;
}

#undef HOME_VERSION
