/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"

  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
 long int diff(long int&,const long int&);

 void prob_calc(prevariable& prob,prevariable& sigg,prevariable& mean_len,
   MY_DOUBLE_TYPE& t);
 void prob_calc2(prevariable& prob,prevariable& sigg,prevariable& mean_len,
   MY_DOUBLE_TYPE& t);
 void prob_calc3(prevariable& prob,prevariable& sigg,prevariable& mean_len,
   MY_DOUBLE_TYPE& t);
 void DF_main_weight_calcs(void);


void dvar_len_fish_stock_history::main_weight_calcs_len_based
  (dvar_vector& tprobb,dvar_vector& propp,dvar_matrix& alis)
{
  tprobb=propp*alis;
  //tprobb/=sum(1.e-20+tprobb);
}


void dvar_len_fish_stock_history::fast_weight_pred_frequency_calc(void)
{
  wtprob.initialize();
  //tprob+=1.e-5;
  int i;
  int j;
  int js;
  dvariable temp;
  dvariable fdiff2;
  dvariable prob;
  dvariable u;
  //double cutoff = 5.e0;
  //double ecut = .0067379469991e0;
  MY_DOUBLE_TYPE a2 = 35.91908881;
  //double cutoff1,u;
  MY_DOUBLE_TYPE a3 = 785.858644;
  //double a = 0.9;
  long int ztmp=gradient_structure::totalbytes();
  long int ztmp1;

  //NMD_24jan2023
  /*
  dmatrix species_wm2;
  dvector wm2=square(wmid);
  if (pmsd)
  {
    species_wm2=square(pmsd->wmid);
  }
  */
  for (int ir=1;ir<=num_regions;ir++)
  {
    if(pmsd) pmsd->current_species=pmsd->region_species_pointer(ir);  //NMD 21Aug2012
    int np=num_real_fish_periods(ir);
    if (projected_simulated_data_flags[1]==1)
    {
      np=num_fish_periods(ir);
    }

    for (int ip=1;ip<=np;ip++)  // Loop over fishing periods
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        if (wght_sample_size(ir,ip,fi)>0)
        {
          js=1;
          dvar_vector& tprobb=wtprob(ir,ip,fi);
          dvar_vector& propp=prop(ir,ip,fi);
          dvar_vector& mean_len=mean_length(ir,ip,fi);
          //dvar_vector& mean_wght=mean_weight(ir,ip,fi);
          dvar_vector& sigg=sdevs(ir,ip,fi);
          int pi=parent(ir,ip,fi);
          if (sum(value(propp))<1.e-80)
          {
            cerr << "This can't happen" << endl;
            ad_exit(1);
          }
          if (fish_flags(pi,57)==3 && fish_flags(pi,26)==3)
          {
            //age_weight_sel(1,12,1,4,1,num_fisheries,1,nage,1,nwint);
            int mn=month(ir,ip);
            int wk=week(ir,ip);
            dvar_matrix& aws=age_weight_fishery_size_dist(mn,wk,pi);
            main_weight_calcs_len_based(tprobb,propp,aws);
            if (pmsd==0 || pmsd->current_species==1)
              tprobb=elem_div(tprobb,wm2);
            else
              tprobb=elem_div(tprobb,pmsd->wm2(pmsd->current_species));
            tprobb/=(1.e-20+sum(tprobb));
          }
          else
          {
            main_weight_calcs(tprobb,propp,mean_len,sigg);
            if (pmsd==0 || pmsd->current_species==1)
              tprobb=elem_div(tprobb,wm2);
            else
              tprobb=elem_div(tprobb,pmsd->wm2(pmsd->current_species));
            tprobb/=sum(tprobb);
          }
        }
      }
    }
  }
  ztmp1=diff(ztmp,gradient_structure::totalbytes());
}

/*
void dvar_len_fish_stock_history::main_weight_calcs(dvar_vector& tprobb,
  dvar_vector& propp,dvar_vector& mean_len,dvar_vector& sigg)
{
  const MY_DOUBLE_TYPE a2 = 35.91908881;
  const MY_DOUBLE_TYPE a3 = 785.858644;
  const MY_DOUBLE_TYPE a22 = 2*35.91908881;
  const MY_DOUBLE_TYPE a33 = 3*785.858644;
  dvariable prob=0.0;
  int js=1;
  int break_flag=0;
  for (int i=1; i<=nwint; i++) 
  {
    MY_DOUBLE_TYPE& fmidd=wmid(i);
   
    for (int j=1; j<=nage; j++) 
    {
      dvariable t=(fmidd-mean_len(j))/sigg(j);
      if ( fabs(t) <= 3.03)
      {
        if ( fabs(t) <= 3.00)
          prob=exp(-t*t/2.e0)/sigg(j);
        else if (t > 0.e0)
        {
          dvariable u=t-3.03;
          prob=(u*u*(a2+a3*u))/sigg(j);
        }
        else if (t < 0.e0)
        {
          dvariable u=t+3.03;
          prob=(u*u*(a2-a3*u))/sigg(j);
        }
        tprobb(i) += propp(j) * prob;
      }
      else
      {
        //if (mean_len(j)>fmidd)
          //break_flag=1;
        //else
          //js=j;
        //if (break_flag ==1) break;
      }
      //if (break_flag ==1) break;
    }
    //break_flag=0;
  }
  
  dvariable ssum=0.0;
  for (i=1; i<=nwint; i++) 
  {
    ssum+=tprobb(i);
  }
  for (i=1; i<=nwint; i++) 
  {
    tprobb(i)/=ssum;
  }
}
*/

void dvar_len_fish_stock_history::main_weight_calcs(dvar_vector& tprobb,
  dvar_vector& propp,dvar_vector& mean_len,dvar_vector& sigg)
{
  int i;
  const MY_DOUBLE_TYPE a2 = 35.91908881;
  const MY_DOUBLE_TYPE a3 = 785.858644;
  const MY_DOUBLE_TYPE a22 = 2*35.91908881;
  const MY_DOUBLE_TYPE a33 = 3*785.858644;
  MY_DOUBLE_TYPE prob=0.0;
  int js=1;
  int break_flag=0;
  int ng=propp.indexmax();
  for (i=1; i<=nwint; i++) 
  {
    MY_DOUBLE_TYPE fmidd=0.0;                  //NMD 21Aug2012
    if(!pmsd)
    {
      fmidd=wmid(i);
    }
    else
    {
      if(pmsd->current_species==1)
      {
        fmidd=wmid(i);
      }
      else
      {
        int is=pmsd->current_species;
        fmidd=pmsd->wmid(is,i);
      }
    }                                 //NMD 21Aug2012
   
    for (int j=1; j<=ng; j++) 
    {
      MY_DOUBLE_TYPE t=(fmidd-value(mean_len(j)))/value(sigg(j));
      if ( fabs(t) <= 3.03)
      {
        if ( fabs(t) <= 3.00)
          prob=exp(-t*t/2.e0)/value(sigg(j));
        else if (t > 0.e0)
        {
          MY_DOUBLE_TYPE u=t-3.03;
          prob=(u*u*(a2+a3*u))/value(sigg(j));
        }
        else if (t < 0.e0)
        {
          MY_DOUBLE_TYPE u=t+3.03;
          prob=(u*u*(a2-a3*u))/value(sigg(j));
        }
        value(tprobb(i)) += value(propp(j)) * prob;
      }
      else
      {
        //if (mean_len(j)>fmidd)
          //break_flag=1;
        //else
          //js=j;
        //if (break_flag ==1) break;
      }
      //if (break_flag ==1) break;
    }
    //break_flag=0;
  }
  
  MY_DOUBLE_TYPE ssum=0.0;
  for (i=1; i<=nwint; i++) 
  {
    ssum+=value(tprobb(i));
  }
  for (i=1; i<=nwint; i++) 
  {
    value(tprobb(i))/=ssum;
  }
  
//    save_identifier_string("uv");
  const char * str1;
  str1="uv";
  char* strx1=const_cast <char*> (str1);
  save_identifier_string(strx1);
  tprobb.save_dvar_vector_position();
  propp.save_dvar_vector_value();
  propp.save_dvar_vector_position();
  mean_len.save_dvar_vector_value();
  mean_len.save_dvar_vector_position();

  if(!pmsd)
  {
    wmid.save_dvector_value();
    wmid.save_dvector_position();
  }
  else
  {
    if(pmsd->current_species==1)
    {
      wmid.save_dvector_value();
      wmid.save_dvector_position();
    }
    else
    {
      int is=pmsd->current_species;
      pmsd->wmid(is).save_dvector_value();
      pmsd->wmid(is).save_dvector_position();
    }
  }                                 //NMD 21Aug2012

  sigg.save_dvar_vector_value();
  sigg.save_dvar_vector_position();
  save_double_value(ng);
  save_double_value(nwint);
//    save_identifier_string("yu");
  const char * str2;
  str2="yu";
  char* strx2=const_cast <char*> (str2);
  save_identifier_string(strx2);
  gradient_structure::GRAD_STACK1->
    set_gradient_stack(DF_main_weight_calcs);
}

void DF_main_weight_calcs(void)
{
  verify_identifier_string("yu");
  int i;
  MY_DOUBLE_TYPE nwint=restore_double_value();
  MY_DOUBLE_TYPE nage=restore_double_value();
  dvar_vector_position siggpos=restore_dvar_vector_position();
  dvector sigg=restore_dvar_vector_value(siggpos);
  dvector_position wmidpos=restore_dvector_position();
  dvector wmid=restore_dvector_value(wmidpos);
  dvar_vector_position mean_lenpos=restore_dvar_vector_position();
  dvector mean_len=restore_dvar_vector_value(mean_lenpos);
  dvar_vector_position proppos=restore_dvar_vector_position();
  dvector propp=restore_dvar_vector_value(proppos);
  dvar_vector_position tprobbpos=restore_dvar_vector_position();
  dvector dftprobb=restore_dvar_vector_derivatives(tprobbpos);
  verify_identifier_string("uv");

  dvector tmp;
  dvector dftmp;
  tmp.allocate(dftprobb);
  dftmp.allocate(dftprobb);

  tmp.initialize();
  dftmp.initialize();

  dvector dfpropp;
  dfpropp.allocate(propp);
  dfpropp.initialize();

  dvector dfsigg;
  dvector dfmean_len;
  dfsigg.allocate(sigg);
  dfmean_len.allocate(mean_len);
  dfsigg.initialize();
  dfmean_len.initialize();
  dmatrix prob(1,static_cast<int>(nwint),1,static_cast<int>(nage));
  dmatrix dfprob(1,static_cast<int>(nwint),1,static_cast<int>(nage));
  MY_DOUBLE_TYPE dft=0.0;
  MY_DOUBLE_TYPE dfssum=0.0;
  const MY_DOUBLE_TYPE a2 = 35.91908881;
  const MY_DOUBLE_TYPE a3 = 785.858644;
  const MY_DOUBLE_TYPE a22 = 2*35.91908881;
  const MY_DOUBLE_TYPE a33 = 3*785.858644;
  dfprob.initialize();

  int break_flag=0;
  int js=1;
  
  for (i=1; i<=nwint; i++) 
  {
    MY_DOUBLE_TYPE& fmidd=wmid(i);
    for (int j=1; j<=nage; j++) 
    {
      MY_DOUBLE_TYPE t=(fmidd-mean_len(j))/sigg(j);
      if ( fabs(t) <= 3.03)
      {
        if ( fabs(t) <= 3.00)
          prob(i,j)=exp(-t*t/2.e0)/sigg(j);
        else if (t > 0.e0)
        {
          MY_DOUBLE_TYPE u=t-3.03;
          prob(i,j)=(u*u*(a2+a3*u))/sigg(j);
        }
        else if (t < 0.e0)
        {
          MY_DOUBLE_TYPE u=t+3.03;
          prob(i,j)=(u*u*(a2-a3*u))/sigg(j);
        }
        tmp(i) += propp(j) * prob(i,j);
      }
      else
      {
        //if (mean_len(j)>fmidd)
          //break_flag=1;
        //else
          //js=j;
        //if (break_flag ==1) break;
      }
      //if (break_flag ==1) break;
    }
    //break_flag=0;
  }

  MY_DOUBLE_TYPE ssum=sum(tmp);
  dvector tprob=tmp/ssum;

  // dvector tprob=tmp/ssum;
  for (i=1; i<=nwint; i++) 
  {
    dfssum-=dftprobb(i)*tprob(i)/ssum;
  }  
  dftmp=dftprobb/ssum;

  //double ssum=sum(tmp);
  for (i=1; i<=nwint; i++) 
  {
    dftmp(i)+=dfssum;
  }  
  dfssum=0.0;
 
  //dftmp=dftprobb;

  break_flag=0;
  js=1;
  for (i=1; i<=nwint; i++) 
  {
    MY_DOUBLE_TYPE& fmidd=wmid(i);
    for (int j=1; j<=nage; j++) 
    {
      MY_DOUBLE_TYPE t=(fmidd-mean_len(j))/sigg(j);
      if ( fabs(t) <= 3.03)
      {
        
        if ( fabs(t) <= 3.00)
          prob(i,j)=exp(-t*t/2.e0)/sigg(j);
        else if (t > 0.e0)
        {
          MY_DOUBLE_TYPE u=t-3.03;
          prob(i,j)=(u*u*(a2+a3*u))/sigg(j);
        }
        else if (t < 0.e0)
        {
          MY_DOUBLE_TYPE u=t+3.03;
          prob(i,j)=(u*u*(a2-a3*u))/sigg(j);
        }
        

        // tmp(i) += propp(j) * prob;
        
        dfpropp(j)+=dftmp(i) * prob(i,j);
        dfprob(i,j)+=dftmp(i) * propp(j);

        if ( fabs(t) <= 3.00)
        {
          //prob=exp(-t*t/2.d0)/sigg;

          dft-=dfprob(i,j)*t*prob(i,j);
          dfsigg(j)-=dfprob(i,j)*prob(i,j)/sigg(j);
          dfprob(i,j)=0.0;
        }
        else if (t > 0.e0)
        {
          MY_DOUBLE_TYPE u=t-3.03;
          //prob=(u*u*(a2+a3*u))/sigg(j);
          dft+=dfprob(i,j)*u*(a22+a33*u)/sigg(j);
          dfsigg(j)-=dfprob(i,j)*prob(i,j)/sigg(j);
          dfprob(i,j)=0.0;
        }
        else if (t < 0.e0)
        {
          MY_DOUBLE_TYPE u=t+3.03;
          dft+=dfprob(i,j)*u*(a22-a33*u)/sigg(j);
          dfsigg(j)-=dfprob(i,j)*prob(i,j)/sigg(j);
          dfprob(i,j)=0.0;
        }
        // t=(fmidd-mean_len(j))/sigg(j);
        dfmean_len(j)-=dft/sigg(j);
        dfsigg(j)-=dft*t/sigg(j);
        dft=0.0;
      }
      else
      {
        //if (mean_len(j)>fmidd)
          //break_flag=1;
        //else
          //js=j;
        //if (break_flag ==1) break;
      }
      //if (break_flag ==1) break;
    }
    //break_flag=0;
  }
  dfmean_len.save_dvector_derivatives(mean_lenpos);
  dfsigg.save_dvector_derivatives(siggpos);
  dfpropp.save_dvector_derivatives(proppos);
}

