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
 void DF_main_length_calcs(void);
 void DF_main_length_calcs2(void);


void dvar_len_fish_stock_history::fast_pred_frequency_calc(void)
{
  //d3_array tprob(1,num_fish_periods, 1,num_fish_incidents, 1, nlint);
  tprob.initialize();
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

  dvar_vector rnuml(1,nlint);
  for (int ir=1;ir<=num_regions;ir++)
  {
    int np=num_real_fish_periods(ir);
    if (projected_simulated_data_flags[1]==1)
    {
      np=num_fish_periods(ir);
    }
      
    for (int ip=1;ip<=np;ip++)  // Loop over fishing periods
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        if (len_sample_size(ir,ip,fi)>0)
        {
          int pi=parent(ir,ip,fi);
          if (fish_flags(pi,57)==3 && fish_flags(pi,26)==3)
          {
            int mn=month(ir,ip);
            int wk=week(ir,ip);
            dvar_vector& tprobb=tprob(ir,ip,fi);
            dvar_vector& propp=prop(ir,ip,fi);
            main_length_calcs_len_based(tprobb,propp,
              //age_len_incident_sel(ir,ip,fi));
              age_length_fishery_size_dist(mn,wk,pi));
          }
          else
          {
            js=1;
            dvar_vector& tprobb=tprob(ir,ip,fi);
            dvar_vector& propp=prop(ir,ip,fi);
            dvar_vector& mean_len=mean_length(ir,ip,fi);
            dvar_vector& sigg=sdevs(ir,ip,fi);
            main_length_calcs(tprobb,propp,mean_len,sigg);
          }
        }
      }
    }
  }
  ztmp1=diff(ztmp,gradient_structure::totalbytes());
}


void dvar_len_fish_stock_history::main_length_calcs(dvar_vector& tprobb,
  dvar_vector& propp,dvar_vector& mean_len,dvar_vector& sigg)
{
  const MY_DOUBLE_TYPE a2 = 35.91908881;
  const MY_DOUBLE_TYPE a3 = 785.858644;
  const MY_DOUBLE_TYPE a22 = 2*35.91908881;
  const MY_DOUBLE_TYPE a33 = 3*785.858644;
  MY_DOUBLE_TYPE prob=0.0;
  int js=1;
  int break_flag=0;
  int i;
  int ng=propp.indexmax();
  for (i=1; i<=nlint; i++) 
  {
    MY_DOUBLE_TYPE& fmidd=fmid(i);
   
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
  for (i=1; i<=nlint; i++) 
  {
    ssum+=value(tprobb(i));
  }
  for (i=1; i<=nlint; i++) 
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
  sigg.save_dvar_vector_value();
  sigg.save_dvar_vector_position();
  save_double_value(shlen);
  save_double_value(filen);
  save_double_value(ng);
  save_double_value(nlint);
//    save_identifier_string("rs");
  const char * str2;
  str2="rs";
  char* strx2=const_cast <char*> (str2);
  save_identifier_string(strx2);
  gradient_structure::GRAD_STACK1->
    set_gradient_stack(DF_main_length_calcs);
}

void DF_main_length_calcs(void)
{

  verify_identifier_string("rs");
  MY_DOUBLE_TYPE dlint=restore_double_value();
  int nlint=static_cast<int>(dlint);
  MY_DOUBLE_TYPE dng=restore_double_value();
  int ng=static_cast<int>(dng);
  MY_DOUBLE_TYPE filen=restore_double_value();
  MY_DOUBLE_TYPE shlen=restore_double_value();
  dvar_vector_position siggpos=restore_dvar_vector_position();
  dvector sigg=restore_dvar_vector_value(siggpos);
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
  dvector fmid(1,nlint);
  fmid.fill_seqadd(shlen+.5*filen,filen);
  dmatrix prob(1,nlint,1,ng);
  dmatrix dfprob(1,nlint,1,ng);
  MY_DOUBLE_TYPE dft=0.0;
  MY_DOUBLE_TYPE dfssum=0.0;
  const MY_DOUBLE_TYPE a2 = 35.91908881;
  const MY_DOUBLE_TYPE a3 = 785.858644;
  const MY_DOUBLE_TYPE a22 = 2*35.91908881;
  const MY_DOUBLE_TYPE a33 = 3*785.858644;
  dfprob.initialize();

  int break_flag=0;
  int js=1;
  
  int i;
  for (i=1; i<=nlint; i++) 
  {
    MY_DOUBLE_TYPE& fmidd=fmid(i);
    for (int j=1; j<=ng; j++) 
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
  for (i=1; i<=nlint; i++) 
  {
    dfssum-=dftprobb(i)*tprob(i)/ssum;
  }  
  dftmp=dftprobb/ssum;

  //double ssum=sum(tmp);
  for (i=1; i<=nlint; i++) 
  {
    dftmp(i)+=dfssum;
  }  
  dfssum=0.0;
 
  //dftmp=dftprobb;

  break_flag=0;
  js=1;
  for (i=1; i<=nlint; i++) 
  {
    MY_DOUBLE_TYPE& fmidd=fmid(i);
    for (int j=1; j<=ng; j++) 
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

void dvar_len_fish_stock_history::main_length_calcs_len_based
  (dvar_vector& tprobb,dvar_vector& propp,dvar_matrix& alis)
{
  tprobb=propp*alis;
  tprobb/=sum(1.e-20+tprobb);
}

