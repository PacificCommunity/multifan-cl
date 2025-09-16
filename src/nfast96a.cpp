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
 void DF_main_length_calcs_len_based(void);

void dvar_len_fish_stock_history::fast_pred_frequency_calc_len_based(void)
{
  len_dist.initialize();
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

  // ***********************************************************************
  MY_DOUBLE_TYPE aa=1.;
  //dvariable cutlength=0.67*vb_coff(1)+0.33*vb_coff(2);
  int cutindex=12;
  // ***********************************************************************

  dvar_vector rnuml(1,nlint);
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_fish_periods(ir);ip++)  // Loop over fishing periods
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        dvar_vector& nnum_fish=num_fish(ir,ip);
        dvar_vector& llendist=len_dist(ir,ip);
        if (len_sample_size(ir,ip,fi)>0 || fi==1)
        {
          js=1;
          dvar_vector& tprobb=tprob(ir,ip,fi);
          dvar_vector& propp=prop(ir,ip,fi);
          dvar_vector& mean_len=mean_length(ir,ip,fi);
          dvar_vector& sigg=sdevs(ir,ip,fi);
          main_length_calcs_len_based(tprobb,propp,mean_len,sigg,
            nnum_fish,llendist,fi);
        }
      } 
    }
  }

  ztmp1=diff(ztmp,gradient_structure::totalbytes());
}

void dvar_len_fish_stock_history::main_length_calcs_len_based
(dvar_vector& tprobb,dvar_vector& propp,dvar_vector& mean_len,
   dvar_vector& sigg,dvar_vector& nnum_fish,dvar_vector& llendist,
   int fi)
{
  const MY_DOUBLE_TYPE a2 = 35.91908881;
  const MY_DOUBLE_TYPE a3 = 785.858644;
  const MY_DOUBLE_TYPE a22 = 2*35.91908881;
  const MY_DOUBLE_TYPE a33 = 3*785.858644;
  MY_DOUBLE_TYPE prob=0.0;
  int js=1;
  int break_flag=0;
  int i;
  for (i=1; i<=nlint; i++) 
  {
    MY_DOUBLE_TYPE& fmidd=fmid(i);
   
    for (int j=js; j<=nage; j++) 
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
        if (fi==1)
        {
          value(llendist(i)) += exp(value(nnum_fish(j))) * prob;
        }  
      }
      else
      {
        if (mean_len(j)>fmidd)
          break_flag=1;
        else
          js=j;
        if (break_flag ==1) break;
      }
      if (break_flag ==1) break;
    }
    break_flag=0;
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
/*
  if (fi==1)
  {
    MY_DOUBLE_TYPE ssum=0.0;
    for (i=1; i<=nlint; i++) 
    {
      ssum+=value(llendist(i));
    }
    for (i=1; i<=nlint; i++) 
    {
      value(llendist(i))/=ssum;
    }
  }
 */
 
//  save_identifier_string("uv");
  const char * str3;
  str3="uv";
  char* strx3=const_cast <char*> (str3);
  save_identifier_string(strx3);
  nnum_fish.save_dvar_vector_value();
  nnum_fish.save_dvar_vector_position();
  llendist.save_dvar_vector_position();
  tprobb.save_dvar_vector_position();
  propp.save_dvar_vector_value();
  propp.save_dvar_vector_position();
  mean_len.save_dvar_vector_value();
  mean_len.save_dvar_vector_position();
  sigg.save_dvar_vector_value();
  sigg.save_dvar_vector_position();
  save_double_value(shlen);
  save_double_value(filen);
  save_double_value(nage);
  save_double_value(nlint);
  save_int_value(fi);
//  save_identifier_string("rs");
  const char * str4;
  str4="rs";
  char* strx4=const_cast <char*> (str4);
  save_identifier_string(strx4);

  gradient_structure::GRAD_STACK1->
    set_gradient_stack(DF_main_length_calcs_len_based);
}

void DF_main_length_calcs_len_based(void)
{


  verify_identifier_string("rs");
  int fi=restore_int_value();
  MY_DOUBLE_TYPE dlint=restore_double_value();
  MY_DOUBLE_TYPE dnage=restore_double_value();
  int nlint=static_cast<int>(dlint);
  int nage=static_cast<int>(dnage);
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
  dvar_vector_position lendistpos=restore_dvar_vector_position();
  dvector dflendist=restore_dvar_vector_der_nozero(lendistpos);
  dvar_vector_position num_fishpos=restore_dvar_vector_position();
  dvector num_fish=restore_dvar_vector_value(num_fishpos);
  verify_identifier_string("uv");

/*
  dvector tmp1;
  dvector dftmp1;
  if (fi==1)
  {
    tmp1.allocate(dftprobb);
    dftmp1.allocate(dftprobb);
    tmp1.initialize();
    dftmp1.initialize();
  }
*/

  dvector tmp;
  dvector dftmp;
  tmp.allocate(dftprobb);
  dftmp.allocate(dftprobb);
  tmp.initialize();
  dftmp.initialize();

  dvector dfnum_fish;
  dfnum_fish.allocate(num_fish);
  dfnum_fish.initialize();

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
  dmatrix prob(1,nlint,1,nage);
  dmatrix dfprob(1,nlint,1,nage);
  MY_DOUBLE_TYPE dft=0.0;
  MY_DOUBLE_TYPE dfssum=0.0;
  MY_DOUBLE_TYPE dfssum1=0.0;
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
    for (int j=js; j<=nage; j++) 
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
/*
        if (fi==1)
        {
          tmp1(i) += exp(num_fish(j)) * prob(i,j);
        }  
*/
      }
      else
      {
        if (mean_len(j)>fmidd)
          break_flag=1;
        else
          js=j;
        if (break_flag ==1) break;
      }
      if (break_flag ==1) break;
    }
    break_flag=0;
  }
// ******************************************************************
/*
  if (fi==1)
  {
    MY_DOUBLE_TYPE ssum1=sum(tmp1);
    dvector lendist=tmp1/ssum1;

    // dvector tprob=tmp/ssum;
    for (i=1; i<=nlint; i++) 
    {
      dfssum1-=dflendist(i)*lendist(i)/ssum1;
    }  
    dftmp1=dflendist/ssum1;

    //double ssum=sum(tmp);
    for (i=1; i<=nlint; i++) 
    {
      dftmp1(i)+=dfssum1;
    }  
    dfssum1=0.0;
  }
*/
// ******************************************************************

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
    for (int j=js; j<=nage; j++) 
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
        if (fi==1)
        {
          //value(lendist(i)) += exp(value(num_fish(j))) * prob;
          dfprob(i,j)+=dflendist(i) * exp(num_fish(j));
          dfnum_fish(j)+=dflendist(i) * exp(num_fish(j))* prob(i,j);
        }  

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
        if (mean_len(j)>fmidd)
          break_flag=1;
        else
          js=j;
        if (break_flag ==1) break;
      }
      if (break_flag ==1) break;
    }
    break_flag=0;
  }
  dfmean_len.save_dvector_derivatives(mean_lenpos);
  dfnum_fish.save_dvector_derivatives(num_fishpos);
  dfsigg.save_dvector_derivatives(siggpos);
  dfpropp.save_dvector_derivatives(proppos);
}
