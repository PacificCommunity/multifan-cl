/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"

  extern dvar_vector * psv;
  extern int _file_version_no;

dvariable bs_sel_pen(dvar_fish_stock_history& fsh)
{
  ivector ff74=column(fsh.fish_flags,74);
  ivector ff71=column(fsh.fish_flags,71);
  ivector ff57=column(fsh.fish_flags,57);
  int i;
  dvariable pen=0.0;
/*
  MY_DOUBLE_TYPE wt=fsh.parest_flags(74);
  for (i=1;i<=fsh.num_fisheries;i++)  // Loop over fisheries
  {
    for (int is=1;is<=ff74(i);is++)
    {
      for (int ib=1;ib<=ff71(i)+1;ib++)
      {
        pen+=wt*norm2(fsh.bs_selcoff(i,is,ib));
      }
    }
  }
*/
  dvector wt(1,fsh.num_fisheries);
  if(fsh.parest_flags(74)>0)
  {
    wt=fsh.parest_flags(74);
  }
  else
  {
    for (i=1;i<=fsh.num_fisheries;i++)
    {
      wt(i)=fsh.fish_flags(i,72);
    }
  }
  for (i=1;i<=fsh.num_fisheries;i++)  // Loop over fisheries
  {
    for (int is=1;is<=ff74(i);is++)
    {
      for (int ib=1;ib<=ff71(i)+1;ib++)
      {
        pen+=wt(i)*norm2(fsh.bs_selcoff(i,is,ib));
      }
    }
  }
  return pen;
}

void logistic_sel_calc(dvar_len_fish_stock_history& fsh)
{
  ivector ff74=column(fsh.fish_flags,74);
  ivector ff75=column(fsh.fish_flags,75);
  ivector ff71=column(fsh.fish_flags,71);
  ivector ff57=column(fsh.fish_flags,57);
  int i;
  for (i=1;i<=fsh.num_fisheries;i++)  // Loop over fisheries
  {
    switch(ff57(i))
    {
    case 1:
      for (int is=1;is<=ff74(i);is++)
      {
        for (int ib=1;ib<=ff71(i)+1;ib++)
        {
          dvariable fp9=fsh.bs_selcoff(i,is,ib,1)/fsh.bl_sel_scaling;
//          dvariable fp10=1.0+fsh.bs_selcoff(i,is,ib,2)/fsh.bl_sel_scaling;
          dvariable fp10;   //NMD_22Nov_2024
          if (fsh.bs_selcoff(i,is,ib,2) < 0.0)
          {
            double tmp=value(-fsh.bs_selcoff(i,is,ib,2)/fsh.bl_sel_scaling);
            if ((1.0-tmp) < 1.0e-06)
            {
              fp10=1.0+
                ((fsh.bs_selcoff(i,is,ib,2)-1.0e-06)/fsh.bl_sel_scaling);
            }
            else
            {
              fp10=1.0+fsh.bs_selcoff(i,is,ib,2)/fsh.bl_sel_scaling;
            }
          }
          else
          {
            fp10=1.0+fsh.bs_selcoff(i,is,ib,2)/fsh.bl_sel_scaling;
          }

          //NMD_18feb2015
          dvar_vector& v=fsh.bstempsel(i,is,ib);
          int j;
          int mmax=v.indexmax();
          MY_DOUBLE_TYPE jmid=0.5*(1+mmax);
          MY_DOUBLE_TYPE jmid1=1.0-jmid;
          //  *******************************************
          //  *******************************************
//          int new_length_flag=1;

//          if (new_length_flag==0)
          if (fsh.age_flags(193) && fsh.pmsd)
          {
            if (fsh.pmsd->num_species !=2)
            {
              cerr << "option for af192 only implemented for two species/sexes"
                   << endl;
              ad_exit(1);
            }
            dvar_vector ts=1.0-2.0*fsh.tlength(i);
            for (j=ff75(i)+1;j<=mmax;j++)
            {
              dvariable ex= (ts(j)-fp9)/fp10;
              v(j)= 1.0/(1+mfexp(ex*log(19.0)));
            }
          }
          else
          {
            for (j=ff75(i)+1;j<=mmax;j++)
            {
              dvariable ex= ((j-jmid)/jmid1-fp9)/fp10;
              //v(j)= 1.0/(1+pow(19,ex));
              v(j)= 1.0/(1+mfexp(ex*log(19.0)));
            }
            for (j=ff75(i)+1;j<=mmax-1;j++)
            {
              v(j)/=(v(mmax)+1.e-12);  //NMD_28Feb.2015
            }
            v(mmax)=1.0;
            if (fsh.fish_flags(i,30)==0)
            {
              v/=mean(v);   //NMD_18feb2015
            }
          }
      //    for (j=ff75(i)+1;j<=mmax-1;j++)
      //    {
      //      v(j)/=(v(mmax)+1.e-12);  //NMD_28Feb.2015
      //    }
      //    v(mmax)=1.0;
      //    if (fsh.fish_flags(i,30)==0)
      //    {
      //      v/=mean(v);   //NMD_18feb2015
      //    }
          v=log(v);
        }
      }
      break;
    case 2:
      for (int is=1;is<=ff74(i);is++)
      {
        for (int ib=1;ib<=ff71(i)+1;ib++)
        {
          dvariable fp9=fsh.bs_selcoff(i,is,ib,1)/fsh.bl_sel_scaling;
          dvariable fp10=1.0+fsh.bs_selcoff(i,is,ib,2)/fsh.bl_sel_scaling;
          dvariable efp11=exp(fsh.bs_selcoff(i,is,ib,3)/fsh.bl_sel_scaling);
          dvar_vector& v=fsh.bstempsel(i,is,ib);
          int mmax=v.indexmax();
          MY_DOUBLE_TYPE jmid=0.5*(1+mmax);
          MY_DOUBLE_TYPE jmid1=1.0-jmid;
          for (int j=ff75(i)+1;j<=mmax;j++)
          {
            MY_DOUBLE_TYPE fj=(j-jmid)/jmid1;
            if (fj<=value(fp9))
            {
              v(j)= pow(2,-square(fj/(fp10*efp11)));
            }
            else
            {
              v(j)= pow(2,-square(fj*efp11/fp10));
            }
          }
          if (fsh.fish_flags(i,30)==0)
          {
            v/=mean(v);   //NMD_18feb2015
          }
          v=log(v);
        }    
      }
      break;
    default:
      break;
    }
  }
}

