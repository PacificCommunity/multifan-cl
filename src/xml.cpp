/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


#include "all.hpp"

void dvar_len_fish_stock_history::mean_weights_calc(int mode_flag)
{
  int ir;
  ivector begin_period(1,num_regions);
  ivector end_period(1,num_regions);
  if (mode_flag==0)   // not projection
  {
    begin_period=1;
    end_period=num_real_fish_periods;
  }
  else
  {
    begin_period=num_real_fish_periods+1;
    end_period=num_fish_periods;
  }
    
  dvar_vector ttmp(1,nage);
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int ip=begin_period(ir);ip<=end_period(ir);ip++)  
    {     // Loop over fishing periods
      // Only need this if there are incidents DF 13-07-05
      if (num_fish_incidents(ir,ip)>0)
      {
        dvariable sv27=get_sv_region(ir,27);  //NMD 12Dec2011
        dvariable sv28=get_sv_region(ir,28);  //NMD 12Dec2011
        if (value(sv28)==3.0)   //NMD 12Dec2011
        {
          pmsd_error();
          ttmp=len_wt_coff*
            (pow(mean_length(ir,ip,1),3)+
             3.0*elem_prod(mean_length(ir,ip,1),vars(ir,ip,1)));
        }
        else
        {
          for (int j=1;j<=nage;j++)
          {  
            ttmp(j)=normal_length_to_weight(0.5,-3.5,
              mean_length(ir,ip,1,j),sqrt(vars(ir,ip,1,j)),
              value(sv27),value(sv28));   //NMD 12Dec2011
          }
        }
      }
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)
      {
        mean_weight(ir,ip,fi)=ttmp;
      }
    }
  }
}

  void dvar_len_fish_stock_history::mean_length_calc(int mode_flag)
  {
    int ir;
    MY_DOUBLE_TYPE tpi=2.e0*3.14159e0;
    dvariable tmppow;
    dvariable temp1,temp2,tmp;
    if (!pmsd)
    {
      dvariable rho=exp(-vb_coff(3));
#if !defined(NO_MY_DOUBLE_TYPE)
      temp2=1.-pow(rho,nage-1.0L);
#else
      temp2=1.-pow(rho,nage-1.0);
#endif
      dvar_vector rel_strength;
      if (mode_flag==0)   // not projection
      {
        rel_strength.allocate(1,nyears);
        rel_strength.initialize();
        if (parest_flags(157)>0)
        {
          if (age_flags(92) || af170q0)
          {
            cerr << "Incompatible flags af92 and pf157 or af170q0" << endl;
            ad_exit(1);
          }
          for (int i=1;i<=nyears;i++)
          {
            for (ir=1;ir<=num_regions;ir++)
            {
              rel_strength(i)+=exp(N(ir,i,1));
            }
          }
          rel_strength=rel_strength-mean(rel_strength);
          rel_strength/=sqrt((.001+norm2(rel_strength))/rel_strength.size());
        }
      }
      ivector begin_period(1,num_regions);
      ivector end_period(1,num_regions);
      if (mode_flag==0)   // not projection
      {
        begin_period=1;
        if(have_projection_periods_flag==0)
        {
          end_period=num_fish_periods;
        }
        else
        {
          end_period=num_real_fish_periods;
        }
      }
      else
      {
        begin_period=num_real_fish_periods+1;
        end_period=num_fish_periods;
      }
      
      for (ir=1;ir<=num_regions;ir++)
      {
        for (int ip=begin_period(ir);ip<=end_period(ir);ip++)  
        {     // Loop over fishing periods
          MY_DOUBLE_TYPE xx=-1+(month(ir,ip)-1.)/12.+(week(ir,ip)-1.)/48.;
          for (int j=1;j<=nage;j++)         // Loop over age classes
          {
            tmppow=j+xx;
            if (parest_flags(21))
            {
              tmppow+=sv(1)/tpi*sin(tpi*(true_month(ir,ip)/12.-sv(2)));
            }
            if (mode_flag==0)   // not projection
            {
              if (parest_flags(157)>0 || parest_flags(163)>0)
              {
                int cohort=year(ir,ip)-j+1;
                tmp=0.0;
                    if (parest_flags(157)>0)
                {
                  if (cohort>0)
                  {
                    tmp=sv(3)*rel_strength(cohort);
                    if (parest_flags(158)>0 && cohort >1)
                    {
                      tmp+=sv(4)*rel_strength(cohort-1);
                    }
                    int numyrs= actual_recruit.indexmax();
                    if (parest_flags(159)>0 && cohort < numyrs)
                    {
                      tmp+=sv(5)*rel_strength(cohort+1);
                    }
                  }
                }
                if (parest_flags(163)>0 )
                {
                  if (cohort>0) tmp+=growth_dev(cohort);
                }
                tmp=1./(1.+exp(-tmp));
#if !defined(NO_MY_DOUBLE_TYPE)
                tmp=1.9*(tmp-.5L);
#else
                tmp=1.9*(tmp-.5);
#endif
                tmppow+=tmp;
              }
            }
            if (parest_flags(160)>0 && j==nage)
            {
              tmppow += sv(6);
            }
            if (parest_flags(168)) 
            {
              if (!parest_flags(169)) 
                tmppow=pow(tmppow/nage,sv(16))*nage;
              else if (!parest_flags(170))
              {
                dvariable tmp=tmppow/nage;
#if !defined(NO_MY_DOUBLE_TYPE)
                tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5L))*nage;
#else
                tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5))*nage;
#endif
              }
              else if (!parest_flags(171))
              {
                dvariable tmp=tmppow/nage;
#if !defined(NO_MY_DOUBLE_TYPE)
                tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5L)
#else
                tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5)
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
                   +sv(18)*square(tmp-0.5L))*nage;
#else
                   +sv(18)*square(tmp-0.5))*nage;
#endif
              }
              else if (!parest_flags(172))
              {
                dvariable tmp=tmppow/nage;
#if !defined(NO_MY_DOUBLE_TYPE)
                tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5L)
#else
                tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5)
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
                   +sv(18)*square(tmp-0.5L)
#else
                   +sv(18)*square(tmp-0.5)
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
                   +sv(19)*cube(tmp-0.5L))*nage;
#else
                   +sv(19)*cube(tmp-0.5))*nage;
#endif
              }
              else
              {
                dvariable tmp=tmppow/nage;
#if !defined(NO_MY_DOUBLE_TYPE)
                tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5L)
#else
                tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5)
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
                   +sv(18)*square(tmp-0.5L)
#else
                   +sv(18)*square(tmp-0.5)
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
                   +sv(20)*square(square(tmp-0.5L))
#else
                   +sv(20)*square(square(tmp-0.5))
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
                   +sv(19)*cube(tmp-0.5L))*nage;
#else
                   +sv(19)*cube(tmp-0.5))*nage;
#endif
              }
            }
            if (parest_flags(174))
            {
              dvariable temp4=sv(19)/(sv(18)*2.506628275)*
                           exp(-0.5*square((j-sv(17))/sv(18)));
              rho=exp(-(vb_coff(3)-temp4));
#if !defined(NO_MY_DOUBLE_TYPE)
              temp2=1.-pow(rho,nage-1.0L);
#else
              temp2=1.-pow(rho,nage-1.0);
#endif
            }
            
            temp1=1.-pow(rho,tmppow);
  
            // precompiled code for ...
            //dvariable tt=vb_coff(1)+(vb_coff(2)-vb_coff(1))*
              //    temp1/temp2;
             dvariable tt;
//            length_calc(tt,vb_coff(1),vb_coff(2),temp1,temp2);
            if (parest_flags(226)==0)
            {
              length_calc(tt,vb_coff(1),vb_coff(2),temp1,temp2);
            }
            else
            {
              dvariable T;
              if (parest_flags(226)==1)
                T=-exp(vb_coff(4));
              else
                T=exp(vb_coff(4));
            
              dvariable c1=pow(vb_coff(1),-1.0/T);
              dvariable cN=pow(vb_coff(2),-1.0/T);

              length_calc(tt,c1,cN,temp1,temp2);
              tt=pow(tt,-T);
            }
            
            for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
            {
              mean_length(ir,ip,fi,j)=tt;
            }
          }
        }
      }
      if (parest_flags(173))
      {
        int num=parest_flags(173);
        for (ir=1;ir<=num_regions;ir++) 
        {
          for (int ip=begin_period(ir);ip<=end_period(ir);ip++)  
          {
            for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
            { 
              for (int j=1;j<num;j++) 
              {
                mean_length(ir,ip,fi,j+1)+= age_pars(3,j);
  	    }
            }
          }
        }
      }
      for (ir=1;ir<=num_regions;ir++)
      {
        for (int ip=begin_period(ir);ip<=end_period(ir);ip++)  
        {
          for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
          {                                         // incidents for this period
            if (fish_flags(parent(ir,ip,fi),11) > 0)
            {
              mean_length(ir,ip,fi,1)+=vb_bias(parent(ir,ip,fi))*
                (12.-month(ir,ip))/12.;
            }
          }
        }
      }
    }
    else
    {
      int nrr=pmsd->num_real_regions;
      int numsp=pmsd->num_species;
    
      for (int isp=1;isp<=numsp;isp++)  // loop over species
      {
        dvariable rho;
        dvariable vbc3;
        dvariable vbc1;
        dvariable vbc2;
        dvar_vector ap3=get_age_pars_species(isp,3);
        if (isp==1)
        {
          vbc1=vb_coff(1);
          vbc2=vb_coff(2);
          vbc3=vb_coff(3);
          rho=exp(-vb_coff(3));
        }
        else
        {
          vbc1=pmsd->vb_coff(isp,1);
          vbc2=pmsd->vb_coff(isp,2);
          vbc3=pmsd->vb_coff(isp,3);
          rho=exp(-vbc3);
        }
#if !defined(NO_MY_DOUBLE_TYPE)
        temp2=1.-pow(rho,nage-1.0L);
#else
        temp2=1.-pow(rho,nage-1.0);
#endif
        dvar_vector rel_strength;
        if (mode_flag==0)   // not projection
        {
          rel_strength.allocate(1,nyears);
          rel_strength.initialize();
          if (parest_flags(157)>0)
          {
            if (age_flags(92) || af170q0)
            {
              cerr << "Incompatible flags af92 and pf157 or af170q0" << endl;
              ad_exit(1);
            }
            for (int i=1;i<=nyears;i++)
            {
              for (ir=1;ir<=num_regions;ir++)
              {
                rel_strength(i)+=exp(N(ir,i,1));
              }
            }
            rel_strength=rel_strength-mean(rel_strength);
            rel_strength/=sqrt((.001+norm2(rel_strength))/rel_strength.size());
          }
        }
        ivector begin_period(1,num_regions);
        ivector end_period(1,num_regions);
        if (mode_flag==0)   // not projection
        {
          begin_period=1;
          if(have_projection_periods_flag==0)
          {
            end_period=num_fish_periods;
          }
          else
          {
            end_period=num_real_fish_periods;
          }
        }
        else
        {
          begin_period=num_real_fish_periods+1;
          end_period=num_fish_periods;
        }
        
        int rmin=pmsd->region_bounds(isp,1);
        int rmax=pmsd->region_bounds(isp,2);
        for (ir=rmin;ir<=rmax;ir++)
        {
          for (int ip=begin_period(ir);ip<=end_period(ir);ip++)  
          {     // Loop over fishing periods
            MY_DOUBLE_TYPE xx=-1+(month(ir,ip)-1.)/12.+(week(ir,ip)-1.)/48.;
            for (int j=1;j<=nage;j++)         // Loop over age classes
            {
              tmppow=j+xx;
              if (parest_flags(21))
              {
                tmppow+=sv(1)/tpi*sin(tpi*(true_month(ir,ip)/12.-sv(2)));
              }
              if (mode_flag==0)   // not projection
              {
                if (parest_flags(157)>0 || parest_flags(163)>0)
                {
                  int cohort=year(ir,ip)-j+1;
                  tmp=0.0;
                      if (parest_flags(157)>0)
                  {
                    if (cohort>0)
                    {
                      tmp=sv(3)*rel_strength(cohort);
                      if (parest_flags(158)>0 && cohort >1)
                      {
                        tmp+=sv(4)*rel_strength(cohort-1);
                      }
                      int numyrs= actual_recruit.indexmax();
                      if (parest_flags(159)>0 && cohort < numyrs)
                      {
                        tmp+=sv(5)*rel_strength(cohort+1);
                      }
                    }
                  }
                  if (parest_flags(163)>0 )
                  {
                    if (cohort>0) tmp+=growth_dev(cohort);
                  }
                  tmp=1./(1.+exp(-tmp));
#if !defined(NO_MY_DOUBLE_TYPE)
                  tmp=1.9*(tmp-.5L);
#else
                  tmp=1.9*(tmp-.5);
#endif
                  tmppow+=tmp;
                }
              }
              if (parest_flags(160)>0 && j==nage)
              {
                tmppow += sv(6);
              }
              if (parest_flags(168)) 
              {
                if (!parest_flags(169)) 
                  tmppow=pow(tmppow/nage,sv(16))*nage;
                else if (!parest_flags(170))
                {
                  dvariable tmp=tmppow/nage;
#if !defined(NO_MY_DOUBLE_TYPE)
                  tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5L))*nage;
#else
                  tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5))*nage;
#endif
                }
                else if (!parest_flags(171))
                {
                  dvariable tmp=tmppow/nage;
#if !defined(NO_MY_DOUBLE_TYPE)
                  tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5L)
#else
                  tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5)
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
                     +sv(18)*square(tmp-0.5L))*nage;
#else
                     +sv(18)*square(tmp-0.5))*nage;
#endif
                }
                else if (!parest_flags(172))
                {
                  dvariable tmp=tmppow/nage;
#if !defined(NO_MY_DOUBLE_TYPE)
                  tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5L)
#else
                  tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5)
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
                     +sv(18)*square(tmp-0.5L)
#else
                     +sv(18)*square(tmp-0.5)
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
                     +sv(19)*cube(tmp-0.5L))*nage;
#else
                     +sv(19)*cube(tmp-0.5))*nage;
#endif
                }
                else
                {
                  dvariable tmp=tmppow/nage;
#if !defined(NO_MY_DOUBLE_TYPE)
                  tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5L)
#else
                  tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5)
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
                     +sv(18)*square(tmp-0.5L)
#else
                     +sv(18)*square(tmp-0.5)
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
                     +sv(20)*square(square(tmp-0.5L))
#else
                     +sv(20)*square(square(tmp-0.5))
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
                     +sv(19)*cube(tmp-0.5L))*nage;
#else
                     +sv(19)*cube(tmp-0.5))*nage;
#endif
                }
              }
              if (parest_flags(174))
              {
                dvariable temp4=sv(19)/(sv(18)*2.506628275)*
                             exp(-0.5*square((j-sv(17))/sv(18)));
                rho=exp(-(vbc3-temp4));
#if !defined(NO_MY_DOUBLE_TYPE)
                temp2=1.-pow(rho,nage-1.0L);
#else
                temp2=1.-pow(rho,nage-1.0);
#endif
              }
              
              temp1=1.-pow(rho,tmppow);
    
              // precompiled code for ...
              //dvariable tt=vb_coff(1)+(vb_coff(2)-vb_coff(1))*
                //    temp1/temp2;
              dvariable tt;
//              length_calc(tt,vbc1,vbc2,temp1,temp2);
              if (parest_flags(226)==0)
              {
                length_calc(tt,vb_coff(1),vb_coff(2),temp1,temp2);
              }
              else
              {
                dvariable T;
                if (parest_flags(226)==1)
                  T=-exp(vb_coff(4));
                else
                  T=exp(vb_coff(4));

                dvariable c1=pow(vb_coff(1),-1.0/T);
                dvariable cN=pow(vb_coff(2),-1.0/T);

                length_calc(tt,c1,cN,temp1,temp2);
                tt=pow(tt,-T);
              }
              
              for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
              {
                mean_length(ir,ip,fi,j)=tt;
              }
            }
          }
        }
        if (parest_flags(173))
        {
          int num=parest_flags(173);
          for (ir=rmin;ir<=rmax;ir++) 
          {
            for (int ip=begin_period(ir);ip<=end_period(ir);ip++)  
            {
              for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
              { 
                for (int j=1;j<num;j++) 
                {
                  mean_length(ir,ip,fi,j+1)+= ap3(j);
    	    }
              }
            }
          }
        }
        for (ir=rmin;ir<=rmax;ir++)
        {
          for (int ip=begin_period(ir);ip<=end_period(ir);ip++)  
          {
            for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
            {                                         // incidents for this period
              if (fish_flags(parent(ir,ip,fi),11) > 0)
              {
                mean_length(ir,ip,fi,1)+=vb_bias(parent(ir,ip,fi))*
                  (12.-month(ir,ip))/12.;
              }
            }
          }
        }
      }
    }
  }
 
