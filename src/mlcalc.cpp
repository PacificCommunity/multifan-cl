/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
#if !defined(ZFEV)
extern dvar_len_fish_stock_history * pcfsh;
#endif

void dvar_len_fish_stock_history::mean_lengths_by_year_calc(void)
{
  if (!pmsd)
  {
    int tmult=1;
    if (age_flags(57)) tmult= age_flags(57);
    MY_DOUBLE_TYPE tpi=2.e0*3.14159e0;
    dvariable tmppow;
    dvariable rho=exp(-vb_coff(3));
    dvariable temp1,temp2,tmp;
    dvar_vector rel_strength(1,nyears);
    rel_strength.initialize();
    int numyrs= actual_recruit.indexmax();
    int ir;
    if (parest_flags(157)>0)
    {
      for (int i=1;i<=nyears;i++)
      {
        for (int ir=1;ir<=num_regions;ir++)
        {
          rel_strength(i)+=exp(N(ir,i,1));
        }
      }
      rel_strength=rel_strength-mean(rel_strength);
      rel_strength/=sqrt((.001+norm2(rel_strength))/rel_strength.size());
    }
    int month_1m=month_1-1;
    //dvar_matrix length_month_yr(1,tmult,1,nage);
    //dvar_matrix sdevs_yr(1,tmult,1,nage);
    dvar_matrix u(1,nage,1,12);
    //dvariable sd=var_coff(1);
    MY_DOUBLE_TYPE v;
    for (int tmonth=1;tmonth<=12;tmonth++)
    {
      for (int j=1;j<=nage;j++)         // Loop over age classes
      {
        if (parest_flags(34) == 0)
        {
          u(j,tmonth)= ( 2.*(j-1.e0+(tmonth-1)/12.)-(nage-1.0))/(nage-1.e0);
        }
        else
        {
          v= j-1.e0+(tmonth-1)/12.;
          u(j,tmonth)= -1.e0+2.e0*(1.0-pow(vb_coff(3),v))/
                 (1.0-pow(vb_coff(3),nage-1.0));
        }
      }
    }
  
    for (int iy=1;iy<=age_flags(57);iy++)  // Loop over fishing periods
    {
      MY_DOUBLE_TYPE xmonth=fmod(month_1m+12.0/tmult*(iy-1),12)+1;
      //for (int j=1;j<=nage;j++)         // Loop over age classes
      {
        dvariable v1=var_coff(1);
        dvariable v2=var_coff(2);
        MY_DOUBLE_TYPE v;
        dvariable u;
        for (int j=1;j<=nage;j++)         // Loop over age classes
        {
          if (parest_flags(34) == 0)
          {
             u= ( 2.*(j-1.e0+(xmonth-1)/12.)-(nage-1.0))/(nage-1.e0);
          }
          else
          {
            v= j-1.e0+(xmonth-1)/12.;
            u= -1.e0+2.e0*(1.0-pow(vb_coff(3),v))/
              (1.0-pow(vb_coff(3),nage-1.0));
          }
          sdevs_yr(iy,j)=v1*exp(v2*u);
        }
      }
    }
  
    temp2=1.-pow(rho,nage-1.0);
  
    for (ir=1;ir<=num_regions;ir++)
    {
      for (int iy=1;iy<=nyears;iy++)  // Loop over fishing periods
      {
        MY_DOUBLE_TYPE xmonth=fmod(month_1m+12.0/tmult*(iy-1),12)+1;
        for (int j=1;j<=nage;j++)         // Loop over age classes
        {
          tmppow=j-1;
          if (parest_flags(21))
          {
            tmppow+=sv(1)/tpi*sin(tpi*(xmonth/12.-sv(2)));
          }
          if (parest_flags(157)>0 || parest_flags(163)>0)
          {
            int cohort=iy-j+1;
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
            tmp=1.9*(tmp-.5);
            tmppow+=tmp;
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
              tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5))*nage;
            }
            else if (!parest_flags(171))
            {
              dvariable tmp=tmppow/nage;
              tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5)
                 +sv(18)*square(tmp-0.5))*nage;
            }
            else if (!parest_flags(172))
            {
              dvariable tmp=tmppow/nage;
              tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5)
                 +sv(18)*square(tmp-0.5)
                 +sv(19)*cube(tmp-0.5))*nage;
            }
            else
            {
              dvariable tmp=tmppow/nage;
              tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5)
                 +sv(18)*square(tmp-0.5)
                 +sv(20)*square(square(tmp-0.5))
                 +sv(19)*cube(tmp-0.5))*nage;
            }
          }
          if (parest_flags(174))
          {
            dvariable temp4=sv(19)/(sv(18)*2.506628275)*
                         exp(-0.5*square((j-sv(17))/sv(18)));
            rho=exp(-(vb_coff(3)-temp4));
            temp2=1.-pow(rho,nage-1.0);
          }
          
          temp1=1.-pow(rho,tmppow);
  
          // precompiled code for ...
          //dvariable tt=vb_coff(1)+(vb_coff(2)-vb_coff(1))*
          //    temp1/temp2;
  
          dvariable tt;
       
  
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
          
          mean_length_yr(ir,iy,j)=tt;
        }
        if (iy<=tmult)
        {
          length_month_yr(iy)=mean_length_yr(ir,iy); 
        }
      }
    }
    if (parest_flags(173))
    {
      int num=parest_flags(173);
      for (ir=1;ir<=num_regions;ir++) 
        for (int iy=1;iy<=nyears;iy++)  // Loop over fishing periods
            for (int j=1;j<num;j++) 
              mean_length_yr(ir,iy,j+1)+= age_pars(3,j);
      
    }
  }
  else
  {
    int nrr=pmsd->num_real_regions;
    int numsp=pmsd->num_species;
    
    for (int isp=1;isp<=numsp;isp++)  // loop over species
    {
      int ng=get_nage_species(isp);
      int tmult=1;
      if (age_flags(57)) tmult= age_flags(57);
      MY_DOUBLE_TYPE tpi=2.e0*3.14159e0;
      dvariable tmppow;
      dvariable rho;
      dvariable vbc1;
      dvariable vbc2;
      dvariable vbc3;
      dvariable vbc4;
      if (isp==1)
      {
        vbc1=vb_coff(1);
        vbc2=vb_coff(2);
        vbc3=vb_coff(3);
//        vbc3=vb_coff(4);
        vbc4=vb_coff(4);    //NMD 7Mar2012
        rho=exp(-vb_coff(3));
      }
      else
      {
        vbc1=pmsd->vb_coff(isp,1);
        vbc2=pmsd->vb_coff(isp,2);
        vbc3=pmsd->vb_coff(isp,3);
//        vbc3=pmsd->vb_coff(isp,4);
        vbc4=pmsd->vb_coff(isp,4);  //NMD 7Mar2012
        rho=exp(-vbc3);
      }
      dvariable temp1,temp2,tmp;
      dvar_vector rel_strength(1,nyears);
      rel_strength.initialize();
      int numyrs= actual_recruit.indexmax();
      int ir;
      if (parest_flags(157)>0)
      {
        for (int i=1;i<=nyears;i++)
        {
          for (int ir=pmsd->region_bounds(isp,1);ir<=pmsd->region_bounds(isp,2);ir++)
          {
            rel_strength(i)+=exp(N(ir,i,1));
          }
        }
        rel_strength=rel_strength-mean(rel_strength);
        rel_strength/=sqrt((.001+norm2(rel_strength))/rel_strength.size());
      }
      int month_1m=month_1-1;
      dvar_matrix u(1,ng,1,12);
      MY_DOUBLE_TYPE v;
      for (int tmonth=1;tmonth<=12;tmonth++)
      {
        for (int j=1;j<=ng;j++)         // Loop over age classes
        {
          if (parest_flags(34) == 0)
          {
            u(j,tmonth)= ( 2.*(j-1.e0+(tmonth-1)/12.)-(ng-1.0))/(ng-1.e0);
          }
          else
          {
            v= j-1.e0+(tmonth-1)/12.;
            u(j,tmonth)= -1.e0+2.e0*(1.0-pow(vbc3,v))/
                   (1.0-pow(vbc3,ng-1.0));
          }
        }
      }
    
      for (int iy=1;iy<=age_flags(57);iy++)  // Loop over fishing periods
      {
        MY_DOUBLE_TYPE xmonth=fmod(month_1m+12.0/tmult*(iy-1),12)+1;
        //for (int j=1;j<=ng;j++)         // Loop over age classes
        {
          dvariable v1=get_var_coff_species(isp)(1);
          dvariable v2=get_var_coff_species(isp)(2);
          MY_DOUBLE_TYPE v;
          dvariable u;
          for (int j=1;j<=ng;j++)         // Loop over age classes
          {
            if (parest_flags(34) == 0)
            {
               u= ( 2.*(j-1.e0+(xmonth-1)/12.)-(ng-1.0))/(ng-1.e0);
            }
            else
            {
              v= j-1.e0+(xmonth-1)/12.;
              u= -1.e0+2.e0*(1.0-pow(vbc3,v))/
                (1.0-pow(vbc3,ng-1.0));
            }
            sdevs_yr(iy,j)=v1*exp(v2*u);
          }
        }
      }
    
      temp2=1.-pow(rho,ng-1.0);
    
      //      for (ir=1;ir<=num_regions;ir++)
      for (int ir=pmsd->region_bounds(isp,1);ir<=pmsd->region_bounds(isp,2);ir++) //NMD16Apr2012
      {
        for (int iy=1;iy<=nyears;iy++)  // Loop over fishing periods
        {
          MY_DOUBLE_TYPE xmonth=fmod(month_1m+12.0/tmult*(iy-1),12)+1;
          for (int j=1;j<=ng;j++)         // Loop over age classes
          {
            tmppow=j-1;
            if (parest_flags(21))
            {
              tmppow+=sv(1)/tpi*sin(tpi*(xmonth/12.-sv(2)));
            }
            if (parest_flags(157)>0 || parest_flags(163)>0)
            {
              int cohort=iy-j+1;
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
              tmp=1.9*(tmp-.5);
              tmppow+=tmp;
            }
            if (parest_flags(160)>0 && j==ng)
            {
              tmppow += sv(6);
            }
            if (parest_flags(168)) 
            {
              if (!parest_flags(169)) 
                tmppow=pow(tmppow/ng,sv(16))*ng;
              else if (!parest_flags(170))
              {
                dvariable tmp=tmppow/ng;
                tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5))*ng;
              }
              else if (!parest_flags(171))
              {
                dvariable tmp=tmppow/ng;
                tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5)
                   +sv(18)*square(tmp-0.5))*ng;
              }
              else if (!parest_flags(172))
              {
                dvariable tmp=tmppow/ng;
                tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5)
                   +sv(18)*square(tmp-0.5)
                   +sv(19)*cube(tmp-0.5))*ng;
              }
              else
              {
                dvariable tmp=tmppow/ng;
                tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5)
                   +sv(18)*square(tmp-0.5)
                   +sv(20)*square(square(tmp-0.5))
                   +sv(19)*cube(tmp-0.5))*ng;
              }
            }
            if (parest_flags(174))
            {
              dvariable temp4=sv(19)/(sv(18)*2.506628275)*
                           exp(-0.5*square((j-sv(17))/sv(18)));
              rho=exp(-(vbc3-temp4));
              temp2=1.-pow(rho,ng-1.0);
            }
            
            temp1=1.-pow(rho,tmppow);
    
            // precompiled code for ...
            //dvariable tt=vb_coff(1)+(vb_coff(2)-vb_coff(1))*
            //    temp1/temp2;
    
            dvariable tt;
         
    
            if (parest_flags(226)==0)
            {
              length_calc(tt,vbc1,vbc2,temp1,temp2);
            }
            else
            {
              dvariable T;
              if (parest_flags(226)==1)
                T=-exp(vbc4);
              else
                T=exp(vbc4);
              
              dvariable c1=pow(vbc1,-1.0/T);
              dvariable cN=pow(vbc2,-1.0/T);
      
              length_calc(tt,c1,cN,temp1,temp2);
              tt=pow(tt,-T);
            }
            
            mean_length_yr(ir,iy,j)=tt;
          }

          if (iy<=tmult)
          {
            int is=pmsd->region_species_pointer(ir);
            if (is==1)
              length_month_yr(iy)=mean_length_yr(ir,iy); 
            else
              pmsd->length_month_yr(is,iy)=mean_length_yr(ir,iy); 
          }
        }
      }
      if (parest_flags(173))
      {
        dvar_vector ap3=get_age_pars_species(isp,3);
        int num=parest_flags(173);
	//        for (ir=1;ir<=num_regions;ir++) {
        for (int ir=pmsd->region_bounds(isp,1);ir<=pmsd->region_bounds(isp,2);ir++) {  //NMD16Apr2012
          for (int iy=1;iy<=nyears;iy++) { // Loop over fishing periods
	    for (int j=1;j<num;j++) {
                mean_length_yr(ir,iy,j+1)+= ap3(j);
		//		cout << " ir " << ir << " iy " << iy << "  j " << j << " mn_len_age " << 
		//	    mean_length_yr(ir,iy,j) << "  ap3 " << ap3(j) << endl; //NMD12Apr2012
            }
          }
        }
      }
    }
  }
}

void dvar_len_fish_stock_history::mean_weights_by_year_calc(void)
{
  //  cout << " Debug mean_length_yr " << endl;  //NMD12Apr2012
  if (parest_flags(380)==2)
  {
    //mean_weight_yr.initialize();
    ifstream ifs("mean_weight_yr");
    ifs >> mean_weight_yr_alternative;
    if (!ifs)
    {
      cerr << "File read of file mean_weight_yr error " << endl;
      ad_exit(1);
    }
  }
  for (int ir=1;ir<=num_regions;ir++)
  {
    dvar_vector vvar=get_global_vars_region(ir);
    dvariable sv27=get_sv_region(ir,27);
    dvariable sv28=get_sv_region(ir,28);
    dvariable lwc=get_len_wt_coff_region(ir);
    for (int iy=1;iy<=nyears;iy++)  // Loop over fishing periods
    {
      if (value(sv28)==3.0)
        mean_weight_yr(ir,iy)= lwc*(cube(mean_length_yr(ir,iy))
          +3.0*elem_prod(mean_length_yr(ir,iy),vvar));
      else
      {
        int ng=get_nage_region(ir);
        for (int j=1;j<=ng;j++)
        {  
          mean_weight_yr(ir,iy,j)= normal_length_to_weight(0.5,-3.5,
            mean_length_yr(ir,iy,j),sqrt(vvar(j)),value(sv27),value(sv28));
        }
      }
      //      cout << " ir " << ir << " iy " << iy << " mn_wt_age " << 
      //	mean_weight_yr(ir,iy) << " sv27 " << value(sv27) << 
      //	"  sv28 " << value(sv28) << endl; //NMD12Apr2012
    }
  }
  if (parest_flags(380)==1)
  {
    ofstream ofs("mean_weight_yr");
    mean_weight_yr_alternative=value(mean_weight_yr);
    ofs << mean_weight_yr_alternative;
    parest_flags(380)=2;
  }
}

dvariable normal_length_to_weight(MY_DOUBLE_TYPE delta,MY_DOUBLE_TYPE xmin,
  const prevariable& _mean,const prevariable & _sigma,const prevariable& _a,
  const prevariable & _b)
{
  MY_DOUBLE_TYPE eps=.01;
  MY_DOUBLE_TYPE mean=value(_mean);
  MY_DOUBLE_TYPE sigma=value(_sigma);
  MY_DOUBLE_TYPE b=value(_b);
  MY_DOUBLE_TYPE a=value(_a);
  int n= int((-2*xmin)/delta)+1;
  dvector x(1,n);
  dvector w(1,n);
  dvector y(1,n);
  dvector z(1,n);
  x[1]=xmin-delta;
  MY_DOUBLE_TYPE f=0;
  int i;
  for (i=1;i<=n;i++)
  {
    w[i]=exp(-0.5*x(i)*x(i));
    y[i]=sigma*x(i)+mean;
    if (y[i]<eps) {
      z[i]=eps/(2-y[i]/eps);
      f+=pow(z[i],b)*w[i];
    }
    else
    {
      f+=pow(y[i],b)*w[i];
    }
    if (i<n) x[i+1]=x[i]+delta;
  }
  MY_DOUBLE_TYPE f1=a*f*0.39894*delta;
  dvariable tmp;
  value(tmp)=f1;
  dvector dfy(1,n);
  dvector dfz(1,n);
  dfy.initialize();
  dfz.initialize();
  MY_DOUBLE_TYPE dfa,dfb,dfsigma,dfmean,df1,dff;
  dfsigma=0.0;
  dfmean=0.0;
  // adjoint code
  df1=1.0;
  //double f1=f*0.39894*delta;
  dff=df1*a*0.39894*delta;
  dfa=df1*f*0.39894*delta;
  for (i=n;i>=1;i--)
  {
    if (y[i]<eps) {
      //f+=pow(z[i],b)*w[i];
      dfz[i]+=dff*b*pow(z[i],b-1.0)*w[i];
      dfb+=dff*pow(z[i],b)*log(z[i])*w[i];
      //z[i]=eps/(2-y[i]/eps);
      dfy[i]+=dfz[i]/square(2-y[i]/eps);
    }
    else
    {
      //f+=pow(y[i],b)*w[i];
      dfy[i]+=dff*b*pow(y[i],b-1.0)*w[i];
      dfb+=dff*pow(y[i],b)*log(y[i])*w[i];
    }
    //y[i]=sigma*x(i)+mean;
    dfsigma+=dfy[i]*x(i);
    dfmean+=dfy[i];
  }
  gradient_structure::GRAD_STACK1->set_gradient_stack(default_evaluation4ind,
    address(tmp),
    address(_mean),dfmean ,
    address(_sigma),dfsigma ,
    address(_a),dfa,
    address(_b),dfb);

  return tmp;
}


dvariable normal_length_to_weight(MY_DOUBLE_TYPE delta,MY_DOUBLE_TYPE xmin,
  const prevariable& _mean,const prevariable & _sigma,MY_DOUBLE_TYPE a,
  MY_DOUBLE_TYPE & b)
{
  MY_DOUBLE_TYPE eps=.01;
  MY_DOUBLE_TYPE mean=value(_mean);
  MY_DOUBLE_TYPE sigma=value(_sigma);
  int n= int((-2*xmin)/delta)+1;
  dvector x(1,n);
  dvector w(1,n);
  dvector y(1,n);
  dvector z(1,n);
  x[1]=xmin-delta;
  MY_DOUBLE_TYPE f=0;
  int i;
  for (i=1;i<=n;i++)
  {
    w[i]=exp(-0.5*x(i)*x(i));
    y[i]=sigma*x(i)+mean;
    if (y[i]<eps) {
      z[i]=eps/(2-y[i]/eps);
      f+=pow(z[i],b)*w[i];
    }
    else
    {
      f+=pow(y[i],b)*w[i];
    }
    if (i<n) x[i+1]=x[i]+delta;
  }
  MY_DOUBLE_TYPE f1=a*f*0.39894*delta;
  dvariable tmp;
  value(tmp)=f1;
  dvector dfy(1,n);
  dvector dfz(1,n);
  dfy.initialize();
  dfz.initialize();
  MY_DOUBLE_TYPE dfsigma,dfmean,df1,dff;
  dfsigma=0.0;
  dfmean=0.0;
  // adjoint code
  df1=1.0;
  //double f1=f*0.39894*delta;
  dff=df1*a*0.39894*delta;
  for (i=n;i>=1;i--)
  {
    if (y[i]<eps) {
      //f+=pow(z[i],b)*w[i];
      dfz[i]+=dff*b*pow(z[i],b-1.0)*w[i];
      //z[i]=eps/(2-y[i]/eps);
      dfy[i]+=dfz[i]/square(2-y[i]/eps);
    }
    else
    {
      //f+=pow(y[i],b)*w[i];
      dfy[i]+=dff*b*pow(y[i],b-1.0)*w[i];
    }
    //y[i]=sigma*x(i)+mean;
    dfsigma+=dfy[i]*x(i);
    dfmean+=dfy[i];
  }
  gradient_structure::GRAD_STACK1->set_gradient_stack(default_evaluation,
    address(tmp),
    address(_mean),dfmean ,
    address(_sigma),dfsigma);

  return tmp;
}

void dvar_len_fish_stock_history::mean_lengths_by_year_for_projections(void)
{
  pmsd_error();
  int iy;
  int tmult=1;
  if (age_flags(57)) tmult= age_flags(57);
  MY_DOUBLE_TYPE tpi=2.e0*3.14159e0;
  dvariable tmppow;
  dvariable rho=exp(-vb_coff(3));
  dvariable temp1,temp2,tmp;
  temp2=1.-pow(rho,nage-1.0);
  int month_1m=month_1-1;
  ivector ttmp(1,nyears);
  for (iy=1;iy<=tmult;iy++)  // Loop over fishing periods
  {
    ttmp(iy)=fmod(month_1m+12.0/tmult*(iy-1),12)+1.0001;
  }
  for (iy=1;iy<=tmult;iy++)  // Loop over fishing periods
  {
    MY_DOUBLE_TYPE xmonth=fmod(month_1m+12.0/tmult*(iy-1),12)+1.0001;
    for (int j=1;j<=nage;j++)         // Loop over age classes
    {
      tmppow=j-1;
      if (parest_flags(21))
      {
        tmppow+=sv(1)/tpi*sin(tpi*(xmonth/12.-sv(2)));
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
          tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5))*nage;
        }
        else if (!parest_flags(171))
        {
          dvariable tmp=tmppow/nage;
          tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5)
             +sv(18)*square(tmp-0.5))*nage;
        }
        else if (!parest_flags(172))
        {
          dvariable tmp=tmppow/nage;
          tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5)
             +sv(18)*square(tmp-0.5)
             +sv(19)*cube(tmp-0.5))*nage;
        }
        else
        {
          dvariable tmp=tmppow/nage;
          tmppow=pow(tmp,sv(16)+sv(17)*(tmp-0.5)
             +sv(18)*square(tmp-0.5)
             +sv(20)*square(square(tmp-0.5))
             +sv(19)*cube(tmp-0.5))*nage;
        }
      }
      if (parest_flags(174))
      {
        dvariable temp4=sv(19)/(sv(18)*2.506628275)*
                     exp(-0.5*square((j-sv(17))/sv(18)));
        rho=exp(-(vb_coff(3)-temp4));
        temp2=1.-pow(rho,nage-1.0);
      }
      
      temp1=1.-pow(rho,tmppow);

      // precompiled code for ...
      //dvariable tt=vb_coff(1)+(vb_coff(2)-vb_coff(1))*
      //    temp1/temp2;
      dvariable tt;
      length_calc(tt,vb_coff(1),vb_coff(2),temp1,temp2);
      
      mean_length_yr_proj(iy,j)=tt;
    }
  }
  if (parest_flags(173))
  {
    int num=parest_flags(173);
    for (iy=1;iy<=tmult;iy++)  // Loop over fishing periods
      for (int j=1;j<num;j++) 
        mean_length_yr(iy,j+1)+= age_pars(3,j);
  }
}

/*       //NMD 12Dec2011 - will need a call to get_sv_region or get_sv_species to work properly
void dvar_len_fish_stock_history::mean_weights_by_year_for_projections(void)
{
  int tmult=1;
  if (age_flags(57)) tmult= age_flags(57);
  dvar_vector vvar=global_vars;
  for (int iy=1;iy<=tmult;iy++)  // Loop over fishing periods
  {
    if (sv(28)==3.0)
      mean_weight_yr_proj(iy)= len_wt_coff*(cube(mean_length_yr_proj(iy))
        +3.0*elem_prod(mean_length_yr_proj(iy),vvar));
    else
    {
      for (int j=1;j<=nage;j++)
      {  
        mean_weight_yr_proj(iy,j)= normal_length_to_weight(0.5,-3.5,
          mean_length_yr_proj(iy,j),sqrt(vvar(j)),value(sv(27)),value(sv(28)));
      }
    }
  }
}
*/   //NMD 12Dec2011


#undef HOME_VERSION

