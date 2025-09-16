/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"

  extern dvar_vector * psv;

  static int lbsel_flag=0;
  static dvar_matrix * plbsel=NULL;
  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
  dvariable mfexp(dvariable& v1);
  dvar_vector mfexp(dvar_vector& );
  dvector mfexp(const dvector& );
  dvariable age_at_length_calcx(dvariable& v,dvar_vector& gml,int nslots);

  dvariable age_at_length_calc(const prevariable& v,const dvar_vector& vb_coff,
    int nage)
  {
    dvariable rho=exp(-vb_coff(3));
    dvariable tmp= (v-vb_coff(1))/(vb_coff(2)-vb_coff(1));
    dvariable tmp1= 1.-tmp*(1-pow(rho,nage-1));
    if (tmp1<=0.0) tmp1=1.e-20;
    dvariable age= 1.-log(tmp1)/vb_coff(3);
    if (age<1) age=1;
    if (age>nage) age=nage;
    return age;
  }

  dvariable age_at_length_calc(MY_DOUBLE_TYPE v,const dvar_vector& vb_coff,int nage)
  {
    if (vb_coff(4)>1.0e-08)
    {
      cerr << "Entered wrong routine for Richards curve option" << endl;
      ad_exit(1);
    }
    dvariable rho=exp(-vb_coff(3));
    dvariable tmp= (v-vb_coff(1))/(vb_coff(2)-vb_coff(1));
    dvariable tmp1= 1.-tmp*(1-pow(rho,nage-1));
    if (tmp1<=0.0) tmp1=1.e-20;
    dvariable age= 1.-log(tmp1)/vb_coff(3);
    if (age<1) age=1;
    if (age>nage) age=nage;
    return age;
  }

  dvariable age_at_length_calc(MY_DOUBLE_TYPE v,const dvar_vector& vb_coff,int nage,
    const ivector& pf)
  {
    if (pf(226)==0)
    {
      dvariable rho=exp(-vb_coff(3));
      dvariable tmp= (v-vb_coff(1))/(vb_coff(2)-vb_coff(1));
      dvariable tmp1= 1.-tmp*(1-pow(rho,nage-1));
      if (tmp1<=0.0) tmp1=1.e-20;
      dvariable age= 1.-log(tmp1)/vb_coff(3);
      if (age<1) age=1;
      if (age>nage) age=nage;
      return age;
    }
    else
    {
      dvariable T=0.0;
      if (pf(226)==1)
        T=-exp(vb_coff(4));
      else
        T=exp(vb_coff(4));
      dvariable pv=pow(v,-1.0/T);
      dvariable c1=pow(vb_coff(1),-1.0/T);
      dvariable cN=pow(vb_coff(2),-1.0/T);
      dvariable rho=exp(-vb_coff(3));
      dvariable tmp= (pv-c1)/(cN-c1);
      dvariable tmp1= 1.-tmp*(1-pow(rho,nage-1));
      if (tmp1<=0.0) tmp1=1.e-20;
      dvariable age= 1.-log(tmp1)/vb_coff(3);
      if (age<1) age=1;
      if (age>nage) age=nage;
      return age;
    }
  }

  dvariable age_at_length_calc(MY_DOUBLE_TYPE v,const dvariable& rho,
    const dvariable& vbdiff,const dvariable& vbtmp,const dvar_vector& vb_coff,int nage)
  {
    //dvariable rho=exp(-vb_coff(3));
    //dvariable tmp= (v-vb_coff(1))/(vb_coff(2)-vb_coff(1));
    dvariable tmp= (v-vb_coff(1))/vbdiff;
    //dvariable tmp1= 1.-tmp*(1-pow(rho,nage-1));
    dvariable tmp1= 1.-tmp*vbtmp;
    if (tmp1<=0.0) tmp1=1.e-20;
    dvariable age= 1.-log(tmp1)/vb_coff(3);
    if (age<1) age=1;
    if (age>nage) age=nage;
    return age;
  }

  dvariable age_at_length_calc(MY_DOUBLE_TYPE v,const dvariable& rho,
    const dvariable& _vbdiff,const dvariable& _vbtmp,const dvar_vector& vb_coff,int nage,const ivector& pf)
  {
    ADUNCONST(dvariable,vbdiff);   //NMD_jan16-19
    ADUNCONST(dvariable,vbtmp);   //NMD_jan16-19
    if (pf(226)==0)      //NMD_27Sep2018
    {
      //dvariable rho=exp(-vb_coff(3));
      //dvariable tmp= (v-vb_coff(1))/(vb_coff(2)-vb_coff(1));
      dvariable tmp= (v-vb_coff(1))/vbdiff;
      //dvariable tmp1= 1.-tmp*(1-pow(rho,nage-1));
      dvariable tmp1= 1.-tmp*vbtmp;
      if (tmp1<=0.0) tmp1=1.e-20;
      dvariable age= 1.-log(tmp1)/vb_coff(3);
      if (age<1) age=1;
      if (age>nage) age=nage;
      return age;
    }
    else
    {
      dvariable T=0.0;
      if (pf(226)==1)
        T=-exp(vb_coff(4));
      else
        T=exp(vb_coff(4));
      dvariable pv=pow(v,-1.0/T);
      dvariable c1=pow(vb_coff(1),-1.0/T);
      dvariable cN=pow(vb_coff(2),-1.0/T);
      vbdiff=cN-c1;
      vbtmp=1-pow(exp(-vb_coff(3)),nage-1);
      dvariable tmp= (pv-c1)/vbdiff;
      dvariable tmp1= 1.-tmp*vbtmp;
      if (tmp1<=0.0) tmp1=1.e-20;
      dvariable age= 1.-log(tmp1)/vb_coff(3);
      if (age<1) age=1;
      if (age>nage) age=nage;
      return age;
    }                //NMD_27Sep2018
    //dvariable rho=exp(-vb_coff(3));
    //dvariable tmp= (v-vb_coff(1))/(vb_coff(2)-vb_coff(1));
//    dvariable tmp= (v-vb_coff(1))/vbdiff;
    //dvariable tmp1= 1.-tmp*(1-pow(rho,nage-1));
//    dvariable tmp1= 1.-tmp*vbtmp;
//    if (tmp1<=0.0) tmp1=1.e-20;
//    dvariable age= 1.-log(tmp1)/vb_coff(3);
//    if (age<1) age=1;
//    if (age>nage) age=nage;
//    return age;
  }

  dvariable age_at_length_calc(const prevariable& v,const dvar_vector& vb_coff,int nage,
    const dvar_vector& sigma)
  {
    dvariable xl=vb_coff(1)-sigma(1);
    dvariable xu=vb_coff(2)+sigma(nage);
    //dvariable xl=vb_coff(1);
    //dvariable xu=vb_coff(2);
    dvariable age= 1.+(nage-1)*(v-xl)/(xu-xl);
    if (age<1)
    {

      //cerr << "need to deal with differentiability in age_at_length_calc" << endl;
      age=1;
    }
    else if (age>nage)
    {
       // cerr << "need to deal with differentiability in age_at_length_calc" << endl;
	age=nage;
    }
    else if (age<1.1)
    {
      dvariable u=age-1.0;
      u=  u*u*( 20. - 100. *u);
      age= 1.+u;
    }
    else if (age>nage-.1)
    {
      dvariable u=nage-age;
      u=  u*u*( 20. - 100. *u);
      age=nage-u;
    }
    return age;
  }

  dvariable age_at_length_calc(MY_DOUBLE_TYPE v,const dvar_vector& vb_coff,int nage,
    const dvar_vector& sigma)
  {
    dvariable xl=vb_coff(1)-sigma(1);
    dvariable xu=vb_coff(2)+sigma(nage);
    //dvariable xl=vb_coff(1);
    //dvariable xu=vb_coff(2);
    dvariable age= 1.+(nage-1)*(v-xl)/(xu-xl);
    if (age<1)
    {

      //cerr << "need to deal with differentiability in age_at_length_calc" << endl;
      age=1;
    }
    else if (age>nage)
    {
       // cerr << "need to deal with differentiability in age_at_length_calc" << endl;
	age=nage;
    }
    else if (age<1.1)
    {
      dvariable u=age-1.0;
      u=  u*u*( 20. - 100. *u);
      age= 1.+u;
    }
    else if (age>nage-.1)
    {
      dvariable u=nage-age;
      u=  u*u*( 20. - 100. *u);
      age=nage-u;
    }
    return age;
  }


  dvariable daves_kludge(prevariable& x)
  {
    MY_DOUBLE_TYPE cx=value(x);
    int i=static_cast<int>(cx);
#if !defined(NO_MY_DOUBLE_TYPE)
    if (cx-i <= 0.5L)
#else
    if (cx-i <= 0.5)
#endif
    {
      dvariable tmp=x-i;
      dvariable tmp2=tmp*tmp;
      dvariable tmp3=tmp*tmp*tmp;
      return (24*tmp3-64*tmp3*tmp+48*tmp3*tmp2);
    }
    else
    {
      dvariable tmp=1-(x-i);
      dvariable tmp2=tmp*tmp;
      dvariable tmp3=tmp*tmp*tmp;
      return (1.-24*tmp3+64*tmp3*tmp-48*tmp3*tmp2);
    } 
  }

void lbsel_calc(dvar_matrix& lbsel,dvar_fish_stock_history& fsh,
  dvar_vector& _vb_coff,dvar_vector& _var_coff)
{
  dvar_matrix ml(1,fsh.num_fisheries,1,fsh.nage);
  dvar_matrix sgml(1,fsh.num_fisheries,1,fsh.nage);  //NMD_aug24_2018
  MY_DOUBLE_TYPE w0=1;
#if !defined(NO_MY_DOUBLE_TYPE)
  MY_DOUBLE_TYPE w1=exp(-.25*.25*.5L);
#else
  MY_DOUBLE_TYPE w1=exp(-.25*.25*.5);
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
  MY_DOUBLE_TYPE w2=exp(-.5*.5*.5L);
#else
  MY_DOUBLE_TYPE w2=exp(-.5*.5*.5);
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
  MY_DOUBLE_TYPE w3=2.*exp(-.5L);
#else
  MY_DOUBLE_TYPE w3=2.*exp(-.5);
#endif
  MY_DOUBLE_TYPE w4=2.*exp(-1.125);
  MY_DOUBLE_TYPE w5=2.*exp(-2.0);
  MY_DOUBLE_TYPE wsum=w0+2*w1+2*w2+2*w3+2.*w4+2.0*w5;
  w0/=wsum;
  w1/=wsum;
  w2/=wsum;
  w3/=wsum;
  w4/=wsum;
  w5/=wsum;
  int ip;
  MY_DOUBLE_TYPE month;
  int i;
  for (i=1;i<=fsh.num_fisheries;i++)  // Loop over fisheries
  {
    int ir=fsh.realization_region(i,1);
    int isp=1;
    if (fsh.pmsd)
    {
      isp=fsh.pmsd->region_species_pointer(ir);
    }
    dvar_vector vb_coff;
    dvar_vector var_coff;
    int nage=0;
    if (isp==1)
    {
      vb_coff=_vb_coff;
      var_coff=_var_coff;
      nage=fsh.nage;     
    }
    else
    {
      vb_coff=fsh.pmsd->vb_coff(isp);
      var_coff=fsh.pmsd->var_coff(isp);
      nage=fsh.pmsd->nage(isp);     
    }
//    dvar_vector sigma(1,nage);     //NMD_aug24_2018
    dvariable rho=exp(-vb_coff(3));
  
    switch(fsh.fish_flags(i,57))
    {
    case 0:
      if (fsh.fish_flags(i,26))
      {
        ip=fsh.realization_period(i,1);
        //month=(fsh.month(i,ip)-1.)/12.;
        month=0.0;
        if (fsh.parest_flags(226)==0)
        {
          for (int j=1;j<=fsh.nage;j++)
          {
            dvariable tmp=(1.-pow(rho,j-1+month))/(1.-pow(rho,fsh.nage-1));
//            sigma(j)=var_coff(1)*exp(var_coff(2)*(-1+2*tmp));    //NMD_aug24_2018
            sgml(i,j)=var_coff(1)*exp(var_coff(2)*(-1+2*tmp));	    
            ml(i,j)=vb_coff(1)+(vb_coff(2)-vb_coff(1))*tmp;
          }
        }
        else
        {
          dvariable T;
          if (fsh.parest_flags(226)==1)
            T=exp(var_coff(4));   // change sign to make simpler
          else
            T=-exp(var_coff(4));
          dvariable c1=pow(var_coff(1),1.0/T);
          dvariable cN=pow(var_coff(2),1.0/T);
          dvariable diff=cN-c1;
          dvariable xn1inv=1.0/(1.-pow(rho,nage-1));
          for (int j=1;j<=fsh.nage;j++)
          {
            dvariable tmp=(1.-pow(rho,j-1))*xn1inv;
            dvariable tt=c1+diff*tmp;
            ml(i,j)=pow(tt,T);
          }
          dvar_vector scaled_ml=setm11(ml(i));   //scale the ml between -1 and 1
//          sigma=var_coff(1)*exp(var_coff(2)*scaled_ml);  //NMD_aug24_2018
          sgml(i)=var_coff(1)*exp(var_coff(2)*scaled_ml);	  
        }
      }
      if (fsh.parest_flags(173))
      {
        fsh.pmsd_error();
        int num=fsh.parest_flags(173);
        for (int j=1;j<num;j++) 
          ml(i,j+1)+= fsh.age_pars(3,j);
      }
      break;
    case 1:
    case 2:
    case 3:
      break;
    default:
      cerr << " Illegal value for fish_flags(i,57) " << endl;
      ad_exit(1);
      break;
    }
  }  
//  dvar_matrix esel=exp(fsh.selcoff);
  dvar4_array esel=exp(fsh.bs_selcoff);  //NMD_25feb2015
    
  ivector ff71=column(fsh.fish_flags,71);
  ivector ff75=column(fsh.fish_flags,75);
  ivector xff3=column(fsh.fish_flags,3);
  ivector ff74=column(fsh.fish_flags,74);
  for (i=1;i<=fsh.num_fisheries;i++)  // Loop over fisheries
  {
    for (int is=1;is<=ff74(i);is++) //NMD_25feb2015    
    {
      for (int ib=1;ib<=ff71(i)+1;ib++) //NMD_25feb2015    
      {
        int ir=fsh.realization_region(i,1);
        int isp=1;
        if (fsh.pmsd)
        {
          isp=fsh.pmsd->region_species_pointer(ir);
        }
        dvar_vector vb_coff;
        if (isp==1)
        {
          vb_coff=_vb_coff;
        }
        else
        {
          vb_coff=fsh.pmsd->vb_coff(isp);
        }
        // modify ff3 to deal with initial 0 selectivities
        int ff3=fsh.fish_flags(i,3)-ff75(i);
        int j;
        switch(fsh.fish_flags(i,57))
        {
        case 0:
          if (fsh.fish_flags(i,26))
          {
            for (j=ff75(i)+1;j<=ff3;j++)
            {
              dvariable lm0=ml(i,j);
//              dvariable lm1=ml(i,j)-.25*sigma(j);  //NMD_aug24_2018
//              dvariable lp1=ml(i,j)+.25*sigma(j);
//              dvariable lm2=ml(i,j)-.5*sigma(j);
//              dvariable lp2=ml(i,j)+.5*sigma(j);
//              dvariable lm3=ml(i,j)-sigma(j);
//              dvariable lp3=ml(i,j)+sigma(j);
//              dvariable lm4=ml(i,j)-1.5*sigma(j);
//              dvariable lp4=ml(i,j)+1.5*sigma(j);
//              dvariable lm5=ml(i,j)-2.0*sigma(j);
//              dvariable lp5=ml(i,j)+2.0*sigma(j);
              dvariable lm1=ml(i,j)-.25*sgml(i,j);
              dvariable lp1=ml(i,j)+.25*sgml(i,j);
              dvariable lm2=ml(i,j)-.5*sgml(i,j);
              dvariable lp2=ml(i,j)+.5*sgml(i,j);
              dvariable lm3=ml(i,j)-sgml(i,j);
              dvariable lp3=ml(i,j)+sgml(i,j);
              dvariable lm4=ml(i,j)-1.5*sgml(i,j);
              dvariable lp4=ml(i,j)+1.5*sgml(i,j);
              dvariable lm5=ml(i,j)-2.0*sgml(i,j);
              dvariable lp5=ml(i,j)+2.0*sgml(i,j);
        
        
              dvariable alm0;
              dvariable alm1;
              dvariable alp1;   
              dvariable alm2;   
              dvariable alp2;   
              dvariable alm3;   
              dvariable alp3;   
              dvariable alm4;   
              dvariable alp4;   
              dvariable alm5;   
              dvariable alp5;   
  
              switch (fsh.fish_flags(i,26))
              {
              case 1:
                alm0=age_at_length_calc(lm0,vb_coff,fsh.nage);
                alm1=age_at_length_calc(lm1,vb_coff,fsh.nage);
                alp1=age_at_length_calc(lp1,vb_coff,fsh.nage);
                alm2=age_at_length_calc(lm2,vb_coff,fsh.nage);
                alp2=age_at_length_calc(lp2,vb_coff,fsh.nage);
                alm3=age_at_length_calc(lm3,vb_coff,fsh.nage);
                alp3=age_at_length_calc(lp3,vb_coff,fsh.nage);
                alm4=age_at_length_calc(lm4,vb_coff,fsh.nage);
                alp4=age_at_length_calc(lp4,vb_coff,fsh.nage);
                alm5=age_at_length_calc(lm5,vb_coff,fsh.nage);
                alp5=age_at_length_calc(lp5,vb_coff,fsh.nage);
        	break;
              case 2:
                if (!fsh.parest_flags(175))
        	{
//                  alm0=age_at_length_calc(lm0,vb_coff,fsh.nage,sigma);  //NMD_aug24_2018
//                  alm1=age_at_length_calc(lm1,vb_coff,fsh.nage,sigma);
//                  alp1=age_at_length_calc(lp1,vb_coff,fsh.nage,sigma);
//                  alm2=age_at_length_calc(lm2,vb_coff,fsh.nage,sigma);
//                  alp2=age_at_length_calc(lp2,vb_coff,fsh.nage,sigma);
//                  alm3=age_at_length_calc(lm3,vb_coff,fsh.nage,sigma);
//                  alp3=age_at_length_calc(lp3,vb_coff,fsh.nage,sigma);
//                  alm4=age_at_length_calc(lm4,vb_coff,fsh.nage,sigma);
//                  alp4=age_at_length_calc(lp4,vb_coff,fsh.nage,sigma);
//                  alm5=age_at_length_calc(lm5,vb_coff,fsh.nage,sigma);
//                  alp5=age_at_length_calc(lp5,vb_coff,fsh.nage,sigma);

                  alm0=age_at_length_calc(lm0,vb_coff,fsh.nage,sgml(i));
                  alm1=age_at_length_calc(lm1,vb_coff,fsh.nage,sgml(i));
                  alp1=age_at_length_calc(lp1,vb_coff,fsh.nage,sgml(i));
                  alm2=age_at_length_calc(lm2,vb_coff,fsh.nage,sgml(i));
                  alp2=age_at_length_calc(lp2,vb_coff,fsh.nage,sgml(i));
                  alm3=age_at_length_calc(lm3,vb_coff,fsh.nage,sgml(i));
                  alp3=age_at_length_calc(lp3,vb_coff,fsh.nage,sgml(i));
                  alm4=age_at_length_calc(lm4,vb_coff,fsh.nage,sgml(i));
                  alp4=age_at_length_calc(lp4,vb_coff,fsh.nage,sgml(i));
                  alm5=age_at_length_calc(lm5,vb_coff,fsh.nage,sgml(i));
                  alp5=age_at_length_calc(lp5,vb_coff,fsh.nage,sgml(i));
		  
        	}
        	else
    	        {
        	  int nage1=fsh.nage+1;
                  alm0=age_at_length_calcx(lm0,fsh.gml,nage1);
                  alm1=age_at_length_calcx(lm1,fsh.gml,nage1);
                  alp1=age_at_length_calcx(lp1,fsh.gml,nage1);
                  alm2=age_at_length_calcx(lm2,fsh.gml,nage1);
                  alp2=age_at_length_calcx(lp2,fsh.gml,nage1);
                  alm3=age_at_length_calcx(lm3,fsh.gml,nage1);
                  alp3=age_at_length_calcx(lp3,fsh.gml,nage1);
                  alm4=age_at_length_calcx(lm4,fsh.gml,nage1);
                  alp4=age_at_length_calcx(lp4,fsh.gml,nage1);
                  alm5=age_at_length_calcx(lm5,fsh.gml,nage1);
                  alp5=age_at_length_calcx(lp5,fsh.gml,nage1);
                }
        	break;
              default:
                cerr << "Illegal value for fish_flag(" << i << ",26)"
                   << endl;
    
              }
    
              int im0=static_cast<int>(value(alm0));   
              int im1=static_cast<int>(value(alm1));   
              int ip1=static_cast<int>(value(alp1));   
              int im2=static_cast<int>(value(alm2));   
              int ip2=static_cast<int>(value(alp2));   
              int im3=static_cast<int>(value(alm3));   
              int ip3=static_cast<int>(value(alp3));   
              int im4=static_cast<int>(value(alm4));   
              int ip4=static_cast<int>(value(alp4));   
              int im5=static_cast<int>(value(alm5));   
              int ip5=static_cast<int>(value(alp5));   
    
              dvariable b=daves_kludge(alm0-im0);
              
              fsh.bstempsel(i,is,ib,j)=((1-b)*esel(i,is,ib,min(ff3,im0))+
				    b*esel(i,is,ib,min(ff3,im0+1)))*w0;
    
              b=daves_kludge(alm1-im1);
//              lbsel(i,j)+=((1-b)*esel(i,min(ff3,im1))+b*esel(i,min(ff3,im1+1)))*w1;
              fsh.bstempsel(i,is,ib,j)+=((1-b)*esel(i,is,ib,min(ff3,im1))+
                                     b*esel(i,is,ib,min(ff3,im1+1)))*w1;    
    
              b=daves_kludge(alp1-ip1);
//              lbsel(i,j)+=((1-b)*esel(i,min(ff3,ip1))+b*esel(i,min(ff3,ip1+1)))*w1;
              fsh.bstempsel(i,is,ib,j)+=((1-b)*esel(i,is,ib,min(ff3,ip1))+
                                     b*esel(i,is,ib,min(ff3,ip1+1)))*w1;    
        
              b=daves_kludge(alm2-im2);
//              lbsel(i,j)+=((1-b)*esel(i,min(ff3,im2))+b*esel(i,min(ff3,im2+1)))*w2;
              fsh.bstempsel(i,is,ib,j)+=((1-b)*esel(i,is,ib,min(ff3,im2))+
                                     b*esel(i,is,ib,min(ff3,im2+1)))*w2;
    
              b=daves_kludge(alp2-ip2);
//              lbsel(i,j)+=((1-b)*esel(i,min(ff3,ip2))+b*esel(i,min(ff3,ip2+1)))*w2;
              fsh.bstempsel(i,is,ib,j)+=((1-b)*esel(i,is,ib,min(ff3,ip2))+
				     b*esel(i,is,ib,min(ff3,ip2+1)))*w2;
    
              b=daves_kludge(alm3-im3);
//              lbsel(i,j)+=((1-b)*esel(i,min(ff3,im3))+b*esel(i,min(ff3,im3+1)))*w3;
              fsh.bstempsel(i,is,ib,j)+=((1-b)*esel(i,is,ib,min(ff3,im3))+
                                     b*esel(i,is,ib,min(ff3,im3+1)))*w3;
    
              b=daves_kludge(alp3-ip3);
//              lbsel(i,j)+=((1-b)*esel(i,min(ff3,ip3))+b*esel(i,min(ff3,ip3+1)))*w3;
              fsh.bstempsel(i,is,ib,j)+=((1-b)*esel(i,is,ib,min(ff3,ip3))+
                                     b*esel(i,is,ib,min(ff3,ip3+1)))*w3;
         
              b=daves_kludge(alm4-im4);
//              lbsel(i,j)+=((1-b)*esel(i,min(ff3,im4))+b*esel(i,min(ff3,im4+1)))*w4;
              fsh.bstempsel(i,is,ib,j)+=((1-b)*esel(i,is,ib,min(ff3,im4))+
                                     b*esel(i,is,ib,min(ff3,im4+1)))*w4;
    
              b=daves_kludge(alp4-ip4);
//              lbsel(i,j)+=((1-b)*esel(i,min(ff3,ip4))+b*esel(i,min(ff3,ip4+1)))*w4;
              fsh.bstempsel(i,is,ib,j)+=((1-b)*esel(i,is,ib,min(ff3,ip4))+
                                     b*esel(i,is,ib,min(ff3,ip4+1)))*w4;
         
              b=daves_kludge(alm5-im5);
//              lbsel(i,j)+=((1-b)*esel(i,min(ff3,im5))+b*esel(i,min(ff3,im5+1)))*w5;
              fsh.bstempsel(i,is,ib,j)+=((1-b)*esel(i,is,ib,min(ff3,im5))+
                                     b*esel(i,is,ib,min(ff3,im5+1)))*w5;
    
              b=daves_kludge(alp5-ip5);
//              lbsel(i,j)+=((1-b)*esel(i,min(ff3,ip5))+b*esel(i,min(ff3,ip5+1)))*w5;
              fsh.bstempsel(i,is,ib,j)+=((1-b)*esel(i,is,ib,min(ff3,ip5))+
                                     b*esel(i,is,ib,min(ff3,ip5+1)))*w5;
         
            }
            //NMD_26feb2015
            dvar_vector tmpsl(1,fsh.nage);
            tmpsl.initialize();
            tmpsl(ff75(i)+1,xff3(i))=fsh.bstempsel(i,is,ib);
              
            for (j=xff3(i)+1;j<=fsh.nage;j++)
            {
              tmpsl(j)=tmpsl(xff3(i));
       	    }
            tmpsl/=mean(tmpsl);
            tmpsl=log(tmpsl);
            fsh.bstempsel(i,is,ib)=tmpsl(ff75(i)+1,xff3(i));
  //NMD_25feb2015
          }
          break;
        case 1:
          {
            const prevariable & fp9=fsh.fish_pars(9,i);
            dvariable fp10=1.0+fsh.fish_pars(10,i);
            MY_DOUBLE_TYPE jmid=0.5*(1+fsh.nage);
            MY_DOUBLE_TYPE jmid1=1.0-jmid;
            int j;
            for (j=1;j<=fsh.nage;j++)
            {
              lbsel(i,j)= 1.0/(1+pow(19,((j-jmid)/jmid1-fp9)/fp10));
            }
            for (j=1;j<=fsh.nage-1;j++)
            {
              lbsel(i,j)/=lbsel(i,fsh.nage);
            }
            lbsel(i,fsh.nage)=1.0;
            if (fsh.fish_flags(i,71))
            {
              dvar3_array & blbsel=fsh.blbsel;
              if (!allocated(blbsel(i)))
              {
                blbsel(i).allocate(1,fsh.fish_flags(i,71),1,fsh.nage);
              }
              for (int ic=1;ic<=fsh.fish_flags(i,71);ic++)
              {
                dvar_matrix & tmp=blbsel(i);
                const dvariable & fp9=fsh.fish_pars(9,i)+fsh.sel_dev_coffs(i,ic,1);
                dvariable fp10=1.0+fsh.fish_pars(10,i)+fsh.sel_dev_coffs(i,ic,2);
                for (j=1;j<=fsh.nage;j++)
                {
                  tmp(ic,j)= 1.0/(1+pow(19,((j-jmid)/jmid1-fp9)/fp10));
                }
                for (j=1;j<=fsh.nage-1;j++)
                {
                  tmp(ic,j)/=tmp(ic,fsh.nage);
                }
                tmp(ic,fsh.nage)=1.0;
              }
            }

          } 
          break;
        case 2:
          {
            const prevariable & fp9=fsh.fish_pars(9,i);
            dvariable fp10=1.0+fsh.fish_pars(10,i);
            const dvariable & efp11=exp(fsh.fish_pars(11,i));
            MY_DOUBLE_TYPE jmid=0.5*(1+fsh.nage);
            MY_DOUBLE_TYPE jmid1=1.0-jmid;
            for (int j=1;j<=fsh.nage;j++)
            {
              MY_DOUBLE_TYPE fj=(j-jmid)/jmid1;
              if (fj<=value(fp9))
              {
                lbsel(i,j)= pow(2,-square(fj/(fp10*efp11)));
              }
              else
              {
                lbsel(i,j)= pow(2,-square(fj*efp11/fp10));
              }
            }
            if (fsh.fish_flags(i,71))
            {
              dvar3_array & blbsel=fsh.blbsel;
              if (!allocated(blbsel(i)))
              {
                blbsel(i).allocate(1,fsh.fish_flags(i,71),1,fsh.nage);
              }
              for (int ic=1;ic<=fsh.fish_flags(i,71);ic++)
              {
                dvar_matrix & tmp=blbsel(i);
                const dvariable & fp9=fsh.fish_pars(9,i)+fsh.sel_dev_coffs(i,ic,1);
                dvariable fp10=1.0+fsh.fish_pars(10,i)+fsh.sel_dev_coffs(i,ic,2);
                const dvariable & efp11=
                  exp(fsh.fish_pars(11,i)+fsh.sel_dev_coffs(i,ic,3));
                for (int j=1;j<=fsh.nage;j++)
                {
                  MY_DOUBLE_TYPE fj=(j-jmid)/jmid1;
                  if (fj<=value(fp9))
                  {
                    tmp(ic,j)= pow(2,-square(fj/(fp10*efp11)));
                  }
                  else
                  {
                    tmp(ic,j)= pow(2,-square(fj*efp11/fp10));
                  }
                }
              }
            }
          }    
          break;
        case 3:

          break;
        default:
          cerr << " Illegal value for fish_flags(i,57) " << endl;
          ad_exit(1);
          break;
        }
      } //NMD_25feb2015    
    } //NMD_25feb2015    
  }
}

void length_basedsel_calc(dvar_matrix& lengthbsel,dvar_fish_stock_history& fsh,
  dvar_vector& vb_coff,dvar_vector& var_coff,int nlint, dvector& fmid)
{
  int icut=fsh.parest_flags(180);
  dvariable rho=exp(-vb_coff(3));
  dvar_matrix esel=exp(fsh.selcoff);
  
  for (int i=1;i<=fsh.num_fisheries;i++)  // Loop over fisheries
  {
    int ff3=fsh.fish_flags(i,3);
    for (int ii=1;ii<=nlint;ii++)
    {
      dvariable lm1=fmid(ii);
      dvariable alm1=age_at_length_calc(lm1,vb_coff,fsh.nage);   
      int im1=static_cast<int>(value(alm1));   
      dvariable b=daves_kludge(alm1-im1);
      lengthbsel(i,ii)=(1-b)*esel(i,min(ff3,im1))+b*esel(i,min(ff3,im1+1));
    }
  }
}

  dvariable age_at_length_calcx(dvariable& v,dvar_vector& gml,int nslots)
  {
    dvariable age;
    MY_DOUBLE_TYPE ub=nslots;
    MY_DOUBLE_TYPE lb=1;
    if (v<=gml(1))
    {
      age=1;
      return age;
    }

    if (v>=gml(nslots))
    {
      age=nslots;
      return age;
    }

    int n1;
    do
    {
      n1=static_cast<int>((ub-lb+1)/2);
      if (v<gml(static_cast<int>(n1+lb)))
        ub=n1+lb;
      else if (v>gml(static_cast<int>(n1+lb)))
        lb=n1+lb;
      else
      {
        lb=n1+lb;
        ub=n1+lb;
	break;
      }
    }
    while (n1>1);	

    dvariable u=(v-gml(static_cast<int>(lb)))/
      (gml(static_cast<int>(ub))-gml(static_cast<int>(lb)));
    if (u<0.0 || u>1.0)
    {
      cerr << "Arithmetic error in age-at-length calc for " << endl;
      cerr << "non-functional selectivity form with growth devs" << endl;
      cerr << "includes sd(length-at-age) - parameter out of bounds" << endl;
      cerr << " Exiting..." << endl;
      exit(1);
    }

    if (u<0.1)
    {
      u=  u*u*( 20. - 100. *u);
    }
    else if (u>0.9)
    {
      u=1.0-u;
      u=  u*u*( 20. - 100. *u);
      u=1.0-u;
    }
    return lb+u;
  }
#ifdef __MSVC32__
  dvariable age_at_length_calcxx(MY_DOUBLE_TYPE& v,dvar_vector& gml,int nslots)
  {
    dvariable age;
    MY_DOUBLE_TYPE ub=nslots;
    MY_DOUBLE_TYPE lb=1;
    if (v<=gml(1))
    {
      age=1;
      return age;
    }

    if (v>=gml(nslots))
    {
      age=nslots;
      return age;
    }

    int n1;
    do
    {
      n1=static_cast<int>((ub-lb+1)/2);
      if (v<gml(static_cast<int>(n1+lb)))
        ub=n1+lb;
      else if (v>gml(static_cast<int>(n1+lb)))
        lb=n1+lb;
      else
      {
        lb=n1+lb;
        ub=n1+lb;
	break;
      }
    }
    while (n1>1);	

    dvariable u=(v-gml(lb))/(gml(ub)-gml(lb));
#if !defined(NO_MY_DOUBLE_TYPE)
    if (u<0.0 || u>1.0L)
#else
    if (u<0.0 || u>1.0)
#endif
    {
      cerr << "Variable u is out of bounds in lbselclc" << endl;
      exit(1);
    }

    return lb+u;
  }
#else
  dvariable age_at_length_calcxx(const dvariable& v,dvar_vector& gml,int nslots)
  {
    dvariable age;
    MY_DOUBLE_TYPE ub=nslots;
    MY_DOUBLE_TYPE lb=1;
    if (v<=gml(1))
    {
      age=1;
      return age;
    }

    if (v>=gml(nslots))
    {
      age=nslots;
      return age;
    }

    int n1;
    do
    {
      n1=static_cast<int>((ub-lb+1)/2);
      if (v<gml(static_cast<int>(n1+lb)))
        ub=n1+lb;
      else if (v>gml(static_cast<int>(n1+lb)))
        lb=n1+lb;
      else
      {
        lb=n1+lb;
        ub=n1+lb;
	break;
      }
    }
    while (n1>1);	

    dvariable u=(v-gml(static_cast<int>(lb)))/
      (gml(static_cast<int>(ub))-gml(static_cast<int>(lb)));
    if (u<0.0 || u>1.0)
    {
      cerr << "Variable u is out of bounds in lbselclc" << endl;
      exit(1);
    }

    return lb+u;
  }
#endif
