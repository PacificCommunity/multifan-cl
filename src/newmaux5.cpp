/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"
#include "scbd.hpp"
void xinit_message(ofstream& xof,int iisave,int ii,const char * msg);

int get_max_bias(_CONST ivector& common_vb_bias,int num_fisheries);

void set_value_inv(dvector& w,dvector& x,int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  ivector& flags,ivector& group);

int num_active(dvector& w,ivector& flags,ivector& group);


static int * ppii = NULL;

dvariable dvar_len_fish_stock_history::reset(dvar_vector& x)
{
  if (parest_flags(385)==1)
  {
    x=elem_div(x,sqrt(pos_hess_scale));
  }
  else if (parest_flags(385)==2)
  {
    x=psqrt_hess_inv*x;
  }
  int ii=1;
  ppii = &ii;
  dvariable pen=0.;

  if (parest_flags(392)==0)
  {
    if (age_flags(30)>0)
    {
      if (parest_flags(155)!=0)
      {
        recinpop_orth_reset(x,ii);
      }
    }
  }

  if (parest_flags(361)>0)
  {
    set_value(testpar,x,ii);
  }

  pen=dvar_fish_stock_history::reset(x,ii,len_sample_size);
  
  if (parest_flags(392)==0)
  {

//    pen=dvar_fish_stock_history::reset(x,ii,len_sample_size);
  
    // von Betalanffy parameters
    //cout << "In reset delta2 =" << endl << delta2 << endl;
    if (!pmsd || sum(column(pmsd->species_flags,2))==0)
    {
      if (age_flags(146)>0  && age_flags(94) !=3)
      {
        set_value(sv(21),x,ii,1.e-6,1000.0,pen,1000.);
        if (pmsd)
        {
          int numsp=pmsd->num_species;
          for (int i=2;i<=numsp;i++)
          {
            set_value(pmsd->sv(i,21),x,ii,1.e-6,1000.0,pen,1000.);
          }
        }
      }
      if (age_flags(162)>0)
      {
        set_value(sv(29),x,ii,0.201,0.999,pen);
        if (pmsd)
        {
          int numsp=pmsd->num_species;
          for (int i=2;i<=numsp;i++)
          {
            set_value(pmsd->sv(i,29),x,ii,0.201,0.999,pen);
          }
        }
      }
    }
    else
    {
      ivector sf2=column(pmsd->species_flags,2);
      if (age_flags(146)>0  && age_flags(94) !=3)
      {
        if (sf2(1))
          set_value(sv(21),x,ii,1.e-6,1000.0,pen,1000.);
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          if (sf2(i))
            set_value(pmsd->sv(i,21),x,ii,1.e-6,1000.0,pen,1000.);
        }
      }
      if (age_flags(162)>0)
      {
        if (sf2(1))
          set_value(sv(29),x,ii,0.201,0.999,pen);
  
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          if (sf2(i))
            set_value(pmsd->sv(i,29),x,ii,0.201,0.999,pen);
        }
      }
    }
    if (age_flags(151)>0)
    {
      set_value(sv(22),x,ii,-.19,10.,pen);
    }
  
    if (age_flags(152)>0)
    {
      set_value(sv(23),x,ii,-.8,5.,pen);
    }
  
    if (age_flags(152)>0)
    {
      set_value(sv(24),x,ii,-0.9,10.0,pen);
    }
  
    set_value(vb_bias,x,ii,0.0,fmax1-fmin1,pen,
      column(fish_flags,11),column(fish_flags,22));
  
    if (age_flags(58)>0)
    {
      MY_DOUBLE_TYPE bnd=age_flags(58)/10.;
      sv(10)=boundp(x(ii++),0.0,bnd,pen,10.);
    }
   
    if (parest_flags(174))
    {
      sv(17)=boundp(x(ii++),0.0,20,pen,10.);
      sv(18)=boundp(x(ii++),0.0,5,pen,10.);
      sv(19)=boundp(x(ii++),-1,1,pen,10.);
    }
  
    if (parest_flags(168)>0)
    {
      sv(16)=boundp(x(ii++),.2,5.0,pen,10.);
      //sv(16)=boundp(x(ii++),-.4,.4,pen,10.);
    }
    if (parest_flags(169)>0)
    {
      sv(17)=boundp(x(ii++),-1.9,1.9,pen,10.);
      //sv(17)=boundp(x(ii++),-.4,.4,pen,10.);
    }
   
    if (parest_flags(170)>0)
    {
      sv(18)=boundp(x(ii++),-6.9,6.9,pen,10.);
      //sv(17)=boundp(x(ii++),-.4,.4,pen,10.);
    }
   
    if (parest_flags(171)>0)
    {
      sv(19)=boundp(x(ii++),-6.9,6.9,pen,10.);
      //sv(17)=boundp(x(ii++),-.4,.4,pen,10.);
    }
   
    if (parest_flags(172)>0)
    {
      sv(20)=boundp(x(ii++),-6.9,6.9,pen,10.);
      //sv(17)=boundp(x(ii++),-.4,.4,pen,10.);
    }
   
    if (age_flags(59)>0)
    {
      sv(11)=boundp(x(ii++),-6.,6.,pen,10.);
      sv(12)=boundp(x(ii++),-6.,6.,pen,10.);
    }
  
    if (age_flags(63)>0)
    {
      sv(13)=boundp(x(ii++),-3.,3.,pen,1000.);
    }
  
    if (age_flags(65)>0)
    {
      sv(14)=boundp(x(ii++),-5.,5.,pen,1000.);
    }
  
    if (age_flags(66)>0)
    {
      sv(15)=boundp(x(ii++),-2.,2.,pen,100000.);
      if (pmsd)
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          pmsd->sv(i,15)=boundp(x(ii++),-2.,2.,pen,100000.);
        }
      }
    }
  
    if (parest_flags(157)>0)
    {
      sv(3)=boundp(x(ii),-5.,5.,pen,100.);
      ii=ii+1;
    }
  
    if (parest_flags(158)>0)
    {
      sv(4)=boundp(x(ii),-5.,5.,pen,100.);
      ii=ii+1;
    }
  
    if (parest_flags(159)>0)
    {
      sv(5)=boundp(x(ii),-5.,5.,pen,100.);
      ii=ii+1;
    }
  
    if (age_flags(102)>0)
    {
      set_value(sv(25),x,ii,-10.,10.,pen);
    }
  
    if (age_flags(153)>0)
    {
      set_value(sv(26),x,ii,-2.0,2.0,pen);
      if (pmsd)                             //NMD 16Aug2012
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          set_value(pmsd->sv(i,26),x,ii,-2.0,2.0,pen);
        }
      }                                     //NMD 16Aug2012
    }
  
    if (age_flags(67)>0)
    {
  #if !defined(NO_MY_DOUBLE_TYPE)
      set_value(age_pars(1),x,ii,.0,1.0L,pen);
  #else
      set_value(age_pars(1),x,ii,.0,1.0,pen);
  #endif
    }
  
    if (parest_flags(175)>1)
    {
      MY_DOUBLE_TYPE onet=1000.;
      MY_DOUBLE_TYPE fivet=5000.;
      set_value(age_pars(4),x,ii,0.001,100.,pen,fivet);
      set_value(sv(20),x,ii,0.001,100.0,pen,onet);
    }
  
    if (age_flags(73)>0)
    {
      MY_DOUBLE_TYPE u=5000.;
      if (age_flags(81)<1)
      {
        set_value(age_pars(2),x,ii,-5.0,5.0,pen,u);
      }
      else
      {
        int num=nage-age_flags(81)+1;
        set_value_partial(age_pars(2),x,ii,num,-5.0,5.0,pen,5000.0);
      }
    }
  
    if (parest_flags(184)>0)
    {
      int num=parest_flags(173)-1;
      set_value_partial(age_pars(3),x,ii,num,-50.0,50.0,pen,1000.0);
      if (pmsd)
      {
        int numsp=pmsd->num_species;
        for (int is=2;is<=numsp;is++)
        {
          set_value_partial(get_age_pars_species(is,3),x,ii,num,
            -50.0,50.0,pen,1000.0);
        }
      }
    }
    if (parest_flags(121)>0)
    {
      int num=parest_flags(121);
      MY_DOUBLE_TYPE nscal=1.0;
      if (parest_flags(387))  //NMD_4sep2020
      {
        nscal=4500.0;
      }
      else
      {
        nscal=200.0;
      }
  //    set_value_partial(age_pars(5),x,ii,num,-20.0,2.0,pen,4500.0);
      set_value_partial(age_pars(5),x,ii,num,-20.0,2.0,pen,nscal);
  //NMD_4sep2020
      if (pmsd)
      {
        int numsp=pmsd->num_species;
        for (int is=2;is<=numsp;is++)
        {
          set_value_partial(get_age_pars_species(is,5),x,ii,num,
            -20.0,2.0,pen,nscal);      //NMD_4sep2020
  //          -20.0,2.0,pen,4500.0);
        }
      }
    }
  
    if (parest_flags(160)>0)
    {
      sv(6)=boundp(x(ii),-0.1,5.,pen,100.);
      ii=ii+1;
    }
    if (parest_flags(163)>0)
    {
      set_value(growth_dev,x,ii,-5.0,5.0,pen);
    }
  
    if (parest_flags(12)>0)
    {

      vb_coff(1)=boundp(x(ii),fmin1,fmax1,pen,1000.);
      //cout << "vb_coff(1)" << ii << endl;
      ii=ii+1;
  
     ///*
      //if (parest_flags(250))
      if (pmsd)
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          pmsd->vb_coff(i,1)=boundp(x(ii),pmsd->fmin1(i),
            pmsd->fmax1(i),pen,1000.);
          ii=ii+1;
        }
      }
     // */
     
    }
  
    if (parest_flags(13) > 0)
    {
      MY_DOUBLE_TYPE nscal=1.0;
      if (parest_flags(387))  //NMD_4sep2020
      {
        nscal=3000.0;
      }
      else
      {
        nscal=1000.0;
      }
  //    vb_coff(2)=boundp(x(ii),fminl,fmaxl,pen,3000.);
      vb_coff(2)=boundp(x(ii),fminl,fmaxl,pen,nscal);
  //NMD_4sep2020
      //cout << "vb_coff(2)" << ii << endl;
      ii=ii+1;
     
      ///*
      //if (parest_flags(251))
      if (pmsd)
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          pmsd->vb_coff(i,2)=boundp(x(ii),pmsd->fminl(i),
            pmsd->fmaxl(i),pen,nscal);  //NMD_4sep2020
  //          pmsd->fmaxl(i),pen,3000.);
          ii=ii+1;
        }
      }
      //*/
    
    }
    if (parest_flags(14) > 0)
    {
      MY_DOUBLE_TYPE nscal=1.0;
      if (parest_flags(387))  //NMD_4sep2020
      {
        nscal=3000.0;
      }
      else
      {
        nscal=1000.0;
      }
      vb_coff(3)=boundp(x(ii),rhomin,rhomax,pen,nscal);    //NMD_4sep2020
  //    vb_coff(3)=boundp(x(ii),rhomin,rhomax,pen,3000.);
      //cout << "vb_coff(3)" << ii << endl;
      ii=ii+1;
      if (pmsd)
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          pmsd->vb_coff(i,3)=boundp(x(ii),pmsd->rhomin(i),
            pmsd->rhomax(i),pen,nscal);      //NMD_4sep2020
  //          pmsd->rhomax(i),pen,3000.);
          ii=ii+1;
        }
      }
    }
  
    if (parest_flags(227) > 0)
    {
      vb_coff(4)=boundp(x(ii),-5.0,5.0,pen,100.);
      //cout << "vb_coff(3)" << ii << endl;
      ii=ii+1;
      if (pmsd)  //NMD_1Apr2020
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          pmsd->vb_coff(i,4)=boundp(x(ii),-5.0,5.0,pen,100.);
          ii=ii+1;
        }
      }  //NMD_1Apr2020
    }
  
    if (parest_flags(21) > 0)
    {
      sv(1)=boundp(x(ii),-1.e0,1.e0,pen);
      ii=ii+1;
      sv(2)=boundp(x(ii),-2.5e0,2.5e0,pen);
      ii=ii+1;
    }
  
    if (age_flags(55) > 0)
    {
      sv(7)=boundp(x(ii),.0001,5.,pen);
      ii=ii+1;
      sv(8)=boundp(x(ii),.01,100.,pen);
      ii=ii+1;
    }
    if (age_flags(56) > 0)
    {
      sv(9)=boundp(x(ii),1.0001,10.,pen);
      ii=ii+1;
    }
  
    if (parest_flags(15) > 0)
    {
      var_coff(1)=boundp(x(ii),vmin1,vmax1,pen,1000.);
      ii=ii+1;
      if (pmsd)
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          pmsd->var_coff(i,1)=boundp(x(ii),pmsd->vmin1(i),pmsd->vmax1(i),
            pen,1000.);
          ii=ii+1;
        }
      }
    }
  
    if (parest_flags(16) > 0)
    {
      var_coff(2)=boundp(x(ii),vminl,vmaxl,pen,1000.);
      ii=ii+1;
      if (pmsd)
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          pmsd->var_coff(i,2)=boundp(x(ii),pmsd->vminl(i),pmsd->vmaxl(i),
            pen,1000.);
          ii=ii+1;
        }
      }
    }
    //  natural mortality deviation for terminal age classes
    if (age_flags(47)>0)
    {
      cout << "In sv(1) reset" << endl;
      sv(1)=boundp(x(ii++),-1.,1.,pen);
    }
  }
  // extra par for debugging
  if (parest_flags(381)>0)
  {
    debug_par=boundp(x(ii++),-1.,1.,pen);
  }
  return pen;
}


int dvar_len_fish_stock_history::nvcal()
{
  int ii=0;

  if (parest_flags(392)==0)
  {
    if (age_flags(30)>0)
    {
      if (parest_flags(155)!=0)
      {
        ii+=recinpop_orth_size_count();
      }
    }
  }
  
  if (parest_flags(361)>0)
  {
    // testpar
    ii+=1;
  }
  ii+=dvar_fish_stock_history::nvcal(len_sample_size);
  if (parest_flags(392)==0)
  {
//    ii+=dvar_fish_stock_history::nvcal(len_sample_size);
  
    if (!pmsd || sum(column(pmsd->species_flags,2))==0)
    {
      if (age_flags(146)>0  && age_flags(94) !=3)
      {
        ii++;
        if (pmsd)
        {
          ii+=pmsd->num_species-1;
        }
      }
      if (age_flags(162)>0)
      {
        ii++;
        if (pmsd)
        {
          ii+=pmsd->num_species-1;
        }
      }
    }
    else
    {
      ivector sf2=column(pmsd->species_flags,2);
      if (age_flags(146)>0  && age_flags(94) !=3)
      {
        if (sf2(1))
          ii++;
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          if (sf2(i))
          ii++;
        }
      }
      if (age_flags(162)>0)
      {
        if (sf2(1))
          ii++;
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          if (sf2(i))
          ii++;
        }
      }
    }
      
    if (age_flags(151)>0)
      ii++;
    if (age_flags(152)>0)
      ii++;
    if (age_flags(152)>0)
      ii++;
    // length bias parameters
    ii+=num_active(vb_bias,column(fish_flags,11),column(fish_flags,22));
  
    if (age_flags(58)>0)
    {
      ii++;
    }
    if (parest_flags(174))
    {
      ii++;
      ii++;
      ii++;
    }
    if (parest_flags(168)>0)
    {
      ii++;
    }
  
    if (parest_flags(169)>0)
    {
      ii++;
    }
  
    if (parest_flags(170)>0)
    {
      ii++;
    }
  
    if (parest_flags(171)>0)
    {
      ii++;
    }
  
    if (parest_flags(172)>0)
    {
      ii++;
    }
  
    if (age_flags(59)>0)
    {
      ii++;
      ii++;
    }
    if (age_flags(63)>0)
    {
      ii++;
    }
    if (age_flags(65)>0)
    {
  	 ii++;
    }
    if (age_flags(66)>0)
    {
      ii++;
      if (pmsd)
      {
        int numsp=pmsd->num_species;
        ii+=numsp-1;
      }
    }
    // von Betalanffy parameters
    if (parest_flags(157)>0)
    {
  	 ii=ii+1;
    }
  
    if (parest_flags(158)>0)
    {
      ii=ii+1;
    }
  
    if (parest_flags(159)>0)
    {
      ii=ii+1;
    }
  
    if (age_flags(102)>0)
    {
      ii++;
    }
    if (age_flags(153)>0)
    {
      ii++;
      if (pmsd)                                      //NMD 16Aug2011
      {
        ii+=pmsd->num_species-1;
      }                                      //NMD 16Aug2011
    }
  
    if (age_flags(67)>0)
    {
      ii+=size_count(age_pars(1));
    }
  
    if (parest_flags(175)>1)
    {
      ii+=size_count(age_pars(4));
      ii++;
    }
    if (age_flags(73)>0)
    {
      if (age_flags(81)<1)
      {
        ii+=size_count(age_pars(2));
      }
      else
      {
        int num=nage-age_flags(81)+1;
        ii+=size_count_partial(age_pars(2),num);
      }
    }
    if (parest_flags(184)>0)
    {
      int num=parest_flags(173)-1;
      ii+=size_count_partial(age_pars(3),num);
      if (pmsd)
      {
        int numsp=pmsd->num_species;
        for (int is=2;is<=numsp;is++)
        {
          ii+=size_count_partial(get_age_pars_species(is,3),num);
        }
      }
    }
    if (parest_flags(121)>0)
    {
      int num=parest_flags(121);
      ii+=size_count_partial(age_pars(5),num);
      if (pmsd)
      {
        int numsp=pmsd->num_species;
        for (int is=2;is<=numsp;is++)
        {
          ii+=size_count_partial(get_age_pars_species(is,5),num);
        }
      }
    }
  
  
    if (parest_flags(160)>0)
    {
      ii=ii+1;
    }
    if (parest_flags(163)>0)
    {
      ii+=size_count(growth_dev);
    }
  
    if (parest_flags(12)>0)
    {
      ii=ii+1;
     
     ///*
      //if (parest_flags(250))
      if (pmsd)
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          ii=ii+1;
        }
      }
      //*/
      
    }
  
    if (parest_flags(13) > 0)
    {
      ii=ii+1;
      
      ///*
      //if (parest_flags(251))
      if (pmsd)
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          ii=ii+1;
        }
      }
      //*/
     
    }
    if (parest_flags(14) > 0)
    {
      ii=ii+1;
      if (pmsd)
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          ii=ii+1;
        }
      }
    }
    if (parest_flags(227) > 0)
    {
      ii=ii+1;
      if (pmsd)  //NMD_1Apr2020
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          ii=ii+1;
        }
      }  //NMD_1Apr2020
    }
  
    if (parest_flags(21) > 0)
    {
      ii=ii+1;
      ii=ii+1;
    }
  
    if (age_flags(55) > 0)
    {
      ii=ii+1;
      ii=ii+1;
    }
  
    if (age_flags(56) > 0)
    {
      ii=ii+1;
    }
  
    if (parest_flags(15) > 0)
    {
      ii=ii+1;
      if (pmsd)
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          ii=ii+1;
        }
      }
    }
  
    if (parest_flags(16) > 0)
    {
      ii=ii+1;
      if (pmsd)
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          ii=ii+1;
        }
      }
    }
    //  natural mortality deviation for terminal age classes
    if (age_flags(47)>0)
    {
     // cout << "In sv(1) nvcal" << endl;
      ii=ii+1;
    }
  }
  if (parest_flags(381)>0)
  {
    ii=ii+1;
  }
  if (ii==0 && parest_flags(197)) ii=1;
  return ii;
}


void dvar_len_fish_stock_history::xinit(dvector& x)
{

  {
    ofstream xof("xinit.rpt");
    //cout << "PPP  " << setprecision(10) << norm(value(diff_coffs)) << endl;
    int af184=historical_age_flags(184);
    rationalize_movement_coffs(af184);
    int ii=1;
    int iisave;

    if (parest_flags(392)==0)
    {
      if (age_flags(30)>0)
      {
        if (parest_flags(155)!=0)
        {
          recinpop_orth_xinit(xof,x,ii);
        }
      }
    }

    if (parest_flags(361)>0)
    {
      iisave=ii;
      set_value_inv(testpar,x,ii);
      xinit_message(xof,iisave,ii,"testpar");
    }
    dvar_fish_stock_history::xinit(x,ii,len_sample_size,xof);

    if (parest_flags(392)==0)
    {
//      dvar_fish_stock_history::xinit(x,ii,len_sample_size,xof);
      // Length bias parameters
      if (!pmsd || sum(column(pmsd->species_flags,2))==0)
      {
        if (age_flags(146)>0  && age_flags(94) !=3)
        {
          iisave=ii;
          set_value_inv(sv(21),x,ii,1.e-6,1000.0,1000.);
          xinit_message(xof,iisave,ii,"sv(21)");
          if (pmsd)
          {
            int numsp=pmsd->num_species;
            for (int i=2;i<=numsp;i++)
            {
              iisave=ii;
              set_value_inv(pmsd->sv(i,21),x,ii,1.e-6,1000.00,1000.);
              xinit_message(xof,iisave,ii,"sv(21)");
            }
          }
        }
        if (age_flags(162)>0)
        {
          iisave=ii;
          set_value_inv(sv(29),x,ii,0.201,0.999);
          xinit_message(xof,iisave,ii,"sv(29)");
          if (pmsd)
          {
            int numsp=pmsd->num_species;
            for (int i=2;i<=numsp;i++)
            {
              iisave=ii;
              set_value_inv(pmsd->sv(i,29),x,ii,0.201,0.999);
              xinit_message(xof,iisave,ii,"sv(29)");
            }
          }
        }
      }
      else
      {
        ivector sf2=column(pmsd->species_flags,2);
        if (age_flags(146)>0  && age_flags(94) !=3)
        {
          if (sf2(1))
          {
            iisave=ii;
            set_value_inv(sv(21),x,ii,1.e-6,1000.0,1000.);
            xinit_message(xof,iisave,ii,"sv(21)");
          }
          int numsp=pmsd->num_species;
          for (int i=2;i<=numsp;i++)
          {
            if (sf2(i))
            {
              iisave=ii;
              set_value_inv(pmsd->sv(i,21),x,ii,1.e-6,1000.00,1000.);
              xinit_message(xof,iisave,ii,"sv(21)");
            }
          }
        }
        if (age_flags(162)>0)
        {
          if (sf2(1))
          {
            iisave=ii;
            set_value_inv(sv(29),x,ii,0.201,0.999);
            xinit_message(xof,iisave,ii,"sv(29)");
          }
          int numsp=pmsd->num_species;
          for (int i=2;i<=numsp;i++)
          {
            if (sf2(i))
            {
              iisave=ii;
              set_value_inv(pmsd->sv(i,29),x,ii,0.201,0.999);
              xinit_message(xof,iisave,ii,"sv(29)");
            }
          }
        }
      }
    
    
      if (age_flags(151)>0)
      {
        iisave=ii;
        set_value_inv(sv(22),x,ii,-.19,10.0);
        xinit_message(xof,iisave,ii,"sv(22)");
      }
    
      if (age_flags(152)>0)
      {
        iisave=ii;
        set_value_inv(sv(23),x,ii,-0.8,5.0);
        xinit_message(xof,iisave,ii,"sv(23)");
      }
    
      if (age_flags(152)>0)
      {
        iisave=ii;
        set_value_inv(sv(24),x,ii,-0.9,10.0);
        xinit_message(xof,iisave,ii,"sv(24)");
      }
    
      iisave=ii;
      set_value_inv(vb_bias,x,ii,0.0,fmax1-fmin1,column(fish_flags,11),
        column(fish_flags,22));
      xinit_message(xof,iisave,ii,"sv(13)");
    
      if (age_flags(58)>0)
      {
        iisave=ii;
        MY_DOUBLE_TYPE bnd=age_flags(58)/10.;
        x(ii++)=boundpin(sv(10),0.0,bnd,10.);
        xinit_message(xof,iisave,ii,"sv(10)");
      }
      if (parest_flags(174))
      {
        iisave=ii;
        x(ii++)=boundpin(sv(17),0.0,20.0,10.);
        xinit_message(xof,iisave,ii,"sv(17)");
        iisave=ii;
        x(ii++)=boundpin(sv(18),0.0,5.0,10.);
        xinit_message(xof,iisave,ii,"sv(18)");
        iisave=ii;
    #if !defined(NO_MY_DOUBLE_TYPE)
        x(ii++)=boundpin(sv(19),-1.0L,1.0L,10.);
    #else
        x(ii++)=boundpin(sv(19),-1.0,1.0,10.);
    #endif
        xinit_message(xof,iisave,ii,"sv(19)");
      }
    
      if (parest_flags(168)>0)
      {
        iisave=ii;
        x(ii++)=boundpin(sv(16),0.2,5.0,10.);
        //x(ii++)=boundpin(sv(16),-.4,.4,10.);
        xinit_message(xof,iisave,ii,"sv(16)");
      }
     
      if (parest_flags(169)>0)
      {
        iisave=ii;
        x(ii++)=boundpin(sv(17),-1.9,1.9,10.);
        //x(ii++)=boundpin(sv(17),-.4,.4,10.);
        xinit_message(xof,iisave,ii,"sv(17)");
      }
     
      if (parest_flags(170)>0)
      {
        iisave=ii;
        x(ii++)=boundpin(sv(18),-6.9,6.9,10.);
        //x(ii++)=boundpin(sv(17),-.4,.4,10.);
        xinit_message(xof,iisave,ii,"sv(18)");
      }
     
      if (parest_flags(171)>0)
      {
        iisave=ii;
        x(ii++)=boundpin(sv(19),-6.9,6.9,10.);
        //x(ii++)=boundpin(sv(17),-.4,.4,10.);
        xinit_message(xof,iisave,ii,"sv(19)");
      }
     
      if (parest_flags(172)>0)
      {
        iisave=ii;
        x(ii++)=boundpin(sv(20),-6.9,6.9,10.);
        //x(ii++)=boundpin(sv(17),-.4,.4,10.);
        xinit_message(xof,iisave,ii,"sv(19)");
      }
     
    
      if (age_flags(59)>0)
      {
        iisave=ii;
        x(ii++)=boundpin(sv(11),-6.,6.,10.);
        xinit_message(xof,iisave,ii,"sv(11)");
        iisave=ii;
        x(ii++)=boundpin(sv(12),-6.,6.,10.);
        xinit_message(xof,iisave,ii,"sv(12)");
      }
    
      if (age_flags(63)>0)
      {
        iisave=ii;
        x(ii++)=boundpin(sv(13),-3.,3.,1000.);
        xinit_message(xof,iisave,ii,"sv(13)");
      }
    
      if (age_flags(65)>0)
      {
        iisave=ii;
        x(ii++)=boundpin(sv(14),-5.,5.,1000.);
        xinit_message(xof,iisave,ii,"sv(14)");
      }
      if (age_flags(66)>0)
      {
        iisave=ii;
        x(ii++)=boundpin(sv(15),-2.,2.,100000.);
        if (pmsd)
        {
          int numsp=pmsd->num_species;
          for (int i=2;i<=numsp;i++)
          {
            x(ii++)=boundpin(pmsd->sv(i,15),-2.,2.,100000.);
          }
        }
      }
    
      if (parest_flags(157)>0)
      {
    	 iisave=ii;
    	 x(ii++)=boundpin(sv(3),-5.,5.,100.);
    	 xinit_message(xof,iisave,ii,"sv(3)");
      }
    
      if (parest_flags(158)>0)
      {
    	 iisave=ii;
    	 x(ii++)=boundpin(sv(4),-5.,5.,100.);
    	 xinit_message(xof,iisave,ii,"sv(4)");
      }
    
      if (parest_flags(159)>0)
      {
        iisave=ii;
        x(ii++)=boundpin(sv(5),-5.,5.,100.);
        xinit_message(xof,iisave,ii,"sv(5)");
      }
    
    
      if (age_flags(102)>0)
      {
        iisave=ii;
        set_value_inv(sv(25),x,ii,-10.,10.);
        xinit_message(xof,iisave,ii,"sv(25)");
      }
      if (age_flags(153)>0)
      {
        iisave=ii;                             // added by JH & PK 21-7-04
        set_value_inv(sv(26),x,ii,-2.0,2.0);
        xinit_message(xof,iisave,ii,"sv(26)"); // added by JH & PK 21-7-04
        if (pmsd)                                      //NMD 16Aug2011
        {
          int numsp=pmsd->num_species;
          for (int i=2;i<=numsp;i++)
          {
            iisave=ii;
            set_value_inv(pmsd->sv(i,26),x,ii,-2.0,2.0);
            xinit_message(xof,iisave,ii,"sv(26)");
          }
        }                                      //NMD 16Aug2011
      }
      if (age_flags(67)>0)
      {
        cout <<"newmaux5.cpp " << age_pars(1) << endl;
        iisave=ii;
    #if !defined(NO_MY_DOUBLE_TYPE)
        set_value_inv(age_pars(1),x,ii,.0,1.0L);
    #else
        set_value_inv(age_pars(1),x,ii,.0,1.0);
    #endif
        xinit_message(xof,iisave,ii,"age_pars(1)");
      }
      if (parest_flags(175)>1)
      {
        iisave=ii;
        MY_DOUBLE_TYPE fivet=5000.;
        set_value_inv(age_pars(4),x,ii,0.001,100.,fivet);
        xinit_message(xof,iisave,ii,"age_pars(4)");
        iisave=ii;
        set_value_inv(sv(20),x,ii,0.001,100.0,1000.0);
        xinit_message(xof,iisave,ii,"sv(20)");
      }
    
      if (age_flags(73)>0)
      {
        iisave=ii;
        if (age_flags(81)<1)
        {
          set_value_inv(age_pars(2),x,ii,-5.0,5.0,5000.0);
        }
        else
        {
          int num=nage-age_flags(81)+1;
          set_value_inv_partial(age_pars(2),x,ii,num,-5.0,5.0,5000.0);
        }
        xinit_message(xof,iisave,ii,"age_pars(2)");
      }
    
      if (parest_flags(184)>0)
      {
        iisave=ii;
        int num=parest_flags(173)-1;
        set_value_inv_partial(age_pars(3),x,ii,num,-50.0,50.0,1000.0);
        xinit_message(xof,iisave,ii,"xage_pars(3)");
        if (pmsd)
        {
          int numsp=pmsd->num_species;
          for (int is=2;is<=numsp;is++)
          {
            iisave=ii;
            set_value_inv_partial(get_age_pars_species(is,3),x,ii,num,
              -50.0,50.0,1000.0);
            xinit_message(xof,iisave,ii,"age_pars_species(3)");
          }
        }
      }
      if (parest_flags(121)>0)
      {
        iisave=ii;
        int num=parest_flags(121);
        MY_DOUBLE_TYPE nscal=1.0;
        if (parest_flags(387))  //NMD_4sep2020
        {
          nscal=4500.0;
        }
        else
        {
//          nscal=100.0;
          nscal=200.0;
        }
        set_value_inv_partial(age_pars(5),x,ii,num,-20.0,2.0,nscal);
    //    set_value_inv_partial(age_pars(5),x,ii,num,-20.0,2.0,4500.0);
    //NMD_4sep2020
        xinit_message(xof,iisave,ii,"age_pars(5)");
        if (pmsd)
        {
          int numsp=pmsd->num_species;
          for (int is=2;is<=numsp;is++)
          {
            iisave=ii;
            set_value_inv_partial(get_age_pars_species(is,5),x,ii,num,
              -20.0,2.0,nscal);    //NMD_4sep2020
    //          -20.0,2.0,4500.0);
            xinit_message(xof,iisave,ii,"age_pars(5)");
          }
        }
      }
    
      if (parest_flags(160)>0)
      {
        iisave=ii;
        x(ii)=boundpin(sv(6),-0.1,5.,100.);
        //cout << "vb_coff(1)" << ii << endl;
        ii=ii+1;
        xinit_message(xof,iisave,ii,"sv(5)");
      }
    
      if (parest_flags(163)>0)
      {
        set_value_inv(growth_dev,x,ii,-5.0,5.0);
      }
    
      if (parest_flags(12)>0)
      {
        iisave=ii;
        x(ii)=boundpin(vb_coff(1),fmin1,fmax1,1000.);
        ii=ii+1;
        xinit_message(xof,iisave,ii,"vb_coff(1)");
        ///*
        //if (parest_flags(250))
        if (pmsd)
        {
          int numsp=pmsd->num_species;
          for (int i=2;i<=numsp;i++)
          {
            iisave=ii;
            x(ii)=boundpin(pmsd->vb_coff(i,1),pmsd->fmin1(i),
              pmsd->fmax1(i),1000.);
            ii=ii+1;
            xinit_message(xof,iisave,ii,"vb_coff(1)");
          }
        }
      // */
      }
    
      if (parest_flags(13) > 0)
      {
        iisave=ii;
        MY_DOUBLE_TYPE nscal=1.0;
        if (parest_flags(387))  //NMD_4sep2020
        {
          nscal=3000.0;
        }
        else
        {
          nscal=1000.0;
        }
        x(ii)=boundpin(vb_coff(2),fminl,fmaxl,nscal);
    //    x(ii)=boundpin(vb_coff(2),fminl,fmaxl,3000.);
    //NMD_4sep2020
        ii=ii+1;
        xinit_message(xof,iisave,ii,"vb_coff(2)");
       // /*
        //if (parest_flags(251))
        if (pmsd)
        {
          int numsp=pmsd->num_species;
          for (int i=2;i<=numsp;i++)
          {
            iisave=ii;
            x(ii)=boundpin(pmsd->vb_coff(i,2),pmsd->fminl(i),
              pmsd->fmaxl(i),nscal);
    //          pmsd->fmaxl(i),3000.);    //NMD_4sep2020
            xinit_message(xof,iisave,ii,"vb_coff(2)");
            ii=ii+1;
          }
        }
       // */
      }
      if (parest_flags(14) > 0)
      {
        iisave=ii;
        MY_DOUBLE_TYPE nscal=1.0;
        if (parest_flags(387))  //NMD_4sep2020
        {
          nscal=3000.0;
        }
        else
        {
          nscal=1000.0;
        }
        x(ii)=boundpin(vb_coff(3),rhomin,rhomax,nscal);
    //    x(ii)=boundpin(vb_coff(3),rhomin,rhomax,3000.);
    //NMD_4sep2020
        ii=ii+1;
        xinit_message(xof,iisave,ii,"vb_coff(3)");
        if (pmsd)
        {
          int numsp=pmsd->num_species;
          for (int i=2;i<=numsp;i++)
          {
            iisave=ii;
            x(ii)=boundpin(pmsd->vb_coff(i,3),pmsd->rhomin(i),
              pmsd->rhomax(i),nscal);
    //          pmsd->rhomax(i),3000.);     //NMD_4sep2020
            xinit_message(xof,iisave,ii,"vb_coff(3)");
            ii=ii+1;
          }
        }
      }
    
      if (parest_flags(227) > 0)
      {
        iisave=ii;
        x(ii)=boundpin(vb_coff(4),-5.0,5.0,100.);
        ii=ii+1;
        xinit_message(xof,iisave,ii,"vb_coff(4)");
        if (pmsd)  //NMD_1Apr2020
        {
          int numsp=pmsd->num_species;
          for (int i=2;i<=numsp;i++)
          {
            iisave=ii;
            x(ii)=boundpin(pmsd->vb_coff(i,4),-5.0,5.0,100.);
            xinit_message(xof,iisave,ii,"vb_coff(4)");
            ii=ii+1;
          }
        }  //NMD_1Apr2020
      }
    
      if (parest_flags(21) > 0)
      {
        iisave=ii;
        x(ii)=boundpin(sv(1),-1.e0,1.e0);
        ii=ii+1;
        x(ii)=boundpin(sv(2),-2.5e0,2.5e0);
        ii=ii+1;
        xinit_message(xof,iisave,ii,"sv");
      }
    
      if (age_flags(55) > 0)
      {
        iisave=ii;
        x(ii)=boundpin(sv(7),.0001,5.);
        ii=ii+1;
        x(ii)=boundpin(sv(8),.01,100.);
        ii=ii+1;
        xinit_message(xof,iisave,ii,"sv");
      }
    
      if (age_flags(56) > 0)
      {
        iisave=ii;
        x(ii)=boundpin(sv(9),1.0001,10.);
        ii=ii+1;
        xinit_message(xof,iisave,ii,"sv");
      }
    
      if (parest_flags(15) > 0)
      {
        iisave=ii;
        x(ii)=boundpin(var_coff(1),vmin1,vmax1,1000.);
        ii=ii+1;
        xinit_message(xof,iisave,ii,"var_coff(1)");
        if (pmsd)
        {
          int numsp=pmsd->num_species;
          for (int i=2;i<=numsp;i++)
          {
            iisave=ii;
            x(ii)=boundpin(pmsd->var_coff(i,1),pmsd->vmin1(i),
              pmsd->vmax1(i),1000.);
            ii=ii+1;
            xinit_message(xof,iisave,ii,"pmsd->var_coff(i,1)");
          }
        }
      }
    
      if (parest_flags(16) > 0)
      {
        iisave=ii;
        x(ii)=boundpin(var_coff(2),vminl,vmaxl,1000.);
        ii=ii+1;
        xinit_message(xof,iisave,ii,"var_coff(2)");
        if (pmsd)
        {
          int numsp=pmsd->num_species;
          for (int i=2;i<=numsp;i++)
          {
            iisave=ii;
            x(ii)=boundpin(pmsd->var_coff(i,2),pmsd->vminl(i),
              pmsd->vmaxl(i),1000.);
            ii=ii+1;
            xinit_message(xof,iisave,ii,"pmsd->var_coff(i,2)");
          }
        }
      }
      //  natural mortality deviation for terminal age classes
      if (age_flags(47)>0)
      {
        iisave=ii;
        cout << "In sv(1) xinit" << endl;
        x(ii++)=boundpin(sv(1),-1.,1.);
        xinit_message(xof,iisave,ii,"sv(1)");
      }
    }
    // extra par for debugging
    if (parest_flags(381)>0)
    {
      iisave=ii;
      // always set it to zero.
      x(ii++)=0.0;
      //x(ii++)=boundpin(debug_par,-1.,1.);
      xinit_message(xof,iisave,ii,"debug_par");
    }
    if (parest_flags(385)==1)
    {
      x=elem_prod(x,sqrt(pos_hess_scale));
    }
    else if (parest_flags(385)==2)
    {
      x=psqrt_hess*x;
    }
  }
  if (parest_flags(390)  && parest_flags(391))
  {
    ifstream ifs("x_vector." + str(parest_flags(390)));
    if (!ifs)
    {
      cerr << "Error trying to open restart file "
        << "x_vector." + str(parest_flags(390)) << endl;
      ad_exit(1);
    }
    int nvar=x.indexmax();;
    int nvar1=0;
    int ifn=0;
    int itn=0;
    dvector xx(1,nvar);
    xx.initialize();
    do
    { 
      ifs >> itn >> ifn >>nvar1;
      if (nvar != nvar1)
      {
        cerr << "Dimension Error  in restart file " <<
           " dimension is " << nvar1 << "For current model is " << nvar << endl;
        ad_exit(1);
        if (!ifs)
        {
          cerr << "Error trying read headers from restart file "
            << "x_vector." + str(parest_flags(390)) << endl;
          ad_exit(1);
        }
      }
      ifs >> xx;
      if (!ifs)
      {
        cerr << "Error trying x_vector from restart file "
          << "x_vector." + str(parest_flags(390)) << endl;
        ad_exit(1);
      }
    }
    while (ifn != parest_flags(391));
    x=xx;
    return;
  }
}



int get_max_bias(_CONST ivector& common_vb_bias,int num_fisheries)
{
  int mmax=common_vb_bias(1);
  for (int i=2;i<=num_fisheries;i++)
  { 
    if (common_vb_bias(i)>mmax)
    {
      mmax=common_vb_bias(i);
    }
  } 
  return mmax;
}
#undef HOME_VERSION
