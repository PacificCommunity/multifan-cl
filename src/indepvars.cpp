/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
const int REGION_REC_SCALE=100.;


void dvar_fish_stock_history::get_indep_vars_for_report(dvar_vector& x, int& ii,
							int& runsetv)
{
  dvariable pen=0.;
  int jj=0;
  jj=ii;
  if (age_flags(68)>0)
  {
    switch(age_flags(184))
    {
    case 0:
      {
        if (runsetv==1)
        {
          MY_DOUBLE_TYPE nscal=1.0;
          nscal=1000.0;
          set_value(diff_coffs,x,ii,0.0,3.0,pen,nscal,age_flags(68),move_map);
        }
        dmatrix idv;
        idv=value(diff_coffs);
        put_in_indep_vars(idv,0.0,3.0,indep_var,indep_var_lo,
                          indep_var_hi,jj);
      }  
      break;
    case 1:
      {
        if (runsetv==1)
        {
          set_value(xdiff_coffs,x,ii,-200.0,200.0,pen,1.e+3,age_flags(68),move_map);
        }
        indepvars_error();
      }  
      break;
    case 2:
      {
        if (runsetv==1)
        {
          set_value(zdiff_coffs,x,ii,-10.0,10.0,pen,1.e+2,age_flags(68),move_map);
        }
        indepvars_error();
      }  
      break;
    default:
      cerr << "Illegal value for age_flags(184)" << endl;
      ad_exit(1);
    }
  }

// - region_rec_diff_coffs
//  jj=ii;
  if (age_flags(70)>0)
  {
    if (parest_flags(155)==0)
    {
      if (num_regions>1)
      {
        if (runsetv==1)
        {
          dvar_matrix tmp(1,num_regions,first_unfixed_year,last_real_year-1);
          set_value(tmp,x,ii,-3.,3.,pen,REGION_REC_SCALE,1,region_flags(3));
        }
        dmatrix idv;
        idv=
          value(trans(region_rec_diff_coffs.sub(first_unfixed_year,last_real_year-1)));
        put_in_indep_vars(idv,-3.0,3.0,indep_var,indep_var_lo,
                        indep_var_hi,jj);
      }
    }
  }

  if (sum_ff48_flag)
  { 
    if (runsetv==1)
    {
      ivector ff48=column(fish_flags,48);
      ivector ffgroup=column(fish_flags,24);
      MY_DOUBLE_TYPE scale=1000.0;
      set_value(bs_selcoff,x,ii,-20.,7.,pen,scale ,ff48,ffgroup);
    }
    d4_array idv;
    imatrix fshgptr;
    fshgptr=fishery_group_ptr(24);
    int mmin=fshgptr.indexmin();
    int mmax=fshgptr.indexmax();
    idv.allocate(mmin,mmax);
    for (int i=mmin;i<=mmax;i++)
    {
      idv(i)=value(bs_selcoff(fshgptr(i,1)));
    }
    put_in_indep_vars(idv,-20.0,7.0,indep_var,indep_var_lo,
                        indep_var_hi,jj);
  }

  // relative ecruitment
  ivector sf=season_flags(1);
  if (age_flags(30)>0)
  {
    if (parest_flags(155)==0)
    {
      if (age_flags(51)==0)
      {
        if (sum(sf))
        {
          if (parest_flags(400)==0)  //NMD_19May2016
          {
            if (runsetv==1)
            {
              set_value(recr(max(2,first_unfixed_year),last_real_year),
                 x,ii,-20.,20,pen,100,year_flags(1),-10.0);
            }
            indepvars_error();
          }
          else
          {
            if (runsetv==1)
            {
              recr(last_real_year-parest_flags(400)+1,last_real_year).
                initialize();
              set_value(recr(max(2,first_unfixed_year),
                last_real_year-parest_flags(400)),
                x,ii,-20.,20,pen,100,year_flags(1),-10.0);
            }
            indepvars_error();
          }  //NMD_19May2016
          if (pmsd && pmsd->num_species>1)
          {
            int ns=pmsd->num_species;
            for (int is=2;is<=ns;is++)
            {
              if (parest_flags(400)==0)   //NMD_19May2016
              {
                if (runsetv==1)
                {
                  set_value(pmsd->recr(is)(max(2,first_unfixed_year),
                    last_real_year),x,ii,-20.,20,pen,
                    100,year_flags(1),-10.0);
                }
                indepvars_error();
              }
              else
              {
                pmsd->recr(is)(last_real_year-parest_flags(400)+1,last_real_year).initialize();
                if (runsetv==1)
                {
                  set_value(pmsd->recr(is)(max(2,first_unfixed_year),
                    last_real_year-parest_flags(400)),
                    x,ii,-20.,20,pen,100,year_flags(1),-10.0);
                }
                indepvars_error();
              }   //NMD_19May2016
            }
          }
        }
        else
        {
          if (pmsd && pmsd->num_species>1)
          {
            int ns=pmsd->num_species;
            int srf3=sum(column(pmsd->species_flags,1));
            if (!srf3)
            {
//              set_value(recr(1,last_real_year),x,ii,-20.,20.,1000.);
              if (parest_flags(400)==0)  //NMD_16May2016
              {
                if (runsetv==1)
                {
                  set_value(recr(max(2,first_unfixed_year),
                    last_real_year),x,ii,-20.,20.,1000.);
                }
                indepvars_error();
              }
              else
              {
                if (runsetv==1)
                {
                  recr(last_real_year-parest_flags(400)+1,last_real_year).
                    initialize();
                  set_value(recr(max(2,first_unfixed_year),
                    last_real_year-parest_flags(400)),x,ii,-20.,20.,1000.);
                }
                indepvars_error();
              }  //NMD_16May2016
              for (int is=2;is<=ns;is++)
              {
                if (parest_flags(400)==0)   //NMD_19May2016
                {
                  if (runsetv==1)
                  {
                    set_value(pmsd->recr(is)(max(2,first_unfixed_year),
                      last_real_year),x,ii,-20.,20.,1000.);
                  }
                  indepvars_error();
                }
                else
                {
                  if (runsetv==1)
                  {
                    pmsd->recr(is)(last_real_year-parest_flags(400)+1,last_real_year).
                      initialize();
                    set_value(pmsd->recr(is)(max(2,first_unfixed_year),
                      last_real_year-parest_flags(400)),x,ii,-20.,20.,1000.);
                  }
                  indepvars_error();
                }   //NMD_19May2016
              }
            }
            else
            {
              dvar_matrix M(1,ns);
              M(1)=recr(max(2,first_unfixed_year),last_real_year);
              if (parest_flags(400)!=0)
                M(1)(last_real_year-parest_flags(400)+1,last_real_year).
                initialize();
              const MY_REAL_DOUBLE th=1000.;
              for (int is=2;is<=ns;is++)
              {
                M(is)=pmsd->recr(is)(max(2,first_unfixed_year),last_real_year);
                if (parest_flags(400)!=0)
                  M(is)(last_real_year-parest_flags(400)+1,last_real_year).
                  initialize();
              }
              if (parest_flags(400)==0)  //NMD_19May2016
              {
                if (runsetv==1)
                {
                  set_value(M,x,ii,-20.,20.,th, 1,column(pmsd->species_flags,1));
                }
                indepvars_error();
              }
              else
              {
                dvar_matrix MM(1,ns);
                for (int is=1;is<=ns;is++)
                {
                  MM(is)=M(is)(2,last_real_year-parest_flags(400));
                }
                if (runsetv==1)
                {
                  set_value(MM,x,ii,-20.,20.,th, 1,column(pmsd->species_flags,1));
                }
                indepvars_error();
              }
            }
          }
          else
          {
            if (parest_flags(400)==0)
            {
              if (runsetv==1)
              {
                set_value(recr(max(2,first_unfixed_year),last_real_year),
                  x,ii,-20.,20.,1000.);
              }
              dvector idv;
              idv=value(recr.sub(max(2,first_unfixed_year),
                    last_real_year).shift(1));
              put_in_indep_vars(idv,-20.0,20.0,indep_var,indep_var_lo,
                indep_var_hi,jj);
            }
            else
            {
              recr(last_real_year-parest_flags(400)+1,last_real_year).
                initialize();
              if (runsetv==1)
              {
                set_value(recr(max(2,first_unfixed_year),
                  last_real_year-parest_flags(400)),
                  x,ii,-20.,20.,1000.);
              }
              dvector idv;
              idv=value(recr.sub(max(2,first_unfixed_year),
                    last_real_year-parest_flags(400)).shift(1));
              put_in_indep_vars(idv,-20.0,20.0,indep_var,indep_var_lo,
                indep_var_hi,jj);
            }
          }
        }
      }
      else
      {
        MY_DOUBLE_TYPE bd=age_flags(51);
        MY_DOUBLE_TYPE bd1=age_flags(46)/10.0;
        if (parest_flags(400)==0)
        {
          mfcl_error();
          if (runsetv==1)
          {
            set_value(recr,x,ii,-bd,bd1,pen);
          }
          indepvars_error();
        }
        else
        {
          if (runsetv==1)
          {
            recr(last_real_year-parest_flags(400)+1,last_real_year).
              initialize();
            set_value(recr(max(2,first_unfixed_year),
              last_real_year-parest_flags(400)),x,ii,-bd,bd1,pen);
          }
          indepvars_error();
        }
        if (pmsd && pmsd->num_species>1)
        {
          int ns=pmsd->num_species;
          for (int is=2;is<=ns;is++)
          {
            if (parest_flags(400)==0)  //NMD_19May2016
            {
              mfcl_error();
              if (runsetv==1)
              {
                set_value(pmsd->recr(is),x,ii,-bd,bd1,pen);
              }
              indepvars_error();
            }
            else
            {
              pmsd->recr(is)(last_real_year-parest_flags(400)+1,last_real_year).
                initialize();
              if (runsetv==1)
              {
                set_value(pmsd->recr(is)(max(2,first_unfixed_year),
                  last_real_year-parest_flags(400)),x,ii,-bd,bd1,pen);
              }
              indepvars_error();
            }
          }
        }
      }
    }
    else 
    {
      //recinpop_orth_reset(x,ii);
    }
  }

// overall population scaling
//  jj=ii;
  if (age_flags(32)>0)
  {
    if (parest_flags(155)==0)
    {
      int pf;
      if (parest_flags(146)>0)
      {
        pf=parest_flags(146);
      }
      else
      {
        pf=1;
      }
      if (age_flags(177)==0 && age_flags(92) ==0)
      {
        if (runsetv==1)
        {
          set_value(totpop_coff,x,ii,-5.0,5.0,pen);
        }
        indepvars_error();
        if (pmsd)
        {
          if (runsetv==1)
          {
            set_value(pmsd->totpop_coff,x,ii,-5.0,5.0,pen);
          }
          indepvars_error();
        }
      }
      else
      {
        MY_DOUBLE_TYPE nscal=1.0;
        if (parest_flags(387))  //NMD_4sep2020
        {
          nscal=100.0;
        }
        else
        {
          nscal=5000.0;
        }
        if (!pmsd || pmsd->num_species==1)
        {
          if (runsetv==1)
          {
            set_value(totpop_coff,x,ii,8.0L,25,pen,nscal);
          }
          MY_DOUBLE_TYPE idv;
          idv=value(totpop_coff);
          put_in_indep_vars(idv,8.0L,25.0,indep_var,indep_var_lo,
                            indep_var_hi,jj);
        }
        else
        {
          if (runsetv==1)
          {
            ivector sf1=column(pmsd->species_flags,1);
            dvar_vector tmp(1,pmsd->num_species);
            ivector v(1,pmsd->num_species);
            v=1;
            set_value(tmp,x,ii,8.0L,25,pen,v,sf1,nscal);  //NMD_4sep2020
            totpop_coff=tmp(1);
            pmsd->totpop_coff=tmp(2,pmsd->num_species);
          }
          indepvars_error();
        }
      }
    }
  }

    // *********************************************
    // logvN length for self scaling multinomial
    {
      MY_DOUBLE_TYPE lbd=1.0;
      MY_DOUBLE_TYPE ubd=7.5;
      if (parest_flags(334)>0)
        ubd=parest_flags(334);
      if (runsetv==1)
      {
        set_value(fish_pars(14),x,ii,lbd,ubd,pen,
          column(fish_flags,67),column(fish_flags,68),500);
      }
      dvector idv;
      idv=value(fish_pars(14));
      put_in_indep_vars(idv,lbd,ubd,indep_var,indep_var_lo,indep_var_hi,jj,
        column(fish_flags,67),column(fish_flags,68));
    }


// variance multiplier coff for length dirichlet multinomial
  MY_DOUBLE_TYPE lbd=-7.0;
  MY_DOUBLE_TYPE ubd=7.0;
  if (runsetv==1)
  {
    set_value(fish_pars(22),x,ii,lbd,ubd,pen,
      column(fish_flags,69),column(fish_flags,68),100.);
  }
  dvector idv;
  idv=value(fish_pars(22));
  put_in_indep_vars(idv,lbd,ubd,indep_var,indep_var_lo,indep_var_hi,jj,
    column(fish_flags,69),column(fish_flags,68));

// sample size covariate for length dirichlet multinomial
  if (runsetv==1)
  {
    set_value(fish_pars(23),x,ii,lbd,ubd,pen,
      column(fish_flags,89),column(fish_flags,68),100.);
  }
  idv=value(fish_pars(23));
  put_in_indep_vars(idv,lbd,ubd,indep_var,indep_var_lo,indep_var_hi,jj,
    column(fish_flags,89),column(fish_flags,68));

// length sample size covariate for self scaling multinomial
  {
    MY_DOUBLE_TYPE bd=1.0;
    if (parest_flags(317)>0)
      bd=parest_flags(317)/100.;
    if (runsetv==1)
    {
      set_value(fish_pars(20),x,ii,-0.05,bd,pen,
        column(fish_flags,85),column(fish_flags,68),100.);
    }
    dvector idv;
    idv=value(fish_pars(20));
    put_in_indep_vars(idv,lbd,ubd,indep_var,indep_var_lo,indep_var_hi,jj,
      column(fish_flags,85),column(fish_flags,68));
  }


  if (parest_flags(33)>0)
  {
    MY_DOUBLE_TYPE ub=parest_flags(33)/100.;  
    if (!age_flags(198))
    {
      if (runsetv==1)
      {
        set_value_exp(fish_pars(3),x,ii,.001,ub,pen,
          column(fish_flags,33),column(fish_flags,34),1000.);
      }
      indepvars_error();
    } 
    else
    {
      MY_DOUBLE_TYPE scale=1000.0;
      if (parest_flags(387))  //NMD_4sep2020
      {
        scale=10000.0;
      }
      else
      {
        scale=1000.0;
      }
      if (runsetv==1)
      {
        set_value_exp(tag_fish_rep,x,ii,.001,ub,pen,
          tag_fish_rep_active_flags,tag_fish_rep_group_flags,scale);
      }
      dmatrix idv;
      idv=value(tag_fish_rep);
      put_in_indep_vars(idv,0.001,ub,indep_var,indep_var_lo,
	indep_var_hi,jj,tag_fish_rep_active_flags,tag_fish_rep_group_flags);
    }
     //xofs1  << ii << " ";
    //    cout   << "AJ1 " << ii << " ";
  }
  else
  {          
    if (!age_flags(198))
    {
      if (runsetv==1)
      {
        set_value(fish_pars(3),x,ii,.001,1.01,pen,
          column(fish_flags,33),column(fish_flags,34));
      }
      ivector ivtmp=column(fish_flags,33);
      if (sum(ivtmp)>0.0)
      {
        indepvars_error();
      }
    }
    else
    {
      if (runsetv==1)
      {
        set_value(tag_fish_rep,x,ii,.001,1.01,pen,
          tag_fish_rep_active_flags,tag_fish_rep_group_flags);
      }
      dmatrix idv;
      idv=value(tag_fish_rep);
      put_in_indep_vars(idv,0.001,1.01,indep_var,indep_var_lo,
	indep_var_hi,jj,tag_fish_rep_active_flags,tag_fish_rep_group_flags);
    }
  }  



  if (parest_flags(111) == 4)
  {
    ivector ff43;
    ff43=column(fish_flags,43);
    ivector ff44;
    ff44=column(fish_flags,44);
    if (parest_flags(305)==0)
    {
      if (runsetv==1)
      {
        set_value(fish_pars(4),x,ii,-50.,50.,pen,
          column(fish_flags,43),column(fish_flags,44),10.0);
      }
      if (sum(ff43))
      {
        dvector idv;
        idv=value(fish_pars(4));
        put_in_indep_vars(idv,-50.,50.0,indep_var,indep_var_lo,indep_var_hi,jj,
          ff43,ff44);
      }


//      indepvars_error();
    }
    else
    {
      MY_DOUBLE_TYPE lb=-5.0;
      MY_DOUBLE_TYPE ub=5.0;
      if (parest_flags(306))
      {
#if !defined(NO_MY_DOUBLE_TYPE)
        lb=log(parest_flags(306)/100.-1.0L);
#else
        lb=log(parest_flags(306)/100.-1.0);
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
        ub=log(50.0*parest_flags(306)/100.-1.0L);
#else
        ub=log(50.0*parest_flags(306)/100.-1.0);
#endif
      }
      if (parest_flags(358))
      {
        ub=log(parest_flags(358)/100.);
      }
      if (runsetv==1)
      {
        set_value(fish_pars(4),x,ii,lb,ub,pen,
          column(fish_flags,43),column(fish_flags,44),10.0);
      }
      if (sum(ff43))
      {
        dvector idv;
        idv=value(fish_pars(4));
        put_in_indep_vars(idv,lb,ub,indep_var,indep_var_lo,indep_var_hi,jj,
          ff43,ff44);
      }
    }
  }

  if (runsetv==1)
  {
    set_value(fish_pars(32),x,ii,0.8,5.0,pen,column(fish_flags,93),
      column(fish_flags,29),10.0);
  }
  ivector ff93;
  ff93=column(fish_flags,93);
  ivector ff29;
  ff29=column(fish_flags,29);
  if (sum(ff93))
  {
    dvector idv;
    idv=value(fish_pars(32));
    put_in_indep_vars(idv,0.8,5.0,indep_var,indep_var_lo,indep_var_hi,jj,
      ff93,ff29);
  }

  if (parest_flags(155)==0)
  {
    if (runsetv==1)
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
      set_value(region_pars(1),x,ii,1.e-6,1.5,pen,region_flags(1),
        region_flags(3),nscal);
    }
    if (sum(region_flags(1)))
    {
      dvector idv;
      idv=value(region_pars(1));
      put_in_indep_vars(idv,1.e-06,1.5,indep_var,indep_var_lo,indep_var_hi,jj,
        region_flags(1),region_flags(3));
    }
  }

  if (parest_flags(393))
  {

    if (runsetv==1)
    {
      set_value_partial(kludged_equilib_coffs,
         x, ii, parest_flags(374), -3.0,0.3,pen,100.0);
    }
    d4_array idv;

    int mmin=kludged_equilib_coffs.indexmin();
    int mmax=kludged_equilib_coffs.indexmax();
    idv.allocate(mmin,mmax);
    for (int i=mmin;i<=mmax;i++)
    {
      idv(i)=value(kludged_equilib_coffs(i));
    }
    put_in_indep_vars(idv,-3.0,3.0,indep_var,indep_var_lo,
      indep_var_hi,jj,parest_flags(374));

    if (runsetv==1)
    {
      MY_DOUBLE_TYPE ten10=10.0;
      MY_DOUBLE_TYPE hun100=100.0;
      MY_DOUBLE_TYPE thou1000=1000.0;
      set_value(kludged_equilib_level_coffs,x,ii,
        -3.0,0.3,pen,thou1000);
    }
    dvector idv2;
    idv2=value(kludged_equilib_level_coffs);
    put_in_indep_vars(idv2,-3.0,3.0,indep_var,indep_var_lo,indep_var_hi,jj);
  }

  if (parest_flags(377))
  {
    ivector num_ifmlrp=implicit_fml_bounds(3);
    if (runsetv==1)
    {
      set_value_partial(implicit_fm_level_regression_pars,x,ii,
        num_ifmlrp,-500.0,500.0,pen,ff29,10000.0);
    }
    dmatrix idv;
    idv=value(implicit_fm_level_regression_pars);
    put_in_indep_vars(idv,-500.0,500.0,indep_var,indep_var_lo,indep_var_hi,
      jj,num_ifmlrp,ff29);
  }
  ii=jj;
}

void dvar_len_fish_stock_history::get_indep_vars_for_report(dvar_vector& x)
{
  int ii=1;
  dvariable pen=0.;
  // flag=1 to include the set_value() calls during indepvar placements
  int runsetv=0;

  if (age_flags(30)>0)
  {
    if (parest_flags(155)!=0)
    {
      recinpop_orth_reset_indepvar(x,ii,runsetv);
    }
  }

  dvar_fish_stock_history::get_indep_vars_for_report(x,ii,runsetv);
  
  int jj;
  jj=ii;
  if (!pmsd || sum(column(pmsd->species_flags,2))==0)
  {
    if (age_flags(146)>0  && age_flags(94) !=3)
    {
      if (runsetv==1)
      {
        set_value(sv(21),x,ii,1.e-6,1000.0,pen,1000.);
      }
      MY_DOUBLE_TYPE idv;
      idv=value(sv(21));	 
      put_in_indep_vars(idv,1.e-6,1000.0,indep_var,indep_var_lo,
                        indep_var_hi,jj);
      if (pmsd)
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          if (runsetv==1)
          {
            set_value(pmsd->sv(i,21),x,ii,1.e-6,1000.0,pen,1000.);
          }
          indepvars_error();
        }
      }
    }
    if (age_flags(162)>0)
    {
      if (runsetv==1)
      {
        set_value(sv(29),x,ii,0.201,0.999,pen);
      }
      indepvars_error();
      if (pmsd)
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          if (runsetv==1)
          {
            set_value(pmsd->sv(i,29),x,ii,0.201,0.999,pen);
          }
          indepvars_error();
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
        if (runsetv==1)
        {
          set_value(sv(21),x,ii,1.e-6,1000.0,pen,1000.);
        }
        indepvars_error();	
      }
      int numsp=pmsd->num_species;
      for (int i=2;i<=numsp;i++)
      {
        if (sf2(i))
        {
          if (runsetv==1)
          {	    
            set_value(pmsd->sv(i,21),x,ii,1.e-6,1000.0,pen,1000.);
          }
          indepvars_error();
        }
      }
    }
    if (age_flags(162)>0)
    {
      if (sf2(1))
      {
        if (runsetv==1)
        {
          set_value(sv(29),x,ii,0.201,0.999,pen);
        }
      }
      indepvars_error();
      int numsp=pmsd->num_species;
      for (int i=2;i<=numsp;i++)
      {
        if (sf2(i))
        {
          if (runsetv==1)
          {
            set_value(pmsd->sv(i,29),x,ii,0.201,0.999,pen);
          }
        }
        indepvars_error();
      }
    }
  }

  if (parest_flags(184)>0)
  {
    int num=parest_flags(173)-1;
    if (runsetv==1)
    {
      set_value_partial(age_pars(3),x,ii,num,-50.0,50.0,pen,1000.0);
    }
    dvector idv;
    idv=value(age_pars(3));
    put_in_indep_vars(idv,-50.0,50.0,indep_var,indep_var_lo,indep_var_hi,
      jj,num);
    if (pmsd)
    {
      int numsp=pmsd->num_species;
      for (int is=2;is<=numsp;is++)
      {
        if (runsetv==1)
        {
          set_value_partial(get_age_pars_species(is,3),x,ii,num,
            -50.0,50.0,pen,1000.0);
        }
        indepvars_error();
      }
    }
  }


  if (parest_flags(121)>0)
  {
    int num=parest_flags(121);
    if (runsetv==1)
    {
      MY_DOUBLE_TYPE nscal=1.0;
      if (parest_flags(387))  //NMD_4sep2020
      {
        nscal=4500.0;
      }
      else
      {
        nscal=200.0;
      }
      set_value_partial(age_pars(5),x,ii,num,-20.0,2.0,pen,nscal);
    }
    dvector idv;
    idv=value(age_pars(5));
    put_in_indep_vars(idv,-20.0,2.0,indep_var,indep_var_lo,indep_var_hi,
      jj,num);

    if (pmsd)
    {
      if (runsetv==1)
      {
        int numsp=pmsd->num_species;
        MY_DOUBLE_TYPE nscal=1.0;
        for (int is=2;is<=numsp;is++)
        {
          set_value_partial(get_age_pars_species(is,5),x,ii,num,
            -20.0,2.0,pen,nscal);      //NMD_4sep2020
        }
      }
      indepvars_error();      
    }
  }

  if (parest_flags(12)>0)
  {
    if (runsetv==1)
    {
      vb_coff(1)=boundp(x(ii),fmin1,fmax1,pen,1000.);
      ii=ii+1;
    }
    MY_DOUBLE_TYPE idv;
    idv=value(vb_coff(1));	 
    put_in_indep_vars(idv,fmin1,fmax1,indep_var,indep_var_lo,
      indep_var_hi,jj);
      
    if (pmsd)
    {
      if (runsetv==1)
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          pmsd->vb_coff(i,1)=boundp(x(ii),pmsd->fmin1(i),
            pmsd->fmax1(i),pen,1000.);
          ii=ii+1;
        }
      }
      indepvars_error();      
    }
  }
  
  if (parest_flags(13) > 0)
  {
    if (runsetv==1)
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
      vb_coff(2)=boundp(x(ii),fminl,fmaxl,pen,nscal);
      ii=ii+1;
    }
    MY_DOUBLE_TYPE idv;
    idv=value(vb_coff(2));	 
    put_in_indep_vars(idv,fminl,fmaxl,indep_var,indep_var_lo,
      indep_var_hi,jj);
     
    if (pmsd)
    {
      if (runsetv==1)
      {
        int numsp=pmsd->num_species;
        MY_DOUBLE_TYPE nscal=1.0;
        for (int i=2;i<=numsp;i++)
        {
          pmsd->vb_coff(i,2)=boundp(x(ii),pmsd->fminl(i),
            pmsd->fmaxl(i),pen,nscal);  //NMD_4sep2020
          ii=ii+1;
        }
      }
      indepvars_error();      
    }
  }

  if (parest_flags(14) > 0)
  {
    if (runsetv==1)
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
      ii=ii+1;
    }
    MY_DOUBLE_TYPE idv;
    idv=value(vb_coff(3));	 
    put_in_indep_vars(idv,rhomin,rhomax,indep_var,indep_var_lo,
      indep_var_hi,jj);
    if (pmsd)
    {
      if (runsetv==1)
      {
        int numsp=pmsd->num_species;
        MY_DOUBLE_TYPE nscal=1.0;
        for (int i=2;i<=numsp;i++)
        {
          pmsd->vb_coff(i,3)=boundp(x(ii),pmsd->rhomin(i),
            pmsd->rhomax(i),pen,nscal);      //NMD_4sep2020

          ii=ii+1;
        }
      }
      indepvars_error();      
    }
  }

  if (parest_flags(227) > 0)
  {
    if (runsetv==1)
    {
      vb_coff(4)=boundp(x(ii),-5.0,5.0,pen,100.);
      ii=ii+1;
    }
    MY_DOUBLE_TYPE idv;
    idv=value(vb_coff(4));	 
    put_in_indep_vars(idv,-5.0,5.0,indep_var,indep_var_lo,
      indep_var_hi,jj);
    if (pmsd)  //NMD_1Apr2020
    {
      if (runsetv==1)
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          pmsd->vb_coff(i,4)=boundp(x(ii),-5.0,5.0,pen,100.);
          ii=ii+1;
        }
      }
      indepvars_error();      
    }  //NMD_1Apr2020
  }
  
  if (parest_flags(15) > 0)
  {
    if (runsetv==1)
    {
      var_coff(1)=boundp(x(ii),vmin1,vmax1,pen,1000.);
      ii=ii+1;
    }
    MY_DOUBLE_TYPE idv;
    idv=value(var_coff(1));	 
    put_in_indep_vars(idv,vmin1,vmax1,indep_var,indep_var_lo,
      indep_var_hi,jj);
    if (pmsd)
    {
      if (runsetv==1)
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          pmsd->var_coff(i,1)=boundp(x(ii),pmsd->vmin1(i),pmsd->vmax1(i),
            pen,1000.);
          ii=ii+1;
        }
      }
      indepvars_error();    
    }
  }
  
  if (parest_flags(16) > 0)
  {
    if (runsetv==1)
    {
      var_coff(2)=boundp(x(ii),vminl,vmaxl,pen,1000.);
      ii=ii+1;
    }
    MY_DOUBLE_TYPE idv;
    idv=value(var_coff(2));	 
    put_in_indep_vars(idv,vminl,vmaxl,indep_var,indep_var_lo,
      indep_var_hi,jj);
    if (pmsd)
    {
      if (runsetv==1)
      {
        int numsp=pmsd->num_species;
        for (int i=2;i<=numsp;i++)
        {
          pmsd->var_coff(i,2)=boundp(x(ii),pmsd->vminl(i),pmsd->vmaxl(i),
            pen,1000.);
          ii=ii+1;
        }
      }
      indepvars_error();      
    }
  }
//  cout << "jj= " << jj << endl;
  int real_nvar;
  real_nvar=indep_var.indexmax();
  if (real_nvar != jj-1)
    {
      cerr << "Indepvars report number of parameters doesn't match xinit.rpt "
           << " real_nvar = " << real_nvar << "  indep_vars max = "
           << jj-1 << "  exiting..." << endl;
      ad_exit(1);
    }
//  ad_exit(1);
}

void dvar_fish_stock_history::recinpop_orth_reset_indepvar(dvar_vector& x,
       int& ii, int& runsetv)
{
  int jj=ii;
  
  if (parest_flags(155)<0)
  {
    indepvars_error();      
    if (!pmsd || pmsd->num_species==1)
    {
      set_value(new_orth_recr,x,ii);
    }
    else 
    {
      ivector sf2=column(pmsd->species_flags,2);
      if (sum(sf2)==0)
      {
        set_value(new_orth_recr,x,ii);
        for (int is=2;is<=pmsd->num_species;is++)
        {
          set_value(pmsd->new_orth_recr(is),x,ii);
        }
      }
      else
      {
        pmsd_error();
        if (sf2(1))
        {
          set_value(new_orth_recr,x,ii);
        }
        for (int is=2;is<=pmsd->num_species;is++)
        {
          if (sf2(is))
            set_value(pmsd->new_orth_recr(is),x,ii);
        }
      }
    }
  }
  else
  {
    int ns=1;
    if (pmsd) ns=pmsd->num_species;
    for (int is=1;is<=ns;is++)
    {
      if (is==1)
      {
        int ir;
        int jjj=0;
        for (int i=1;i<=numcomp(1);i++)
        {
          if (degree_yr>0)
          {
            if (runsetv==1)
            {
              set_value_partial(orth_recr_all(++jjj),x,ii,degree_yr,100.);
            }else{
              ++jjj;
            }
            int num=degree_yr;
            dvector tmp;
            tmp=value(orth_recr_all(jjj));
            dvector idv;
            idv=tmp.shift(1);
            put_in_indep_vars(idv,-999.0,999.0,indep_var,indep_var_lo,
              indep_var_hi,jj,num);
          }
          else
          {
            indepvars_error();
            int by=1+ny_begin_yr-1;
            int ey=num_real_years-ny_end_yr+1;
            if (runsetv==1)
            {
              set_value_partial(yearly_recr_all(++jjj),x,ii,ey-by+2,100.);
            }else{
              ++jjj;
            }
          }
        }
  
        for (int i=1;i<=numcomp(2);i++)
        {
          if (runsetv==1)
          {
            set_value_partial(orth_recr_all(++jjj),x,ii,degree_ses,100.);
          }else{
            ++jjj;
          }
	  int num=degree_ses;
          dvector tmp;
          tmp=value(orth_recr_all(jjj));
          dvector idv;
          idv=tmp.shift(1);
          put_in_indep_vars(idv,-999.0,999.0,indep_var,indep_var_lo,
            indep_var_hi,jj,num);
        }
  
        for (int i=1;i<=numcomp(3);i++)
        {
          if (runsetv==1)
          {
            set_value_partial(orth_recr_all(++jjj),x,ii,degree_reg,100.);
          }else{
            ++jjj;
          }
	  int num=degree_reg;
          dvector tmp;
          tmp=value(orth_recr_all(jjj));
          dvector idv;
          idv=tmp.shift(1);
          put_in_indep_vars(idv,-999.0,999.0,indep_var,indep_var_lo,
            indep_var_hi,jj,num);
        }
  
        for (int i=1;i<=numcomp(4);i++)
        {
          if (runsetv==1)
          {
            set_value_partial(orth_recr_all(++jjj),x,ii,degree_ses_reg,100.);
          }else{
            ++jjj;
          }
	  int num=degree_ses_reg;
          dvector tmp;
          tmp=value(orth_recr_all(jjj));
          dvector idv;
          idv=tmp.shift(1);
          put_in_indep_vars(idv,-999.0,999.0,indep_var,indep_var_lo,
            indep_var_hi,jj,num);
        }
      }
      else
      {
        ivector numcomp=pmsd->numcomp(is);
        ivector sf2=column(pmsd->species_flags,2);
        int degree_yr=pmsd->degree_yr(is);
        int degree_ses=pmsd->degree_ses(is);
        int degree_reg=pmsd->degree_reg(is);
        int degree_ses_reg=pmsd->degree_ses_reg(is);
        if (sum(sf2)==0)
        {
          int ir;
          int jjj=0;
          for (int i=1;i<=numcomp(1);i++)
          {
            if (degree_yr>0)
            {
              if (runsetv==1)
              {
                set_value_partial(pmsd->orth_recr_all(is)(++jjj),x,ii,
                  degree_yr,100.);
              }else{
                ++jjj;
              }
              int num=degree_yr;
              dvector tmp;
              tmp=value(orth_recr_all(is)(jjj));
              dvector idv;
              idv=tmp.shift(1);
              put_in_indep_vars(idv,-999.0,999.0,indep_var,indep_var_lo,
                indep_var_hi,jj,num);
            }
            else
            {
              pmsd_error();
              int by=1+ny_begin_yr-1;
              int ey=num_real_years-ny_end_yr+1;
              set_value_partial(yearly_recr_all(++jjj),x,ii,ey-by+2,10.);
            }
          }
    
          for (int i=1;i<=numcomp(2);i++)
          {
            if (runsetv==1)
            {
              set_value_partial(pmsd->orth_recr_all(is)(++jjj),x,ii,degree_ses,
                100);
            }else{
              ++jjj;
            }
            int num=degree_ses;
            dvector tmp;
            tmp=value(orth_recr_all(is)(jjj));
            dvector idv;
            idv=tmp.shift(1);
            put_in_indep_vars(idv,-999.0,999.0,indep_var,indep_var_lo,
              indep_var_hi,jj,num);
          }
    
          for (int i=1;i<=numcomp(3);i++)
          {
            if (runsetv==1)
            {
              set_value_partial(pmsd->orth_recr_all(is)(++jjj),x,ii,degree_reg,
                100);
            }else{
              ++jjj;
            }
            int num=degree_reg;
            dvector tmp;
            tmp=value(orth_recr_all(is)(jjj));
            dvector idv;
            idv=tmp.shift(1);
            put_in_indep_vars(idv,-999.0,999.0,indep_var,indep_var_lo,
              indep_var_hi,jj,num);
          }
    
          for (int i=1;i<=numcomp(4);i++)
          {
            if (runsetv==1)
            {
              set_value_partial(pmsd->orth_recr_all(is)(++jjj),
                x,ii,degree_ses_reg,100.);
            }else{
              ++jjj;
            }
            int num=degree_ses_reg;
            dvector tmp;
            tmp=value(orth_recr_all(is)(jjj));
            dvector idv;
            idv=tmp.shift(1);
            put_in_indep_vars(idv,-999.0,999.0,indep_var,indep_var_lo,
              indep_var_hi,jj,num);
          }
        }
        else
        {
          // do nothng here recruitments are assigned in another routine
          //pmsd_error();
        }
      }
    }
  }
  ii=jj;
}




void put_in_indep_vars(MY_DOUBLE_TYPE idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
  dvector indep_var_hi, int& jjsave)
{
  int i1=jjsave;
  indep_var(i1)=idv;
  indep_var_lo(i1)=fmin;
  indep_var_hi(i1)=fmax;
  i1++;
  jjsave=i1;
}

void put_in_indep_vars(dvector idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
  dvector indep_var_hi, int& jjsave)
{
//  int jjsave=jj;
  int i1=jjsave;
  int lmin=idv.indexmin();
  int lmax=idv.indexmax();
  int i2=i1+lmax-1;
  dvector tmp(lmin,lmax);
  tmp.initialize();
  tmp=idv(lmin,lmax);
  indep_var.sub(i1,i2).shift(1)=tmp;
  indep_var_lo.sub(i1,i2).shift(1)=fmin;
  indep_var_hi.sub(i1,i2).shift(1)=fmax;
  i1=i2+1;

  int itot;
  itot=i1-jjsave;
  jjsave=i1;
}



void put_in_indep_vars(dvector idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
  dvector indep_var_hi, int& jjsave, int num)
{
//  int jjsave=jj;
  int i1=jjsave;
  int lmin=idv.indexmin();
  int lmax=lmin+num-1;
  int i2=i1+lmax-1;
  dvector tmp(lmin,lmax);
  tmp.initialize();
  tmp=idv(lmin,lmax);
  indep_var.sub(i1,i2).shift(1)=tmp;
  indep_var_lo.sub(i1,i2).shift(1)=fmin;
  indep_var_hi.sub(i1,i2).shift(1)=fmax;
  i1=i2+1;

  int itot;
  itot=i1-jjsave;
  jjsave=i1;
}


void put_in_indep_vars(dvector idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
  dvector indep_var_hi, int& jjsave, imatrix fshgptr)
{
//  int jjsave=jj;
  int i1=jjsave;

  int mmin=fshgptr.indexmin();
  int mmax=fshgptr.indexmax();
  for (int i=mmin;i<=mmax;i++)    // loop over the groups
  {
    indep_var(i1)=idv(fshgptr(i,1));
    indep_var_lo(i1)=fmin;
    indep_var_hi(i1)=fmax;
    i1=i1+1;
  }
  
  int itot;
  itot=i1-jjsave;
  jjsave=i1;
}

void put_in_indep_vars(dvector idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
  dvector indep_var_hi, int& jjsave, ivector fflgs, ivector fgrps)
{
//  int jjsave=jj;
  dvector w;
  w=idv;
  int i1=jjsave;
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  if (sum(fflgs))
  {
    if (sum(fgrps))
    {
      ivector key(mmin,mmax);
      sort(fgrps,key);
      int kin=key(mmin);
      int flag_value=fflgs(kin);
      if (fflgs(kin))
      {
        indep_var(i1)=w(kin);
        indep_var_lo(i1)=fmin;
        indep_var_hi(i1)=fmax;
        i1=i1+1;
      }
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (fgrps(key(i))==fgrps(key(i-1)))
        {
          if (fflgs(key(i))!=flag_value)
          {
            cerr << "Error -- grouped initial parameters have unequal flags"
                 << endl;
            ad_exit(1);
          }
          if (fflgs(key(i))) w(key(i))=w(kin);
        }
        else
        {
          kin=key(i);
          flag_value=fflgs(kin);
          if (fflgs(kin))
          {
            indep_var(i1)=w(kin);
            indep_var_lo(i1)=fmin;
            indep_var_hi(i1)=fmax;
            i1=i1+1;
          }
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (fflgs(i))
        {
          indep_var(i1)=w(i);
          indep_var_lo(i1)=fmin;
          indep_var_hi(i1)=fmax;
          i1=i1+1;
        }
      }
    }
  }

  int itot;
  itot=i1-jjsave;
  jjsave=i1;
}


void put_in_indep_vars(dmatrix idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
		       dvector indep_var_hi, int& jjsave)
{
//  int jjsave=jj;
  cout << jjsave << endl;
  int mmin=idv.indexmin();
  int mmax=idv.indexmax();
  int i1=jjsave;
  for (int i=mmin;i<=mmax;i++)
  {
    int lmin=idv(i).indexmin();
    int lmax=idv(i).indexmax();
    int i2=i1+lmax-1;
    dvector tmp(lmin,lmax);
    tmp.initialize();
    tmp=idv(i)(lmin,lmax);
    indep_var.sub(i1,i2).shift(1)=tmp;
    indep_var_lo.sub(i1,i2).shift(1)=fmin;
    indep_var_hi.sub(i1,i2).shift(1)=fmax;
    i1=i2+1;
  }
  int itot;
  itot=i1-jjsave;
  jjsave=i1;
  cout << jjsave << endl;
}

void put_in_indep_vars(dmatrix idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
  dvector indep_var_hi, int& jjsave, imatrix fflgs, imatrix fgrps)
{
//  int jjsave=jj;

  dvector w=rowstack(idv);
  ivector flags=rowstack(fflgs);
  ivector group=rowstack(fgrps);
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  int i1=jjsave;
  if (sum(flags))
  {
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      int flag_value=flags(kin);
      if (flags(kin))
      {
        indep_var(i1)=w(kin);
        indep_var_lo(i1)=fmin;
        indep_var_hi(i1)=fmax;
        i1=i1+1;
      }
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          if (flags(key(i))) w(key(i))=w(kin);
          if (flags(key(i))!=flag_value)
          {
            cerr << "Error -- grouped initial parameters have unequal flags"
                 << endl;
            exit(1);
          }
        }
        else
        {
          kin=key(i);
          flag_value=flags(kin);
          if (flags(kin))
          {
            indep_var(i1)=w(kin);
            indep_var_lo(i1)=fmin;
            indep_var_hi(i1)=fmax;
            i1=i1+1;
          }
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i))
        {
          indep_var(i1)=w(i);
          indep_var_lo(i1)=fmin;
          indep_var_hi(i1)=fmax;
          i1=i1+1;
        }
      }
    }
  }
  int itot;
  itot=i1-jjsave;
  jjsave=i1;
}

void put_in_indep_vars(dmatrix idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
  dvector indep_var_hi, int& jjsave,ivector n,ivector fgrps)
{
//  int jjsave=jj;
  int mmin=idv.indexmin();
  int mmax=idv.indexmax();
  int i1=jjsave;
  ivector ig(mmin,mmax);
  if (sum(fgrps))
  {
    int brkflg=0;
    ig(1)=fgrps(1);
    if (n(1)>0)
    {
      int lmin=idv(1).indexmin();
      int lmax=lmin+n(1)-1;
      int i2=i1+lmax-1;
      dvector tmp(lmin,lmax);
      tmp.initialize();
      tmp=idv(1)(lmin,lmax);
      indep_var.sub(i1,i2).shift(1)=tmp;
      indep_var_lo.sub(i1,i2).shift(1)=fmin;
      indep_var_hi.sub(i1,i2).shift(1)=fmax;
      i1=i2+1;
    }
    for (int i=mmin+1;i<=mmax;i++)
    {
      brkflg=0;
      ig(i)=fgrps(i);
      for (int j=mmin;j<=(i-1);j++)
      {
          if (ig(i)==ig(j)) brkflg=1;
      }
      if (!brkflg)
      {
        if (n(i)>0)
        {
          int lmin=idv(i).indexmin();
          int lmax=lmin+n(i)-1;
          int i2=i1+lmax-1;
          dvector tmp(lmin,lmax);
          tmp.initialize();
          tmp=idv(i)(lmin,lmax);
          indep_var.sub(i1,i2).shift(1)=tmp;
          indep_var_lo.sub(i1,i2).shift(1)=fmin;
          indep_var_hi.sub(i1,i2).shift(1)=fmax;
          i1=i2+1;
        }
      }
    }
  }
  else
  {
    for (int i=mmin;i<=mmax;i++)
    {
      if (n(i)>0)
      {
        int lmin=idv(i).indexmin();
        int lmax=lmin+n(i)-1;
        int i2=i1+lmax-1;
        dvector tmp(lmin,lmax);
        tmp.initialize();
        tmp=idv(i)(lmin,lmax);
        indep_var.sub(i1,i2).shift(1)=tmp;
        indep_var_lo.sub(i1,i2).shift(1)=fmin;
        indep_var_hi.sub(i1,i2).shift(1)=fmax;
        i1=i2+1;
      }
    }
  }
  int itot;
  itot=i1-jjsave;
  jjsave=i1;
}


void put_in_indep_vars(dmatrix idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
  dvector indep_var_hi, int& jjsave, ivector ffgroup)
{
  cerr << "code not finished - exiting" << endl;
  ad_exit(1);
//  int jjsave=jj;
  int mmin=idv.indexmin();
  int mmax=idv.indexmax();
  int i1=jjsave;
  for (int i=mmin;i<=mmax;i++)
  {
    int lmin=idv(i).indexmin();
    int lmax=idv(i).indexmax();
    int i2=i1+lmax-1;
    dvector tmp(lmin,lmax);
    tmp.initialize();
    tmp=idv(i)(lmin,lmax);
    indep_var.sub(i1,i2).shift(1)=tmp;
    indep_var_lo.sub(i1,i2).shift(1)=fmin;
    indep_var_hi.sub(i1,i2).shift(1)=fmax;
    i1=i2+1;
  }
  int itot;
  itot=i1-jjsave;
  jjsave=i1;
}


void put_in_indep_vars(d4_array idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
		       dvector indep_var_hi, int& jjsave)
{
//  int jjsave=jj;
  int mmin=idv.indexmin();
  int mmax=idv.indexmax();
  int i1=jjsave;
  for (int i=mmin;i<=mmax;i++)
  {
    int jmin=idv(i).indexmin();
    int jmax=idv(i).indexmax();
    for (int j=jmin;j<=jmax;j++)
    {
      int kmin=idv(i,j).indexmin();
      int kmax=idv(i,j).indexmax();
      for (int k=kmin;k<=kmax;k++)
      {
        int lmin=idv(i,j,k).indexmin();
        int lmax=idv(i,j,k).indexmax();

        int i2=i1+lmax-1;
        dvector tmp(lmin,lmax);
        tmp.initialize();
        tmp=idv(i,j,k)(lmin,lmax);
        indep_var.sub(i1,i2).shift(1)=tmp;
        indep_var_lo.sub(i1,i2).shift(1)=fmin;
        indep_var_hi.sub(i1,i2).shift(1)=fmax;
        i1=i2+1;
      }
    }
  }
  int itot;
  itot=i1-jjsave;
  jjsave=i1;
}

void put_in_indep_vars(d4_array idv,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,dvector indep_var,dvector indep_var_lo,
  dvector indep_var_hi, int& jjsave,MY_DOUBLE_TYPE n)
{
//  int jjsave=jj;
  int mmin=idv.indexmin();
  int mmax=idv.indexmax();
  int i1=jjsave;
  for (int i=mmin;i<=mmax;i++)
  {
    int jmin=idv(i).indexmin();
    int jmax=idv(i).indexmax();
    for (int j=jmin;j<=jmax;j++)
    {
      int kmin=idv(i,j).indexmin();
      int kmax=idv(i,j).indexmax();
      for (int k=kmin;k<=kmax;k++)
      {
        int lmin=idv(i,j,k).indexmin();
        int lmax=lmin+n-1;

        int i2=i1+lmax-1;
        dvector tmp(lmin,lmax);
        tmp.initialize();
        tmp=idv(i,j,k)(lmin,lmax);
        indep_var.sub(i1,i2).shift(1)=tmp;
        indep_var_lo.sub(i1,i2).shift(1)=fmin;
        indep_var_hi.sub(i1,i2).shift(1)=fmax;
        i1=i2+1;
      }
    }
  }
  int itot;
  itot=i1-jjsave;
  jjsave=i1;
}
