/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
#include "scbd.hpp"

const int REGION_REC_SCALE=100.;

static void scrout_total(int newtotal)
{
  cout <<"newmau5a.cpp " << newtotal << endl;
}
int getsize(const dvar3_array& s);
void mykludge(int ii1){ii1=0;}

dvariable dvar_fish_stock_history::reset(dvar_vector& x, int& ii,
  d3_array& len_sample_size)
{
  //ofstream xofs1("counts");
  int * pii = &ii;
  dvariable pen=0.;
  //ii=1;
 //if (parest_flags(392)==0)
 //{
 

  //if (parest_flags(361)>0)
  //  set_value(testpar,x,ii);

  if (parest_flags(392)==0)   //NMD_7jun2024
  {

    if (age_flags(68)>0)
    {
      switch(age_flags(184))
      {
      case 0:
        {
          MY_DOUBLE_TYPE nscal=1.0;
          /*
          if (!parest_flags(387))
          {
            nscal=300.0;
          }
          else
          {
            nscal=1000.0;
          }
          */
           nscal=1000.0;
          set_value(diff_coffs,x,ii,0.0,3.0,pen,nscal,age_flags(68),move_map);
        }  
        break;
      case 1:
        {
          set_value(xdiff_coffs,x,ii,-200.0,200.0,pen,1.e+3,age_flags(68),move_map);
        }  
        break;
      case 2:
        set_value(zdiff_coffs,x,ii,-10.0,10.0,pen,1.e+2,age_flags(68),move_map);
        break;
      default:
        cerr << "Illegal value for age_flags(184)" << endl;
        ad_exit(1);
      }
    }

    if (age_flags(88)>0)
    {
      set_value(diff_coffs2,x,ii,-4.0,4.0,pen,1000.,age_flags(88),move_map);
    }

    if (age_flags(90)>0)
    {
      set_value(diff_coffs3,x,ii,-4.0,4.0,pen,1000.,age_flags(90),move_map);
    }

//  if (parest_flags(392)==0)
//  {

    if (age_flags(70)>0)
    {
      if (parest_flags(155)==0)
      {
        if (num_regions>1)
        {
          ivector iv(1,num_regions);
          iv=1;
          dvar_matrix tmp(1,num_regions,first_unfixed_year,last_real_year-1);
        //set_value(region_rec_diff_coffs.sub(1,last_real_year-1),x,ii,-3.,3.,pen,
          set_value(tmp,x,ii,-3.,3.,pen,REGION_REC_SCALE,1,region_flags(3));
          region_rec_diff_coffs.sub(first_unfixed_year,last_real_year-1)=
            trans(tmp);
        }
      }
    }

    //ofstream xxofs("effdev.rpt");
    int i;
    if (age_flags(34)>0 && !age_flags(143) && !age_flags(92))
    {
      MY_DOUBLE_TYPE max_dev;
      if (age_flags(35)==0)
      {
        max_dev=1;
      }
      else
      {
        max_dev=age_flags(35);
      }
      if (!sum_ff79_flag)
      {
        for (i=1;i<=num_fisheries;i++)
        {
          if (fish_flags(i,4)>1)
          {
            dvar_vector edc=effort_dev_coffs(i)(first_unfixed_fish_time(i),
              num_real_fish_times(i));
            edc.initialize();
            set_value(edc,zero_effdev_flag(i)(first_unfixed_fish_time(i),
              num_real_fish_times(i)),x,ii,-max_dev,max_dev,pen,100.);
          }
        }
      }
      else
      {
        ivector ff79=column(fish_flags,79);
        ivector ff4=column(fish_flags,4);
        int projflag=0;
        for (i=1;i<=num_fisheries;i++)
        {
          if (num_real_fish_times(i)!=num_fish_times(i))
          {
            projflag=1;
            break;
          }
        }
        if (projflag)
        {
          dvar_matrix sub_effort_dev_coffs(1,num_fisheries);
          for (i=1;i<=num_fisheries;i++)
          {
            sub_effort_dev_coffs(i)=effort_dev_coffs(i)(1,num_real_fish_times(i));
          }
          set_value(sub_effort_dev_coffs,x,ii,-max_dev,max_dev,pen,100.,
             ff4,zero_effdev_flag,ff79);
        }
        else
        {
          int ii1=ii;
          // this was called to check that ii1=ii using new pgroup manager
          set_value(effort_dev_coffs,x,ii1,-max_dev,max_dev,pen,100.,
             ff4,
             //zero_effdev_flag,
             missing_effort_by_realization_flag,
             ff79);

          set_value(effort_dev_coffs,x,ii,-max_dev,max_dev,pen,100.,
            pgroup_manager_1);

          //set_value(effort_dev_coffs,x,ii,-max_dev,max_dev,pen,100.,
            // ff4,zero_effdev_flag,ff79,pgroup_manager_1);

        }
      }
    }
    //xxofs.close();
    //exit(1);

    //removed this line feb07 05 DF
    if (!age_flags(92))
    {
      set_value(q0,x,ii,-25.,1.,pen,column(fish_flags,1),
        column(fish_flags,60),10000.0);
    }
    else if (missing_catch_flag)
    {
      // need to do this by fishery group
      set_value(q0_miss,x,ii,-25.,1.,pen,column(fish_flags,1),
        column(fish_flags,60),10000.0);

      if (age_flags(104)>0)
      {
        set_value(fm_level_devs,x,ii,-2.,2.,pen,100.0);
      }
    }
    // TTTTTTTTTTTTTTTTTT DIST
    if (parest_flags(290)>0)
      set_value(log_length_variance,x,ii,-2.,5.,pen,100.0);

    if (parest_flags(310)>0)
       set_value(length_tot_exp,x,ii,-5.0,5.0,pen,100.);

    if (parest_flags(291)>0)
    {
      switch (parest_flags(297))
      {
      case 0:
      case 1:  //LN3m
        set_value(length_rho,x,ii,-0.85,0.85,pen,100.0);
        break;
      case 2:  //LN3
        set_value(length_rho,x,ii,-1.98,1.98,pen,100.0);
        break;
      default:
        cerr << "error in pf297" << endl;
        ad_exit(1);
      }
    }


    if (parest_flags(292)>0)
      set_value(log_length_dof,x,ii,0.0,10,pen,100.0);

    if (parest_flags(296)>0)
    {
      switch (parest_flags(297))
      {
      case 0:
        set_value(length_psi,x,ii,-0.0001,0.99,pen,100.0);
        break;
      case 1: // LN3m
        set_value(length_psi,x,ii,100.0);
        break;
      case 2: // LN3
        set_value(length_psi,x,ii,0.01,0.99,pen,100.0);
        break;
      default:
        cerr << "error in pf297" << endl;
        ad_exit(1);
      }
    }

    if (parest_flags(298)>0)
      set_value(length_exp,x,ii,-1.98,1.98,pen,100.0);

    // TTTTTTTTTTTTTTTTTT DIST
    if (parest_flags(280)>0)
      set_value(log_weight_variance,x,ii,-10.,15.,pen,100.0);

    if (parest_flags(300)>0)
      set_value(weight_tot_exp,x,ii,-5.0,5.0,pen,100.);

    if (age_flags(134)>0)
    {
      MY_DOUBLE_TYPE ten=10.;
      set_value(species_pars(1),x,ii,-0.975,0.975,pen,ten);
    }
    if (age_flags(136)>0)
    {
      MY_DOUBLE_TYPE ten=10.;
      set_value(species_pars(2),x,ii,-0.975,0.975,pen,ten);
    }


    if (parest_flags(281)>0)
    {
      switch (parest_flags(287))
      {
      case 0:
      case 1:  //LN3m
        set_value(weight_rho,x,ii,-0.90,0.90,pen,100.0);
        break;
      case 2:  //LN3
        set_value(weight_rho,x,ii,-1.98,1.98,pen,100.0);
        break;
      default:
        cerr << "error in pf287" << endl;
        ad_exit(1);
      }
    }

    if (parest_flags(282)>0)
      set_value(log_weight_dof,x,ii,0.0,10,pen,100.0);

    if (parest_flags(286)>0)
    {
      switch (parest_flags(287))
      {
      case 0:
        set_value(weight_psi,x,ii,-0.0001,0.99,pen,100.0);
        break;
      case 1: // LN3m
        set_value(weight_psi,x,ii,100.0);
        break;
      case 2:
        set_value(weight_psi,x,ii,0.01,0.99,pen,100.0);
        break;
      default:
        cerr << "error in pf287" << endl;
        ad_exit(1);
      }
    }
    if (parest_flags(288)>0)
      set_value(weight_exp,x,ii,-1.98,1.98,pen,100.0);

    for (i=1;i<=num_fisheries;i++)
    {
      // time trend in catchability
      if (fish_flags(i,2)>1)
      {
        q1(i)=boundp(x(ii++),-.99,.99,pen);
      }
    }

   // Don't think we need this any more with new bs_selcoff stuff
   // fishing_selectivity_interface(x,ii,pen);

     //xofs1  << ii << " ";
    //  cout   << "AA " << ii << " ";
    dvector tmp(1,num_fisheries);
    tmp.initialize();
  //  cout << " ZZZZ0 " << ii << endl;
    if (parest_flags(323)==0)
    {
      new_fishing_selectivity_interface_reset(ii,x,pen,-3.0,3.0,100.);
    }

  //  cout << " ZZZZ2 " << ii << " " << getsize(sel_dev_coffs) << endl;

    for (i=1;i<=num_fisheries;i++)
    {
      if (fish_flags(i,19)>0 and fish_flags(i,71)>0)
      {
        cerr << "cant have fish_flags(19)>0 and fish_flags(71)>0" << endl;
//        ad_exit(1);     /////////////////CHECK THIS
      }
    }

    if (sum_ff48_flag)
    { 
      ivector ff48=column(fish_flags,48);
      ivector ffgroup=column(fish_flags,24);
      MY_DOUBLE_TYPE scale=1000.0;
      set_value(bs_selcoff,x,ii,-20.,7.,pen,scale ,ff48,ffgroup);
    }

     //xofs1  << ii << " ";
    //  cout   << "AB " << ii << " ";
    for (i=1;i<=num_fisheries;i++)
    {
      // selectivity deviations
      if (fish_flags(i,5)>0)
      {
        MY_DOUBLE_TYPE max_dev;
        if (age_flags(36)==0)
        {
          max_dev=1;
        }
        else
        {
          max_dev=age_flags(36);
        }
        if (fish_flags(i,3)>0)
        {
          MY_DOUBLE_TYPE delt_scale=10.;
          if (parest_flags(151)>0)
          {
            delt_scale=parest_flags(151);
          }
          if (parest_flags(151)<0)
          {
            delt_scale=-1./parest_flags(151);
          }
          for (int j=2;j<=num_fish_times(i);j++)
          {
            if (len_sample_size(realization_region(i,j),realization_period(i,j),
              realization_incident(i,j))>0)
            {
              set_value(delta2(i,j),x,ii,-max_dev,max_dev,pen,delt_scale);
            }
          }
        }
      }
    }
    for (i=1;i<=num_fisheries;i++)
    {
      if (fish_flags(i,6)>0)
      {
        corr_wy(i)=boundp(x(ii++),-.95,.95,pen);
      }
      if (fish_flags(i,7)>0)
      {
        corr_by(i)=boundp(x(ii++),-.95,.95,pen);
      }
      if (fish_flags(i,8)>0)
      {
        corr_wc(i)=boundp(x(ii++),-.95,.95,pen);
      }
      if (fish_flags(i,9)>0)
      {
        corr_eff(i)=boundp(x(ii++),-.95,.95,pen);
      }
    }

     //xofs1  << ii << " ";
    //  cout   << "AC " << ii << " ";
    if (!sum(column(fish_flags,29)))
    {
      for (i=1;i<=num_fisheries;i++)
      {
        if (!age_flags(92))
        {
          if (fish_flags(i,10)>0)
          {
            if (fish_flags(i,23)==0)
            { 
              MY_DOUBLE_TYPE oneh=100.;
              set_value(catch_dev_coffs(i),x,ii,-.8,.8,pen,oneh);
            }
            else
            {
              MY_DOUBLE_TYPE oneh=100.;
              set_value(catch_dev_coffs(i),between_times(i),x,ii,-.8,.8,pen,oneh);
            }
          }
        }
      }
    }
    else
    {
      for (i=1;i<=ngroups;i++)
      {
        if (!age_flags(92))
        {
          if (fish_flags(gfish_index(i,1),10)>0)
          {
            {
              set_value(grouped_catch_dev_coffs(i),grouped_between_times(i),
                 x,ii,-.8,.8,pen,1000.);
            }
          }
        }
      }
    }

     //xofs1  << ii << " ";
    //  cout   << "AD " << ii << " ";

    for (i=1;i<=num_fisheries;i++)
    {
      {
        if (fish_flags(i,37)>0)
        {
          if (fish_flags(i,45)==0)
          { 
            set_value(rep_dev_coffs(i),x,ii,-.5,.5,pen);
          }
          else
          {
            set_value(rep_dev_coffs(i),between_times(i),x,ii,-.5,.5,pen);
          }
        }
      }
    }

    // relative ecruitment
    ivector sf=season_flags(1);
    if (age_flags(30)>0)
    {
      if (parest_flags(155)==0)
      {
        if (age_flags(51)==0)
        {
          //set_value(recr,x,ii);
          if (sum(sf))
          {
  //          set_value(recr(1,last_real_year),x,ii,-20.,20,pen,100,year_flags(1),-10.0);
            if (parest_flags(400)==0)  //NMD_19May2016
            {
              set_value(recr(max(2,first_unfixed_year),last_real_year),
                 x,ii,-20.,20,pen,100,year_flags(1),-10.0);
            }
            else
            {
              recr(last_real_year-parest_flags(400)+1,last_real_year).
                initialize();
              set_value(recr(max(2,first_unfixed_year),
                last_real_year-parest_flags(400)),
                x,ii,-20.,20,pen,100,year_flags(1),-10.0);
            }  //NMD_19May2016
            if (pmsd && pmsd->num_species>1)
            {
              int ns=pmsd->num_species;
              for (int is=2;is<=ns;is++)
              {
  //              set_value(pmsd->recr(is)(1,last_real_year),x,ii,-20.,20,pen,
  //                100,year_flags(1),-10.0);
                if (parest_flags(400)==0)   //NMD_19May2016
                {
                  set_value(pmsd->recr(is)(max(2,first_unfixed_year),
                   last_real_year),x,ii,-20.,20,pen,
                    100,year_flags(1),-10.0);
                }
                else
                {
                  pmsd->recr(is)(last_real_year-parest_flags(400)+1,last_real_year).
                    initialize();
                  set_value(pmsd->recr(is)(max(2,first_unfixed_year),
                    last_real_year-parest_flags(400)),
                    x,ii,-20.,20,pen,100,year_flags(1),-10.0);
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
                  set_value(recr(max(2,first_unfixed_year),
                    last_real_year),x,ii,-20.,20.,1000.);
                }
                else
                {
                  recr(last_real_year-parest_flags(400)+1,last_real_year).
                    initialize();
                  set_value(recr(max(2,first_unfixed_year),
                    last_real_year-parest_flags(400)),x,ii,-20.,20.,1000.);
                }  //NMD_16May2016
                for (int is=2;is<=ns;is++)
                {
  //                set_value(pmsd->recr(is)(1,last_real_year),x,ii,
  //                  -20.,20.,1000.);
                  if (parest_flags(400)==0)   //NMD_19May2016
                  {
                    set_value(pmsd->recr(is)(max(2,first_unfixed_year),
                      last_real_year),x,ii,-20.,20.,1000.);
                  }
                  else
                  {
                    pmsd->recr(is)(last_real_year-parest_flags(400)+1,last_real_year).
                      initialize();
                    set_value(pmsd->recr(is)(max(2,first_unfixed_year),
                      last_real_year-parest_flags(400)),x,ii,-20.,20.,1000.);
                  }   //NMD_19May2016
                }
              }
              else
              {
                dvar_matrix M(1,ns);
  //              M(1)=recr(1,last_real_year);
                M(1)=recr(max(2,first_unfixed_year),last_real_year);  //NMD_19May2016
                if (parest_flags(400)!=0)
                  M(1)(last_real_year-parest_flags(400)+1,last_real_year).
                  initialize();  //NMD_19May2016
                const MY_REAL_DOUBLE th=1000.;
                for (int is=2;is<=ns;is++)
                {
  //                M(is)=pmsd->recr(is)(1,last_real_year);
                  M(is)=pmsd->recr(is)(max(2,first_unfixed_year),last_real_year);  //NMD_19May2016
                  if (parest_flags(400)!=0)
                    M(is)(last_real_year-parest_flags(400)+1,last_real_year).
                    initialize();  //NMD_19May2016
                }
  //              set_value(M,x,ii,-20.,20.,th, 1,column(pmsd->species_flags,1));
                if (parest_flags(400)==0)  //NMD_19May2016
                {
                  set_value(M,x,ii,-20.,20.,th, 1,column(pmsd->species_flags,1));
                }
                else
                {
                  dvar_matrix MM(1,ns);
                  for (int is=1;is<=ns;is++)
                  {
                    MM(is)=M(is)(2,last_real_year-parest_flags(400));
                  }
                  set_value(MM,x,ii,-20.,20.,th, 1,column(pmsd->species_flags,1));
                }  //NMD_19May2016
              }
            }
            else
            {
              if (parest_flags(400)==0)
  //              set_value(recr(1,last_real_year),x,ii,-20.,20.,1000.);
              set_value(recr(max(2,first_unfixed_year),last_real_year),
                x,ii,-20.,20.,1000.);  //NMD_19May2016
              else
              {
                recr(last_real_year-parest_flags(400)+1,last_real_year).
                  initialize();
  //              set_value(recr(1,last_real_year-parest_flags(400)),
  //                x,ii,-20.,20.,1000.);
                set_value(recr(max(2,first_unfixed_year),
                  last_real_year-parest_flags(400)),
                  x,ii,-20.,20.,1000.);   //NMD_19May2016
              }
            }
          }
        }
        else
        {
          MY_DOUBLE_TYPE bd=age_flags(51);
          MY_DOUBLE_TYPE bd1=age_flags(46)/10.0;
  //        set_value(recr,x,ii,-bd,bd1,pen);
          if (parest_flags(400)==0)  //NMD_19May2016
          {
            mfcl_error();
            set_value(recr,x,ii,-bd,bd1,pen);
          }
          else
          {
            recr(last_real_year-parest_flags(400)+1,last_real_year).
              initialize();
            set_value(recr(max(2,first_unfixed_year),
              last_real_year-parest_flags(400)),x,ii,-bd,bd1,pen);
          }  //NMD_19May2016
          if (pmsd && pmsd->num_species>1)
          {
            int ns=pmsd->num_species;
            for (int is=2;is<=ns;is++)
            {
  //            set_value(pmsd->recr(is),x,ii,-bd,bd1,pen);
              if (parest_flags(400)==0)  //NMD_19May2016
              {
                mfcl_error();
                set_value(pmsd->recr(is),x,ii,-bd,bd1,pen);
              }
              else
              {
                pmsd->recr(is)(last_real_year-parest_flags(400)+1,last_real_year).
                  initialize();
                set_value(pmsd->recr(is)(max(2,first_unfixed_year),
                  last_real_year-parest_flags(400)),x,ii,-bd,bd1,pen);
              }  //NMD_19May2016
            }
          }
        }
      }
      else 
      {
        //recinpop_orth_reset(x,ii);
      }
    }
    // relative initial population
     //xofs1  << ii << " ";
    //  cout   << "AE " << ii << " ";
    if (age_flags(31)>0 && age_flags(94)==0)
    {
      //cout << "initpop ii = " << ii << endl;
      if (age_flags(51)==0)
      {
        set_value(initpop,x,ii,-20.,20.,pen,100.0);
      }
      else
      {
        MY_DOUBLE_TYPE bd=age_flags(51);
        MY_DOUBLE_TYPE bd1=age_flags(46)/10.0;
        set_value(initpop,x,ii,-bd,bd1,pen);
      }
    }

     //xofs1  << ii << " ";
    //  cout   << "AF " << ii << " ";
    // overall population scaling
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
          set_value(totpop_coff,x,ii,-5.0,5.0,pen);
          if (pmsd)
            set_value(pmsd->totpop_coff,x,ii,-5.0,5.0,pen);
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
  //          set_value(totpop_coff,x,ii,1,25,pen,100.0);
            set_value(totpop_coff,x,ii,8.0L,25,pen,nscal);
  //NMD_4sep2020
          }
          else
          {
            ivector sf1=column(pmsd->species_flags,1);
            dvar_vector tmp(1,pmsd->num_species);
            ivector v(1,pmsd->num_species);
            v=1;
            //set_value(pmsd->totpop_coff,x,ii,1,25,pen,500);
  //          set_value(tmp,x,ii,1,25,pen,v,sf1,100.0);
            set_value(tmp,x,ii,8.0L,25,pen,v,sf1,nscal);  //NMD_4sep2020
            totpop_coff=tmp(1);
            pmsd->totpop_coff=tmp(2,pmsd->num_species);
          }
        }
        //totpop=x(ii++)/double(pf);
      }
    }
    if (age_flags(113)>0)
    {
      MY_DOUBLE_TYPE oneh=100.;
      set_value(rec_init_diff,x,ii,-12.0,12.0,pen,oneh);
      if (pmsd)
        set_value(pmsd->rec_init_diff,x,ii,-12.0,12.0,pen,oneh);
    }
      
     //xofs1  << ii << " ";
    //  cout   << "AG " << ii << " ";
    //  average natural mortality
    if (age_flags(33)>0)
    {
      nat_mort_coff=x(ii++)/1000.;
    }
    if (pmsd && pmsd->num_species>1)
    {
      for (int is=2;is<=pmsd->num_species;is++)
      {
        if (pmsd->age_flags(is,33)>0)
        {
          pmsd->nat_mort_coff(is)=x(ii++)/1000.;
        }
      }
    }
    if (age_flags(54)>1)
    {
      set_value_partial(avail_coff,x,ii,age_flags(54)-1,-15.,.0,pen);
    }

     //xofs1  << ii << " ";
    //  cout   << "AH " << ii << " ";
    if (sum(column(fish_flags,47)))
    {           
      set_value_partial(seasonal_catchability_pars,x,ii,column(fish_flags,47),
        -4.,4.,pen,column(fish_flags,47) ,column(fish_flags,28),100.0);
    }  

    if (!age_flags(92)) //NMD_29may2020
    {
      set_value(fish_pars(1),x,ii,-4.,4.,pen,
        column(fish_flags,27),column(fish_flags,28),1000.0);

      set_value(fish_pars(2),x,ii,-3.,3.,pen,
        column(fish_flags,27),column(fish_flags,28),10000.0);
    }
   
    set_value(fish_pars(7),x,ii,-2.,2.,pen,
      column(fish_flags,51),column(fish_flags,52),100.0);

    set_value(fish_pars(8),x,ii,-3.,3.,pen,
      column(fish_flags,53),column(fish_flags,54));

    set_value(fish_pars(12),x,ii,-5.,5.,pen,
      column(fish_flags,58),column(fish_flags,59));

    set_value(fish_pars(13),x,ii,-5.,5.,pen,
      column(fish_flags,63),column(fish_flags,64));

    // logvN length for self scaling multinomial
    //set_value(fish_pars(14),x,ii,0.0,9.0,pen,
    //  column(fish_flags,67),column(fish_flags,68));

    // *********************************************
    // logvN length for self scaling multinomial
    {
      MY_DOUBLE_TYPE lbd=1.0;
      MY_DOUBLE_TYPE ubd=7.5;
      if (parest_flags(334)>0)
        ubd=parest_flags(334);
    
      set_value(fish_pars(14),x,ii,lbd,ubd,pen,
        column(fish_flags,67),column(fish_flags,68),500);

    }
    {
      MY_DOUBLE_TYPE lbd=-7.0;
      MY_DOUBLE_TYPE ubd=7.0;
      // variance multiplier coff for length dirichlet multinomial
      set_value(fish_pars(22),x,ii,lbd,ubd,pen,
        column(fish_flags,69),column(fish_flags,68),100.);
      set_value(fish_pars(23),x,ii,lbd,ubd,pen,
        column(fish_flags,89),column(fish_flags,68),100.);
    }
    {
      MY_DOUBLE_TYPE lbd=-7.0;
      MY_DOUBLE_TYPE ubd=7.0;
      // variance multiplier coff for weight dirichlet multinomial
      set_value(fish_pars(24),x,ii,lbd,ubd,pen,
        column(fish_flags,87),column(fish_flags,77),100.);
      set_value(fish_pars(25),x,ii,lbd,ubd,pen,
        column(fish_flags,88),column(fish_flags,77),100.);
    }

    {
      MY_DOUBLE_TYPE lbd=-0.5;
      MY_DOUBLE_TYPE ubd=0.5;
      // variance multiplier coff for RE length hetergeneity
      set_value(fish_pars(26),x,ii,lbd,ubd,pen,
        column(fish_flags,90),column(fish_flags,68),100.); //NMD_5feb2025
//        column(fish_flags,90),column(fish_flags,77),100.);
    }

    {
      MY_DOUBLE_TYPE lbd=-0.5;
      MY_DOUBLE_TYPE ubd=0.5;
      // variance multiplier coff for RE weight hetergeneity
      set_value(fish_pars(27),x,ii,lbd,ubd,pen,
        column(fish_flags,91),column(fish_flags,77),100.);
    }

    // logvN weight for self scaling multinomial
    //set_value(fish_pars(15),x,ii,0.0,9.0,pen,
    //  column(fish_flags,76),column(fish_flags,77));

    // *********************************************
    // logvN weight for self scaling multinomial
    {
      MY_DOUBLE_TYPE lbd=1.0;
      MY_DOUBLE_TYPE ubd=7.5;
      if (parest_flags(335)>0)
        ubd=parest_flags(335);
    
      set_value(fish_pars(15),x,ii,lbd,ubd,pen,
        column(fish_flags,76),column(fish_flags,77),500);
    }

    // *********************************************
    // *********************************************

    // rho length for self scaling multinomial
    {
      MY_DOUBLE_TYPE bd=0.85;
      if (parest_flags(315)>0)
        bd=parest_flags(315)/100.;
    
      set_value(fish_pars(16),x,ii,0.05,bd,pen,
        column(fish_flags,78),column(fish_flags,68),100.0);
    }

    // rho weight for self scaling multinomial
    {
      MY_DOUBLE_TYPE bd=0.85;
      if (parest_flags(370)>0)
        bd=parest_flags(370)/100.;
    
      set_value(fish_pars(17),x,ii,0.05,bd,pen,
        column(fish_flags,80),column(fish_flags,77),100.0);
    }

    // log var length for self scaling multinomial
    {
      MY_DOUBLE_TYPE bd=4.0;
      MY_DOUBLE_TYPE ubd=-0.2;
      if (parest_flags(318)!=0) 
       ubd=parest_flags(318)/100.;
      //if (parest_flags(316)>0)
      //  bd=parest_flags(316)/100.;
      cout << "AA ubd " << ubd << endl;
      set_value(fish_pars(18),x,ii,-bd,ubd,pen,
        column(fish_flags,82),column(fish_flags,68),100.);
    }

    // log var weight for self scaling multinomial
    {
      MY_DOUBLE_TYPE bd=4.0;
      MY_DOUBLE_TYPE ubd=-0.2;
      if (parest_flags(318)!=0) 
       ubd=parest_flags(318)/100.;
      //if (parest_flags(316)>0)
      //  bd=parest_flags(316)/100.;
      cout << "BB ubd " << ubd << endl;
      set_value(fish_pars(19),x,ii,-bd,ubd,pen,
        column(fish_flags,84),column(fish_flags,77),100.);
    }

    // length sample size covariate for self scaling multinomial
    {
      MY_DOUBLE_TYPE bd=1.0;
      if (parest_flags(317)>0)
        bd=parest_flags(317)/100.;
      set_value(fish_pars(20),x,ii,-0.05,bd,pen,
        column(fish_flags,85),column(fish_flags,68),100.);
    }

    // weight sample size covariate for self scaling multinomial
    {
      MY_DOUBLE_TYPE bd=1.0;
      if (parest_flags(317)>0)
        bd=parest_flags(317)/100.;
      set_value(fish_pars(21),x,ii,-0.05,bd,pen,
        column(fish_flags,86),column(fish_flags,77),100.);
    }


     //xofs1  << ii << " ";
    //  cout   << "AJ " << ii << " ";
    if (parest_flags(33)>0)
    {
      MY_DOUBLE_TYPE ub=parest_flags(33)/100.;  
      if (!age_flags(198))
      {
        set_value_exp(fish_pars(3),x,ii,.001,ub,pen,
          column(fish_flags,33),column(fish_flags,34),1000.);
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
        set_value_exp(tag_fish_rep,x,ii,.001,ub,pen,
          tag_fish_rep_active_flags,tag_fish_rep_group_flags,scale);
      }
       //xofs1  << ii << " ";
      //    cout   << "AJ1 " << ii << " ";
    }
    else
    {          
      if (!age_flags(198))
      {
        set_value(fish_pars(3),x,ii,.001,1.01,pen,
          column(fish_flags,33),column(fish_flags,34));
      }
      else
      {
        set_value(tag_fish_rep,x,ii,.001,1.01,pen,
          tag_fish_rep_active_flags,tag_fish_rep_group_flags);
      }
       //xofs1  << ii << " ";
      //    cout   << "AJ2 " << ii << " ";
    }  

     //xofs1  << ii << " ";
    //  cout   << "AK " << ii << " ";
    //  int iiold=ii;
    if (parest_flags(111) == 4)    //NMD_24Sep2019
    {
      if (parest_flags(305)==0)
      {
        set_value(fish_pars(4),x,ii,-50.,50.,pen,
          column(fish_flags,43),column(fish_flags,44),10.0);
        //set_value(fish_pars(28),x,ii,0.001,0.3,pen,
        //  column(fish_flags,43),column(fish_flags,44),100.0);
      }
      else if(parest_flags(305)>0)
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
        set_value(fish_pars(4),x,ii,lb,ub,pen,
        column(fish_flags,43),column(fish_flags,44),10.0);
        //set_value(fish_pars(28),x,ii,0.001,0.3,pen,
        //  column(fish_flags,43),column(fish_flags,44),100.0);
      }
    }
    // !!! XXXX  ****************************************************
    // ****************************************************
    // ****************************************************
    // put in preliminary fish_pars(30) for gamma tags
    // only do it for pf111 from 5 to 8
    // use ff43 to make it active
    // use same grouping flags as negbin ?
    if (parest_flags(111) >= 5 && parest_flags(111) <=8)
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
      set_value(fish_pars(30),x,ii,lb,ub,pen,
      column(fish_flags,43),column(fish_flags,44),10.0);

    }
    /*
    set_value(fish_pars(32),x,ii,0.8,5.0,pen,column(fish_flags,93),
      column(fish_flags,29),10.0);

    set_value(fish_pars(31),x,ii,-2.0,2.0,pen,column(fish_flags,81),
      column(fish_flags,29),100.0);
    */

  //  set_value(fish_pars(33),x,ii,-10.0,2.0,pen,column(fish_flags,93),
  //    column(fish_flags,29),10.0);

    // ****************************************************
    // ****************************************************
    // ****************************************************
     //xofs1  << ii << " ";
    //  cout   << "AK1 " << ii << " ";
    set_value(fish_pars(5),x,ii,.001,.999,pen,
      column(fish_flags,46),column(fish_flags,44));
     //xofs1  << ii << " ";
    //  cout   << "AK2 " << ii << " ";

    set_value(fish_pars(6),x,ii,.001,10.,pen,
      column(fish_flags,46),column(fish_flags,44));

     //xofs1  << ii << " ";
    //  cout   << "AL " << ii << " ";

    if (parest_flags(155)==0)
    {
      MY_DOUBLE_TYPE nscal=1.0;
      if (parest_flags(387))  //NMD_4sep2020
      {
        nscal=3000.0;
      }
      else
      {
        nscal=1000.0;
      }  //NMD_4sep2020
      set_value(region_pars(1),x,ii,1.e-6,1.5,pen,region_flags(1),
        region_flags(3),nscal);
  //      region_flags(3),3000.0);
    }
     //xofs1  << ii << " ";
    // VVVVVVVv
    if (age_flags(125) && age_flags(126))
    {
      int nya=age_flags(95);
      if (nya<1 || nya>95)
      {
        cerr << "values of af95 must be between 1 and 95 when af125"
                " not equal 0" << endl;
        ad_exit(1);
      }
      set_value(region_pars.sub(2,nya),x,ii,-5.0,5.0,pen);
    }
     //xofs1  << ii << " ";
    //  cout   << "AM " << ii << " ";
    if (parest_flags(393))
    {
      set_value_partial(kludged_equilib_coffs,
         x, ii, parest_flags(374), -3.0,0.3,pen,100.0);

      //double one1=1.0;
      MY_DOUBLE_TYPE ten10=10.0;
      MY_DOUBLE_TYPE hun100=100.0;
      MY_DOUBLE_TYPE thou1000=1000.0;
      set_value(kludged_equilib_level_coffs,x,ii,
        -3.0,0.3,pen,thou1000);
    }

    /*
    if (parest_flags(377))
    {
      ivector ff29=column(fish_flags,29);
      ivector num_ifmlrp=implicit_fml_bounds(3);
//      set_value_partial(implicit_fm_level_regression_pars,x,ii,
//        num_ifmlrp,-250.0,250.0,pen,ff29,10000.0);
      set_value_partial(implicit_fm_level_regression_pars,x,ii,
        num_ifmlrp,-500.0,500.0,pen,ff29,10000.0);
    }
    */

    MY_DOUBLE_TYPE one1=1.0;
    //set_value(kludged_equilib_level_coffs,x,ii,-3.0,0.3,pen,one1);
  } // if (parest_flags(392)==0)
  
  set_value(fish_pars(32),x,ii,0.8,5.0,pen,column(fish_flags,93),
    column(fish_flags,29),10.0);

  set_value(fish_pars(31),x,ii,-2.0,2.0,pen,column(fish_flags,81),
    column(fish_flags,29),100.0);

  if (parest_flags(377))
  {
    ivector ff29=column(fish_flags,29);
    ivector num_ifmlrp=implicit_fml_bounds(3);
    set_value_partial(implicit_fm_level_regression_pars,x,ii,
      num_ifmlrp,-500.0,500.0,pen,ff29,10000.0);
  }


  return pen;
}

void dvar_fish_stock_history::xinit(dvector& x,int& ii,d3_array& len_sample_size,
  ofstream& xof)
{
  //ii=1;
  int iisave;
  //if (parest_flags(392)==0)
  //{
  /*
  if (parest_flags(361)>0)
  {
    iisave=ii;
    set_value_inv(testpar,x,ii);
    xinit_message(xof,iisave,ii,"testpar");
  }
  */

  if (parest_flags(392)==0)
  {

    if (age_flags(68)>0)
    {
      switch(age_flags(184))
      {
      case 0:
        {
          MY_DOUBLE_TYPE nscal=1.0;
          /*
          if (!parest_flags(387))
          {
            nscal=300.0;
          }
          else
          {
            nscal=1000.0;
          }
          */
          iisave=ii;
          nscal=1000.0;
          set_value_inv(diff_coffs,x,ii,0.0,3.0,nscal,age_flags(68),move_map);
//          xinit_message(xof,iisave,ii,"diff_coffs");
          xinit_message(xof,iisave,ii,"diff_coffs",diff_coffs);
        }  
        break;
      case 1:
        iisave=ii;
        set_value_inv(xdiff_coffs,x,ii,-200.0,200.0,1.e+3,age_flags(68),move_map);
        xinit_message(xof,iisave,ii,"xdiff_coffs");
        break;
      case 2:
        iisave=ii;
        set_value_inv(zdiff_coffs,x,ii,-10.0,10.0,1.e+2,age_flags(68),move_map);
        xinit_message(xof,iisave,ii,"zdiff_coffs");
        break;
      default:
        cerr << "Illegal value for age_flags(184)" << endl;
        ad_exit(1);
      }
    }

    if (age_flags(88)>0)
    {
      iisave=ii;
      set_value_inv(diff_coffs2,x,ii,-4.0,4.0,1000.,age_flags(88),move_map);
      xinit_message(xof,iisave,ii,"diff_coffs2");
    }

    if (age_flags(90)>0)
    {
      iisave=ii;
      set_value_inv(diff_coffs3,x,ii,-4.0,4.0,1000.,age_flags(90),move_map);
      xinit_message(xof,iisave,ii,"diff_coffs3");
    }

//  if (parest_flags(392)==0)
//  {
    if (age_flags(70)>0)
    {
      if (parest_flags(155)==0)
      {
        iisave=ii;
        if (num_regions>1)
        {
          ivector iv(1,num_regions);
          iv=1;
          const dvar_matrix & tmp
           =trans(region_rec_diff_coffs.sub(first_unfixed_year,last_real_year-1));
        //set_value_inv(region_rec_diff_coffs.sub(1,last_real_year-1),x,ii,-3.,3.,
          set_value_inv(tmp,x,ii,-3.,3.,REGION_REC_SCALE,1,region_flags(3));
//          xinit_message(xof,iisave,ii,"region_rec_diffs");
          xinit_message(xof,iisave,ii,"region_rec_diffs",tmp); //NMD_18Oct2024
        }
      }
    }

    int i;
    if (age_flags(34)>0 && !age_flags(143) && !age_flags(92))
    {
      MY_DOUBLE_TYPE max_dev;
      if (age_flags(35)==0)
      {
        max_dev=1;
      }
      else
      {
        max_dev=age_flags(35);
      }
      if (!sum_ff79_flag)
      {
        for (i=1;i<=num_fisheries;i++)
        {
          if (fish_flags(i,4)>1)
          {
            iisave=ii;
            dvar_vector edc=effort_dev_coffs(i)(first_unfixed_fish_time(i),
              num_real_fish_times(i));
            set_value_inv(edc,zero_effdev_flag(i)(first_unfixed_fish_time(i),
              num_real_fish_times(i)),
              x,ii,-max_dev,max_dev,100.);
  //          xinit_message(xof,iisave,ii,"effort_dev_coffs");
            xinit_message(xof,iisave,ii,"effort_dev_coffs",edc,i); //NMD_5Aug2021
          }
        }
      }
      else
      {
        ivector ff79=column(fish_flags,79);
        ivector ff4=column(fish_flags,4);
        iisave=ii;
        int projflag=0;
        for (i=1;i<=num_fisheries;i++)
        {
          if (num_real_fish_times(i)!=num_fish_times(i))
          {
            projflag=1;
            break;
          }
        }
        if (projflag)
        {
          dvar_matrix sub_effort_dev_coffs(1,num_fisheries);
          for (i=1;i<=num_fisheries;i++)
          {
            sub_effort_dev_coffs(i)=effort_dev_coffs(i)(1,num_real_fish_times(i));
          }
          set_value_inv(sub_effort_dev_coffs,x,ii,-max_dev,max_dev,100.,
             ff4,zero_effdev_flag,ff79);
        }
        else
        {
          int ii1=ii;
          set_value_inv(effort_dev_coffs,x,ii1,-max_dev,max_dev,100.,
             ff4,zero_effdev_flag,ff79);

          set_value_inv(effort_dev_coffs,x,ii,-max_dev,max_dev,100.,
             pgroup_manager_1);

          //set_value_inv(effort_dev_coffs,x,ii,-max_dev,max_dev,100.,
            // ff4,zero_effdev_flag,ff79,pgroup_manager_1);

        }
  //      xinit_message(xof,iisave,ii,"effort_dev_coffs");
        xinit_message(xof,iisave,ii,"effort_dev_coffs",effort_dev_coffs);  //NMD_5Aug2021
      }
    }
    if (!age_flags(92))
    {
      iisave=ii;
      set_value_inv(q0,x,ii,-25.,1.,column(fish_flags,1),
        column(fish_flags,60),10000.0);
      xinit_message(xof,iisave,ii,"   q0");
    }
    else if (missing_catch_flag)
    {
      iisave=ii;
      set_value_inv(q0_miss,x,ii,-25.,1.,column(fish_flags,1),
        column(fish_flags,60),10000.0);
      xinit_message(xof,iisave,ii,"   q0_missing");
      if (age_flags(104)>0)
      {
        iisave=ii;
        set_value_inv(fm_level_devs,x,ii,-2.,2.,100.0);
        xinit_message(xof,iisave,ii,"   fm_level_devs");
      }
    }

    if (parest_flags(290)>0)
    {
      iisave=ii;
      set_value_inv(log_length_variance,x,ii,-2.,5.,100.0);
      xinit_message(xof,iisave,ii,"   log_length_variance");
    }

    if (parest_flags(310)>0)
    {
      iisave=ii;
      set_value_inv(length_tot_exp,x,ii,-5.0,5.0,100.);
      xinit_message(xof,iisave,ii,"   length__tot_exp");
    }


    if (parest_flags(291)>0)
    {
      iisave=ii;
      switch (parest_flags(297))
      {
      case 0:
      case 1:  //LN3m
        set_value_inv(length_rho,x,ii,-0.85,0.85,100.0);
        break;
      case 2:  //LN3
        set_value_inv(length_rho,x,ii,-1.98,1.98,100.0);
        break;
      default:
        cerr << "error in pf297" << endl;
        ad_exit(1);
      }
      xinit_message(xof,iisave,ii,"   length_rho");
    }

    if (parest_flags(292)>0)
    {
      iisave=ii;
      set_value_inv(log_length_dof,x,ii,0.0,10,100.0);
      xinit_message(xof,iisave,ii,"   log_length_dof");
    }

    if (parest_flags(296)>0)
    {
      iisave=ii;
      switch (parest_flags(297))
      {
      case 0:
        set_value_inv(length_psi,x,ii,-0.0001,0.99,100.0);
        break;
      case 1: // LN3m
        set_value_inv(length_psi,x,ii,100.0);
        break;
      case 2:
        set_value_inv(length_psi,x,ii,0.01,0.99,100.0);
        break;
      default:
        cerr << "error in pf297" << endl;
        ad_exit(1);
      }
      xinit_message(xof,iisave,ii,"   length_psi");
    }

    if (parest_flags(298)>0)
    {
      iisave=ii;
      set_value_inv(length_exp,x,ii,-1.98,1.98,100.0);
      xinit_message(xof,iisave,ii,"   length_exp");
    }

    // TTTTTTTTTTTTTTTTTT DIST
    if (parest_flags(280)>0)
    {
      iisave=ii;
      set_value_inv(log_weight_variance,x,ii,-10.,15.,100.0);
      xinit_message(xof,iisave,ii,"   log_weight_variance");
    }

    if (parest_flags(300)>0)
    {
      iisave=ii;
      set_value_inv(weight_tot_exp,x,ii,-5.0,5.0,100.);
      xinit_message(xof,iisave,ii,"   weight_tot_exp");
    }
    if (age_flags(134)>0)
    {
      iisave=ii;
      set_value_inv(species_pars(1),x,ii,-0.975,0.975,10.);
      xinit_message(xof,iisave,ii,"   species pars(1)");
    }
    if (age_flags(136)>0)
    {
      iisave=ii;
      set_value_inv(species_pars(2),x,ii,-0.975,0.975,10.);
      xinit_message(xof,iisave,ii,"   species pars(2)");
    }

    if (parest_flags(281)>0)
    {
      iisave=ii;
      switch (parest_flags(287))
      {
      case 0:
      case 1:  //LN3m
        set_value_inv(weight_rho,x,ii,-0.90,0.90,100.0);
        break;
      case 2:  //LN3
        set_value_inv(weight_rho,x,ii,-1.98,1.98,100.0);
        break;
      default:
        cerr << "error in pf287" << endl;
        ad_exit(1);
      }
      xinit_message(xof,iisave,ii,"   weight_rho");
    }

    if (parest_flags(282)>0)
    {
      iisave=ii;
      set_value_inv(log_weight_dof,x,ii,0.0,10,100.0);
      xinit_message(xof,iisave,ii,"   log_weight_dof");
    }

    if (parest_flags(286)>0)
    {
      iisave=ii;
      switch (parest_flags(287))
      {
      case 0:
        set_value_inv(weight_psi,x,ii,-0.0001,0.99,100.0);
        break;
      case 1: // LN3m
        set_value_inv(weight_psi,x,ii,100.0);
        break;
      case 2:
        set_value_inv(weight_psi,x,ii,0.01,0.99,100.0);
        break;
      default:
        cerr << "error in pf287" << endl;
        ad_exit(1);
      }
      xinit_message(xof,iisave,ii,"   weight_psi");
    }

    if (parest_flags(288)>0)
    {
      iisave=ii;
      set_value_inv(weight_exp,x,ii,-1.98,1.98,100.0);
      xinit_message(xof,iisave,ii,"   weight_exp");
    }

    for (i=1;i<=num_fisheries;i++)
    {
      // time trend in catchability
      if (fish_flags(i,2)>1)
      {
        xof << ii << "   q1" << endl; 
        x(ii++)=boundpin(q1(i),-.99,.99);
      }
    }
    iisave=ii;
    dvector tmp(1,num_fisheries);
    tmp.initialize();
  //  cout << " UUUU " << ii << endl;
    //set_value_inv_partial(value(sel_dev_coffs),x,ii,
    //  column(fish_flags,19),-3.0,3.0,
    //  column(fish_flags,19),column(fish_flags,51),1000.);
    ivector ff19=column(fish_flags,19);
    if (parest_flags(323)==0)
    {
      if (sum(ff19))
      {
        iisave=ii;
        new_fishing_selectivity_interface_xinit(ii,x,-3.0,3.0,100.);
        xinit_message(xof,iisave,ii,"sel_dev_coffs");
      }
    }


    if (sum_ff48_flag)
    { 
      iisave=ii;
      ivector ff48=column(fish_flags,48);
      ivector ffgroup=column(fish_flags,24);
      MY_DOUBLE_TYPE scale=1000.0;
      set_value_inv(bs_selcoff,x,ii,-20.,7.,scale ,ff48,ffgroup);
      if (!sum(ffgroup))
      {
        xinit_message(xof,iisave,ii,"bs_selcoff",bs_selcoff);
      }
      else
      {
//      new call for xinit_message with grouping capability
        imatrix fshgptr;
        fshgptr=fishery_group_ptr(24);
        xinit_message(xof,iisave,ii,"bs_selcoff_gp:",fshgptr,bs_selcoff);
      }
    }

    for (i=1;i<=num_fisheries;i++)
    {
      // selectivity deviations
      if (fish_flags(i,5)>0)
      {
        MY_DOUBLE_TYPE max_dev;
        if (age_flags(36)==0)
        {
          max_dev=1;
        }
        else
        {
          max_dev=age_flags(36);
        }
        if (fish_flags(i,3)>0)
        { 
          MY_DOUBLE_TYPE delt_scale=10.;
          if (parest_flags(151)>0)
          {
            delt_scale=parest_flags(151);
          }
          if (parest_flags(151)<0)
          {
            delt_scale=-1./parest_flags(151);
          }
          for (int j=2;j<=num_fish_times(i);j++)
          {
            if (len_sample_size(realization_region(i,j),realization_period(i,j),
              realization_incident(i,j))>0)
            {
              iisave=ii;
              set_value_inv(delta2(i,j),x,ii,-max_dev,max_dev,delt_scale);
              xinit_message(xof,iisave,ii,"delta2");
            }
          }
        }
      }
    }
    for (i=1;i<=num_fisheries;i++)
    {

      if (fish_flags(i,6)>0)
      {
        iisave=ii;
        x(ii++)=boundpin(corr_wy(i),-.95,.95);
        xinit_message(xof,iisave,ii,"corr_wy");
      }
      if (fish_flags(i,7)>0)
      {
        iisave=ii;
        x(ii++)=boundpin(corr_by(i),-.95,.95);
        xinit_message(xof,iisave,ii,"corr_by");
      }
      if (fish_flags(i,8)>0)
      {
        iisave=ii;
        x(ii++)=boundpin(corr_wc(i),-.95,.95);
        xinit_message(xof,iisave,ii,"corr_wc");
      }
      if (fish_flags(i,9)>0)
      {
        iisave=ii;
        x(ii++)=boundpin(corr_eff(i),-.95,.95);
        xinit_message(xof,iisave,ii,"corr_eff");
      }
    }

    if (!sum(column(fish_flags,29)))
    {
      for (i=1;i<=num_fisheries;i++)
      {
        if (!age_flags(92))
        {
          if (fish_flags(i,10)>0)
          {
            if (fish_flags(i,23)==0)
            {
              iisave=ii;
            //  set_value_inv(catch_dev_coffs(i),x,ii,-.5,.5,100.);
              set_value_inv(catch_dev_coffs(i),x,ii,-.8,.8,100.);
              xinit_message(xof,iisave,ii,"catch_dev_coffs");
            }
            else
            {
              iisave=ii;
            //  set_value_inv(catch_dev_coffs(i),between_times(i),x,ii,-.5,.5,100.);
              set_value_inv(catch_dev_coffs(i),between_times(i),x,ii,-.8,.8,100.);
              xinit_message(xof,iisave,ii,"catch_dev_coffs");
            }
          }
        }
      }
    }
    else
    {
      for (i=1;i<=ngroups;i++)
      {
        if (!age_flags(92))
        {
          if (fish_flags(gfish_index(i,1),10)>0)
          {
            {
              iisave=ii;   // added by JH & PK 21-07-04
              set_value_inv(grouped_catch_dev_coffs(i),grouped_between_times(i),
                 x,ii,-.8,.8,1000.);
              xinit_message(xof,iisave,ii,"grouped catch_dev_coffs");
            }
          }
        }
      }
    }



    for (i=1;i<=num_fisheries;i++)
    {
      {
        if (fish_flags(i,37)>0)
        {
          if (fish_flags(i,45)==0)
          {
            iisave=ii;
            set_value_inv(rep_dev_coffs(i),x,ii,-.5,.5,1.);
            xinit_message(xof,iisave,ii,"rep_dev_coffs");
          }
          else
          {
            iisave=ii;
            set_value_inv(rep_dev_coffs(i),between_times(i),x,ii,-.5,.5,1.);
            xinit_message(xof,iisave,ii,"rep_dev_coffs");
          }
        }
      }
    }

    // relative recruitment
    ivector sf=season_flags(1);
    if (age_flags(30)>0)
    {
      iisave=ii;
      if (parest_flags(155)==0)
      {
        if (age_flags(51)==0)
        {
          //set_value_inv(recr,x,ii);
          
          if (sum(sf))
          {
  //          set_value_inv(recr(1,last_real_year),x,ii,-20.,20.,100,year_flags(1));
            if (parest_flags(400)==0) //NMD_19May2016
            {
              set_value_inv(recr(max(2,first_unfixed_year),
                last_real_year),x,ii,-20.,20.,100,year_flags(1));
            }
            else
            {
              set_value_inv(recr(max(2,first_unfixed_year),
                last_real_year-parest_flags(400)),
                x,ii,-20.,20.,100,year_flags(1));
            } //NMD_19May2016
            if (pmsd && pmsd->num_species>1)
            {
              int ns=pmsd->num_species;
              for (int is=2;is<=ns;is++)
              {
  //              set_value_inv(pmsd->recr(is)(1,last_real_year),x,ii,-20.,20.,100,
  //                year_flags(1));
                if (parest_flags(400)==0)  //NMD_19May2016
                {
                  set_value_inv(pmsd->recr(is)(max(2,first_unfixed_year),
                    last_real_year),x,ii,-20.,20.,100,year_flags(1));
                }
                else
                {
                  set_value_inv(pmsd->recr(is)(max(2,first_unfixed_year),
                    last_real_year-parest_flags(400)),x,ii,-20.,20.,100,
                    year_flags(1));
                }  //NMD_19May2016
              }
            }
          }
          else
          {
            if (pmsd && pmsd->num_species>1)
            {
              int ns=pmsd->num_species;
              int srf3=sum(region_flags(3));
              if (!srf3)
              {
  //              set_value_inv(recr(1,last_real_year),x,ii,-20.,20.,1000.);
                if (parest_flags(400)==0)  //NMD_19May2016
                {
                  set_value_inv(recr(max(2,first_unfixed_year),
                    last_real_year),x,ii,-20.,20.,1000.);
                }
                else
                {
                  set_value_inv(recr(max(2,first_unfixed_year),
                    last_real_year-parest_flags(400)),
                    x,ii,-20.,20.,1000.);
                }  //NMD_19May2016
                for (int is=2;is<=ns;is++)
                {
  //                set_value_inv(pmsd->recr(is)(1,last_real_year),x,ii,
  //                  -20.,20.,1000.);
                  if (parest_flags(400)==0)  //NMD_19May2016
                  {
                    set_value_inv(pmsd->recr(is)(max(2,first_unfixed_year),
                      last_real_year),x,ii,-20.,20.,1000.);
                  }
                  else
                  {
                    set_value_inv(pmsd->recr(is)(max(2,first_unfixed_year),
                      last_real_year-parest_flags(400)),x,ii,-20.,20.,1000.);
                  }  //NMD_19May2016
                }
              }
              else
              {
                dvar_matrix M(1,ns);
  //              M(1)=recr(1,last_real_year);
                M(1)=recr(max(2,first_unfixed_year),last_real_year);  //NMD_19May2016
                for (int is=2;is<=ns;is++)
                {
  //                M(is)=pmsd->recr(is)(1,last_real_year);
                  M(is)=pmsd->recr(is)(max(2,first_unfixed_year),last_real_year);    //NMD_19May2016
                }
                if (parest_flags(400)==0)  //NMD_19May2016
                {
                  set_value_inv(M,x,ii,-20.,20.,1000.0, 1,
                    column(pmsd->species_flags,1));
                }
                else
                {
                  dvar_matrix MM(1,ns);
                  for (int is=1;is<=ns;is++)
                  {
                    MM(is)=M(is)(max(2,first_unfixed_year),
                      last_real_year-parest_flags(400));
                  }
                  set_value_inv(MM,x,ii,-20.,20.,1000.0, 1,
                    column(pmsd->species_flags,1));
                }  //NMD_19May2016
  //              set_value_inv(M,x,ii,-20.,20.,1000.0, 1,
  //                column(pmsd->species_flags,1));
              }
            }
            else
            {
              if (parest_flags(400)==0)
  //              set_value_inv(recr(1,last_real_year),x,ii,-20.,20.,1000.);
                set_value_inv(recr(max(2,first_unfixed_year),
                  last_real_year),x,ii,-20.,20.,1000.);  //NMD_19May2016
              else
              {
  //              set_value_inv(recr(1,last_real_year-parest_flags(400)),
  //                x,ii,-20.,20.,1000.);
                set_value_inv(recr(max(2,first_unfixed_year),
                  last_real_year-parest_flags(400)),x,ii,-20.,20.,1000.);  //NMD_19May2016
              }
            }
          }
        }
        else
        {
          MY_DOUBLE_TYPE bd=age_flags(51);
          MY_DOUBLE_TYPE bd1=age_flags(46)/10.0;
  //        set_value_inv(recr,x,ii,-bd,bd1);
          if (parest_flags(400)==0)  //NMD_19May2016
          {
            mfcl_error();
            set_value_inv(recr,x,ii,-bd,bd1);
          }
          else
            set_value_inv(recr(max(2,first_unfixed_year),last_real_year-parest_flags(400)),x,ii,-bd,bd1);
          if (pmsd && pmsd->num_species>1)   // XXXYYY
          {
            int ns=pmsd->num_species;
            for (int is=2;is<=ns;is++)
            {
  //          set_value_inv(pmsd->recr(is),x,ii,-bd,bd1);
              if (parest_flags(400)==0)  //NMD_19May2016
              {
                mfcl_error();
                set_value_inv(pmsd->recr(is),x,ii,-bd,bd1);
              }
              else
              {
                set_value_inv(pmsd->recr(is)(max(2,first_unfixed_year),
                  last_real_year-parest_flags(400)),x,ii,-bd,bd1);
              }  //NMD_19May2016
            }
          }
        }
//        xinit_message(xof,iisave,ii,"recr");   //NMD_22feb2018
        if (parest_flags(400)==0)  //NMD_18oct2024
        {
          xinit_message(xof,iisave,ii,"recr",recr);   //NMD_18oct2024
        }
        else
        {
          xinit_message(xof,iisave,ii,"recr",recr(max(2,first_unfixed_year),
                  last_real_year-parest_flags(400)));
        } 
      }
      else 
      {
        //recinpop_orth_xinit(xof,x,ii);
      }
    }
    // relative initial population
    if (age_flags(31)>0 && age_flags(94)==0)
    {
      //cout << " starting initpop ii = " << ii << endl;
      iisave=ii;
      if (age_flags(51)==0)
      {
        set_value_inv(initpop,x,ii,-20.,20.,100.);
      }
      else
      {
        MY_DOUBLE_TYPE bd=age_flags(51);
        MY_DOUBLE_TYPE bd1=age_flags(46)/10.0;
        set_value_inv(initpop,x,ii,-bd,bd1);
      }
      xinit_message(xof,iisave,ii,"initpop");
      //cout << " after initpop ii = " << ii << endl;
    }

    // overall population scaling
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
        //cout << " starting totpop ii = " << ii << endl;
        iisave=ii;
        if (age_flags(177)==0 && age_flags(92) ==0)
        {
          set_value_inv(totpop_coff,x,ii,-5.0,5.0);
          
          if (pmsd)
            set_value_inv(pmsd->totpop_coff,x,ii,-5.0,5.0);
          
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
          }  //NMD_4sep2020
          if (!pmsd || pmsd->num_species==1)
          {
  //          set_value_inv(totpop_coff,x,ii,1,25,100.0);
            set_value_inv(totpop_coff,x,ii,8.0L,25,nscal);
          }
          else
          {
            ivector sf1=column(pmsd->species_flags,1);
            dvar_vector tmp(1,pmsd->num_species);
            ivector v(1,pmsd->num_species);
            v=1;
            tmp(1)=totpop_coff;
            tmp(2,pmsd->num_species)=pmsd->totpop_coff;
#if !defined(NO_MY_DOUBLE_TYPE)
  //          set_value_inv(tmp,x,ii,1.0L,25.0,v,sf1,100.0);
#else
  //          set_value_inv(tmp,x,ii,1.0,25.0,v,sf1,100.0);
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
            set_value_inv(tmp,x,ii,8.0L,25.0,v,sf1,nscal);
#else
            set_value_inv(tmp,x,ii,8.0,25.0,v,sf1,nscal);
#endif
          }
        }
        //x(ii++)=pf*totpop;
        xinit_message(xof,iisave,ii,"totpop");
      }
    }
    if (age_flags(113)>0)
    {
      iisave=ii;
      set_value_inv(rec_init_diff,x,ii,-12.0,12.0,100.);
      if (pmsd)
        set_value_inv(pmsd->rec_init_diff,x,ii,-12.0,12.0,100.);

      xinit_message(xof,iisave,ii,"rec_init_diff");
    }
    //  average natural mortality
    if (age_flags(33)>0)
    {
      iisave=ii;
      x(ii++)=value(nat_mort_coff)*1000;;
      xinit_message(xof,iisave,ii,"nat_mort_coff");
    }
    if (pmsd && pmsd->num_species>1)
    {
      for (int is=2;is<=pmsd->num_species;is++)
      {
        if (pmsd->age_flags(is,33)>0)
        {
          iisave=ii;
          x(ii++)=value(pmsd->nat_mort_coff(is))*1000;;
          xinit_message(xof,iisave,ii,"nat_mort_coff");
        }
      }
    }
    if (age_flags(54)>1)
    {
      iisave=ii;
      set_value_inv_partial(avail_coff,x,ii,age_flags(54)-1,-15.,.0);
      xinit_message(xof,iisave,ii,"avail_coff");
    }

    if (sum(column(fish_flags,47)))
    {           
      iisave=ii;
    set_value_inv_partial(value(seasonal_catchability_pars),x,ii,column(fish_flags,47),-4.,4.,
        column(fish_flags,47) ,column(fish_flags,28),100.0);
      xinit_message(xof,iisave,ii,"seasonal_catchability_pars");
    }  

    if (!age_flags(92))  //NMD_29may2020
    {
      iisave=ii;

      set_value_inv(fish_pars(1),x,ii,-4.,4.,column(fish_flags,27),
        column(fish_flags,28),1000.0);
      xinit_message(xof,iisave,ii,"fish_pars(1)");
      iisave=ii;

      set_value_inv(fish_pars(2),x,ii,-3.,3.,column(fish_flags,27),
        column(fish_flags,28),10000.);
      xinit_message(xof,iisave,ii,"fish_pars(2)");
    }
    iisave=ii;

    set_value_inv(fish_pars(7),x,ii,-3.,3.,column(fish_flags,51),
      column(fish_flags,52),100.);
    xinit_message(xof,iisave,ii,"fish_pars(7)");
    iisave=ii;

    set_value_inv(fish_pars(8),x,ii,-3.,3.,column(fish_flags,53),
      column(fish_flags,54));
    xinit_message(xof,iisave,ii,"fish_pars(8)");
    iisave=ii;

    set_value_inv(fish_pars(12),x,ii,-5.,5.,column(fish_flags,58),
      column(fish_flags,59));
    xinit_message(xof,iisave,ii,"fish_pars(12)");
    iisave=ii;

    set_value_inv(fish_pars(13),x,ii,-5.,5.,column(fish_flags,63),
      column(fish_flags,64));
    xinit_message(xof,iisave,ii,"fish_pars(13)");
    iisave=ii;

    // logvN length for self scaling multinomial
    //set_value_inv(fish_pars(14),x,ii,0.0,9.0,column(fish_flags,67),
    //  column(fish_flags,68));
    //xinit_message(xof,iisave,ii,"fish_pars(14) logvN len");
    //iisave=ii;

    // *************************************************
    // *************************************************
    // logvN length for self scaling multinomial
    {
      MY_DOUBLE_TYPE lbd=1.0;
      MY_DOUBLE_TYPE ubd=7.5;
      if (parest_flags(334)>0)
        ubd=parest_flags(334);
      set_value_inv(fish_pars(14),x,ii,lbd,ubd,column(fish_flags,67),
        column(fish_flags,68),500);
      xinit_message(xof,iisave,ii,"fish_pars(14)_logvN_len");
      iisave=ii;
    }
    {
      MY_DOUBLE_TYPE lbd=-7.0;
      MY_DOUBLE_TYPE ubd=7.0;
      // variance multiplier coff for length dirichlet multinomial
      set_value_inv(fish_pars(22),x,ii,lbd,ubd,column(fish_flags,69),
        column(fish_flags,68),100.);
      xinit_message(xof,iisave,ii,"fish_pars(22)_length_dirichlet_mult_var_mult_coff");
      iisave=ii;
      set_value_inv(fish_pars(23),x,ii,lbd,ubd,column(fish_flags,89),
        column(fish_flags,68),100.);
      xinit_message(xof,iisave,ii,"fish_pars(23)_length_dirichlet_mult_var_mult_covariate");
      iisave=ii;
    }
    {
      MY_DOUBLE_TYPE lbd=-7.0;
      MY_DOUBLE_TYPE ubd=7.0;
      // variance multiplier coff for weight dirichlet multinomial
      set_value_inv(fish_pars(24),x,ii,lbd,ubd,column(fish_flags,87),
        column(fish_flags,77),100.);
      xinit_message(xof,iisave,ii,"fish_pars(24) weight dirichlet mult var mult coff");
      iisave=ii;
      set_value_inv(fish_pars(25),x,ii,lbd,ubd,column(fish_flags,88),
        column(fish_flags,77),100.);
      xinit_message(xof,iisave,ii,"fish_pars(25) weight dirichlet mult var mult covariate");
      iisave=ii;
    }


    {
      MY_DOUBLE_TYPE lbd=-0.5;
      MY_DOUBLE_TYPE ubd=0.5;
//      xinit_message(xof,iisave,ii,"fish_pars(26) len vexp ");
//      iisave=ii;
      set_value_inv(fish_pars(26),x,ii,lbd,ubd,
        column(fish_flags,90),column(fish_flags,68),100.); //NMD_5feb2025
//        column(fish_flags,90),column(fish_flags,77),100.);
      xinit_message(xof,iisave,ii,"fish_pars(26) len vexp ");  //NMD_5feb2025
      iisave=ii;  //NMD_5feb2025
    }

    {
      MY_DOUBLE_TYPE lbd=-0.5;
      MY_DOUBLE_TYPE ubd=0.5;
//      xinit_message(xof,iisave,ii,"fish_pars(27) wght vexp ");
//      iisave=ii;
      set_value_inv(fish_pars(27),x,ii,lbd,ubd,
        column(fish_flags,91),column(fish_flags,77),100.);
      xinit_message(xof,iisave,ii,"fish_pars(27) wght vexp "); //NMD_5feb2025
      iisave=ii;  //NMD_5feb2025
    }

    // logvN weight for self scaling multinomial
    //set_value_inv(fish_pars(15),x,ii,0.0,9.0,column(fish_flags,76),
    //  column(fish_flags,77));
    //xinit_message(xof,iisave,ii,"fish_pars(15) logvN wt ");
    //iisave=ii;

    // logvN weight for self scaling multinomial
    {
    
      MY_DOUBLE_TYPE lbd=1.0;
      MY_DOUBLE_TYPE ubd=7.5;
      if (parest_flags(335)>0)
        ubd=parest_flags(335);
    
      set_value_inv(fish_pars(15),x,ii,lbd,ubd,column(fish_flags,76),
        column(fish_flags,77),500);
      xinit_message(xof,iisave,ii,"fish_pars(15) logvN wt ");
      iisave=ii;
    }

    // *************************************************
    // *************************************************
    // rho length for self scaling multinomial
    {
      MY_DOUBLE_TYPE bd=0.85;
      if (parest_flags(315)>0)
        bd=parest_flags(315)/100.;
    
      set_value_inv(fish_pars(16),x,ii,0.05,bd,column(fish_flags,78),
        column(fish_flags,68),100.0);
      xinit_message(xof,iisave,ii,"fish_pars(16) rho len");
      iisave=ii;
    }

    // rho weight for self scaling multinomial
    {
      MY_DOUBLE_TYPE bd=0.85;
      if (parest_flags(370)>0)
        bd=parest_flags(370)/100.;
    
      set_value_inv(fish_pars(17),x,ii,0.05,bd,column(fish_flags,80),
        column(fish_flags,77),100.0);
      xinit_message(xof,iisave,ii,"fish_pars(17) rho wt");
      iisave=ii;
    }


    // var length for self scaling multinomial
    {
      MY_DOUBLE_TYPE bd=4.0;
      MY_DOUBLE_TYPE ubd=-0.2;
      if (parest_flags(318)!=0) 
       ubd=parest_flags(318)/100.;
      //if (parest_flags(316)>0)
      //  bd=parest_flags(316)/100.;
      cout << "A ubd " << ubd << endl;
      set_value_inv(fish_pars(18),x,ii,-bd,ubd,column(fish_flags,82),
        column(fish_flags,68),100.);
      xinit_message(xof,iisave,ii,"fish_pars(18) var len");
      iisave=ii;
    }

    // var weight for self scaling multinomial
    {
      MY_DOUBLE_TYPE bd=4.0;
      MY_DOUBLE_TYPE ubd=-0.2;
      if (parest_flags(318)!=0) 
       ubd=parest_flags(318)/100.;
      //if (parest_flags(316)>0)
      //  bd=parest_flags(316)/100.;
      cout << "B ubd " << ubd << endl;
      set_value_inv(fish_pars(19),x,ii,-bd,ubd,column(fish_flags,84),
        column(fish_flags,77),100.);
      xinit_message(xof,iisave,ii,"fish_pars(19) var weight");
      iisave=ii;
    }

    // length sample size covariate for self scaling multinomial
    {
      MY_DOUBLE_TYPE bd=1.0;
      if (parest_flags(317)>0)
        bd=parest_flags(317)/100.;
      set_value_inv(fish_pars(20),x,ii,-0.05,bd,
        column(fish_flags,85),column(fish_flags,68),100.);
      xinit_message(xof,iisave,ii,"fish_pars(20)_len_sample_size_covariate");
      iisave=ii;
    }

    // weight sample size covariate for self scaling multinomial
    {
      MY_DOUBLE_TYPE bd=1.0;
      if (parest_flags(317)>0)
        bd=parest_flags(317)/100.;
      set_value_inv(fish_pars(21),x,ii,-0.05,bd,
        column(fish_flags,86),column(fish_flags,77),100.);
      xinit_message(xof,iisave,ii,"fish_pars(21)_weight_sample_size_covariate");
      iisave=ii;
    }

    if (parest_flags(33)>0)
    {
      MY_DOUBLE_TYPE ub=parest_flags(33)/100.;  
      if (!age_flags(198))
      {
        set_value_inv_exp(fish_pars(3),x,ii,.001,ub,
        column(fish_flags,33),column(fish_flags,34),1000.);
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
        set_value_inv_exp(tag_fish_rep,x,ii,.001,ub,
          tag_fish_rep_active_flags,tag_fish_rep_group_flags,scale);
      }
    }
    else
    {
      if (!age_flags(198))
      {
        set_value_inv(fish_pars(3),x,ii,.001,1.01,
          column(fish_flags,33),column(fish_flags,34));
      }
      else
      {
        set_value_inv(tag_fish_rep,x,ii,.001,1.01,
          tag_fish_rep_active_flags,tag_fish_rep_group_flags);
      }
    }  
    if (!age_flags(198))
    {
      xinit_message(xof,iisave,ii,"fish_pars(3)");
    }
    else
    {
//      xinit_message(xof,iisave,ii,"tag_fish_rep");
      xinit_message(xof,iisave,ii,"tag_fish_rep",tag_fish_rep_active_flags,
        tag_fish_rep_group_flags,tag_fish_rep); //NMD_18Oct_2024
    }

    if (parest_flags(111) == 4)   //NMD_24Sep2019
    {
      if (parest_flags(305)==0)
      {
        iisave=ii;
        set_value_inv(fish_pars(4),x,ii,-50.,50.,column(fish_flags,43),
          column(fish_flags,44),10.0);
        xinit_message(xof,iisave,ii,"fish_pars(4)");
        //iisave=ii;
        //set_value_inv(fish_pars(28),x,ii,0.001,0.3,column(fish_flags,43),
        //  column(fish_flags,44),100.0);
        //xinit_message(xof,iisave,ii,"fish_pars(28)");
      }
      else if(parest_flags(305)>0)
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
        iisave=ii;
        set_value_inv(fish_pars(4),x,ii,lb,ub,column(fish_flags,43),
          column(fish_flags,44),10.0);
        xinit_message(xof,iisave,ii,"fish_pars(4)");
        //iisave=ii;
        //set_value_inv(fish_pars(28),x,ii,0.001,0.3,column(fish_flags,43),
        //  column(fish_flags,44),100.0);
        //xinit_message(xof,iisave,ii,"fish_pars(28)");
      }
    }
    // !!! XXXX  ****************************************************
    // ****************************************************
    // ****************************************************
    // put in preliminary fish_pars(30) for gamma tags
    // only do it for pf111 from 5 to 8
    // use ff43 to make it active
    // use same grouping flags as negbin ?
    if (parest_flags(111) >= 5 && parest_flags(111) <=8)
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
      iisave=ii;
      set_value_inv(fish_pars(30),x,ii,lb,ub,column(fish_flags,43),
        column(fish_flags,44),10.0);
      xinit_message(xof,iisave,ii,"fish_pars(30)");  //NMD_1Oct2019

    }

    /*
    iisave=ii;
    set_value_inv(fish_pars(32),x,ii,0.8,5.0,column(fish_flags,93),
      column(fish_flags,29),10.0);
    xinit_message(xof,iisave,ii,"fish_pars(32)");  //NMD_1Oct2019
    iisave=ii;
    set_value_inv(fish_pars(31),x,ii,-2.0,2.0,column(fish_flags,81),
      column(fish_flags,29),100.0);
    xinit_message(xof,iisave,ii,"fish_pars(31)");  //NMD_1Oct2019
    */

    /*
    iisave=ii;
    set_value_inv(fish_pars(33),x,ii,-10.0,2.0,column(fish_flags,93),
      column(fish_flags,29),10.0);
    xinit_message(xof,iisave,ii,"fish_pars(33)");  //NMD_1Oct2019
    */
    
    iisave=ii;
    set_value_inv(fish_pars(5),x,ii,.001,.999,column(fish_flags,46),
      column(fish_flags,44));
    xinit_message(xof,iisave,ii,"fish_pars(5)");

    iisave=ii;
    set_value_inv(fish_pars(6),x,ii,.001,10.,column(fish_flags,46),
      column(fish_flags,44));
    xinit_message(xof,iisave,ii,"fish_pars(6)");

    if (parest_flags(155)==0)
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
      set_value_inv(region_pars(1),x,ii,1.e-6,1.5,region_flags(1),
        region_flags(3),nscal);
  //      region_flags(3),3000.0); //SDH revert back 10/12/2007
  //NMD_4sep2020
      xinit_message(xof,iisave,ii,"region_pars(1)");
    }
    if (age_flags(125) && age_flags(126))
    {
      int nya=age_flags(95);
      if (nya<1 || nya>95)
      {
        cerr << "values of af95 must be between 1 and 95 when af125"
                " not equal 0" << endl;
        ad_exit(1);
      }
      iisave=ii;
      set_value_inv(region_pars.sub(2,nya),x,ii,-5.0,5.0);
      xinit_message(xof,iisave,ii,"region_pars(1)");
    }
    if (parest_flags(393))
    {
      iisave=ii;
      set_value_inv_partial(kludged_equilib_coffs,
         x, ii, parest_flags(374),-3.0,0.3,100.0);
//      xinit_message(xof,iisave,ii,"kludged_equilib_coffs");
      xinit_message(xof,iisave,ii,"kludged_equilib_coffs",parest_flags(374),
                    kludged_equilib_coffs); //NMD_18Oct_2024

      iisave=ii;
      MY_DOUBLE_TYPE thou1000=1000.0;
      set_value_inv(kludged_equilib_level_coffs,x,ii,-3.0,0.3,thou1000);
      xinit_message(xof,iisave,ii,"kludged_equilib_level_coffs");
    }

    /*
    if (parest_flags(377))
    {
      ivector ff29=column(fish_flags,29);
      ivector num_ifmlrp=implicit_fml_bounds(3);
      iisave=ii;
//      set_value_inv_partial(implicit_fm_level_regression_pars,x,ii,
//        num_ifmlrp,-250.0,250.0,ff29,10000.0);
      set_value_inv_partial(implicit_fm_level_regression_pars,x,ii,
        num_ifmlrp,-500.0,500.0,ff29,10000.0);
      xinit_message(xof,iisave,ii,"implicit_fm_level_regression_pars",
        implicit_fm_level_regression_pars,num_ifmlrp);
    }
    */
    
  } // if (parest_flags(392)==0)
  
  iisave=ii;
  set_value_inv(fish_pars(32),x,ii,0.8,5.0,column(fish_flags,93),
    column(fish_flags,29),10.0);
  xinit_message(xof,iisave,ii,"fish_pars(32)");  //NMD_1Oct2019
  iisave=ii;
  set_value_inv(fish_pars(31),x,ii,-2.0,2.0,column(fish_flags,81),
    column(fish_flags,29),100.0);
  xinit_message(xof,iisave,ii,"fish_pars(31)");  //NMD_1Oct2019

  if (parest_flags(377))
  {
    ivector ff29=column(fish_flags,29);
    ivector num_ifmlrp=implicit_fml_bounds(3);
    iisave=ii;
    set_value_inv_partial(implicit_fm_level_regression_pars,x,ii,
      num_ifmlrp,-500.0,500.0,ff29,10000.0);
    xinit_message(xof,iisave,ii,"implicit_fm_level_regression_pars",
      implicit_fm_level_regression_pars,num_ifmlrp);
  }

}
int dvar_fish_stock_history::nvcal(d3_array& len_sample_size)
{
  int ii=0;
  //if (parest_flags(392)==0)
  //{
  /*
  if (parest_flags(361)>0)
  {
    ii+=1;
  }
  */
  if (parest_flags(392)==0)
  {
    if (age_flags(68)>0)
    {
      switch(age_flags(184))
      {
      case 0:
        ii+=num_active(diff_coffs,age_flags(68),move_map);
        break;
      case 1:
        ii+=num_active(xdiff_coffs,age_flags(68),move_map);
        break;
      case 2:
        ii+=num_active(zdiff_coffs,age_flags(68),move_map);
        break;
      default:
        cerr << "Illegal value for age_flags(184)" << endl;
        ad_exit(1);
      }
    }

    if (age_flags(88)>0)
    {
      ii+=num_active(diff_coffs2,age_flags(88),move_map);
    }

    if (age_flags(90)>0)
    {
      ii+=num_active(diff_coffs3,age_flags(90),move_map);
    }

//  if (parest_flags(392)==0)
//  {
    if (age_flags(70)>0)
    {
      if (parest_flags(155)==0)
      {
        if (num_regions>1)
        {
          dvar_matrix tmp(1,num_regions,first_unfixed_year,last_real_year-1);
          ii+=size_count(tmp,1,region_flags(3));
        }
      }
    }

    int i;
    if (age_flags(34)>0 && !age_flags(143) && !age_flags(92))
    {
      MY_DOUBLE_TYPE max_dev;
      if (age_flags(35)==0)
      {
        max_dev=1;
      }
      else
      {
        max_dev=age_flags(35);
      }
      if (!sum_ff79_flag)
      {
        for (i=1;i<=num_fisheries;i++)
        {
          if (fish_flags(i,4)>1)
          {
            ii+=sum(zero_effdev_flag(i)(first_unfixed_fish_time(i),
              num_real_fish_times(i)));
          }
        }
      }
      else
      {
        ivector ff79=column(fish_flags,79);
        ivector ff4=column(fish_flags,4);
        int projflag=0;
        for (i=1;i<=num_fisheries;i++)
        {
          if (num_real_fish_times(i)!=num_fish_times(i))
          {
            projflag=1;
            break;
          }
        }
        if (projflag)
        {
          dvar_matrix sub_effort_dev_coffs(1,num_fisheries);
          for (i=1;i<=num_fisheries;i++)
          {
            sub_effort_dev_coffs(i)=effort_dev_coffs(i)(1,num_real_fish_times(i));
          }
          ii+=num_active(sub_effort_dev_coffs,ff4,zero_effdev_flag,ff79);
        }
        else
        {
          ivector oldff77=column(old_fish_flags,77);
          
          int ii1=num_active(effort_dev_coffs,ff4,zero_effdev_flag,ff79);
          mykludge(ii1);
          
          ii+=num_active(effort_dev_coffs,ff4,zero_effdev_flag,ff79,
            pgroup_manager_1,oldff77);
        }
      }
    }

    //removed this line feb07 05 DF
    if (!age_flags(92))
    {
      ii+=num_active(q0,column(fish_flags,1),column(fish_flags,60));
    }
    else if (missing_catch_flag)
    {
      ii+=num_active(q0_miss,column(fish_flags,1),column(fish_flags,60));
      if (age_flags(104)>0)
      {
        ii+=size_count(fm_level_devs);
      }
    }
    if (parest_flags(290)>0)
      ii+=1;    //   log_length_variance;

    // TTTTTTTTTTTTTTTTTT DIST
    if (parest_flags(310)>0)
      ii+=1;    //  length_tot_exp;


    // TTTTTTTTTTTTTTTTTT DIST
    if (parest_flags(291)>0)
      ii+=1;    //   length_rho;

    // TTTTTTTTTTTTTTTTTT DIST
    if (parest_flags(292)>0)
      ii+=1;    //   length_rho;

    // TTTTTTTTTTTTTTTTTT DIST
    if (parest_flags(296)>0)
      ii+=1;    //   length_psi;

    // TTTTTTTTTTTTTTTTTT DIST
    if (parest_flags(298)>0)
      ii+=1;    //   length_exp;

    // TTTTTTTTTTTTTTTTTT DIST
    if (parest_flags(280)>0)
      ii+=1;    //   log_length_variance;

    // TTTTTTTTTTTTTTTTTT DIST
    if (parest_flags(300)>0)
      ii+=1;    //   weight_tot_exp;;
   if (age_flags(134)>0)
    {
      ii+=size_count(species_pars(1));
    }
    if (age_flags(136)>0)
    {
      ii+=size_count(species_pars(2));
    }

    // TTTTTTTTTTTTTTTTTT DIST
    if (parest_flags(281)>0)
      ii+=1;    //   length_rho;

    // TTTTTTTTTTTTTTTTTT DIST
    if (parest_flags(282)>0)
      ii+=1;    //   length_rho;

    // TTTTTTTTTTTTTTTTTT DIST
    if (parest_flags(286)>0)
      ii+=1;    //   length_psi;

    if (parest_flags(288)>0)
      ii+=1;    //   weight_exp;

    for (i=1;i<=num_fisheries;i++)
    {
      // time trend in catchability
      if (fish_flags(i,2)>1)
      {
        ii++;
      }
    }

    if (parest_flags(323)==0)
    {
      if (sum(column(fish_flags,19)))
      {
        //int itmp=0;
        int itmp1=0;
        //itmp=num_active_partial(value(sel_dev_coffs),column(fish_flags,19),
        //  column(fish_flags,51),column(fish_flags,19));
        itmp1+=new_fishing_selectivity_interface_nvar();
        cout <<  itmp1 << endl;
        ii+=itmp1;
      }
    }


    // Don't think we need this any more with new bs_selcoff stuff
    //ii+=fishing_selectivity_interface_nvar();

    if (sum_ff48_flag)
    { 
      ivector ff48=column(fish_flags,48);
      ivector ffgroup=column(fish_flags,24);
      ii+=size_count(bs_selcoff,ff48,ffgroup);
    }

    for (i=1;i<=num_fisheries;i++)
    {
      // selectivity deviations
      if (fish_flags(i,5)>0)
      {
        if (fish_flags(i,3)>0)
        {
          for (int j=2;j<=num_fish_times(i);j++)
          {
            if (len_sample_size(realization_region(i,j),realization_period(i,j),
              realization_incident(i,j))>0)
            {
              ii+=size_count(delta2(i,j));
            }
          }
        }
      }
    }
    for (i=1;i<=num_fisheries;i++)
    {
      if (fish_flags(i,6)>0)
      {
        ii++;
      }
      if (fish_flags(i,7)>0)
      {
        ii++;
      }
      if (fish_flags(i,8)>0)
      {
        ii++;
      }
      if (fish_flags(i,9)>0)
      {
        ii++;
      }
    }

    if (!sum(column(fish_flags,29)))
    {
      for (i=1;i<=num_fisheries;i++)
      {
        if (!age_flags(92))
        {
          if (fish_flags(i,10)>0)
          {
            if (fish_flags(i,23)==0)
            {
              ii+=size_count(catch_dev_coffs(i));
            }
            else
            {
              int mmin=between_times(i).indexmin();
              int mmax=between_times(i).indexmax();
              for (int j=mmin;j<=mmax;j++)
              {
                if (between_times(i,j)) ii++;
              }
            }
          }
        }
      }
    }
    else
    {
      for (i=1;i<=ngroups;i++)
      {
        if (!age_flags(92))
        {
          if (fish_flags(gfish_index(i,1),10)>0)
          {
            int mmin=grouped_between_times(i).indexmin();
            int mmax=grouped_between_times(i).indexmax();
            for (int j=mmin;j<=mmax;j++)
            {
              if (grouped_between_times(i,j)) ii++;
            }
          }
        }
      }
    }


    for (i=1;i<=num_fisheries;i++)
    {
      {
        if (fish_flags(i,37)>0)
        {
          if (fish_flags(i,45)==0)
          {
            ii+=size_count(rep_dev_coffs(i));
          }
          else
          {
            int mmin=between_times(i).indexmin();
            int mmax=between_times(i).indexmax();
            for (int j=mmin;j<=mmax;j++)
            {
              if (between_times(i,j)) ii++;
            }
          }
        }
      }
    }

    // relative recruitment
    if (age_flags(30)>0)
    {
      if (parest_flags(155)==0)
      {
        if (sum(season_flags(1))==0)
        {
          int srf3=0;
          if (pmsd)
                srf3=sum(pmsd->species_flags(1));
          if (!srf3)
          {
            if (parest_flags(400)==0)
  //            ii+=size_count(recr(1,last_real_year));
              ii+=size_count(recr(max(2,first_unfixed_year),
                last_real_year));  //NMD_19May2016
            else
            {
  //            ii+=size_count(recr(1,last_real_year-parest_flags(400)));
              ii+=size_count(recr(max(2,first_unfixed_year),
               last_real_year-parest_flags(400))); //NMD_19May2016
            }
            if (pmsd && pmsd->num_species>1)
            {
              int ns=pmsd->num_species;
              for (int is=2;is<=ns;is++)
              {
  //              ii+=size_count(pmsd->recr(is)(1,last_real_year));
                if (parest_flags(400)==0)   //NMD_19May2016
                {
                  ii+=size_count(pmsd->recr(is)(max(2,first_unfixed_year),
                    last_real_year));
                }
                else
                {
                  ii+=size_count(pmsd->recr(is)(max(2,first_unfixed_year),
                    last_real_year-parest_flags(400)));
                }   //NMD_19May2016
              }
            }
          }
          else
          {
            int ns=pmsd->num_species;
            dvar_matrix M(1,ns);
  //          M(1)=recr(1,last_real_year);
            M(1)=recr(2,last_real_year);  //NMD_19May2016
            for (int is=2;is<=ns;is++)
            {
  //            M(is)=pmsd->recr(is)(1,last_real_year);
              M(is)=pmsd->recr(is)(max(2,first_unfixed_year),last_real_year);  //NMD_19May2016
            }
  //          ii+=size_count(M,1,column(pmsd->species_flags,1));
            if (parest_flags(400)==0)  //NMD_19May2016
            {
              ii+=size_count(M,1,column(pmsd->species_flags,1));
            }
            else
            {
              dvar_matrix MM(1,ns);
              for (int is=1;is<=ns;is++)
              {
                MM(is)=M(is)(max(2,first_unfixed_year),
                  last_real_year-parest_flags(400));
              }
              ii+=size_count(MM,1,column(pmsd->species_flags,1));
            }  //NMD_19May2016
          }
        }
        else
        {
          ii+=sum(year_flags(1));
          if (pmsd && pmsd->num_species>1)
          {
            int ns=pmsd->num_species;
            for (int is=2;is<=ns;is++)
            {
              ii+=sum(year_flags(1));
            }
          }
        }
      }
      else 
      {
        //ii+=recinpop_orth_size_count();
      }
    }
    // relative initial population
    if (age_flags(31)>0 && age_flags(94)==0)
    {
      ii+=size_count(initpop);
    }

    // overall population scaling
    if (age_flags(32)>0)
    {
      if (parest_flags(155)==0)
      {
        if (!pmsd)
          ii++;
        else
        {
          ivector v(1,pmsd->num_species);
          v=1;
          dvar_vector tmp(1,pmsd->num_species);
        
          //ii+=size_count(pmsd->totpop_coff);
          ii+=
             num_active(tmp,v,column(pmsd->species_flags,1));
        }
      }
    }
    if (age_flags(113)>0)
    {
      ii++;
      if (pmsd)
        ii+=size_count(pmsd->rec_init_diff);
    }
    //  average natural mortality
    if (age_flags(33)>0)
    {
      ii++;
    }
    if (pmsd && pmsd->num_species>1)
    {
      for (int is=2;is<=pmsd->num_species;is++)
      {
        if (pmsd->age_flags(is,33)>0)
        {
          ii++;;
        }
      }
    }
    if (age_flags(54)>1)
    {
      ii+=size_count_partial(avail_coff,age_flags(54)-1);
    }

    if (sum(column(fish_flags,47)))
    {           
      ii+=num_active_partial(seasonal_catchability_pars,column(fish_flags,47),
      column(fish_flags,28),column(fish_flags,47));
    }  

    if (!age_flags(92)) //NMD_29may2020
    {
      ii+=num_active(fish_pars(1),column(fish_flags,27),column(fish_flags,28));
      ii+=num_active(fish_pars(2),column(fish_flags,27),column(fish_flags,28));
    }
    ii+=num_active(fish_pars(7),column(fish_flags,51),column(fish_flags,52));
    ii+=num_active(fish_pars(8),column(fish_flags,53),column(fish_flags,54));
    ii+=num_active(fish_pars(12),column(fish_flags,58),column(fish_flags,59));
    ii+=num_active(fish_pars(13),column(fish_flags,63),column(fish_flags,64));
    ii+=num_active(fish_pars(14),column(fish_flags,67),column(fish_flags,68));
    ii+=num_active(fish_pars(22),column(fish_flags,69),column(fish_flags,68));
    ii+=num_active(fish_pars(23),column(fish_flags,89),column(fish_flags,68));
    ii+=num_active(fish_pars(24),column(fish_flags,87),column(fish_flags,77));
    ii+=num_active(fish_pars(25),column(fish_flags,88),column(fish_flags,77));
    ii+=num_active(fish_pars(26),column(fish_flags,90),column(fish_flags,68));
//    ii+=num_active(fish_pars(26),column(fish_flags,90),column(fish_flags,77));
    ii+=num_active(fish_pars(27),column(fish_flags,91),column(fish_flags,77));
    ii+=num_active(fish_pars(15),column(fish_flags,76),column(fish_flags,77));
    ii+=num_active(fish_pars(16),column(fish_flags,78),column(fish_flags,68));
    ii+=num_active(fish_pars(17),column(fish_flags,80),column(fish_flags,77));
    ii+=num_active(fish_pars(18),column(fish_flags,82),column(fish_flags,68));
    ii+=num_active(fish_pars(19),column(fish_flags,84),column(fish_flags,77));
    ii+=num_active(fish_pars(20),column(fish_flags,85),column(fish_flags,68));
    ii+=num_active(fish_pars(21),column(fish_flags,86),column(fish_flags,77));

    if (!age_flags(198))
    {
      ii+=num_active(fish_pars(3),
        column(fish_flags,33),column(fish_flags,34));
    }
    else
    {
      ii+=num_active(tag_fish_rep,
          tag_fish_rep_active_flags,tag_fish_rep_group_flags);
    }

    //  ii+=num_active(fish_pars(4),column(fish_flags,43),column(fish_flags,44));
    if (parest_flags(111) == 4)   //NMD_24Sep2019
    {
      if (parest_flags(305)>=0)
      ii+=num_active(fish_pars(4),column(fish_flags,43),column(fish_flags,44));
    }
    if (parest_flags(111) >= 5 && parest_flags(111) <=8)
    {
      ii+=num_active(fish_pars(30),column(fish_flags,43),column(fish_flags,44));
    }
//    ii+=num_active(fish_pars(32),column(fish_flags,93),column(fish_flags,29));
//    ii+=num_active(fish_pars(31),column(fish_flags,81),column(fish_flags,29));
  //  ii+=num_active(fish_pars(33),column(fish_flags,93),column(fish_flags,29));
    //ii+=num_active(fish_pars(28),column(fish_flags,43),column(fish_flags,44));
    
    ii+=num_active(fish_pars(5),column(fish_flags,46),column(fish_flags,44));
    ii+=num_active(fish_pars(6),column(fish_flags,46),column(fish_flags,44));
    if (parest_flags(155)==0)
    {
      ii+=num_active(region_pars(1),region_flags(1),region_flags(3));
    }

    if (age_flags(125) && age_flags(126))
    {
      int nya=age_flags(95);
      if (nya<1 || nya>95)
      {
        cerr << "values of af95 must be between 1 and 95 when af125"
                " not equal 0" << endl;
        ad_exit(1);
      }
      ii+=(nya-1)*num_regions;
    }
    if (parest_flags(393))
    {
      ii+=num_active_partial(kludged_equilib_coffs,parest_flags(374));

      ii+=size_count(kludged_equilib_level_coffs);
    }

    /*
    if (parest_flags(377))
    {
      ivector ff29=column(fish_flags,29);
      ivector num_ifmlrp=implicit_fml_bounds(3);
      ii+=num_active_partial(implicit_fm_level_regression_pars,ff29,
            num_ifmlrp);
    }
    */
    
  } //if (parest_flags(392)==0)
  
  ii+=num_active(fish_pars(32),column(fish_flags,93),column(fish_flags,29));
  ii+=num_active(fish_pars(31),column(fish_flags,81),column(fish_flags,29));
  if (parest_flags(377))
  {
    ivector ff29=column(fish_flags,29);
    ivector num_ifmlrp=implicit_fml_bounds(3);
    ii+=num_active_partial(implicit_fm_level_regression_pars,ff29,
          num_ifmlrp);
  }

  return ii;
}

void xinit_message(ofstream& xof,int iisave,int ii,const char * msg)
{
  for (int i=iisave;i<ii;i++)
  {
    xof << i << "  " << msg << endl;
  }
}

void xinit_message(ofstream& xof,int iisave,int ii,const char * msg,
  const dvar_vector& M )
{
  int mmin=M.indexmin();
  int mmax=M.indexmax();
  int i1=iisave;
  for (int i=mmin;i<=mmax;i++)
  {
    xof << i1++ << "  " << msg << "(" << i << ")" << endl;
  }
}

void xinit_message(ofstream& xof,int iisave,int ii,const char * msg,
  const dvar_vector& M, int fi)
{
  int mmin=M.indexmin();
  int mmax=M.indexmax();
  int i1=iisave;
  for (int i=mmin;i<=mmax;i++)
  {
    xof << i1++ << "  " << msg << "(" << fi << "," << i << ")" << endl;
  }
}

void xinit_message(ofstream& xof,int iisave,int ii,const char * msg,
  const dvar_matrix& M )
{
  int mmin=M.indexmin();
  int mmax=M.indexmax();
  int i1=iisave;
  for (int i=mmin;i<=mmax;i++)
  {
    int jmin=M(i).indexmin();
    int jmax=M(i).indexmax();
    for (int j=jmin;j<=jmax;j++)
    {
      xof << i1++ << "  " << msg << "(" << i << "," << j << ")" << endl;
    }
  }
}

void xinit_message(ofstream& xof,int iisave,int ii,const char * msg,
  const dvar_matrix& M ,ivector& jmax)
{
  int mmin=M.indexmin();
  int mmax=M.indexmax();
  int i1=iisave;
  for (int i=mmin;i<=mmax;i++)
  {
    int jmin=M(i).indexmin();
    for (int j=jmin;j<=jmax(i);j++)
    {
      xof << i1++ << "  " << msg << "(" << i << "," << j << ")" << endl;
    }
  }
}

void xinit_message(ofstream& xof,int iisave,int ii,const char * msg,
  const imatrix& mflags,const imatrix& mgroup,const dvar_matrix& M )
{
  int i1=iisave;
  dvar_vector w=rowstack(M);
  ivector flags=rowstack(mflags);
  ivector group=rowstack(mgroup);

  int mmin=w.indexmin();
  int mmax=w.indexmax();
  if (sum(flags))
  {
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      int flag_value=flags(kin);
      MY_DOUBLE_TYPE w_value=value(w(kin));
      if (flags(kin)) xof << i1++ << "  " << msg << "(" << group(kin)
                          << ")" << endl;
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          if (value(w(key(i)))!=w_value)
          {
            cerr << "Error -- grouped initial parameters have unequal values"
                 << endl;
            ad_exit(1);
          }
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
          w_value=value(w(kin));
          if (flags(kin)) xof << i1++ << "  " << msg << "(" << group(kin)
                              << ")" << endl;
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) xof << i1++ << "  " << msg << "(" << i << ")" << endl;
      }
    }
  }
  rowunstack(w,M);
}



void xinit_message(ofstream& xof,int iisave,int ii,const char * msg,
  const dvar3_array& M )
{
  int mmin=M.indexmin();
  int mmax=M.indexmax();
  int i1=iisave;
  for (int i=mmin;i<=mmax;i++)
  {
    int jmin=M(i).indexmin();
    int jmax=M(i).indexmax();
    for (int j=jmin;j<=jmax;j++)
    {
      int kmin=M(i,j).indexmin();
      int kmax=M(i,j).indexmax();
      for (int k=kmin;k<=kmax;k++)
      {
        xof << i1++ << "  " << msg << "(" << i << "," << j 
            << "," << k << ")" << endl;
      }
    }
  }
}


void xinit_message(ofstream& xof,int iisave,int ii,const char * msg,
  const dvar4_array& M )
{
  int mmin=M.indexmin();
  int mmax=M.indexmax();
  int i1=iisave;
  for (int i=mmin;i<=mmax;i++)
  {
    int jmin=M(i).indexmin();
    int jmax=M(i).indexmax();
    for (int j=jmin;j<=jmax;j++)
    {
      int kmin=M(i,j).indexmin();
      int kmax=M(i,j).indexmax();
      for (int k=kmin;k<=kmax;k++)
      {
        int lmin=M(i,j,k).indexmin();
        int lmax=M(i,j,k).indexmax();
        for (int l=lmin;l<=lmax;l++)
        {
          xof << i1++ << "  " << msg << "(" << i << "," << j 
              << "," << k << "," << l << ")" << endl;
        }
      }
    }
  }
}


void xinit_message(ofstream& xof,int iisave,int ii,const char * msg,
  int vmax,const dvar4_array& M )
{
  int mmin=M.indexmin();
  int mmax=M.indexmax();
  int i1=iisave;
  for (int i=mmin;i<=mmax;i++)
  {
    int jmin=M(i).indexmin();
    int jmax=M(i).indexmax();
    for (int j=jmin;j<=jmax;j++)
    {
      int kmin=M(i,j).indexmin();
      int kmax=M(i,j).indexmax();
      for (int k=kmin;k<=kmax;k++)
      {
        int lmin=M(i,j,k).indexmin();
        int lmax=vmax;
        for (int l=lmin;l<=lmax;l++)
        {
          xof << i1++ << "  " << msg << "(" << i << "," << j 
              << "," << k << "," << l << ")" << endl;
        }
      }
    }
  }
}


void xinit_message(ofstream& xof,int iisave,int ii,const char * msg,
  const imatrix& fshgptr, const dvar4_array& M )
{
  int mmin=fshgptr.indexmin();
  int mmax=fshgptr.indexmax();
  int i1=iisave;
  for (int i=mmin;i<=mmax;i++)
  {
    int jmin=M(fshgptr(i,1)).indexmin();
    int jmax=M(fshgptr(i,1)).indexmax();
    adstring suffix=itoa(i,10);
    for (int j=jmin;j<=jmax;j++)
    {
      int kmin=M(fshgptr(i,1),j).indexmin();
      int kmax=M(fshgptr(i,1),j).indexmax();
      for (int k=kmin;k<=kmax;k++)
      {
        int lmin=M(fshgptr(i,1),j,k).indexmin();
        int lmax=M(fshgptr(i,1),j,k).indexmax();
        for (int l=lmin;l<=lmax;l++)
        {
          xof << i1++ << "  " << msg+suffix << "(" << i << "," << j 
              << "," << k << "," << l << ")" << endl;
        }
      }
    }
  }
}


void print_big_effort(dvar_vector& edc,ivector& zcf,ofstream& ofs)
{
  int mmin=edc.indexmin();
  int mmax=edc.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    if (!zcf(i) && fabs(value(edc(i)))>1.e-10)
      ofs << zcf(i) << " " << edc(i) << endl;
  }
}

#undef HOME_VERSION
