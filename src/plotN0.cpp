/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE



#define HOME_VERSION
#include "all.hpp"

#include "variable.hpp"

 void F0biomass_calcs(adstring& filename,dvar_len_fish_stock_history& fsh,
                      dvar_matrix& sel) // JH 11-Nov-03
 {
   ofstream ofs((char*) filename,ios::app);
  //  Calculate population biomass in absence of fishing
  dmatrix bio_F0(1,fsh.nyears,1,fsh.num_regions);
  dmatrix sbio_F0(1,fsh.nyears,1,fsh.num_regions);
  bio_F0.initialize();
  sbio_F0.initialize();
  int nr, na, ny, i;
  for (ny=1;ny<=fsh.nyears;ny++)
  {
    for (nr=1;nr<=fsh.num_regions;nr++)
    {
      dvar_vector pmature = fsh.get_pmature_region(nr);   //NMD_16jul2021     
      for (na=1;na<=fsh.nage;na++)
      {
        bio_F0(ny,nr)+=exp(value(fsh.N_q0(nr,ny,na)))
                       *value(fsh.mean_weight_yr(nr,ny,na))/1000.;
        sbio_F0(ny,nr)+=exp(value(fsh.N_q0(nr,ny,na)))*value(pmature(na))
                        *value(fsh.mean_weight_yr(nr,ny,na))/1000.; //NMD_16jul2021
      }
    }
  }
  ofs << "# Total biomass in absence of fishing" << endl;
  ofs << setprecision(5)  << bio_F0 << endl; 
  ofs << "# Adult biomass in absence of fishing" << endl;
  ofs << setprecision(5)  << sbio_F0 << endl; 
  //--------------------
  // JH 11-Nov-03
  ofs << "# Exploitable populations in absence of fishing " << endl;
  dmatrix exp_pop(1,fsh.nyears,1,fsh.num_fisheries);
  exp_pop.initialize();
  int iy, j;
  for (iy=1;iy<=fsh.nyears;iy++)
  {
    for (int ifish=1;ifish<=fsh.num_fisheries;ifish++)
    {
      int rr=fsh.realization_region(ifish,1);
      for (j=1;j<=fsh.nage;j++)
      {
        exp_pop(iy,ifish)+=value(sel(ifish,j))*exp(value(fsh.N_q0(rr,iy,j)))*
                           value(fsh.mean_weight_yr(rr,iy,j))/1000.;
      }
    }
    ofs <<  setscientific() << setprecision(3) << exp_pop(iy) << endl;
  }
  //--------------------

  ofs << "# Predicted catch for interaction analysis by fishery (down) and time (across)" << endl;
  for (i=1;i<=fsh.num_fisheries;i++)
  {
    for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
    {
      int rr=fsh.realization_region(i,nt);
      int rp=fsh.realization_period(i,nt);
      int ri=fsh.realization_incident(i,nt);
      if (!fsh.data_fish_flags(1,fsh.parent(rr,rp,ri)))  // JH jan20-06
      //if (fsh.data_fish_flags(1,fsh.parent(rr,rp,ri)))     //  PK jan14-05
      {
       // ofs << setprecision(3) << setw(12) << value(fsh.tot_catch(rr,rp,ri));
       // JH jan20-06
        ofs << setprecision(3) << setw(12) << sum(exp(value(fsh.catch_q0(rr,rp,ri))));
      }
      else
      {
        dvariable totcatch_q0=0.0;
        dvariable sv27=fsh.get_sv_region(rr,27);  //NMD 12Dec2011
        dvariable sv28=fsh.get_sv_region(rr,28);  //NMD 12Dec2011
        if (value(sv28)==3.0)  //NMD 12Dec2011
        {
          totcatch_q0=fsh.len_wt_coff/1000.* (exp(fsh.catch_q0(rr,rp,ri))*
            (pow(fsh.mean_length(rr,rp,ri),3)+
            3.0*elem_prod(fsh.mean_length(rr,rp,ri),fsh.vars(rr,rp,ri))));
        }
        else
        {
          for (int j=1;j<=fsh.nage;j++)
          {  
            totcatch_q0+=exp(fsh.catch_q0(rr,rp,ri,j))*
              normal_length_to_weight(0.5,-3.5,
              fsh.mean_length(rr,rp,ri,j),sqrt(fsh.vars(rr,rp,ri,j)),
              value(sv27),value(sv28));  //NMD 12Dec2011
          }
          totcatch_q0/=1000.;
        }
        if (totcatch_q0<=0.0) 
        {
          //cout << "negative totcatch for ir " << rr << " ip " << rp 
            //   << "  fi " << ri << endl;
        }
        ofs << setprecision(3) << setw(12) << value(totcatch_q0);
      }
    }
    ofs << endl;
  }
/*  for (i=1;i<=fsh.num_fisheries;i++)
  {
    for (int nt=1;nt<=fsh.num_fish_times(i);nt++)
    {
      int rr=fsh.realization_region(i,nt);
      int rp=fsh.realization_period(i,nt);
      int ri=fsh.realization_incident(i,nt);
      MY_DOUBLE_TYPE tmp=fsh.obs_tot_catch(rr,rp,ri);
#if !defined(NO_MY_DOUBLE_TYPE)
      if (tmp==-1.0L) tmp=0.0;
#else
      if (tmp==-1.0) tmp=0.0;
#endif
      ofs << setprecision(3) << setw(12) << tmp;
    }
    ofs << endl;
  }*/
 }
