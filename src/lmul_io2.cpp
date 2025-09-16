/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#define HOME_VERSION
#include "all.hpp"
    extern int  _version_no;
    extern int  _file_version_no;
    extern int  _ini_version_no;
    extern ofstream * clogf;

#ifndef __GNUC__
//#include <mf_menu.h>
#endif
//#include <ztc_writ.h>
#ifdef HERE
  #undef HERE
#endif
//#define HERE 
#define HERE *clogf << "reached line " << __LINE__ << " in " << __FILE__ << endl;

        static void yyy(int x){;}
extern int _data_file_version_no;
    cifstream& operator >> (cifstream& , dvar_fish_stock_history&);

    par_cifstream& operator >> (par_cifstream& pif, dvar_len_fish_stock_history& fsh)
    {
#ifndef __GNUC__
      HERE
#endif
      cifstream& cif= *(cifstream *) &pif;
      test_the_pointer();
      cif >> * (dvar_fish_stock_history *) &fsh;
      test_the_pointer();

      int version_no=fsh.parest_flags(200);
      _file_version_no=version_no;
      
      if (!cif)
      {
        cerr << "Error reading PAR file"<<endl;
        ad_exit(1);
      }
      if (!fsh.parest_flags(197))
      {
        if (version_no > 1046)
        {
          cif >> fsh.log_length_variance;
          cif >> fsh.length_rho;
          cif >> fsh.log_length_dof;
          cif >> fsh.length_psi;
          cif >> fsh.length_exp;

          cif >> fsh.log_weight_variance;
          cif >> fsh.weight_rho;
          cif >> fsh.log_weight_dof;
          cif >> fsh.weight_psi;
          cif >> fsh.weight_exp;
        }
        else
        {
          fsh.log_length_variance=0.0;
          fsh.length_rho=0.0;
          fsh.length_psi=0.0;
          fsh.length_exp=0.0;

          fsh.log_weight_variance=0.0;
          fsh.weight_rho=0.0;
          fsh.weight_psi=0.0;
          fsh.weight_exp=0.0;
        }
        if (version_no >1046 )
        {
          cif >> fsh.length_tot_exp;
          cif >> fsh.weight_tot_exp;
        }
        else
        {
          fsh.length_tot_exp=0.0;
          fsh.weight_tot_exp=0.0;
        }
        //if (fsh.parest_flags(197) || _file_version_no > 1055)
        if (_file_version_no > 1055)
        {
          cif >> fsh.cpmature_at_length;
          for (int i=fsh.cpmature_at_length.indexmin();i<=fsh.cpmature_at_length.indexmax();i++)
          {
            if (fsh.cpmature_at_length(i)<1.0e-08)
            {
              cout << "Bounds error reading cpmature_at_length(" << i
                   << ") in par file "
                   << endl << " value is " << fsh.cpmature_at_length(i) << endl;
              fsh.cpmature_at_length(i)=1.0e-08;
              cout << " Has been reset to: " << fsh.cpmature_at_length(i) << endl;
            }
            if (fsh.cpmature_at_length(i)>1.0)
            {
              cout << "Bounds error reading cpmature_at_length(" << i
                   << ") in par file "
                   << endl << " value is " << fsh.cpmature_at_length(i)<< endl;
              ad_exit(1);
            }
          }
          if (fsh.pmsd)   //NMD_Aug28_2018
          {
            int numsp=fsh.pmsd->num_species;
            for (int is=2;is<=numsp;is++)
            {
              cif >> fsh.pmsd->cpmature_at_length(is);
              for (int i=fsh.pmsd->cpmature_at_length(is).indexmin();
                       i<=fsh.pmsd->cpmature_at_length(is).indexmax();i++)
              {
                if (fsh.pmsd->cpmature_at_length(is,i)<1.0e-08)
                {
                  cout << "Bounds error reading cpmature_at_length(" << i
                       << ") in par file "
                       << endl << " value is " << fsh.pmsd->cpmature_at_length(is,i) << endl;
                  fsh.pmsd->cpmature_at_length(is,i)=1.0e-08;
                  cout << " Has been reset to: " << fsh.pmsd->cpmature_at_length(is,i) << endl;
                }
                if (fsh.pmsd->cpmature_at_length(is,i)>1.0)
                {
                  cout << "Bounds error reading multi-species cpmature_at_length(" << i 
                   << ") in par file "
                   << endl << " value is " << fsh.pmsd->cpmature_at_length(is,i)<< endl;
                   ad_exit(1);
                }
              }
            }
          }
        }
        else
        {
//          fsh.cpmature_at_length.initialize();
          fsh.cpmature_at_length=1.e-08;
          if (fsh.pmsd)
          {
//            fsh.pmsd->cpmature_at_length.initialize();
            fsh.pmsd->cpmature_at_length=1.e-08;
          }
        }   //NMD_Aug28_2018
        cif >> fsh.vb_coff(1) >> fsh.fmin1 >> fsh.fmax1;
//        cout <<  fsh.vb_coff(1) <<  fsh.fmin1 <<  fsh.fmax1 << endl;
        cif >> fsh.vb_coff(2) >> fsh.fminl >> fsh.fmaxl;
        cif >> fsh.vb_coff(3) >> fsh.rhomin >> fsh.rhomax;
	if (fsh.vb_coff(3)==0.0)
        {
          cerr << "Error reading par file vb_coff(3)=0.0 "  << endl;
          ad_exit(1);
        }
	if (fsh.vb_coff(3)==0.0)
        {
          cerr << "Error reading par file vb_coff(3)=0.0 "  << endl;
          ad_exit(1);
        }
        if (version_no > 1039)
        {
          cif >> fsh.vb_coff(4);
        }
        if (fsh.pmsd)
        {
          int numsp=fsh.pmsd->num_species;
          for (int i=2;i<=numsp;i++)
          {
            cif >> fsh.pmsd->vb_coff(i,1) >> fsh.pmsd->fmin1(i) 
                >> fsh.pmsd->fmax1(i);
            cif >> fsh.pmsd->vb_coff(i,2) >> fsh.pmsd->fminl(i) 
                >> fsh.pmsd->fmaxl(i);
            cif >> fsh.pmsd->vb_coff(i,3) >> fsh.pmsd->rhomin(i) 
                >> fsh.pmsd->rhomax(i);
            if (version_no > 1039)
            {
              cif >> fsh.pmsd->vb_coff(i,4);
            }
            
          }
        }
        if (!cif)
        {
    	  cerr << "Error reading PAR file"<<endl;
	  ad_exit(1);
        }
        cif >> fsh.vb_bias;
        if (version_no > 1001)
        {
          cif >> fsh.common_vb_bias;
          cif >> fsh.common_vb_bias_coffs;
        }
        else
        {
          fsh.common_vb_bias.initialize();
          fsh.common_vb_bias_coffs.initialize();
        }

        if (version_no>1063) 
        {
          int nr=0;
          int lry=0;
          cif >> nr >> lry;
          if (!allocated(fsh.Ninit_standard))
          {
            fsh.Ninit_standard.allocate(1,nr,1,lry);
          }
          cif >> fsh.Ninit_standard;
          if (!cif)
          {
            cerr << "Error reading PAR file"<<endl;
  	    ad_exit(1);
          }
          if (!allocated(fsh.Ninit_orthogonal))
          {
            fsh.Ninit_orthogonal.allocate(1,nr,1,lry);
          }
          cif >> fsh.Ninit_orthogonal;
          if (!cif)
          {
            cerr << "Error reading PAR file"<<endl;
  	    ad_exit(1);
          }
        }

        if (!cif)
        {
          cerr << "Error reading PAR file"<<endl;
          ad_exit(1);
        }
        if (version_no > 1020)
        {
          pif >> fsh.sv;
        }
        else if (version_no > 1005)
        {
          pif >> fsh.sv(1,20);
          fsh.sv(21,fsh.sv.indexmax()).initialize();
        }
        else if (version_no > 1002)
        {
          fsh.sv.initialize();
          for (int i=1;i<=8;i++)
          {
            pif >> fsh.sv(i);
          }
        }
        else
        {
          fsh.sv.initialize();
          pif >> fsh.sv(1) >> fsh.sv(2);
        }
        if (!cif)
        {
    	  cerr << "Error reading PAR file"<<endl;
	  ad_exit(1);
        }

        if (_data_file_version_no<3)
          fsh.sv(27)=fsh.len_wt_coff;
        else
          fsh.len_wt_coff=value(fsh.sv(27));
        if (fsh.sv(28)==0.0) fsh.sv(28)=3.0;

        if (fsh.sv(29)==0.0) fsh.sv(29)=0.9;

//        cout << "VVV " << fsh.sv << endl;
        if (fsh.pmsd)
        {
          pif >> fsh.pmsd->sv;
		  //NMD 31Oct2011
          for(int i=2;i<=fsh.pmsd->num_species; i++){
            if (fsh.pmsd->sv(i,28)==0.0) fsh.pmsd->sv(i,28)==3.0;
            if (fsh.pmsd->sv(i,29)==0.0) fsh.pmsd->sv(i,29)=0.9;
            // DF 30NOV2011
            fsh.pmsd->len_wt_coff(i)=value(fsh.pmsd->sv(i,27));
            // DF 30NOV2011
          }
          //NMD 31Oct2011
	}
        if (!cif)
        {
    	  cerr << "Error reading PAR file"<<endl;
	  ad_exit(1);
        }

        if (version_no > 1004)
        {
          pif >> fsh.growth_dev;
        }
        else
        {
          fsh.growth_dev.initialize();
        }

        cif >> fsh.var_coff(1) >> fsh.vmin1 >> fsh.vmax1;
        cif >> fsh.var_coff(2) >> fsh.vminl >> fsh.vmaxl;
        if (fsh.pmsd)
        {
          int numsp=fsh.pmsd->num_species;
          for (int i=2;i<=numsp;i++)
          {
            cif >> fsh.pmsd->var_coff(i,1) >> fsh.pmsd->vmin1(i) 
                >> fsh.pmsd->vmax1(i);
            cif >> fsh.pmsd->var_coff(i,2) >> fsh.pmsd->vminl(i) 
                >> fsh.pmsd->vmaxl(i);
          }
        }
        ivector mps = fsh.kludged_get_equilibrium_movements_per_season();
        if (!allocated(fsh.kludged_equilib_coffs))
        {
          fsh.kludged_equilib_coffs.allocate(1,fsh.age_flags(57),1,mps,
            1,fsh.num_regions,1,fsh.nage);
        }
        if (_file_version_no<1060)
        {
          fsh.kludged_equilib_coffs.initialize();
        }
        else
        {
          cif >> fsh.kludged_equilib_coffs;
          if (!cif) error_msg2("kludged_equilib_coffs");
        }
        // kludged_equilib_level_coffs;
        if (!allocated(fsh.kludged_equilib_level_coffs))
        {
          fsh.kludged_equilib_level_coffs.allocate(1,fsh.num_regions);
        }
        if (_file_version_no<1060)
        {
          fsh.kludged_equilib_level_coffs.initialize();
        }
        else
        {
          cif >> fsh.kludged_equilib_level_coffs;
          if (!cif) error_msg2("kludged_equilib_level_coffs");
        }
        // *********************************************************
        int ndeg=0; //NMD 17Nov2017
        if (version_no > 1054) //NMD 17Nov2017
        {
          cif >> ndeg;
        } //NMD 17Nov2017
        if (ndeg>0)
        {
          if (!allocated(fsh.new_orth_recr))
          {
            fsh.new_orth_recr.allocate(1,ndeg);
          }
          else 
          {
            if (fsh.new_orth_recr.indexmax() != ndeg)
            {
              cerr << "size error" << endl;
              ad_exit(1);
            }
          }
          cif >>  fsh.new_orth_recr;
        }
        if (fsh.pmsd)
        {
          int num_species=fsh.pmsd->num_species;
          for (int is=2;is<=num_species;is++)
          {
            int ndeg=0; //NMD 17Nov2017
            if (version_no > 1054) //NMD 17Nov2017
            {
              cif >> ndeg;
            } //NMD 17Nov2017
//            int ndeg;
//            cif >> ndeg;
            if (ndeg>0)
            {
              cif >> fsh.pmsd->num_new_weights(is);
              if (!allocated(fsh.pmsd->new_orth_recr(is)))
              {
                fsh.pmsd->new_orth_recr(is).allocate(1,ndeg);
              }
              else
              {
                if (fsh.pmsd->new_orth_recr(is).indexmax() != ndeg)
                {
                  cerr << "size error" << endl;
                  ad_exit(1);
                }
              }
              cif >> fsh.pmsd->new_orth_recr(is);
            }
          }
        }
        if (!cif)
        {
    	  cerr << "Error reading PAR file"<<endl;
	  ad_exit(1);
        }
        // *********************************************************
        cif >> fsh.nfmbound;
        yyy(fsh.nfmbound);
        if (fsh.nfmbound > 30)
        {
          cout <<" need to increase the dimension  of arrays"
		     " used for bounding means"<< endl;
	 ad_exit(1);
        }
        for (int i=1; i<=fsh.nfmbound; i++)
        {
          cif >> fsh.ifper(i) >> fsh.ifinc(i) >> fsh.iageclass(i)
	    >>fsh.fmmin(i) >> fsh.fmmax(i) >> fsh.pdown(i)
	    >> fsh.pup(i);
        }
        cif >> fsh.D;

        if (version_no > 1065)
        {
	  int tmpyr1;
          cif >> tmpyr1;
//	  cout << " Read in year 1 from par version 1066: " << tmpyr1 << endl;
        }
	
        if (!cif)
        {
    	  cerr << "Error reading PAR file"<<endl;
	  ad_exit(1);
        }

        if (version_no > 1066)
        {
          cif >> fsh.historical_version;
        }
	
        if (!cif)
        {
    	  cerr << "Error reading PAR file"<<endl;
	  ad_exit(1);
        }
	
      }
      else
      {
        if (!allocated(fsh.kludged_equilib_coffs))
        {
          ivector mps=
            fsh.kludged_get_equilibrium_movements_per_season();
          fsh.kludged_equilib_coffs.allocate(1,fsh.age_flags(57),1,mps,
             1,fsh.num_regions,1,fsh.nage);
          fsh.kludged_equilib_coffs.initialize();
        }
        //cifstream cif("generic.par");
        if (!allocated(fsh.kludged_equilib_level_coffs))
        {
          fsh.kludged_equilib_level_coffs.allocate(1,fsh.num_regions);
          fsh.kludged_equilib_level_coffs.initialize();
        }
        if (_ini_version_no>2)
        {
          cif >> fsh.cpmature_at_length;
          if (fsh.pmsd)   //NMD_Aug30_2018
          {
            int numsp=fsh.pmsd->num_species;
            for (int i=2;i<=numsp;i++)
            {
              cif >> fsh.pmsd->cpmature_at_length(i);
            }
          }	  
        }
        else
        {
          fsh.cpmature_at_length.initialize();
          if (fsh.pmsd)   //NMD_Aug30_2018
          {
            fsh.pmsd->cpmature_at_length.initialize();
          }
        }

        cif >> fsh.vb_coff(1) >> fsh.fmin1 >> fsh.fmax1;
        cif >> fsh.vb_coff(2) >> fsh.fminl >> fsh.fmaxl;
        cif >> fsh.vb_coff(3) >> fsh.rhomin >> fsh.rhomax;
	if (fsh.vb_coff(3)==0.0)
        {
          cerr << "Error reading par file vb_coff(3)=0.0 "  << endl;
          ad_exit(1);
        }
        if (_ini_version_no>4)
        {
          cif >> fsh.vb_coff(4);
        }
        else
        {
          fsh.vb_coff(4)=0.0; 
        }
//        fsh.vb_coff(4)=0.0; 
        if (fsh.pmsd)
        {
          int numsp=fsh.pmsd->num_species;
          for (int i=2;i<=numsp;i++)
          {
            cif >> fsh.pmsd->vb_coff(i,1) >> fsh.pmsd->fmin1(i) //NMD 3Dec2011
                >> fsh.pmsd->fmax1(i);
            cif >> fsh.pmsd->vb_coff(i,2) >> fsh.pmsd->fminl(i) //NMD 3Dec2011
                >> fsh.pmsd->fmaxl(i);
            cif >> fsh.pmsd->vb_coff(i,3) >> fsh.pmsd->rhomin(i) //NMD 3Dec2011
                >> fsh.pmsd->rhomax(i);
            if (_ini_version_no>4)
            {
              cif >> fsh.pmsd->vb_coff(i,4);
            }
            else
            {
              fsh.pmsd->vb_coff(i,4)=0.0; 
            }
//            fsh.pmsd->vb_coff(i,4)=0.0;                           //NMD 3Dec2011
          }
        }

        if (!cif)
        {
    	  cerr << "Error reading PAR file"<<endl;
	  ad_exit(1);
        }

        fsh.vb_bias.initialize();
        fsh.common_vb_bias.initialize();
        fsh.common_vb_bias_coffs.initialize();
        fsh.sv.initialize();
        //       a             b       in      w = al^b
        cif >> fsh.sv(27) >> fsh.sv(28);
        // steepness
        if (_ini_version_no>0) cif >> fsh.sv(29);   //NMD_5Mar2013
        if (fsh.pmsd)
        {
          for (int i=2;i<=fsh.pmsd->num_species;i++)
          {
            cif >> fsh.pmsd->sv(i,27);
            cif >> fsh.pmsd->sv(i,28);
            if (_ini_version_no>0) cif >> fsh.pmsd->sv(i,29);  //NMD_5Mar2013
          }
        }
        if (!cif)
        {
    	  cerr << "Error reading PAR file"<<endl;
	  ad_exit(1);
        }
        
        fsh.growth_dev.initialize();
        cif >> fsh.var_coff(1) >> fsh.vmin1 >> fsh.vmax1;
        cif >> fsh.var_coff(2) >> fsh.vminl >> fsh.vmaxl;
        if (fsh.pmsd)
        {
          int numsp=fsh.pmsd->num_species;
          for (int i=2;i<=numsp;i++)
          {
            cif >> fsh.pmsd->var_coff(i,1) >> fsh.pmsd->vmin1(i) //NMD 3Dec2011
                >> fsh.pmsd->vmax1(i);
            cif >> fsh.pmsd->var_coff(i,2) >> fsh.pmsd->vminl(i) //NMD 3Dec2011
                >> fsh.pmsd->vmaxl(i);
          }
        }
        if (!cif)
        {
    	  cerr << "Error reading PAR file"<<endl;
	  ad_exit(1);
        }


        cif >> fsh.nfmbound;
        if (!cif)
        {
    	  cerr << "Error reading PAR file"<<endl;
	  ad_exit(1);
        }
        if (fsh.nfmbound > 30)
        {
          cout <<" need to increase the dimension  of arrays"
		     " used for bounding means"<< endl;
	 ad_exit(1);
        }
        for (int i=1; i<=fsh.nfmbound; i++)
        {
          cif >> fsh.ifper(i) >> fsh.ifinc(i) >> fsh.iageclass(i)
	    >>fsh.fmmin(i) >> fsh.fmmax(i) >> fsh.pdown(i)
	    >> fsh.pup(i);
        }
        //cif >> fsh.D;
        if (!cif)
        {
    	cerr << "Error reading PAR file"<<endl;
	ad_exit(1);
        }
      }
      if (fsh.age_flags(61))
      {
        ifstream ifs("relbio.dat");
        if (!ifs)
        {
          cerr << "Error trying to open file relbio.dat" << endl;
        }
        else
        {
          ifs >> fsh.parest_flags(180); 
          ifs >> fsh.biomass_index; 
          int i;
          for ( i=1;i<=fsh.num_fisheries;i++)
          {
            MY_DOUBLE_TYPE sum=0.0;
            int ii=0;
            int j;
            for (j=1;j<=fsh.num_fish_times(i);j++)
            {
              int rr=fsh.realization_region(i,j);
              int rp=fsh.realization_period(i,j);
              int ri=fsh.realization_incident(i,j);
              
              if (fsh.biomass_index(i,j)>0)
              {
                ii++;
                sum+=fsh.biomass_index(i,j);
              }
            }
            sum/=ii;
            for (j=1;j<=fsh.num_fish_times(i);j++)
            {
              int rr=fsh.realization_region(i,j);
              int rp=fsh.realization_period(i,j);
              int ri=fsh.realization_incident(i,j);
              
              if (fsh.biomass_index(i,j)>0)
              {
                fsh.biomass_index(i,j)/=sum;
              }
            }
          }
          for (i=1;i<=fsh.num_fisheries;i++)
          {
            for (int j=1;j<=fsh.num_fish_times(i);j++)
            {
              int rr=fsh.realization_region(i,j);
              int rp=fsh.realization_period(i,j);
              int ri=fsh.realization_incident(i,j);
              fsh.ORB(rr,rp,ri)=fsh.biomass_index(i,j);
            }
          }
        }
        if (!ifs)
        {
          cerr << "Error trying to read biomass data from file relbio.dat" 
               << endl;
        }
      }

      ivector ff71=column(fsh.fish_flags,71);
      //ivector ff72=column(fsh.fish_flags,72);
      //ivector ff73=column(fsh.fish_flags,73);
      //ivector ff74=column(fsh.fish_flags,74);
      //int sff71=sum(ff71);
      //if (sff71)
      {
        read_selectivity_blocks(fsh,ff71);
      }

      test_the_pointer();
      return pif;
    }

#undef HOME_VERSION
