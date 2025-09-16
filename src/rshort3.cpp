/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
#ifdef __MSC__
  MY_DOUBLE_TYPE mfexp(const MY_DOUBLE_TYPE& );
  dvariable mfexp(const prevariable& );
  dvar_vector mfexp(const dvar_vector& );
  dvector mfexp(const dvector& );
#else
  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
  dvariable mfexp(prevariable& );
  //dvar_vector mfexp(dvar_vector& );
  //dvector mfexp(dvector& );
#endif

#if !defined(ZFEV)
extern dvar_len_fish_stock_history * pcfsh;
#endif
  extern dvar_vector * psv;

   void   uxxs(int u){;}

MY_DOUBLE_TYPE rmax_penalty_weight=1.e+5; 
MY_DOUBLE_TYPE rmax_multiplier= 0.96;


void check_sanity(ivector &t)
{
  int mmin=t.indexmin();
  int mmax=t.indexmax();
  for (int i=mmin;i<mmax;i++)
  {
    if (t(i) != t(i+1))
    {
      cerr << "sanity error in move_periods" << endl;
      ad_exit(1);
    }
  }
}
void dvar_fish_stock_history::set_totpop(void)
{
  totpop=totpop_coff;
}

void setup_diffusion(dvar3_array& Dad,int nage,int num_regions,
  imatrix& Dflags, dvar_vector& diff_coffs,dvar_vector& diff_coffs2,
  dvar_vector& diff_coffs3,ivector& age_flags,pmulti_species_data & pmsd);

int dvar_fish_stock_history::get_numyears_average(void)
{ 
  int nya=1; // num_years_for_averaging total mortality
  /*
  //if (age_flags(125))
  {
    nya=age_flags(95); // num_years_for_average
    if (nya==0)
    {
      cerr << "Error age_flags(95) must be >0 " << endl;
      ad_exit(1);
    }
  }
  */
//NMD_25May2021
  nya=age_flags(95); // num_years_for_average
  if (age_flags(94)==2)
  {
    if (nya==0)
    {
      cerr << "Error age_flags(95) must be >0 " << endl;
      cerr << "when age_flags(94)==2 option is used " << endl;
      ad_exit(1);
    }
    if (age_flags(125))
    {
      if (nya==0)
      {
        cerr << "Error age_flags(95) must be >0 " << endl;
        cerr << "when age_flags(94)==2 option and " << endl;
        cerr << "when age_flags(125) option is used " << endl;
        ad_exit(1);
      }
    }
  }
//NMD_25May2021
  return nya;
}

void check_sanity(ivector &t, dvar3_array& tg, ivector& rip,int it,
  dvar_fish_stock_history& fsh)
{
  ivector& itp=fsh.initial_tag_period(it);
  check_sanity(t);
# if !defined(OPT_LIB)
    
    int mmin=tg.indexmin();
    int mmax=tg.indexmax();
    int err=0;
    int i;
    for (i=mmin;i<mmax;i++)
    {
      if (tg(i,rip(i)).indexmin() != tg(i+1,rip(i+1)).indexmin())
      {
        err=1;
        break;
      }
    }
    if (err)
    {
      cerr << "sanity error in youngest age class for tag group " 
           << it << endl
           << "for regions " << i << " and " << i+1
           << " mimum ages are "
           <<  tg(i,rip(i)).indexmin() << " and "
           <<  tg(i+1,rip(i+1)).indexmin() 
           << " fishing periods are " << endl
           <<  rip << endl
           << " initial tag periods for this tag group are "
           << itp << endl;
   
      cerr << "years are" << endl;
      for (i=mmin;i<=mmax;i++)
      {
        cerr << fsh.year(i,rip(i)) << " ";
      }
      cerr << endl;
      cerr << "months are" << endl;
      for (i=mmin;i<=mmax;i++)
      {
        cerr << fsh.month(i,rip(i)) << " ";
      }
      cerr << endl;

      ad_exit(1);
    }
     /*
      cerr << " fishing periods are " << endl
           <<  rip << endl;

      cerr << "years are" << endl;
      for (i=mmin;i<=mmax;i++)
      {
        cerr << fsh.year(i,rip(i)) << " ";
      }
      cerr << endl;
      cerr << "months are" << endl;
      for (i=mmin;i<=mmax;i++)
      {
        cerr << fsh.month(i,rip(i)) << " ";
      }
      cerr << endl;
     */

#  endif
}
  int sip;
  
void dvar_len_fish_stock_history::catch_equations_calc(dvar_vector& sv,
  ivector* pq_flag,const prevariable & _fzpen)
{
  ADUNCONST(prevariable,fzpen)
  //greport("beginning catch_equations_calc");
  tmprecr.initialize();
  tmpinitpop.initialize();
  MY_DOUBLE_TYPE one_plus=1.e0+1.e-10;
  if (af170q0==0)
  {
    tot_mort.initialize();
    N.initialize();
  }
  else
  {
    tot_mort_q0.initialize();
    N_q0.initialize();
  }
  // calculate the totalfishing mortality and survival rates for each
  // fishing period in each region
  //cout <<"rshort3.cpp " <<  num_fish_periods << endl;
  int ir;
  dvar_vector * ps=0;
  dvar_vector * ptm=0;
  dvar_matrix * pfm=0;

  int nya=get_numyears_average(); 
  /*
  int nya=1; // num_years_for_average
  if (age_flags(125))
  {
    nya=age_flags(95); // num_years_for_average
    if (nya==0)
    {
      cerr << "Error age_flags(95) must be >0 " << endl;
      ad_exit(1);
    }
  }
  */

  for (ir=1;ir<=num_regions;ir++)
  {
    //cout <<"rshort3.cpp " << fraction(ir) << endl;
    int tmp_nfp;
    if (do_fishery_projections_flag==0)
      tmp_nfp=num_real_fish_periods(ir);
    else
      tmp_nfp=num_fish_periods(ir);

    for (int ip=1;ip<=tmp_nfp;ip++)  // Loop over fishing periods
    {
      int yr=year(ir,ip);
      // VVVVVV
      if (age_flags(125) && yr>nya) break;
      if (af170q0==0)
      {
        ptm=&(tot_mort(ir,ip));
        pfm=&(fish_mort(ir,ip));
        ps=&(survival(ir,ip));
      }
      else
      {
        ptm=&(tot_mort_q0(ir,ip));
        pfm=&(fish_mort_q0(ir,ip));
        ps=&(survival_q0(ir,ip));
      }
      

      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) // Loop over fishing
      {                                         // incidents for this period
        (*ptm)+=mfexp((*pfm)(fi));
      }
      if (age_flags(115))   // put bounds on Z (acutally total F)
      {                     // for SS2 parameterization
        bounds_for_ss2(*ptm,*pfm,fzpen);
      }
      
      int af115=get_age_flags_region(ir,115);
      if (!pmsd)
      {
        if (af115==0)   
        {                     
          (*ptm)+=mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
        }
        else
        {  // for SS2 parameterization
           // this is Z=-log(S)=-log(1-sum F)+M
           (*ptm)=-log(1-*ptm)+mfexp(nat_mort(year(ir,ip))+fraction(ir,ip));
        }
      }
      else
      {
        if (af115==0)   
        {                     
          (*ptm)+=
            mfexp(get_nat_mort_region(ir)(year(ir,ip))+fraction(ir,ip));
        }
        else
        {  // for SS2 parameterization
           // this is Z=-log(S)=-log(1-sum F)+M
           (*ptm)=-log(1-*ptm)
             +mfexp(get_nat_mort_region(ir)(year(ir,ip))+fraction(ir,ip));
        }
      }
    
      (*ps)=mfexp(-(*ptm));
    }
  }

  do_fish_mort_intermediate_calcs();

  if (num_regions>1 && !age_flags(114))
  {
    setup_diffusion();
  }
     if (ff263flag==0)
     {
#if !defined(ZFEV)
       if (!parest_flags(143))
#endif
       {
         if (!parest_flags(175))
           mean_length_calc(0);
         else
           mean_length_calcx();
       }
     }
#if !defined(ZFEV)
     else
     {
       // do this once at the beginning of the run
       if (parest_flags(144))
       {
         pcfsh->mean_length_calc();
         pcfsh->variance_calc();
       }
     }
#endif

  //mean_weight=mean_weight_calcs();

  if (age_flags(177))
  {
    set_totpop();
    if (pmsd)
    {
      pmsd->totpop=pmsd->totpop_coff;
    }

    if (pq_flag)
    {
      age_flags(192)=1;
    }
    fzpen+=get_initial_population(sv,0,pq_flag);
    age_flags(192)=0;
  }
  else
  {
    if (af170q0==0)
    {
      get_population_multipliers(sv,pq_flag);
      fzpen+=get_initial_population(sv,0,pq_flag);
    }
    else
    {
      if (age_flags(171)==0)
      {
        get_population_multipliers(sv,pq_flag);
      }
      // get the initial age structure
      if (!age_flags(94))
      {
        get_initial_age_structure(totpop,sv);
      }
      else
      {  
        xget_initial_age_structure_equilibrium();
      }
    }
  }
  if (projection_sim_index>0 && do_fishery_projections_flag)
  {
    //ofstream ofs("projected_randomized_catches");
    int num_times=age_flags(20);
    random_number_generator rng(217+projection_sim_index);
    randrec.allocate(num_real_years,nyears);
    MY_DOUBLE_TYPE sd=age_flags(21)/100.;
    gradient_structure::set_NO_DERIVATIVES(); 
    {
      //sd=0.0;
      randrec.fill_randn(rng);
      randrec=exp(sd*randrec);
  
      get_numbers_at_age(sv,pq_flag);

      if(!pq_flag && parest_flags(191))   //NMD_8May2018
      {
        *proj_output_files[0] << "#simulation " << projection_sim_index << endl;
        *proj_output_files[0] << "#catch " << endl;
        for (int ir=1;ir<=num_regions;ir++)
        {
          *proj_output_files[0]  << "#region " << ir << endl;
          for (int ip=1;ip<=num_fish_periods(ir);ip++)
          {
            if (ip==num_real_fish_periods(ir)+1)
              *(proj_output_files[0]) << "#projected catch" << endl;
            for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++)
            {
              *(proj_output_files[0]) << exp(catch(ir,ip,fi)) << endl;
            }
          }
        }
        *proj_output_files[0] << "#num fish in pop " << endl;
        for (int ir=1;ir<=num_regions;ir++)
        {
          *(proj_output_files[0]) << "#region " << ir << endl;
          for (int ip=1;ip<=num_fish_periods(ir);ip++)
          {
            if (ip==num_real_fish_periods(ir)+1)
              *(proj_output_files[0]) << "#projected numbers" << endl;
            *(proj_output_files[0])  << exp(num_fish(ir,ip)) << endl;
          }
        }
      }    //NMD_8May2018
    }
  }

  if (age_flags(94)!=3  && !age_flags(182) )
    get_annual_recruitment_flag=1;

  get_numbers_at_age(sv,pq_flag);

  get_annual_recruitment_flag=0;
      
  MY_DOUBLE_TYPE diff=0.0;
  int icount=0;
  dmatrix oldvtmp;  //NMD_13Apr2021
  if (!pmsd)
  {
    oldvtmp.allocate(1,1,1,last_real_year+1);   //NMD_13Apr2021
  }
  else
  {
    oldvtmp.allocate(1,pmsd->num_species,1,last_real_year+1);   //NMD_13Apr2021
  }
  oldvtmp.initialize();
  do
  {
    if ( af170q0==1 && age_flags(171)==1)
    {
      dmatrix vtmp;      
      if (age_flags(94)==3  || age_flags(182) )
      {
        get_annual_recruitment_flag=1;
        vtmp=value(get_numbers_at_age(sv,pq_flag));
        get_annual_recruitment_flag=0;
        if (icount>0)
        {
          diff=norm2(vtmp-oldvtmp);
        }
        else
        {
          diff=1.0;
        }
        icount ++;
        oldvtmp=vtmp;
      }
    }
  }    
  while (diff>1.e-8);

  if (age_flags(180))
  {   
    (*clogf)  <<  "num_fish" << endl;
    (*clogf)  <<  num_fish << endl;
  }
  if (age_flags(180))
  {   
    (*clogf)  <<  "catch" << endl;
    (*clogf)  <<  catch << endl;
  }
}

void dvar_fish_stock_history::get_initial_age_structure(
  const dvariable& totpop,dvar_vector& sv)
{
  dvar3_array * pN=0;
  if (af170q0==0)
  {
    pN=&N;
  }
  else
  {
    pN=&N_q0;
  }
  //pmsd_error();
  // get the initial age structure as independent parameters
  int ng=nage;
  for (int ir=1;ir<=num_regions;ir++)
  {
    //if (ir>1) N(ir,1,1)=0.0;
    if (pmsd)
    {
//      ng=pmsd->region_species_pointer(ir);
      ng=pmsd->nage_by_region(ir);  //NMD_jan17-19
    }
    for (int j=2;j<=ng;j++)
    {
      (*pN)(ir,1,j)=initpop(ir,j)-initmean-rec_init_diff+totpop;
    }
  }
  // here is where we "diffuse " the fish to move them between the regions
  if (num_regions>1) 
  {
    do_the_diffusion(1,sv,*pN);
  }
}
dvar_matrix fast_diffusion_calcs(int nage,int num_regions,dvar3_array& N,
  dvar3_array& Dad,int year,pmulti_species_data & pmsd);

void dvar_fish_stock_history::do_the_diffusion(int year,dvar_vector& sv,dvar3_array& N,ofstream * pofs)
{
  int i=1;
  if(parest_flags(357)) i=2;
//  if(fish_nonmv_flag) i=2;    //NMD_14Sep2018
  //if(parest_flags(357)) i=move_map(1);
  dvar_matrix tmp=fast_diffusion_calcs(nage,num_regions,N,Dad(i),year,pmsd);
  //check_derivative_values("epos1");
  //for (int ir=1;ir<=num_regions;ir++)
  int it=0;
  if (pmsd) it=pmsd->tag_index;
  if (it)
  {
    for (int ir=tag_region_bounds(1,it);ir<=tag_region_bounds(2,it);ir++)
    {
      N(ir,year)=log( 5.e-10+tmp(ir) );
    }
  }
  else
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      N(ir,year)=log( 5.e-10+tmp(ir) );
    }
  }
  //check_derivative_values("fpos1");
}

dvariable setup_alpha(int j,dvar_vector& sv,int nage,ivector& age_flags)
{
  dvariable fj;
  MY_DOUBLE_TYPE tmp=double(j-1)/(nage-1);
  if (!age_flags(66))
  {
   fj=2.*tmp-1.;
  }
  else
  {
   if (j>1)
   {
    fj=2.*pow(tmp,exp(sv(15))) - 1.;
   }
   else
   {
    fj=-1.0;
    }
  }
  if (age_flags(64))
  {
    if (fj<-1.e-8)
    {
      fj=-pow(-fj,exp(sv(14)));
    }
    else if (fj>1.e-8)
    {
      fj=pow(fj,exp(sv(14)));
    }
  }
  dvariable alpha=sv(10)*exp(sv(13)*fj);
  return alpha;
}

void dvar_fish_stock_history::get_initial_population_test(dvar_vector& sv)
{
  N.initialize();
  int ir;
  for (ir=1;ir<=num_regions;ir++)
  {
    for (int iy=1;iy<=nyears;iy++)
    {
      N(ir,iy)(2,nage)=-20.0;
      if (ir==1)
        N(ir,iy,1)=log(1000.);
      else
        N(ir,iy,1)=-20.0;
    }
  }
  for (ir=1;ir<=num_regions;ir++)
  {
    num_fish(ir,1)(2,nage)=-20.0;
    if (ir==1)
      num_fish(ir,1,1)=log(1000.0);
    else
      num_fish(ir,1,1)=-20.0;
  }
}

void dvar_fish_stock_history::bounds_for_ss2(dvar_vector & TF,
  const dvar_matrix & _LF,const prevariable & _fzpen)
{
  ADUNCONST(prevariable,fzpen)
  ADUNCONST(dvar_matrix,LF)
  MY_DOUBLE_TYPE rmax=0.7;
  if (age_flags(116)!=0)
  {
    rmax=age_flags(116)/100.;
  }
  int nf=LF.indexmax();

  for (int j=1;j<=nage;j++)
  {
    if (value(TF(j))>rmax_multiplier*rmax)
    {
      fzpen+=rmax_penalty_weight*square(square(log(TF(j)/(rmax_multiplier*rmax))));
    }
    if (value(TF(j))>rmax)
    {
      age_flags(181)=2;
      //fpen+=100.*square(TF(j)-rmax);
      dvariable dd=TF(j)-rmax;
      dvariable ppen=0.0;
      dvariable TFstar=rmax-posfun(0.2-dd,0.1,ppen)+0.2;
  
      dvariable lambda=TFstar/TF(j);
  
      for (int fi=1;fi<=nf;fi++)
      {
        LF(fi,j)*=lambda;
      }
      TF(j)*=lambda; 
    }
  }
}

void setup_diffusion(dvar3_array& Dad,int nage,int num_regions,
  imatrix& Dflags, dvar_vector& diff_coffs,dvar_vector& diff_coffs2,
  dvar_vector& diff_coffs3,ivector& age_flags,pmulti_species_data & pmsd)
{
  if (age_flags(184)) 
  {
    if(pmsd) cerr << "not implemented for multi-species" << endl;
    ad_exit(1);
  }
  
  Dad.initialize();
  int ns=pmsd->num_species;
  int nrr=pmsd->num_real_regions;
  if (!age_flags(89))
  {
    int ii=1;
    for (int is=1;is<=ns;is++)
    {
      int rmin=pmsd->region_bounds(is,1);
      int rmax=pmsd->region_bounds(is,2);
      dvar_matrix td=Dad(1).sub(rmin,rmax).shift(1);
      
      for (int i=1;i<=nrr;i++)
      {
        td(i,i)=1;
        for (int j=1;j<=nrr;j++)
        {
          if (Dflags(i,j)) 
          {
            td(i,i)+=diff_coffs(ii);
            td(j,i)-=diff_coffs(ii++);
          }
        }
      }
      Dad(1).sub(rmin,rmax).shift(1)=inv(td);
      //cout <<"rshort3.cpp " << colsum(Dad(1)) << endl;
      for (int j=2;j<=nage;j++)
      {
        Dad(j).sub(rmin,rmax).shift(1)=
          Dad(1).sub(rmin,rmax).shift(1);
      }
    }
  }
  else
  {
    cout << "Not done yet" << endl;
    ad_exit(1);
    for (int is=1;is<=ns;is++)
    {
    //  const MY_DOUBLE_TYPE ninv=1.0/nage;
      const MY_DOUBLE_TYPE ninv=1.0/(nage-1);
  
      for (int j=1;j<=nage;j++)
      {
        int ii=1;
        dvariable fj=-1.+2.*(j-1)*ninv;
        dvar_matrix& td=Dad(j);
        for (int i=1;i<=num_regions;i++)
        {
          td(i,i)=1;
          for (int jj=1;jj<=num_regions;jj++)
          {
            if (age_flags(91))
            {
              if (fj<-1.e-8)
              {
                fj=-pow(-fj,exp(diff_coffs3(jj)));
              }
              else if (fj>1.e-8)
              {
                fj=pow(fj,exp(diff_coffs3(jj)));
              }
  	  }
            if (Dflags(i,jj)) 
            {
              dvariable tmp;
              if (!age_flags(91))
  	    {  
                tmp=diff_coffs(ii)*exp(diff_coffs2(ii)*fj);
  	    }
  	    else
  	    {
                if (fj<-1.e-8)
                {
                  fj=-pow(-fj,exp(diff_coffs3(ii)));
                }
                else if (fj>1.e-8)
                {
                  fj=pow(fj,exp(diff_coffs3(ii)));
                }
                tmp=diff_coffs(ii)*exp(diff_coffs2(ii)*fj);
  	    }
              td(i,i)+=tmp;
              td(jj,i)-=tmp;
              ii++;
            }
          }
        }
        td=inv(td);
        {
          dvar_matrix tdtrans=trans(td);
          int mmin=tdtrans.indexmin();
          int mmax=tdtrans.indexmax();
          for (int i=mmin;i<=mmax;i++)
          {
            tdtrans(i)+=1.e-7;
            tdtrans(i)/=sum(tdtrans(i));
          }
          td=trans(tdtrans);
        }
        //cout <<"rshort3.cpp " << colsum(td) << endl;
      }
    }
  }
}


void setup_diffusion(dvar3_array& Dad,int nage,int num_regions,
  imatrix& Dflags, dvar_vector& diff_coffs,dvar_vector& diff_coffs2,
  dvar_vector& diff_coffs3,ivector& age_flags)
{
  
  cout << "Error"<< endl;
  ad_exit(1);
  Dad.initialize();
  if (!age_flags(89))
  {
    int ii=1;
    dvar_matrix& td=Dad(1);
    for (int i=1;i<=num_regions;i++)
    {
      td(i,i)=1;
      for (int j=1;j<=num_regions;j++)
      {
        if (Dflags(i,j)) 
        {
          td(i,i)+=diff_coffs(ii);
          td(j,i)-=diff_coffs(ii++);
        }
      }
    }
    Dad(1)=inv(td);
    //cout <<"rshort3.cpp " << colsum(Dad(1)) << endl;
    for (int j=2;j<=nage;j++)
    {
      Dad(j)=Dad(1);
    }
  }
  else
  {
  //  const MY_DOUBLE_TYPE ninv=1.0/nage;
    const MY_DOUBLE_TYPE ninv=1.0/(nage-1);

    for (int j=1;j<=nage;j++)
    {
      int ii=1;
      dvariable fj=-1.+2.*(j-1)*ninv;
      dvar_matrix& td=Dad(j);
      for (int i=1;i<=num_regions;i++)
      {
        td(i,i)=1;
        for (int jj=1;jj<=num_regions;jj++)
        {
          if (age_flags(91))
          {
            if (fj<-1.e-8)
            {
              fj=-pow(-fj,exp(diff_coffs3(jj)));
            }
            else if (fj>1.e-8)
            {
              fj=pow(fj,exp(diff_coffs3(jj)));
            }
	  }
          if (Dflags(i,jj)) 
          {
            dvariable tmp;
            if (!age_flags(91))
	    {  
              tmp=diff_coffs(ii)*exp(diff_coffs2(ii)*fj);
	    }
	    else
	    {
              if (fj<-1.e-8)
              {
                fj=-pow(-fj,exp(diff_coffs3(ii)));
              }
              else if (fj>1.e-8)
              {
                fj=pow(fj,exp(diff_coffs3(ii)));
              }
              tmp=diff_coffs(ii)*exp(diff_coffs2(ii)*fj);
	    }
            td(i,i)+=tmp;
            td(jj,i)-=tmp;
            ii++;
          }
        }
      }
      td=inv(td);
      {
        dvar_matrix tdtrans=trans(td);
        int mmin=tdtrans.indexmin();
        int mmax=tdtrans.indexmax();
        for (int i=mmin;i<=mmax;i++)
        {
          tdtrans(i)+=1.e-7;
          tdtrans(i)/=sum(tdtrans(i));
        }
        td=trans(tdtrans);
      }
      //cout <<"rshort3.cpp " << colsum(td) << endl;
    }
  }
}
