/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
#include "here.h"
#include "globals.h"

  void  cblas_matrix_prod(int n,dvector& u,dvector & v,dvector & prod);
  dvar_matrix reg_vbiocalc(dvar_len_fish_stock_history& fsh);
  dvar_matrix unnormalized_reg_vbiocalc(dvar_len_fish_stock_history& fsh);
  dvar_matrix unnormalized_adult_reg_vbiocalc(dvar_len_fish_stock_history& fsh);

  dvariable get_dd_agedif(dvar_vector& sv,dvar_vector& rs);
  dvar_vector sbcalc(dvar_len_fish_stock_history& fsh);
  void do_average_exploitation(dvar_len_fish_stock_history& fsh,int& ii,
    int& num_grads,dvar_vector& dep_vars,ofstream& ofl);

   void depvars_bounds_check(int i,int n);

  int dep_gradients_calc2_for_projections(dvar_len_fish_stock_history& fsh, int noeff_sw)   //NMD 22Feb2012
  {
// NMD 22Feb2012
    ofstream ofl(noeff_sw == 1 ? "deplabel_noeff.tmp" : "deplabel.tmp");
    ofstream tmpout(noeff_sw == 1 ? "tmp_noeff.out"  :  "tmp.out");
// NMD 22Feb2012
    int num_grads=500;
    int ii=1;
    dvar_vector dep_vars(1,num_grads);
    //  new stuff for projections Df 16/02/2010
    
    for (int ir=1;ir<=fsh.num_regions;ir++) 
    {
      int ip=fsh.num_real_fish_periods(ir)+1;
      for (int j=2;j<=fsh.nage;j++)
      {
        if (ip > fsh.num_fish(ir).indexmax())
        {
          cerr << "index too big" << endl;
          ad_exit(1);
        }
// NMD_8Sep2015  - place no-fished population numbers
        if (noeff_sw==1)
        {
          dep_vars(ii) << fsh.num_fish_q0(ir,ip,j);
        } 
        else
        {
          dep_vars(ii) << fsh.num_fish(ir,ip,j);
        }
// NMD_8Sep2015
//        dep_vars(ii) << fsh.num_fish(ir,ip,j);
        ofl << dep_vars(ii++) << endl; 
        depvars_bounds_check(ii,num_grads);
        ofl << adstring("num_fish(")  
           + str(ir) +adstring(",")
           + str(ip) +adstring(",")
           + str(j)  + adstring(")") << endl;
      }
    }
   
    // save values of dependent variables
    int ndep=ii-1;
    mfglobals.dep_vars_values.allocate(1,ndep);
    int ij;
// NMD 22Feb2012
    uostream uos(noeff_sw == 1 ? "proparams_noeff.val" :  "proparams.val");
// NMD 22Feb2012
    uos << ndep;
    for (ij=1;ij<=ndep;ij++)
    {
      uos << dep_vars(ij);
      mfglobals.dep_vars_values(ij) = value(dep_vars(ij));
    }
    // save the historical recruitment in case it is needed
// NMD 22Feb2012
    if (noeff_sw == 0)
    {
      ofstream ofs("histrec"); 
      for (int ir=1;ir<=fsh.num_regions;ir++)
      {
        for (int iy=1;iy<=fsh.last_real_year;iy++)
        {
          ofs << fsh.N(ir,iy,1) << " ";
        }
      }
    }
// NMD 22Feb2012
    return ndep;
  }

void make_covariance_report_for_projections(const dvar_fish_stock_history& fsh)
{
//NMD 22Feb2012
  int noeff_sw = 0;
  ivector q_flag=column(fsh.fish_flags,55);
  if (sum(q_flag)) noeff_sw = 1;
//NMD 22Feb2012
  int nvar;
  int nvar1;
  int num_deps;
  int i,j;

  adstring current_path;
  adstring root;

  //getpaths(current_path,root);

  adstring fname= ad_root+".hes";

  uistream ifs(fname);
//NMD30Mar2020  -skip Hessian calcs needed for dependent variables
// but create temporary declaration of Hessinv
  dmatrix Hessinv;   //NMD30Mar2020
  if (fsh.parest_flags(237)==1)
  {
    if (!ifs)
    {
      cerr << "Error trying to open file " << fname << endl;
    }
    ifs >> nvar;
    if (!ifs)
    {
      cerr << "Error reading nvar from " << fname << endl;
    }
    int test=0;
    if (test)
    {
      nvar=100;
    }
    dvector vHess;
    dmatrix Hess=make_dmatrix(1,nvar,1,nvar,vHess);
    dvector row(1,nvar);

    if (!test)
    {
      ifs >> Hess;
    }
    else
    {
      random_number_generator rng(101);
      Hess.fill_randn(rng);
      Hess=Hess*trans(Hess);
    }
    if (!ifs)
    {
      cerr << "Error reading hessian from " << fname << endl;
    }
    for (i=1;i<=nvar;i++)
    {
      for (j=1;j<=i;j++)
      {
	Hess(i,j)=0.5*(Hess(i,j)+Hess(j,i));
	Hess(j,i)=Hess(i,j);
      }
    }
    for (i=1;i<=nvar;i++)
    {
      for (j=1;j<=nvar;j++)
      {
	if (Hess(i,j)!=0) row(i)=1;
      }
    }

    for (i=1;i<=nvar;i++)
    {
      if (row(i)==0)
      {
	cout << " Error there was a zero row in row " << i << endl;
	exit(1);
      }
    }
    //}

    dvector veigvecs;
    dmatrix eigvecs=make_dmatrix(1,nvar,1,nvar,veigvecs);
    dvector eigvals;
    int ierr=0;
    lapack_symmetric_eigen(nvar,vHess,veigvecs,eigvals,ierr);

    dvector vHessinv;
    //    dmatrix Hessinv=make_dmatrix(1,nvar,1,nvar,vHessinv);
    Hessinv=make_dmatrix(1,nvar,1,nvar,vHessinv);  //NMD_30Mar2020

    if (ierr)
    {
      cerr << "Error calculating eigenvectors in "
	   << " lapack_symmetric_eigen(Hess,eigvecs,eigvals,ierr)"
	   << endl;
      cerr << "If you want to proceed anyway enter a 1  Otherwise enter 0"
	   << endl;
      int nn=0;
      cin  >>  nn;
      if (nn !=1) ad_exit(1); 
    }

    MY_DOUBLE_TYPE mmin=min(eigvals);
    if (mmin<=0.0)
    {      
      cerr << "Error -- hessian not positive definite in " 
       "lapack_choleski_inverse(M,ierr)" << endl
       << "calculating eigenvectors and eigenvalues " << endl;
    }	
    dvector dinv=1.0/eigvals;
    dvector vtmp;
    dmatrix tmp=make_dmatrix(1,nvar,1,nvar,vtmp);
    cout<<"starting to invert hessian nvar = " << nvar <<endl;
    for (int i=1;i<=nvar;i++)
    {
      for (int j=1;j<=nvar;j++)
      {
	tmp(i,j)=eigvecs(j,i)*dinv(i);
      }
    }
    // This puts inverse into vHessinv
    cblas_matrix_prod(nvar,veigvecs,vtmp,vHessinv);
    if (test)
    {
      ofstream ofs("hesstest");
      ofs << Hess*Hessinv << endl;
    }

    if (mmin<=0.0)
    {      
      ofstream ofs3("neigenvalues");
      ofstream ofs4("sorted eigenvectors");
      ofstream ofs4x("xsorted eigenvectors");
      int izero=0;
      for (int i=1;i<=nvar;i++)
      {
	if (eigvals(i)<=0.) izero++;
      }
      cout << " There are " << izero << " nonpositive eigenvalues" << endl;
      ofs3 << izero << " " << nvar << endl;
      izero=0;
      dmatrix unsorted_eigen_vectors(1,nvar,0,nvar);
       dmatrix bin(1,2,1,nvar);
       bin(2)=eigvals;
       bin.fill_seqadd(1,1);
       dmatrix bb=sort(trans(bin),2);
       //cout << column(bb,1) << endl << endl;
       for (int i=1;i<=nvar;i++)
       {
	 unsorted_eigen_vectors(i)(1,nvar)=eigvecs(i);
	 unsorted_eigen_vectors(i,0)=eigvals(i);
       }
       dmatrix sorted_eigen_vectors=sort(unsorted_eigen_vectors,0);
       for (int i=nvar;i>=1;i--)
       {
	 ofs4 << sorted_eigen_vectors(i,0) << "    ";
	 ofs4 << sorted_eigen_vectors(i)(1,nvar) << endl;
	 ofs4x << setfixed() << setprecision(5) << sorted_eigen_vectors(i,0)
	       << "     ";
	 ofs4x << setfixed() << setprecision(2)
	       << sorted_eigen_vectors(i)(1,nvar) << endl;
       }
       {
	 ofstream ofs("new_cor_report");
	 MY_DOUBLE_TYPE cut=0.0;
	 ofs << " Smallest eigenvalues" << endl;
	 int nrep=min(nvar/2,30);
	 for(int ind=1;ind<=nrep;ind++)
	 {
	   //cut=max(0.8*max(fabs(sorted_eigen_vectors(ind))),0.1);
	   //cut=min(0.4,cut);
#if !defined(NO_MY_DOUBLE_TYPE)
           float128 p1=0.1;
#else
           double p1=0.1;
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
           float128 c1=0.4;
#else
           double c1=0.4;
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
	   float128 d=0.8*max(fabs(sorted_eigen_vectors(ind)));
#else
	   double d=0.8*max(fabs(sorted_eigen_vectors(ind)));
#endif
	   cut=min(max(d,p1),c1);
	   ofs << setprecision(5) << setw(12) 
	       << sorted_eigen_vectors(ind,0) << "   ";
	   int ip=0;
	   cut*=1.4;
	   do
	   {
	     cut/=1.4;
	     for (int i=1;i<=nvar;i++)
	     {
	       if (fabs(sorted_eigen_vectors(ind,i))>cut)
	       {
		 ip++;
	       }
	     }
	   }
	   while (ip==0 || (ip ==1 && cut> 0.08) || (ip==2 && cut>0.12));
	   for (int i=1;i<=nvar;i++)
	   {
	     if (fabs(sorted_eigen_vectors(ind,i))>cut)
	     {
	       ofs << "(" << i << " " << setprecision(3) 
		   << sorted_eigen_vectors(ind,i) << ") ";
	     }
	   }
	   ofs << endl;
	 }
	 ofs << endl << " Largest eigenvalues" << endl;
	 for(int ind=nvar;ind>nvar-nrep;ind--)
	 {
#if !defined(NO_MY_DOUBLE_TYPE)
	   //cut=max(0.8*max(fabs(sorted_eigen_vectors(ind))),0.1L);
#else
	   //cut=max(0.8*max(fabs(sorted_eigen_vectors(ind))),0.1);
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
	   //cut=min(0.4L,cut);
#else
	   //cut=min(0.4,cut);
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
           float128 p1=0.1;
#else
           double p1=0.1;
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
           float128 c1=0.4;
#else
           double c1=0.4;
#endif
#if !defined(NO_MY_DOUBLE_TYPE)
	   float128 d=0.8*max(fabs(sorted_eigen_vectors(ind)));
#else
	   double d=0.8*max(fabs(sorted_eigen_vectors(ind)));
#endif
	   cut=min(max(d,p1),c1);
	   ofs << setprecision(5) << setw(12) 
	       << sorted_eigen_vectors(ind,0) << "   ";
	   int ip=0;
	   cut*=1.4;
	   do
	   {
	     cut/=1.4;
	     for (int i=1;i<=nvar;i++)
	     {
	       if (fabs(sorted_eigen_vectors(ind,i))>cut)
	       {
		 ip++;
	       }
	     }
	   }
	   while (ip==0 || (ip ==1 && cut> 0.08) || (ip==2 && cut>0.12));
	   for (int i=1;i<=nvar;i++)
	   {
	     if (fabs(sorted_eigen_vectors(ind,i))>cut)
	     {
	       ofs << "(" << i << " " << setprecision(3) 
		   << sorted_eigen_vectors(ind,i) << ") ";
	     }
	   }
	   ofs << endl;
	 }
       }

       for (int i=1;i<=nvar;i++)
       {
	 if (eigvals(i)<=0.)
	 {
	   ofs3 << setprecision(15) << setscientific() << eigvals(i) 
		<< "  "  << eigvecs(i) << endl;
	 }
       }
    }
  }   // -skips Hessian calcs needed for dependent variables
  fname=ad_root+".dep";
  uistream ifs1(fname);
  if (!ifs1)
  {
    cerr << "Error trying to open file " << fname << endl;
  }
  ifs1 >> nvar1 >> num_deps;
  if (!ifs1)
  {
    cerr << "Error reading nvar from depgrad.rpt" << endl;
  }

  if (fsh.parest_flags(237)==1)  //NMD16Mar2018 skip *.hes input unless needed
  {
    if (nvar!=nvar1)
    {
      cerr << "number of independent variables is not consistent"
  	    " between files" << "Hessinv.rpt and depgrad.rpt" << endl;
      ad_exit(1);
    }
  }
  else
  {
    nvar=nvar1;  //NMD20Mar2018 - set nvar to nvar1
  }
  dmatrix depgrad(1,num_deps,1,nvar);
  dmatrix covar(1,num_deps,1,num_deps);
  dvector stddev(1,num_deps);
  fname=ad_root+".cov";
  ofstream ofs(fname);
  ifstream ifs2("deplabel.tmp");
  adstring_array deplabel(1,num_deps);
  dvector depvalue(1,num_deps);
  depvalue.initialize();

  if (fsh.parest_flags(237)==1)  //NMD26Aug2015
  {
    for (i=1;i<=num_deps;i++)
    {
      ifs1 >> depgrad(i);
      ifs2 >> depvalue(i);
      ifs2 >> deplabel(i);
      if (!ifs1)
      {
        cerr << "Error reading gradient for dependent variable " << i
	         << "  from depgrad.rpt" << endl;
        ad_exit(1);
      }
    }
    for (i=1;i<=num_deps;i++)
    {
      dvector tmp=Hessinv*depgrad(i);
      for (j=1;j<=i;j++)
      {
        covar(i,j)=depgrad(j)*tmp;
        covar(j,i)=covar(i,j);
      }
      if (covar(i,i)>=0)
      {
        stddev(i)=sqrt(covar(i,i));
      }
      else
      {
        cerr << "non positive variance for dependent variable " << i 
             << " value is " << covar(i,i) << endl;
        cerr << "No covariance report will be produced" << endl;
        return;
      }
    }
  }
  for (i=1;i<=num_deps;i++)
  {
    //ofs << "                  "  << setw(9) << i << "  " ;
  }
//  - might need this section for nage
  int tndeps;
  uistream uis("proparams.val");
  uis >> tndeps;
  if (num_deps != tndeps)
  {
    cerr << "Bad value for the number of dependent variables"
            " in proparams.val" << endl;
    ad_exit(1);
  }
  dvector params(1,num_deps);
  uis>>params;
  if (!uis)
  {
    cerr << "Error reading from file proparams.val" << endl;
    ad_exit(1);
  }
//NMD 22Feb2012
  dvector params_noeff(1,num_deps);  //NMD 22Feb2012
  if(noeff_sw == 1)
  {
    uistream uis("proparams_noeff.val");
    if(!uis)
    {
      cerr << "Error trying to open file: proparams_noeff.val " << endl;
      cerr << "Need to run Hessian option 7 with no-effort option" << endl;
      ad_exit(1);
    }
    uis >> tndeps;
    if (num_deps != tndeps)
    {
      cerr << "Bad value for the number of dependent variables"
              " in proparams_noeff.val" << endl;
      ad_exit(1);
    }
    uis>>params_noeff;
  }
  //NMD 22Feb2012

  dmatrix chk;    //NMD26Aug2015
  dvector va;     //NMD26Aug2015
  if (fsh.parest_flags(237)==1)  //NMD26Aug2015
  {
    int ierr=0;
    do
    {
      chk=choleski_decomp_error(covar,ierr);   //NMD26Aug2015
      if (ierr)
      {
        cerr << "corr matrix not pos def -- tough luck"
             << endl;
        for (int i=1;i<=num_deps;i++)
        {
          for (int j=1;j<=num_deps;j++)
          {
            if (i!=j) 
              covar(i,j)*=0.99;
          }
        }
      }
    }
    while (ierr>0);
    chk=choleski_decomp_error(covar,ierr);    //NMD26Aug2015
    va=diagonal(covar);                       //NMD26Aug2015
  }

  int num_sims=fsh.age_flags(20);
  const int& num_regions=fsh.num_regions;
  const int& nage=fsh.nage;
  const int& last_real_year=fsh.last_real_year;
  const int& nyears=fsh.nyears;
  const ivector& num_real_fish_periods=fsh.num_real_fish_periods;
  d3_array numbers_at_age;
  d3_array numbers_at_age_noeff;  //NMD 22Feb2012
  //NMD_26Aug2015
  int ns=fsh.age_flags(57);
  int rem=fsh.last_real_year%ns;
  int na=(fsh.last_real_year-1)/ns+1;
  int na_proj=nyears/ns;

  // check to see if first projection period is in a new year
  if (fsh.year(1,num_real_fish_periods(1)+1)>
    fsh.year(1,num_real_fish_periods(1)))
  {
    numbers_at_age.allocate(1,num_sims,1,num_regions,2,nage);
    numbers_at_age_noeff.allocate(1,num_sims,1,num_regions,2,nage);  //NMD 22Feb2012
  }
  else
  {
    cerr << "Don't have code to deal with this yet" << endl;
    ad_exit(1);
    numbers_at_age.allocate(1,num_sims,1,num_regions,1,nage);
    numbers_at_age_noeff.allocate(1,num_sims,1,num_regions,2,nage);  //NMD 22Feb2012
  }
  


  d3_array recruitments(1,num_sims,1,num_regions,last_real_year+1,nyears);
  int iseed=1801;
  if (fsh.parest_flags(231)>0)
    iseed=fsh.parest_flags(231);
  random_number_generator rng(iseed);
  dmatrix eps;
  int pf237=fsh.parest_flags(237);
  if (pf237)
  {
    eps.allocate(1,num_sims,1,num_deps);
    eps.fill_randn(rng);
  }
  for (i=1;i<=num_sims;i++)
  {
    dvector tmp;
    dvector tmp_noeff;  //NMD 22Feb2012
    if (pf237)
    {
      tmp=exp(params+chk*eps(i)-0.5*va);
      if(noeff_sw == 1) tmp_noeff=exp(params_noeff+chk*eps(i)-0.5*va);  //NMD 22Feb2012
    }
    else
    {
      tmp=exp(params);
      if(noeff_sw == 1) tmp_noeff=exp(params_noeff);  //NMD 22Feb2012
    }
    int offset=1;
    for (int ir=1;ir<=num_regions;ir++)
    {
      numbers_at_age(i,ir)=tmp(offset,offset+nage-2).shift(2);
      if(noeff_sw == 1) numbers_at_age_noeff(i,ir)=tmp_noeff(offset,offset+nage-2).shift(2);  //NMD 22Feb2012
      offset+=nage-1;
    }
  }
  // now generate the recruitment by resampling 
  // historical recuitment
  int pf232=fsh.parest_flags(232);
  int pf233=fsh.parest_flags(233);
  if (pf232==0) pf232=1;
  if (pf232==-1) pf232=last_real_year;
  if (pf233==0) pf233=last_real_year;

  if (fsh.age_flags(182)==1 && fsh.age_flags(57)>1)  //NMD_26Aug2015
  {
    pf232=fsh.parest_flags(232)/fsh.age_flags(57)+1;
    pf233=(fsh.parest_flags(233)-1)/fsh.age_flags(57)+1;
  }

  int nslots=pf233-pf232+1;
  dvector p(1,nslots);
  p=1.0;
  p/=sum(p);

//NMD_26Aug2015
  int year1;
  int yearlst;
  if (fsh.age_flags(182)==1 && fsh.age_flags(57)>1)  //NMD_26Aug2015
  {
    year1=na+1;
    yearlst=na_proj;
  } 
  else
  {
    year1=last_real_year+1;
    yearlst=nyears;

  }
  const imatrix & _simyears=fsh.simyears;
  ADUNCONST(imatrix,simyears)
  if (allocated(simyears))
    simyears.deallocate();
  simyears.allocate(1,num_sims,year1,yearlst);
//NMD_26Aug2015


  for (i=1;i<=num_sims;i++)
  {
    simyears(i).fill_multinomial(rng,p);
    simyears(i)+=pf232-1;
  }
  {
    ofstream ofs("simyears");
    ofs << num_sims << " " << year1 << " " << yearlst << endl;
    ofs << simyears << endl;
  }
  
  dmatrix N(1,num_regions,1,last_real_year);
  {
    ifstream ifs("histrec"); 
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int iy=1;iy<=last_real_year;iy++)
      {
        ifs >> N(ir,iy);
      }
    }
    //cout << max(N) << endl;
    if (!ifs)
    {
      cerr << "Error reading from file histrec" << endl;
      ad_exit(1);
    }
  }
  //NMD_26Aug2015
  dmatrix Nann(1,num_regions,1,na);
  if (fsh.age_flags(182)==1 && fsh.age_flags(57)>1)  //NMD_26Aug2015
  {
    for (int ir=1;ir<=num_regions;ir++)
    {
      Nann(ir).initialize();
      for (int iy=1;iy<=fsh.last_real_year;iy++)
      {
        int ia=(iy-1)/ns+1;
        Nann(ir,ia)+=exp(N(ir,iy));
      }
    }
    for (int ir=1;ir<=num_regions;ir++)
    {
      for (int ia=1;ia<=na;ia++)
      {
        Nann(ir,ia)=log(Nann(ir,ia));
      }
    }
  }     //NMD_26Aug2015

  int pf238=fsh.parest_flags(238);
  if (pf238==0) 
  {
    for (int isim=1;isim<=num_sims;isim++)
    {
      for (int ir=1;ir<=num_regions;ir++)
      {
        if (fsh.age_flags(182)==1 && fsh.age_flags(57)>1)  //NMD_26Aug2015
        {
          for (int iy=last_real_year+1;iy<=nyears;iy++)
          {
            int ia=(iy-1)/ns+1;
            recruitments(isim,ir,iy)=Nann(ir,simyears(isim,ia));
          }
        }
        else
        {
          for (int iy=last_real_year+1;iy<=nyears;iy++)
            recruitments(isim,ir,iy)=N(ir,simyears(isim,iy));
        }  //NMD_27Aug2015
      }
    }
  }
  else
  {
    MY_DOUBLE_TYPE alpha=value(fsh.alpha);
    MY_DOUBLE_TYPE beta=value(fsh.beta);
    MY_DOUBLE_TYPE phi=value(fsh.phi);
    MY_DOUBLE_TYPE steepness=value(fsh.steepness);
  }
   
  ofstream ofs2("simulated_numbers_at_age");
  ofs2 <<"# simulated_numbers_at_age"<< endl;
  ofs2 << num_sims << endl;
  ofs2 << setprecision(3) << setfixed() << numbers_at_age  << endl;
  ofs2 <<"# simulated_recruitments"<< endl;
  ofs2 << recruitments << endl;
//NMD 22Feb2012
  if(noeff_sw == 1)
  {
    ofstream ofs3("simulated_numbers_at_age_noeff");
    ofs3 <<"# simulated_numbers_at_age"<< endl;
    ofs3 << num_sims << endl;
    ofs3 << setprecision(3) << setfixed() << numbers_at_age_noeff  << endl;
    ofs3 <<"# simulated_recruitments"<< endl;
    ofs3 << recruitments << endl;
  }
//NMD 22Feb2012

  ofs << "# number of dependent variables "  << endl;
  ofs << num_deps << endl;
  for (i=1;i<=num_deps;i++)
  {
    cout<<"calc covariance of dep vars. Got to "<<i<<endl;
    //ofs << "# " << setw(3) << i; 
    //ofs << "  "  << setw(9) << setprecision(4) 
     //    << depvalue(i) << "  "  << deplabel(i) << endl;
    for (j=1;j<=i;j++)
    {
      MY_DOUBLE_TYPE tmp=stddev(i)*stddev(j);
      if (tmp)
        ofs << "  "  << setw(8) << setprecision(3) << covar(i,j);
      else
        ofs << "  "  << setw(8) << setprecision(3) << tmp;
    }
    ofs << endl;
  }
}
