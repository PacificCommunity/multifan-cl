/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

void dvar_len_fish_stock_history::all_age_length_calcs(void)
{
  for (int ir=1;ir<=num_regions;ir++)
  {
    for (int ip=1;ip<=num_real_fish_periods(ir);ip++)  
    {
      for (int fi=1;fi<=num_fish_incidents(ir,ip);fi++) 
      {                                         
        if (age_len_sample_size(ir,ip,fi)>0)
        {
//          cout << age_len_sample_size(ir,ip,fi) << endl;
          int pi=parent(ir,ip,fi);
          dvar_matrix& q=qij(ir,ip,fi);
          if (!allocated(q))
          {
            q.allocate(1,nlint,1,nage);
          }
          dvar_vector& mean_len=mean_length(ir,ip,fi);
          dvar_vector& sigg=sdevs(ir,ip,fi);
          age_length_calcs(q,mean_len,sigg);
//          cout << q(1) << endl;
        }
      }
    }
  }
}

void dvar_len_fish_stock_history::age_length_calcs(dvar_matrix& q,
  dvar_vector& mean_len,dvar_vector& sigg)
{
  const MY_DOUBLE_TYPE a2 = 35.91908881;
  const MY_DOUBLE_TYPE a3 = 785.858644;
  const MY_DOUBLE_TYPE a22 = 2*35.91908881;
  const MY_DOUBLE_TYPE a33 = 3*785.858644;
  int js=1;
  int break_flag=0;
  int i;
  q.initialize();
  for (i=1; i<=nlint; i++) 
  {
    MY_DOUBLE_TYPE& fmidd=fmid(i);
   
    for (int j=1; j<=nage; j++) 
    {
      dvariable t=(fmidd-mean_len(j))/sigg(j);
      if ( fabs(t) <= 3.03)
      {
        if ( fabs(t) <= 3.00)
          q(i,j)=exp(-t*t/2.e0)/sigg(j);
        else if (t > 0.e0)
        {
          dvariable u=t-3.03;
          q(i,j)=(u*u*(a2+a3*u))/sigg(j);
        }
        else if (t < 0.e0)
        {
          dvariable u=t+3.03;
          q(i,j)=(u*u*(a2-a3*u))/sigg(j);
        }
      }
    }
  }
}
cifstream & operator >> (cifstream & cif, age_length_record& alr)
{
  cif >> alr.year >> alr.month >> alr.fishery >> alr.species;
  cif >> alr.age_length;
  return cif;
}

age_length_record_array::age_length_record_array(const age_length_record_array & alarray)
{
  mmin=alarray.mmin;
  mmax=alarray.mmax;
  ptr=alarray.ptr;
  ncopies=alarray.ncopies;
  if (ncopies==0)
  {
    cerr << "making a copy of an unallocated age_length_record_array" << endl;
  }
  else
  {
    *ncopies++;
  }
}
  
void age_length_record_array::allocate(int _mmin,int _mmax,
  int nlint,int nage)  
{
  mmin=_mmin;
  mmax=_mmax;
  int sz=mmax-mmin+1;
  sumwght=0.0;
  effective_sample_size.allocate(mmin,mmax);
  ptr=new age_length_record[sz];
  ptr-=mmin;
  ncopies=new int;
  obs_sample_size_by_length.allocate(mmin,mmax,1,nlint);
  *ncopies=0;
  for (int i=mmin;i<=mmax;i++)
  {
    ptr[i].allocate(nlint,nage);
  }
}

age_length_record_array::age_length_record_array(void)
{
  mmin=0;
  mmax=-1;
  ptr=0;
  ncopies=0;
}

void dvar_len_fish_stock_history::get_age_length_periods(void)
{
  int mmin=alarray.indexmin();
  int mmax=alarray.indexmax();

  for (int i=mmin;i<=mmax;i++)
  {
    int fi=alarray(i).fishery;
    int year=alarray(i).year;
    int month=alarray(i).month;
    int no_match_flag=1;
    for (int nt=1;nt<=num_fish_times(fi);nt++)
    {
      int rr=realization_region(fi,nt);
      int rp=realization_period(fi,nt);
      int ty=really_true_year(rr,rp)+year1-1;
      int tm=really_true_month(rr,rp);
      int inc=realization_incident(fi,nt);
      if (ty==year && tm==month) 
      {
        alarray(i).region=rr;
        alarray(i).period=rp;
        alarray(i).incident=inc;
        age_len_sample_size(rr,rp,inc)=sum(alarray(i).age_length);
        no_match_flag=0;
      }
    }
    if (no_match_flag==1)
    {
      cerr << "Error -- age-length data year and/or month for sample "
           << i << "  do not match  any fishing incident " << endl;
      ad_exit(1);
    }
  }
}

dvariable dvar_len_fish_stock_history::fit_age_length_data(void)
{
  int mmin=alarray.indexmin();
  int mmax=alarray.indexmax();
 
  dvariable loglik=0.0;
  ofstream * pofs=0;
  if (generate_report)
    pofs =new ofstream("agelengthresids.dat");
  ppstf->age_length_like.initialize();
  MY_DOUBLE_TYPE totclike=0.0;
  for (int i=mmin;i<=mmax;i++)
  {
    MY_DOUBLE_TYPE clike=0.0;
    dmatrix & obs_sample_size_by_length=alarray.obs_sample_size_by_length;
    dvector & effective_sample_size=alarray.effective_sample_size;
    int year=alarray(i).year;
    int fishery=alarray(i).fishery;
    int month=alarray(i).month;
    int rr=alarray(i).region;
    int rp=alarray(i).period;
    int inc=alarray(i).incident;

    if (age_len_sample_size(rr,rp,inc)>0)
    {
      dvar_matrix& q=qij(rr,rp,inc);
      dvar_vector& p=prop(rr,rp,inc);
      int lmin=q.indexmin();
      int lmax=q.indexmax();
      dmatrix& age_length=alarray(i).age_length;
      dmatrix res;
      if (generate_report)
        res.allocate(lmin,lmax);
      if (generate_report)
      {
        (*pofs) << "Fishery " << fishery << "  Year " << year 
                << "  Month " << month << endl;
      }
      for (int ii=lmin;ii<=lmax;ii++)
      {
        dvar_vector r=elem_prod(prop(rr,rp,inc),q(ii));
        r/=(1.e-10+sum(r));  // condition on length interval ii
        dvariable ll=0.0;
        if (alarray.sumwght)
        {
          ll=effective_sample_size(i)*obs_sample_size_by_length(i,ii)
               *(age_length(ii)*log(1.e-8+r)); 
        }      
        else
        {
          ll=age_length(ii)*log(1.e-8+r);   
        }    
        loglik+=ll;       
        clike+=value(ll);
        // for age length fit report
        if (generate_report)
        {
          dvector obsp=age_length(ii)/(1.e-10+sum(age_length(ii)));
          dvector cr=value(r);
          res(ii)=elem_div(obsp-cr,sqrt(1.e-20+elem_prod(cr,1.0-cr)));
          (*pofs) << sum(age_length(ii)) << endl;
          (*pofs) << obsp << endl;
          (*pofs) << cr << endl;
        }
      }
    }
    ppstf->age_length_like(i)=-1.0*clike; //NMD_25May2015
    totclike+=clike;
  }
  if (pofs)
  {
    delete pofs;
    pofs=0;
  }
  cout << "Total age length data " << -1.0*totclike << endl; //NMD_25May2015
  return -1.0*loglik;  // return negatvie log-likelihood
}

void read_age_length_records(cifstream & cif,int nlint,int nage,
  dvar_len_fish_stock_history& fsh)
{
  int num_age_length_records;
  cif >> num_age_length_records;
  fsh.alarray.allocate(1,num_age_length_records,nlint,nage);
  cif >> fsh.alarray.effective_sample_size;
  fsh.alarray.sumwght=sum(fsh.alarray.effective_sample_size);
  for (int i=1;i<=num_age_length_records;i++)
  {
    cif >> fsh.alarray.ptr[i];
  }

  if (fsh.pmsd)
  {
    for (int i=1;i<=num_age_length_records;i++)
    {
//      cout << "Record: " << i << " " << "Fshry: " << fsh.alarray(i).fishery
//           << "   Species: " << fsh.alarray(i).species << endl;
      if (fsh.alarray(i).species > 1)
      {
	int isp=fsh.alarray(i).species;
        int offset=(isp-1)*fsh.pmsd->num_real_fisheries;
        fsh.alarray(i).fishery=fsh.alarray(i).fishery+offset;
      }
//      cout << "Record: " << i << " " << "Fshry: " << fsh.alarray(i).fishery
//           << "   Species: " << fsh.alarray(i).species << endl;
    }
  }
  
  int lmin=fsh.alarray(1).age_length.indexmin();
  int lmax=fsh.alarray(1).age_length.indexmax();
  for (int i=1;i<=num_age_length_records;i++)
  {
    for (int ii=lmin;ii<=lmax;ii++)
    {
      fsh.alarray.obs_sample_size_by_length(i,ii)
        =sum(fsh.alarray(i).age_length(ii));
    }
  }
  if (fsh.alarray.sumwght)
  {
    int lmin=fsh.alarray(1).age_length.indexmin();
    int lmax=fsh.alarray(1).age_length.indexmax();
    for (int i=1;i<=num_age_length_records;i++)
    {
      for (int ii=lmin;ii<=lmax;ii++)
      {
        fsh.alarray(i).age_length(ii)
          /=(1.e-10+fsh.alarray.obs_sample_size_by_length(i,ii));
      }
    }
  }
  
  if (fsh.ppstf)
  {
    __SAFE_DEALLOCATE__(fsh.ppstf->age_length_like);
    fsh.ppstf->age_length_like.allocate(1,num_age_length_records);
  }
}

age_length_record_array::~age_length_record_array()
{
  if (ncopies)
  { 
    if (*ncopies)
    {
      *ncopies--;
    }
    else
    {
      ptr+=mmin;
      delete [] ptr;
      ptr=0;
      delete ncopies;
      ncopies=0;
      mmin=0;
      mmax=-1;
    }
  }
}

