/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#define USE_NO_LAPACKE
#pragma implementation "newmult.hpp"
#pragma implementation "variable.hpp"
#include "version.h"
#include "all.hpp"
#define _INC_ERRNO
#if defined(__MSVC32__)
//#include <error.h>
//  int _errno=0;
extern "C" {
 //int * errno(){ return 0;}
}

extern int need_this_int;   // to ensure this gets linked in

//using namespace std;   //NMD_23Mar2015

#endif
#if !defined(__MINGW__ ) && !defined(_WIN32) && !defined(_MAC)
#  include <fpu_control.h>
#endif
#include <fenv.h>

//#define USE_DD_NOT
#if defined(USE_DD)
#  include <qd/fpu.h>
#endif

#include "safe_mem.h"
//#include <ctrl87.h>
#include <stdio.h>

#if defined(WIN32)
#  include <windows.h>
#endif
#if defined(USE_ADPVM)
extern "C" {
#include "pvm3.h"
}
#include <adpvm2.h>
#endif //#if defined(USE_ADPVM)

//#define  __declspec(dllexport) 
/*
 std::basic_ostream<char>& operator << (std::basic_ostream<char>& s,
    _Quad& x);
  
  std::basic_istream<char>& operator >> (std::basic_istream<char>& s,
    _Quad& x);
*/
  
void print_flfra_for_debug(fishery_freq_record_array&  flfra, int I)
{
  ofstream ofs("flfra_debug" +str(I));
  int mmin=flfra.indexmin();
  int mmax=flfra.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    ofs << flfra[i] << endl;
  }
}


extern "C" {
  #include <signal.h>

#if !defined(__MINGW__)  && !defined(_WIN64)
static void handle_fpe_signal (int sig, siginfo_t *sip, ucontext_t *uap) 
{
  signal(SIGFPE, handle_fpe_signal);
  ad_exit(1); 
} 
#endif
}

int TURNOFF=1;
adtimer * profile_timer_1=0;
adtimer * profile_timer_2=0;

void print (fishery_freq_record_array&  flfra)
{
  ofstream ofs("ffiles");
  int mmin=flfra.indexmin();
  int mmax=flfra.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    ofs << flfra(i) << endl;
  }
}

void dvar_fish_stock_history::calculate_survey_indices(void)
{
  ivector ff92=column(fish_flags,92);

  first_survey_time.initialize();
  for (int fi=1;fi<=num_fisheries;fi++)
  {
    if (ff92(fi))
    {
      for (int i=1;i<=num_fish_times(fi);i++)
      {
        int rr=realization_region(fi,i);
        int rp=realization_period(fi,i);
        int ri=realization_incident(fi,i);
        if (missing_effort_by_realization_flag(fi,i))
        {
          first_survey_time(fi)=i;
          break;
        }
      }
    }
  }
  ivector nft;    //NMD_28jun2022
  if (sum(data_fish_flags(2)))
  {
    nft=num_real_fish_times;
  }
  else
  {
    nft=num_fish_times;
  }
  if (!allocated(survey_index))
  {
//     survey_index.allocate(1,num_fisheries,first_survey_time,num_fish_times);
     survey_index.allocate(1,num_fisheries,first_survey_time,nft);
  }
  if (!allocated(survey_cpue_obs))   //NMD_8mar2022
  {
//     survey_cpue_obs.allocate(1,num_fisheries,first_survey_time,
//       num_fish_times);
     survey_cpue_obs.allocate(1,num_fisheries,first_survey_time,
       nft);
  }
  survey_cpue_obs.initialize();
  if (!allocated(survey_cpue_pred))   //NMD_8mar2022
  {
//     survey_cpue_pred.allocate(1,num_fisheries,first_survey_time,
//       num_fish_times);
     survey_cpue_pred.allocate(1,num_fisheries,first_survey_time,
       nft);
  }
  survey_cpue_pred.initialize();

  for (int fi=1;fi<=num_fisheries;fi++)
  {
    if (ff92(fi))
    {
//      for (int i=first_survey_time(fi);i<=num_fish_times(fi);i++)
      for (int i=first_survey_time(fi);i<=nft(fi);i++)
      {
        int rr=realization_region(fi,i);
        int rp=realization_period(fi,i);
        int ri=realization_incident(fi,i);
        if (missing_effort_by_realization_flag(fi,i))
        {
          survey_index(fi,i)=obs_tot_catch(rr,rp,ri)/
            really_true_effort(rr,rp,ri);
//            effort(rr,rp,ri);   //NMD_1dec2021
        }
        else
        {
          cerr << "Can not have missing effort in survey fishery"
             " indices " << endl;
          cerr << "Missing effort for fishery " << fi << " fish time " 
            << i << " I give up " << endl;
          ad_exit(1);
        }
      }
    }
  }
}



void set_option_flag(const char * s,int& ss,int argc,char * argv[])
{
  int nopt=0;
  int on1=0;
  if ( (on1=option_match(argc,argv,s,nopt))>-1)
  {
    if (!nopt)
    {
      cerr << "Usage " << s << " option needs integer  -- ignored" << endl;
    }
    else
    {   
      int jj=atoi(argv[on1+1]);
      ss=jj;
    }
  }
}

extern "C" {
 void mf_system(const char * s) { system(s);}
 //void mf_system(const char *){ cerr << " Not implemented " << endl; 
 // ad_exit(1);}
 }

//grad_stack_entry * gse= (grad_stack_entry*)(284164092);

gradient_structure * pgs = 0;
int gss = 0;
// shape printing
void ms(ivector & v,int i,int d)
{
  v(i)=d;  cout << i << " " << v(i) << endl;
}

void ms(dvector & v,int i,MY_DOUBLE_TYPE d)
{
  v(i)=d; cout << i << " " << v(i) << endl;
}

void ms(dvar_vector & v,int i,MY_DOUBLE_TYPE d)
{
  v(i)=d; cout << i << " " << v(i) << endl;
}

/*
void mp(const fishery_header_record * fhr,int i)
{
  cout << fhr[i] << endl;
}
*/

void mp(const banded_lower_triangular_dmatrix & v)
{
  cout << make_dmatrix(v) << endl;
}

void mp(const banded_lower_triangular_dvar_matrix & v)
{
  cout << make_dvar_matrix(v) << endl;
}

void mp(const dvar5_array & v)
{
  cout << v << endl;
}

void mpc(const imatrix & v,int i)
{
  ivector iv=column(v,i);
  cout << iv << endl;
}

void mpt(const dmatrix & v,int i)
{
  cout << (trans(v))(i) << endl;
}

void mpt(const dvar_matrix & v,int i)
{
  cout << (trans(value(v)))(i) << endl;
}

void mpinit(const dvar_matrix & _v,int i)
{
  ADUNCONST(dvar_matrix,v)
  v(i).initialize();
}
void mpset(dvar_matrix & v,int i,MY_DOUBLE_TYPE x)
{
  v(i)=x;
}
void mpset(dvar_matrix & v,int i,int j,MY_DOUBLE_TYPE x)
{
  v(i,j)=x;
  cout << " v(i,j) = " << x << endl;
}
void mpset(dvar_matrix & v,MY_DOUBLE_TYPE x)
{
  v=x;
}

void mpinit(const dmatrix & _v,int i)
{
  ADUNCONST(dmatrix,v)
  v(i).initialize();
}
void mpset(const dmatrix & _v,int i,MY_DOUBLE_TYPE x)
{
  ADUNCONST(dmatrix,v)
  v(i)=x;
}
void mpset(const dmatrix & _v,MY_DOUBLE_TYPE x)
{
  ADUNCONST(dmatrix,v)
  v=x;
}

void mpc(const dmatrix & v,int i)
{
  dvector iv=column(v,i);
  cout << iv << endl;
}

void mpc(const dvar_matrix & v,int i)
{
  dvar_vector iv=column(v,i);
  cout << iv << endl;
}

void mp(const dvar5_array & v,int i)
{
  cout << v(i) << endl;
}
void mp(const dvar5_array & v,int i,int j)
{
  cout << v(i,j) << endl;
}
void mp(const dvar5_array & v,int i,int j,int k)
{
  cout << v(i,j,k) << endl;
}
void mp(const dvar5_array & v,int i,int j,int k,int l)
{
  cout << v(i,j,k,l) << endl;
}

void mp(const adstring_array & v,int i)
{
  cout << v(i) << endl;
}

void mp(const adstring & v)
{
  cout << v << endl;
}

void mp(const dvar3_array & v)
{
  cout << v << endl;
}
void mp(const dvar3_array & v,int i)
{
  cout << v(i) << endl;
}
void mp(const dvar3_array & v,int i,int j)
{
  cout << v(i,j) << endl;
}
void mp(const dvar4_array & v)
{
  cout << v << endl;
}
void mp(const dvar4_array & v,int i)
{
  cout << v(i) << endl;
}
void mp(const dvar4_array & v,int i,int j)
{
  cout << v(i,j) << endl;
}

void mpsum(const imatrix & v,int i)
{
  cout << sum(v(i)) << endl;
}

void mpsum(const imatrix & v)
{
  cout << sum(v) << endl;
}


void mpsumexp(const dvar_matrix & v,int i)
{
  cout << sum(exp(v(i))) << endl;
}

void mpsumexp(const dmatrix & v,int i)
{
  cout << sum(exp(v(i))) << endl;
}

void mpsum(const dmatrix & v,int i)
{
  cout << sum(v(i)) << endl;
}

void mpsum(const dvar_matrix & v,int i)
{
  cout << sum(v(i)) << endl;
}

void mpsumexp(const dvar3_array & v,int i,int j)
{
  cout << sum(exp(v(i,j))) << endl;
}

void mpsum(const dvar3_array & v,int i,int j)
{
  cout << sum(v(i,j)) << endl;
}

void mpsumexp(const dvar4_array & v,int i,int j,int k)
{
  cout << sum(exp(v(i,j,k))) << endl;
}

void mpsum(const dvar4_array & v,int i,int j,int k)
{
  cout << sum(v(i,j,k)) << endl;
}

void mpexp(const dvar4_array & v,int i,int j,int k)
{
  cout << exp(v(i,j,k)) << endl;
}

void mpexp(const dvar3_array & v,int i,int j)
{
  cout << exp(v(i,j)) << endl;
}

void mpexp(const dvar_matrix & v,int i)
{
  cout << exp(v(i)) << endl;
}


void mp(const d5_array & v,int i)
{
  cout << v(i) << endl;
}

void mp(const d5_array & v,int i,int j)
{
  cout << v(i,j) << endl;
}

void mp(const d5_array & v,int i,int j,int k)
{
  cout << v(i,j,k) << endl;
}

void mp(const d5_array & v,int i,int j,int k,int l)
{
  cout << v(i,j,k,l) << endl;
}


void mp(const dvar4_array & v,int i,int j,int k)
{
  cout << v(i,j,k) << endl;
}
void mp(const d4_array & v,int i,int j,int k)
{
  cout << v(i,j,k) << endl;
}
void mp(const d4_array & v)
{
  cout << v << endl;
}
void mp(const d4_array & v,int i)
{
  cout << v(i) << endl;
}
void mp(const d4_array & v,int i,int j)
{
  cout << v(i,j) << endl;
}
void mp(const dvar3_array & v,int i,int j,int k)
{
  cout << v(i,j,k) << endl;
}
void mp(const d3_array & v)
{
  cout << v << endl;
}
void mp(const d3_array & v,int i)
{
  cout << v(i) << endl;
}
void mp(const d3_array & v,int i,int j)
{
  cout << v(i,j) << endl;
}

void mp(const d3_array & v,int i,int j,int k)
{
  cout << v(i,j,k) << endl;
}


void mp(const i3_array & v)
{
  cout << v << endl;
}

void mp(const i3_array & v,int i)
{
  cout << v(i) << endl;
}
void mp(const i3_array & v,int i,int j)
{
  cout << v(i,j) << endl;
}

void mp(const i3_array & v,int i,int j,int k)
{
  cout << v(i,j,k) << endl;
}

void mp(const i4_array & v)
{
  cout << v << endl;
}

void mp(const i4_array & v,int i)
{
  cout << v(i) << endl;
}
void mp(const i4_array & v,int i,int j)
{
  cout << v(i,j) << endl;
}

void mp(const i4_array & v,int i,int j,int k)
{
  cout << v(i,j,k) << endl;
}
void mp(const i4_array & v,int i,int j,int k,int l)
{
  cout << v(i,j,k,l) << endl;
}

void mp(const dvar_vector & v,const char * s)
{
  if (!strcmp(s,"s"))
  {
    cout << sum(v) << endl;
  }
}

void emp(const dvar_vector & v )
{
  cout << exp(v) << endl;
}

void emp(const dvar_vector & v,int i)
{
  cout << exp(v[i]) << endl;
}

void emp(const dvector & v )
{
  cout << exp(v) << endl;
}

void emp(const dvector & v,int i)
{
  cout << exp(v[i]) << endl;
}

void emp(const dvar_matrix & v )
{
  cout << exp(v) << endl;
}

void emp(const dvar_matrix & v,int i)
{
  cout << exp(v[i]) << endl;
}

void emp(const dvar3_array & v,int i)
{
  cout << exp(v[i]) << endl;
}

void emp(const dvar3_array & v,int i,int j)
{
  cout << exp(v(i,j)) << endl;
}

void emp(const d3_array & v,int i)
{
  cout << exp(v[i]) << endl;
}

void emp(const d3_array & v,int i,int j)
{
  cout << exp(v(i,j)) << endl;
}

void emp(const dvar4_array & v,int i)
{
  cout << exp(v[i]) << endl;
}

void emp(const dvar4_array & v,int i,int j)
{
  cout << exp(v(i,j)) << endl;
}
void emp(const dvar4_array & v,int i,int j,int k)
{
  cout << exp(v(i,j,k)) << endl;
}

void emp(const d4_array & v,int i)
{
  cout << exp(v[i]) << endl;
}

void emp(const d4_array & v,int i,int j,int k)
{
  cout << exp(v(i,j,k)) << endl;
}
void emp(const d4_array & v,int i,int j)
{
  cout << exp(v(i,j)) << endl;
}

void emp(const dmatrix & v )
{
  cout << exp(v) << endl;
}

void semp(const dmatrix & v )
{
  cout << sum(exp(v)) << endl;
}

void semp(const dmatrix & v ,int i)
{
  cout << sum(exp(v(i))) << endl;
}

void semp(const dvar_matrix & v )
{
  cout << sum(exp(v)) << endl;
}

void semp(const dvar_matrix & v ,int i)
{
  cout << sum(exp(v(i))) << endl;
}

void emp(const dmatrix & v,int i)
{
  cout << exp(v[i]) << endl;
}


void mp(const dvar_vector & v,int i)
{
  cout << v[i] << endl;
}

void mp(const dvar_vector & _v,int i,int j)
{
  ADUNCONST(dvar_vector,v)
  cout << v(i,j) << endl;
}

void mp(const dvector & _v,int i,int j)
{
  ADUNCONST(dvector,v)
  cout << v(i,j) << endl;
}

void setmp(const int & _v,int x)
{
  ADUNCONST(int,v)
  v=x;
}

void setmp(const ivector & _v,int x,int i)
{
  ADUNCONST(ivector,v)
  if (i<=v.indexmax() && i>=v.indexmin())
    v(i)=x;
  else
    cout << "illegal index" << endl;
}

void setmp(const dvector & _v,MY_DOUBLE_TYPE x,int i)
{
  ADUNCONST(dvector,v)
  if (i<=v.indexmax() && i>=v.indexmin())
    v(i)=x;
  else
    cout << "illegal index" << endl;
}

void setmp(const dvar_vector & _v,MY_DOUBLE_TYPE x,int i)
{
  ADUNCONST(dvar_vector,v)
  if (i<=v.indexmax() && i>=v.indexmin())
    v(i)=x;
  else
    cout << "illegal index" << endl;
}



void mp(const dvector & _v,int i)
{
  ADUNCONST(dvector,v)
  cout << v(i) << endl;
}

void mps(const dvariable & _v)
{
  ADUNCONST(dvariable,v)
  cout << setscientific() << setprecision(5) << v << endl;
}

void mps(const dvar_vector & _v,int i)
{
  ADUNCONST(dvar_vector,v)
  cout << setscientific() << setprecision(5) << v[i] << endl;
}

void mps(const dvar_vector & _v)
{
  ADUNCONST(dvar_vector,v)
  cout << setscientific() << setprecision(5) << v << endl;
}

void mpss(const dvariable & _v)
{
  ADUNCONST(dvariable,v)
  cout << setscientific() << setprecision(10) << v << endl;
}

void mpss(const dvar_vector & _v,int i)
{
  ADUNCONST(dvar_vector,v)
  cout << setscientific() << setprecision(10) << v[i] << endl;
}

void mpss(const dvar_vector & _v)
{
  ADUNCONST(dvar_vector,v)
  cout << setscientific() << setprecision(10) << v << endl;
}

void mp(const dvar_vector & _v)
{
  ADUNCONST(dvar_vector,v)
  cout << v << endl;
}

void mpnz(const dvar_vector & _v)
{
  ADUNCONST(dvar_vector,v)
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  int print_flag=0;
  for (int i=mmin;i<=mmax;i++)
  {
    if (v[i] != 0.0)
    {
      print_flag=1;
      cout << setscientific() << v[i] << "  ";
    }
  }
  if (print_flag) cout << endl;
}

void mpnz(const dvar_matrix & _v)
{
  ADUNCONST(dvar_matrix,v)
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    mpnz(v[i]);
  }
}
 
void mpnz(const dvar3_array & _v)
{
  ADUNCONST(dvar3_array,v)
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    mpnz(v[i]);
  }
}
 
void mpnz(const dvector & _v)
{
  ADUNCONST(dvector,v)
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    if (v[i] != 0.0)
    {
      cout << setscientific() << v[i] << "  ";
    }
  }
  cout << endl;
}

void mpnz(const dmatrix & _v)
{
  ADUNCONST(dmatrix,v)
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    mpnz(v[i]);
  }
}
 
void mpnz(const d3_array & _v)
{
  ADUNCONST(d3_array,v)
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    mpnz(v[i]);
  }
}

void mpgt(const dvar_vector & _v)
{
  ADUNCONST(dvar_vector,v)
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    if (fabs(value(v[i])) > 1.e-50)
    {
      cout << setscientific() << v[i] << "  ";
    }
  }
  cout << endl;
}

void mpgt(const dvar_matrix & _v)
{
  ADUNCONST(dvar_matrix,v)
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    mpgt(v[i]);
  }
}
 
void mpgt(const dvar3_array & _v)
{
  ADUNCONST(dvar3_array,v)
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    mpgt(v[i]);
  }
}
 
void mpgt(const dvector & _v)
{
  ADUNCONST(dvector,v)
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  int print_flag=0;
  for (int i=mmin;i<=mmax;i++)
  {
    if (fabs(v[i]) > 1.e-50)
    {
      print_flag=1;
      cout << setscientific() << v[i] << "  ";
    }
  }
  if (print_flag) cout << endl;
}

void mpgt(const dmatrix & _v)
{
  ADUNCONST(dmatrix,v)
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    mpgt(v[i]);
  }
}
 
void mpgt(const d3_array & _v)
{
  ADUNCONST(d3_array,v)
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    mpgt(v[i]);
  }
}

void mp(MY_DOUBLE_TYPE v)
{
  cout << v << endl;
}

void mp(const prevariable & v)
{
  cout << v << endl;
}

void mp(const dvector & v)
{
  cout << v << endl;
}

void mp(const ivector & v,int i)
{
  if (i<=v.indexmax() && i>=v.indexmin())
    cout << v(i) << endl;
  else
    cout << "illegal index" << endl;
}

void mp(const ivector & v)
{
  cout << v << endl;
}

void mpp(const dvariable & v)
{
  cout << setprecision(12) << v << endl;
}

void mpp(const dmatrix & v,int i)
{
  cout << setprecision(12) << v(i) << endl;
}

void mpp(const dvar_matrix & v,int i)
{
  cout << setprecision(12) << v(i) << endl;
}

void mp(const dvar_matrix & v,int i)
{
  cout << v(i) << endl;
}

void mp(const dvar_matrix & v,int i,int j)
{
  cout << v(i,j) << endl;
}

void mp(const dvar_matrix & v)
{
  cout << v << endl;
}

void mp(const dmatrix & v)
{
  cout << v << endl;
}
void mp(const dmatrix & v,int i)
{
  cout << v(i) << endl;
}
void mp(const dmatrix & v,int i,int j)
{
  cout << v(i,j) << endl;
}


void mp(const imatrix & v)
{
  cout << v << endl;
}

void mp(const imatrix & v,int i)
{
  cout << v(i) << endl;
}
void mp(const imatrix & v,int i,int j)
{
  cout << v(i,j) << endl;
}

int getsize(const dvar_vector & v)
{
  if (allocated(v))
    return v.indexmax()-v.indexmin()+1;
  else
    return 0;
}


void sumprint(const dvar_vector & v)
{
  cout << sum(v) << endl;
}

void sumprint(const dvector & v)
{
  cout << sum(v) << endl;
}

void sumprint(const ivector & v)
{
  cout << sum(v) << endl;
}

void sumprint(const dvar_matrix & v,int i)
{
  cout << v(i) << endl;
}

void sumprint(const dvar_matrix & v,int i,int j)
{
  cout << v(i,j) << endl;
}

void sumprint(const dvar_matrix & v)
{
  cout << sum(v) << endl;
}

void sumprint(const dmatrix & v)
{
  cout << sum(v) << endl;
}

void mpa(const dmatrix& M,int i,int j)
{
  cout << &(M(i,j)) << endl;
}

void mpa(const dvector& M,int i)
{
  cout << &(M(i)) << endl;
}

void mpa(const d3_array& M,int i,int j,int k)
{
  cout << &(M(i,j,k)) << endl;
}

void mpa(const dvar3_array& M,int i,int j,int k)
{
  cout << &(value(M(i,j,k))) << endl;
}

void mpa(const dvar_matrix& M,int i,int j)
{
  cout << &(value(M(i,j))) << endl;
}

void mpa(const dvar_vector& M,int i)
{
  cout << &(value(M(i))) << endl;
}

void sumprint(const imatrix & v)
{
  cout << sum(v) << endl;
}
void mp(const dvar3_array & v,char * s,int j,int k)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    cout << v(i,j,k) << " ";
  }
  cout << endl;
}

void emp(const dvar3_array & v,char * s,int j,int k)
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    cout << exp(v(i,j,k)) << " ";
  }
  cout << endl;
}

  ofstream * clogf=0;
  int ad_argc=0;
  char ** ad_argv=NULL;
  extern unsigned int _stack=36000;
  int fishery_freq_record::real_nwint;
  int fishery_freq_record::wgroup_factor;
  int fishery_freq_record::real_nlint;
  int fishery_freq_record::age_nage;
  int fishery_freq_record::age_age1;
  int fishery_freq_record::group_factor;
  adstring ad_root;
  long int ad_array_mbl_size=0L;

#ifdef __BCPLUSPLUS__
  unsigned _stklen=36000;
#endif

#include "multjunk.h"
#include <float.h>

//ofstream DERCH_ofs("derch.rpt");
dvector mean_weights_kludge;
adstring * _current_path=NULL;
//int _sd_flag=0;

int _ini_version_no=0;        // acts as ini file version for building the
                            // par file vi -makepar option
int _data_file_version_no;  // Set this for new versions. This is read in
int _data_version_no=9;    // Set this for new versions
int _version_no=1067; //  This is the multifan-cl major version number of the par file.
//int _version_no=1059; //  This is the multifan-cl major version number.
//int _version_no=1056; //  This is the multifan-cl major version number.
int _file_version_no;    // Set this for new versions this is read in from the
                         // from the par file to test for incompatibility
                         
int _minor_version_number=1;  //  This is the multifan-cl minor version number.
                         // It is used to distinguish between changes that do not
                       // affect par file compatibility
int fish_nonmv_flag=0;  // for backward compatibility of fishing in non-move periods
int _NUMSV=30;
int makeparflag=0;
  
_mf_cl_globals mfglobals;

class multi_species_data;

typedef multi_species_data * pmulti_species_data;

void process_multi_species_data(pmulti_species_data & pmsd, 
  cifstream& infile,int num_species,int num_regions,
  ivector& fishery_regions);

void set_region_incidence_matrix(int nswitch,cifstream& infile,
  ivector& dataswitch,imatrix& Dflags,int& num_fisheries,int& num_regions,
  ivector& fishery_regions,int& tag_group_flag,dvector& _region_area,int & _data_file_version_no,pmulti_species_data &  pmsd);

int age_data_flag=0;
int sex_flag=0;
adstring movement_coffs_file;
adstring full_datafile_path;
adstring full_input_parfile_path;
adstring full_output_parfile_path;
adstring directory_path;
imatrix cl_switches;
void mfsend_x_to_slaves(const imatrix&  x);
imatrix mfget_imatrix_from_master(void);
void mfget_file_from_master(const adstring& s);
void mfsend_file_to_slaves(const adstring& s);

void read_option_file(const adstring& full_input_option_path,
  imatrix& cls);


void allocate_species(fishery_freq_record_array&  flfra,int num_species);

void do_grouped_between_time_calculations(dvar_fish_stock_history& fsh);

 void test_fpu(void)
 {
#if defined(USE_DD)
   unsigned int cw;
   fpu_fix_start(&cw);
   unsigned int cw2;
   fpu_fix_start(&cw2);
   
  /*
   //if (cw != cw2)
   //{
   //  ivector itmp(1,2);
   //  itmp(1)=cw;
   //  itmp(2)=cw2;
   //  cout << "fpu status changed" << endl;
   //  for (int i=1;i<=16;i++)
   //  {
   //    cout << setw(2) << i-1 << " ";
   //  }
   //  //cout <<"nnewlan.cpp " << endl;
   //  for (int ii=1;ii<=2;ii++)
   //  {
   //    int tcw=itmp[ii];
   //    unsigned int mask=1;
   //    for (int i=1;i<=16;i++)
   //    {
   //      if (mask & tcw)  
   //        cout << " 1 ";
   //      else
   //        cout << " 0 ";
   //      mask*=2;
   //      //cout << (mask & cw) << " " << (mask & cw2) << endl;
   //    }
   //    cout <<"nnewlan.cpp " << endl;
   //  }
   //} 
  */
#endif
 }

//extern ad_kill_flag=1;
#if !defined(__BORLANDC__)
  int heapcheck(void){return 0;}
#endif

/*
#if defined(__MSVC32__)
#  include <crtdbg.h>
#endif
*/
#if defined(USE_PTHREADS)
pthread_key_t admb_pthread_key=0;
#endif //USE_PTHREADS
void show_version(void);
int show_version2(void);
void print_version(ofstream& ofs);

void usage(const char* progName);
#if !defined(__MINGW__) && !defined(_WIN32)
#  include <sys/resource.h>
#endif
int tag_thread_tester(int nslaves, dvar_len_fish_stock_history& fsh);

dvar_vector choleski_solve(_CONST dvar_matrix& MM,const dvar_vector& vv,
  const prevariable& det,const int& sgn);

dmatrix stupid_test(void)
{
  int n=2;
  int m=2;
  dvector a(1,n*m);
  a.fill_seqadd(1,1);
  dmatrix A(1,n);
  int offset=0;
  for (int i=1;i<=n;i++)
  {
    A(i)=a(1+offset,m+offset).shift(1);
    offset+=m;
  }
  return A;
}
void bs(dmatrix M)
{
  cout << M << endl;
}
 void print_not_done(const char * s)
 {
  cerr << s << endl;
  cerr << " is not finished -- You ned to deal with this" << endl;
  ad_exit(1);
}
/*
 void calculate_variance_by_delta_method_noeff(const char * _s)
 {
   print_not_done("calculate_variance_by_delta_method_noeff");
 }
 void calculate_variance_by_delta_method(const char * _s)
 {
   print_not_done("calculate_variance_by_delta_method");
 }
 void make_covariance_report_for_projections(dvar_len_fish_stock_history& f)
 {
   print_not_done("make_covariance_report_for_projections");
 }
void make_correlation_report(dvar_len_fish_stock_history& fsh)
{
  print_not_done("make_correlation_report");
}

void cblas_matrix_prod(int n, dvector& x, dvector& y, dvector& z)
{
  print_not_done("cblas_matrix_prod");
}
void lapack_symmetric_eigen(int n, dvector& x, dvector& y, dvector& z, int& u)
{
  print_not_done("lapack_symmetric_eigen");
}
void check_hessian_pd(char const*)
{
  print_not_done("check_hessian_pd");
}
*/
void orth_test(int n);
  

dvariable beta(const  prevariable& a,const prevariable& b )
{
  return exp(gammln(a)+gammln(b)-gammln(a+b));
}
dvariable betaln(double a,const prevariable& b )
{
  return exp(gammln(a)+gammln(b)-gammln(a+b));
}

dvariable beta(const double a,const prevariable& b )
{
  return exp(gammln(a)+gammln(b)-gammln(a+b));
}

dvar_vector betaln(const double a,const dvar_vector& b )
{
  int mmin=b.indexmin();
  int mmax=b.indexmax();
  dvar_vector tmp(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    tmp(i)=gammln(a)+gammln(b(i))-gammln(a+b(i));
  }
  return tmp;
}

  
extern int dertest_flag;
extern "C"  {
  void ad_boundf(int i)
  {
    // so we can stop here for breakpoint
    exit(i);
  }
}

int Lapack_Choleski_Inverse(int n,int lda,MY_DOUBLE_TYPE * v);

void lapack_symmetric_eigen(dmatrix & _MM,dmatrix& Eigenvectors,
  dvector& eigenvalues,int & ierr);

dmatrix make_dmatrix(int cmin,int cmax,int rmin,int rmax,dvector & v);

int main(int argc,char * argv[])
{

  ad_exit=&ad_boundf;
  
  
  ofstream ofs1("argv1");
  ofs1 << "Arg 0" << endl;
  //exit(1);
  ofs1 << argv[0] << endl;
  cout << "Starting 1" << endl;
  cout << argc << endl;
  cout << argv[0] << endl;
  {
    ofstream ofs("argv");
    if (!ofs)
    {
      cout << "error opening file argv" << endl;
    }
    else
    {
      cout << "opened file argv" << endl;
    }
    ofs << "argc = " << argc << endl;
    ofs << "argv = " << endl;
    for (int i=0;i<argc;i++)
    {
      ofs << argv[i] << endl;
    }
  } 
  mytestcin();
  cout << "Starting 2" << endl;
  fenv_t newval;
  fegetenv(&newval);
#if !defined(__MINGW__) && ! defined(_WIN32) && !defined(_MAC)
  newval.__control_word &= ~(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif
  //newval.__control_word &= ~(FE_ALL_EXCEPT);
#if !defined(__MINGW__) && ! defined(_WIN32) && !defined(_MAC)
  //fesetenv(&newval);
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );
#elif defined(_MAC)
  feraiseexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);

#endif

  mytestcin();
  if (option_match(argc, argv, "--version") > 0) 
  {
    show_version();
    exit(EXIT_SUCCESS);
  }
  else if (option_match(argc, argv, "--help") > 0) 
  {
    usage(argv[0]);
    exit(EXIT_SUCCESS);
  }
#if !defined(__MINGW__) && !defined(_WIN32)
  signal(SIGFPE, handle_fpe_signal); 
#endif
#if defined(USE_PTHREADS)
  //int rc=pthread_key_create(&admb_pthread_key, my_do_nothing);
  int * pidentifier = new int;
  *pidentifier = 0;
  //rc = pthread_setspecific(admb_pthread_key, (void*)(pidentifier));
#endif //USE_PTHREADS
#if defined(USE_DD)
  unsigned int old_cw;
  fpu_fix_start(&old_cw);
#endif
  test_fpu();
  //cout <<"nnewlan.cpp " << pow(10.1,25.2) << endl;
  test_fpu();
/*
//#  if defined(__MSVC32__)
//   int tmpflag=_CrtSetDbgFlag(_CRTDBG_REPORT_FLAG);
//
//   tmpflag |= _CRTDBG_CHECK_ALWAYS_DF;
//   _CrtSetDbgFlag(tmpflag);
//#  endif
*/
  clogf = new ofstream("error.log");

 system("pwd");
# if (defined(_PRINT_WD_))
#   if (defined(linux))
    system("pwd > /tmp/the_wd");
#   else
    system("pwd > c:\\the_wd");
#    endif
# endif
    ad_exit=&ad_boundf;
    ad_argc=argc;
    ad_argv=argv;
#  if !defined(linux)
     unsigned cw= _control87(0, 0);
#  endif
  check_working_directory_option(argc,argv);

  mytestcin();
  deal_with_command_line(argc,argv,full_datafile_path,
    full_input_parfile_path,full_output_parfile_path);


  mytestcin();
  int ret_val=sub_main(directory_path,full_input_parfile_path,
    full_datafile_path,full_output_parfile_path);

#if defined(USE_DD)
  fpu_fix_end(&old_cw);
#  endif
  exit(ret_val);
}

newstyle_frq_data newfrqdata;

#define HERE cout << "reached line " << __LINE__ << " in " << __FILE__ << endl;

/*
// struct pvmhostinfo *ad_hostp=NULL;
// int ad_nhost=0;
// ivector * ad_stid=NULL;
// int pvm_switch=0;
// int pvm_save_switch=0;
*/

// Print out program's version number.
void show_version() {
  cout <<"Version number: " << VERSION << endl;
}
int show_version2() {
//  cout <<"Version number: " << VERSION << endl;
//  cout <<"Major Version number: " << MAJOR_VERSION << endl;
//  cout <<"Minor Version number: " << MINOR_VERSION << endl;
//  cout <<"Revision Version number: " << REVISION_VERSION << endl;
//  cout <<"Build Version number: " << BUILD_VERSION << endl;
  adstring vsnstrng;
  vsnstrng=VERSION;
//  cout <<"Version number: " << vsnstrng << endl;
//  cout <<"HERE " << endl;
//  vsnstrng=MAJOR_VERSION;
  int itmp;
  itmp=MAJOR_VERSION;  
//  cout <<"Integer major version number: " << itmp << endl;
  adstring tmpstrng2=itoa(itmp,10);
//  cout <<"adstring major version number: " << tmpstrng2 << endl;

  itmp=MINOR_VERSION;  
//  cout <<"Integer minor version number: " << itmp << endl;
  adstring tmpstrng3=itoa(itmp,10);
//  cout <<"adstring minor version number: " << tmpstrng3 << endl;

  itmp=REVISION_VERSION;  
//  cout <<"Integer revision version number: " << itmp << endl;
  adstring tmpstrng4=itoa(itmp,10);
//  cout <<"adstring revision version number: " << tmpstrng4 << endl;

  itmp=BUILD_VERSION;  
//  cout <<"Integer build version number: " << itmp << endl;
  adstring tmpstrng5=itoa(itmp,10);
//  cout <<"adstring build version number: " << tmpstrng5 << endl;

  adstring vsnagg;
  vsnagg=tmpstrng2+tmpstrng3+tmpstrng4+tmpstrng5;
//  cout <<"Concatenated adstring version number: " << vsnagg << endl;
  
  itmp=atoi(vsnagg);
  return itmp;
}
void print_version(ofstream& ofs) 
{
  ofs <<"# MULTIFAN-CL version number: " << VERSION << endl;
}


// Print out program's usage.
void usage(const char* progName) {
  cout << "MFCL " << progName << " " << VERSION << " Usage:" << endl;
  //* @todo Put usage here. *
}

mf_pvm_manager * mf_pvm = 0;

#include "svd.h";
svd svd_decomp1(dmatrix& M);

int sub_main(adstring& directory_path,adstring& full_parfile_path,
    adstring& full_datafile_path,adstring& full_output_parfile_path)
{

  /*
  ofstream of1("randoms.txt");
  random_number_generator rng(901);
  MY_DOUBLE_TYPE eps;
  MY_DOUBLE_TYPE sd=0.3;
  MY_DOUBLE_TYPE bias=0.5*sd*sd;
  for (int i=1;i<=500;i++)
  {
    eps=sd*randn(rng)-bias;
    of1 << exp(eps) << endl;
  }
  cout << "done the random numbers " << endl;
  ad_exit(1);
  */

  /*
  int ny=20;
  int a=5;
  int b=1;
  double tmp;
  tmp=ny-a+b;
  cout << tmp << endl;
  ad_exit(1);
  
  dvector aaa(1,5);
  dvector bbb(1,5);
  for (int i=1;i<=5;i++)
  {
    aaa(i)=2+i;
    bbb(i)=5+i;
  }
  cout << aaa << endl;
  cout << bbb << endl;
  for (int i=1;i<=5;i++)
  {
    bbb(i)-=aaa(i);
  }
  cout << bbb << endl;
  ad_exit(1);

  {
    int n=1000;
    int m=120;
    int k=700;
    random_number_generator rng(901);
    dmatrix M(1,n,1,k);
    dmatrix B(1,k,1,m);
    //dmatrix B(1,k,1,m);
    M.fill_randn(rng);
    B.fill_randn(rng);
    // cout << M << endl << endl;
    dmatrix C=lapack_matrix_multiplcation(M,B);
    cout << C << endl << endl;
    cout << M*B << endl << endl;
    cout << norm2(C-M*B) << endl;
    ad_exit(1);
 }
 { 
    cout << "starting luinv " << endl;
    dmatrix Minv=lapack_luinv(M);
    cout << "ending luinv " << endl;
    //cout << Minv << endl << endl;
    //cout << M*Minv << endl << endl;
    //cout << M*trans(Minv) << endl << endl;
    ad_exit(1);
  }
  {
    int nvar=n;
    random_number_generator rng(901);
    dmatrix S(1,n,1,n);
    dmatrix U(1,n,1,n);
    dmatrix VT(1,n,1,n);
    U.fill_randn(rng);
    S.initialize();
    S(1,1)=-1.0;
    for (int i=2;i<=n;i++)
    {
      S(i,i)=i/10.0;
    }
    for (int i=1;i<=n;i++)
    {
      U(i)/=norm(U(i));
      for (int j=i+1;j<=n;j++)
      {
        U(j)-=(U(j)*U(i))*U(i);
      }
    }
    VT=trans(U);
    dvector veigvecs;
    dvector vH;
    dmatrix H=make_dmatrix(1,nvar,1,nvar,vH);
    dmatrix M(1,nvar,1,nvar);
    H=U*S*VT;
    M=H;
    svd SVD=svd_decomp1(M);
    dmatrix eigvecs=make_dmatrix(1,nvar,1,nvar,veigvecs);
    dvector eigvals;
    int ierr=0;
    lapack_symmetric_eigen(nvar,vH,veigvecs,eigvals,ierr);
    cout << "S"  << endl;
    cout << S << endl << endl;
    cout << "U"  << endl;
    cout << U << endl << endl;
    cout << eigvals << endl << endl;
    cout << eigvecs << endl << endl;
  }
  */
  //ad_exit(1);
  profile_timer_1=new adtimer;
  profile_timer_2=new adtimer;
  mytestcin();
  char ** nptr=0;
  char * tmpp=new char[1000];
  strcpy(tmpp,"1");
  MY_DOUBLE_TYPE ff=strtod(tmpp,0);
  delete [] tmpp;
  tmpp=0;
  cout << "input (.frq) file: "<<full_datafile_path << endl;
  mf_pvm = new mf_pvm_manager;

  if (option_match(ad_argc,ad_argv,"-master")>0)
  {
#   if !defined(USE_ADPVM)
      cerr << "can not use parallel processing unless compiled with"
        " USE_ADPVM option" << endl;
      *clogf << "can not use parallel processing unless compiled with"
        " USE_ADPVM option" << endl;
      ad_exit(1);
#   endif
  if (mf_pvm)
    mf_pvm->pvm_switch=1;
  }
  if (option_match(ad_argc,ad_argv,"-slave")>0)
  {
#if !defined(USE_ADPVM)
    *clogf << "can not use parallel processing unless compiled with"
      " USE_ADPVM option" << endl;
    cerr << "can not use parallel processing unless compiled with"
      " USE_ADPVM option" << endl;
    ad_exit(1);
#endif
    if (mf_pvm)
      mf_pvm->pvm_switch=2;
    cout << "slave mode:  options are"<< endl;
    for (int i=0;i<ad_argc;i++)
      cout <<"nnewlan, ad_argv: " << ad_argv << endl;
  }
#if defined(USE_ADPVM)
  if (mf_pvm->pvm_switch)
  {
    //if (mf_pvm->pvm_switch==1)
    {
      FILE * flog = 0;
#    if (defined(linux))
        flog = fopen("/tmp/adpvmlog","w");
#    else
        flog = fopen("c:\\temp\\adpvmlog","w");
#    endif
      if (flog==0) 
      {
        cerr << "Error trying to open pvm log file" << endl;
      }
      pvm_catchout(flog);
    }
  } 
#endif //#if defined(USE_ADPVM)
  int result, check, i, narch;

// ********************************************************
// ********************************************************
// set up pvm stuff

#if defined(USE_ADPVM)
  if (mf_pvm->pvm_switch)
  {
    mf_pvm->setup_pvm();
  }
  if (mf_pvm->pvm_switch==1) // master
  {
    if (allocated(cl_switches))
    {
      mfsend_int_to_slaves(cl_switches.indexmax());
      mfsend_x_to_slaves(cl_switches);
    }
    else
    {
      mfsend_int_to_slaves(0);
    }
  }
  if (mf_pvm->pvm_switch==2) // slave
  {
    int nopt=mfget_int_from_master();
    if (nopt)
    {
      imatrix m=mfget_imatrix_from_master();
      if (allocated(cl_switches)) 
        cl_switches.deallocate();
      cl_switches.allocate(1,nopt,1,3);
      cl_switches=m;
    }
  }
#endif //#if defined(USE_ADPVM)

// ********************************************************
// ********************************************************
  result=0;
#if defined(USE_ADPVM)
    if (mf_pvm->pvm_switch)
    {
      if (mf_pvm->pvm_switch==1) // master
        mfsend_file_to_slaves(full_input_parfile_path);
      else if (mf_pvm->pvm_switch==2) // slave
        mfget_file_from_master(full_input_parfile_path);
    }
#endif //#if defined(USE_ADPVM)

  int tot_fi,nlint,nage; 
  int nwint=0;
  MY_DOUBLE_TYPE wshlen,wfilen;
  MY_DOUBLE_TYPE shlen,filen;
  MY_DOUBLE_TYPE tmp_len_wt_coff=0.0;

  // here is where you set up a lot of array sizes etc.
  set_gtradient_structure_globals();
  long int itmp;
  if (!ad_array_mbl_size)
    itmp=25000000L;
  else
    itmp=ad_array_mbl_size;

  pgs = new gradient_structure(itmp);
  {
    dvar_vector x(1,10);
  }
  gss=itmp;


  /*
  independent_variables XX(1,3);
  XX(1)=0.9;
  XX(2)=1.1
  XX(2)=2.1
  dvariable xx=XX(1);
  dvariable mmu=XX(2);
  dvariable ttau=XX(3);

  dvariable utmp=my_gamma_density(xx,mmu,ttau);
  */


  getroot(ad_root,full_datafile_path);
  cifstream infile(full_datafile_path);
  if (!infile)
  {
    cerr << "Error opening file " << full_datafile_path << endl;
    exit(1);
  }
  imatrix Dflags;
  int nswitch=10;
  //infile >> nswitch;
  ivector dataswitch;
  ivector _fishery_regions;
  int num_regions=1;
  int num_fisheries=1;
  int tag_group_flag=0;
  dvector _region_area;

  mytestcin();
  MY_DOUBLE_TYPE y;

  //  new c++ smart pointer
  std::shared_ptr<group_manager_1> pgroup_manager_1(nullptr);
  //group_manager_1 * pgroup_manager_1=0;

  multi_species_data * pmsd=0;
  // read in info about region relationships if any
  set_region_incidence_matrix(nswitch,infile,dataswitch,Dflags,
    num_fisheries,num_regions,_fishery_regions,tag_group_flag,_region_area,
    _data_file_version_no,pmsd);
  
  int _mfactor=dataswitch(8);
  imatrix _season_region_flags(1,_mfactor,1,num_regions);
  if (pmsd && pmsd->num_species>1)   // using region kludge
  {
    //pmsd->season_region_flags.deallocate();
    ivector nr=pmsd->num_region_by_species;
    int ns=pmsd->num_species;
    pmsd->ses_reg_recr_flags.allocate(2,ns);
    for (int is=2;is<=ns;is++)
    {
      pmsd->ses_reg_recr_flags(is).allocate(1,_mfactor,1,nr(is));
    }
  }

  if (_data_version_no< _data_file_version_no)
  {
    cerr << "Error -- Data file version is newer than executable version"
         << ".frq file version is " << _data_file_version_no << endl;
    cerr << " supported executable version is " << _data_version_no << endl;
    exit(1);
  }
  ivector move_weeks;

  imatrix dff=read_datafile_flaginfo(infile,dataswitch(2),
    _data_file_version_no,_season_region_flags,pmsd);
  if (pmsd)   // using region kludge
  {
    cout << "Need to deal with season_region_flags for"
      " multiple species" << endl;
  }

  int month_1=set_month1_flag(dataswitch);

  read_move_weeks(move_weeks,dataswitch,infile,month_1);
  

  infile >> tot_fi >> fishery_freq_record::real_nlint >> shlen >> filen
    >>  fishery_freq_record::group_factor;

  mytestcin();
  if (_data_file_version_no >=3) 
  {
    infile >> fishery_freq_record::real_nwint >> wshlen >> wfilen
      >>  fishery_freq_record::wgroup_factor;
  }
  else
  {
    fishery_freq_record::real_nwint=0.0;
    nwint=0;
    wshlen=0.0;
    wfilen=0.0;
    fishery_freq_record::wgroup_factor=0.0;
  }
  if (_data_file_version_no >=5) 
  {
    infile >> fishery_freq_record::age_nage;
    infile >> fishery_freq_record::age_age1;
  }
  if (_data_file_version_no <3) 
  {
    infile >> tmp_len_wt_coff;
  }

  filen*=fishery_freq_record::group_factor;

  if (!age_data_flag)
  {
    nlint=fishery_freq_record::real_nlint/fishery_freq_record::group_factor;
    if(fishery_freq_record::real_nlint%fishery_freq_record::group_factor)
	  nlint++;
    if (fishery_freq_record::real_nwint)
    {
      nwint=fishery_freq_record::real_nwint/fishery_freq_record::wgroup_factor;
      if(fishery_freq_record::real_nwint%fishery_freq_record::wgroup_factor)
	  nwint++;
    }
  }
  else
  {
    nlint=fishery_freq_record::real_nlint;
  }

  fishery_freq_record_array  flfra(1,tot_fi,nlint,shlen,filen,
      nwint,wshlen,wfilen);

 
 
  if (pmsd)   // using region kludge
  {
    allocate_species(flfra,pmsd->num_species);
    pmsd->num_real_fisheries=num_fisheries;
  }
  {
   
   /*
    for (int i=1;i<=tot_fi;i++)
    {
      flfra(i).freq.deallocate();
    }
   
    MY_DOUBLE_TYPE * pd;
    //pd = (MY_DOUBLE_TYPE *)  0x14920000;
    pd = (MY_DOUBLE_TYPE *)  0x14a60000;
    // *pd = 100.0;  
    */
   
  }

 
  int num_real_regions=num_regions;
  if (pmsd)
  {
    num_real_regions=pmsd->num_real_regions;
  }
  
  ivector regmin(1,num_real_regions);
  ivector regmax(1,num_real_regions);
  fishery_freq_record_array * pkludged_flfra=0;
  {
    infile >> flfra; // read in the fishery_header_records;
  
    if (!pmsd)
    {
      fishery_freq_record_array  flfra1(1,tot_fi,nlint,shlen,filen,
         nwint,wshlen,wfilen);

      sort(flfra,_fishery_regions,flfra1,regmin,regmax);
    //*************************************************
      // For DEBUG purposes
//      print_flfra_for_debug(flfra,1);
    //*************************************************

    }
    else
    {
      regmin.deallocate();
      regmax.deallocate();
      regmin.allocate(1,pmsd->num_kludged_regions);
      regmax.allocate(1,pmsd->num_kludged_regions);
      for (i=flfra.indexmin();i<=flfra.indexmax();i++)
      {
        int sp=flfra(i).get_species();
        flfra(i).fishery=(sp-1)*num_fisheries+flfra(i).fishery;
      }

      fishery_freq_record_array  flfra1(1,tot_fi,nlint,shlen,filen,
         nwint,wshlen,wfilen);
  
      ivector kludged_fishery_regions(1,pmsd->num_species*num_fisheries);
      
      kludged_fishery_regions(1,num_fisheries)=_fishery_regions;
      for (int i=2;i<=pmsd->num_species;i++) 
      {
        int offset=(i-1)*num_fisheries;
        kludged_fishery_regions(1+offset,num_fisheries+offset).shift(1)=
          _fishery_regions+(i-1)*num_real_regions;
      }
      _fishery_regions.deallocate();
      _fishery_regions=kludged_fishery_regions;

      sort(flfra,kludged_fishery_regions,flfra1,regmin,regmax);
      num_regions=pmsd->num_kludged_regions;
      dvector dtmp(1,pmsd->num_real_regions);
      dtmp=_region_area;
      _region_area.deallocate();
      _region_area.allocate(1,pmsd->num_kludged_regions);
      int offset=0;
      for (int i=1;i<=pmsd->num_species;i++)
      {
        _region_area(1+offset,pmsd->num_real_regions+offset).shift(1)
          =dtmp;
        offset+=pmsd->num_real_regions;
      }
    }
  }

  // make sure the first year is 1
  int direction_flag=0;
  int min_year=0;
  int first_data_month=0;
  //if (pmsd==0)
  imatrix month_detector=
      normalize_year(flfra,month_1,direction_flag,min_year,first_data_month);

    imatrix months_used=months_analyzer(month_detector);
  //else
  //  normalize_year(*pkludged_flfra,month_1,direction_flag,min_year,first_data_month);
  int _first_time=0;

  //if (_mfactor!=0 && _mfactor!=1)
  //{
  //  month_doubling_kludge(flfra,_mfactor,regmin,_first_time);
  //}

  if (!infile)
  {
    cerr << "Error reading headers from file "
      << (char*)full_datafile_path <<endl;
    return 1;
  }
#if defined(close)
#  undef close
  infile.close();
#  define close _close
#else
  infile.close();
#endif
      test_fpu();
  *clogf << " readfrq file" << endl;
  par_cifstream * pinfile=NULL;
  ivector parest_flags(1,400);
  parest_flags.initialize();
  pinfile= new par_cifstream ((char*)full_input_parfile_path);
  if (!pinfile || !(*pinfile))
  {
    cerr << " Error trying to open par file " << full_input_parfile_path
         << endl;
    exit(1);
  }
  ivector historical_age_flags(1,200);
  imatrix historical_multspp_age_flags; //NMD_23Apr2015tmp;
  imatrix historical_multspp_parest_flags; 
  if(pmsd)
  {
    historical_multspp_age_flags.allocate(2,pmsd->num_species,1,200);
    historical_multspp_parest_flags.allocate(2,pmsd->num_species,1,400);
  }
  ivector historical_parest_flags(1,400);

  int isp=1; //NMD_23Apr2015
  if(pmsd)
  {
    isp*=pmsd->num_species;
  }
  int nfsh=isp*num_fisheries;
  imatrix historical_fish_flags(1,nfsh,1,100);   //NMD_23Apr2015
  historical_age_flags.initialize();
  historical_multspp_age_flags.initialize(); //NMD_23Apr2015
  historical_parest_flags.initialize();
  historical_fish_flags.initialize();
  // PAR file.
  if(!makeparflag)
  { 
    ivector itmp(1,200);
    // Input historical flags
    cifstream cif((char*)full_input_parfile_path);
    if(!cif) cout << "ERROR !cif" << endl;
    cif >> itmp;
  mytestcin();
    if(itmp(200) > 1047)
    {
      // search for historical flags NMD_24Mar2015
      if (search_for_string(cif, "Historical_flags","End_historical_flags"))
      {
        cif >> historical_parest_flags(1,200);
        int hist_file_version_no = historical_parest_flags(200);
        if (hist_file_version_no > 1047)
        {
          cif >> historical_parest_flags(201,400);
        }
        if(pmsd && itmp(200) > 1054)
        {
          if (pmsd->num_species>1)
            cif >>  historical_multspp_parest_flags;
        }
        cif >> historical_age_flags;
        if (pmsd)    //NMD_23Apr2015
        {
          if (pmsd->num_species>1)
            cif >>  historical_multspp_age_flags;
        }
        cif >> historical_fish_flags;
      }      
//      for (int i=1;i<=num_fisheries;i++)
      for (int i=1;i<=nfsh;i++)   //NMD_23Apr2015
      {
        if (historical_fish_flags(i,74)==0)
          historical_fish_flags(i,74)=1;
      }
      if (search_for_string(cif, "Historical_zero_effdev_flag",
        "End_Historical_zero_effdev_flag"))
      {
        if (allocated(pgroup_manager_1->old_zero_effdev_flags))
        {
          pgroup_manager_1->old_zero_effdev_flags.deallocate();
        }
        //pgroup_manager_1->old_zero_effdev_flags.allocate(1,nfsh,
        //  1,num_fish_times);
        //cif >> pgroup_manager_1->old_zero_effdev_flags;
      }
      // End of historical flags input
    } //NMD_24Mar2015

    *pinfile >> parest_flags(1,200);

    _file_version_no = parest_flags(200);
    if (_file_version_no>=1046) 
    {
      *pinfile >> parest_flags(201,400);
    }
    else if (_file_version_no>=1036) 
    {
      *pinfile >> parest_flags(201,300);
    }
    if (_version_no < _file_version_no )
    {
      cerr << "Error -- PAR file has higher version number than program"
         << endl;
      cerr << "file version is " << _file_version_no << endl
         << " program version is " << _version_no << endl;
      return 1;
    }
    if (_file_version_no>=1055) 
    {
      if (pmsd && pmsd->num_species>1)
      {
        *pinfile >> pmsd->parest_flags;
        pmsd->historical_age_flags=historical_multspp_age_flags;
        pmsd->historical_parest_flags=historical_multspp_parest_flags;
      }
    }
    else
    {
      if (pmsd && pmsd->num_species>1)
      {
        pmsd->parest_flags.initialize();
      }
    }
    *pinfile >> nage;
    if (pmsd)
    {
      *pinfile >> pmsd->nage;
    }
  }
  // INI file.
  else
  {
    _file_version_no=_version_no;  
	  // read in the ini version number
    *pinfile >> _ini_version_no;
    parest_flags(197)=1;
      // Old format.
    if (_ini_version_no < 1000) 
    {
      nage = _ini_version_no;
      _ini_version_no = 0;
    }
    else // New format.
    {
		  _ini_version_no = _ini_version_no - 1000;
      /*
      if (_ini_version_no < 0 || _ini_version_no > 5)
      {
        cerr << "Illegal ini version number -- assuming old format"
          << endl;
        cout << "Illegal ini version number -- assuming old format"  //NMD
          << endl;
        ad_exit(1);
      }	
      */
      *pinfile >> nage;
      if (pmsd)
      {
        *pinfile >> pmsd->nage;
      }
    }
  }
  if (pmsd)
  {
    if (!(*pinfile))
    {
      cerr << "Error reading parameters from cifstream file " 
        << (char*)full_input_parfile_path <<endl;
      ad_exit(1);
    }
//NMD 4Nov2011
    if (allocated(pmsd->rec_init_diff)) pmsd->rec_init_diff.deallocate();
    pmsd->rec_init_diff.allocate(2,pmsd->num_species);
    if (allocated(pmsd->totpop)) pmsd->totpop.deallocate();
    pmsd->totpop.allocate(2,pmsd->num_species);
    int tmp = max(pmsd->nage);
    if (allocated(pmsd->pmature)) pmsd->pmature.deallocate();
    pmsd->pmature.allocate(2,pmsd->num_species,1,tmp);
    if (allocated(pmsd->cpmature_at_length))    //NMD_Aug28_2018
      pmsd->cpmature_at_length.deallocate();
    pmsd->cpmature_at_length.allocate(2,pmsd->num_species,1,nlint);
    //NMD_Aug28_2018
    if (allocated(pmsd->age_pars)) pmsd->age_pars.deallocate();
    for (int is=2;is<=pmsd->num_species;is++)
    {
      int tmp = pmsd->nage(is);
      pmsd->age_pars.allocate(2,pmsd->num_species,1,10,1,tmp);
    }
//NMD 4Nov2011
  }
  

    //if (_sd_flag)
    //{
    //  cerr << " setting parest_flags(30)=0 " << endl;
    //  parest_flags(30)=0;
    //}
    
    if (parest_flags(41)>0)
    {
      mean_weights_kludge.allocate(1,nage);
      cifstream cif("meanwght.dat");
      cif >> mean_weights_kludge;
      if (!cif)
      {
        cerr << "Error reading mean weight data from file meanwght.dat"
             << endl;
        return 1;
      }
    }

    if (!(*pinfile))
    {
      cerr << "Error reading parameters from cifstream file " 
        << (char*)full_input_parfile_path <<endl;
      return 1;
    }


    
      const dvar_len_fish_stock_history & rfsh= 
      flfra.get_history_data(tag_group_flag,num_regions,nage,
            parest_flags,regmin,regmax,dataswitch,Dflags,_region_area,
            month_1,_mfactor,_first_time,move_weeks,direction_flag,
            _season_region_flags,pmsd);
      //else
      //  const dvar_len_fish_stock_history * rfsh= 
      //    &(pkludged_flfra->get_history_data(tag_group_flag,num_regions,nage,
      //    parest_flags,regmin,regmax,dataswitch,Dflags,_region_area,
      //    month_1,_mfactor,_first_time,move_weeks,direction_flag,
      //    _season_region_flags,pmsd));

    dvar_len_fish_stock_history fsh= (dvar_len_fish_stock_history&) rfsh;

    fsh.setup_some_stuff(tmp_len_wt_coff,wshlen,wfilen,nwint,num_fisheries,
      dff,pinfile,month_1,_first_time, _mfactor,_fishery_regions,
      first_data_month,pmsd);

   // get current compilation version number
    fsh.current_version=show_version2();
  
   //  *************************************************************
   //    read in fish_flags before doing switch changes
   //  *************************************************************

    fsh.read_fish_flags(pinfile);

    fsh.months_used=months_used;
   //  *************************************************************
   //  *************************************************************

    fsh.age_flags_sanity_check();

    // set some flags
    fsh.set_some_flags(dataswitch,min_year);

    if(!(!cl_switches))
    {
      do_switch_changes(cl_switches,fsh);
    }
    if (tag_group_flag>0)
    {
      fsh.read_tagging_data(ad_root,pmsd);
    }
    else    //NMD_jan17-19
    {
      if (pmsd) pmsd->tag_index=0;
    }    //NMD_jan17-19

    fsh.age_flags_sanity_check();

#if defined(USE_ADPVM)
    if (mf_pvm->pvm_switch)
    {
      *clogf << " begin tag assignments" << endl;
      set_tag_group_assignments(fsh);
      *clogf << " finsihed tag assignments" << endl;
    }
#endif //#if defined(USE_ADPVM)

    pcfsh=&fsh;

      test_fpu();
  //cout << "calling heapcheck" << endl;
  //cout << "heapcheck = " << heapcheck() << endl;
  *clogf << " begin readpar file" << endl;
    test_the_pointer();

    fsh.set_missing_totcatch_flags();
    fsh.historical_parest_flags=historical_parest_flags;
    fsh.historical_age_flags=historical_age_flags;
    fsh.historical_fish_flags=historical_fish_flags;
    if (pmsd)    //NMD_23Apr2015
    {
      int nry=(fsh.nyears-1)/fsh.num_seasons+1;
      ivector nr=fsh.pmsd->num_region_by_species;
      for(int is=2;is<=fsh.pmsd->num_species;is++)
      { 
        fsh.pmsd->historical_age_flags(is)=historical_multspp_age_flags(is);
        fsh.pmsd->orth_recr_all(is).allocate(1,fsh.num_seasons*nr(is),0,nry-1);
        fsh.pmsd->yearly_recr_all(is).allocate
          (1,fsh.num_seasons*nr(is),0,nry);
        fsh.pmsd->orth_recr_all(is).initialize();
        fsh.pmsd->yearly_recr_all(is).initialize();
      }
    }
  
  mytestcin();
    *pinfile >> fsh;
  mytestcin();
    test_the_pointer();
  *clogf << " end readpar file" << endl;
    fsh.last_real_year=fsh.nyears;
    if (sum(fsh.data_fish_flags(2)))
    {
      normalize_projection_year(month_1,direction_flag,min_year,
        fsh.data_fish_flags,fsh);
    }
    else
    {
      fsh.projection_year=fsh.nyears;
      fsh.projection_month=13;
    }

    if (fsh.age_flags(101)>0)
    {
      fsh.read_recruitment_env_data(ad_root);
    }
    int nopt,nnopt;
    if ((nopt=option_match(ad_argc,ad_argv,"-move",nnopt))>0)
    {
      movement_coffs_file=ad_argv[nopt+1];
      fsh.age_flags(114)=1;
      fsh.read_movement_coffs(movement_coffs_file);
    }
    test_the_pointer();

    if (option_match(ad_argc,ad_argv,"-dbgio",nnopt)>0)
    {
      par_ofstream ofs(adstring("tmpfile"));
      ofs << fsh;
      exit(0);
    }
    // *******************************************************
    // *******************************************************
    // *******************************************************
    // *******************************************************
    if (fsh.parest_flags(363))
    {
        cout << "converting bet_128.hes " << endl;
        cout << " sizeof(double)" << sizeof(double) << endl;
        cout << " sizeof(MY_REAL_DOUBLE)" << sizeof(MY_REAL_DOUBLE) << endl;
        adstring fname="bet_128.hes";
        uistream ifs(fname);
        if (!ifs)
        {
          cerr << "Error trying to open file " << fname << endl;
        }
        int nvar;
        ifs >> nvar;
        if (!ifs)
        {
          cerr << "Error reading nvar from " << fname << endl;
        }
        dvector vH;
        //nvar=1000;
        dmatrix Hess=make_dmatrix(1,nvar,1,nvar,vH);
        dvector vhessinv;
        dmatrix hessinv=make_dmatrix(1,nvar,1,nvar,vhessinv);
        dvector row(1,nvar);

        ifs >> Hess;
        if (!ifs)
        {
          cerr << "Error reading hessian from " << fname << endl;
        }
        {
          ofstream ofs("bet.hes_converted_to_64_bit",std::ios::binary);
          ofs.write((const char *)(&nvar),sizeof(int));
          for (int i=1;i<=nvar*nvar;i++)
          {
            MY_REAL_DOUBLE tmp=(MY_REAL_DOUBLE)(vH(i));
            char ctmp[sizeof(MY_REAL_DOUBLE)];
            if (isnan(tmp))
            {
              cout << "got nan at i = " << endl;
            }
            ofs.write(ctmp,sizeof(MY_REAL_DOUBLE));
          }
          ofs.flush();
        }
        {
          ifstream ifs("bet.hes_converted_to_64_bit",std::ios::binary);
          int nvar=0;
          char ctmp[sizeof(int)];
          ifs.read(ctmp,sizeof(int));
          for (int i=1;i<=nvar*nvar;i++)
          {
            MY_REAL_DOUBLE tmp=0.0;
            char ctmp1[sizeof(MY_REAL_DOUBLE)];
            ifs.read(ctmp1,sizeof(MY_REAL_DOUBLE));
            if (isnan(tmp))
            {
              cout << "got nan at i = " << endl;
            }
          }
        }
        ad_exit(1);
      }
    // *******************************************************
    // *******************************************************
    // *******************************************************
    // *******************************************************
    if(!(!cl_switches))
    {
      do_switch_changes(cl_switches,fsh);
    }

    fsh.set_year_flags_for_recruitment();
    fsh.allocate_fishery_projection_flags();
    fsh.allocate_grouped_fishery_projection_flags();

    fsh.resize_if_necessary();

    if (sum(fsh.data_fish_flags(2)))
    {
      fsh.set_fishery_projection_flags();
      if (colsum(fsh.fish_flags,32))
        fsh.set_grouped_fishery_projection_flags();
    }

    // for grouped fishery catchability calculations
    if (sum(column(fsh.fish_flags,29)))
    {
      fsh.do_grouped_between_time_calculations();
      if (_file_version_no>=1022)
      {
        int ntmp;
        *pinfile >> ntmp;
        if (ntmp) 
        {
          *pinfile >> fsh.grouped_catch_dev_coffs;
          if (!(*pinfile)) cerr << "Error trying to read "
            "grouped_catch_dev_coffs" << endl;
        }
        else
          fsh.grouped_catch_dev_coffs.initialize();
      }
      else
      {
        fsh.grouped_catch_dev_coffs.initialize();
      }
    }       
    /*   //NMD_7jun-19
    if (fsh.parest_flags(311))
    {
      fsh.make_tail_compressed_samples();
    }
    if (fsh.parest_flags(301))
    {
      fsh.make_tail_compressed_weight_samples();
    }
    */   //NMD_7jun-19


    // BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
    // check historical versus current version number
    if (fsh.historical_version == 0)
    {
      cout << "Input solution .par was generated by an unknown version" << endl;
      cout << "Check the version number in the plot.rep file" << endl;
    }
    else if (fsh.historical_version != fsh.current_version)
    {
      cout << "Input solution .par was generated by an earlier version" << endl;
    }
    else if (fsh.historical_version == fsh.current_version)
    {
      cout << "Input solution .par was generated by the current version" << endl;
    }

    fsh.allocate_plotstuff();

    fsh.set_effective_length_and_weight_sizes();
    fsh.allocate_optional_stuff();
    
    if (fsh.parest_flags(311))   //NMD_7jun-19
    {
      fsh.make_tail_compressed_samples();
    }
    if (fsh.parest_flags(301))
    {
      fsh.make_tail_compressed_weight_samples();
    }                            //NMD_7jun-19

//    fsh.allocate_some_tag_stuff();
    if (tag_group_flag>0)    //NMD_Jan14-18
    {
      fsh.allocate_some_tag_stuff();
    }
    
    set_true_months(fsh);

    fsh.setup_sel_ptr();
    fsh.new_set_selectivity_breaks();
   
    if (sum(column(fsh.fish_flags,29)))
    {
      fsh.set_grouped_true_months();
    }

    fsh.get_first_unfixed_year();

    if (sum(column(fsh.fish_flags,47))) get_seasonal_catchability_pars_index(fsh);
      test_fpu();
    get_global_fish_periods(fsh);
      test_fpu();

    fsh.effort_multiplier_for_cobb_douglas_production_function();

    normalize_effort_data(fsh);
    fsh.process_effort_weight_data();
   
    // new code for effort weights   DF 22/04/08
    get_effort_data_by_fishery(fsh);

    //fsh.do_grouped_tags_between_time_calculations();
    fsh.set_sel_seasons();

    if (!(*pinfile))
    {
      cerr << "Error reading parameters from cifstream file " 
        << (char*)full_input_parfile_path <<endl;
      return 3;
    }
      test_fpu();

    delete pinfile;
    pinfile=NULL;

    if (tag_group_flag>0)
    {
      fsh.var_convert_tag_lengths_to_age();
      fsh.observed_tag_catch_calc(direction_flag);
      fsh.observed_tag_catch_by_length_calc();
      fsh.tot_tags_catch();
    }

    fsh.set_zero_effdev_flag();

    // get the environmental stuff
    read_environmental_data(full_datafile_path,fsh);

    if(!(!cl_switches))
    {
      do_switch_changes(cl_switches,fsh);
    }

    fsh.get_fishery_group_pointers();

    do_between_time_calculations(fsh);

    test_the_pointer();
      test_fpu();

    fsh.num_real_years=(fsh.last_real_year-1)/_mfactor+1;


    if (fsh.parest_flags(155)>0)
    {
      //fsh.get_orth_recr_info();
      
      if (!pmsd || pmsd->num_species==1)
      {
        // This routine returns an orthonormal basis for the vector space of num seasons x num regions matrices
        // (think of them as an  num seasons x num regions dimensional  vector space)
        // Don't know if turning off recruitment for certain regions/seasons has been implemented within this calculation
        // This is a new test routine which deals with the season-region matrices only
        //fsh.get_new orthogonal_season_region_recruitment_poly_info();
        // This is the old routine which deals appears to have a bunch of 
        // uneccessary code in it. 
        fsh.get_orthogonal_recruitment_poly_info();
        fsh.get_orthogonal_polynomial_weights();
      }
      else
      {
        for (int is=1;is<=pmsd->num_species;is++)
        {
          pmsd->current_species=is;
          fsh.get_orthogonal_recruitment_poly_info();
          fsh.get_orthogonal_polynomial_weights();
        }
      }
    }
    else if (fsh.parest_flags(155)<0)
    {
      if (fsh.historical_parest_flags(155)<0)
      {
        fsh.initial_orthp_estimate_flag=1;
        fsh.get_corresponding_orthogonal_coefficients();
        fsh.initial_orthp_estimate_flag=0;
      }
      fsh.new_get_orthogonal_recruitment_poly_info();
    }

    double chkninit=value(sum(fsh.Ninit_standard));  //NMD_24aug_2021
    if (fsh.parest_flags(155)!=0 && fsh.historical_parest_flags(155)==0
        && chkninit>0.0) 
    {
      fsh.initial_orthp_estimate_flag=1;
      fsh.get_orthogonal_recruitment_poly_info();
      fsh.get_orthogonal_polynomial_weights();
      fsh.get_corresponding_orthogonal_coefficients();
      fsh.initial_orthp_estimate_flag=0;
    }
    

    if (fsh.parest_flags(240))
    {
      fsh.age_len_sample_size.allocate(1,fsh.num_regions,
        1,fsh.num_fish_periods,1,fsh.num_fish_incidents);
      fsh.age_len_sample_size.initialize();
      adstring alname= ad_root + ".age_length";
      cifstream cif(alname);
      if (!cif)
      {
        cerr << "Error trying to open file " << alname << endl;
        ad_exit(1);
      }
      read_age_length_records(cif,fsh.nlint,fsh.nage,fsh);
      if (!cif)
      {
        cerr << "Error trying to read age length records from file  " 
             << alname << endl;
        ad_exit(1);
      }
      fsh.get_age_length_periods();

    }

    if (fsh.age_flags(92)) do_robust_mean_spread_calcs(fsh);
    int quit_flag=0;
    // estimate some parameters
    // parest_flags(30) is the switch for automatic parameter relaxation
    cout << "testing threaded tag_catch_equations" << endl;
    cout << "enter 1 for test or 0 for no test" << endl;
    int iflag=parest_flags(239);
    fsh.threaded_tag_flag=0;
    //cin >> iflag;
    if (iflag==1)
    {
      fsh.threaded_tag_flag=1;
      // ?????????????
      tag_thread_tester(1,fsh);
      cout << "Slave number fo master is "
           << ad_comm::pthread_manager->get_slave_number() << endl;
      //ad_exit(1);
    }

    fsh.check_for_no_catch_or_effort();

    if (fsh.sum_ff79_flag)
    {
      ivector ff79=column(fsh.fish_flags,79);
      ivector ff4=column(fsh.fish_flags,4);
      zero_effdev_sanity_check(fsh.effort_dev_coffs,ff4,fsh.zero_effdev_flag,
        ff79);
    }
    {
    do
    {
      //fsh.parest_flags(30)=0;
      delayed_infile="";
      delayed_outfile="";
      {

        test_fpu();
        //cout <<"nnewlan, pflg30: " << fsh.parest_flags(30) << endl;
        //cout <<"nnewlan, pflg32: " << fsh.parest_flags(32) << endl;
      
        char ch;
        if (fsh.parest_flags(30))
        {
          switch (fsh.parest_flags(32))
          {
          case 1: 
            fsh.albacore_control_switches();
           cout << "calling albacore" << endl;
            break;
          case 2:
            cout << "calling kleiber control " << endl;
            fsh.set_shark_control_switches();
            break;
	        
          case 3:
            cout << "calling grid control " << endl;
            fsh.set_gridsearch_control_switches();
            break;
	        
          case 4:
            cout << "calling Yukio control " << endl;
            fsh.set_Yukio_control_switches();
            break;
          case 5:
            cout << "calling Takeuchi_san control " << endl;
            fsh.set_Takeuchi_san_control_switches();
            break;
          case 6:
            cout << "calling fsh.set_shark_no_vonb_control_switches" << endl;
            fsh.set_shark_no_vonb_control_switches();
            break;
          case 7:
            cout << "calling fsh.set_shark_no_vonb_var_control_switches" 
             << endl;
            fsh.set_shark_no_vonb_var_control_switches();
            break;
          default:
	    fsh.set_control_switches();
          }
        }
        // add parameters for cub-spline interpolated selectivity
        // if desired
        ivector ff62=column(fsh.fish_flags,62);
        if (sum(ff62)>0)
        {
          fsh.add_cs_selectivity_coffs();
        }
        if (fsh.parest_flags(376)>1)
        {
          fsh.add_equilibrium_selectivity_coffs();
        }


        test_the_pointer();

        fsh.check_flag_stuff();

        fsh.fishing_incident_times=fsh.get_time();
        if (sum(fsh.data_fish_flags(2)))   //NMD_27jun2022
        {
          fsh.fishing_incident_real_times=fsh.get_real_time();
        }
        else
        {
          fsh.fishing_incident_real_times=fsh.fishing_incident_times;
        }
        if (fsh.parest_flags(378))
        {
          fsh.new_build_implicit_catch_effort_design_matrix();
          //fsh.implicit_fml_bounds=
          //  fsh.get_bounds_for_implicit_catch_effort_design_matrix();
        }
        if (sum(fsh.data_fish_flags(2)) && !fsh.parest_flags(378)
           && !fsh.parest_flags(377))   //NMD_7may2024
        {
          fsh.new_do_build_part_for_projections1();
          for (int fi=1;fi<=fsh.num_fisheries;fi++)
          {    
            fsh.q_level_finyr(fi).allocate(1,fsh.fshtms_finyr(fi));
          }
          fsh.q_level_finyr.initialize();
        }

        int bdflag=0;
        do
        { 
          {
            // calculate the number of active parameters
            int nvar=fsh.nvcal();
            gradient_structure::set_USE_FOR_HESSIAN(fsh.parest_flags(198));
#if defined(USE_ADPVM)
            switch (mf_pvm->pvm_switch)
            {
            case 0:
#endif //#if defined(USE_ADPVM)

  if (parest_flags(360)>0)
  {
    ivector tf2=column(fsh.tag_flags,2);
    ivector tf3=column(fsh.tag_flags,3);
    fsh.tag_mortality_group.safe_allocate(tf2,tf3);
  }

  int nopt,nnopt;
  nopt=0;
  if ((nopt=option_match(ad_argc,ad_argv,"-maxfn",nnopt))>0)
  {
    fsh.parest_flags(1)=atoi(ad_argv[nopt+1]);
  }
  // do tail compression for length data. this is not
  // always needed
  if (parest_flags(320))  //NMD_15Aug2016
  {
    len_tail_compress(fsh);
  }  //NMD_15Aug2016
   
  if (parest_flags(330))  //NMD_15Aug2016
  {
    wght_tail_compress(fsh);
  }  //NMD_15Aug2016
   
  if (fsh.parest_flags(241))
  {
    fsh.do_simulation_stuff();
  }

  fsh.fishing_incident_times=fsh.get_time();
  if (sum(fsh.data_fish_flags(2)))   //NMD_27jun2022
  {
    fsh.fishing_incident_real_times=fsh.get_real_time();
  }
  else
  {
    fsh.fishing_incident_real_times=fsh.fishing_incident_times;
  }


  fsh.ff92sum=sum(column(fsh.fish_flags,92));
  if (fsh.ff92sum)
  {
    fsh.calculate_survey_indices();
    remove("grouped_cpue_obs");
    remove("grouped_cpue_pred");
    remove("cpue_obs");
    remove("cpue_pred");
    remove("cpue_obs_mls");
    remove("cpue_pred_mlv");
  }

  mytestcin();
  fsh.do_all_sanity_checks();
  fitting_procedure(fsh,quit_flag);
#if defined(USE_ADPVM)
              break;
            case 1:
 // cout << "calling heapcheck" << endl;
  //cout << "heapcheck = " << heapcheck() << endl;
              *clogf << " starting master fit" << endl;
	      fitting_procedure(fsh,quit_flag);
	      //master_fitting_procedure(fsh,quit_flag);
              break;
            case 2:
              *clogf << " starting slave fit" << endl;
	      slave_fitting_procedure(fsh);
              break;
            default:
              cerr << "Illegal value for pvm_switch = " << mf_pvm->pvm_switch <<endl;
              ad_exit(1);
            }
#endif //#if defined(USE_ADPVM)

            // update parameter values in len_fish_stock_history 
#         if defined(USE_ADPVM)
            if (mf_pvm->pvm_switch)
            {
              if (mf_pvm->pvm_switch==1) // master
                mfsend_file_to_slaves(full_output_parfile_path);
              else if (mf_pvm->pvm_switch==2) // slave
                mfget_file_from_master(full_output_parfile_path);
            }
#         endif
            read_again(fsh,full_output_parfile_path);
            // check bounds on the mean lengths
#         if defined(USE_ADPVM)
            if (mf_pvm->pvm_switch==0)
            {
#         endif
              bdflag=check_mean_length_bounds(fsh);
#         if defined(USE_ADPVM)
            }
            if (mf_pvm->pvm_switch==1) // master
            {
              bdflag=check_mean_length_bounds(fsh);
              mfsend_int_to_slaves(bdflag);
            }
            if (mf_pvm->pvm_switch==2) // master
            {
              bdflag=mfget_int_from_master();
            }
#         endif
            char ch;
            if (bdflag)
            {
              cout << "relaxing bounds on mean length" << endl;
              fsh.parest_flags(199)=0;
            }
          }
        } while (bdflag);
#     if defined(USE_ADPVM)
        if (mf_pvm->pvm_switch==1) // master
        {
          mfsend_int_to_slaves(quit_flag);
        }
        if (mf_pvm->pvm_switch==2) //slave 
        {
          quit_flag=mfget_int_from_master();
        }
#     endif

        if (quit_flag)
        {
          if (quit_flag=='Q' && fsh.parest_flags(30) &&
            fsh.parest_flags(20)>1)
          {
            fsh.parest_flags(20)--;
          }
          cout <<"executing break" << endl;
          break;         // either quit or go on to next search
        }
      }
    }
    while (fsh.parest_flags(30));
  }
  if (quit_flag ==1)
  {
    return 3;
  }
  exit(0);
}

int check_mean_length_bounds(dvar_len_fish_stock_history& fsh)
{
  int bdflag=0;
  MY_DOUBLE_TYPE diff = fsh.fmaxl-fsh.fminl;
  if ((fsh.vb_coff(2)-fsh.fminl) < .01*diff)
  {
    bdflag=1;
    fsh.fminl-=.3*diff;
    fsh.fmaxl-=.2*diff;
  }
  if ((fsh.fmaxl-fsh.vb_coff(2)) < .01*diff)
  {
    bdflag=1;
    fsh.fminl+=.2*diff;
    fsh.fmaxl+=.3*diff;
  }
  return bdflag;
}

  void read_again(dvar_len_fish_stock_history& fsh,adstring& file_root_par)
  {
    char ch;

    par_cifstream infile1((char*)full_output_parfile_path);
    ivector tmp;
    tmp.allocate(fsh.parest_flags);

    ivector itmp(1,200);
    // Input historical flags
    cifstream cif((char*)full_output_parfile_path);
    if(!cif) cout << "ERROR !cif" << endl;
    cif >> itmp;
    if(itmp(200) > 1047){
// search for historical flags NMD_30Mar2015
      if (search_for_string(cif, "Historical_flags","End_historical_flags"))
      {
        cif >> fsh.historical_parest_flags(1,200);
        int hist_file_version_no = fsh.historical_parest_flags(200);
        if (hist_file_version_no > 1047)
        {
          cif >> fsh.historical_parest_flags(201,400);
        }
        if( hist_file_version_no > 1054)
        {
          if (fsh.pmsd && fsh.pmsd->num_species>1)
          {
            cif >>  fsh.pmsd->historical_parest_flags;
          }
        }
        cif >> fsh.historical_age_flags;
        if (fsh.pmsd)    //NMD_23Apr2015
        {
          if (fsh.pmsd->num_species>1)
            cif >>  fsh.pmsd->historical_age_flags;
        }
        cif >> fsh.historical_fish_flags;
      }      
      for (int i=1;i<=fsh.num_fisheries;i++)
      {
        if (fsh.historical_fish_flags(i,74)==0)
          fsh.historical_fish_flags(i,74)=1;
      }
      // End of historical flags input
    } //NMD_30Mar2015

    infile1 >> tmp(1,200);
    fsh.parest_flags(200) =tmp(200);
    _file_version_no = fsh.parest_flags(200);
    if (_file_version_no>=1046) 
    {
      infile1 >> tmp(201,400);
    }
    else if (_file_version_no>=1036) 
    {
      infile1 >> tmp(201,300);
    }
    if(fsh.pmsd && itmp(200) > 1054)
    {
      if (fsh.pmsd->num_species>1)
        infile1 >>   fsh.pmsd->parest_flags;
    }
    int n;
    infile1 >> n;
    if (fsh.pmsd)
    {
       infile1 >> fsh.pmsd->nage;
    }
    infile1 >> fsh.age_flags;
    if (fsh.pmsd)
    {
      if (fsh.pmsd->num_species>1)
        infile1 >>   fsh.pmsd->age_flags;
    }
    check_number(infile1,871,fsh.parest_flags);
    fsh.read_fish_flags(&infile1);
  
    infile1 >> fsh;
    // save current phase so you can put it back
    int tmp_ipar=fsh.parest_flags(20);
    fsh.parest_flags=tmp;
    fsh.parest_flags(20)=tmp_ipar;
    if (!infile1)
    {
      cerr << "Error re-reading parameters from file " 
           << (char*)full_output_parfile_path <<endl;
      exit(1);
    }
#if defined(close)
#  undef close
    infile1.close();
#  define close _close
#else
    infile1.close();
#endif
  }
  void set_command_line_switches(imatrix& cls,int argc,
    char * argv[], int opos)
  {
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
    //ad_exit(11);
    int ii=opos+2;
    if (ii> (argc-1) )
    {
      cerr << "Error not enough command line arguments"
       " for  command line switch changes" << endl;
      exit(1);
    }
    int mmin=cls.rowmin();
    int mmax=cls.rowmax();
    for (int i=mmin;i<=mmax;i++)
    {
      if (ii+2> (argc-1) )
      {
        cerr << "Error not enough command line arguments"
         " for  command line switch changes" << endl;
        exit(1);
      }
      cls(i,1)=atoi(argv[ii++]);
      cls(i,2)=atoi(argv[ii++]);
      cls(i,3)=atoi(argv[ii++]);
    }
  }


  void do_switch_changes(imatrix& cls,dvar_len_fish_stock_history& fsh)
  {
    //    cout << "B1" << endl;
    int mmin=cls.rowmin();
    int mmax=cls.rowmax();
    for (int i=mmin;i<=mmax;i++)
    {
      //    cout << "B2" << endl;
      switch (cls(i,1))
      {
      case 1:
        if (fsh.pmsd && cls(i,4)==1)
        {
          if (cls(i,5)>fsh.pmsd->num_species || cls(i,5)< 1)
          {
            cerr << " cls(i,5) num_species out of range "  
                 << cls(i,5) << endl; 
            ad_exit(1);
          }
          fsh.pmsd->parest_flags(cls(i,5),cls(i,2))=cls(i,3);
        }
        else
        {
          fsh.parest_flags(cls(i,2))=cls(i,3);
        }   
        break;
      case 2:
        if (fsh.pmsd && cls(i,4)==1)
        {
          if (cls(i,5)>fsh.pmsd->num_species || cls(i,5)< 1)
          {
            cerr << " cls(i,5) num_species out of range "  
                 << cls(i,5) << endl; 
            ad_exit(1);
          }
          fsh.pmsd->age_flags(cls(i,5),cls(i,2))=cls(i,3);
        }
        else
        {
          fsh.age_flags(cls(i,2))=cls(i,3);
        }   
        break;
      case 5:
        fsh.season_flags(1,cls(i,2))=cls(i,3);
        break;
      default:
        if (cls(i,1) < 0)
        {
          int snum=-cls(i,1);
          if (snum <= fsh.num_fisheries)
          {
            if (cls(i,2)!=3 || cls(i,3) !=0)
              fsh.fish_flags(snum,cls(i,2))=cls(i,3);
            else
              fsh.fish_flags(snum,cls(i,2))=fsh.nage-1;
          }
          else if (snum == 999)
          {
            if (cls(i,2)!=3 || cls(i,3) !=0)
            {
              for (int jj=1;jj<=fsh.num_fisheries;jj++)
              {
                fsh.fish_flags(jj,cls(i,2))=cls(i,3);
              }
            }
            else
            {
              for (int jj=1;jj<=fsh.num_fisheries;jj++)
              {
                fsh.fish_flags(jj,cls(i,2))=fsh.nage-1;
              }
            }
          }
          else if (snum == 9999)
          {
            if (fsh.pmsd==0)
            {
              for (int jj=1;jj<=fsh.num_tag_releases;jj++)
              {
                fsh.tag_flags(jj,cls(i,2))=cls(i,3);
              }
            }
            else
            {
              for (int it=1;it<=fsh.old_num_tag_releases;it++)
              {
                fsh.true_tag_flags(it,cls(i,2))=cls(i,3);    //NMD_5Nov13
              }
	      /*  //NMD_16Dec2013
              ivector iind(1,fsh.pmsd->num_species);
              iind.initialize();
              for (int it=1;it<=fsh.old_num_tag_releases;it++)
              {
                for (int is=1;is<=fsh.pmsd->num_species;is++)
                {
                  if (fsh.pmsd->tag_species_flag(it,is))
                  {
                    int ioff=fsh.offset(is)+iind(is)+1;
                    fsh.tag_flags(ioff,cls(i,2))
                      =fsh.true_tag_flags(it,cls(i,2));
                    iind(is)++;
                  }
                }
              }
	      */
            }
          }
          else if (snum>99999)
          {
            int reg=snum-99999;
            if (reg>fsh.num_regions)
              cerr << "Illegal command line option for region_flags"<< endl;
            else
              fsh.region_flags(reg,cls(i,2))=cls(i,3);
          }
          else if (snum>9999)
          {
            int it=snum-9999;
            if (fsh.pmsd==0)
            {
              if (it>fsh.num_tag_releases)
                cerr << "Illegal command line option for tag_flags"<< endl;
              else
                fsh.tag_flags(it,cls(i,2))=cls(i,3);
            }
            else
            {
              if (it>fsh.old_num_tag_releases)
              {
                cerr << "Illegal command line option for tag_flags"<< endl;
              }
              else
              {
// NMD_7Dec2013
                fsh.true_tag_flags(it,cls(i,2))=cls(i,3);

		/*  //NMD_16Dec2013
                ivector iind(1,fsh.pmsd->num_species);
                iind.initialize();
                for (int itt=1;itt<=it;itt++)
                {
                  for (int is=1;is<=fsh.pmsd->num_species;is++)
                  {
                    if (fsh.pmsd->tag_species_flag(itt,is))
                    {
                      int ioff=fsh.offset(is)+iind(is)+1;
                      if(itt==it)
		      {
                        fsh.tag_flags(ioff,cls(i,2))
                          =fsh.true_tag_flags(it,cls(i,2));
                      }
                      iind(is)++;
                    }
                  }
                }
// NMD_7Dec2013
*/  //NMD_16Dec2013
              }
            }
          }
          else if (snum>999)
          {
            int reg=snum-999;
            if (reg>fsh.num_fisheries)
              cerr << "Illegal command line option for fish_flags"<< endl;
            else
              fsh.fish_flags(reg,cls(i,2))=cls(i,3);
          }
          else
          {
            cerr << "Unrecognized flag option in do_switch_changes"
                 << endl;
            ad_exit(1);
          }
        }
        else
        {
          cerr << "Unrecognized flag option in do_switch_changes"
               << endl;
          ad_exit(1);
        }
      }
    }
    // do this in case these switches have been changed in this routine
    fsh.sum_ff71_flag=sum(column(fsh.fish_flags,71));
    fsh.sum_ff48_flag=sum(column(fsh.fish_flags,48));
    fsh.sum_ff79_flag=sum(column(fsh.fish_flags,79));
  }
  dvar_vector dddfun(int n)
  {
    dvar_matrix v(1,3,1,n);
    dmatrix w(1,3,1,n);
    w.fill_seqadd(1,1);
    v=w;
    return v(3);
  }

  void getroot(adstring& adr,adstring& fdp)
  {
    int sz=fdp.size();
    adstring tmp=fdp;
    int i;
    for (i=sz;i>=1;i--)
    {
      if (fdp(i)=='/' || fdp(i)=='\\')
      {
	tmp=fdp(i+1,sz);
	break;
      }	   	
    }
    sz=tmp.size();
    for (i=1;i<=sz;i++)
    {
      if (tmp(i)=='.') break;
    }
    adr=tmp(1,i-1);
  }  
#undef HOME_VERSION

int ddcheck_flag=-1;

void deal_with_command_line(int argc,char * argv[],
  adstring& full_datafile_path,adstring& full_input_parfile_path,
  adstring& full_output_parfile_path)
{
  int on1;
  int nopt=0;
  mytestcin();
  if ( (on1=option_match(argc,argv,"-dd",nopt))>-1)
  {
    if (!nopt)
    {
      cerr << "Usage -dd option needs integer  -- ignored" << endl;
    }
    else
    {   
      cout <<"nnewlan, on1: " << on1 << endl;
      cout <<"nnewlan, argc: " << argc << endl;
      int jj=atoi(argv[on1+1]);
      ddcheck_flag=jj;
    }
  }
  mytestcin();
  int nnopt=0;
  if (option_match(ad_argc,ad_argv,"-makepar",nnopt)>0)
  {
    makeparflag=1;
    full_datafile_path= argv[1];
    full_input_parfile_path = argv[2];
    full_output_parfile_path=argv[3];
    if ( argc < 4 ) 
    {
      cerr << "Not enough parameters in call to analyzer" << endl;
      cerr << "Usage  progname datafile inparfile outparfile [-switch options]"
       << endl;
      exit(1);
    }
    return;
  }

  mytestcin();
  if ((nopt=option_match(ad_argc,ad_argv,"-file",nnopt))>0)
  {
    full_datafile_path= argv[1];
    full_input_parfile_path = argv[2];
    full_output_parfile_path=argv[3];
	#ifdef __MSVC32__
		read_option_file((adstring)argv[nopt+1],cl_switches);
	#else
		read_option_file(argv[nopt+1],cl_switches);
	#endif
  mytestcin();
  }
  else
  {
  mytestcin();
    if ( argc < 4 ) 
    {
      cerr << "Not enough parameters in call to analyzer" << endl;
      cerr << "Usage  progname datafile inparfile outparfile [-switch options]"
         << endl;
      exit(1);
    }
    if (argc>=4) 
    {
      full_datafile_path= argv[1];
      full_input_parfile_path = argv[2];
      full_output_parfile_path=argv[3];
    }
  mytestcin();
    if (argc>=5)
    {
      int nopt=0;
      int opos=option_match(ad_argc,ad_argv,"-switch",nopt);
      if (opos>0)
      {
        int nswitch=atoi(argv[opos+1]);
        cl_switches.allocate(1,nswitch,1,5);
        set_command_line_switches(cl_switches,argc,argv,opos);
      }
    }
  }
}
  
void set_gtradient_structure_globals(void)
{
  long int gs_size=0;
  long int cp_size=0;
# if (defined(_PRINT_WD_))
  system("pwd > c:\\the_wd");
#endif
  ifstream ifs("mfcl.cfg",ios::in);
  if (!(!ifs))
  {
    ifs >> gs_size;
    if (gs_size)
    {
      cout << " GRADSTACK_BUFFER_SIZE read in and set equal to " << gs_size << endl;
      gradient_structure::set_GRADSTACK_BUFFER_SIZE(gs_size);
    }
    else
    {
      gradient_structure::set_GRADSTACK_BUFFER_SIZE(1000000);
    }

    ifs >> cp_size;
    if (cp_size)
    {
      cout << " CMPDIF_BUFFER_SIZE read in and set equal to " << cp_size
        << endl;
      gradient_structure::set_CMPDIF_BUFFER_SIZE(cp_size);
    }
    else
    {
      gradient_structure::set_CMPDIF_BUFFER_SIZE(56000000L);
    }
    ifs >> ad_array_mbl_size;
    if (ad_array_mbl_size)
    {
      cout << " array memblock size read in and set equal to " 
        << ad_array_mbl_size << endl;
    }
  }
  else
  {
    gradient_structure::set_GRADSTACK_BUFFER_SIZE(1000);
    gradient_structure::set_CMPDIF_BUFFER_SIZE(5600L);
    //gradient_structure::set_GRADSTACK_BUFFER_SIZE(1000000);
    //gradient_structure::set_CMPDIF_BUFFER_SIZE(56000000L);
  }
  gradient_structure::set_MAX_NVAR_OFFSET(90000);   //NMD 1June2015
  //gradient_structure::set_CMPDIF_BUFFER_SIZE(56000000L);
  gradient_structure::set_NO_DERIVATIVES();
  gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(50000);
}
void datafilereaderror(const char * s)
{
  cerr << "error reading " << s << " from frq file" << endl;
  exit(1);
}

void set_region_incidence_matrix(int nswitch,cifstream& infile,
  ivector& dataswitch,imatrix& Dflags,int& num_fisheries,int& num_regions,
  ivector& fishery_regions,int& tag_group_flag,dvector& _region_area,int & _data_file_version_no,pmulti_species_data & pmsd)
{
  // check for new format frq file 
  int num_species=1;
  if (!nswitch) 
  {
    cerr << "does this still work for nswitch ==0" << endl;
    exit(1);
  }
  else
  {
    dataswitch.allocate(1,nswitch);
    infile >> dataswitch;
    if (!infile) datafilereaderror("dataswitch");
    if (nswitch>7)
      _data_file_version_no=dataswitch(10);
    if (dataswitch(1)) num_regions=dataswitch(1);
    if (dataswitch(2)) num_fisheries=dataswitch(2);
    cout << nswitch << endl;
    if (nswitch>=4)
    {
      if (dataswitch(6))
        num_species=dataswitch(6);
      else
        num_species=1;
      //if (_data_file_version_no<7 || num_species==1)
        if (dataswitch(4)) tag_group_flag=dataswitch(4);
    }

    if (nswitch>=7)
    {
      if (dataswitch(4)) age_data_flag=dataswitch(7);
    }

    _region_area.allocate(1,num_regions);
    fishery_regions.allocate(1,num_fisheries);

    if (num_species>1)
    {
      // create structure to hold multi species imnformation
      process_multi_species_data(pmsd,infile,num_species,num_regions,
        fishery_regions);
    }
    infile >> _region_area;
    infile >> fishery_regions;
    // read in region connectivity description
    if (num_species<=1)
    {
      if (nswitch>=3)
      {
        if (dataswitch(3)) 
        {
          Dflags.allocate(1,num_regions,1,num_regions);
          for (int i=1;i<=num_regions;i++)
          { 
            Dflags(i,i)=0;
            for (int j=i+1;j<=num_regions;j++)
            {
              infile >> Dflags(i,j);
              Dflags(j,i)=Dflags(i,j);
            }
          }
        }
      }
    }
    else
    {
      if (nswitch>=3)
      {
        if (dataswitch(3)) 
        {
          Dflags.allocate(1,pmsd->num_kludged_regions);
          for (int ii=1;ii<=num_species;ii++)
          {
            for (int k=1;k<=pmsd->num_region_by_species(ii);k++)
            {
              Dflags(k+pmsd->rowoffset(ii)).allocate(1,pmsd->num_region_by_species(ii));
            }
            imatrix tmp=Dflags.sub(1+pmsd->rowoffset(ii),
              pmsd->num_region_by_species(ii)+pmsd->rowoffset(ii));
            tmp.rowshift(1);
            for (int i=1;i<=pmsd->num_region_by_species(ii);i++)
            { 
              tmp(i,i)=0;
              for (int j=i+1;j<=pmsd->num_region_by_species(ii);j++)
              {
                infile >> tmp(i,j);
                tmp(j,i)=tmp(i,j);
              }
            }
          }
        }
      }
    }
    
    if (!infile)
    {
      cerr << "Error trying to read in initial fishery and region"
               " configuration data" << endl;
      exit(1);
    }
  }
  if (!Dflags) 
  {
    Dflags.allocate(1,1,1,1);
    Dflags(1,1)=1;
  }
}

void read_move_weeks(ivector& move_weeks,ivector& dataswitch,
  cifstream & infile,int month_1)
{
  int nmove;
  int data_file_version_no = dataswitch(10);
  if (data_file_version_no>=4)
  {
    infile >> nmove;
    if (!infile) datafilereaderror("nmove");
    if (allocated(move_weeks)) move_weeks.deallocate();
    move_weeks.allocate(1,nmove);
    infile >> move_weeks;
    if (!infile) datafilereaderror("move_weeks");
  }
  else
  {
    nmove=1;
    if (allocated(move_weeks)) move_weeks.deallocate();
    move_weeks.allocate(1,nmove);
    move_weeks(1)=4*(month_1-1)+1;
  }
}

imatrix read_datafile_flaginfo(cifstream& infile,
  int num_fish,int data_file_version_no,imatrix& _season_region_flags,
  pmulti_species_data & pmsd)
{
  imatrix dff(1,5,1,num_fish);
  if (_data_file_version_no > 0)
    infile >> dff;
  else
    dff.initialize();
  if (!infile) {
    cerr << "error reading in data fish flags" << endl;
    exit(1);
  }

  if (_data_file_version_no > 4)
  {
    infile >> _season_region_flags;
    if (pmsd && pmsd->num_species>1)
    {
      ivector ub=pmsd->num_region_by_species;
      infile >> pmsd->ses_reg_recr_flags;
    }
  }
  else
  {
    _season_region_flags=1;
  }

  if (!infile) {
    cerr << "error reading in season_region_flags" << endl;
    ad_exit(1);
  }
  return dff;
}

int set_month1_flag(ivector& dataswitch)
{
  int m=1;
  int& ds9=dataswitch(9);
  if (ds9)
  {
    if (ds9<1 || ds9>12)
    {
      cerr << "Illegal month 1 value of " << ds9 << " ignored" << endl;
    }
    else
    {
      m=ds9;
    }
  }
  return m;
}

void dvar_len_fish_stock_history::setup_some_stuff(MY_DOUBLE_TYPE tmp_len_wt_coff,
  MY_DOUBLE_TYPE _wshlen,MY_DOUBLE_TYPE _wfilen,int _nwint,int _num_fisheries,imatrix _dff,
  par_cifstream * _pinfile,int _month_1,int _first_time,int _mfactor,
  ivector & _fishery_regions,int _first_data_month, 
  pmulti_species_data & pmsd )
{
  length_month_yr.allocate(1,_mfactor,1,nage);
  sdevs_yr.allocate(1,_mfactor,1,nage);
  first_data_month=_first_data_month;
  fishery_regions=_fishery_regions;
  if (_data_file_version_no >=3) 
    len_wt_coff=0.0;
  else
    len_wt_coff=tmp_len_wt_coff;
  nwint=_nwint;
  wfilen=_wfilen;
  wshlen=_wshlen;
  if (!pmsd)
    data_fish_flags.allocate(1,5,1,_num_fisheries);
  else
    data_fish_flags.allocate(1,5,1,pmsd->num_species*_num_fisheries);
  age_flags.initialize();
  if (!pmsd)
    data_fish_flags=_dff;
  else
  {
    for (int i=1;i<=5;i++)
    {
      for (int j=1;j<=pmsd->num_species;j++)
      {
        int offset=(j-1)*_num_fisheries;
        data_fish_flags(i)(1+offset,_num_fisheries+offset).shift(1)=_dff(i);
      }
    }
  }
  if (!  parest_flags(197)) 
  {
    *_pinfile >>   age_flags;
    if (pmsd)
    {
      if (pmsd->num_species>1)
        *_pinfile >>   pmsd->age_flags;
    }
  }
  check_number(*_pinfile,871,parest_flags);
  month_1=_month_1;
  first_time=_first_time;
  month_factor=_mfactor;
}
    

void dvar_fish_stock_history::set_year_flags_for_recruitment(void)
{
  // Determine what season grouping the first freq record is in
 
  int month_offset;
  
  if (month_1<=first_data_month)
    month_offset=first_data_month-month_1+1;
  else
    month_offset=first_data_month-month_1+13;

  int season1=(month_offset-1)/(12/num_seasons);
  ivector& sf=season_flags(1);
  ivector& yf=year_flags(1);
  seasonal_recruitment_flag=0;
  if (sum(sf) !=0)
  {
    seasonal_recruitment_flag=1;
    for (int i=1;i<=nyears;i++)
    {
      int imod=(season1+i-1)%num_seasons+1;
      yf(i)=sf(imod);
    }
  }
}

void dvar_len_fish_stock_history::set_effective_length_and_weight_sizes(void)
{
  test_fpu();
  int i;
  ivector w=column(fish_flags,49);
  effective_len_size=10;
  test_fpu();
  for (i=1;i<=num_fisheries;i++)
  {
  test_fpu();
    if (w(i) !=0)  effective_len_size(i)=w(i);
  test_fpu();
  }
  test_fpu();

  w=column(fish_flags,50);
  effective_weight_size=10;
  test_fpu();
  for (i=1;i<=num_fisheries;i++)
  {
    if (w(i) !=0)  effective_weight_size(i)=w(i);
  }
  if (nwint)
  {
    if (!allocated(wmid)) wmid.allocate(1,nwint);
    if (!allocated(realwmid)) realwmid.allocate(1,nwint);
    realwmid.fill_seqadd(wshlen+.5*wfilen,wfilen);
    wmid=pow(realwmid/value(sv(27)),1.0/value(sv(28)));
    if (pmsd)
      pmsd->wmid=get_multi_wmid();
#if !defined(NO_MY_DOUBLE_TYPE)
    wm2=pow(wmid,value(sv(28)-1.0L));
#else
    wm2=pow(wmid,value(sv(28)-1.0));
#endif
    if (pmsd)
      pmsd->wm2=get_multi_wm2();
  }
}

int check_working_directory_option(int argc,char *argv[])
{
  int nopt,nnopt;
  if ((nopt=option_match(ad_argc,ad_argv,"-wd",nnopt))>0)
  {
    chdir(argv[nopt+1]);
    return 1;
  }
  return 0;
}

#if defined(USE_ADPVM)
void mfsend_x_to_slaves(const dvar_vector&  x)
{
  // *********  begin variable send block  *************
  for (int i=mf_pvm->minspawn; i<=mf_pvm->maxspawn; i++)
  {
    int bufid = adpvm_master_vinitsend( PvmDataDefault );
    adpvm_pack(x);
    adpvm_master_vsend((mf_pvm->ad_stid)(i));
  }
  // *********  end variable send block  *************
}

void mfsend_x_to_slaves(const dvector&  x)
{
  // *********  begin variable send block  *************
  for (int i=mf_pvm->minspawn; i<=mf_pvm->maxspawn; i++)
  {
    int bufid = adpvm_master_cinitsend( PvmDataDefault );
    adpvm_pack(x);
    adpvm_master_csend((mf_pvm->ad_stid)(i));
  }
  // *********  end variable send block  *************
}

void mfsend_x_to_slaves(const dmatrix&  x)
{
  // *********  begin variable send block  *************
  for (int i=mf_pvm->minspawn; i<=mf_pvm->maxspawn; i++)
  {
    int bufid = adpvm_master_cinitsend( PvmDataDefault );
    adpvm_pack(x);
    adpvm_master_csend((mf_pvm->ad_stid)(i));
  }
  // *********  end variable send block  *************
}

void mfsend_x_to_slaves(const imatrix&  x)
{
  // *********  begin variable send block  *************
  for (int i=mf_pvm->minspawn; i<=mf_pvm->maxspawn; i++)
  {
    int bufid = adpvm_master_cinitsend( PvmDataDefault );
    adpvm_pack(x);
    adpvm_master_csend((mf_pvm->ad_stid)(i));
  }
  // *********  end variable send block  *************
}
dvector mfget_dvector_from_slave(int i)
{
  dvector tmp;
  // *********  begin variable send block  *************
  {
    // *********  begin variable receive block  *************
    adpvm_master_crecv((mf_pvm->ad_stid)(i)); //get the values from slave
    adpvm_unpack(tmp);
    adpvm_master_end_creceive();  // have got all the values
  }
  // *********  end variable receive block  *************
  return tmp;
}


dvar_vector mfget_f_from_slaves(void)
{
 // cout << "calling heapcheck" << endl;
  //cout << "heapcheck = " << heapcheck() << endl;
  dvar_vector fslave(mf_pvm->minspawn,mf_pvm->maxspawn);
  // *********  begin variable send block  *************
  for (int i=mf_pvm->minspawn; i<=mf_pvm->maxspawn; i++)
  {
    {
      // *********  begin variable receive block  *************
      adpvm_master_vrecv((mf_pvm->ad_stid)(i)); //get the values from slave
      adpvm_unpack(fslave(i));
      adpvm_master_end_vreceive();  // have got all the values
    }
  }
  // *********  end variable receive block  *************
  //cout << "calling heapcheck" << endl;
  //cout << "heapcheck = " << heapcheck() << endl;
  return fslave;
}

void  mfget_dv3_sum_from_slaves(const dvar3_array& _w)
{
  dvar3_array& w=(dvar3_array&)(_w);
  // *********  begin variable send block  *************
  for (int i=mf_pvm->minspawn; i<=mf_pvm->maxspawn; i++)
  {
    dvar3_array v;
    {
      // *********  begin variable receive block  *************
      adpvm_master_vrecv((mf_pvm->ad_stid)(i)); //get the values from slave
      adpvm_unpack(v);
      adpvm_master_end_vreceive();  // have got all the values
    }
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (mmin !=v.indexmin() || mmax !=v.indexmax() )
    {
      cerr << "Incompatible dvar3_arrays in "
        "void  mfget_dv3_sum_from_slaves(const dvar3_array& _w)"
        << endl;
      ad_exit(1);
    }
    for (int i=mmin;i<=mmax;i++) w(i)+=v(i);
  }
  // *********  end variable receive block  *************
}

void send_taggroup_assignments(const ivector& min_group,
  const ivector& max_group,dvar_len_fish_stock_history& fsh)
{
  for (int i=mf_pvm->minspawn; i<=mf_pvm->maxspawn; i++)
  {
    int bufid = adpvm_master_cinitsend( PvmDataDefault );
    adpvm_pack(min_group(i));
    adpvm_pack(max_group(i));
    adpvm_master_csend((mf_pvm->ad_stid)(i));
  }
  // *********  end constant send block  *************
}

void get_taggroup_assignments(dvar_len_fish_stock_history& fsh)
{
  int ptid=pvm_parent();
  // ***************  begin constant receive block *********************
  adpvm_slave_crecv(ptid);
  adpvm_unpack(fsh.min_tag_group);  
  adpvm_unpack(fsh.max_tag_group);  
  adpvm_slave_end_creceive();
  // ***************  end receive block ***********************************
}

void mfsend_file_to_slaves(const adstring& s)
{
  const int BUFSIZE=4096;
  adstring zipped_file = s + ".zip"; 
  adstring command_string = " zip " + zipped_file + " " + s;
#if defined(__BORLANDC__)
  mf_system((char*)(command_string));
#else
  system((char*)(command_string));
#endif
  FILE * file=fopen((char*) zipped_file,"rb");
  char buffer[BUFSIZE]; 
  do
  {
    int nread=fread(buffer,1,BUFSIZE,file);
    if (nread<BUFSIZE)
    {
      if (feof(file)!=0)
      {
         cout << "end of file reached OK in " << s << endl;
      }
      else if (ferror(file))
      {
        cout << "Error reading file " << s << " in mfsend_file_to_slaves"
             << endl;
        ad_exit(1);
      }
    }
    mfsend_int_to_slaves(nread);
  
    if (nread>0)
    {
      mfsend_bytes_to_slaves(buffer,nread);
    }
    else
    {
      break;
    }
  }
  while(1);   
  fclose(file);
  file=0;
}  

void mfget_bytes_from_master(char * buffer,int nbytes)
{
  int ptid=pvm_parent();
  // ***************  begin constant receive block *********************
  int atag,atid,alen;
  pvm_precv(ptid,71,buffer,nbytes,PVM_BYTE,&atid,&atag,&alen);
 /*
  adpvm_slave_crecv(ptid);
  adpvm_unpack(buffer,nbytes);  
  adpvm_slave_end_creceive();
  */
  // ***************  end receive block ***********************************
}

void mfget_file_from_master(const adstring& s)
{
  const int BUFSIZE=4096;
  adstring zipped_file = s + ".zip"; 
  FILE * file=fopen((char*)zipped_file,"wb");
  char buffer[BUFSIZE]; 
  do
  {
    int nread=mfget_int_from_master();
  
    if (nread>0)
    {
      mfget_bytes_from_master(buffer,nread);
      int nr=fwrite(buffer,1,nread,file);
    }
    else
    {
      break;
    }
  }
  while(1);
  fclose(file);
  adstring command_string = " unzip -o " + zipped_file;
#if defined(__BORLANDC__)
  mf_system((char*)(command_string));
#else
  system((char*)(command_string));
#endif
  file=0;
}  

void mfsend_bytes_to_slaves(char * buffer,int n)
{
  for (int i=mf_pvm->minspawn; i<=mf_pvm->maxspawn; i++)
  {
    pvm_psend((mf_pvm->ad_stid)(i),71,buffer,n,PVM_BYTE);
    
    /*
    int bufid = adpvm_master_cinitsend( PvmDataDefault );
    adpvm_pack(buffer,n);
    adpvm_master_csend((mf_pvm->ad_stid)(i));
    */
  }
  // *********  end constant send block  *************
}

  
void mfsend_int_to_slaves(int n)
{
  for (int i=mf_pvm->minspawn; i<=mf_pvm->maxspawn; i++)
  {
    //int bufid = adpvm_master_cinitsend( PvmDataDefault );
    int bufid = adpvm_master_cinitsend( PvmDataDefault );
    adpvm_pack(n);
    adpvm_master_csend((mf_pvm->ad_stid)(i));
  }
  // *********  end constant send block  *************
}

imatrix mfget_imatrix_from_master(void)
{
  int ptid=pvm_parent();
  // ***************  begin constant receive block *********************
  adpvm_slave_crecv(ptid);
  imatrix m;
  adpvm_unpack(m);  
  adpvm_slave_end_creceive();
  // ***************  end receive block ***********************************
  return m;
}
int mfget_int_from_master(void)
{
  int ptid=pvm_parent();
  // ***************  begin constant receive block *********************
  adpvm_slave_crecv(ptid);
  int n;
  adpvm_unpack(n);  
  adpvm_slave_end_creceive();
  // ***************  end receive block ***********************************
  return n;
}
#endif //#if defined(USE_ADPVM)

void dvar_fish_stock_history::read_movement_coffs(adstring & movement_coffs_file)
{
  cifstream cif(movement_coffs_file);

  if (!cif)
  {
    cerr << "Error opening movement coffs file" << endl;
    ad_exit(1);
  }
  cif >> Dad;
  if (!cif)
  {
    cerr << "Error readin movement coffs from file" << endl;
    ad_exit(1);
  }
  int mmin=Dad.indexmin();
  int mmax=Dad.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    int jmin=Dad(i).indexmin();
    int jmax=Dad(i).indexmax();
    for (int j=jmin;j<=jmax;j++)
    {
      if (!pmsd)
      {
        int kmin=Dad(i,j).indexmin();
        int kmax=Dad(i,j).indexmax();
        Dad(i,j)=-1.0*Dad(i,j);
        int k;
        for (k=kmin;k<=kmax;k++)
        {
          Dad(i,j)(k,k)=0.0;
        }
        dvector tmp=value(colsum(Dad(i,j)));
        for (k=kmin;k<=kmax;k++)
        {
          Dad(i,j)(k,k)=1-tmp(k);
        }
        Dad(i,j)=inv(Dad(i,j));
      }
      else
      {
        int ns=pmsd->num_species;
	int nregs=pmsd->num_real_regions;
        for (int is=1;is<=ns;is++)
        {
          int tmpmin=((is-1)*nregs)+1;
          int tmpmax=tmpmin+nregs-1;
	  dmatrix tmp_Dad(1,nregs,1,nregs);
	  tmp_Dad.initialize();
	  tmp_Dad=value(Dad(i,j).sub(tmpmin,tmpmax).shift(1));
          int kmin=tmp_Dad.indexmin();
          int kmax=tmp_Dad.indexmax();
          tmp_Dad=-1.0*tmp_Dad;
          int k;
          for (k=kmin;k<=kmax;k++)
          {
            tmp_Dad(k,k)=0.0;
          }
          dvector tmp=colsum(tmp_Dad);
          for (k=kmin;k<=kmax;k++)
          {
            tmp_Dad(k,k)=1-tmp(k);
          }
          tmp_Dad=inv(tmp_Dad);
	  int icnt;
	  icnt=1;
          for (k=tmpmin;k<=tmpmax;k++)
          {
	    Dad(i,j,k)=tmp_Dad(icnt);
	    icnt++;
          }
        }
      }


    }
  }
  cout <<"nnewlan, Dad(1,2): " << Dad(1,2) << endl;
  cout <<"nnewlan,colsum(Dad(1,2)): " << colsum(Dad(1,2)) << endl;
}
  
void  dvar_len_fish_stock_history::resize_if_necessary(void)
{
  if (parest_flags(155)>0)
  {
    if (allocated(recr_polys))
    {
      recr_polys.deallocate();
    }
    int ntr=(nyears-1)/num_seasons+1;
    recr_polys.allocate(1,ntr,0,parest_flags(155));
    recr_polys=orthpoly(ntr,parest_flags(155));
  }
}

void colp(const dvar_matrix& m,int col)
{
  cout <<"nnewlan, column(m,col): " << column(m,col) << endl;
}

void process_multi_species_data(pmulti_species_data & pmsd, 
  cifstream& infile,int num_species,int num_regions,
  ivector& fishery_regions)
{
  // create structure to hold multi species information
  pmsd=new multi_species_data(num_species,num_regions);
  // read in the number of tag groups for each species
  infile >> pmsd->num_tag_groups;
  infile >> pmsd->species_by_region;
  cout <<  pmsd->species_by_region << endl;
  cout << rowsum(pmsd->species_by_region) << endl;
  pmsd->num_region_by_species=rowsum(pmsd->species_by_region);
  pmsd->num_kludged_regions=sum(pmsd->num_region_by_species);
  pmsd->region_species_pointer.allocate(1,pmsd->num_kludged_regions);
  pmsd->tag_species_pointer.allocate(1,sum(pmsd->num_tag_groups));
  
  int i,j;
  int ii=0;
  for (i=1;i<=num_species;i++)
  { 
    for (j=1;j<=pmsd->num_tag_groups(i);j++)
    { 
      pmsd->tag_species_pointer(++ii)=i;
    }
  }
  ii=0;
  for (i=1;i<=num_species;i++)
  {
    for (int j=1;j<=pmsd->num_region_by_species(i);j++)
    {
      pmsd->region_species_pointer(++ii)=i;
    }
  }
  cout << pmsd->region_species_pointer << endl;
  pmsd->region_bounds(1,1)=1;
  pmsd->region_bounds(1,2)=pmsd->num_region_by_species(1);
  pmsd->rowoffset(1)=0;
  for (i=2;i<=num_species;i++)
  {
    pmsd->region_bounds(i,1)=pmsd->region_bounds(i-1,2)+1;
    pmsd->region_bounds(i,2)=pmsd->region_bounds(i,1)
      +pmsd->num_region_by_species(i)-1;
    pmsd->rowoffset(i)=pmsd->rowoffset(i-1)+pmsd->num_region_by_species(i-1);
  }
}
void allocate_species(fishery_freq_record_array&  flfra,int num_species)
{
  int mmin=flfra.indexmin();
  int mmax=flfra.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    flfra[i].species.allocate(1,num_species);
    flfra[i].species1.allocate(1,num_species);
    flfra[i].species2.allocate(1,num_species);
    flfra[i].species3.allocate(1,num_species);
  }
}

void dvar_fish_stock_history::check_flag_stuff(void)
{
  if (age_flags(94)==1 && (parest_flags(373)>0  
    && parest_flags(374)>0 ))
  {
    cerr << " if age_flags(94)==1 then parest_flags(373)"
        " parest_flags(374) be zero" << endl;
    ad_exit(1);
  }
  if (age_flags(94)==3 && age_flags(95) )
  {
    cerr << " if age_flags(94)==3 then age_flags(95) must be zero" << endl;
    ad_exit(1);
  }
}


int search_for_string(cifstream & in,
                             const string& s,const string& es)
{
  string line;
  int posfound=0;
  size_t pos = 0;
  while (getline(in, line))
  {
    pos = line.find(es);
    if (pos != string::npos) 
    {
      break;
    }
    pos = line.find(s);
    if (pos != string::npos) 
    {
      posfound=1;
      break;
    }
  }
  return posfound;
}

int getsize(const dvar_matrix & _m)
{
  ADUNCONST(dvar_matrix,m)
  int ns=0;
  if (allocated(m))
  {
    int mmin=m.indexmin();
    int mmax=m.indexmax();
    for (int i=mmin;i<=mmax;i++)
    {
      ns+=getsize(m(i));
    }
    return ns;
  }
  else
    return 0;
}


int getsize(const dvar3_array & _m)
{
  ADUNCONST(dvar3_array,m)
  int ns=0;
  if (allocated(m))
  {
    int mmin=m.indexmin();
    int mmax=m.indexmax();
    for (int i=mmin;i<=mmax;i++)
    {
      ns+=getsize(m(i));
    }
    return ns;
  }
  else
    return 0;
}



dvar_fish_stock_history::~dvar_fish_stock_history()
{
  if (ppstf)
  {
    delete ppstf;
    ppstf=0;
  }
}      

void mfcl_error(void)
{
  cerr << "Check this bfore using" << endl;
  ad_exit(1);
}

void df_print_identifier_string_30(void)
{
  adstring str=get_string_marker_30();
  cout << "GS = " << str << endl;
  //gradlog << "GS = " << str << endl;
}

void insert_identifier_string_30(const char * _s)
{
  char * s = (char *) _s;
  int n=strlen(s);
  char ss[31];
  ss[30]='\0';
  int n1=min(n,30);
  for (int i=0;i<n1;i++)
  {
    ss[i]=_s[i];
  }
  for (int i=n1+1;i<30;i++)
  {
    ss[i]=' ';
  }
  //save_identifier_string(ss);
  save_identifier_string_dbg_and_opt(ss);
}

void df_grad_profiler_30(void);
void set_grad_profiler_30(const char * s)
{
  insert_identifier_string_30(s);

  gradient_structure::GRAD_STACK1->
            set_gradient_stack(df_grad_profiler_30);
  //extern adtimer * profile_timer_1;
  //extern adtimer * profile_timer_2;
}

void df_print_identifier_string(void);
void df_grad_profiler_30(void)
{
  //verify_identifer_string("fff1");
  //adstring get_identifier_string(10);
  df_print_identifier_string_30();
  if (profile_timer_1)
  {
    cout << profile_timer_1->get_elapsed_time_and_reset()/1000. << endl;
    cout << profile_timer_2->get_elapsed_time()/1000. << endl;
  }
}
    

adstring get_string_marker_30(void)
{
  adstring str1;
  AD_LONG_INT num_bytes=30;
  char str[31];
  str[num_bytes]='\0';
  gradient_structure::get_fp()->fread(str,num_bytes);
  //clogf << "in verify_id_string " << str1 << endl;
  str1=str;
 return str1;
}

 
int save_identifier_string_dbg_and_opt(char* str)
{
 int wsize=sizeof(char);
 int length=strlen(str);
 gradient_structure::get_fp()->fwrite(str,length);
 return 0;
}

ostream& operator << (ostream& s,const hemholtz_contrast & _H)
{
  ADUNCONST(hemholtz_contrast,H)
  s << H.get_C();
  return s;
}

  hemholtz_contrast::hemholtz_contrast(int n) 
     // Hemholtz contrast, I think;
  {
    allocate(n); // Hemholtz contrast, I think;
  }


  void hemholtz_contrast::allocate(int n)  // Hemholtz contrast, I think;
  {
    if (!allocated(C))
    {
      C.allocate(1,n,1,n);
    }
    else
    {
      cerr << "C already allocated" << endl;
      ad_exit(1);
    }
 
    C.initialize();
    C(1)=1.0;
    C(1)/=norm(C(1));
    for (int i=2;i<=n;i++)
    {
      for (int j=1;j<i;j++)
      {
        C(i,j)=1.0;
      }
      C(i,i)= -double(i-1);
      C(i)/=norm(C(i));
    }
    cout << setw(6) << setprecision(3) << setfixed() << C << endl;
    dmatrix M(1,n,1,n);
    M.initialize();
    for (int i=1;i<=n;i++)
    {
      M(i,i)=C(i)*C(i);
      for (int j=1;j<i;j++)
      {
        M(i,j)=C(i)*C(j);
        M(j,i)=M(i,j);
      }
    }
    cout << setw(6) << setprecision(3) << setfixed() << M << endl;
    C=trans(C);
  }

ostream& operator << (ostream& s,const grouped_df_orthogonal_contrast & _H)
{
  ADUNCONST(grouped_df_orthogonal_contrast,H)
  s << H.get_ungrouped_v();
  return s;
}

ostream& operator << (ostream& s,const df_orthogonal_contrast & _H)
{
  ADUNCONST(df_orthogonal_contrast,H)
  s << H.get_v();
  return s;
}

cifstream& operator >> (cifstream& cif,const grouped_df_orthogonal_contrast & _H)
{
  ADUNCONST(grouped_df_orthogonal_contrast,H)
  cif >> H.get_ungrouped_v();
  return cif;
}

cifstream& operator >> (cifstream& cif,const df_orthogonal_contrast & _H)
{
  ADUNCONST(df_orthogonal_contrast,H)
  cif >> H.get_v();
  return cif;
}

df_orthogonal_contrast::df_orthogonal_contrast(int n) 
{
  allocate(n); // Hemholtz contrast, I think;
}

grouped_df_orthogonal_contrast::grouped_df_orthogonal_contrast(int n) :
  ungrouped_v(1,n) 
{

}

grouped_df_orthogonal_contrast::grouped_df_orthogonal_contrast(ivector& v ) :
  ungrouped_v(v.indexmin(),v.indexmax()), gmgr(v)
{ 
  df_orthogonal_contrast::allocate(get_ngroups()); 
}

int df_orthogonal_contrast::nvcal(void) 
{
  return v.indexmax()-v.indexmin()+1;
}

void df_orthogonal_contrast::xinit(dvector& x,int & ii) 
{
  int n=C.indexmax()-C.indexmin()+1;
  dvector xc=value(v)*C;    // C is orthog so this is inv(C)*v
  x(ii,ii+n-1).shift(1)=xc; 
  ii+=n;
}

void grouped_df_orthogonal_contrast::xinit(dvector& x,int & ii, 
  ivector& groups) 
{
  if (!allocated())   // check for partial allocation of ungrouped vector
  {                   // of parmeters
    if (!::allocated(ungrouped_v))
    {
      cerr << "Error in xinit ungouped vector not allocated" << endl;
      ad_exit(1);
    }
    // late allocation of the group manager
    gmgr.allocate(groups);
    int ngroups=gmgr.get_ngroups();
    df_orthogonal_contrast::allocate(ngroups);
    if (!::allocated(v))
    {
      v.allocate(1,ngroups);
    }
  }
  int ngroups=gmgr.get_ngroups();
  for (int i=1;i<=ngroups;i++)
  {
    get_v()(i)=ungrouped_v(gmgr.get_grouped_vector_ptr()(i));
  }
  df_orthogonal_contrast::xinit(x,ii);
}

int grouped_df_orthogonal_contrast::nvcal(ivector& groups) 
{
  // late allocation of the group manager
  if (!gmgr.allocated())
  {
    gmgr.allocate(groups);
  }
  int ngroups=gmgr.get_ngroups();
  df_orthogonal_contrast::allocate(ngroups);
  return df_orthogonal_contrast::nvcal(); 
}

void df_orthogonal_contrast::reset(dvar_vector& x,int & ii) 
{
  int n=C.indexmax()-C.indexmin()+1;
  v=C*x(ii,ii+n-1).shift(1);   
  ii+=n;
}

void grouped_df_orthogonal_contrast::reset(dvar_vector& x,int & ii) 
{
  int n=C.indexmax()-C.indexmin()+1;
  v=C*x(ii,ii+n-1).shift(1);         // there are n groups and v holds the par
  imatrix groups=gmgr.get_groups();  // value for each group
  ivector ng=gmgr.get_ng();
  for (int i=1;i<=n;i++)
  {
    for (int j=1;j<=ng(i);j++)
    {
      ungrouped_v(groups(i,j))=v(i);  // ungrouped_v holds the ungrouped values
    }
  }
  ii+=n;
}


void df_orthogonal_contrast::allocate(int n)  
{
  if (!allocated(v))
  {
    v.allocate(1,n);
  }
  else
  {
    cerr << "v already allocated" << endl;
    ad_exit(1);
  }
  if (!allocated(C))
  {
    C.allocate(1,n,1,n);
  }
  else
  {
    cerr << "C already allocated" << endl;
    ad_exit(1);
  }
  C(1)=1.0;
  C(2).fill_seqadd(-1.0,2.0/(n-1.0));
  for (int i=2;i<n;i++)     // make unnormalized orthogonal polynomials
  {
    C(i+1)=elem_prod(C(2),C(i));
  }
  for (int i=1;i<=n;i++)    // modified gram-schmidt
  {
    C(i)/=norm(C(i));
    for (int j=i+1;j<=n;j++)
    {
      C(j)-=(C(j)*C(i))*C(i);
    }
  }
  C=trans(C);
}

void hem_test(int n)
{
  hemholtz_contrast C(n);

  dvector x(1,n);

  x.initialize();
  x(1)=1.0;
  cout << C << endl;

  cout << C*x << endl;

  x.initialize();
  x(2)=1.0;
  cout << C*x << endl;


  ad_exit(1);
}    

void make_dave_contrast(int n);
void make_hemholtz_contrast(int n);
void orth_test(int _n)
{
  make_hemholtz_contrast(200);

  ad_exit(1);


  df_orthogonal_contrast CC(200);
  {
    ofstream ofs("contrast");
    dmatrix C=value(CC.get_C());
    ofs << setw(10) << setprecision(3) << setscientific() << trans(C)*C << endl;
  }

  cifstream cif("orthtest");
  int nn=0;
  cif >> nn;
  ivector v(1,40);
  v.initialize();
  v(1)=5;
  for (int i=2;i<=40;i++)
  {
    v(i)=(i-1)%10+1;
    v(i)+=v(i-1);
    v(i)=v(i)/2;
  }
  int vsize=v.indexmax()-v.indexmin()+1;
  grouped_df_orthogonal_contrast C(vsize);
  cif >> C;
  cout << C << endl;
  int n=C.nvcal(v);
  dvector x(1,n);
  int ii=1;
  C.xinit(x,ii,v);
  cout << x << endl;
  dvar_vector vx=x;
  ii=1;
  C.reset(vx,ii);
  cout << C << endl;
  

  x.initialize();
  x(1)=1.0;
  cout << C << endl;

  cout << C*x << endl;

  x.initialize();
  x(2)=1.0;
  cout << C*x << endl;


  ad_exit(1);
}    

  void new_group_manager::allocate(ivector& v)
  {
    int mmin=v.indexmin();
    int mmax=v.indexmax();
    imatrix tmp(1,2,mmin,mmax);
    tmp(1).fill_seqadd(1,1);
    tmp(2)=v;
    imatrix tmp2=trans(sort(trans(tmp),2));
    ivector ind=tmp2(1);
    ivector w=tmp2(2);
    ngroups=1;
    for (int i=mmin+1;i<=mmax;i++)
    {
      if (w(i)>w(i-1)) ngroups++;
    }  
    ng.allocate(1,ngroups);
    ng.initialize();
    int ii=1;
    ng(ii)=1;
    for (int i=mmin+1;i<=mmax;i++)
    {
      if (w(i)>w(i-1)) ii++;
      ng(ii)+=1;
    }  
    groups.allocate(1,ngroups,1,ng);
    ii=1;
    int jj=0; 
    groups(ii,++jj)=ind(mmin);
    for (int i=mmin+1;i<=mmax;i++)
    {
      if (w(i)>w(i-1)) 
      {
        ii++;
        jj=0;
      }
      groups(ii,++jj)=ind(i);
    }  
    if (!::allocated(grouped_vector_ptr))
    {
      grouped_vector_ptr.allocate(1,ngroups); // picks an ungrouped vector 
                                              // value for each group
    }
    for (int i=1;i<=ngroups;i++)
    {
      grouped_vector_ptr(i)=min(groups(i));
    }
  }


  void make_dave_contrast_line(dvector v)
  {
    int n=v.indexmax()-v.indexmin()+1;
    int n2=n/2;
    int offset=v.indexmin()-1;
    v(1+offset,n2+offset)=1.0;
    if (n%2==0)
    {
      v(1+n2+offset,n+offset)=-1.0;
    }
    else
    {
      v(1+n2+offset,n+offset)=-1.0*sqrt(n2/(n2+1.0));
    }
  }   

  void make_dave_contrast(int n)
  {
    dmatrix C(1,n,1,n);
    C.initialize();
    C(1)=1.0;
    make_dave_contrast_line(C(2));

    int n2=n/2;
    make_dave_contrast_line(C(3)(1,n2));
    make_dave_contrast_line(C(4)(n2+1,n));

    int n4=n/4;
    make_dave_contrast_line(C(5)(1,n4));
    make_dave_contrast_line(C(6)(n4+1,2*n4));
    make_dave_contrast_line(C(7)(2*n4+1,3*n4));
    make_dave_contrast_line(C(8)(3*n4+1,n));



    cout << setw(6) << setfixed() << setprecision(3) << C << endl;
  }
    
  void make_hemholtz_contrast(int n)  // Hemholtz contrast, I think;
  {
    dmatrix C;
    C.allocate(1,n,1,n);
 
    C.initialize();
    C(1)=1.0;
    C(1)/=norm(C(1));
    for (int i=2;i<=n;i++)
    {
      for (int j=1;j<i;j++)
      {
        C(i,j)=1.0;
      }
      C(i,i)= -double(i-1);
      C(i)/=norm(C(i));
    }
    //cout << setw(6) << setprecision(3) << setfixed() << C << endl;
    dmatrix M(1,n,1,n);
    M.initialize();
    for (int i=1;i<=n;i++)
    {
      M(i,i)=C(i)*C(i);
      for (int j=1;j<i;j++)
      {
        M(i,j)=C(i)*C(j);
        M(j,i)=M(i,j);
      }
    }
    cout << setw(9) << setprecision(6) << setfixed() << "Orth Test" << endl;
    cout << setw(9) << setprecision(6) << setfixed() << M << endl;
    C=trans(C);
  }
ostream & operator << (const ostream& _s,
  const fishery_header_record& hr)
{
  ADUNCONST(ostream,s)
  s << "year " << hr.year << " movement_period " << hr.movement_period
    << " month " << hr.month << " week " << hr.week << endl;
  return s;
}

void  mp(fishery_freq_record_array& ffra,int i)
{
  cout << ffra(i);
}
void  mp(fishery_header_record_array& fhra,int i)
{
  cout << fhra(i);
}
