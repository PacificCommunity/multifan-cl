/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#if !defined(__TRIDIAG__)
#define __TRIDIAG__
  dvector tridag(const dmatrix & S,dvector & r);
  class  symmetric_tridiagonal_dmatrix 
  {
    int mmin;
    int mmax;
    dvector diag;
    dvector subdiag;
  public:
    void initialize() { diag.initialize(); subdiag.initialize();}
    int indexmin() { return mmin;}
    int indexmax() { return mmax;}
    symmetric_tridiagonal_dmatrix(int _mmin,int _mmax) : mmin(_mmin),
      mmax(_mmax),diag(_mmin,_mmax),subdiag(_mmin+1,_mmax){;}
    MY_DOUBLE_TYPE & get_diag(int i) { return diag(i); }
    MY_DOUBLE_TYPE & get_subdiag(int i) { return subdiag(i); }
    dvector & get_diag(void) { return diag; }
    dvector & get_subdiag(void) { return subdiag; }
    MY_DOUBLE_TYPE & operator () (int i,int j) 
    {
      if (i==j)
       return diag(i);
      else if (i==(j+1))
      {
        return subdiag(i);
      }
      else if ((i+1)==j)
       return subdiag(j);
      else
      {
        cerr << "index error in symmetric_tridiagonal_dmatrix ()"
             << endl;
        ad_exit(1);
      }
    }
         
  };
  dvector solve(symmetric_tridiagonal_dmatrix& S,const dvector& r);
  dvector operator * (const symmetric_tridiagonal_dmatrix& S,const dvector& r);
  dvector operator * (const dvector& r,const symmetric_tridiagonal_dmatrix& S);
  
  MY_DOUBLE_TYPE ln_det_choleski(const symmetric_tridiagonal_dmatrix& S,const int& ierr);

dmatrix make_dmatrix(const symmetric_tridiagonal_dmatrix& S);

dvector solve(const symmetric_tridiagonal_dmatrix& _S,const dmatrix& u,
  const dmatrix& v,const dvector& x);
MY_DOUBLE_TYPE ln_det(const symmetric_tridiagonal_dmatrix& _S,const dmatrix& u,
  const dmatrix& v,const int& ierr);

#endif  //if !defined(__TRIDIAG__)
