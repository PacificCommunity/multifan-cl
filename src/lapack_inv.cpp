/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#include "all.hpp"
#if defined(NO_MY_DOUBLE_TYPE)
  #define lapack_complex_float std::complex<float>
  #define lapack_complex_double std::complex<double>
#endif
#include <complex>
#include <lapacke.h>
lapack_int matInv(double *A, unsigned n);

dmatrix lapack_luinv(const dmatrix& M)
{
  int n=M.indexmax();
  dvector v;
  dmatrix M1=make_dmatrix(1,n,1,n,v);
  M1=M;
  int iret=matInv(&(v[1]),n);
  return M1;
}




lapack_int matInv(double *A, unsigned n)
{
  cerr << "Starting" << endl;
    ivector ipiv(1,n);
    ipiv.initialize();
    //int ipiv[n+1];
    lapack_int ret;

    ret =  LAPACKE_dgetrf(LAPACK_ROW_MAJOR,
                          n,
                          n,
                          A,
                          n,
                          &(ipiv[1]));

    if (ret !=0)
    {
      cerr << "Error" << endl;
      return ret;
    }

    ret = LAPACKE_dgetri(LAPACK_ROW_MAJOR,
                       n,
                       A,
                       n,
                       &(ipiv[1]));
    if (ret !=0)
    {
      cerr << "Error 2" << ret << endl;
      return ret;
    }
    return ret;
}

