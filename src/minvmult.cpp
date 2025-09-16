/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include <admodel.h>
#include <cblas.h>
#include "svd.h"


MY_DOUBLE_TYPE * getpointer(const dvector& _v)
{
  ADUNCONST(dvector,v)
  return &(v(v.indexmin()));
}
dvector minvmult(svd& SVD,dvector& G,int n)
{
#if !defined(USE_NO_LAPACKE)
  //double *s = &(SVD.get_S()(1));
  //double *u = &(SVD.get_u()(1));
  int m=SVD.get_n();
  int k=m;
  std::unique_ptr<MY_REAL_DOUBLE[]>UC=
    std::unique_ptr<MY_REAL_DOUBLE[]>(new MY_REAL_DOUBLE[n]);
  MY_REAL_DOUBLE * C = UC.get();
  MY_REAL_DOUBLE* A=SVD.get_pvt(void)
  MY_REAL_DOUBLE*B get_double_pointer(G);
  //MY_DOUBLE_TYPE *A = getpointer(SVD.get_vt());
  //MY_DOUBLE_TYPE *C = getpointer(c);
  //MY_DOUBLE_TYPE *B = getpointer(G);
  MY_REAL_DOUBLE_TYPE alpha=1.0;
  MY_REAL_DOUBLE beta=0.0;
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                m, n, k, alpha,A, k, B, n, beta, C, n);
  dvector c(1,m*n);
  for (int i=1;i<=n*m;i++)
  {
    c(i)=C(i-1);
  }
  return c;
#else
  dvector c(1,5);
  c=0.0;
   cerr <<  "Need to implement this for long MY_DOUBLE_TYPE " << endl;
   ad_exit(1);
  return c;
#endif
}
dmatrix inverse2(svd& SVD)
{
#if !defined(USE_NO_LAPACKE)
  int n=SVD.get_n();
  int m=n;
  int k=m;
  dvector c(1,m*n);
  dmatrix CC(1,m);
  int offset=0;
  for (int i=1;i<=m;i++)
  {
    CC(i)=c(1+offset,n+offset).shift(1);
    offset+=n;
  }
  dmatrix VT=SVD.get_VT();
  dvector S=SVD.get_S();
  for (int i=1;i<=n;i++)
  {
    VT(i)*=(1.0/S(i));
  }
  MY_DOUBLE_TYPE *B = getpointer(SVD.get_vt());
  MY_DOUBLE_TYPE *C = getpointer(c);
  MY_DOUBLE_TYPE *A = getpointer(SVD.get_u());
  MY_DOUBLE_TYPE alpha=1.0;
  MY_DOUBLE_TYPE beta=0.0;
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                m, n, k, alpha,A, k, B, n, beta, C, n);
  return CC;
#else
  dmatrix CC(1,5);
  CC=0.0;
   cerr <<  "Need to implement this for long MY_DOUBLE_TYPE " << endl;
   ad_exit(1);
  return CC;
#endif
}
dmatrix multiply(dmatrix & AA,dmatrix & BB)
{
//#if !defined(USE_NO_LAPACKE)
  int m=AA.indexmax();
  int k=AA(1).indexmax();
  int k1=BB.indexmax();
  int n=BB(1).indexmax();
  //dvector a(1,m*k);
  MY_REAL_DOUBLE *a = new MY_REAL_DOUBLE[m*k];
  int ii=0;
  for (int i=1;i<=m;i++)
  {
    for (int j=1;j<=k;j++)
    {
      a[ii++]=(double)(AA(i,j));
    }
  }
 
  MY_REAL_DOUBLE * b = new MY_REAL_DOUBLE[m*k];
  //dvector b(1,k*n);
  ii=0;
  for (int i=1;i<=k;i++)
  {
    for (int j=1;j<=n;j++)
    {
      b[ii++]=double(BB(i,j));
    }
  }
 
  MY_REAL_DOUBLE * c = new MY_REAL_DOUBLE[m*n];
  //MY_DOUBLE_TYPE *A = getpointer(a);
  //MY_DOUBLE_TYPE *B = getpointer(b);
  //MY_DOUBLE_TYPE *C = getpointer(c);
  MY_REAL_DOUBLE alpha=1.0;
  MY_REAL_DOUBLE beta=0.0;
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                m, n, k, alpha,a, k, b, n, beta, c, n);
  dmatrix CC(1,m,1,n);
  ii=0;
  for (int i=1;i<=m;i++)
  {
    for (int j=1;j<=n;j++)
    {
      CC(i,j)=c[ii++];
    }
  }
  //cout << c << endl;
  return CC;
//#else
//  dmatrix CC(1,5);
//  CC=0.0;
//   cerr <<  "Need to implement this for long MY_DOUBLE_TYPE " << endl;
//   ad_exit(1);
//  return CC;
//#endif
}
