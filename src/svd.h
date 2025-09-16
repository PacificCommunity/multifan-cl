/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#if !defined(__SVD_H__)
# define __SVD_H__
#include <admodel.h>
//#pragma message "MY_DOUBLE_TYPE=" BOOST_PP_STRINGIZE(MY_DOUBLE_TYPE)

MY_REAL_DOUBLE* get_double_pointer(const dvector& _v);

class svd
{
  int n;
  int m;
  dmatrix U;
  dvector S;
  dmatrix VT;
  dvector u;
  dvector vt;
public:
  svd(int _n)  
  {
    n=_n;
    m=_n;
    u.allocate(1,n*n);
    vt.allocate(1,m*m);
    S.allocate(1,max(n,m));
    U.allocate(1,n);
    int offset=0;
    for (int i=1;i<=n;i++)
    {
      U(i)=u(1+offset,n+offset).shift(1);
      offset+=n;
    }
    VT.allocate(1,m);
    offset=0;
    for (int i=1;i<=m;i++)
    {
      VT(i)=vt(1+offset,m+offset).shift(1);
      offset+=m;
    }
  }
  int & get_n(void){ return n;}
  int & get_m(void){ return m;}
  MY_REAL_DOUBLE* get_pu(void)
  { 
# if (defined(BOOST_MP_USE_FLOAT128))
    return get_double_pointer(u);
# else
    return &(u[1]);
# endif
  }
  MY_REAL_DOUBLE* get_pvt(void)
  { 
#  if (defined(BOOST_MP_USE_FLOAT128))
    return get_double_pointer(vt);
# else
    return &(vt[1]);
# endif
  }
  MY_REAL_DOUBLE* get_pS(void)
  { 
#  if (defined(BOOST_MP_USE_FLOAT128))
    return get_double_pointer(S);
# else
    return &(S[1]);
# endif
  }
  dvector & get_u(void){ return u;}
  dvector & get_vt(void){ return vt;}
  dvector & get_S(void){ return S;}
  dmatrix & get_U(void){ return U; }
  dmatrix & get_VT(void){ return VT; }
  dmatrix mult(void)
  {
    dmatrix _U=get_U();
    dmatrix _VT=get_VT();
    dvector _S=get_S();
    dmatrix SS(1,n,1,n);
    SS.initialize();
    for (int i=1;i<=n;i++)
    {
      SS(i,i)=_S(i);
    }
    cout << _U*SS*_VT << endl << endl;
    dmatrix M(1,n,1,n);
    return M;
  }
};
void lapack_symmetric_eigen(dmatrix & MM,dmatrix& E,dvector& e,int& ierr);
dmatrix lapack_choleski_inverse(dmatrix& M,int & ierr);
svd svd_decomp(dmatrix& M);
int Lapack_Choleski_Inverse(int n,int lda,double *a);
#endif
