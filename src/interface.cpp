/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include <fvar.hpp>
//dvector dd_newton_raphson2(dvector& _g){}

dvector dd_newton_raphson2(dvector& _g,banded_symmetric_dmatrix& _H,dvector& _x)
{
  int bw=_H.bandwidth();
  int n=_g.indexmax();
  banded_symmetric_dmatrix H(1,n,bw);
  dvector x(1,n);
  dvector g(1,n);

  int i;
  for (i=1;i<=n;i++)
  {
    x(i)=_x(i);
    g(i)=_g(i);
  }
    
  for (i=0;i<=bw-1;i++)
  {
    H(i)=_H(i);
  }
  dvector y=solve(H,g);
  dvector w=x-y;
  return w;
}

#if (defined(linux) || __MSVC32__>=8)
#  define _export  
#endif

extern "C" _export  void dd_newton_raphson(int n,MY_DOUBLE_TYPE * v,MY_DOUBLE_TYPE * diag,
    MY_DOUBLE_TYPE * ldiag, MY_DOUBLE_TYPE * yy)
 {
   banded_symmetric_dmatrix H(1,n,2);
   dvector x(1,n);

   MY_DOUBLE_TYPE * w=v;
   MY_DOUBLE_TYPE * z=yy;
   MY_DOUBLE_TYPE * d=diag;
   MY_DOUBLE_TYPE * ld=ldiag;
   w--;
   z--;
   d--;
   ld--;
   
   int i;
   for (i=1;i<=n;i++)
   {
     x(i)=w[i];
     H(i,i)=d[i];
   }
   for (i=1;i<n;i++)
   {
     H(i+1,i)=ld[i];
   }
   dvector y=solve(H,x);
   for (i=1;i<=n;i++)
   {
     z[i]=y(i);
   }
}
