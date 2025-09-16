/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#include "all.hpp"


#if defined(NO_MY_DOUBLE_TYPE)
#  if defined(MY_DOUBLE_TYPE)
#    undef MY_DOUBLE_TYPE
#  endif
#  define MY_DOUBLE_TYPE double
#  if defined(MY_REAL_DOUBLE)
#    undef MY_REAL_DOUBLE double
#  endif
#  define MY_REAL_DOUBLE double
#endif

#if defined(NO_MY_DOUBLE_TYPE)
  #define lapack_complex_float std::complex<float>
  #define lapack_complex_double std::complex<double>
#endif
#include <complex>
#include <memory>
#include <lapacke.h>
void print_identifier_stuff(ofstream& ofs);

extern "C" {
 void goto_set_num_threads(int num_threads);
}

#include "svd.h"
  
  MY_REAL_DOUBLE* get_double_pointer(const dvector& _v)
  {
    ADUNCONST(dvector,v)
    if (!allocated(v))
    {
      return  (MY_REAL_DOUBLE *)(0);
    }
    int mmin=v.indexmin();
    int mmax=v.indexmax();
    int n=mmax-mmin+1;
    if (n<=0)
    {
      return  (MY_REAL_DOUBLE *)(0);
    }
    MY_REAL_DOUBLE * pv=new MY_REAL_DOUBLE[n];
    MY_DOUBLE_TYPE * dv = &(v(mmin));  
    for (int i=0;i<n;i++)
    {
      pv[i] = static_cast<MY_REAL_DOUBLE>(dv[i]);
    }
    return pv;
  }
    
  MY_REAL_DOUBLE * get_double_pointer(const dmatrix& M)
  {
    if (!allocated(M))
    {
      return  (MY_REAL_DOUBLE *)(0);
    }
    int mmin=M.indexmin();
    int mmax=M.indexmax();
    int nc=mmax-mmin+1;
    if (nc<=0)
    {
      return  (MY_REAL_DOUBLE *)(0);
    }
    ivector lb(mmin,mmax);
    ivector ub(mmin,mmax);
    int n=0;
    for (int i=mmin;i<=mmax;i++)
    {
      if (!allocated(M(i)))
      {
        lb(i)=0;
        ub(i)=-1;
      }
      else
      {
        lb(i)=M(i).indexmin();
        ub(i)=M(i).indexmax();
        n+=ub(i)-lb(i)+1;
      }
    }
    MY_REAL_DOUBLE * pv=new MY_REAL_DOUBLE[n];
    int ii=0;
    for (int i=mmin;i<=mmax;i++)
    {
      for (int j=lb(i);j<=ub(i);j++)
      {
        pv[ii++] = static_cast<MY_REAL_DOUBLE>(M(i,j));
      }
    }
    return pv;
  }
    

  void double_pointer_values_return(dvector& v,MY_REAL_DOUBLE * pv)
  {
    if (allocated(v))
    {
      int mmin=v.indexmin();
      int mmax=v.indexmax();
      int n=mmax-mmin+1;
      if (n>0)
      {
        for (int i=0;i<n;i++)
        {
          v[i+mmin] = static_cast<MY_DOUBLE_TYPE>(pv[i]);
        }
      }
    }
  }

    
  void double_pointer_values_return(dmatrix& M,MY_REAL_DOUBLE * pv)
  {
    if (allocated(M))
    {
      int mmin=M.indexmin();
      int mmax=M.indexmax();
      int nc=mmax-mmin+1;
      if (nc>0)
      {
        ivector lb(mmin,mmax);
        ivector ub(mmin,mmax);
        int n=0;
        for (int i=mmin;i<=mmax;i++)
        {
          if (!allocated(M(i)))
          {
            lb(i)=0;
            ub(i)=-1;
          }
          else
          {
            lb(i)=M(i).indexmin();
            ub(i)=M(i).indexmax();
            n+=ub(i)-lb(i)+1;
          }
        }
        int ii=0;
        for (int i=mmin;i<=mmax;i++)
        {
          for (int j=lb(i);j<=ub(i);j++)
          {
            M(i,j) = static_cast<MY_DOUBLE_TYPE>(pv[ii++]);
          }
        }
      }
    }
  }


dmatrix make_dmatrix(int cmin,int cmax,int rmin,int rmax,dvector & v);


void double_pointer_values_get(const dmatrix &M,
  std::unique_ptr<MY_REAL_DOUBLE[]>& up) 
{
  int mmin=M.colmin();
  int mmax=M.colmax();
  ivector lb(mmin,mmax);
  ivector ub(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    lb=M(i).indexmin();
    ub=M(i).indexmax();
  }
  MY_REAL_DOUBLE * p = up.get();
  if (p==0)
  {
    int n=0;
    for (int i=mmin;i<=mmax;i++)
    {
      n+=max(ub(i)-lb(i)+1,0);
    }
    up=std::unique_ptr<MY_REAL_DOUBLE[]>(new MY_REAL_DOUBLE[n]);
  }

  int ii=0;
  for (int i=mmin;i<=mmax;i++)
  {
    for (int j=lb(i);j<=ub(i);j++)
    {
      up[ii++]=static_cast<double>(M(i,j));
    }
  }
}  

std::unique_ptr<MY_REAL_DOUBLE[]> double_pointer_values_get(const dmatrix &M)
{
  int mmin=M.colmin();
  int mmax=M.colmax();
  ivector lb(mmin,mmax);
  ivector ub(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    lb=M(i).indexmin();
    ub=M(i).indexmax();
  }
  int n=0;
  for (int i=mmin;i<=mmax;i++)
  {
    n+=max(ub(i)-lb(i)+1,0);
  }
  std::unique_ptr<MY_REAL_DOUBLE[]> up(new MY_REAL_DOUBLE[n]);

  int ii=0;
  for (int i=mmin;i<=mmax;i++)
  {
    for (int j=lb(i);j<=ub(i);j++)
    {
      up[ii++]=static_cast<double>(M(i,j));
    }
  }
  return up;
}  

void double_pointer_values_return(dmatrix &M,
  std::unique_ptr<MY_REAL_DOUBLE[]>& up) 
{
  int mmin=M.colmin();
  int mmax=M.colmax();
  ivector lb(mmin,mmax);
  ivector ub(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    lb=M(i).indexmin();
    ub=M(i).indexmax();
  }
  MY_REAL_DOUBLE * p = up.get();
  if (p==0)
  {
    cerr << "This can not happen" << endl;
    ad_exit(1);
  }

  int ii=0;
  for (int i=mmin;i<=mmax;i++)
  {
    for (int j=lb(i);j<=ub(i);j++)
    {
#if !defined(NO_MY_DOUBLE_TYPE)
      //M(i,j)=static_cast<boost::multiprecision::float128>(up[ii++]);
#else
      //M(i,j)=static_cast<boost::multiprecision::double>(up[ii++]);
#endif
      M(i,j)=static_cast<MY_DOUBLE_TYPE>(up[ii++]);
    }
  }
}  

void double_pointer_values_return(dvector &v,
  std::unique_ptr<MY_REAL_DOUBLE[]>& up) 
{
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  MY_REAL_DOUBLE * p = up.get();
  if (p==0)
  {
    cerr << "This can not happen" << endl;
    ad_exit(1);
  }

  int ii=0;
  for (int i=mmin;i<=mmax;i++)
  {
#if !defined(NO_MY_DOUBLE_TYPE)
    //M(i,j)=static_cast<boost::multiprecision::float128>(up[ii++]);
#else
    //M(i,j)=static_cast<boost::multiprecision::double>(up[ii++]);
#endif
    v(i)=static_cast<MY_DOUBLE_TYPE>(up[ii++]);
  }
}  


void lapack_symmetric_eigen(dmatrix & _MM,dmatrix& eigenvectors,
  dvector& eigenvalues,int & ierr)
{
  int n=_MM.indexmax();
  int lda=n;
  dmatrix MM(1,n,1,n);
  MM=_MM;
  if (!allocated(eigenvectors))
  {
    eigenvectors.allocate(1,n,1,n);
  }
  if (!allocated(eigenvalues))
  {
    eigenvalues.allocate(1,n);
  }
  //std::unique_ptr<MY_REAL_DOUBLE[]> pa (new MY_REAL_DOUBLE[n*n]);
  std::unique_ptr<MY_REAL_DOUBLE[]> pe (new MY_REAL_DOUBLE[n]);

  /*
  std::unique_ptr<MY_REAL_DOUBLE[]> pa;
  double_pointer_values_get(MM,pa);
  */

  std::unique_ptr<MY_REAL_DOUBLE[]> pa=double_pointer_values_get(MM);
  /*
  int ii=0;
  for (int i=1;i<=n;i++)
  {
    for (int j=1;j<=n;j++)
    {
      pa[ii++]=static_cast<double>(MM(i,j));
    }
  }
  */
  int ii=0;
  cout << endl;
  for (int i=1;i<=n;i++)
  {
    for (int j=1;j<=n;j++)
    {
      cout << pa[ii++] << "  ";
    }
    cout << endl;
  }
  cout << endl;
  cout << endl;
    
  int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n,pa.get(),n,pe.get());
  double_pointer_values_return(eigenvectors,pa);
  double_pointer_values_return(eigenvalues,pe);
}

dmatrix lapack_choleski_inverse(dmatrix& M,int & ierr)
{
  int n=M.indexmax();
  int n2=n*n;
  int lda=n;
  cout << "In lapack_choleski_inverse(dmatrix& M) n = " << n << endl;
  MY_REAL_DOUBLE * pa = get_double_pointer(M);
  //dvector a(1,n2);
  cout << "Calling Lapack_Choleski_Inverse" << endl;     
  ierr=Lapack_Choleski_Inverse(n,lda,pa); 
  dmatrix N(1,n,1,n);
  double_pointer_values_return(N,pa);
  return N;
}

ostream& operator << (ostream & os,const svd& _SVD)
{
  ADUNCONST(svd,SVD)
  os << SVD.get_n() << " " << SVD.get_m() << endl;
  os << SVD.get_U() << endl << endl;
  os << SVD.get_U() << endl << endl;
  os << SVD.get_VT()<< endl << endl;
  os << SVD.get_S();
  return os;
}
  

svd SingularValueDecomposition( int n,     // number of columns in matrix
                                int lda,   // leading dimension of matrix
                                double *a) // pointer to top-left corner
{
  int m=n;
    svd SVD(n);
    double *s= SVD.get_pS();
    double *u= SVD.get_pu();
    double *vt= SVD.get_pvt();
//#if !defined(USE_NO_LAPACKE)
    int info = LAPACKE_dgesdd( LAPACK_ROW_MAJOR, 'S', m, n, a, n, s,
            u, m, vt, n );
    if (info) // handle error conditions here
    {
      cout << "svd failed" << endl;
      ad_exit(1);
    }
//#else
//   cerr <<  "Need to implement this for long double " << endl;
//   ad_exit(1);
//#endif
    return SVD;
}
svd SingularValueDecomposition1( int n,     // number of columns in matrix
                                int lda,   // leading dimension of matrix
                                double *a) // pointer to top-left corner
{
  int m=n;

    svd SVD(n);
    //dvector S(1,n);
    //dvector U(1,n*n);
    //dvector VT(1,n*n);
    double *s= SVD.get_pS();
    double *u= SVD.get_pu();
    double *vt= SVD.get_pvt();
    //double *s =&(S(1)); 
    //double *u = &(U(1));
    //double *vt = &(VT(1));
#if !defined(USE_NO_LAPACKE)
    int info = LAPACKE_dgesdd( LAPACK_ROW_MAJOR, 'S', m, n, a, n, s,
            u, m, vt, n );
    if (info) // handle error conditions here
    {
      cout << "svd failed" << endl;
      ad_exit(1);
    }
#else
   cerr <<  "Need to implement this for long double " << endl;
   ad_exit(1);
#endif
    cout << SVD.get_S() << endl << endl;
    cout << SVD.get_U() << endl << endl;
    cout << SVD.get_VT() << endl << endl;
    return SVD;
}
dvector minvmult(svd& SVD,dvector& G,int n);
dvector makeexample(int n)
{
  dvector a(1,n*n);
  dmatrix A(1,n);
  int ii=0;
  for (int i=1;i<=n;i++)
  {
    A(i)=a(1+ii,n+ii).shift(1);
    ii+=n;
  }
  random_number_generator rng(679);
  A.fill_randn(rng);
  A=A*trans(A);
  //cout << A << endl;
  //cout << a << endl;
  return a;
}

dmatrix inverse(svd& SVD)
{
  int n=SVD.get_n();
  dvector S=SVD.get_S();
  dmatrix VTS(1,n,1,n);
  dmatrix D(1,n,1,n);
  D.initialize();
  for (int i=1;i<=n;i++)
  {
    VTS(i)=(1.0/S(i))*SVD.get_VT()(i);
    D(i,i)=S(i);
  }
  return SVD.get_U()*VTS;
}  
dmatrix multiply(dmatrix & AA,dmatrix & BB);


dmatrix lapack_lusolve(const dmatrix& M,const dmatrix& G)
{
  dmatrix Minv=lapack_luinv(M);
  //dmatrix C=Minv*G;
  dmatrix tmp=lapack_matrix_multiplcation(Minv,G);
  //cout << norm2(C-tmp) << endl;
  //ad_exit(1);
  return tmp;
}
/*
dmatrix svd_solve(dmatrix& M,dmatrix& G)
{
  int n=M.indexmax();
  int n1=G.indexmax();
  int m=n;
  int lda=n;
  cout << "In svd_solve n = " << n << endl;
  if (n !=n1 )
  {
    cerr << "size error" << endl;
    ad_exit(1);
  }
  //dvector a(1,n*n);
  dmatrix Msave(1,n,1,n);
  Msave=M;

  MY_REAL_DOUBLE * pa=get_double_pointer(M);
  cout << "Calling SingularValueDecomposition(" << endl;     
  svd SVD=SingularValueDecomposition(     // number of rows in matrix
    n,     // number of columns in matrix
    lda,   // leading dimension of matrix
    pa);    // pointer to top-left corner
  double_pointer_values_return(M,pa);
  cout << "Finished SingularValueDecomposition(" << endl;     
  dvector S=SVD.get_S();
  dmatrix VTS(1,n,1,n);
  dmatrix UU(1,n,1,n);
  dmatrix D(1,n,1,n);
  D.initialize();
  for (int i=1;i<=n;i++)
  {
    D(i,i)=1/(S(i));
  }
  dmatrix Test=D*SVD.get_VT();
  
  UU=SVD.get_U();
  cout << "Finished SingularValueDecomposition A" << endl;     
  for (int i=1;i<=n;i++)
  {
    VTS(i)=(1.0/S(i))*SVD.get_VT()(i);
  }
  {
    ofstream ofs("Sorted_ S");
    ofs << sort(S) << endl;
  }
  cout << "norm2(Test-VTS)" << endl;
  cout << norm2(Test-VTS) << endl;
  dmatrix NNN=SVD.get_U()*VTS;
  cout << "Finished SingularValueDecomposition B" << endl;     
  cout << "Starting Multiply" << endl;     

  dmatrix MMM=multiply(SVD.get_U(),VTS);
  cout << "norm2(NNN-MMM)" << endl;
  cout << norm2(NNN-MMM) << endl;
  dmatrix QQ=NNN*G;
  cout << "Finished Multiply" << endl;     
  uostream uos("hessinv");
  cout << "writing MMM" << endl;     
  uos << n << MMM;
  cout << "Starting second Multiply" << endl;     
  dmatrix PP=multiply(MMM,G);
  cout << "norm2(QQ-PP)" << endl;
  cout << norm2(QQ-PP) << endl;
  cout << "norm2(G-Msave*QQ)" << endl;
  cout << norm2(G-Msave*QQ) << endl;
  //ad_exit(1);
  cout << "finished second Multiply" << endl;     
  return PP;
}
*/


svd svd_decomp(dmatrix& M)
{
  int n=M.indexmax();
  int m=n;
  int lda=n;
  cout << "In svd_decomp n = " << n << endl;
  dvector a(1,n*n);
  int ii=1;
  for (int i=1;i<=n;i++)
  {
    for (int j=1;j<=n;j++)
    {
      a(ii++)=M(i,j);
    }
  }
  MY_REAL_DOUBLE * pa=get_double_pointer(a);
  cout << "Calling SingularValueDecomposition" << endl;     
  svd SVD=SingularValueDecomposition(n,lda,pa); 
  double_pointer_values_return(a,pa);
  cout << "Finished SingularValueDecomposition" << endl;     
  return SVD;
}
svd svd_decomp1(dmatrix& M)
{
  int n=M.indexmax();
  int m=n;
  int lda=n;
  cout << "In svd_decomp n = " << n << endl;
  dvector a(1,n*n);
  int ii=1;
  for (int i=1;i<=n;i++)
  {
    for (int j=1;j<=n;j++)
    {
      a(ii++)=M(i,j);
    }
  }
# if (defined(BOOST_MP_USE_FLOAT128))
  double * pa=get_double_pointer(a);
# else
    double * pa=&(a(1));
# endif
  //MY_REAL_DOUBLE * pa=get_double_pointer(a);
  cout << "Calling SingularValueDecomposition" << endl;     
  svd SVD=SingularValueDecomposition1(n,lda,pa); 
  double_pointer_values_return(a,pa);
  cout << "Finished SingularValueDecomposition" << endl;     
  return SVD;
}

void makedebugdata(int n,int m,int m2)
{
  random_number_generator rng(821);
  dmatrix H(1,n,1,n);
  H.fill_randn(rng);
  H=H*trans(H);
  uostream uos("debug.hes");
  uos << n;
  uos << H;

  dmatrix G(1,m,1,n);
  G.fill_randn(rng);
  uostream uos1("debug.dep");
  uos1 << n << m;
  uos1 << G;


  dmatrix G1(1,m2,1,n);
  dmatrix Gall(1,m+m2,1,n);
  G1.fill_randn(rng);
  uostream uos2("debug.dp2");
  uos2 << n << m2;
  uos2 << G1;

  for (int i=1;i<=m;i++)
  {
    Gall(i)=G(i);
  }

  for (int i=1;i<=m2;i++)
  {
    Gall(m+i)=G(i)-G1(i);
  }

  dmatrix Hinv=inv(H);
  dvector var(1,m+m2);
  for (int i=1;i<=m+m2;i++)
  {
    var(i)=Gall(i)*(Hinv*Gall(i));
  }
  cout << "True variances" << endl;
  cout << var << endl;
}

int fileSize(const char *add)
{
    ifstream mySource;
    mySource.open(add, ios_base::binary);
    mySource.seekg(0,ios_base::end);
    int size = mySource.tellg();
    mySource.close();
    return size;
}
dvector * global_v_pointer=0;

void calculate_variance_by_delta_method_noeff(const char * _s)
{
  adstring s=_s;
  goto_set_num_threads(4);
  int debugflag=0; 
  int n=0;
  adstring fname;
  adstring fname1;
  adstring fname2;
  if (debugflag==0) 
  {
    fname=s+".hes";
    fname1=s+".dep";
    fname2=s+".dp2";
  }
  else
  {
    makedebugdata(5,7,2);
    fname="debug.hes";
    fname1="debug.dep";
    fname2="debug.dp2";
  }
  uistream uis(fname);
  uis >> n;
  dmatrix M(1,n,1,n);
  uis >> M;
  if(!uis)
  {
    cerr << "error reading hessian from file " << fname << endl;
    ad_exit(1);
  }
  for (int i=1;i<=n;i++)
  {
    for (int j=1;j<=n;j++)
    {
      M(i,j)=0.5*(M(i,j)+M(j,i));
      M(j,i)=M(i,j);
    }
  }
  adtimer adt;
  int n1=0;
  int ndep=0;
  int n2=0;
  int ndep2=0;
  int nbytes1=fileSize(fname);
  uistream uis1(fname1);
  if (!uis1)
  {
    cerr << "Error trying to open file " << fname1 << endl;
  }
  uis1 >> n1 >> ndep;
  cout << 8*n1*ndep+8 << endl;
  int nbytes2=fileSize(fname2);
  uistream uis2(fname2);
  if (!uis2)
  {
    cerr << "Error trying to open file " << fname2 << endl;
  }
  uis2 >> n2 >> ndep2;
  cout << 8*n2*ndep2+8 << endl;
  if (!uis2) {
    cerr << "Error reading n2 from file " << fname1 << endl;
    ad_exit(1);
  }
  if (n1!=n || n2 !=n) {
    cerr << "number of independent variables is not consistent"
            " between files" << endl;
    ad_exit(1);
  }
  dmatrix tdepgrad(1,ndep+ndep2,1,n);
    
  {
    //read in as doubles and convert
    MY_REAL_DOUBLE tmp1;
    //char tmp[8];
    for (int i=1;i<=ndep;i++)
    {
      for (int j=1;j<=n;j++)
      {
        // horrrible cast to char * to read in
        uis1.read((char*)(&tmp1),8);
        tdepgrad(i,j)=tmp1;
      }
    }
  }
  if (!uis1) {
    cerr << "Error reading gradients from " << fname << endl;
    ad_exit(1);
  }
  {
    //read in as doubles and convert
    MY_REAL_DOUBLE tmp1;
    //char tmp[8];
    for (int i=1;i<=ndep2;i++)
    {
      for (int j=1;j<=n;j++)
      {
        // horrrible cast to char * to read in
        uis2.read((char*)(&tmp1),8);
        tdepgrad(i+ndep,j)=tmp1;
      }
    }
  }
 // uis2 >> tdepgrad.sub(1+ndep,ndep+ndep2);
  if (!uis2) {
    cerr << "Error reading gradients from " << fname2 << endl;
   ad_exit(1);
  }
  for (int i=1;i<=ndep2;i++)
  {
//    tdepgrad(ndep+i)-=tdepgrad(i);
    tdepgrad(ndep+i)=tdepgrad(i)-tdepgrad(ndep+i);  //NMD_31may2022
  }
  
  //fname="skj.var";
  fname= s + ".var";
  ofstream ofs(fname);
  //ofstream ofsz(fname1);

  cout << "calling print identifier stuff" << endl;
  print_identifier_stuff(ofs);
  //print_identifier_stuff(ofsz);
  cout << "finished print identifier stuff" << endl;

  cout << "starting trans -- elapsed time = " 
       << adt.get_elapsed_time_and_reset()/1000 
       << " seconds " << endl;
  dmatrix G=trans(tdepgrad);
    
  cout << "did trans -- elapsed time = " 
       << adt.get_elapsed_time_and_reset()/1000 
       << " seconds " << endl;

  cout << "starting lu solve -- elapsed time = " 
       << adt.get_elapsed_time_and_reset()/1000 
       << " seconds " << endl;
  dmatrix PP=lapack_lusolve(M,G);
    
  int nl=M.indexmin();
  int nu=M.indexmax();
  cout << "number of parameters " << n << " number of dependent variables "
       << ndep << endl;
  cout << " elapsed time = " << adt.get_elapsed_time_and_reset()/1000 
       << " seconds " << endl;
  cout << "calculating variances" << endl;
  dvector var(1,ndep+ndep2);
  for (int i=1;i<=ndep+ndep2;i++)
  {
    var(i)=tdepgrad(i)*column(PP,i);
  }
  
  /*
  {
    ofstream ofs1("lu_realvar");
    for (int i=1;i<=ndep+ndep2;i++)
    {
      ofs1 << var(i)  << endl;
    }
  }
  */
  

  ifstream ifs2("deplabel.tmp");
  if (!ifs2)
  {
    cerr << "Error trying to open file deplabel.tmp" << endl;
  }
  line_adstring deplabel;
  dvector depvalue(1,ndep);
  dvector ndepvalue(1,ndep2);

  adstring_array dstrings(1,ndep);
  for (int i=1;i<=ndep;i++)
  {
    ifs2 >> depvalue(i);
    ifs2 >> deplabel;
    dstrings(i)=deplabel;
    ofs << setw(3) << i; 
    //ofsz << setw(3) << i; 
    if (var(i)>=0)
    {
      ofs << "  "  << setw(11) << setprecision(8) << sqrt(var(i));
    }
    else
    {
      cout << "variance = " << var(i) << endl;
      ofs << "  ********";
    }
    ofs << "  "  << setw(9) << setprecision(4) 
         << depvalue(i) << "  "  << deplabel << endl;
  }



  ifstream nifs2("deplabel_noeff.tmp");
  if (!nifs2)
  {
    cerr << "Error trying to open file deplabel_noeff.tmp" << endl;
  }
  dmatrix ndepgrad(1,ndep2,1,n);
  adstring_array ndstrings(1,ndep2);
  for (int i=1;i<=ndep2;i++)
  {
    nifs2 >> ndepvalue(i);
    nifs2 >> deplabel;
    ndstrings(i)=deplabel;
  }
  for (int i=1;i<=ndep2;i++)
  {
    if (i>ndep) break;
    if (i>ndep2) break;
    ofs << setw(3) << i; 
    //ofsz << setw(3) << i; 
    if (var(ndep+i)>0)
    {
      ofs << "  "  << setscientific() << setw(11) << setprecision(8) << sqrt(var(ndep+i));
    }
    else
    {
      cout << "variance = " << var(ndep+i) << endl;
      ofs << "  ********";
    }
    ofs << "   "  << setscientific() << setw(13) << setprecision(8) 
         << depvalue(i) - ndepvalue(i) << "  "  
         << dstrings(i) << " - " << ndstrings(i) << " " << depvalue(i) << " " 
         << ndepvalue(i) << endl;
  }
  cout << " elapsed time = " << adt.get_elapsed_time_and_reset()/1000 
       << " seconds " << endl;
}

void calculate_variance_by_delta_method(const char * _s)
{
  adstring s=_s;
  int debugflag=0; 
  //if (argc>1)
  //  debugflag=1;
  goto_set_num_threads(4);
  int n=0;

  adstring fname;
  adstring fname1;
  if (debugflag==0) 
  {
    fname=s+".hes";
    fname1=s+".dep";
  }
  else
  {
    makedebugdata(5,7,2);
    fname="debug.hes";
    fname1="debug.dep";
  }
  
 
  uistream uis(fname);
  uis >> n;
  dvector a(1,n*n);
  uis >> a;
  if(!uis)
  {
    cerr << "error reading hessian from file " << fname << endl;
    exit(1);
  }
  
 
  /*
  n=100;
  dvector a=makeexample(n);
  */
  dmatrix M(1,n,1,n);
  int offset=0;
  for (int i=1;i<=n;i++)
  {
    M(i)=a(1+offset,n+offset).shift(1);
    offset+=n;
  }
  for (int i=1;i<=n;i++)
  {
    for (int j=1;j<=n;j++)
    {
      M(i,j)=0.5*(M(i,j)+M(j,i));
      M(j,i)=M(i,j);
    }
  }
  adtimer adt;
  int n1=0;
  int ndep=0;
  uistream uis1(fname1);
  if (!uis1)
  {
    cerr << "Error trying to open file " << fname1 << endl;
  }
  uis1 >> n1 >> ndep;
  cout << 8*n1*ndep+8 << endl;
  if (n1!=n) {
    cerr << "number of independent variables is not consistent"
            " between files" << endl;
    exit(1);
  }
  dmatrix tdepgrad(1,ndep,1,n);
  uis1 >> tdepgrad.sub(1,ndep);
  if (!uis1) {
    cerr << "Error reading gradients from " << fname << endl;
    exit(1);
  }
  //fname="skj.var";
  fname= s + ".var";
  ofstream ofs(fname);

  cout << "calling print identifiter stuff" << endl;
  print_identifier_stuff(ofs);
  cout << "finished print identifiter stuff" << endl;

  cout << "starting trans -- elapsed time = " 
       << adt.get_elapsed_time_and_reset()/1000 
       << " seconds " << endl;
  dmatrix G=trans(tdepgrad);
    
  cout << "did trans -- elapsed time = " 
       << adt.get_elapsed_time_and_reset()/1000 
       << " seconds " << endl;

  //cout << "starting svd solve -- elapsed time = " 
  cout << "starting lapack_lusolve  -- elapsed time = " 
       << adt.get_elapsed_time_and_reset()/1000 
       << " seconds " << endl;
  dmatrix PP=lapack_lusolve(M,G);
  //dmatrix PP=svd_solve(M,G);
  //cout << "Done" << endl;
  cout << "number of parameters " << n << " number of dependent variables "
       << ndep << endl;
  cout << " elapsed time = " << adt.get_elapsed_time_and_reset()/1000 
       << " seconds " << endl;
  cout << "calculating variances" << endl;
  dvector var(1,ndep);
  for (int i=1;i<=ndep;i++)
  {
    var(i)=tdepgrad(i)*column(PP,i);
  }

  ifstream ifs2("deplabel.tmp");
  if (!ifs2)
  {
    cerr << "Error trying to open file deplabel.tmp" << endl;
  }
  line_adstring deplabel;
  dvector depvalue(1,ndep);

  adstring_array dstrings(1,ndep);
  for (int i=1;i<=ndep;i++)
  {
    ifs2 >> depvalue(i);
    ifs2 >> deplabel;
    dstrings(i)=deplabel;
    ofs << setw(3) << i; 
    if (var(i)>=0)
    {
      ofs << "  "  << setw(8) << setprecision(3) << sqrt(var(i));
    }
    else
    {
      cout << "variance = " << var(i) << endl;
      ofs << "  ********";
    }
    ofs << "  "  << setw(9) << setprecision(4) 
         << depvalue(i) << "  "  << deplabel << endl;
  }

  if (debugflag)
    cout << var << endl;
  cout << " elapsed time = " << adt.get_elapsed_time_and_reset()/1000 
       << " seconds " << endl;
   ofstream ofs2("vars");
   ofs2 << var << endl;
}

void check_hessian_pd(const char * _s)
{
  adstring s=_s;
  int debugflag=0; 
  //if (argc>1)
  //  debugflag=1;
  goto_set_num_threads(4);
  int n=0;

  adstring fname;
  adstring fname1;
  fname=s+".hes";
  fname1=s+".svd_report";
 
  uistream uis(fname);
  uis >> n;
  dvector a(1,n*n);
  uis >> a;
  if(!uis)
  {
    cerr << "error reading hessian from file " << fname << endl;
    ad_exit(1);
  }
  
 
  dmatrix M(1,n,1,n);
  int offset=0;
  for (int i=1;i<=n;i++)
  {
    M(i)=a(1+offset,n+offset).shift(1);
    offset+=n;
  }
  for (int i=1;i<=n;i++)
  {
    for (int j=1;j<=n;j++)
    {
      M(i,j)=0.5*(M(i,j)+M(j,i));
      M(j,i)=M(i,j);
    }
  }
  adtimer adt;
  svd SVD=svd_decomp(M);
  ofstream ofs(fname1);
  ofs << "Eigenvalues " << endl << SVD.get_S() << endl;
  ofs << endl << "Eigenvectors " << endl << SVD.get_U() << endl;
  cout << "Done" << endl;
  ad_exit(1);
}
int Lapack_Choleski_Inverse( int n,     // number of columns in matrix
                              int lda,   // leading dimension of matrix
                              double *a) // pointer to top-left corner
{
  int info=0;
#if !defined(USE_NO_LAPACKE)
  info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L',  n, a, n);
#else
   cerr <<  "Need to implement this for long double " << endl;
   ad_exit(1);
#endif
  if (info) // handle error conditions here
  {
    cout << "cholesky factorization failed" << endl;
    return info;
    //ad_exit(1);
  }
#if !defined(USE_NO_LAPACKE)
  info = LAPACKE_dpotri( LAPACK_ROW_MAJOR, 'L',  n, a, n);
  if (info) // handle error conditions here
  {
    cout << "cholesky inverse failed" << endl;
    return info;
    //ad_exit(1);
  }
#else
   cerr <<  "Need to implement this for long double " << endl;
   ad_exit(1);
#endif
  return info;
}

void matrix_multiplication(double *A, int A_width, int A_height,
         double *B, int B_width, int B_height,
         double *C);
// routines using lapack/blas called from mfcl have lapack_ prefix
dmatrix lapack_matrix_multiplcation(const dmatrix& A,const dmatrix& B)
{
  int acmin=A.indexmin();
  int acmax=A.indexmax();
  int armin=A(acmin).indexmin();
  int armax=A(acmin).indexmax();
  int bcmin=B.indexmin();
  int bcmax=B.indexmax();
  int brmin=B(acmin).indexmin();
  int brmax=B(acmin).indexmax();
  dvector va;
  dmatrix MA=make_dmatrix(acmin,acmax,armin,armax,va);
  dvector vb;
  dmatrix MB=make_dmatrix(bcmin,bcmax,brmin,brmax,vb);
  dvector vc;
  dmatrix MC=make_dmatrix(acmin,acmax,brmin,brmax,vc);

  MA=A;
  MB=B;
  MC.initialize();

  int A_width=armax-armin+1;
  int A_height=acmax-acmin+1;
  int B_width=brmax-brmin+1;
  int B_height=bcmax-bcmin+1;
  
  matrix_multiplication(&(va[1]),A_width,A_height,
         &(vb[1]),B_width,B_height,&(vc[1]));
  return MC;
}
#include <cblas.h>
void matrix_multiplication(double *A, int A_width, int A_height,
         double *B, int B_width, int B_height,
         double *C)
{
    int m = A_height;
    int n = B_width;
    int k = A_width;
    if (A_width != B_height)
    {
      cerr << "Matrix shape error in  matrix_multiplcation"
           << endl;
    }
    // This corresponds to 
    // C = 1.0 * A * B + 0.0 * C
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                m, n, k, 1.0, 
                A, k,
                B, n,
                0.0, 
                C, n);
}                

