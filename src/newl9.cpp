/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#define HOME_VERSION
#include "all.hpp"
#include "svd.h"
#include "variable.hpp"
#include <complex>
#if defined(NO_MY_DOUBLE_TYPE)
  #define lapack_complex_float std::complex<float>
  #define lapack_complex_double std::complex<double>
#endif
#include <lapacke.h>
#include <cblas.h>

#if defined(_WIN32)
#  if defined(close)
#    undef close
#  endif
#endif
const int  local_num_times=12;
MY_DOUBLE_TYPE local_min_bound[local_num_times]={1.e-13,1.e-12,1.e-11,1.e-10,1.e-9,1.e-8,1.e-7,1.e-5,1.e-3,
  1.e-2,1.e-1,1.0};

MY_REAL_DOUBLE * get_double_pointer(const dvector& v);
MY_REAL_DOUBLE * get_double_pointer(const dmatrix& M);
void double_pointer_values_return(dmatrix& m,MY_REAL_DOUBLE * pv);
void double_pointer_values_return(dvector& v,MY_REAL_DOUBLE * pv);

void  cblas_matrix_prod(int n,dvector& u,dvector & v,dvector & prod)
{
  MY_REAL_DOUBLE * pu=get_double_pointer(u);
  MY_REAL_DOUBLE * pv=get_double_pointer(v);
  MY_REAL_DOUBLE * pprod=get_double_pointer(prod);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    n, n, n, 1.0, pu, n,pv , n, 0.0,pprod,n);
  double_pointer_values_return(u,pu);
  double_pointer_values_return(v,pv);
  double_pointer_values_return(prod,pprod);
}

MY_DOUBLE_TYPE min(double x,MY_DOUBLE_TYPE y)
{
  MY_DOUBLE_TYPE z=static_cast<MY_DOUBLE_TYPE>(x);
  if (z<y)
    return z;
  else
   return y;
} 

MY_DOUBLE_TYPE max(double x,MY_DOUBLE_TYPE y)
{
  MY_DOUBLE_TYPE z=static_cast<MY_DOUBLE_TYPE>(x);
  if (z>y)
    return z;
  else
   return y;
} 

dmatrix get_pHess_sqrt(int nvar,dvar_len_fish_stock_history& fsh,dmatrix& eigvecs,dvector& veigvecs,
  dvector & eigvals)
{
  dvector vtmp;
  dmatrix tmp=make_dmatrix(1,nvar,1,nvar,vtmp);
  //dvector dinv=1.0/eigvals;
  double minval=0.001;
  if (fsh.parest_flags(386)>0)
  {
    minval=fsh.parest_flags(386)/1.e+9;
  }
  for (int i=1;i<=nvar;i++)
  {
    MY_DOUBLE_TYPE d=sqrt(max(minval,eigvals(i)));
    for (int j=1;j<=nvar;j++)
    {
      tmp(i,j)=eigvecs(j,i)*d;
    }
  }
  dvector v;
  dmatrix m=make_dmatrix(1,nvar,1,nvar,v);
  cblas_matrix_prod(nvar,veigvecs,vtmp,v);
  return m;
}

dmatrix get_pHess(int nvar,dvar_len_fish_stock_history& fsh,dmatrix& eigvecs,dvector& veigvecs,
  dvector & eigvals)
{
  dvector vtmp;
  dmatrix tmp=make_dmatrix(1,nvar,1,nvar,vtmp);
  //dvector dinv=1.0/eigvals;
  double minval=0.001;
  if (fsh.parest_flags(386)>0)
  {
    minval=fsh.parest_flags(386)/1.e+9;
  }
  for (int i=1;i<=nvar;i++)
  {
    MY_DOUBLE_TYPE d=max(minval,eigvals(i));
    for (int j=1;j<=nvar;j++)
    {
      tmp(i,j)=eigvecs(j,i)*d;
    }
  }
  dvector v;
  dmatrix m=make_dmatrix(1,nvar,1,nvar,v);
  cblas_matrix_prod(nvar,veigvecs,vtmp,v);
  return m;
}

dmatrix get_pHessinv_sqrt(int nvar,dvar_len_fish_stock_history& fsh,dmatrix& eigvecs,dvector& veigvecs,
  dvector & eigvals)
{
  dvector vtmp;
  dmatrix tmp=make_dmatrix(1,nvar,1,nvar,vtmp);
  //dvector dinv=1.0/eigvals;
  double minval=0.001;
  if (fsh.parest_flags(386)>0)
  {
    minval=fsh.parest_flags(386)/1.e+9;
  }
  for (int i=1;i<=nvar;i++)
  {
    MY_DOUBLE_TYPE d=sqrt(max(minval,eigvals(i)));
    for (int j=1;j<=nvar;j++)
    {
      tmp(i,j)=eigvecs(j,i)/d;
    }
  }
  dvector v;
  dmatrix m=make_dmatrix(1,nvar,1,nvar,v);
  cblas_matrix_prod(nvar,veigvecs,vtmp,v);
  return m;
}

dmatrix get_pHessinv(int nvar,dvar_len_fish_stock_history& fsh,dmatrix& eigvecs,dvector& veigvecs,
  dvector & eigvals)
{
  dvector vtmp;
  dmatrix tmp=make_dmatrix(1,nvar,1,nvar,vtmp);
  //dvector dinv=1.0/eigvals;
  double minval=0.001;
  if (fsh.parest_flags(386)>0)
  {
    minval=fsh.parest_flags(386)/1.e+9;
  }
  for (int i=1;i<=nvar;i++)
  {
    MY_DOUBLE_TYPE d=max(minval,eigvals(i));
    for (int j=1;j<=nvar;j++)
    {
      tmp(i,j)=eigvecs(j,i)/d;
    }
  }
  dvector v;
  dmatrix m=make_dmatrix(1,nvar,1,nvar,v);
  cblas_matrix_prod(nvar,veigvecs,vtmp,v);
  return m;
}
dmatrix get_Hessinv(int nvar,dvar_len_fish_stock_history& fsh,dmatrix& eigvecs,dvector& veigvecs,
  dvector & eigvals)
{
  dvector vtmp;
  dmatrix tmp=make_dmatrix(1,nvar,1,nvar,vtmp);
  //dvector dinv=1.0/eigvals;
  for (int i=1;i<=nvar;i++)
  {
    MY_DOUBLE_TYPE d=eigvals(i);
    for (int j=1;j<=nvar;j++)
    {
      if (fabs(d)<1.e-20)
      {
        tmp(i,j)=1.e+20;
      }
      else
      {
        tmp(i,j)=eigvecs(j,i)/d;
      }
    }
  }
  dvector v;
  dmatrix m=make_dmatrix(1,nvar,1,nvar,v);
  cblas_matrix_prod(nvar,veigvecs,vtmp,v);
  return m;
}


void make_correlation_report(dvar_len_fish_stock_history& fsh)
{
  int nvar;
  int nvar1;
  int num_deps;
  int i,j;

  adstring current_path;
  adstring root;

  adstring fname= ad_root+".hes";

  uistream ifs(fname);
  if (!ifs)
  {
    cerr << "Error trying to open file " << fname << endl;
  }
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
  // check for zero row
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
      cout << " Error there was a zero Hessian row - exiting " << endl;
      ad_exit(1);
    }
  }
  dvector veigvecs;
  dmatrix eigvecs=make_dmatrix(1,nvar,1,nvar,veigvecs);
  dvector eigvals;
  int ierr=0;
  dvector gscale(1,nvar);
  gscale.initialize();
  dvector pos_gscale(1,nvar);
  pos_gscale.initialize();
  dvector pos_hess_diag(1,nvar);
  pos_hess_diag.initialize();
  dvector hess_diag(1,nvar);
  hess_diag.initialize();
  {
    cout<<"starting to calculate hessian eigenvectors and eigenvalues"<<endl;
    {
      lapack_symmetric_eigen(nvar,vH,veigvecs,eigvals,ierr);
      //dmatrix pH=get_pHess(nvar,fsh,eigvecs,veigvecs,eigvals);
      dmatrix pHs=get_pHess_sqrt(nvar,fsh,eigvecs,veigvecs,eigvals);
      dmatrix pHsinv=get_pHessinv_sqrt(nvar,fsh,eigvecs,veigvecs,eigvals);
      dmatrix pHinv=get_pHessinv(nvar,fsh,eigvecs,veigvecs,eigvals);
      dmatrix Hinv=get_Hessinv(nvar,fsh,eigvecs,veigvecs,eigvals);
      int mmin1=pHsinv.indexmin();
      int mmax1=pHsinv.indexmax();
      dvector diag1(mmin1,mmax1);
      dvector diag2(mmin1,mmax1);
      dvector diag3(mmin1,mmax1);
      dmatrix corr(mmin1,mmax1,mmin1,mmax1);
      corr.initialize();
      //cout << pH*pHinv << endl;
      //ad_exit(1);
      {
        for (int i=mmin1;i<=mmax1;i++)
        {
          diag1(i)=pHsinv(i,i);
          diag2(i)=pHinv(i,i);
          diag3(i)=Hinv(i,i);
        }
        adstring pfname1= ad_root+"_pos_sqrt_hess_diag2";
        adstring pfname2= ad_root+"_pos_hess_inv_diag2";
        adstring pfname3= ad_root+"_pos_hess_cov";
        adstring pfname4= ad_root+"_pos_hess_cor";
        adstring pfname5= ad_root+"_hess_inv_diag";
        ofstream ofs1(pfname1);
        ofstream ofs2(pfname2);
        ofstream ofs3(pfname3);
        ofstream ofs4(pfname4);
        ofstream ofs5(pfname5);
        dmatrix output1(1,2,mmin1,mmax1);
        dmatrix output2(1,2,mmin1,mmax1);
        dmatrix output3(1,2,mmin1,mmax1);
        output1(2)=diag1;
        output2(2)=diag2;
        output3(2)=diag3;
        output1(1).fill_seqadd(1.0,1.0);
        output2(1).fill_seqadd(1.0,1.0);
        output3(1).fill_seqadd(1.0,1.0);
        ofs1 << trans(output1) << endl;
        ofs2 << trans(output2) << endl;
        ofs5 << trans(output3) << endl;
        ofs3 << pHinv << endl;
        for (int i=mmin1;i<=mmax1;i++)
        {
          for (int j=mmin1;j<=mmax1;j++)
          {
            corr(i,j)=pHinv(i,j)/sqrt((1.e-20+pHinv(i,i))*(1.e-20+pHinv(j,j)));
          }
        }
        ofs4 << "    " << setw(6) << output1(1) << endl;
        for (int i=mmin1;i<=mmax1;i++)
        {
          ofs4 << setw(6) << i << setw(6) << setprecision(3) << setfixed() 
               << corr(i) << endl;
        }
      }
      //lapack_symmetric_eigen(Hess,eigvecs,eigvals,ierr);
      { 
        dmatrix teigvecs=trans(eigvecs);   // now the rows of teigvecs are the eigenvectors

        dvector ek(1,nvar);
        double minval=0.001;
        if (fsh.parest_flags(386)>0)
        {
          minval=fsh.parest_flags(386)/1.e+9;
        }
        for (int k=1;k<=nvar;k++)
        {
          ek.initialize();
          ek(k)=1;
          MY_DOUBLE_TYPE tmp=0.0;
          MY_DOUBLE_TYPE tmp1=0.0;
          MY_DOUBLE_TYPE ptmp=0.0;
          MY_DOUBLE_TYPE ptmp1=0.0;
          for (int i=1;i<=nvar;i++)
          {
            tmp1+=square(teigvecs(i,k))*eigvals(i);
            tmp+=square(teigvecs(i,k))/eigvals(i);
            ptmp+=square(teigvecs(i,k))/max(minval,eigvals(i));
            ptmp1+=square(teigvecs(i,k))*max(minval,eigvals(i));
          }
          
          gscale(k)=tmp;         
          pos_gscale(k)=ptmp;         
          pos_hess_diag(k)=ptmp1;
          hess_diag(k)=tmp1;
        }
      }
      {
        adstring fname= ad_root+"_new.std";
        adstring pfname= ad_root+"_pos_new.std";
        adstring pfname3= ad_root+"_pos_sqrt_hess";
        adstring pfname4= ad_root+"_pos_sqrt_hess_inv";
        adstring pfname1= ad_root+"_pos_hess_diag";
        adstring pfname2= ad_root+"_hess_diag";
        {
          uostream p1(pfname3);
          p1 << nvar;
          p1 << pHs;
          p1.close();
          uostream p2(pfname4);
          p2 << nvar;
          p2 << pHsinv;
          p2.close();
        }
        ofstream ofs(fname);
        ofstream pofs(pfname);
        ofstream pofs1(pfname1);
        ofstream pofs2(pfname2);
        ofs << gscale << endl;
        pofs << pos_gscale << endl;
        pofs1 << nvar;
        pofs1 << pos_hess_diag << endl;
        pofs2 << hess_diag << endl;
      }
      if (ierr)
      {
        cerr << "Error calculating eigenvectors in "
             << " lapack_symmetric_eigen(Hess,eigvecs,eigvals,ierr)"
             << endl;
        ad_exit(1);
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
      for (int i=1;i<=nvar;i++)
      {
        for (int j=1;j<=nvar;j++)
        {
          tmp(i,j)=eigvecs(j,i)*dinv(i);
        }
      }
      cout<<"starting to invert hessian nvar = " << nvar <<endl;
      cblas_matrix_prod(nvar,veigvecs,vtmp,vhessinv);

      //uostream uos1("bin_eigenvectors");
     // uostream uos2("bin_eigenvalues");
     // ofstream ofs1("eigenvectors");
      //ofstream ofs2("eigenvalues");
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

      cout.flush();
      {
         ofstream ofs("new_cor_report");
         MY_DOUBLE_TYPE cut=0.0;
         ofs << " Smallest eigenvalues" << endl;
         int nrep=min(nvar/2,80);
         MY_DOUBLE_TYPE p10=0.10;
         for(int ind=1;ind<=nrep;ind++)
         {
           cut=max(0.8*max(fabs(sorted_eigen_vectors(ind))),p10);
           cut=min(0.4,cut);
           ofs << setprecision(5) << setw(12) 
               << sorted_eigen_vectors(ind,0) << "   ";
           int ip=0;
           cut*=1.4;
           do
           {
             cut/=1.4;   //more detail: cut/=1.8
             for (int i=1;i<=nvar;i++)
             {
               if (fabs(sorted_eigen_vectors(ind,i))>cut) //more detail cut*0.1
               {
                 ip++;
               }
             }
           }
           while (ip==0 || (ip ==1 && cut> 0.08) || (ip==2 && cut>0.12));
           //for (int i=1;i<=nvar;i++)
           dmatrix tmpmat(1,nvar,1,2);
           tmpmat.initialize();
           int ii=0;
           for (int i=1;i<=nvar;i++)
           {
             if (fabs(unsorted_eigen_vectors(ind,i))>0.01)
             {
               tmpmat(++ii,1)=i;
               tmpmat(ii,2)=-fabs(unsorted_eigen_vectors(ind,i));
             }
           }
           dmatrix tmpmat2=sort(tmpmat.sub(1,ii),2);
           ivector index=ivector(column(tmpmat2,1));
             
           
           //for (int i=1;i<=nvar;i++)
           for (int i=1;i<=ii;i++)
           {
             //if (fabs(sorted_eigen_vectors(ind,i))>0.01)
             {
               //ofs << "(" << i << " " << setprecision(3) 
               //    << sorted_eigen_vectors(ind,i) << ") ";
               ofs << "(" << index(i) << " " << setprecision(3) 
                   << sorted_eigen_vectors(ind,index(i)) << ") ";
             }
           }
           ofs << endl;
         }
         ofs << endl << " Largest eigenvalues" << endl;
         for(int ind=nvar;ind>nvar-nrep;ind--)
         {
           cut=max(0.8*max(fabs(sorted_eigen_vectors(ind))),p10);
           cut=min(0.4,cut);
           cout  << setprecision(5) << setw(12) 
               << sorted_eigen_vectors(ind,0) << "   ";
           cout.flush();
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
      //do

      for (int i=1;i<=nvar;i++)
      {
        if (eigvals(i)<=0.)
        {
          ofs3 << setprecision(15) << setscientific() << eigvals(i) 
               << "  "  << eigvecs(i) << endl;
        }
      }
    
      //uos1 << nvar;
      //uos1 << eigvals;
      //uos2 << nvar;
      //uos2 << eigvecs;
      //ofs1 << nvar << endl;
      //ofs1 << setprecision(15) << setscientific() << eigvals << endl;
      //ofs2 << nvar << endl;
      //ofs2 << setprecision(15) << setscientific() << eigvecs << endl;
      //ad_exit(1);
    }
  }

// Code for including no-effort dep_vars in correlation report NMD_6Apr2022
  adstring fname2;
  fname2=ad_root+".dp2";
  uistream uis2(fname2);
  if (!uis2)
  {
    cerr << "Error trying to open file " << fname2 << endl;
  }
  int zero_fflag=0;
  if (!uis2)
  {
    cout << " Correlation report excludes zero-fishing dep_vars" << endl;
  }
  else
  {
    cout << " Correlation report includes zero-fishing dep_vars" << endl;
    zero_fflag=1;
  }

  fname=ad_root+".dep";
  uistream uis1(fname);
  if (!uis1)
  {
    cerr << "Error trying to open file " << fname << endl;
  }
  uis1 >> nvar1 >> num_deps;
  if (!uis1)
  {
    cerr << "Error reading nvar from depgrad.rpt" << endl;
  }
  if (nvar!=nvar1)
  {
    cerr << "number of independent variables is not consistent"
            " between files" << "Hess.rpt and depgrad.rpt" << endl;
    ad_exit(1);
  }
  dmatrix depgrad;
  dmatrix covar;
  dvector stddev;
  line_adstring deplabel;
  line_adstring ndeplabel;
  dvector depvalue;
  dvector ndepvalue;

  
  if (zero_fflag==1)
  {
    int n2=0;
    int ndep2=0;
    uis2 >> n2 >> ndep2;
    if (!uis2)
    {
      cerr << "Error trying to open file " << fname2 << endl;
      cerr << "Error reading n2 from file " << fname2 << endl;
      ad_exit(1);
    }
    if (nvar1!=n2) {
      cerr << "number of independent variables is not consistent"
              " between depgrad files" << endl;
      ad_exit(1);
    }
    dmatrix tdepgrad(1,num_deps+ndep2,1,nvar1);

    covar.allocate(1,num_deps+ndep2,1,num_deps+ndep2);
    covar.initialize();
    stddev.allocate(1,num_deps+ndep2);
    fname=ad_root+".cor";
    ofstream ofs(fname);
    ifstream ifs2("deplabel.tmp");
    ifstream ifs3("deplabel_noeff.tmp");
    depvalue.allocate(1,num_deps);
    depvalue.initialize();
    adstring_array dstrings(1,num_deps);
    ndepvalue.allocate(1,ndep2);
    ndepvalue.initialize();
    adstring_array ndstrings(1,ndep2);
    
    {
      //read in as doubles and convert
      MY_REAL_DOUBLE tmp1;
      //char tmp[8];
      for (int i=1;i<=num_deps;i++)
      {
        for (int j=1;j<=nvar1;j++)
        {
          // horrrible cast to char * to read in
          uis1.read((char*)(&tmp1),8);
          tdepgrad(i,j)=tmp1;
          if (!uis1)
          {
            cerr << "Error reading gradients from " << fname << endl;
            ad_exit(1);
          }
        }
        ifs2 >> depvalue(i);
        ifs2 >> deplabel;
        dstrings(i)=deplabel;
      }
    }
    {
      //read in as doubles and convert
      MY_REAL_DOUBLE tmp1;
      //char tmp[8];
      for (int i=1;i<=ndep2;i++)
      {
        for (int j=1;j<=nvar1;j++)
        {
          // horrrible cast to char * to read in
          uis2.read((char*)(&tmp1),8);
          tdepgrad(i+num_deps,j)=tmp1;
          if (!uis2)
          {
            cerr << "Error reading gradients from " << fname2 << endl;
            ad_exit(1);
          }
        }
        ifs3 >> ndepvalue(i);
        ifs3 >> ndeplabel;
        ndstrings(i)=ndeplabel;
      }
    }
    for (int i=1;i<=ndep2;i++)
    {
//      tdepgrad(num_deps+i)-=tdepgrad(i);
      tdepgrad(num_deps+i)=tdepgrad(i)-tdepgrad(num_deps+i);  //NMD_31may2022
    }

    for (i=1;i<=num_deps+ndep2;i++)
    {
      dvector tmp=hessinv*tdepgrad(i);
      for (j=1;j<=i;j++)
      {
        covar(i,j)=tdepgrad(j)*tmp;
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
    
    for (i=1;i<=num_deps+ndep2;i++)
    {
      ofs << "                  "  << setw(9) << i << "  " ;
    }
    ofs << endl;


    for (i=1;i<=num_deps+ndep2;i++)
    {
      cout<<"calc correlation of dep vars. Got to "<<i<<endl;
      if (i<=num_deps)
      {
        ofs << setw(3) << i; 
        ofs << "  "  << setw(9) << setprecision(4) 
            << depvalue(i) << "  "  << dstrings(i) << endl;
      }
      else
      {
        ofs << setw(3) << i; 
        ofs << "  "  << setw(9) << setprecision(4) 
            << depvalue(i-num_deps) - ndepvalue(i-num_deps) << "  "  
            << dstrings(i-num_deps) << " - " << ndstrings(i-num_deps) << endl;
      }
      for (j=1;j<=i;j++)
      {
        MY_DOUBLE_TYPE tmp=stddev(i)*stddev(j);
        if (tmp)
          ofs << "  "  << setw(8) << setprecision(3) << covar(i,j)/tmp;
        else
          ofs << "  "  << setw(8) << setprecision(3) << tmp;
      }
      ofs << endl;
    }

    cout << " Here" << endl;
  }
  else
  {
    // Existing code
    depgrad.allocate(1,num_deps,1,nvar);
    covar.allocate(1,num_deps,1,num_deps);
    covar.initialize();
    stddev.allocate(1,num_deps);
    fname=ad_root+".cor";
    ofstream ofs(fname);
    ifstream ifs2("deplabel.tmp");
    depvalue.allocate(1,num_deps);
    depvalue.initialize();
    adstring_array dstrings(1,num_deps);

    for (i=1;i<=num_deps;i++)
    {
      uis1 >> depgrad(i);
      ifs2 >> depvalue(i);
      ifs2 >> deplabel;
      dstrings(i)=deplabel;
      
      if (!uis1)
      {
        cerr << "Error reading gradient for dependent variable " << i
             << "  from depgrad.rpt" << endl;
        ad_exit(1);
      }
    }
    for (i=1;i<=num_deps;i++)
    {
      dvector tmp=hessinv*depgrad(i);
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
    for (i=1;i<=num_deps;i++) {
      ofs << "                  "  << setw(9) << i << "  " ;
    }
    ofs << endl;

    for (i=1;i<=num_deps;i++)
    {
      cout<<"calc correlation of dep vars. Got to "<<i<<endl;
      ofs << setw(3) << i; 
      ofs << "  "  << setw(9) << setprecision(4) 
           << depvalue(i) << "  "  << dstrings(i) << endl;
      for (j=1;j<=i;j++)
      {
        MY_DOUBLE_TYPE tmp=stddev(i)*stddev(j);
        if (tmp)
          ofs << "  "  << setw(8) << setprecision(3) << covar(i,j)/tmp;
        else
          ofs << "  "  << setw(8) << setprecision(3) << tmp;
      }
      ofs << endl;
    }
  }

  int negflag=0;
  int num_times=1;
  if (min(eigvals)<1.e-10)
  {
    negflag=1;
    num_times=local_num_times;
  }
  
//////////////////////  Positivized std-dev calculations  /////////////////
  if (fsh.parest_flags(384)==1)
  {   //NMD_9sep2020
  adstring stdfname=ad_root+"_positivized" + ".std";
  {
    ofstream ofsstd(stdfname);
    for (int ii=1;ii<=num_times;ii++)
    {
      dvector vhessinv1;
      dmatrix hessinv1=make_dmatrix(1,nvar,1,nvar,vhessinv1);
      if (negflag)
      {
        fname=ad_root+"_" + str(ii)+".cor";
      }
      else
      {
        fname=ad_root+".cor";
      }
      ofstream ofs(fname);
      dvector ppos(eigvals.indexmin(),eigvals.indexmax());
      ppos=eigvals;
      if (negflag)
      {
        for (int i=eigvals.indexmin();i<=eigvals.indexmax();i++)
        {
          ppos(i)=max(ppos(i),local_min_bound[ii-1]);
        }
      }
      dvector dinv=1.0/ppos;
      //dvector dinv=1.0/eigvals;
      dvector vtmp;
      dmatrix tmp=make_dmatrix(1,nvar,1,nvar,vtmp);
      for (int i=1;i<=nvar;i++)
      {
        for (int j=1;j<=nvar;j++)
        {
          tmp(i,j)=eigvecs(j,i)*dinv(i);
        }
      }
      cout<<"starting to invert hessian nvar = " << nvar <<endl;
      cblas_matrix_prod(nvar,veigvecs,vtmp,vhessinv1);
      cout<<"finished inverting hessian nvar = " << nvar <<endl;
      //depgrad(1)=0.5*(eigvecs(1)+depgrad1);

      int jmin=1;
      for (int i=1;i<=num_deps;i++)
      {
        dvector tmp=hessinv1*depgrad(i);
        for (j=jmin;j<=i;j++)
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
      cout << "FFF" << ii << " " << eigvecs(1)*(hessinv1*eigvecs(1)) << endl;
      for (int i=1;i<=num_deps;i++) {
        ofs << "                  "  << setw(9) << i << "  " ;
      }
      ofs << endl;

      for (int i=1;i<=num_deps;i++)
      {
        ofsstd << setw(9) << stddev(i) << " ";
        cout<<"calc correlation of dep vars. Got to "<<i<<endl;
        ofs << setw(3) << i;
        ofs << "  "  << setw(9) << setprecision(4)
             << depgrad(i) << "  "
             << "  "  << setw(10) << setprecision(6)
             << depvalue(i) << "  "
             << "  "  << setw(10) << setprecision(6)
             << stddev(i)
             << "  "  << setw(9) << setprecision(4)
             << "  " << deplabel(i) << endl;
        for (j=1;j<=i;j++)
        {
          MY_DOUBLE_TYPE tmp=stddev(i)*stddev(j);
          if (tmp)
            ofs << "  "  << setw(8) << setprecision(3) << covar(i,j)/tmp;
          else
            ofs << "  "  << setw(8) << setprecision(3) << tmp;
        }
        ofs << endl;
      }
      ofsstd << endl;
    }
  }
  {
    ifstream ifsstd(stdfname);
    dmatrix stds(1,num_times,1,num_deps);
    dmatrix stds1(1,num_times,1,num_deps);
    ifsstd >> stds;
    if (!ifsstd)
    {
      cerr << "Error reading from file " << stdfname << endl;
      ad_exit(1);
    }
    for (int ii=2;ii<=num_times;ii++)
    {
      for (int j=1;j<=num_deps;j++)
      {
        if (stds(1,j)==0.0)
          stds1(ii,j)=0.0;
        else
          stds1(ii,j)=stds(ii,j)/stds(1,j);
      }
    }
    for (int j=1;j<=num_deps;j++)
    {
      if (stds(1,j)==0.0)
        stds1(1,j)=0.0;
      else
        stds1(1,j)=1.0;
    }
    //stds(1)=1.0;
    adstring cstdfname=ad_root+"_scaled_positivized" + ".std";
    adstring cstdfname1=ad_root+"scaled_positivized1" + ".std";
    adstring cstdfname2=ad_root+"positivized" + ".std";
    ofstream ofsstd(cstdfname);
    ofstream ofsstd1(cstdfname1);
    ofstream ofsstd2(cstdfname2);
    for (int j=1;j<=num_deps;j++)
    {
      ofsstd << " " << setw(4) << j << " " << setfixed()
             << setprecision(4) << setw(8) << column(stds1,j) << endl;
      ofsstd1 << " " << setw(4) << j << " " << setscientific()
             << setprecision(6) << setw(12) << column(stds1,j) << endl;
      ofsstd2 << " " << setw(4) << j << " " << setscientific()
             << setprecision(6) << setw(12) << column(stds,j) << endl;
    }
  }
  }     //NMD_9sep2020
}

void lapack_symmetric_eigen(int n,dvector & vM,dvector& a,
  dvector& eigenvalues,int & ierr)
{
  int lda=n;
  if (!allocated(a))
  {
    a.allocate(1,n*n);
  }
  a=vM;
  if (!allocated(eigenvalues))
  {
    eigenvalues.allocate(1,n);
  }

  /* Locals */
  int info=0;
  MY_REAL_DOUBLE * pa=get_double_pointer(a);
  MY_REAL_DOUBLE * pe=get_double_pointer(eigenvalues);
  info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n,pa,
    lda,pe);
  double_pointer_values_return(a,pa);
  double_pointer_values_return(eigenvalues,pe);
  ierr=info;
  /* Check for convergence */
  if( info != 0 )
  {
    cout <<  "The algorithm in LAPACKE_dsyev failed to compute eigenvalues."
         << " info = " << info << endl;
  }
}

dmatrix make_dmatrix(int cmin,int cmax,int rmin,int rmax,dvector & v)
{
  if (!allocated(v)) v.allocate(1,(cmax-cmin+1)*(rmax-rmin+1));
  dmatrix M(cmin,cmax);
  int offset=0;
  for (int i=cmin;i<=cmax;i++)
  {
    M(i)=v(offset+rmin,offset+rmax).shift(rmin);
    offset+=rmax-rmin+1;
  }
  return M;
}

void stitch_parallel_hessian(void)
{
  cifstream ifs1("parall_hess");
  if (!ifs1)
  {
    cerr << "error opening file parall_hess" << endl;
    ad_exit(1);
  }
  // num_parallel_processes
  int np=0;
  ifs1 >>  np;
  cout << "  num_parallel_processes = " << np << endl;
  int nvars;
  ifs1 >> nvars;
  cout << "  num_indep vars = " << nvars << endl;  
  ivector nvstart(1,np);
  ifs1 >> nvstart;
  cout << " start positions = " << nvstart << endl;

  // Calculate rows in each *.hess file
  ivector nrws(1,np);
  for (int i=1; i<=np-1; i++)
  {
    nrws(i)=nvstart(i+1)-nvstart(i);
  }
  nrws(np)=nvars-(nvstart(np)-1);

  cout << " No. rows in each file: " << nrws << endl;
  cout << " Total rows: " << sum(nrws) << "  nvars: " << nvars << endl;
  cout << " Start positions : "  << nvstart << endl;

  int nvarps;
  dvector vH;
  dmatrix Hess=make_dmatrix(1,nvars,1,nvars,vH);

  for (int i=1; i<=np; i++)
  {
    adstring iternum= itoa(i,10);
    adstring alname= ad_root + ".hes_" + iternum;
    cout << " Filename for :  " << i << " " << alname << endl;
    uistream ifs(alname);
    if (!ifs)
    {
      cerr << "Error trying to open file " << alname << endl;
      ad_exit(1);
    }
    int hmin;
    int hmax; 
    ifs >> nvarps >> hmin >> hmax;

    if (!ifs)
    {
      cerr << "Error reading nvar from " << alname << endl;
    }
    cout << " nvar for file: " << alname << " = " << nvarps
         << " hmin:  " <<  hmin << " hmax:  " << hmax << endl;

    dvector vH;
    dmatrix tmpHess=make_dmatrix(1,nrws(i),1,nvars,vH);
    ifs >> tmpHess;
    if (!ifs)
    {
      cerr << "Error reading hessian from " << alname << endl;
    }
    int jj=0;
    //    int mmin=nvstart(i);
    //    int mmax=nvend(i);
    for (int j=hmin; j<=hmax; j++)
    {
      jj++;
      Hess(j)=tmpHess(jj);
    }

  }

  // Ouput stitched Hessian file
  adstring hessfile=ad_root + adstring(".hes");
  cout << " Writing the Hessian file:  " << hessfile << endl;
  uostream ofs(hessfile);
  if (!ofs)
  {
    cerr << "Error opening file " << (char*) hessfile << endl;
    exit(1);
  }
  ofs <<  "  " << nvars;
  ofs << Hess;

  ofs.CLOSE();

  cout << " Completed stitch of Hessian file: " << hessfile << endl;
  {
    uistream ifs(hessfile);
    int nvar;
    ifs >> nvar;
    dmatrix h(1,nvar,1,nvar);
    for (int i=1;i<=nvar;i++)
    {
      ifs >> h(i); 
    }
    MY_DOUBLE_TYPE fmax=0.0;
    for (int i=1;i<=nvar;i++)
    {
      for (int j=1;j<=nvar;j++)
      {
        if (fabs(h(i,j)-h(j,i))> fmax)
          fmax=fabs(h(i,j)-h(j,i));
      }
    }
    cout << "fmax = " << setprecision(15) << fmax << endl;
  }

}

#undef HOME_VERSION
