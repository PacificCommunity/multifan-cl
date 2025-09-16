/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
#ifdef __MSC__
  MY_DOUBLE_TYPE mfexp(const MY_DOUBLE_TYPE& );
  dvariable mfexp(const prevariable& );
  dvar_vector mfexp(const dvar_vector& );
  dvector mfexp(const dvector& );
#else
  MY_DOUBLE_TYPE mfexp(MY_DOUBLE_TYPE& );
  dvariable mfexp(prevariable& );
  //dvar_vector mfexp(dvar_vector& );
  //dvector mfexp(dvector& );
#endif

void dvar_fish_stock_history::get_orthogonal_recruitment_poly_info(void)
{
  int local_print_switch=0;
  if (!pmsd || pmsd->current_species==1)
  {
    ofstream ofs("matrices");

    // These are not degrees -- they are the number of parameters
    // which is the degree -1 
    degree_reg=0;
    degree_yr=0;
    degree_ses=0;
    degree_ses_reg=0;
    // use degree  since that is the number of parameters
    if (parest_flags(155)>0)
    {
      degree_yr=parest_flags(155);
    }
    if (parest_flags(217)>0)
    {
      degree_ses=parest_flags(217);
    }
    if (parest_flags(216)>0)
    {
      degree_reg=parest_flags(216);
    }
    if (parest_flags(218)>0)
    {
      degree_ses_reg=parest_flags(218);
    }
   /*
   // pf(155) is the default for degree of orthogonal polys for all these
   if (parest_flags(221)!=0)
     degree_yr=parest_flags(221);
   degree_reg=parest_flags(155);
   if (parest_flags(216)>0)
   {
     degree_reg=parest_flags(216)-1;
   }
   else if(parest_flags(216)<0)    //NMD19dec2017
   {
     degree_reg=0;
   }
   degree_ses=parest_flags(155);
   if (parest_flags(217)>0)
   {
     degree_ses=parest_flags(217);
   }
   else if(parest_flags(217)<0)
   {
     degree_ses=0;
   }
   degree_ses_reg=parest_flags(155);
   if (parest_flags(218)>0)
   {
     degree_ses_reg=parest_flags(218);
   }
   else if(parest_flags(218)<0)
   {
     degree_ses_reg=0;
   }    //NMD19dec2017
   */
   
    int i,j; 
    int N=num_seasons;
    int M=num_regions;
    //N=4;
    //M=3;
    if (pmsd) M=pmsd->get_num_real_regions();
    dmatrix a(1,N,1,M);
    d3_array b(1,N-1,1,N,1,M);
    d3_array c(1,M-1,1,N,1,M);
    d4_array d(1,N-1,1,M-1,1,N,1,M);
  
    a.initialize();
    b.initialize();
    c.initialize();
    d.initialize();
  
    imatrix & Active = ses_reg_recr_flags;
    a=1.0;
    ivector csum=colsum(Active);
    ivector rsum=rowsum(Active);
    for (i=1;i<=N-1;i++)
    {
      b(i,i)=1.0;
    }
  
    for (j=1;j<=M-1;j++)
    {
      for (i=1;i<=N;i++)
      {
        c(j,i,j)=1.0;
      }
    }
    int Nrownzero=N-1;
    int Ncolnzero=M-1;
    int Nrowcolnzero=(N-1)*(M-1);
  
  
    for (i=1;i<=N-1;i++)
    {
      for (j=1;j<=M-1;j++)
      {
        d(i,j,i,j)=1.0;
      }
    }
    
    d3_array v(1,N*M,1,N,1,M);
    d3_array w(1,N*M,1,N,1,M);
    int ii=0;
    v(++ii)=a;
    int N1=N-1;
    int M1=M-1;
    int NM1=(N-1)*(M-1);
    for (i=1;i<=N-1;i++)
    {
      v(++ii)=b(i);
    }
     
    for (j=1;j<=M-1;j++)
    {
      v(++ii)=c(j);
    }
    for (i=1;i<=N-1;i++)
    {
      for (j=1;j<=M-1;j++)
      {
        v(++ii)=d(i,j);
      }
    } 
    int nbasis=ii;
  
    // zero out the no recruits
    for (i=1;i<=N;i++)
    {
      for (j=1;j<=M;j++)
      {
        if (Active(i,j)==0)
        { 
          for (ii=1;ii<=nbasis;ii++)
          { 
            v(ii)(i,j)=0.0;
          }
        }
      }
    }
    dvector nv(1,nbasis);
    for (ii=1;ii<=nbasis;ii++)
    { 
      nv(ii)=norm(v(ii)); 
    }
    for (ii=1;ii<=nbasis;ii++)
    { 
      if (nv(ii)>1.e-8)
        v(ii)/=norm(v(ii));
      else
        v(ii)=0.0;
    }
      
    
    // gram-schmtdt 
    ii=0;
  
    for (i=1;i<=nbasis;i++)
    { 
      MY_DOUBLE_TYPE nm=norm(v(i));
      if (nm < 1.e-10)
      {
        if (i<=N)
        {
          int ioffset=i-1;
          N1--;
          Nrownzero--;
        }
        else if (i<=N+M-1)
        {
          int ioffset=i-1-(N-1);
          Ncolnzero--;
          M1--;
        }
        else 
        {
          int ioffset=i-1-(M-1)-(N-1);
          int iindex=(ioffset-1)/(M-1)+1;
          int jindex=(ioffset-1)%(M-1)+1;
          NM1--;
        }
      }
      else
      {
        w(++ii)=v(i)/=norm(v(i));
        for (j=i+1;j<=nbasis;j++)
        { 
          v(j)-=dot(v(j),w(ii))*w(ii);
        }
      }
    }
    if (allocated(numcomp))
    {
      numcomp.deallocate();
    }
    numcomp.allocate(1,4);
  
    numcomp(1)=1;
    numcomp(2)=N1;
    numcomp(3)=M1;
    numcomp(4)=NM1;
  
    if (allocated(orth_recr_basis))
    {
      orth_recr_basis.deallocate();
    }
    orth_recr_basis.allocate(1,ii,1,N,1,M);
    orth_recr_basis=w.sub(1,ii);
    if (local_print_switch)
    {
      for (int i=1;i<=ii;i++)
      {
         ofs << endl << " i = " << i << endl;
         ofs << orth_recr_basis(i) << endl;
      }
      ad_exit(1);
    }
  }
  else
  {
    int is=pmsd->current_species;
    ivector & numcomp=pmsd->numcomp(is);
    ivector pf=pmsd->parest_flags(is);

    // These are not degrees -- they are the number of parameters
    // which is the degree -1 
    degree_reg=0;
    degree_yr=0;
    degree_ses=0;
    degree_ses_reg=0;
    // use degree  since that is the number of parameters
    if (pf(155)>0)
    {
      pmsd->degree_yr(is)=pf(155);
    }
    if (pf(217)>0)
    {
      pmsd->degree_ses(is)=pf(217);
    }
    if (pf(216)>0)
    {
      pmsd->degree_reg(is)=pf(216);
    }
    if (pf(218)>0)
    {
      pmsd->degree_ses_reg(is)=pf(218);
    }
   /*
   // pf(155) is the default for degree of orthogonal polys for all these
   pmsd->degree_yr(is)=pf(155);
   if (pf(221)!=0)
     pmsd->degree_yr(is)=pf(221);
   pmsd->degree_reg(is)=pf(155);
   if (pf(216)>0)
   {
     pmsd->degree_reg(is)=pf(216);
   }
   else if(pf(216)<0)    //NMD19dec2017
   {
     pmsd->degree_reg(is)=0;
   }
   pmsd->degree_ses(is)=pf(155);
   if (pf(217)>0)
   {
     pmsd->degree_ses(is)=pf(217);
   }
   else if(pf(217)<0)
   {
     pmsd->degree_ses(is)=0;
   }
   pmsd->degree_ses_reg(is)=pf(155);
   if (pf(218)>0)
   {
     pmsd->degree_ses_reg(is)=pf(218);
   }
   else if(pf(218)<0)
   {
     pmsd->degree_ses_reg(is)=0;
   }    //NMD19dec2017
   */
    int i,j; 
    int N=num_seasons;
    int M=num_regions;
    if (pmsd) M=pmsd->get_num_real_regions();
    dmatrix a(1,N,1,M);
    d3_array b(1,N-1,1,N,1,M);
    d3_array c(1,M-1,1,N,1,M);
    d4_array d(1,N-1,1,M-1,1,N,1,M);
  
    a.initialize();
    b.initialize();
    c.initialize();
    d.initialize();
  
    imatrix & Active = pmsd->ses_reg_recr_flags(is);
    a=1.0;
    ivector csum=colsum(Active);
    ivector rsum=rowsum(Active);
    for (i=1;i<=N-1;i++)
    {
      b(i,i)=1.0;
      /*
      for (j=1;j<=M;j++)
      {
        b(i,i,j)=1.0;
        for (int ii=N;ii>=1;ii--)
        {
          b(i,N,j)-=1.0;
        }
      }
      */
    }
  
    for (j=1;j<=M-1;j++)
    {
      for (i=1;i<=N;i++)
      {
        c(j,i,j)=1.0;
        /*
        for (int jj=M;jj>=1;jj--)
        {
          c(j,i,M)-=1.0;
        }
        */
      }
    }
    int Nrownzero=N-1;
    int Ncolnzero=M-1;
    int Nrowcolnzero=(N-1)*(M-1);
  
  
    for (i=1;i<=N-1;i++)
    {
      for (j=1;j<=M-1;j++)
      {
        d(i,j,i,j)=1.0;
        /*
        d(i,j,i,M)=-1.0;
        d(i,j,N,j)=-1.0;
        d(i,j,N,M)=+1.0;
        */
      }
    }
    
    d3_array v(1,N*M,1,N,1,M);
    d3_array w(1,N*M,1,N,1,M);
    int ii=0;
    v(++ii)=a;
    int N1=N-1;
    int M1=M-1;
    int NM1=(N-1)*(M-1);
    for (i=1;i<=N-1;i++)
    {
      v(++ii)=b(i);
    }
     
    for (j=1;j<=M-1;j++)
    {
      v(++ii)=c(j);
    }
    for (i=1;i<=N-1;i++)
    {
      for (j=1;j<=M-1;j++)
      {
        v(++ii)=d(i,j);
      }
    } 
    int nbasis=ii;
  
    // zero out the no recruits
    for (i=1;i<=N;i++)
    {
      for (j=1;j<=M;j++)
      {
        if (Active(i,j)==0)
        { 
          for (ii=1;ii<=nbasis;ii++)
          { 
            v(ii)(i,j)=0.0;
          }
        }
      }
    }
    dvector nv(1,nbasis);
    for (ii=1;ii<=nbasis;ii++)
    { 
      nv(ii)=norm(v(ii)); 
    }
    for (ii=1;ii<=nbasis;ii++)
    { 
      if (nv(ii)>1.e-8)
        v(ii)/=norm(v(ii));
      else
        v(ii)=0.0;
    }
      
    
    // gram-schmtdt 
    ii=0;
  
    for (i=1;i<=nbasis;i++)
    { 
      MY_DOUBLE_TYPE nm=norm(v(i));
      if (nm < 1.e-10)
      {
        if (i<=N)
        {
          int ioffset=i-1;
          N1--;
          Nrownzero--;
        }
        else if (i<=N+M-1)
        {
          int ioffset=i-1-(N-1);
          Ncolnzero--;
          M1--;
        }
        else 
        {
          int ioffset=i-1-(M-1)-(N-1);
          int iindex=(ioffset-1)/(M-1)+1;
          int jindex=(ioffset-1)%(M-1)+1;
          NM1--;
        }
      }
      else
      {
        w(++ii)=v(i)/=norm(v(i));
        for (j=i+1;j<=nbasis;j++)
        { 
          v(j)-=dot(v(j),w(ii))*w(ii);
        }
      }
    }
    if (allocated(numcomp))
    {
      numcomp.deallocate();
    }
    numcomp.allocate(1,4);
  
    numcomp(1)=1;
    numcomp(2)=N1;
    numcomp(3)=M1;
    numcomp(4)=NM1;
  
    if (allocated(pmsd->orth_recr_basis(is)))
    {
      pmsd->orth_recr_basis(is).deallocate();
    }
    pmsd->orth_recr_basis(is).allocate(1,ii,1,N,1,M);
    pmsd->orth_recr_basis(is)=w.sub(1,ii);
  }
}

  
#undef HOME_VERSION


