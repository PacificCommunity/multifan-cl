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

#if !defined(ZFEV)
extern dvar_len_fish_stock_history * pcfsh;
#endif
double *  __gdouble__=0;

dvar_vector lineup(dvar_matrix& m)
{
  int mmin=m.indexmin();
  int mmax=m.indexmax();
  int cmin=m(mmin).indexmin();
  int n=0;
  ivector sizes(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    sizes(i)=m(i).indexmax()-m(i).indexmin()+1;
    n+=sizes(i);
  }
  dvar_vector v(cmin,n+cmin-1);
  int offset=cmin;
  for (int i=mmin;i<=mmax;i++)
  {
    v(offset,offset+sizes(i)-1).shift(m(i).indexmin())=m(i);
    offset+=sizes(i);
  }
  return v;
}
dvar_matrix stack(dvar_vector& v,int nr)  // break up vector into nc pieces
{                                         // and stack them up
  int mmin=v.indexmin();
  int mmax=v.indexmax();
  int sz=mmax-mmin+1;
  if (sz%nr)
  {
    cerr << "Incompatible dimension in stack" << endl;
    ad_exit(1);
  }
  int ncols=sz/nr;
  dvar_matrix M(1,nr,mmin,mmin+ncols-1);
  int offset=mmin;
  for (int i=1;i<=nr;i++)
  {
    M(i)=v(offset,offset+ncols-1).shift(M(i).indexmin()); // offset to deal  
                                                 // with multi-species kludge
    //M(i)=v(offset,offset+ncols-1).shift(mmin); // this is the "correct" way
    offset+=ncols;
  }
  return M;
}


void dvar_fish_stock_history::setup_diffusion(void)
{
  if (age_flags(184)==2)
  {
    movement_coffs_option_2_to_0();
    //cout << "PPP2  " << norm(value(diff_coffs)) << endl;
  }
  if (age_flags(184)==1)
  {
    // this as for the orthogonal parameterization
    //get_diff_parameterization();
    // this is the new xdiff parameterization
    get_old_diff_parameterization();
  }
  setup_diffusion2();
}

void dvar_fish_stock_history::setup_diffusion2(void)
{
  int mmin=xdiff_coffs.indexmin();
  int mmax=xdiff_coffs.indexmax();

  int ns=1;
  int nrr=num_regions;
  if (pmsd)
  {
    ns=pmsd->num_species;
    nrr=pmsd->num_real_regions;
  }
  Dad.initialize();
  if (!age_flags(89))
  {
    int offset=0;
    int rmin=1;
    int rmax=num_regions;
    //ofstream aofs("new_MM");
    for (int k=mmin;k<=mmax;k++)   // loop over movement periods
    {
      //for (int is=1;is<=1;is++)
      int ii=1;
      for (int is=1;is<=ns;is++)
      {
        if (pmsd)
        {
          offset=pmsd->rowoffset(is);
          rmin=pmsd->region_bounds(is,1);
          rmax=pmsd->region_bounds(is,2);
        }
        dvar_matrix td=Dad(k,1).sub(rmin,rmax).shift(1);

        for (int i=1;i<=nrr;i++)
        {
          td(i,i)=1;
          for (int j=1;j<=nrr;j++)
          {
            if (Dflags(i+offset,j)) 
            {
              dvariable tmp=0.0;
              tmp=diff_coffs(k,ii);
              td(i,i)+=tmp;
              td(j,i)-=tmp;
              ii++;
            }
          }
        }

        Dad(k,1).sub(rmin,rmax).shift(1)=inv(td);
        for (int j=2;j<=nage;j++)
        {
          Dad(k,j).sub(rmin,rmax).shift(1)=
            Dad(k,1).sub(rmin,rmax).shift(1);
        }
      }
    }
  }  //NMD_20Nov2018   - insert case for age-specific movement with pmsd_error
  else
  {
    if (pmsd)
    {
      pmsd_error();
      ad_exit(1);
    }
    for (int k=mmin;k<=mmax;k++)   // loop over movement periods
    {

  //  const MY_DOUBLE_TYPE ninv=1.0/nage;
      const MY_DOUBLE_TYPE ninv=1.0/(nage-1);

      for (int j=1;j<=nage;j++)
      {
        int ii=1;
        dvariable fj=-1.+2.*(j-1)*ninv;

        dvar_matrix& td=Dad(k,j); //NMD_20Nov2018 - indexed by period
        for (int i=1;i<=num_regions;i++)
        {
          td(i,i)=1;
          for (int jj=1;jj<=num_regions;jj++)
          {
            if (age_flags(91))
            {
              if (fj<-1.e-8)
              {
                fj=-pow(-fj,exp(diff_coffs3(k,jj)));
              }
              else if (fj>1.e-8)
              {
                fj=pow(fj,exp(diff_coffs3(k,jj)));
              }
	    }
            if (Dflags(i,jj)) 
            {
              dvariable tmp;
              if (!age_flags(91))
	      {  
                tmp=diff_coffs(k,ii)*exp(diff_coffs2(k,ii)*fj);
	      }
	      else
	      {
                if (fj<-1.e-8)
                {
                  fj=-pow(-fj,exp(diff_coffs3(k,ii)));
                }
                else if (fj>1.e-8)
                {
                  fj=pow(fj,exp(diff_coffs3(k,ii)));
                }
                tmp=diff_coffs(k,ii)*exp(diff_coffs2(k,ii)*fj);
	      }
              td(i,i)+=tmp;
              td(jj,i)-=tmp;
              ii++;
            }
          }
        }
        td=inv(td);
        {
          dvar_matrix tdtrans=trans(td);
          int mmin=tdtrans.indexmin();
          int mmax=tdtrans.indexmax();
          for (int i=mmin;i<=mmax;i++)
          {
            tdtrans(i)+=1.e-7;
            tdtrans(i)/=sum(tdtrans(i));
          }
          td=trans(tdtrans);
        }
        //cout <<"rshort3.cpp " << colsum(td) << endl;
      }
    }
  }  //NMD_20Nov2018
}

// this is the new method
void dvar_fish_stock_history::get_xdiff_parameterization(void)
{
  if (sum(Dflags)<=0)   //NMD_jan14-19
  {
    xdiff_coffs=0.0;
  }   //NMD_jan14-19
}


void dvar_fish_stock_history::get_diff_parameterization(void)
{
  if (sum(Dflags)>0)   //NMD_jan14-19
  {
    int mmin=xdiff_coffs.indexmin();
    int mmax=xdiff_coffs.indexmax();
    int nrows=mmax-mmin+1;
    dvar_vector v1=lineup(xdiff_coffs);  // lineup the rows of diff_coffs into
    dvar_vector tmp1=exp(v1*new_orthogonal_diffusion_matrix);
    dvar_matrix M=stack(tmp1,nrows);
    diff_coffs=M;
  }
  else
  {
    diff_coffs=0.0;
  }   //NMD_jan14-19
}

dmatrix get_orthogonal_diffusion_matrix(int n)
{
  dmatrix M(1,n,1,n);
  M.initialize();

  M(1)=1.0;
  for (int i=2;i<=n;i++)
  {
    M(i,i-1)=1.0;
  }

  //cout << setfixed() << setprecision(3) << M << endl;

  for (int i=1;i<=n;i++)
  {
    M(i)/=norm(1.e-30+M(i));
    for (int j=i+1;j<=n;j++)
    {
      M(j)-= (M(j)*M(i)) *M(i);
    }
  }
  return M;
}

void dvar_fish_stock_history::
  get_new_xdiff_parameterization(void)
{
  int mmin=xdiff_coffs.indexmin();
  int mmax=xdiff_coffs.indexmax();
//  __gdouble__=(double *) (&xdiff_coffs(1,1));  //NMD_24mar2023

  int ns=1;
  int nrr=num_regions;
  if (pmsd)
  {
    ns=pmsd->num_species;
    nrr=pmsd->num_real_regions;
  }
  if (!age_flags(89))
  {
    int offset=0;
    int rmin=1;
    int rmax=num_regions;
    //ofstream aofs("new_MM");
    for (int k=mmin;k<=mmax;k++)   // loop over movement periods
    {
      int ii=1;
      for (int is=1;is<=ns;is++)
      {
        if (pmsd)
        {
          offset=pmsd->rowoffset(is);
          rmin=pmsd->region_bounds(is,1);
          rmax=pmsd->region_bounds(is,2);
        }
        dvar_matrix MM(1,nrr,1,nrr);
        MM.initialize();
        dvar_matrix MM1(1,nrr,1,nrr);
        MM1.initialize();
        dvar_matrix MM2(1,nrr,1,nrr);
        MM2.initialize();
        dvar_matrix MM3(1,nrr,1,nrr);
        MM3.initialize();
//        cout << " BEFORE " << " norm2(xdiff_coffs)= "  << norm2(xdiff_coffs) 
//             << " norm(diff_coffs) = " << norm2(diff_coffs) << endl;
        // put the xdiff_coffs into the right place to use them
        // for this movement matrix
        int ii1=ii;   // need to use ii twice here 
                      // so first time use ii1
        for (int i=1;i<=nrr;i++)
        {
          for (int j=1;j<=nrr;j++)
          {
            if (Dflags(i+offset,j)) 
            {
              MM(j,i)=diff_coffs(k,ii1);
              ii1++;
            }
          }
        }
        for (int i=1;i<=nrr;i++)
        {
          for (int j=1;j<i;j++)
          {
            if (Dflags(i+offset,j)) 
            {
              MM1(i,j)=log(MM(i,j))+log(MM(j,i));
              MM1(j,i)=log(MM(i,j))-log(MM(j,i));
            }
          }
        }
        ii1=ii;   // need to use ii twice 
                      //here so first time use ii1
        for (int i=1;i<=nrr;i++)
        {
          for (int j=1;j<=i;j++)
          {
            if (Dflags(i+offset,j))
            { 
              y1diff_coffs(k,ii1)=MM(j,i)+MM(i,j);
              y2diff_coffs(k,ii1)=MM(j,i)/(MM(j,i)+MM(i,j));
              ii1++;
            }
          }
        }
        ii1=ii;   // need to use ii twice 
                      //here so first time use ii1
        for (int i=1;i<=nrr;i++)
        {
          for (int j=1;j<=nrr;j++)
          {
            if (Dflags(i+offset,j))
            { 
              xdiff_coffs(k,ii1)=MM1(j,i);
              ii1++;
            }
          }
        }
        ii1=ii;   // need to use ii twice here 
                     //so first time use ii1
        //iflag.initialize();
        for (int i=1;i<=nrr;i++)
        {
          for (int j=1;j<=nrr;j++)
          {
            if (Dflags(i+offset,j))
            {
              MM2(j,i)=xdiff_coffs(k,ii1);
            //  iflag(i,j)=1;
              ii1++;
            }
          }
        }

        for (int i=1;i<=nrr;i++)
        {
          MM3(i,i)=1.0;
          for (int j=1;j<=nrr;j++)
          {
            //if (iflag(i,j))
            if (Dflags(i+offset,j))
            {
              if (i<j)
              {
                MM3(i,i)+=exp(0.5*(MM2(i,j)+MM2(j,i)));
                MM3(j,i)-=exp(0.5*(MM2(i,j)+MM2(j,i)));
              }
              else
              {
                MM3(i,i)+=exp(0.5*(MM2(i,j)-MM2(j,i)));
                MM3(j,i)-=exp(0.5*(MM2(i,j)-MM2(j,i)));
              }
            }
          }
        }

        /*
        aofs << "MM3"  << k << endl;
        aofs << setprecision(3) << setw(7) 
             << MM3 << endl << endl;
        */

      }
    }
//    cout << " AFTER " << " norm(xdiff_coffs)= "  << norm2(xdiff_coffs) 
//         << " norm(diff_coffs) = " << norm2(diff_coffs) << endl;
  }  //NMD_20Nov2018   - insert case for age-specific movement with pmsd_error
  else
  {
    if (pmsd)
    {
      pmsd_error();
      ad_exit(1);
    }
  }
}

void dvar_fish_stock_history::
  get_old_diff_parameterization(void)
{
  if (age_flags(184)==1)
  {
    int mmin=xdiff_coffs.indexmin();
    int mmax=xdiff_coffs.indexmax();
  
    int ns=1;
    int nrr=num_regions;
    if (pmsd)
    {
      ns=pmsd->num_species;
      nrr=pmsd->num_real_regions;
    }
    if (!age_flags(89))
    {
      int offset=0;
      int rmin=1;
      int rmax=num_regions;
      //ofstream aofs("new_MM");
      for (int k=mmin;k<=mmax;k++)   // loop over movement periods
      {
        int ii=1;
        for (int is=1;is<=ns;is++)
        {
          if (pmsd)
          {
            offset=pmsd->rowoffset(is);
            rmin=pmsd->region_bounds(is,1);
            rmax=pmsd->region_bounds(is,2);
          }
          dvar_matrix MMy(1,nrr,1,nrr);
          dvar_matrix MM2(1,nrr,1,nrr);
          MM2.initialize();
          MMy.initialize();
          dvar_matrix MM3(1,nrr,1,nrr);
          MM3.initialize();
          // put the xdiff_coffs into the right place to use them
          // for this movement matrix
          int ii1=ii;   // need to use ii twice here 
                        // so first time use ii1
          for (int i=1;i<=nrr;i++)
          {
            for (int j=1;j<=nrr;j++)
            {
              if (Dflags(i+offset,j))
              {
                MM2(j,i)=xdiff_coffs(k,ii1);
                ii1++;
              }
            }
          }
  
          ii1=ii;   // need to use ii twice here 
                        // so first time use ii1
          for (int i=1;i<=nrr;i++)
          {
            for (int j=1;j<i;j++)
            {
              if (Dflags(i+offset,j))
              {
                MMy(j,i)=y1diff_coffs(k,ii1)*y2diff_coffs(k,ii1);
                MMy(i,j)=y1diff_coffs(k,ii1)*(1.0-y2diff_coffs(k,ii1));
                ii1++;
              }
            }
          }
  
          for (int i=1;i<=nrr;i++)
          {
            MM3(i,i)=1.0;
            for (int j=1;j<=nrr;j++)
            {
              //if (iflag(i,j))
              if (Dflags(i+offset,j))
              {
                if (i<j)
                {
                  MM3(i,i)+=exp(0.5*(MM2(i,j)+MM2(j,i)));
                  MM3(j,i)-=exp(0.5*(MM2(i,j)+MM2(j,i)));
                }
                else
                {
                  MM3(i,i)+=exp(0.5*(MM2(i,j)-MM2(j,i)));
                  MM3(j,i)-=exp(0.5*(MM2(i,j)-MM2(j,i)));
                }
              }
            }
          }
          ii1=ii;   // need to use ii twice 
                        //here so first time use ii1
          for (int i=1;i<=nrr;i++)
          {
            for (int j=1;j<=nrr;j++)
            {
              if (Dflags(i+offset,j))
              { 
                diff_coffs(k,ii1)=-MM3(j,i);
                ii1++;
              }
            }
          }
        }
      }
    }  //NMD_20Nov2018   - insert case for age-specific movement with pmsd_error
    else
    {
      if (pmsd)
      {
        pmsd_error();
        ad_exit(1);
      }
    }
  }
  else if (age_flags(184)==2)
  {

    // for af184==2 active parameters are zdiff_coffs. 
    // Convert them to diff_coffs for calculation of the movement matrices
    movement_coffs_option_2_to_0();
    /*
    int mmin=xdiff_coffs.indexmin();
    int mmax=xdiff_coffs.indexmax();
  
    int ns=1;
    int nrr=num_regions;
    if (pmsd)
    {
      ns=pmsd->num_species;
      nrr=pmsd->num_real_regions;
    }
    if (!age_flags(89))
    {
      int offset=0;
      int rmin=1;
      int rmax=num_regions;
      //ofstream aofs("new_MM");
      for (int k=mmin;k<=mmax;k++)   // loop over movement periods
      {
        int ii=1;
        for (int is=1;is<=ns;is++)
        {
          if (pmsd)
          {
            offset=pmsd->rowoffset(is);
            rmin=pmsd->region_bounds(is,1);
            rmax=pmsd->region_bounds(is,2);
          }
          dvar_matrix MMy(1,nrr,1,nrr);
          dvar_matrix MM2(1,nrr,1,nrr);
          MM2.initialize();
          MMy.initialize();
          dvar_matrix MM3(1,nrr,1,nrr);
          MM3.initialize();
          // put the xdiff_coffs into the right place to use them
          // for this movement matrix
          int ii1=ii;   // need to use ii twice here 
                        // so first time use ii1
          for (int i=1;i<=nrr;i++)
          {
            for (int j=1;j<=nrr;j++)
            {
              if (Dflags(i+offset,j))
              {
                MM2(j,i)=xdiff_coffs(k,ii1);
                ii1++;
              }
            }
          }
  
          ii1=ii;   // need to use ii twice here 
                        // so first time use ii1
          for (int i=1;i<=nrr;i++)
          {
            for (int j=1;j<i;j++)
            {
              if (Dflags(i+offset,j))
              {
                MMy(j,i)=y1diff_coffs(k,ii1)*y2diff_coffs(k,ii1);
                MMy(i,j)=y1diff_coffs(k,ii1)*(1.0-y2diff_coffs(k,ii1));
                ii1++;
              }
            }
          }
  
          for (int i=1;i<=nrr;i++)
          {
            MM3(i,i)=1.0;
            for (int j=1;j<=nrr;j++)
            {
              //if (iflag(i,j))
              if (Dflags(i+offset,j))
              {
                if (i<j)
                {
                  MM3(j,i)-=MMy(j,i)*MMy(i,j);
                  MM3(i,i)+=MM3(j,i);
                }
                else
                {
                  MM3(j,i)-=MMy(j,i)*(1.0-MMy(i,j));
                  MM3(i,i)+=MM3(j,i);
                }
              }
            }
          }
          ii1=ii;   // need to use ii twice 
                        //here so first time use ii1
          for (int i=1;i<=nrr;i++)
          {
            for (int j=1;j<=nrr;j++)
            {
              if (Dflags(i+offset,j))
              { 
                diff_coffs(k,ii1)=-MM3(j,i);
                ii1++;
              }
            }
          }
        }
      }
    }  //NMD_20Nov2018   - insert case for age-specific movement with pmsd_error
    else
    {
      if (pmsd)
      {
        pmsd_error();
        ad_exit(1);
      }
    }
  */
  }
}
