/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include <fvar.hpp>
#include "scbd.hpp"

int getdim(int ff57,int ff61,int nage)
{
  switch (ff57)
  {
  case 0:
    return nage;
  case 1:
    return 2;
  case 2:
    return 3;
  case 3:
    return ff61;
  default:
    cerr << "this cant happen" << endl;
    ad_exit(1);
  }
}
    
int logic_equal(int flag,int flag1)
{
  if (  (flag && flag1) || (!flag && !flag1) )
    return 1;
  else
    return 0;
}

  dvar_matrix colsub(const dvar_matrix& _M,int lb,int ub);

void set_value_partial(const dvar3_array& _w,const dvar_vector& x,
  const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,const ivector& flags,const ivector& group,
  MY_DOUBLE_TYPE scale,ivector& block_ptr,ivector& ff57,int nage,ivector& ff61,
  ivector& sel_block_index)
  {
    ADUNCONST(dvar3_array,w)
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      int flag_value=flags(kin);
      int block_flag=block_ptr(kin);  // is this fishery blocked?
      int dim=getdim(ff57(kin),ff61(kin),nage);
      if (flags(kin) && block_flag) 
      {
        const dvar_matrix& M=
          colsub(w(kin).sub(1,sel_block_index(block_flag)),1,dim);
        set_value(M,x,ii,fmin,fmax,pen,scale);
      }
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          int dim1=getdim(ff57(key(i)),ff61(kin),nage);
          int dim2=getdim(ff57(key(i-1)),ff61(kin),nage);
          if (dim1!=dim2)
          {
            cerr << "Error -- grouped initial parameters have unequal ff57"
                 << endl;
            ad_exit(1);
          }
          if (!logic_equal(block_ptr(key(i)),block_flag))
          {
            cerr << "Error -- grouped initial parameters have unequal"
                " block flags (ff71) "
                " for fishery " << key(i) << " and " << kin <<  endl <<
                " block flags are " << block_ptr(key(i)) << " and " 
                << block_ptr(kin) << endl;
            ad_exit(1);
          }
          if (flags(key(i))!=flag_value)
          {
            cerr << "Error -- grouped initial parameters have unequal flags"
                 << endl;
            ad_exit(1);
          }
          if (block_flag)
          {
            const dvar_matrix& _M1=
              colsub(w(key(i)).sub(1,sel_block_index(block_flag)),1,dim1);
            const dvar_matrix& M2=
              colsub(w(key(i-1)).sub(1,sel_block_index(block_flag)),1,dim2);
            ADUNCONST(dvar_matrix,M1)
            M1=M2;
          }
        }
        else
        {
          kin=key(i);
          flag_value=flags(kin);
          block_flag=block_ptr(kin);  // is this fishery blocked?
          dim=getdim(ff57(kin),ff61(kin),nage);
          flag_value=flags(kin);
          if (flags(kin) && block_flag) 
          {
            const dvar_matrix& M=
              colsub(w(kin).sub(1,sel_block_index(block_flag)),1,dim);
            set_value(M,x,ii,fmin,fmax,pen,scale);
          }
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        int block_flag=block_ptr(i);  // is this fishery blocked?
        int dim=getdim(ff57(i),ff61(i),nage);
        if (flags(i) && block_flag) 
        {
          const dvar_matrix& M=
            colsub(w(i).sub(1,sel_block_index(block_flag)),1,dim);
          set_value(M,x,ii,fmin,fmax,pen,scale);
        }
      }
    }
  }

void set_value_inv_partial(const dvar3_array& _w,const dvector& x,
  const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const ivector& flags,const ivector& group,
  MY_DOUBLE_TYPE scale,ivector& block_ptr,ivector& ff57,int nage,ivector& ff61,
  ivector& sel_block_index)
  {
    ADUNCONST(dvar3_array,w)
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      int flag_value=flags(kin);
      int block_flag=block_ptr(kin);  // is this fishery blocked?
      int dim=getdim(ff57(kin),ff61(kin),nage);
      if (flags(kin) && block_flag) 
      {
        const dvar_matrix& M=
          colsub(w(kin).sub(1,sel_block_index(block_flag)),1,dim);

        set_value_inv(M,x,ii,fmin,fmax,scale);
      }
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          int dim1=getdim(ff57(key(i)),ff61(kin),nage);
          int dim2=getdim(ff57(key(i-1)),ff61(kin),nage);
          if (dim1!=dim2)
          {
            cerr << "Error -- grouped initial parameters have unequal ff57"
                 << endl;
            ad_exit(1);
          }
          if (flags(key(i))!=flag_value)
          {
            cerr << "Error -- grouped initial parameters have unequal flags"
                 << endl;
            ad_exit(1);
          }
        }
        else
        {
          kin=key(i);
          flag_value=flags(kin);
          block_flag=block_ptr(kin);  // is this fishery blocked?
          dim=getdim(ff57(kin),ff61(kin),nage);
          flag_value=flags(kin);
          if (flags(kin) && block_flag) 
          {
            const dvar_matrix& M=
              colsub(w(kin).sub(1,sel_block_index(block_flag)),1,dim);
            set_value_inv(M,x,ii,fmin,fmax,scale);
          }
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        int block_flag=block_ptr(i);  // is this fishery blocked?
        int dim=getdim(ff57(i),ff61(i),nage);
        if (flags(i) && block_flag) 
        {
          const dvar_matrix& M=
            colsub(w(i).sub(1,sel_block_index(block_flag)),1,dim);
          set_value_inv(M,x,ii,fmin,fmax,scale);
        }
      }
    }
  }

int size_count_partial(const dvar3_array& _w,
  const ivector& flags,const ivector& group,
  ivector& block_ptr,ivector& ff57,int nage,ivector& ff61,
  ivector& sel_block_index)
  {
    int itmp=0;
    ADUNCONST(dvar3_array,w)
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      int flag_value=flags(kin);
      int block_flag=block_ptr(kin);  // is this fishery blocked?
      int dim=getdim(ff57(kin),ff61(kin),nage);
      if (flags(kin) && block_flag) 
      {
        const dvar_matrix& M=
          colsub(w(kin).sub(1,sel_block_index(block_flag)),1,dim);

        itmp+=size_count(M);
      }
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          int dim1=getdim(ff57(key(i)),ff61(kin),nage);
          int dim2=getdim(ff57(key(i-1)),ff61(kin),nage);
          if (dim1!=dim2)
          {
            cerr << "Error -- grouped initial parameters have unequal ff57"
                 << endl;
            ad_exit(1);
          }
          if (flags(key(i))!=flag_value)
          {
            cerr << "Error -- grouped initial parameters have unequal flags"
                 << endl;
            ad_exit(1);
          }
        }
        else
        {
          kin=key(i);
          flag_value=flags(kin);
          block_flag=block_ptr(kin);  // is this fishery blocked?
          dim=getdim(ff57(kin),ff61(kin),nage);
          flag_value=flags(kin);
          if (flags(kin) && block_flag) 
          {
            const dvar_matrix& M=
              colsub(w(kin).sub(1,sel_block_index(block_flag)),1,dim);
            itmp+=size_count(M);
          }
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        int block_flag=block_ptr(i);  // is this fishery blocked?
        int dim=getdim(ff57(i),ff61(i),nage);
        if (flags(i) && block_flag) 
        {
          const dvar_matrix& M=
            colsub(w(i).sub(1,sel_block_index(block_flag)),1,dim);
          itmp+=size_count(M);
        }
      }
    }
    return itmp;
  }

