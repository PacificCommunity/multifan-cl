/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include <fvar.hpp>

  int size_count(const dvar_matrix& w,int flags,const ivector& group)
  {
    int n=0;
    int mmin=w.rowmin();
    int mmax=w.rowmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      //if (flags) set_value_inv(w(kin),x,ii,fmin,fmax,s);
      if (flags) n+=size_count(w(kin));
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))!=group(key(i-1)))
        {
          kin=key(i);
          //if (flags) set_value_inv(w(kin),x,ii,fmin,fmax,s);
          if (flags) n+=size_count(w(kin));
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        //if (flags) set_value_inv(w(i),x,ii,fmin,fmax,s);
        if (flags) n+=size_count(w(i));
      }
    }
    return n;
  }

  void set_value_inv(const dmatrix& w,const dvector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    int flags,const ivector& group)
  {
    int mmin=w.rowmin();
    int mmax=w.rowmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      if (flags) set_value_inv(w(kin),x,ii,fmin,fmax);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))!=group(key(i-1)))
        {
          kin=key(i);
          if (flags) set_value_inv(w(kin),x,ii,fmin,fmax);
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags) set_value_inv(w(i),x,ii,fmin,fmax);
      }
    }
  }


  void set_value_inv(const dmatrix& w,const dvector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,
    int flags,const ivector& group)
  {
    int mmin=w.rowmin();
    int mmax=w.rowmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      if (flags) set_value_inv(w(kin),x,ii,fmin,fmax,s);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))!=group(key(i-1)))
        {
          kin=key(i);
          if (flags) set_value_inv(w(kin),x,ii,fmin,fmax,s);
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags) set_value_inv(w(i),x,ii,fmin,fmax,s);
      }
    }
  }


  void set_value_inv(const dvar_matrix& w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,int flags,const ivector& group)
  {
    int mmin=w.rowmin();
    int mmax=w.rowmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      if (flags) set_value_inv(w(kin),x,ii,fmin,fmax,s);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))!=group(key(i-1)))
        {
          kin=key(i);
          if (flags) set_value_inv(w(kin),x,ii,fmin,fmax,s);
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags) set_value_inv(w(i),x,ii,fmin,fmax,s);
      }
    }
  }

  void set_value(const dvar_matrix& _w,const dvar_vector& x,
    const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    const dvariable& pen,int flags,const ivector& group)
  {
    dvar_matrix& w=(dvar_matrix&) _w;
    int mmin=w.rowmin();
    int mmax=w.rowmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      if (flags) set_value(w(kin),x,ii,fmin,fmax,pen);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          w(key(i))=w(key(i-1));
        }
        else
        {
          kin=key(i);
          if (flags) set_value(w(kin),x,ii,fmin,fmax,pen);
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags) set_value(w(i),x,ii,fmin,fmax,pen);
      }
    }
  }

dvariable boundp(const prevariable& x, MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax)
{
  if (gradient_structure::Hybrid_bounded_flag==0)
  {
    dvariable t,y;
    MY_DOUBLE_TYPE diff=fmax-fmin;
    const MY_DOUBLE_TYPE l4=log(4.0);
    dvariable ss=0.4999999999999999*sin(x*1.57079632679489661)+0.50;
    t=fmin + diff*ss;
    return(t);
  }
  else
  {
    MY_DOUBLE_TYPE diff=fmax-fmin;
    dvariable t,y;
    if (x>-20)
    {
      y=1.0/(1+exp(-x));
    }
    else
    {
      dvariable u=exp(x);
      y=u/(1.0+u);
    }
    t=fmin + diff*y;
    return(t);
  }
}
dvariable boundp(const prevariable& x, MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s)
{
  return boundp(x/s,fmin,fmax);
}

void set_value(const dvar_vector& x,_CONST dvar_vector& v, const int& _ii,
  MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const MY_DOUBLE_TYPE s)
{
  int& ii = (int&) _ii;
  if (!(!(x)))
  {
    int min=x.indexmin();
    int max=x.indexmax();
    for (int i=min;i<=max;i++)
    {
      ((dvar_vector&)(x))(i)=boundp(v(ii++),fmin,fmax,s);
    }
  }
} 
  void set_value(const dvar_matrix& _w,const dvar_vector& x,
    const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    const MY_DOUBLE_TYPE s,int flags,const ivector& group)
  {
    dvar_matrix& w=(dvar_matrix&) _w;
    int mmin=w.rowmin();
    int mmax=w.rowmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      if (flags) set_value(w(kin),x,ii,fmin,fmax,s);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          w(key(i))=w(key(i-1));
        }
        else
        {
          kin=key(i);
          if (flags) set_value(w(kin),x,ii,fmin,fmax,s);
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags) set_value(w(i),x,ii,fmin,fmax,s);
      }
    }
  }

  int num_active(const dmatrix& w,int flags,const ivector& group)
  {
    int mmin=w.rowmin();
    int mmax=w.rowmax();
    if (mmin != group.indexmin() || mmax != group.indexmax())
    {
      cerr << "F incompatible array bounds in num_active() " << endl;
      exit(1);
    }
    int nv=0;
    int maxg=max(group);
    if (maxg)
    {
      ivector key(mmin,mmax);
      sort(group,key); 
      int kin=key(mmin);
      if (flags) nv+=w(kin).size();
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))!=group(key(i-1)))
        {
          kin=key(i);
          if (flags) nv+=w(kin).size();
        } 
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags) nv+=w(i).size();
      }
    }
    return nv;
  }


  void set_value(const dvar_matrix& _w,const dvar_vector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    const dvariable& pen,MY_DOUBLE_TYPE s,int flags,const ivector& group)
  {
    dvar_matrix& w=(dvar_matrix&) _w;
    int mmin=w.rowmin();
    int mmax=w.rowmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      if (flags) set_value(w(kin),x,ii,fmin,fmax,pen,s);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          w(key(i))=w(key(i-1));
        }
        else
        {
          kin=key(i);
          if (flags) set_value(w(kin),x,ii,fmin,fmax,pen,s);
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags) set_value(w(i),x,ii,fmin,fmax,pen,s);
      }
    }
  }
  int num_active(const dvar_matrix& w,int flags,const ivector& group)
  {
    int mmin=w.rowmin();
    int mmax=w.rowmax();
    if (mmin != group.indexmin() || mmax != group.indexmax())
    {
      cerr << "G incompatible array bounds in num_active() " << endl;
      exit(1);
    }
    int nv=0;
    if (flags)
    {
      int maxg=max(group);
      if (maxg)
      {
        ivector key(mmin,mmax);
        sort(group,key); 
        int kin=key(mmin);
        nv+=w(kin).size();
        for (int i=mmin+1;i<=mmax;i++)
        {
          if (group(key(i))!=group(key(i-1)))
          {
            kin=key(i);
            nv+=w(kin).size();
          } 
        }
      }
      else
      {
        for (int i=mmin;i<=mmax;i++)
        {
          nv+=w(i).size();
        }
      }
    }
    return nv;
  }

