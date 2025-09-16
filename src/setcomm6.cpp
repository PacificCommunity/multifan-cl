/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include <fvar.hpp>
  int size(const dmatrix& M)
  {
    int tmp=0;
    int mmin=M.indexmin();
    int mmax=M.indexmax();
    for (int i=mmin;i<=mmax;i++)
    {
      tmp+=M(i).size();
    }
    return tmp;
  }

  int size(const dvar_matrix& M)
  {
    int tmp=0;
    int mmin=M.indexmin();
    int mmax=M.indexmax();
    for (int i=mmin;i<=mmax;i++)
    {
      tmp+=M(i).size();
    }
    return tmp;
  }


  int size_count(const dvar3_array& w,int flags,const ivector& group)
  {
    int n=0;
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      if (flags) n+=size_count(w(kin));
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))!=group(key(i-1)))
        {
          kin=key(i);
          if (flags) n+=size_count(w(kin));
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags) n+=size_count(w(i));
      }
    }
    return n;
  }

  void set_value_inv(const d3_array& w,const dvector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    int flags,const ivector& group)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
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


  void set_value_inv(const d3_array& w,const dvector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,
    int flags,const ivector& group)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
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


  void set_value_inv(const dvar3_array& w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,int flags,const ivector& group)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
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

  void set_value(const dvar3_array& _w,const dvar_vector& x,
    const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    const dvariable& pen,int flags,const ivector& group)
  {
    dvar3_array& w=(dvar3_array&) _w;
    int mmin=w.indexmin();
    int mmax=w.indexmax();
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

  void set_value(const dvar4_array& _w,const dvar_vector& x,
    const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    const dvariable& pen,int flags,const ivector& group)
  {
    dvar4_array& w=(dvar4_array&) _w;
    int mmin=w.indexmin();
    int mmax=w.indexmax();
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


  void set_value(const dvar3_array& _w,const dvar_vector& x,
    const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    const MY_DOUBLE_TYPE s,int flags,const ivector& group)
  {
    dvar3_array& w=(dvar3_array&) _w;
    int mmin=w.indexmin();
    int mmax=w.indexmax();
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

  int num_active(const d3_array& w,int flags,const ivector& group)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
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
      if (flags) nv+=size(w(kin));
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))!=group(key(i-1)))
        {
          kin=key(i);
          if (flags) nv+=size(w(kin));
        } 
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags) nv+=size(w(i));
      }
    }
    return nv;
  }


  void set_value(const dvar3_array& _w,const dvar_vector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    const dvariable& pen,MY_DOUBLE_TYPE s,int flags,const ivector& group)
  {
    dvar3_array& w=(dvar3_array&) _w;
    int mmin=w.indexmin();
    int mmax=w.indexmax();
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
  int num_active(const dvar3_array& w,int flags,const ivector& group)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
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
        nv+=size(w(kin));
        for (int i=mmin+1;i<=mmax;i++)
        {
          if (group(key(i))!=group(key(i-1)))
          {
            kin=key(i);
            nv+=size(w(kin));
          } 
        }
      }
      else
      {
        for (int i=mmin;i<=mmax;i++)
        {
          nv+=size(w(i));
        }
      }
    }
    return nv;
  }

  void set_value(const dvar3_array& _w,const dvar_vector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    const dvariable& pen,MY_DOUBLE_TYPE s,ivector&  flags,const ivector& group)
  {
    dvar3_array& w=(dvar3_array&) _w;
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      if (flags(kin)) 
      {
        set_value(w(kin),x,ii,fmin,fmax,pen,s);
      }
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          w(key(i))=w(key(i-1));
        }
        else
        {
          kin=key(i);
          if (flags(i)) 
          {
            set_value(w(kin),x,ii,fmin,fmax,pen,s);
          }
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) set_value(w(i),x,ii,fmin,fmax,pen,s);
      }
    }
  }

  void set_value(const dvar4_array& _w,const dvar_vector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    const dvariable& pen,MY_DOUBLE_TYPE s,ivector&  flags,const ivector& group)
  {
    dvar4_array& w=(dvar4_array&) _w;
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      if (flags(kin)) 
         set_value(w(kin),x,ii,fmin,fmax,pen,s);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          w(key(i))=w(key(i-1));
        }
        else
        {
          kin=key(i);
          if (flags(kin)) 
            set_value(w(kin),x,ii,fmin,fmax,pen,s);
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) 
          set_value(w(i),x,ii,fmin,fmax,pen,s);
      }
    }
  }


  int num_active(const dvar3_array& w,ivector& flags,const ivector& group)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (mmin != group.indexmin() || mmax != group.indexmax())
    {
      cerr << "G incompatible array bounds in num_active() " << endl;
      exit(1);
    }
    int nv=0;

    int maxg=max(group);
    if (maxg)
    {
      ivector key(mmin,mmax);
      sort(group,key); 
      int kin=key(mmin);
      if (flags(kin))
        nv+=size(w(kin));
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))!=group(key(i-1)))
        {
          kin=key(i);
          if (flags(kin))
            nv+=size(w(kin));
        } 
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i))
          nv+=size(w(i));
      }
    }
    return nv;
  }

  void set_value_inv(const dvar3_array& w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,ivector&  flags,const ivector& group)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      if (flags(kin)) set_value_inv(w(kin),x,ii,fmin,fmax,s);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))!=group(key(i-1)))
        {
          kin=key(i);
          if (flags(kin)) set_value_inv(w(kin),x,ii,fmin,fmax,s);
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) set_value_inv(w(i),x,ii,fmin,fmax,s);
      }
    }
  }

  void set_value_inv(const dvar4_array& w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,ivector&  flags,const ivector& group)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      if (flags(kin)) 
          set_value_inv(w(kin),x,ii,fmin,fmax,s);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))!=group(key(i-1)))
        {
          kin=key(i);
          if (flags(kin)) 
            set_value_inv(w(kin),x,ii,fmin,fmax,s);
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) 
          set_value_inv(w(i),x,ii,fmin,fmax,s);
      }
    }
  }

  int size_count(const dvar3_array& w,ivector&  flags,const ivector& group)
  {
    int n=0;
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      if (flags(kin)) n+=size_count(w(kin));
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))!=group(key(i-1)))
        {
          kin=key(i);
          if (flags(kin)) n+=size_count(w(kin));
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) n+=size_count(w(i));
      }
    }
    return n;
  }

  int size_count(const dvar4_array& w,ivector&  flags,const ivector& group)
  {
    int n=0;
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      if (flags(kin)) n+=size_count(w(kin));
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))!=group(key(i-1)))
        {
          kin=key(i);
          if (flags(kin)) n+=size_count(w(kin));
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) n+=size_count(w(i));
      }
    }
    return n;
  }

