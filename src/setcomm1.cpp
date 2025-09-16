/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include <fvar.hpp>

  int num_active(const dmatrix& w,const ivector& flags,const ivector& group)
  {
    int mmin=w.rowmin();
    int mmax=w.rowmax();
    if (mmin != flags.indexmin() || mmax != flags.indexmax()
      ||  mmin != group.indexmin() || mmax != group.indexmax())
    {
      cerr << "D incompatible array bounds in num_active() " << endl;
      ad_exit(1);
    }
    int nv=0;
    int maxg=max(group);
    if (maxg)
    {
      int maxf=max(flags);
      if (maxf)
      {
        for (int i=mmin;i<=mmax;i++)
        {
          if (flags(i)==0)
          {
            cerr << "At present for common (grouped) parameters"
                 " all active flags must be on or off " << endl;
          }
        }  
      }
      
      ivector key(mmin,mmax);
      sort(group,key); 
      int kin=key(mmin);
      int flag_value=flags(kin);
      if (flags(kin)) nv+=w(kin).size();
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
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
          if (flags(kin)) nv+=w(kin).size();
        } 
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) nv+=w(i).size();
      }
    }
    return nv;
  }

  int num_active_partial(const dvar_matrix& w,const ivector& flags,
    const ivector& group,const ivector& range)
  {
    int mmin=w.rowmin();
    int mmax=w.rowmax();
    if (mmin != flags.indexmin() || mmax != flags.indexmax()
      ||  mmin != group.indexmin() || mmax != group.indexmax())
    {
      cerr << "E incompatible array bounds in num_active() " << endl;
      ad_exit(1);
    }
    int nv=0;
    int maxg=max(group);
    if (maxg)
    {
      int maxf=max(flags);
      if (maxf)
      {
        for (int i=mmin;i<=mmax;i++)
        {
          if (flags(i)==0)
          {
            cerr << "At present for common (grouped) parameters"
                 " all active flags must be on or off " << endl;
          }
        }  
      }
      
      ivector key(mmin,mmax);
      sort(group,key); 
      int kin=key(mmin);
      int flag_value=flags(kin);
      if (flags(kin)) nv+=size_count_partial(w(kin),range(kin));
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
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
          if (flags(kin)) nv+=size_count_partial(w(kin),range(kin));
        } 
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) nv+=size_count_partial(w(i),range(i));
      }
    }
    return nv;
  }

  int num_active_partial(const dvar_matrix& w,const ivector& group,
    const ivector& range)
  {
    int mmin=w.rowmin();
    int mmax=w.rowmax();
    if (mmin != group.indexmin() || mmax != group.indexmax())
    {
      cerr << "E incompatible array bounds in num_active() " << endl;
      ad_exit(1);
    }
    int nv=0;
    int maxg=max(group);
    if (maxg)
    {
      ivector key(mmin,mmax);
      sort(group,key); 
      int kin=key(mmin);
      nv+=size_count_partial(w(kin),range(kin));
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))!=group(key(i-1)))
        {
          kin=key(i);
          nv+=size_count_partial(w(kin),range(kin));
        } 
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        nv+=size_count_partial(w(i),range(i));
      }
    }
    return nv;
  }

  int num_active_partial(const dvar3_array& w,const ivector& flags,
    const ivector& group,const ivector& range)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (mmin != flags.indexmin() || mmax != flags.indexmax()
      ||  mmin != group.indexmin() || mmax != group.indexmax())
    {
      cerr << "E incompatible array bounds in num_active() " << endl;
      ad_exit(1);
    }
    int nv=0;
    int maxg=max(group);
    if (maxg)
    {
      int maxf=max(flags);
      if (maxf)
      {
        for (int i=mmin;i<=mmax;i++)
        {
          if (flags(i)==0)
          {
            cerr << "At present for common (grouped) parameters"
                 " all active flags must be on or off " << endl;
          }
        }  
      }
      
      ivector key(mmin,mmax);
      sort(group,key); 
      int kin=key(mmin);
      int flag_value=flags(kin);
      if (flags(kin)) nv+=size_count_partial(w(kin),range(kin));
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
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
          if (flags(kin)) nv+=size_count_partial(w(kin),range(kin));
        } 
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) nv+=size_count_partial(w(i),range(i));
      }
    }
    return nv;
  }

  int num_active_partial(const dvar3_array& w,int flags,int range)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    int nv=0;
    if (flags) 
    {
      for (int i=mmin;i<=mmax;i++)
      {
        nv+=size_count_partial(w(i),range);
      }
    }
    return nv;
  }
  /*
  int num_active_partial(const dvar4_array& w,int flags,int range)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    int nv=0;
    if (flags) 
    {
      for (int i=mmin;i<=mmax;i++)
      {
        nv+=size_count_partial(w(i),range);
      }
    }
    return nv;
  }
  */
  int num_active_partial(const dvar4_array& w,int range)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    int nv=0;
    for (int i=mmin;i<=mmax;i++)
    {
      nv+=size_count_partial(w(i),range);
    }
    return nv;
  }
