/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include <fvar.hpp>

void set_value_inv(const dvar_vector& x,const dvector& _v,const int& _ii,
  ivector mflags,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s)
{
  int& ii=(int&) _ii;
  dvector& v=(dvector&) _v;
  int min=x.indexmin();
  int max=x.indexmax();
  for (int i=min;i<=max;i++)
  {
    if (mflags(i))
    {
      v(ii++)=boundpin(x(i),fmin,fmax)*s;
    }
  }
}
 
  void set_value_inv(const dvar_matrix & _w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,ivector&  flags,imatrix& mflags,
    const ivector& group)
  {
    ADUNCONST(dvar_matrix,w)  //NMD_27apr2018
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      if (flags(kin)) 
          set_value_inv(w(kin),x,ii,mflags(kin),fmin,fmax,s);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          w(key(i))=w(key(i-1));
        }
        else
        {
          kin=key(i);
//          if (flags(i))
          if (flags(kin))    //DF_9jan2019
          {
            set_value_inv(w(kin),x,ii,mflags(kin),fmin,fmax,s);
          }
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) 
          set_value_inv(w(i),x,ii,mflags(i),fmin,fmax,s);
      }
    }
  }

void set_value(const dvar_vector& _v,const dvar_vector& x,const int& _ii,
  ivector& mflags,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const prevariable& pen,MY_DOUBLE_TYPE s)
{
  RETURN_ARRAYS_INCREMENT();
  int& ii=(int&) _ii;
  dvar_vector& v=(dvar_vector&) _v;
  int min=v.indexmin();
  int max=v.indexmax();
  for (int i=min;i<=max;i++)
  {
    if (mflags(i))
    {
      v(i)=boundp(x(ii++),fmin,fmax,pen,s);
    }
  }
  RETURN_ARRAYS_DECREMENT();
}


  void set_value(const dvar_matrix & _w,const dvar_vector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const prevariable& pen,
    MY_DOUBLE_TYPE s,ivector&  flags,imatrix& mflags,
    const ivector& group)
  {
    ADUNCONST(dvar_matrix,w)  //NMD_27apr2018
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      if (flags(kin)) 
          set_value(w(kin),x,ii,mflags(kin),fmin,fmax,pen,s);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          w(key(i))=w(key(i-1));
        }
        else
        {
          kin=key(i);
//          if (flags(i))
          if (flags(kin))    //DF_9jan2019
          {
            set_value(w(kin),x,ii,mflags(kin),fmin,fmax,pen,s);
          }
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) 
          set_value(w(i),x,ii,mflags(i),fmin,fmax,pen,s);
      }
    }
  }
  
  void set_value_inv(const dvar_matrix & _w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,ivector&  flags,
    const ivector& group)
  {
    ADUNCONST(dvar_matrix,w)  //NMD_27apr2018
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
        if (group(key(i))==group(key(i-1)))
        {
          w(key(i))=w(key(i-1));
        }
        else
        {
          kin=key(i);
//          if (flags(i))
          if (flags(kin))    //DF_9jan2019
          {
            set_value_inv(w(kin),x,ii,fmin,fmax,s);
          }
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
  int size(const dvar_vector& v,const ivector& flags)
  {
    int ii=0;
    int mmin=v.indexmin();
    int mmax=v.indexmax(); 
    if (mmin != flags.indexmin() || 
      mmax != flags.indexmax() )
    {
      cerr << "shape miismatch " << endl;
      ad_exit(1);
    }
    for (int i=mmin;i<=mmax;i++)
    {
      if(flags(i)) ii++;
    }
    return ii;
  }
  int num_active(const dvar_matrix& w,ivector& flags,imatrix& mflags,
    const ivector& group)
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
      {
        nv+=size(w(kin),mflags(kin));
      }
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))!=group(key(i-1)))
        {
          kin=key(i);
          if (flags(kin))
            nv+=size(w(kin),mflags(kin));
        } 
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i))
          nv+=size_count(w(i));
      }
    }
    return nv;
  }
  int operator != (const ivector& _v,const ivector& _w)
  {
    ADUNCONST(ivector,v)
    ADUNCONST(ivector,w)
    int mmin=v.indexmin();
    int mmax=v.indexmax();
    
    if (mmin != w.indexmin() || mmax != w.indexmax())
    {
      cerr << "size error in ivector compare" << endl;
      ad_exit(1);
    }
    int ir=0;
    for (int i=mmin;i<=mmax;i++)
    {
      if (v(i) != w(i))
      {
        ir=1;
        break;
      }
    }
    return ir;
  }

  void zero_effdev_sanity_check(dvar_matrix& w,ivector& flags,imatrix& mflags,
    const ivector& group)
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
      int fl=flags(kin);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))!=group(key(i-1)))
        {
          kin=key(i);
          fl=flags(kin);
        } 
        else
        {
          int kinold=kin; //YT 2018-06-19
          kin=key(i);
          if (flags(kin) != fl)
          {
            cerr << "fish flags(4) do not satisfy grouping coherence" 
                 << endl;
            ad_exit(1);
          }
          if (mflags(kin) != mflags(kinold)) 
          {
            cerr << "zero fdev flags do not satisfy grouping coherence" 
                 << endl;
            ad_exit(1);
          }
        }
      }
    }
  }
