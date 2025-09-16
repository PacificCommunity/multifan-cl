/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include <fvar.hpp>
#include "scbd.hpp"


void set_value_inv_partial(const dmatrix& w,const dvector& x,const int& ii,
  const ivector& range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const ivector& flags,
  const ivector& group,MY_DOUBLE_TYPE scale)
  {
    int mmin=w.rowmin();
    int mmax=w.rowmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      int flag_value=flags(kin);
      if (flags(kin)) set_value_inv_partial(w(kin),x,ii,range(kin),fmin,fmax,
        scale);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          if (flags(key(i))!=flag_value)
          {
            cerr << "Error -- grouped initial parameters have unequal flags"
                 << endl;
            exit(1);
          }
        }
        else
        {
          kin=key(i);
          flag_value=flags(kin);
          if (flags(kin)) set_value_inv_partial(w(kin),x,ii,range(kin),fmin,
            fmax,scale);
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) set_value_inv_partial(w(i),x,ii,range(i),fmin,fmax,scale);
      }
    }
  }

void set_value_inv_partial(const dmatrix& w,const dvector& x,const int& ii,
  int range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE scale)
{
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    set_value_inv_partial(w(i),x,ii,range,fmin,fmax,scale);
  }
}



void set_value_inv_partial(const d3_array& w,const dvector& x,const int& ii,
  const ivector& range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const ivector& flags,
  const ivector& group,MY_DOUBLE_TYPE scale)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      int flag_value=flags(kin);
      if (flags(kin)) set_value_inv_partial(w(kin),x,ii,range(kin),fmin,fmax,
        scale);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          if (flags(key(i))!=flag_value)
          {
            cerr << "Error -- grouped initial parameters have unequal flags"
                 << endl;
            exit(1);
          }
        }
        else
        {
          kin=key(i);
          flag_value=flags(kin);
          if (flags(kin)) set_value_inv_partial(w(kin),x,ii,range(kin),fmin,
            fmax,scale);
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) set_value_inv_partial(w(i),x,ii,range(i),fmin,fmax,scale);
      }
    }
  }

void set_value_inv_partial(const d3_array& w,const dvector& x,const int& ii,
  const int range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,int flags, MY_DOUBLE_TYPE scale)
{
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  if (flags) 
  {
    for (int i=mmin;i<=mmax;i++)
    {
      set_value_inv_partial(w(i),x,ii,range,fmin,fmax,scale);
    }
  }
}
void set_value_inv_partial(const dvar_matrix& w,const dvector& x,const int& ii,
  int range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE scale);

void set_value_inv_partial(const dvar3_array& w,const dvector& x,const int& ii,
  const int range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE scale)
{
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    set_value_inv_partial(w(i),x,ii,range,fmin,fmax,scale);
  }
}

void set_value_inv_partial(const dvar4_array& w,const dvector& x,const int& ii,
  const int range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,int flags,MY_DOUBLE_TYPE scale)
{
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  if (flags) 
  {
    for (int i=mmin;i<=mmax;i++)
    {
      set_value_inv_partial(w(i),x,ii,range,fmin,fmax,scale);
    }
  }
}

void set_value_inv_partial(const dvar4_array& w,const dvector& x,const int& ii,
  const int range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE scale)
{
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    set_value_inv_partial(w(i),x,ii,range,fmin,fmax,scale);
  }
}


void set_value_partial(const dvar_matrix& _w,const dvar_vector& x,
  const int& ii,const ivector& range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,const ivector& group,MY_DOUBLE_TYPE scale)
  {
    dvar_matrix& w=(dvar_matrix&) _w;
    int mmin=w.rowmin();
    int mmax=w.rowmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      set_value_partial(w(kin),x,ii,range(kin),fmin,fmax,pen,
        scale);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          w(key(i))=w(key(i-1));
        }
        else
        {
          kin=key(i);
          set_value_partial(w(kin),x,ii,range(kin),fmin,fmax,pen,scale);
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        set_value_partial(w(i),x,ii,range(i),fmin,fmax,pen,scale);
      }
    }
  }

void set_value_inv_partial(const dvar_matrix& _w,const dvector& x,
  const int& ii,const ivector& range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const ivector& group,MY_DOUBLE_TYPE scale)
  {
    dvar_matrix& w=(dvar_matrix&) _w;
    int mmin=w.rowmin();
    int mmax=w.rowmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      set_value_inv_partial(w(kin),x,ii,range(kin),fmin,fmax,
        scale);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          w(key(i))=w(key(i-1));
        }
        else
        {
          kin=key(i);
          set_value_inv_partial(w(kin),x,ii,range(kin),fmin,fmax,scale);
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        set_value_inv_partial(w(i),x,ii,range(i),fmin,fmax,scale);
      }
    }
  }
  

void set_value_partial(const dvar_matrix& _w,const dvar_vector& x,
  const int& ii,const ivector& range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,const ivector& flags,const ivector& group,MY_DOUBLE_TYPE scale)
  {
    dvar_matrix& w=(dvar_matrix&) _w;
    int mmin=w.rowmin();
    int mmax=w.rowmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      int flag_value=flags(kin);
      if (flags(kin)) set_value_partial(w(kin),x,ii,range(kin),fmin,fmax,pen,
        scale);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          w(key(i))=w(key(i-1));
          if (flags(key(i))!=flag_value)
          {
            cerr << "Error -- grouped initial parameters have unequal flags"
                 << endl;
            exit(1);
          }
        }
        else
        {
          kin=key(i);
          flag_value=flags(kin);
          if (flags(kin)) set_value_partial(w(kin),x,ii,range(kin),fmin,
            fmax,pen,scale);
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) set_value_partial(w(i),x,ii,range(i),fmin,fmax,pen,
          scale);
      }
    }
  }
void set_value_partial(const dvar3_array& _w,const dvar_vector& x,
  const int& ii,const ivector& range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,const ivector& flags,const ivector& group,MY_DOUBLE_TYPE scale)
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
      if (flags(kin)) set_value_partial(w(kin),x,ii,range(kin),fmin,fmax,pen,
        scale);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          w(key(i))=w(key(i-1));
          if (flags(key(i))!=flag_value)
          {
            cerr << "Error -- grouped initial parameters have unequal flags"
                 << endl;
            exit(1);
          }
        }
        else
        {
          kin=key(i);
          flag_value=flags(kin);
          if (flags(kin)) set_value_partial(w(kin),x,ii,range(kin),fmin,
            fmax,pen,scale);
        }
      }
    }
    else
    {
      cout << mmin << " " << mmax << endl;
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) set_value_partial(w(i),x,ii,range(i),fmin,fmax,pen,
          scale);
      }
    }
  }

void set_value_partial(const dvar_matrix& _w,const dvar_vector& x,
  const int& ii,int range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,const ivector& flags,const ivector& group,MY_DOUBLE_TYPE scale)
  {
    dvar_matrix& w=(dvar_matrix&) _w;
    int mmin=w.rowmin();
    int mmax=w.rowmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      int flag_value=flags(kin);
      if (flags(kin)) set_value_partial(w(kin),x,ii,range,fmin,fmax,pen,
        scale);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          w(key(i))=w(key(i-1));
          if (flags(key(i))!=flag_value)
          {
            cerr << "Error -- grouped initial parameters have unequal flags"
                 << endl;
            exit(1);
          }
        }
        else
        {
          kin=key(i);
          flag_value=flags(kin);
          if (flags(kin)) set_value_partial(w(kin),x,ii,range,fmin,
            fmax,pen,scale);
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) set_value_partial(w(i),x,ii,range,fmin,fmax,pen,
          scale);
      }
    }
  }

void set_value_partial(const dvar3_array& _w,const dvar_vector& x,
  const int& ii,int range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,const ivector& flags,const ivector& group,MY_DOUBLE_TYPE scale)
  {
    ADUNCONST(dvar3_array,w)
    int mmin=w.rowmin();
    int mmax=w.rowmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      int flag_value=flags(kin);
      if (flags(kin)) set_value_partial(w(kin),x,ii,range,fmin,fmax,pen,
        scale);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          w(key(i))=w(key(i-1));
          if (flags(key(i))!=flag_value)
          {
            cerr << "Error -- grouped initial parameters have unequal flags"
                 << endl;
            exit(1);
          }
        }
        else
        {
          kin=key(i);
          flag_value=flags(kin);
          if (flags(kin)) set_value_partial(w(kin),x,ii,range,fmin,
            fmax,pen,scale);
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) set_value_partial(w(i),x,ii,range,fmin,fmax,pen,
          scale);
      }
    }
  }


void set_value_partial(const dvar3_array& _w,const dvar_vector& x,
  const int& ii,const ivector& range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,const ivector& flags,MY_DOUBLE_TYPE scale)
{
  ADUNCONST(dvar3_array,w)
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  cout << mmin << " " << mmax << endl;
  for (int i=mmin;i<=mmax;i++)
  {
    if (flags(i)) set_value_partial(w(i),x,ii,range(i),fmin,fmax,pen,
      scale);
  }
}

void set_value_partial(const dvar3_array& _w,const dvar_vector& x,
  const int& ii,const ivector& range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,MY_DOUBLE_TYPE scale)
{
  ADUNCONST(dvar3_array,w)
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  cout << mmin << " " << mmax << endl;
  for (int i=mmin;i<=mmax;i++)
  {
    set_value_partial(w(i),x,ii,range(i),fmin,fmax,pen,scale);
  }
}

void set_value_partial(const dvar3_array& _w,const dvar_vector& x,
  const int& ii,const int range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,MY_DOUBLE_TYPE scale)
{
  ADUNCONST(dvar3_array,w)
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  cout << mmin << " " << mmax << endl;
  for (int i=mmin;i<=mmax;i++)
  {
    set_value_partial(w(i),x,ii,range,fmin,fmax,pen,scale);
  }
}

/*
void set_value_partial(const dvar4_array& _w,const dvar_vector& x,
  const int& ii,const ivector& range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,const ivector& flags,MY_DOUBLE_TYPE scale)
{
  ADUNCONST(dvar4_array,w)
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  cout << mmin << " " << mmax << endl;
  for (int i=mmin;i<=mmax;i++)
  {
    if (flags(i)) set_value_partial(w(i),x,ii,range,fmin,fmax,pen,
      scale);
  }
}
*/
void set_value_partial(const dvar4_array& _w,const dvar_vector& x,
  const int& ii,const int range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,MY_DOUBLE_TYPE scale)
{
  ADUNCONST(dvar4_array,w)
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  cout << mmin << " " << mmax << endl;
  for (int i=mmin;i<=mmax;i++)
  {
    set_value_partial(w(i),x,ii,range,fmin,fmax,pen,scale);
  }
}

void set_value_partial(const dvar4_array& _w,const dvar_vector& x,
  const int& ii,const int range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,int flags,MY_DOUBLE_TYPE scale)
{
  ADUNCONST(dvar4_array,w)
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  cout << mmin << " " << mmax << endl;
  if (flags)
  { 
    for (int i=mmin;i<=mmax;i++)
    {
      set_value_partial(w(i),x,ii,range,fmin,fmax,pen,scale);
    }
  }
}

/*
void set_value_partial(const dvar3_array& _w,const dvar_vector& x,
  const int& ii,const ivector& range,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
  const dvariable& pen,const ivector& flags,const ivector& group,
    MY_DOUBLE_TYPE scale)
{
  ADUNCONST(dvar3_array,w)
  int mmin=w.rowmin();
  int mmax=w.rowmax();
  if (sum(group))
  {
    ivector key(mmin,mmax);
    sort(group,key);
    int kin=key(mmin);
    int flag_value=flags(kin);
    if (flags(kin)) set_value_partial(w(kin),x,ii,range(kin),fmin,fmax,pen,
      scale);
    for (int i=mmin+1;i<=mmax;i++)
    {
      if (group(key(i))==group(key(i-1)))
      {
        w(key(i))=w(key(i-1));
        if (flags(key(i))!=flag_value)
        {
          cerr << "Error -- grouped initial parameters have unequal flags"
               << endl;
          exit(1);
        }
      }
      else
      {
        kin=key(i);
        flag_value=flags(kin);
        if (flags(kin)) set_value_partial(w(kin),x,ii,range(kin),fmin,
          fmax,pen,scale);
      }
    }
  }
  else
  {
    for (int i=mmin;i<=mmax;i++)
    {
      if (flags(i)) set_value_partial(w(i),x,ii,range(i),fmin,fmax,pen,
        scale);
    }
  }
}
*/

