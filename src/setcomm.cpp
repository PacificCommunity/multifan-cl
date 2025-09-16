/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include <fvar.hpp>
#include "newmult.hpp"
#include "scbd.hpp"

dvar_vector rowstack(const dvar_matrix& M)
{
  int mmin=M.indexmin();
  int mmax=M.indexmax();
  int ns=0;
  int i;
  for (i=mmin;i<=mmax;i++)
  {
    ns+=M(i).indexmax()-M(i).indexmin()+1;
  }
  dvar_vector tmp(1,ns);
  int offset=0;
  ns=0;
  for (i=mmin;i<=mmax;i++)
  {
    ns=M(i).indexmax()-M(i).indexmin()+1;
    dvar_vector v=M(i);
    tmp(1+offset,ns+offset)=v.shift(1+offset);
    offset+=ns;
  }
  return tmp;
}
   
void rowunstack(const dvar_vector& _w,const dvar_matrix& _M)
{
  ADUNCONST(dvar_vector,w)
  ADUNCONST(dvar_matrix,M)
  int mmin=M.indexmin();
  int mmax=M.indexmax();
  int ns=0;
  int i;
  int offset=0;
  for (i=mmin;i<=mmax;i++)
  {
    ns=M(i).indexmax()-M(i).indexmin()+1;
    M(i)=w(1+offset,ns+offset).shift(M(i).indexmin());
    offset+=ns;
  }
}

dvector rowstack(const dmatrix& M)
{
  int mmin=M.indexmin();
  int mmax=M.indexmax();
  int ns=0;
  int i;
  for (i=mmin;i<=mmax;i++)
  {
    ns+=M(i).indexmax()-M(i).indexmin()+1;
  }
  dvector tmp(1,ns);
  int offset=0;
  ns=0;
  for (i=mmin;i<=mmax;i++)
  {
    ns=M(i).indexmax()-M(i).indexmin()+1;
    dvector v=M(i);
    tmp(1+offset,ns+offset)=v.shift(1+offset);
    offset+=ns;
  }
  return tmp;
}

ivector rowstack(const imatrix& M)
{
  int mmin=M.indexmin();
  int mmax=M.indexmax();
  int ns=0;
  int i;
  for (i=mmin;i<=mmax;i++)
  {
    ns+=M(i).indexmax()-M(i).indexmin()+1;
  }
  ivector tmp(1,ns);
  int offset=0;
  ns=0;
  for (i=mmin;i<=mmax;i++)
  {
    ns=M(i).indexmax()-M(i).indexmin()+1;
    ivector v=M(i);
    tmp(1+offset,ns+offset)=v.shift(1+offset);
    offset+=ns;
  }
  return tmp;
}
   

/*
void set_value(const prevariable& _x,const dvar_vector& v,const int& _ii, 
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax,const dvariable& _fpen,MY_DOUBLE_TYPE s)
{
  int& ii=(int&)_ii; 
  prevariable& x=(prevariable&) _x;
  dvariable& fpen=(dvariable&) _fpen;
  x=boundp(v(ii++),log(fmin),log(fmax),fpen,s);
}
*/


void set_value_exp(const prevariable& _x,const dvar_vector& v,const int& _ii, 
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax,const dvariable& _fpen,MY_DOUBLE_TYPE s)
{
  int& ii=(int&)_ii; 
  prevariable& x=(prevariable&) _x;
  dvariable& fpen=(dvariable&) _fpen;
  x=exp(boundp(v(ii++),log(fmin),log(fmax),fpen,s));
}

  void check_group_for_holes(const ivector& v)
  {
    ivector w=sort(v);
    int mmin=v.indexmin();
    int mmax=v.indexmax();
    for (int i=mmin;i<mmax;i++)
    {
      if (w(i+1)>w(i)+1)
      {
        cerr << "Error in grouping vector -- there is a hole" << endl;
        cerr << "The values for the grouping vector are" << endl
           << setw(3) << v << endl;
        ad_exit(1);
      }
    }
  }
 
  void check_group_flag_sanity(const ivector& flag,const ivector& group)
  {
    int mmin=group.indexmin();
    int mmax=group.indexmax();
    if (flag.indexmin() != mmin || flag.indexmax() != mmax ) 
    {
      cerr << "Error -- grouiping and active flags have differnent shape" 
           << endl;
        ad_exit(1);
    }
    ivector tmp(min(group),max(group));
    tmp.initialize();
    for (int i=mmin;i<mmax;i++)
    {
      if (!tmp(group(i)))
      {
        if (flag(i)==1)
          tmp(group(i))=1;
        else   
          tmp(group(i))=2;
      }
      else
      {
        if (flag(i)==1)
        {
          if (tmp(group(i))==2)
          {
            cerr << "Error -- Different active flag values for grouped"
              " parameters" << endl;
            ad_exit(1);
            
          }
        }
        else   
        {
          if (tmp(group(i))==1)
          {
            cerr << "Error -- Different active flag values for grouped"
              " parameters" << endl;
            ad_exit(1);
            
          }
        }
      }
    }
  } 
 //   void set_value_inv(const dvar_vector& w,const dvector& x,const int& ii,
 //     MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const ivector& flags,const ivector& group)
 //   {
 //     int mmin=w.indexmin();
 //     int mmax=w.indexmax();
 //     if (sum(flags))
 //     {
 //       if (sum(group))
 //       {
 //         ivector key(mmin,mmax);
 //         sort(group,key);
 //         int kin=key(mmin);
 //         int flag_value=flags(kin);
 //         if (flags(kin)) set_value_inv(w(kin),x,ii,fmin,fmax);
 //         for (int i=mmin+1;i<=mmax;i++)
 //         {
 //           if (group(key(i))==group(key(i-1)))
 //           {
 //             if (flags(key(i))!=flag_value)
 //             {
 //               cerr << "Error -- grouped initial parameters have unequal flags"
 //                    << endl;
 //               exit(1);
 //             }
 //           }
 //           else
 //           {
 //             kin=key(i);
 //             flag_value=flags(kin);
 //             if (flags(kin)) set_value_inv(w(kin),x,ii,fmin,fmax);
 //           }
 //         }
 //       }
 //       else
 //       {
 //         for (int i=mmin;i<=mmax;i++)
 //         {
 //           if (flags(i)) set_value_inv(w(i),x,ii,fmin,fmax);
 //         }
 //       }
 //     }
 //   }
 // 

  int num_active(const dvar_vector& w,const ivector& flags,const ivector& group)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (mmin != flags.indexmin() || mmax != flags.indexmax()
      ||  mmin != group.indexmin() || mmax != group.indexmax())
    {
      cerr << "B incompatible array bounds in num_active() " << endl;
      cerr << w.indexmin() << " " << w.indexmax() << endl;
      cerr << flags.indexmin() << " " << flags.indexmax() << endl;
      cerr << group.indexmin() << " " << group.indexmax() << endl;
      exit(1);
    }
    int nv=0;
    if (sum(flags))
    {
      int maxg=max(group);
      if (maxg)
      {
        ivector key(mmin,mmax);
        sort(group,key); 
        /*
        {
          for (int i=mmin;i<=mmax;i++) 
          {
            cout <<"setcomm.cpp " << group(key(i)) << endl;
          }
        }
        */
        int kin=key(mmin);
        int flag_value=flags(kin);
        if (flags(kin)) nv++;
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
            if (flags(kin)) nv++;
          } 
        }
      }
      else
      {
        for (int i=mmin;i<=mmax;i++)
        {
          if (flags(i)) nv++;
        }
      }
    }
    return nv;
  }

  int num_active(const dvar_matrix& _mw,const imatrix& mflags,const imatrix& mgroup)
  {
    dvar_vector w=rowstack(_mw);
    ivector flags=rowstack(mflags);
    ivector group=rowstack(mgroup);
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (mmin != flags.indexmin() || mmax != flags.indexmax()
      ||  mmin != group.indexmin() || mmax != group.indexmax())
    {
      cerr << "B incompatible array bounds in num_active() " << endl;
      cerr << w.indexmin() << " " << w.indexmax() << endl;
      cerr << flags.indexmin() << " " << flags.indexmax() << endl;
      cerr << group.indexmin() << " " << group.indexmax() << endl;
      exit(1);
    }
    int nv=0;
    if (sum(flags))
    {
      int maxg=max(group);
      if (maxg)
      {
        ivector key(mmin,mmax);
        sort(group,key); 
        /*
        {
          for (int i=mmin;i<=mmax;i++) 
          {
            cout <<"setcomm.cpp " << group(key(i)) << endl;
          }
        }
        */
        int kin=key(mmin);
        int flag_value=flags(kin);
        if (flags(kin)) nv++;
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
            if (flags(kin)) nv++;
          } 
        }
      }
      else
      {
        for (int i=mmin;i<=mmax;i++)
        {
          if (flags(i)) nv++;
        }
      }
    }
    return nv;
  }

  int num_active(const dvar_matrix& _mw,int flag)
  {
    ADUNCONST(dvar_matrix,mw)
    int nv=0;
    if (flag)
    {
      int mmin=mw.indexmin();
      int mmax=mw.indexmax();
      for (int i=mmin;i<=mmax;i++)
      {
        nv+=mw(i).indexmax()-mw(i).indexmin()+1;
      }
    }
    return nv;
  }

  int num_active(const dvar_vector& w,const ivector& flags)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (mmin != flags.indexmin() || mmax != flags.indexmax())
    {
      cerr << "C incompatible array bounds in num_active() " << endl;
      exit(1);
    }
    int nv=0;
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) nv++;
      }
    }
    return nv;
  }
  int num_active(const dvar_vector& w,int flags)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    int nv=0;
    if (mmax>=mmin && flags) nv+=(mmax-mmin+1);
    return nv;
  }




  void set_value(const dvar_vector& _w,const dvar_vector& x,const int& _ii,const ivector& onsw,
    const ivector& gpsw)
  {
    dvar_vector& w=(dvar_vector&) _w;
    int& ii=(int&) _ii;
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (!sum(gpsw))
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (onsw(i))
        {
          w(i)=x(ii++);
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (onsw(i))
        {
          w(i)=x(ii+gpsw(i));
        }
      }
      ii+=max(gpsw);
    }
  }

  void set_value_inv(const dvar_vector& _w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,int flags,MY_DOUBLE_TYPE& scale)
  {
    dvar_vector& w=(dvar_vector&) _w;
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags) set_value_inv(w(i),x,ii,fmin,fmax,scale);
      }
    }
  }

  void set_value_inv(const dvar_vector& w,const dvector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    const ivector& flags,MY_DOUBLE_TYPE s)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) set_value_inv(w(i),x,ii,fmin,fmax,s);
      }
    }
  }

void set_value(const dvar_vector& _w,const dvar_vector& x,const int& _ii,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,const dvariable& pen,const ivector& flags)
{
  dvar_vector& w=(dvar_vector&) _w;
  int& ii= (int&) _ii;
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  {
    for (int i=mmin;i<=mmax;i++)
    {
      if (flags(i)) w(i)=boundp(x(ii++),fmin,fmax,pen);
    }
  }
}

void set_value(const dvar_vector& _w,const dvar_vector& x,const int& _ii,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,const dvariable& pen,const ivector& flags,MY_DOUBLE_TYPE s)
{
  dvar_vector& w=(dvar_vector&) _w;
  int& ii= (int&) _ii;
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  {
    for (int i=mmin;i<=mmax;i++)
    {
      if (flags(i)) w(i)=boundp(x(ii++),fmin,fmax,pen,s);
    }
  }
}

/*
//void set_value_inv(const prevariable& ,const dvector& x,const int& ii,MY_DOUBLE_TYPE fmin,
//  MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE scale);

  void set_value_inv(const dvector& w,const dvector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    const ivector& flags,const ivector& group,MY_DOUBLE_TYPE scale)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (sum(group))
    {
      ivector key(mmin,mmax);
      sort(group,key);
      int kin=key(mmin);
      MY_DOUBLE_TYPE w_value=value(w(kin));
      int flag_value=flags(kin);
      if (flags(kin)) set_value_inv(w(kin),x,ii,fmin,fmax,scale);
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))==group(key(i-1)))
        {
          if (value(w(key(i)))!=w_value)
          {
            cerr << "Error -- grouped initial parameters have unequal values"
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
          w_value=value(w(kin));
          if (flags(kin)) set_value_inv(w(kin),x,ii,fmin,fmax);
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) set_value_inv(w(i),x,ii,fmin,fmax);
      }
    }
  }

void set_value(const dvar_vector& w,const dvar_vector& x,const int& ii,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,const dvariable& pen,const ivector& flags,const ivector& group,MY_DOUBLE_TYPE scale)
{
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  int maxg=max(group);
  if (maxg)
  {
    ii--;
    for (int i=mmin;i<=mmax;i++)
    {
      if (flags(i)) w(i)=boundp(x(ii+group(i)),fmin,fmax,pen,scale);
    }
    ii+=max(group)+1;
  }
  else
  {
    for (int i=mmin;i<=mmax;i++)
    {
      if (flags(i)) w(i)=boundp(x(ii++),fmin,fmax,pen,scale);
    }
  }
}
  void set_value(const dvar_vector& w,const dvar_vector& x,const int& ii,const ivector& onsw,
    const ivector& gpsw,MY_DOUBLE_TYPE scale)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (!sum(gpsw))
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (onsw(i))
        {
          w(i)=x(ii++)/scale;
        }
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (onsw(i))
        {
          w(i)=x(ii+gpsw(i))/scale;
        }
      }
      ii+=max(gpsw);
    }
  }

  void set_value_inv(const dvector& w,const dvector& x,const int& ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,
    const ivector& flags,MY_DOUBLE_TYPE scale)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (flags(i)) set_value_inv(w(i),x,ii,fmin,fmax);
      }
    }
  }

void set_value(const dvar_vector& w,const dvar_vector& x,const int& ii,MY_DOUBLE_TYPE fmin,
  MY_DOUBLE_TYPE fmax,const dvariable& pen,const ivector& flags,MY_DOUBLE_TYPE scale)
{
  int mmin=w.indexmin();
  int mmax=w.indexmax();
  {
    for (int i=mmin;i<=mmax;i++)
    {
      if (flags(i)) w(i)=boundp(x(ii++),fmin,fmax,pen,scale);
    }
  }
}
/* ********************************************************************* */
  void set_value_inv(const dvar_vector& w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const ivector& flags,const ivector& group,MY_DOUBLE_TYPE s)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    
    if (sum(flags))
    {
      if (sum(group))
      {
        check_group_for_holes(group);
        check_group_flag_sanity(flags,group);
        ivector key(mmin,mmax);
        sort(group,key);
        int kin=key(mmin);
        MY_DOUBLE_TYPE w_value=value(w(kin));
        int flag_value=flags(kin);
        if (flags(kin)) set_value_inv(w(kin),x,ii,fmin,fmax,s);
        for (int i=mmin+1;i<=mmax;i++)
        {
          if (group(key(i))==group(key(i-1)))
          {
            if (value(w(key(i)))!=w_value)
            {
              cerr << "Error -- grouped initial parameters have unequal values"
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
            w_value=value(w(kin));
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
  }

void set_value_inv_exp(const prevariable& x,const dvector& _v,const int& _ii,
    MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s)
{
  dvector& v=(dvector&)(_v);
  int& ii=(int&)(_ii);
  v(ii++)=boundpin(log(x),log(fmin),log(fmax),s);
}


  void set_value_inv_exp(const dvar_matrix& _mw,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax, const imatrix& mflags,const imatrix& mgroup,
    MY_DOUBLE_TYPE s)
  {
    dvar_vector w=rowstack(_mw);
    ivector flags=rowstack(mflags);
    ivector group=rowstack(mgroup);

    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (sum(flags))
    {
      if (sum(group))
      {
        check_group_for_holes(group);
        check_group_flag_sanity(flags,group);
        ivector key(mmin,mmax);
        sort(group,key);
        int kin=key(mmin);
        int flag_value=flags(kin);
        MY_DOUBLE_TYPE w_value=value(w(kin));
        if (flags(kin)) set_value_inv_exp(w(kin),x,ii,fmin,fmax,s);
        for (int i=mmin+1;i<=mmax;i++)
        {
          if (group(key(i))==group(key(i-1)))
          {
            if (value(w(key(i)))!=w_value)
            {
              cerr << "Error -- grouped initial parameters have unequal values"
                   << endl;
              ad_exit(1);
            }
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
            w_value=value(w(kin));
            if (flags(kin)) set_value_inv_exp(w(kin),x,ii,fmin,fmax,s);
          }
        }
      }
      else
      {
        for (int i=mmin;i<=mmax;i++)
        {
          if (flags(i)) set_value_inv_exp(w(i),x,ii,fmin,fmax,s);
        }
      }
    }
    rowunstack(w,_mw);
  }

  void set_value_inv_exp(const dvar_vector& w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax, const ivector& flags,const ivector& group,
    MY_DOUBLE_TYPE s)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (sum(flags))
    {
      if (sum(group))
      {
        check_group_for_holes(group);
        check_group_flag_sanity(flags,group);
        ivector key(mmin,mmax);
        sort(group,key);
        int kin=key(mmin);
        int flag_value=flags(kin);
        MY_DOUBLE_TYPE w_value=value(w(kin));
        if (flags(kin)) set_value_inv(w(kin),x,ii,fmin,fmax,s);
        for (int i=mmin+1;i<=mmax;i++)
        {
          if (group(key(i))==group(key(i-1)))
          {
            if (value(w(key(i)))!=w_value)
            {
              cerr << "Error -- grouped initial parameters have unequal values"
                   << endl;
              ad_exit(1);
            }
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
            w_value=value(w(kin));
            if (flags(kin)) set_value_inv_exp(w(kin),x,ii,fmin,fmax,s);
          }
        }
      }
      else
      {
        for (int i=mmin;i<=mmax;i++)
        {
          if (flags(i)) set_value_inv_exp(w(i),x,ii,fmin,fmax,s);
        }
      }
    }
  }

  void set_value(const dvar_vector& _w,const dvar_vector& x,const int& _ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const dvariable& pen,const ivector& flags,
    const ivector& group,MY_DOUBLE_TYPE s)
  {
    dvar_vector& w=(dvar_vector&) _w;
    int& ii= (int&) _ii;
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (sum(flags))
    {
      if (sum(group))
      {
        ivector key(mmin,mmax);
        sort(group,key);
        int kin=key(mmin);
        int flag_value=flags(kin);
        if (flags(kin)) set_value(w(kin),x,ii,fmin,fmax,pen,s);
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
            if (flags(key(i))) w(key(i))=w(kin);
          }
          else
          {
            kin=key(i);
            flag_value=flags(kin);
            if (flags(kin)) set_value(w(kin),x,ii,fmin,fmax,pen,s);
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
  }

  void set_value_exp(const dvar_vector& _w,const dvar_vector& x,const int& _ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const dvariable& pen,const ivector& flags,
    const ivector& group,MY_DOUBLE_TYPE s)
  {
    dvar_vector& w=(dvar_vector&) _w;
    int& ii= (int&) _ii;
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (sum(flags))
    {
      if (sum(group))
      {
        ivector key(mmin,mmax);
        sort(group,key);
        int kin=key(mmin);
        int flag_value=flags(kin);
        if (flags(kin)) set_value(w(kin),x,ii,fmin,fmax,pen,s);
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
            if (flags(key(i))) w(key(i))=w(kin);
          }
          else
          {
            kin=key(i);
            flag_value=flags(kin);
            if (flags(kin)) set_value_exp(w(kin),x,ii,fmin,fmax,pen,s);
          }
        }
      }
      else
      {
        for (int i=mmin;i<=mmax;i++)
        {
          if (flags(i)) set_value_exp(w(i),x,ii,fmin,fmax,pen,s);
        }
      }
    }
  }


void xgroup_manager::safe_allocate(const ivector& flags,
  const ivector& groups)
{
  /*
  if (!allocated())
  {
    check_group_for_holes(groups);
    check_group_flag_sanity(flags,groups);
    int mmin=flags.indexmin();
    int mmax=flags.indexmax();
    ivector key(mmin,mmax);
    sort(groups,key);
    int ig=1;
    int kin=key(mmin);
    // get the number of groups
    for (int i=mmin+1;i<=mmax;i++)
    {
      if (groups(key(i))!=groups(key(i-1)))
      {
        ig++;
      }
    }
    ngroups=ig;
    ng.allocate(1,ngroups); 
    group_ptr.allocate(mmin,mmax); 
    inv_group_ptr.allocate(1,ngroups); 
    ng.initialize();
    group_ptr.initialize(); 
    inv_group_ptr.initialize(); 
    ig=1;
    group_ptr(key(mmin))=ig;
    inv_group_ptr(ig)=key(mmin); 
    for (int i=mmin+1;i<=mmax;i++)
    {
      if (groups(key(i))!=groups(key(i-1)))
      {
        ig++;
        inv_group_ptr(ig)=key(i); 
      } 
      cout << "UUUU " << ig << endl;
      group_ptr(key(i))=ig;
      ng(ig)++;
    }
  }
  */
  cout << "Leaving UUU " << endl;
}
  
