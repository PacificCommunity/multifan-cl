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

  void check_group_for_holes(const ivector& v);
  void check_group_flag_sanity(const ivector& flag,const ivector& group);


void set_value_exp(const prevariable& _x,const dvar_vector& v,const int& _ii, 
  MY_DOUBLE_TYPE fmin, MY_DOUBLE_TYPE fmax,const dvariable& _fpen)
{
  int& ii=(int&)_ii; 
  prevariable& x=(prevariable&) _x;
  dvariable& fpen=(dvariable&) _fpen;
  x=exp(boundp(v(ii++),log(fmin),log(fmax),fpen));
}


  void set_value_inv(const dvar_vector& w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const ivector& flags,const ivector& group)
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
        if (flags(kin)) set_value_inv(w(kin),x,ii,fmin,fmax);
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
  }
  
  void set_value_inv(const dvar_matrix& _mw,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const imatrix& mflags,const imatrix& mgroup)
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
        if (flags(kin)) set_value_inv(w(kin),x,ii,fmin,fmax);
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
    rowunstack(w,_mw);
  }
  
    void set_value(const dvar_vector& _w,const dvar_vector& x,const int& _ii,
      MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const dvariable& pen,const ivector& flags,
      const ivector& group)
    {
      int i;
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
          /*
          for (i=mmin+1;i<=mmax;i++)
          {
            cout << group(key(i)) << " ";
          }
          */
          int kin=key(mmin);
          int flag_value=flags(kin);
          if (flags(kin)) set_value(w(kin),x,ii,fmin,fmax,pen);
          for (i=mmin+1;i<=mmax;i++)
          {
            if (group(key(i))==group(key(i-1)))
            {
              //if (flags(i)) w(key(i))=w(kin);
              //  DF  1Mar19
              if (flags(key(i))) w(key(i))=w(kin);
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
              if (flags(kin)) set_value(w(kin),x,ii,fmin,fmax,pen);
            }
          }
        }
        else
        {
          for (int i=mmin;i<=mmax;i++)
          {
            if (flags(i)) set_value(w(i),x,ii,fmin,fmax,pen);
          }
        }
      }
    }

    void set_value(const dvar_matrix& _mw,const dvar_vector& x,const int& _ii,
      MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const dvariable& pen,const imatrix& mflags,
      const imatrix& mgroup)
    {
      int i;
      //dvar_matrix& mw=(dvar_matrix&) _mw;
      dvar_vector w=rowstack(_mw);
      ivector flags=rowstack(mflags);
      ivector group=rowstack(mgroup);
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
          if (flags(kin)) set_value(w(kin),x,ii,fmin,fmax,pen);
          for (i=mmin+1;i<=mmax;i++)
          {
            if (group(key(i))==group(key(i-1)))
            {
              //if (flags(i)) w(key(i))=w(kin);
              //  DF  1Mar19
              if (flags(key(i))) w(key(i))=w(kin);
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
              if (flags(kin)) set_value(w(kin),x,ii,fmin,fmax,pen);
            }
          }
        }
        else
        {
          for (int i=mmin;i<=mmax;i++)
          {
            if (flags(i)) set_value(w(i),x,ii,fmin,fmax,pen);
          }
        }
      }
      rowunstack(w,_mw);
    }

    void set_value_exp(const dvar_matrix& _mw,const dvar_vector& x,
      const int& _ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const dvariable& pen,
      const imatrix& mflags,const imatrix& mgroup)
    {
      int i;
      //dvar_matrix& mw=(dvar_matrix&) _mw;
      dvar_vector w=rowstack(_mw);
      ivector flags=rowstack(mflags);
      ivector group=rowstack(mgroup);
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
          if (flags(kin)) set_value_exp(w(kin),x,ii,fmin,fmax,pen);
          for (i=mmin+1;i<=mmax;i++)
          {
            if (group(key(i))==group(key(i-1)))
            {
              //if (flags(i)) w(key(i))=w(kin);
              //  DF  1Mar19
              if (flags(key(i))) w(key(i))=w(kin);
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
              if (flags(kin)) set_value_exp(w(kin),x,ii,fmin,fmax,pen);
            }
          }
        }
        else
        {
          for (int i=mmin;i<=mmax;i++)
          {
            if (flags(i)) set_value_exp(w(i),x,ii,fmin,fmax,pen);
          }
        }
      }
      rowunstack(w,_mw);
    }

    void set_value_exp(const dvar_matrix& _mw,const dvar_vector& x,
      const int& _ii,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const dvariable& pen,
      const imatrix& mflags,const imatrix& mgroup,MY_DOUBLE_TYPE s)
    {
      int i;
      //dvar_matrix& mw=(dvar_matrix&) _mw;
      dvar_vector w=rowstack(_mw);
      ivector flags=rowstack(mflags);
      ivector group=rowstack(mgroup);
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
          if (flags(kin)) set_value_exp(w(kin),x,ii,fmin,fmax,pen,s);
          for (i=mmin+1;i<=mmax;i++)
          {
            if (group(key(i))==group(key(i-1)))
            {
              //if (flags(i)) w(key(i))=w(kin);
              //  DF  1Mar19
              if (flags(key(i))) w(key(i))=w(kin);
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
      rowunstack(w,_mw);
    }

 //   void set_value(const dvar_vector& _w,const dvar_vector& x,const int& _ii,
 //     MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const dvariable& pen,const ivector& flags,
 //     const ivector& group)
 //   {
 //     dvar_vector& w=(dvar_vector&) _w;
 //     int& ii= (int&) _ii;
 //     int mmin=w.indexmin();
 //     int mmax=w.indexmax();
 //     int maxg=max(group);
 //     if (sum(flags))
 //     {
 //       cout << "flags = " << flags << endl;
 //       if (maxg)
 //       {
 //         int icount=0;
 //         ivector group_flag(1,maxg);
 //         group_flag.initialize();
 //         ii--;
 //         for (int i=mmin;i<=mmax;i++)
 //         {
 //           if (flags(i)) 
 //           {
 //             w(i)=boundp(x(ii+group(i)),fmin,fmax,pen);
 //             if (!group_flag(group(i)))
 //             {
 //               icount++;
 //               group_flag(group(i))=1;
 //             }
 //           }
 //         }
 //         ii+=icount+1;
 //         //ii+=max(group)+1;
 //       }
 //       else
 //       {
 //         for (int i=mmin;i<=mmax;i++)
 //         {
 //           if (flags(i)) w(i)=boundp(x(ii++),fmin,fmax,pen);
 //         }
 //       }
 //     }
 //   }
 //   
