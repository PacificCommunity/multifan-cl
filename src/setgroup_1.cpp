/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
  group_manager_1::~group_manager_1(void) 
  {
    std::cout << "Calling ~group_manager_1"  << std::endl;
  }

  group_manager_1::group_manager_1(void) : maxg(0), groupsum(0) {}

  group_manager_1::group_manager_1(ivector & ff4,imatrix& _mflags,
    const ivector& _group) : active_flags(ff4),
    mflags(_mflags),
    group(_group)
  {
    int mmin=group.indexmin();
    int mmax=group.indexmax();
    active_group.allocate(mmin,mmax);
    active_group=group;
    /*
    for (int i=mmin;i<=mmax;i++)
    {
      if (active_flags(i)>0) 
      {
        active_group(i)=group(i);
      }
      else
      {
        active_group(i)=0;
      }
    }
    */
    key.allocate(mmin,mmax);
    sort(active_group,key); 
    maxg=max(active_group);
    groupsum=sum(group);
  }


  void set_value(const dvar_vector& _v,const dvar_vector& x,const int& _ii,
    ivector& mflags,MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const prevariable& pen,MY_DOUBLE_TYPE s);

  int size(const dvar_vector& v,const ivector& flags);

  void set_value_inv(const dvar_matrix & _w,const dvector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,MY_DOUBLE_TYPE s,
    std::shared_ptr<group_manager_1>& pgroup_manager_1)
  {
    ADUNCONST(dvar_matrix,w)  //NMD_27apr2018
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    const int groupsum=pgroup_manager_1->groupsum;
    imatrix& mflags=(imatrix&)(pgroup_manager_1->mflags);
    const ivector& flags=pgroup_manager_1->active_flags;
    if (groupsum)
    {
      ivector& key=pgroup_manager_1->key;
      const ivector& flags=pgroup_manager_1->active_flags;
      ivector& group=pgroup_manager_1->active_group;
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
 

  void set_value(const dvar_matrix & _w,const dvar_vector& x,const int& ii,
    MY_DOUBLE_TYPE fmin,MY_DOUBLE_TYPE fmax,const prevariable& pen,
    MY_DOUBLE_TYPE s, std::shared_ptr<group_manager_1>& pgroup_manager_1)
  {
    ADUNCONST(dvar_matrix,w)  //NMD_27apr2018
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    const int groupsum=pgroup_manager_1->groupsum;
    imatrix& mflags=(imatrix&)(pgroup_manager_1->mflags);
    const ivector& flags=pgroup_manager_1->active_flags;
    ivector& group=pgroup_manager_1->active_group;
    if (groupsum)
    {
      ivector& key=pgroup_manager_1->key;
      const ivector& flags=pgroup_manager_1->active_flags;
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

  
  int num_active(const dvar_matrix& w,ivector& flags,imatrix& mflags,
    const ivector& group,std::shared_ptr<group_manager_1>& pgroup_manager_1,
    const ivector & oldfflags)
  {
    int mmin=w.indexmin();
    int mmax=w.indexmax();
    if (mmin != group.indexmin() || mmax != group.indexmax())
    {
      cerr << "G incompatible array bounds in num_active() " << endl;
      exit(1);
    }
    if (pgroup_manager_1.get()!=0)
    {
      // group_manager_1 is already allocated so deallocate it with reset()
      pgroup_manager_1.reset();
    }
      pgroup_manager_1=std::make_shared<group_manager_1>  \
        (flags,mflags,group);
    int nv=0;

    if (pgroup_manager_1->maxg)
    {
      int i=mmin;
      ivector& key=pgroup_manager_1->key;
      const ivector& active_flags=pgroup_manager_1->active_flags;
        
      int kin=key(i);
      if (pgroup_manager_1->active_flags(kin))
      {
        nv+=size(w(kin),mflags(kin));
      }
      for (int i=mmin+1;i<=mmax;i++)
      {
        if (group(key(i))!=group(key(i-1)))
        {
          kin=key(i);
          if (pgroup_manager_1->active_flags(key(i)))
            nv+=size(w(kin),mflags(kin));
        } 
      }
    }
    else
    {
      for (int i=mmin;i<=mmax;i++)
      {
        if (pgroup_manager_1->active_flags(i))
          nv+=size_count(w(i));
      }
    }
    return nv;
  }
  
