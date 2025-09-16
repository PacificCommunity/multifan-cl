/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"

int greater_than(fishery_freq_record **tmp,int j, ivector& fr);

  void print(fishery_freq_record_array& flfra1);
    //int relaxation_sequence_control(void);
  void sort(fishery_freq_record_array& flfra1,ivector& fr,
    fishery_freq_record_array& flfra,ivector& regmin,ivector& regmax)
  {
    int num_regions=max(fr);
    int mmax=flfra.indexmax();
    int i;
    for (i=1;i<=mmax;i++)
    {
      flfra[i]=flfra1[i];
    }
    fishery_freq_record **tmp=new fishery_freq_record * [ mmax];
    // *** I like to use 1 offset pointers rather than 0 offset pointers
    tmp-=1;
    fishery_freq_record *rec;
    for (i=1;i<=mmax;i++)
    {
      tmp[i]=flfra.ptr+i;
    }
    for (i=1;i<=mmax-1;i++)
    {
      for (int j=1;j<=mmax-i;j++)
      {
        if (greater_than(tmp,j,fr ))
        {
          rec= tmp[j];
          tmp[j]=tmp[j+1];
          tmp[j+1]=rec;
        }
      }
    }
    int ii=1;
    regmin(ii)=1;
    for (i=2;i<=mmax-1;i++)
    {
      if ( fr(tmp[i+1]->fishery)> fr(tmp[i]->fishery))
      {
        regmax(ii++)=i;
        if (ii<=num_regions)
        {
          regmin(ii)=i+1;
        }
      }
    }
    regmax(num_regions)=mmax;
    //cout << "regmin " << regmin << endl;
    //cout << "regmax " << regmax << endl;
    for (i=1;i<=mmax;i++)
    {
      flfra1[i]=*(tmp[i]);
    }
  //  cout << "flfra1 " << endl << flfra1 << endl;
    delete [] (tmp+1);
  }

static void junk(int i1,int i2){;}

int greater_than(fishery_freq_record **tmp,int j, ivector& fr)
{
  int i1=tmp[j]->fishery;
  int i2=tmp[j+1]->fishery;
  junk(i1,i2);
  if ( fr(tmp[j]->fishery)< fr(tmp[j+1]->fishery))
  {
    return 0;
  }
  else if ( fr(tmp[j]->fishery)> fr(tmp[j+1]->fishery))
  {
    return 1;
  }
  else if ( tmp[j]->year < tmp[j+1]->year)
  {
    return 0;
  }
  else if ( tmp[j]->year > tmp[j+1]->year)
  {
    return 1;
  }
  else if ( tmp[j]->month < tmp[j+1]->month)
  {
    return 0;
  }
  else if ( tmp[j]->month > tmp[j+1]->month)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}
#undef HOME_VERSION

