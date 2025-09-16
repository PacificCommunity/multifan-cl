/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"
#pragma hdrstop

#ifdef __BCPLUSPLUS__
  #include <iomanip.h>
#endif
#ifdef __ZTC__
  #include <iomanip.hpp>
  #include <fstream.hpp>
  #include <strstrea.hpp>
#endif
#ifdef __NDPX__
  #include <fstream.h>
  #include <sstream.h>
  extern "C" {
  #include <dos.h>
  }
#endif

    void fishery_catch_at_age_record_array::
	     allocate(fishery_header_record_array& fhra,int ng)
    {
      for (int i=indexmin();i<=indexmax();i++)
      {
	(*this)(i).allocate(fhra(i),ng);
      }
    }

    fishery_catch_at_age_record_array::fishery_catch_at_age_record_array
    (fishery_header_record_array& fhra,int ng)
    {
      nage=ng;
      index_min=fhra.indexmin();
      index_max=fhra.indexmax();
      int sz=size();
      ptr=new fishery_catch_at_age_record [sz];
      if (ptr==NULL)
      {
        cerr << "Error allocating memory in fishery_freq_record_array"<<endl;
        exit(21);
      }
      ptr-=indexmin();
      allocate(fhra,ng);
    }

    ostream& operator << (ostream& ofs,header_record& fhr)
    {
      ofs << "#     year          month      week\n#" << endl;
      ofs << "      " << fhr.year  << "          " << fhr.month
          << "          " << fhr.week << endl;
      return ofs;
    }

    ostream& operator << (ostream& ofs,fishery_header_record& fhr)
    {
      ofs << * (header_record *) (&fhr);
      ofs << "#  fishery      catch       effort" << endl;
      ofs << "       " << fhr.fishery << "       " << setprecision(4)
	  << fhr.obs_tot_catch  << "       " << setprecision(4)
          << fhr.effort << endl;
      return ofs;
    }

    istream& operator >> (istream& ifs,header_record& fhr)
    {
      ifs >> fhr.year >> fhr.month >> fhr.week;
      return ifs;
    }

    istream& operator >> (istream& ifs,fishery_header_record& fhr)
    {
      ifs >> * (header_record *) (&fhr);
      ifs >> fhr.fishery >> fhr.obs_tot_catch >> fhr.effort;
      return ifs;
    }

    cifstream& operator >> (cifstream& ifs,header_record& fhr)
    {
      ifs >> fhr.year >> fhr.month >> fhr.week;
      return ifs;
    }

    cifstream& operator >> (cifstream& ifs,fishery_header_record& fhr)
    {
      ifs >> * (header_record *) (&fhr);
      ifs >> fhr.fishery;
      if (allocated(fhr.species) && _data_file_version_no>6)
        ifs >> fhr.species;
      if (allocated(fhr.species1) && _data_file_version_no>7)
        ifs >> fhr.species1;
      if (allocated(fhr.species2) && _data_file_version_no>8)
        ifs >> fhr.species2;
      if (allocated(fhr.species3) && _data_file_version_no>8)
        ifs >> fhr.species3;
      ifs >> fhr.obs_tot_catch >> fhr.effort;

      if (_data_file_version_no>5)
        ifs >> fhr.effort_weight;
      else
        fhr.effort_weight=1.0;

      //if (fhr.obs_tot_catch==0)
      //{
      //  cout <<   " fhr.fishery " << fhr.fishery 
      //       << " fhr.obs_tot_catch " << fhr.obs_tot_catch 
      //       << " fhr.effort " << fhr.effort << endl;
      //  fhr.obs_tot_catch =-1.0;
      //}
      return ifs;
    }


ostream& operator << (ostream& ofs,
                      fishery_catch_at_age_record_array& fcara)
{
  for (int i=fcara.indexmin();i<=fcara.indexmax();i++)
  {
    ofs << fcara(i);
  }
  return ofs;
}

istream& operator >> (istream& ifs,
                      fishery_catch_at_age_record_array& fcara)
{
  for (int i=fcara.indexmin();i<=fcara.indexmax();i++)
  {
    ifs >> fcara(i);
  }
  return ifs;
}

cifstream& operator >> (cifstream& ifs,
                      fishery_catch_at_age_record_array& fcara)
{
  for (int i=fcara.indexmin();i<=fcara.indexmax();i++)
  {
    ifs >> fcara(i);
  }
  return ifs;
}

    ostream& operator <<
           (ostream& ofs,fishery_catch_at_age_record& fhr)
    {
      ofs << * (fishery_header_record *) (&fhr);
      ofs << "#\n#  The proportions at age in the catch\n#" << endl;
      ofs <<  setprecision(3) << setw(4) << fhr.props << endl;
      return ofs;
    }

    istream& operator >>
           (istream& ifs,fishery_catch_at_age_record& fhr)
    {
      ifs >> * (fishery_header_record *) (&fhr);
      ifs >> fhr.props;
      return ifs;
    }

    cifstream& operator >>
           (cifstream& ifs,fishery_catch_at_age_record& fhr)
    {
      ifs >> * (fishery_header_record *) (&fhr);
      ifs >> fhr.props;
      return ifs;
    }


#undef HOME_VERSION
