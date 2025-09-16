/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#include "all.hpp"

#include <adstring.hpp>
#include <string.h>
#include <stdlib.h>
#include "safe_mem.h"

int read_newstyle_file(adstring& full_datafile_path)
{
  int new_version=0;
  ifstream ifs(full_datafile_path);
  char ch;
  adstring label;
  int version_no;

  ifs >> ch;
  if (!ifs) return new_version;
  ifs >> label;
  if (!ifs) return new_version;
  ifs >> version_no;
  if (!ifs) return new_version;

  if (ch != '#') return new_version;
  if (!(label == adstring("version")) && !(label == adstring("VERSION")) ) 
    return new_version;
  new_version = version_no;
  return new_version;
}

void error_rpt(const char *s)
{
  cerr << s << " in read_newstyle_data " << endl;
  exit(1);
}

void read_newstyle_data(cifstream& infile,newstyle_frq_data& tmpfrq)
{
  int nflags=0;
  infile >> nflags;
  if (!infile) 
  {
    error_rpt("nflags");
  }
  tmpfrq.new_frq_flags.allocate(1,nflags);
  infile >> tmpfrq.new_frq_flags;
  if (!infile) 
  {
    error_rpt("new_frq_flags");
  }
  int num_fisheries=tmpfrq.new_frq_flags(1);
   
  if(tmpfrq.new_frq_flags(2))
  {
    for (int i=1;i<=num_fisheries;i++)
    {
      line_adstring tstring;
      const line_adstring& rl = tstring;
      infile >> rl;
      tmpfrq.fishery_names+=tstring;
      if (!infile) 
      {
        error_rpt("fishery names");
      }
    }
  }
  else
  {
    for (int i=1;i<=num_fisheries;i++)
    {
      adstring("Fishery: " + str(i));
      tmpfrq.fishery_names+= adstring("Fishery: " + str(i));
      if (!infile) 
      {
        error_rpt("fishery names");
      }
    }
  }
  
  int num_fish_flags=0;
  infile >> num_fish_flags; 
  if (!infile) 
  {
    error_rpt("num_fish_flags");
  }
  tmpfrq.new_fish_flags.allocate(1,num_fisheries,1,num_fish_flags);
  infile >> tmpfrq.new_fish_flags;
  if (!infile) 
  {
    error_rpt("new_fish_flags");
  }
}

