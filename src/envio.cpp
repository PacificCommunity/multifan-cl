/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE


#include "all.hpp"

   extern int _file_version_no;

void read_environmental_data(adstring& full_datafile_path,
  dvar_len_fish_stock_history& fsh)
{
  if (fsh.age_flags(72))
  {
    int p=0; 
    adstring root;
    p=pos(".",full_datafile_path);
    if (p>1)
    {
      root=full_datafile_path(1,p-1);
    }
    else
    {
      root=full_datafile_path; 
    }
    adstring envfilename=root + adstring(".env");
    cout <<"envio.cpp " << envfilename << endl;
    cifstream cif((char*)(envfilename));
    if (!cif)
    {
      cerr << "Error trying to open environmental data file " << (char*) envfilename
           << endl;
      exit(1);
    }
    int ny;
    cif >> ny;
    if (ny != fsh.rec_covars.indexmax())
    {
      cerr << "size of environmental time series in file "
           << envfilename << " is wrong " << endl;
      exit(1);
    }
    
    cif >> fsh.rec_covars;
    if (!cif)
    {
      cerr << "Error trying to read environmental data from file " << (char*) envfilename
           << endl;
      exit(1);
    }
    cout << "In *env*cpp " << endl << fsh.rec_covars << endl;
    fsh.rec_covars-=mean(fsh.rec_covars);
    fsh.rec_covars/=norm(fsh.rec_covars);
    cout << "In *env*cpp " << endl << fsh.rec_covars << endl;
  }
}
#undef HOME_VERSION

