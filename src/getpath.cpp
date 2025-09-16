/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
// read the current path into adstring current_path and the
// root name for the hyphothesis into the adstring root
//
#include <fvar.hpp>

#include <adstring.hpp>
void getpaths(adstring& current_path,adstring& root)
{ 
  {
    ifstream ifs("do_sd.tmp");
    ifs >> current_path;
    if (!ifs)
    {
      cerr << "Error reading current_path in getpaths()"
           << endl;
      exit(1);
    }
  }
  adstring parfilename;
  {
    adstring fname= current_path + adstring("\\") + adstring("file_nam.tmp");
    ifstream ifs(fname);
    if (!ifs)
    {
      cerr << "Error trying to open file " << fname << endl;
    }
    ifs >> parfilename;
    cout << "Line 1 of file+nam.tmp" << endl;
    cout <<"getpath.cpp " << parfilename << endl;
    ifs >> parfilename;
    cout << "Line 1 of file+nam.tmp" << endl;
    cout <<"getpath.cpp " << parfilename << endl;
    if (!ifs)
    {
      cerr << "Error trying to read parfilename from file " << fname << endl;
    }
    int llen=length(parfilename);
    int p=llen;
    for (int i=1;i<=llen;i++)
    {
      if (parfilename(i)=='.')
      {
        p=i;
        break;
      }
    } 
    cout << " p = " << p << endl; 
    //p=pos(".",parfilename);
    if (p>1)
    {
      root=parfilename(1,p-1);
    }
    else
    {
      cerr << "Root of parfilename is empty in routine getpath()" << endl;
    }
  }
}
