/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include "all.hpp"
void ad_eat_white_space(istream & s)
{
  int comment_flag=0;
  int ch;
  while (1)
  {
    ch=s.get();
    if (ch==-1) 
    {
      break;
    }
    if (s.eof()) break;
    if (ch=='#') comment_flag=1;
    if (ch == '\n' || ch == '\r') comment_flag=0;
    if (!comment_flag  && 
      (
          isdigit(ch) || 
          isalpha(ch) || 
          ch=='-' || 
          ch=='+' || 
          ch=='.' || 
          ch=='e' || 
          ch=='E'  
      )
      ) break;
  }
  #ifdef __MSVC32__
  s.putback(ch);
  #else
  // s.putback(ch);
  s.unget();
  #endif
//#endif
}
void strip_microsoft_junk(const adstring& _s)
{
  adstring & s = (adstring &) (_s);
  int n=strlen(s);
  if (s[n]=='\r')
  {
    s[n]='\0';
    cout << "stripped microsoft junk" << endl;
  }
}

#if !defined(__BORLANDC__)
void read_option_file(const adstring& full_input_option_path,
  imatrix& cls)
{
  mytestcin();
  const int MAX_NSWITCH=3000;
  imatrix tmp(1,MAX_NSWITCH,1,5);
  tmp.initialize();
  int ii=0;
  //cifstream* pcif = NULL;
  // If input path is regular file, connect it to cin and
  // treat it as if standard input.  --- PK: 2002-4-10. 

  strip_microsoft_junk(full_input_option_path);

  mytestcin();
  if (!(full_input_option_path==adstring("-")))
  {
    cifstream cif((const char*)(full_input_option_path));
    if (!cif) {
      cerr << " Error trying to open command line file "
           << full_input_option_path << endl;
      exit(1);
    }
  
    while(1) {
      if (cif.eof()) break;
      if (ii++>=MAX_NSWITCH){
        cerr << "Maximum number of command line switches of "
             << MAX_NSWITCH << " exceeded" << endl;
        exit(1);
      }
      ad_eat_white_space(cif);
      adstring tmpstring;
      cif >> tmpstring;
      if (cif.fail() || cif.eof()) {
        ii--;
        break;
      }
      if (isalpha(tmpstring(1)))
      {
        // multi-species situation
        if (tmpstring=="sp")
        {
          tmp(ii,4)=1;
          cif >> tmp(ii,5);
          cif >> tmp(ii,1);
        }
        else
        {
          cerr << "unknown option in option file" << endl;
          ad_exit(1);
        }
      }
      else
      {
        tmp(ii,1)=atoi(tmpstring);
      }
          
      if (cif.eof()) {
        ii--;
        break;
      }
      ad_eat_white_space(cif);
      cif >> tmp(ii,2);
      ad_eat_white_space(cif);
      cif >> tmp(ii,3);
      if (!cif) {
        cerr << " Error trying to read option line " << ii 
           << " from script file " << endl;
        ad_exit(1);
      }
    }
    int nswitch=ii;
    //if (allocated(cls)) cls.deallocate();
    cls.allocate(1,nswitch,1,5);
    for (int i=1;i<=nswitch;i++)
    {
      cls(i)=tmp(i);
      cout <<"optfile.cpp " << cls(i) << endl;
    }
    cout << "number of options was " << nswitch << endl;
  }
  else
  {
    mytestcin();
    while(1) {
      if (cin.eof()) break;
      if (ii++>=MAX_NSWITCH){
        cerr << "Maximum number of command line switches of "
             << MAX_NSWITCH << " exceeded" << endl;
        exit(1);
      }
      ad_eat_white_space(cin);
      adstring tmpstring;
      cin >> tmpstring;
      if (cin.fail() || cin.eof()) {
        ii--;
        break;
      }
      if (isalpha(tmpstring(1)))
      {
        // multi-species situation
        if (tmpstring=="sp")
        {
          tmp(ii,4)=1;
          cin >> tmp(ii,5);
          cin >> tmp(ii,1);
        }
        else
        {
          cerr << "unknown option in option file" << endl;
          ad_exit(1);
        }
      }
      else
      {
        tmp(ii,1)=atoi(tmpstring);
      }
          
    mytestcin();
      if (cin.eof()) {
        ii--;
        break;
      }
      ad_eat_white_space(cin);
      cin >> tmp(ii,2);
      ad_eat_white_space(cin);
      cin >> tmp(ii,3);
      if (!cin) {
        cerr << " Error trying to read option line " << ii 
           << " from script file " << endl;
        ad_exit(1);
      }
    }
    int nswitch=ii;
    //if (allocated(cls)) cls.deallocate();
    cls.allocate(1,nswitch,1,5);
    for (int i=1;i<=nswitch;i++)
    {
      cls(i)=tmp(i);
      cout <<"optfile.cpp " << cls(i) << endl;
    }
    cout << "number of options was " << nswitch << endl;
  }
  cin.clear();
}
#else
void read_option_file(const adstring& full_input_option_path,
  imatrix& cls)
{
  mytestcin();
  const int MAX_NSWITCH=3000;
  imatrix tmp(1,MAX_NSWITCH,1,3);
  int ii=0;
  if (full_input_option_path==adstring("-"))
  {
    while(1) {
      if (cin.eof()) break;
      if (ii++>MAX_NSWITCH){
        cerr << "Maximum number of command line switches of "
             << MAX_NSWITCH << " exceeded" << endl;
        exit(1);
      }
      ad_eat_white_space(cin);
      cin >> tmp(ii,1); 
      if (cin.eof()) {
        ii--;
        break;
      }
      ad_eat_white_space(cin);
      cin >> tmp(ii,2);
      ad_eat_white_space(cin);
      cin >> tmp(ii,3);
      if (!cin) {
        cerr << " Error trying to read option line " << ii 
           << " from script file " << endl;
        exit(1);
      }
    }
    int nswitch=ii;
    //if (allocated(cls)) cls.deallocate();
    cls.allocate(1,nswitch,1,3);
    for (int i=1;i<=nswitch;i++)
    {
      cls(i)=tmp(i);
      cout <<"optfile.cpp " << cls(i) << endl;
    }
    cout << "number of options was " << nswitch << endl;
  }
  else
  {
    cifstream cif((char*)(full_input_option_path));
    int nswitch;
  
    if (!cif) {
      cerr << " Error trying to open command line file "
           << full_input_option_path << endl;
      exit(1);
    }
    ii=0;
    MY_DOUBLE_TYPE xtmp=0.0;
    while(1)
    {
      if (!cif) {
        cerr << " Error trying to read number of options from "
             << full_input_option_path << endl;
        exit(1);
      }
      ad_eat_white_space(cif);
      cif >> xtmp;
      if (cif.eof()) break;
      if (ii++>=MAX_NSWITCH)
      {
        cerr << "Maximum number of command line switches of "
             << MAX_NSWITCH << " exceeded" << endl;
        exit(1);
      }
      
      tmp(ii,1)=xtmp;
      ad_eat_white_space(cif);
      cif >>  tmp(ii,2);
      if (!cif) {
        cerr << " Error trying to read option line " << ii 
           << " from full_input_option_path " << endl;
        exit(1);
      }
      ad_eat_white_space(cif);
      cif  >> tmp(ii,3);
      if (!cif) {
        cerr << " Error trying to read option line " << ii 
           << " from full_input_option_path " << endl;
        exit(1);
      }
    }
    nswitch=ii;
    cls.allocate(1,nswitch,1,3);
    for (int i=1;i<=nswitch;i++)
    {
      cls(i)=tmp(i);
      //cout <<"optfile.cpp " << cls(i) << endl;
    }
  }
  mytestcin();
  cout <<"optfile.cpp " << cls << endl;
}
#endif
