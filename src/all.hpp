/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#ifndef _ALL_HPP
# define AD_FPEXCEPTION
  #define _ALL_HPP
 class par_uostream;
 class fishery_catch_at_age_record_array;
#define HOME_VERSION
#if defined(__WIN32)
#  include <windows.h>
#endif
#include <fvar.hpp>
#include <ciso646>
#include "adthread.h"
#include "safe_mem.h"
 #ifdef __BORLANDC__
  #define catch Catch
 #endif
#define __SAFE_DEALLOCATE__(df) \
    if (allocated(df)) df.deallocate();
#include "newmult.hpp"
//#include <tt.hpp>
#include <cifstrem.h>
#include <stdlib.h>
#include <memory>
#define catch Catch
//#include <iostream>

#include "ils_qr.hpp"
#include "variable.hpp"
#include "globals.h"
#ifdef __MSVC32__
#  include <direct.h>
#endif
#ifdef __BORLANDC__
#  include <dir.h>
#endif
#if defined(linux)
#  include <unistd.h>
#endif

  #if defined(__MSVC32__) || defined(__ADMING__)
    #define lseek _lseek
//    #define  read _read
//    #define write _write 
//    #define open _open
//    #define close _close 
    #define VCCLOSE close 
  #endif

#if !defined(__BORLANDC__)
  int heapcheck(void);
#endif
extern ofstream check_file;
//  #define HERE 
//#ifdef HERE
 //#undef HERE
 //#define HERE cout <<"reached line "<< __LINE__ << " in " << __FILE__ << endl;
//#endif
#define COUT_TRACE(object) cout <<  #object " = " << endl << object << endl;
#define SHOW_ALL(object) check_file  << endl <<  #object " = " << endl << setfixed << setprecision(8) << object << endl;
#define VSHOW_ALL(object) vcheck_file  << endl <<  #object " = " << endl << setfixed << setprecision(8) << object << endl;
#include "makebig2.hpp"

#endif
