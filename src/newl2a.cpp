/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#define HOME_VERSION
#include "all.hpp"

   void check_index(int index_min,int index_max,int i,char * s)
   {
     if (i<index_min)
     {
       cerr << "Index too low in " << s
       << " index_value was "<< i << "Minimum valid index is "
       << index_min<<endl;
       exit(1);
     }
     if (i>index_max)
     {
       cerr << "Index too high in " << s
       << " index_value was "<< i << "Maximum valid index is "
       << index_max<<endl;
       exit(1);
     }
   }

#undef HOME_VERSION

