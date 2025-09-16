/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#include <admodel.h>


 int option_match(int argc,char * argv[], char * string)
 {
   int rval=-1;
   for (int i=0;i<argc;i++)
   {
     if (!strcmp(argv[i],string))
     {
       rval=i;
       break;
     }
   }
   return rval;
 }

 int option_match(int argc,char * argv[], char * string, const int& _nopt)
 { 
   int & nopt=(int&)(_nopt);
   int iflag=0;
   if (!strcmp(string,"-file"))
     iflag=0;
   int rval=-1;
   int i;
   nopt=0;
   for (i=0;i<argc;i++)
   {
     if (iflag)
     {
        cout <<"optmatch.cpp " << argv[i] << endl;
        cout << "i = " << i << endl;
        cout <<"optmatch.cpp " << string << endl;
        cout << "strcmp = " << strcmp(argv[i],string) << endl << endl;
     }
     if (!strcmp(argv[i],string))
     {
       rval=i;
       break;
     }
   }
   for (i=rval+1;i<argc;i++)
   {
     if (argv[i][0] == '-') break;
     nopt++;
   }
   if (iflag)
     cout <<"optmatch.cpp " << rval << endl;
   return rval;
 }


#undef HOME_VERSION
