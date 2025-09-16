/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE

#define HOME_VERSION
#include "all.hpp"

extern adstring delayed_infile;
extern adstring delayed_outfile;

adstring set_outfile_name(dvar_len_fish_stock_history& fsh);

 void dvar_len_fish_stock_history::set_Takeuchi_san_control_switches(void)
 {
   switch(++parest_flags(20))
   {
     case 1:
     {
       // Maximum number of function evaluations
       age_flags(144)=1;
       parest_flags(1)=  60;
       parest_flags(23)=1000;
       age_flags(41)=1;
       age_flags(60)=50;
       age_flags(32)=1;
       parest_flags(149)=10;
       parest_flags(153)=10;
       parest_flags(197)=0;
       if (parest_flags(142))
       { 
         age_flags(143)=(age_flags(143)+1)%3;
         cout << " age_flags(143) " << age_flags(143) << endl;
         if (age_flags(143)==1)
         { 
           parest_flags(1)=  25;
           parest_flags(142)+=10;
           cout << " parest_flags(142) " << parest_flags(142) << endl;
         }
         if (age_flags(143)==2)
         { 
           parest_flags(1)=  75;
           cout << " parest_flags(142) " << parest_flags(142) << endl;
         }
       }
 
	 int err_flag=0;
       for (int i=1;i<=num_fisheries;i++)
       {
 	fish_flags(i,1)=1;
 	fish_flags(i,48)=1;
         if (fish_flags(i,12)==0)
         {
 	  if (fish_flags(i,3)!=nage-1)
           {
             cerr << "Warning -- Error in input parameter file -- " 
                     "fish_flags(" << i << "," << 3 
                  << ") is not equal to " << nage-1 << endl;
             //err_flag=1;
           }
 	  //fish_flags(i,3)=nage-1;
         }
         else
         {
 	  if (fish_flags(i,3)!=max(1,nage-2))
           {
             cerr << "Error in input parameter file -- " 
		       "fish_flags(" << i << "," << 3
                  << ") must be set equal to " << max(1,nage-2) << endl;
             err_flag=1;
           }
 	  fish_flags(i,3)=max(1,nage-2);
         }
       }
       if (err_flag) exit(1);
 
       if (parest_flags(142)>10)
       {
         if (age_flags(143))
         {
           for (int i=1;i<=num_fisheries;i++)
           {
             // effort devs
 	    fish_flags(i,4)=1;
           }
         }
         else
	   {
           for (int i=1;i<=num_fisheries;i++)
           {
             // effort devs
 	    fish_flags(i,4)=2;
           }
         }
         // effort devs
       }
 
       // let recr vary
       age_flags(30)=1;
       // let init_pop vary
       age_flags(31)=1;
       // set initial average fishing mortality rate penalty
	 // age_flags(37)=3*sqrt(10./nage);
       if (!age_flags(37)) age_flags(37)=200;
	 //age_flags(37)=4;
       break;
     }
     case 2:
     {
       char ch;
       // Maximum number of function evaluations
       parest_flags(1)= 70; 
       age_flags(113)=1;
       age_flags(41)=1;
       // let totpop vary
       age_flags(144)=5;
       parest_flags(146)=5;
       parest_flags(149)=1;
       parest_flags(153)=1;
       age_flags(143)=0;
       cout << "Starting phase 2 ... enter a character" << endl;
       // cin >> ch;
       for (int i=1;i<=num_fisheries;i++)
       {
         // effort devs
	   fish_flags(i,4)=2;
	 }
	 // increase maximum effort deviations
	 age_flags(35)=6;
	 break;
     }
     case 3:
     {
	 parest_flags(1)=   60; //50;
       parest_flags(149)=0;
       age_flags(144)=50;
	 for (int i=1;i<=num_fisheries;i++)
	 {
	   // catchability devs
           /*
	   if (parest_flags(150)==0)
	   {
	     fish_flags(i,10)=2;
	   }
           */
	 }
	 break;
     }
     case 4:
     {
       age_flags(144)=200;
     }
     case 5:
     {
       parest_flags(1)=  60; //50;
	 age_flags(144)=500;
       parest_flags(153)=0;
       break;
     }
     case 6:
     {
       age_flags(60)=20;
       cout << "Starting phase 6 " << endl;
       parest_flags(1)=  60; //50;
       parest_flags(146)=25;
       age_flags(144)=1500;
       //age_flags(46)=20;
       //age_flags(51)=10;
       break;
     }
     case 7:
     {
       cout << "Starting phase 7 " << endl;
       parest_flags(1)=  60; //50;
       age_flags(144)=2500;
       break;
     }
     case 8:
     {
       cout << "Starting phase 8 " << endl;
       parest_flags(1)=  65; //50;
       age_flags(144)=5000;
       break;
     }
     case 9:
     {
       char ch;
       cout << "Starting phase 9 " << endl;
       parest_flags(1)=150;
       // mean length of first age class
       parest_flags(12)=1;
       // mean length of last age class
       parest_flags(13)=1;
       // increase maximum selectivity deviations
       age_flags(36)=2;
       //age_flags(51)=0;
       age_flags(144)=10000;
       break;
     }
     case 10:
     {
       char ch;
       cout << "Starting phase 10 " << endl;
       age_flags(60)=0;
       age_flags(46)=0;
       age_flags(51)=0;
       parest_flags(1)=100;
       // average standard deviations 
       parest_flags(15)=1;
       // control initial step size
       parest_flags(23)=10;
       // robust estimation procedure
       //age_flags(38)=1;
       // smaller total fishing mortality penalty
       //age_flags(39)=1;
       age_flags(144)=10000;
       delayed_outfile=set_outfile_name(*this);
       break;
     }
     case 11:
     {
       char ch;
       cout << "Starting phase 11 " << endl;
       parest_flags(1)=60; //500;
       // von Bertalanffy K           
       //parest_flags(14)=1;
       break;
     }
     case 12:
     {
       char ch;
       cout << "Starting phase 12 " << endl;
       parest_flags(1)=60; //500;
       //delayed_outfile="10011000.fit";
       delayed_outfile=set_outfile_name(*this);
       // end constraints on overall average fishing mortality
       //age_flags(37)=0;
       if (parest_flags(150)!=0)
       {
         for (int i=1;i<=num_fisheries;i++)
         {
           fish_flags(i,5)=2;
         }
       }
       // parest_flags(199)=1; //500;
       break;
     }
     case 13:
     {
       char ch;
       parest_flags(1)=70; //500;
       delayed_infile="00011000.fit";
       cout << "Starting phase 13 " << endl;
       // end constraints on overall average fishing mortality
       age_flags(37)=0;
       //delayed_outfile="01011000.fit";
       delayed_outfile=set_outfile_name(*this);
       parest_flags(199)=0; //500;
       break;
     }
     case 14:
     {
       char ch;
       cout << "Starting phase 14 " << endl;
       parest_flags(1)=100; //500;
       age_flags(50)=1;
       delayed_outfile=set_outfile_name(*this);
       break;
     }
     case 15:
     {
       char ch;
       parest_flags(1)=100;
       cout << "Starting phase 15 " << endl;
       // estimate Von Bertananffy  K
       //parest_flags(14)=1;
       delayed_outfile=set_outfile_name(*this);
       age_flags(50)=1;
       break;
     }
     case 16:
     {
       char ch;
       cout << "Starting phase 16 " << endl;
       parest_flags(1)=100; //500;
       // end constraints on overall average fishing mortality
       //age_flags(37)=0;
       // von Bertalanffy K           
       parest_flags(30)=0;
       delayed_outfile=set_outfile_name(*this);
       age_flags(50)=1;
       break;
     }
     default:
     {
       cerr << "switch 20 in parameter relaxation sequence"
             " out of bounds" << endl;
       exit(1);
     }
   }
   if (parest_flags(142))
   {
     parest_flags(20)--;
   }
 } 
