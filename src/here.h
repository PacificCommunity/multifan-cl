/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#if !(defined HERE_H)
# define HERE_H

#ifdef HERE
 #undef HERE
#endif
#define HERE cout << "reached line " << __LINE__ << " in " << __FILE__ << endl; 

#endif
