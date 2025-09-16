/*
Copyright (C) Pacific Community (SPC)
Refer to COPYRIGHT and LICENSE files in the root of the repository
*/
#define MY_REAL_DOUBLE double
#include "extras.h"
#define USE_NO_LAPACKE
#if !defined(__ADMSTHREAD__)
#  define __ADMSTHREAD__
class ad_master_slave_synchro
{
public:
  int threadid;
  int flag;
  pthread_mutex_t  mutex;
  pthread_cond_t   cond1;
  pthread_cond_t   cond2;
  ad_master_slave_synchro(void);
  void master_lock(void); 
  void slave_lock(void); 
  void master_send_signal(void); 
  void master_get_signal(void); 
  void slave_get_signal(int i); 
  void slave_send_signal(void); 
};

#endif   // #if !defined(__ADMSTHREAD__)
